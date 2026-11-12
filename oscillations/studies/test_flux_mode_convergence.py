#!/usr/bin/env python3
"""v4.3.1 full-P8 mode convergence. Does not change the standard fitter."""
import argparse
import csv
import os
from pathlib import Path
import sys

import numpy as np

from compare_normalized_flux_profile import TAGS, mean_center, profile, root_inputs, vec, write_csv
from test_electron_holdouts import fit_nested
from test_electron_mixing import MAX, PARAMS, parameter_row

MASKED_CASES = ("ratios_dropped", "only_ratios")


def current_root_inputs(args):
    osc = Path(__file__).resolve().parents[1]
    sys.path[:0] = [str(osc), str(osc.parent)]
    os.environ.setdefault("CCNUEROOT", str(osc.parent))
    sys.argv = [sys.argv[0]]
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    from fitData_normalizedFluxPrior import NormalizedFluxPrior
    from tools.Helper import GetSliceIndices
    from tools.StitchedHistogram import StitchedHistogram

    for case in args.cases:
        tag = TAGS["ratio_all" if case in MASKED_CASES else case] + args.tag_suffix
        cfg = osc / f"HIST_CONFIG_{tag}.json"
        path = osc / "rootfiles" / f"NuE_stitched_hists_{tag}.root"
        for source in (cfg, path):
            if not source.is_file():
                raise FileNotFoundError(source)
        h = StitchedHistogram("sample", include_shared_ratio_stats=True)
        h.SetHistConfig(str(cfg))
        h.Load(str(path))
        dh = h.GetDataHistogram()
        d, m = vec(dh), vec(h.GetMCHistogram())
        C = h.GetCovarianceMatrix(True)
        A = np.asarray(h.GetAMatrix(), dtype=float)
        band = dh.GetVertErrorBand("Flux")
        if A.shape != (1000, 54) or len(d) != 54 or not band or band.GetNHists() != len(A):
            raise ValueError(f"{case}: expected aligned 1,000-universe Flux bands across 54 bins")
        A = np.array([vec(band.GetHist(i)) - d for i in range(len(A))]) - A
        prior = NormalizedFluxPrior(C, A)
        prior.assert_flux_covariance(h.GetCovarianceMatrix(False) - C)
        if prior.K != len(d):
            raise ValueError(f"{case}: expected 54 numerical modes, found {prior.K}")
        keep = np.arange(len(d))
        profile_keep = keep
        first_ratio_bin = None
        if case != "direct_all":
            ratio = vec(h.ratio_id)
            if not np.isin(ratio, [0., 1.]).all() or np.count_nonzero(ratio) != 24:
                raise ValueError(f"{case}: expected a binary mask with 24 ratio bins")
            if case == "ratios_dropped":
                profile_keep = keep[ratio == 0]
                if args.sterile_k is None:
                    keep = profile_keep
            elif case == "only_ratios":
                keep = keep[ratio == 1]
                profile_keep = keep
        if case in ("ratio_all", "only_ratios") or (
            case == "ratios_dropped" and args.sterile_k is not None
        ):
            bins = GetSliceIndices(str(cfg), "", ["fhc_ratio"])
            if not bins or not ratio[bins[0]]:
                raise ValueError("Could not identify the first FHC ratio bin")
            first_ratio_bin = int(np.flatnonzero(keep == bins[0])[0])
        print(f"{case}: ROOT covariance closure passed for {path}; "
              f"profiling {len(profile_keep)} bins, scoring {len(keep)} bins", flush=True)
        yield case, A[:, keep], C[np.ix_(keep, keep)], d[keep], m[keep], h, first_ratio_bin, keep, profile_keep


def mode_basis(A, C):
    """Order flux modes by variance relative to the non-flux covariance."""
    L = np.linalg.cholesky(C)
    B = mean_center(A) / np.sqrt(len(A))
    U, s, _ = np.linalg.svd(np.linalg.solve(L, B.T), full_matrices=False)
    F = s[:, None] * (L @ U).T
    np.testing.assert_allclose(F.T @ F, B.T @ B, rtol=1e-8, atol=1e-8)
    return L, U, s*s, F


def scores_by_k(r, L, U, eigen):
    """Return scores for K=0,...,full, using the whitened eigencoordinates."""
    r = np.linalg.solve(L, r)
    projection = (U.T @ r)**2
    # All dropped modes contribute their full residual power.
    loss = projection * eigen / (1+eigen)
    full = np.sum(projection/(1+eigen))
    return full + np.r_[np.cumsum(loss[::-1])[::-1], 0.]


def fixed_checks(case, A, C, d, m, predictions, tolerance):
    L, U, eigen, F = mode_basis(A, C)
    kmax = len(F)
    q0 = scores_by_k(d-m, L, U, eigen)
    rows = []
    for k, eps in enumerate(np.r_[eigen, 0.]):
        null_excess = q0[k]-q0[-1]
        # For every mass, 0 <= min chi2_K - min chi2_full <= eps*q0_full,
        # since zero mixing is feasible at every mass.
        bound = eps*q0[-1]
        delta_bound = max(null_excess, abs(null_excess-bound))
        rows.append(dict(case=case, k=k, next_eigenvalue=eps, null_chi2=q0[k],
                         null_excess=null_excess, profiled_chi2_error_bound=bound,
                         null_relative_profile_error_bound=delta_bound,
                         certified_within_tolerance=bound <= tolerance,
                         tolerance=tolerance))
    points = []
    for name, mc in predictions.items():
        values = scores_by_k(d-mc, L, U, eigen)
        _, _, _, full_shift = profile(C, F, 1., d-mc)
        for k, score in enumerate(values):
            _, _, _, shift = profile(C, F[:k], 1., d-mc)
            np.testing.assert_allclose(score, (d-mc) @ np.linalg.solve(C+F[:k].T@F[:k],d-mc),
                                       rtol=1e-9, atol=1e-9)
            points.append(dict(case=case, k=k, point=name, chi2=score, delta_vs_null=q0[k]-score,
                               score_excess=score-values[-1],
                               flux_shift_difference_nonflux_norm=np.linalg.norm(
                                   np.linalg.solve(L, shift-full_shift))))
    certified = next(row["k"] for row in rows if row["certified_within_tolerance"])
    print(f"{case}: full null={q0[-1]:.6f}; first K certified at {tolerance:g} chi2 is "
          f"{certified}/{kmax}", flush=True)
    return rows, points, F, certified


def mass_profile(score, best, reference, masses, full_points=None, report=None):
    """Constrained multistart fit of the three mixings at each fixed mass."""
    from scipy.optimize import LinearConstraint, minimize

    cons = LinearConstraint([[1., 1., 1.]], -np.inf, 1.)
    last = np.zeros(3)
    output = []
    for i, mass in enumerate(masses):
        def objective(u):
            return score(np.r_[mass/100., u])
        tries = []
        starts = [last, best[1:], reference[1:], [.08, .08, .0], [.02, .02, .50], [.2, .4, 0.]]
        if full_points is not None:
            starts.append(full_points[i][1][1:])
        for start in starts:
            start = np.asarray(start)
            tries.append((objective(start), start))
            opt = minimize(objective, start, method="SLSQP", bounds=list(zip(np.zeros(3), MAX[1:])),
                           constraints=(cons,), options={"ftol": 1e-8, "maxiter": 120})
            if np.isfinite(opt.fun) and np.all(opt.x >= -1e-9) and np.all(opt.x <= MAX[1:]+1e-9) and opt.x.sum() <= 1+1e-9:
                x = opt.x.clip(0., MAX[1:])
                x /= max(1., x.sum())
                tries.append((objective(x), x))
        value, last = min(tries, key=lambda entry: entry[0])
        output.append((value, np.r_[mass/100., last]))
        if report and ((i+1) % 20 == 0 or i+1 == len(masses)):
            report(f"fixed-mass profile {i+1}/{len(masses)}")
    return output


def root_study(case, A, C, d, m, templates, h, refs, args, first_ratio_bin=None, keep=None):
    from scipy.linalg import cho_factor, cho_solve

    ref = refs.get(case)
    if ref is None:
        # Seed the search near the current full-basis solution. All K values are refitted.
        mixing = [.014, .039, 0.] if case == "direct_all" else [.026, .044, 0.]
        x_ref = np.array([.1927, *mixing])
    else:
        x_ref = np.array([float(ref[name]) for name in PARAMS])
        x_ref[0] /= 100.
    if keep is None:
        keep = np.arange(len(d))

    def prediction(x):
        h.OscillateHistogram(100*x[0], *x[1:], False, False)
        return vec(h.GetOscillatedHistogram())[keep]

    np.testing.assert_allclose(prediction(np.zeros(4)), m, rtol=1e-7, atol=1e-8)
    rows, points, F, certified = fixed_checks(
        case, A, C, d, m, {"null": m, "saved_best": prediction(x_ref), **templates}, args.tolerance)
    kmax = len(F)
    requested = set(args.ks)
    if not args.current_root:
        requested |= {max(0, certified-1), certified}
    ks = [kmax] + sorted(requested - {kmax})
    if max(ks) > kmax:
        raise ValueError(f"{case}: cannot retain more than {kmax} modes.")
    masses = np.unique(np.r_[np.linspace(0, 100, args.mass_steps), 100*x_ref[0], 19.26, 96.21])
    fits, curves, spectra, ratio_rows, pulls = [], [], [], [], []
    full_points = None
    for k in ks:
        factor = cho_factor(C+F[:k].T@F[:k])
        def score(x):
            r = d-prediction(x)
            return float(r@cho_solve(factor, r))

        def report(msg):
            print(f"{case} K={k}: {msg}", flush=True)
        null = score(np.zeros(4))
        if k == kmax and ref is not None:
            np.testing.assert_allclose([null, score(x_ref)],
                                       [float(ref[n]) for n in ("chi2_null","chi2_bf")],
                                       atol=1e-5, rtol=0.)
        report("global refit")
        seeds = [np.zeros(4), x_ref]
        if full_points is not None:
            seeds.append(xfull)
        (q, x, _, _), _, _ = fit_nested(score, "full", seeds, 3, args.seed, report)
        report("fixed-mass profile")
        curve = mass_profile(score, x, x_ref, masses, full_points, report)
        # A grid fit may find a better basin; feed that point back into the global search.
        grid_best = min(curve, key=lambda item: item[0])
        if grid_best[0] < q-1e-6:
            (q, x, _, _), _, _ = fit_nested(score, "full", [x, grid_best[1]], 3, args.seed, report)
        if q > grid_best[0]+1e-6:
            raise RuntimeError(f"{case} K={k}: global fit missed a lower grid point.")
        if k == kmax:
            full_points, xfull = curve, x
            full_curve = np.array([v for v, _ in curve])
            qfull, nullfull = q, null
        qcurve = np.array([v for v, _ in curve])
        bound = rows[k]["profiled_chi2_error_bound"]
        excess = np.r_[q-qfull, qcurve-full_curve]
        # Violations signal inconsistent optimization, not a failure of
        # the covariance bound. Preserve the profile output for inspection.
        worst_violation = max(0., -excess.min(), excess.max()-bound)
        delta_curve_diff = np.max(np.abs((null-qcurve)-(nullfull-full_curve)))
        profile_diff = np.max(np.abs((qcurve-q)-(full_curve-qfull)))
        fits.append(dict(case=case, k=k, chi2_null=null, chi2_bf=q, delta_chi2=null-q,
                         null_chi2_shift=null-nullfull, best_chi2_shift=q-qfull,
                         delta_chi2_shift=(null-q)-(nullfull-qfull),
                         max_null_relative_grid_difference=delta_curve_diff,
                         max_best_relative_grid_difference=profile_diff,
                         max_optimizer_bound_violation=worst_violation,
                         **parameter_row(x)))
        curves += [dict(case=case,k=k,chi2=v,delta_vs_null=null-v,delta_from_best=v-q,
                        **parameter_row(theta)) for v,theta in curve]
        if first_ratio_bin is not None:
            unit = np.eye(len(d))[first_ratio_bin]
            precision_row = cho_solve(factor, unit)
            conditional_sigma = 1 / np.sqrt(precision_row[first_ratio_bin])
        for name, theta in (("null", np.zeros(4)), ("best", x)):
            mc = prediction(theta)
            total, residual, pull_norm_sq, shift = profile(C, F[:k], 1., d-mc)
            pulls.append(pull_row(k, name, total, residual, pull_norm_sq))
            spectra += [dict(case=case,k=k,point=name,bin=int(keep[i])+1,data=d[i],mc=mc[i],
                             flux_shift=shift[i],profiled_mc=mc[i]+shift[i]) for i in range(len(d))]
            if first_ratio_bin is not None:
                i = first_ratio_bin
                pull = (precision_row @ (d-mc)) * conditional_sigma
                ratio_rows.append(dict(k=k, point=name, root_bin=int(keep[i])+1, data=d[i], raw_mc=mc[i],
                                       profiled_mc=mc[i]+shift[i],
                                       profiled_residual=d[i]-mc[i]-shift[i],
                                       conditional_prediction=d[i]-pull*conditional_sigma,
                                       conditional_sigma=conditional_sigma, conditional_pull=pull))
        report(f"delta chi2={null-q:.5f}; largest observed profile difference={profile_diff:.5f}")
        if worst_violation > 1e-5:
            report(f"WARNING: optimizer solutions violate the covariance bound by {worst_violation:.5g}.")
        for name, entries in (("refits",fits),("mass_profiles",curves),("profiled_spectra",spectra),
                              ("pull_vs_residual",pulls)):
            write_csv(args.out/f"{case}_{name}.csv",entries)
        if ratio_rows:
            write_csv(args.out/f"{case}_first_fhc_ratio_bin.csv", ratio_rows)
    if args.plots:
        plot_study(case, fits, curves, ratio_rows, args.out)
        plot_pulls(case, pulls, args.out)
    return rows, points


def pull_row(k, point, total, residual, pull_norm_sq):
    return dict(k=k, point=point, chi2=total, residual_chi2=residual,
                flux_pull_norm_sq=pull_norm_sq)


def fixed_k_objective(A, C, profile_bins, k):
    """Solve pulls on profile_bins and evaluate the penalized score on all supplied bins."""
    from scipy.linalg import solve_triangular

    B = mean_center(A) / np.sqrt(len(A))
    L = np.linalg.cholesky(C)
    L_profile = np.linalg.cholesky(C[np.ix_(profile_bins, profile_bins)])
    white = solve_triangular(L_profile, B[:, profile_bins].T, lower=True)
    U, s, Vh = np.linalg.svd(white, full_matrices=False)
    if not 0 <= k <= len(s):
        raise ValueError(f"K={k} exceeds the profiling basis of {len(s)} modes")
    # Keep the same right singular vectors when mapping the shifts to ratio bins.
    modes = B.T @ Vh[:k].T
    pull_from_white = (s[:k] / (1. + s[:k]**2))[:, None] * U[:, :k].T

    def evaluate(r):
        y = solve_triangular(L_profile, r[profile_bins], lower=True, check_finite=False)
        pulls = pull_from_white @ y
        shift = modes @ pulls
        after = solve_triangular(L, r-shift, lower=True, check_finite=False)
        residual, penalty = float(after @ after), float(pulls @ pulls)
        return residual+penalty, residual, penalty, shift

    return evaluate


def run_fixed_k(case, A, C, d, m, h, first_ratio_bin, keep, profile_keep, args, k, out):
    profile_bins = np.flatnonzero(np.isin(keep, profile_keep))
    evaluate = fixed_k_objective(A, C, profile_bins, k)
    def predict(x):
        h.OscillateHistogram(100*x[0], *x[1:], False, False)
        return vec(h.GetOscillatedHistogram())[keep]

    def score(x):
        return evaluate(d-predict(x))[0]

    def report(message):
        print(f"{case} K={k}: {message}", flush=True)

    zero = np.zeros(4)
    np.testing.assert_allclose(predict(zero), m, rtol=1e-7, atol=1e-8)
    null = evaluate(d-m)[0]
    if args.optimizer == "production":
        from tools.Fitters import OscillationFitter

        class FixedModeFitter(OscillationFitter):
            def CalChi2(self, x):
                return score(x)

        report("running production OscillationFitter.DoFit")
        fitter = FixedModeFitter(h, exclude="", marginalize_flux=False, hist_config=h.hist_config)
        best, params = fitter.DoFit()
        x = np.array([params["m"]/100., params["ue4"], params["umu4"], params["utau4"]])
        np.testing.assert_allclose(best, score(x), rtol=1e-9, atol=1e-8)
    else:
        seed = np.array([.1927, .014, .039, 0.] if case == "direct_all" else [.1927, .026, .044, 0.])
        (best, x, _, _), _, _ = fit_nested(score, "full", [zero, seed], 3, args.seed, report)
        masses = np.unique(np.r_[np.linspace(0, 100, args.mass_steps), 100*seed[0], 100*x[0], 96.21])
        curve = mass_profile(score, x, seed, masses, report=report)
        grid_best = min(curve, key=lambda item: item[0])
        if grid_best[0] < best-1e-6:
            (best, x, _, _), _, _ = fit_nested(score, "full", [x, grid_best[1]], 3, args.seed, report)
        if best > grid_best[0]+1e-6:
            raise RuntimeError(f"{case} K={k}: global fit missed a lower grid point")
        write_csv(out / f"{case}_mass_profiles.csv",
                  [dict(case=case, k=k, chi2=q, delta_from_best=q-best, **parameter_row(point))
                   for q, point in curve])

    result = dict(case=case, k=k, optimizer=args.optimizer,
                  n_profile_bins=len(profile_bins), n_score_bins=len(d),
                  chi2_null=null, chi2_bf=best, delta_chi2=null-best, **parameter_row(x))
    spectra, pulls = [], []
    profile_set = set(profile_bins)
    for name, point in (("null", zero), ("best", x)):
        mc = predict(point)
        q, residual, penalty, shift = evaluate(d-mc)
        pulls.append(pull_row(k, name, q, residual, penalty))
        spectra.extend(dict(case=case, k=k, point=name, bin=int(root_bin)+1,
                            used_for_profile=i in profile_set, data=d[i], mc=mc[i],
                            flux_shift=shift[i], profiled_mc=mc[i]+shift[i],
                            is_first_fhc_ratio_bin=i == first_ratio_bin)
                       for i, root_bin in enumerate(keep))
    write_csv(out / f"{case}_profiled_spectra.csv", spectra)
    write_csv(out / f"{case}_pull_vs_residual.csv", pulls)
    report(f"null={null:.6f}, best={best:.6f}, sterile gain={null-best:.6f}")
    if args.plots:
        from types import SimpleNamespace
        from fitData_normalizedFluxPrior import plot_fit_spectra

        # Reuse the exact fixed-K objective, including the 30-bin solve / 54-bin
        # score of ratios_dropped. The plotting code must not solve another prior.
        plot_prior = SimpleNamespace(evaluate=evaluate, L=np.linalg.cholesky(C), K=k)
        params = dict(zip(("m", "ue4", "umu4", "utau4"), (100*x[0], *x[1:])))
        tag = TAGS["ratio_all" if case in MASKED_CASES else case] + args.tag_suffix
        path = plot_fit_spectra(h, plot_prior, d, m, predict(x), params,
                               h.hist_config, out / "plots", f"{tag}_{case}",
                               args.ntuple_tag, "residual", 1.)
        report(f"post-fit plot: {path}")
    return result


def null_pull_curve(case, A, C, d, m, args):
    L, U, eigen, _ = mode_basis(A, C)
    z = U.T @ np.linalg.solve(L, d-m)
    # A retained mode leaves z/(1+eigen) in the residual and takes a pull
    # of sqrt(eigen)*z/(1+eigen). Dropped modes have zero pull.
    residual_by_mode = (z / (1 + eigen))**2
    pulls_by_mode = eigen * residual_by_mode
    residual = np.r_[0., np.cumsum(residual_by_mode)] + np.r_[
        np.cumsum((z*z)[::-1])[::-1], 0.]
    pulls = np.r_[0., np.cumsum(pulls_by_mode)]
    total = residual + pulls
    np.testing.assert_allclose(total, scores_by_k(d-m, L, U, eigen),
                               rtol=1e-9, atol=1e-8)
    rows = [pull_row(k, "null", q, r, p)
            for k, (q, r, p) in enumerate(zip(total, residual, pulls))]
    write_csv(args.out / f"{case}_pull_vs_residual.csv", rows)
    if args.plots:
        plot_pulls(case, rows, args.out)
    print(f"{case}: scanned K=0 to {len(eigen)}; null chi2 {total[0]:.6f} to {total[-1]:.6f}",
          flush=True)
    return rows


def plot_combined_null_curves(curves, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = {
        "direct_all": "Direct CC electron",
        "ratio_all": "Ratios and other samples",
        "ratios_dropped": "Ratios dropped",
        "only_ratios": "Only ratios",
    }
    fig, ax = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
    for case, rows in curves.items():
        positive = [row for row in rows
                    if row["residual_chi2"] > 0 and row["flux_pull_norm_sq"] > 0]
        if not positive:
            print(f"{case}: no positive points for the log–log plot", flush=True)
            continue
        x = [row["residual_chi2"] for row in positive]
        y = [row["flux_pull_norm_sq"] for row in positive]
        line, = ax.loglog(x, y, ".-", label=labels[case], markersize=4)
        for index in sorted({0, len(positive)-1, *(i for i, row in enumerate(positive)
                                                   if row["k"] in (10, 20, 30, 40))}):
            ax.annotate(str(positive[index]["k"]), (x[index], y[index]),
                        xytext=(4, 4), textcoords="offset points", fontsize=8,
                        color=line.get_color())
    ax.set(title="Null flux profiling by number of retained modes",
           xlabel=r"Residual $\chi^2$", ylabel="Squared flux pull norm")
    if ax.lines:
        ax.legend()
    ax.grid(True, which="both", alpha=.2)
    fig.savefig(out / "null_lcurves_combined_loglog.png", dpi=180)
    plt.close(fig)


def replay_pulls(case, A, C, d, keep, args):
    """Recover the two chi2 terms from a completed scan without refitting."""
    with (args.pulls_from / f"{case}_refits.csv").open(newline="") as f:
        fits = list(csv.DictReader(f))
    with (args.pulls_from / f"{case}_profiled_spectra.csv").open(newline="") as f:
        spectra = list(csv.DictReader(f))
    _, _, _, F = mode_basis(A, C)
    if not any(int(row["k"]) == len(F) for row in fits):
        raise ValueError(f"{case}: saved scan is missing the full-basis fit")
    pulls = []
    for fit in fits:
        k = int(fit["k"])
        if not 0 <= k <= len(F):
            raise ValueError(f"{case}: invalid saved K={k}")
        for point, field in (("null", "chi2_null"), ("best", "chi2_bf")):
            rows = sorted((row for row in spectra if int(row["k"]) == k and row["point"] == point),
                          key=lambda row: int(row["bin"]))
            np.testing.assert_array_equal([int(row["bin"]) for row in rows], keep+1)
            np.testing.assert_allclose([float(row["data"]) for row in rows], d, rtol=1e-10, atol=1e-10)
            r = d - np.array([float(row["mc"]) for row in rows])
            total, residual, pull_norm_sq, shift = profile(C, F[:k], 1., r)
            np.testing.assert_allclose(total, float(fit[field]), rtol=1e-9, atol=1e-7)
            np.testing.assert_allclose(shift, [float(row["flux_shift"]) for row in rows],
                                       rtol=1e-8, atol=1e-8)
            pulls.append(pull_row(k, point, total, residual, pull_norm_sq))
    write_csv(args.out / f"{case}_pull_vs_residual.csv", pulls)
    if args.plots:
        plot_pulls(case, pulls, args.out)


def plot_pulls(case, pulls, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    points = list(dict.fromkeys(row["point"] for row in pulls))
    fig, axes = plt.subplots(1, len(points), figsize=(5.5*len(points), 4),
                             constrained_layout=True, squeeze=False)
    for ax, point in zip(axes[0], points):
        rows = sorted((row for row in pulls if row["point"] == point), key=lambda row: row["k"])
        ax.plot([row["residual_chi2"] for row in rows], [row["flux_pull_norm_sq"] for row in rows], ".-")
        labels = rows if len(rows) <= 15 else rows[:-1:5] + rows[-1:]
        for row in labels:
            ax.annotate(str(row["k"]), (row["residual_chi2"], row["flux_pull_norm_sq"]),
                        xytext=(3, 4), textcoords="offset points")
        ax.set(title=point.capitalize(), xlabel=r"Residual $\chi^2$", ylabel="Squared flux pull norm")
    fig.savefig(out / f"{case}_pull_vs_residual.png", dpi=160)
    plt.close(fig)


def plot_study(case, fits, curves, ratio_rows, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fits = sorted(fits, key=lambda row: row["k"])
    ks = [row["k"] for row in fits]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    ax = axes[0, 0]
    for field, label in (("null_chi2_shift", "Null"), ("best_chi2_shift", "Best"),
                         ("delta_chi2_shift", "Sterile gain")):
        ax.plot(ks, [row[field] for row in fits], ".-", label=label)
    ax.set(xlabel="Modes retained", ylabel=r"Change in $\chi^2$ from full basis")
    ax.legend()
    ax = axes[0, 1]
    full = [row for row in curves if row["k"] == ks[-1]]
    for k in ks[:-1]:
        subset = [row for row in curves if row["k"] == k]
        ax.plot([row["dm2"] for row in subset],
                [a["delta_from_best"] - b["delta_from_best"] for a, b in zip(subset, full)],
                label=f"K={k}")
    ax.set(xlabel=r"$\Delta m^2$ [eV$^2$]", ylabel=r"Mass profile difference from full basis")
    if len(ks) > 1:
        ax.legend()
    axes[1, 0].plot(ks, [row["dm2"] for row in fits], ".-")
    axes[1, 0].set(xlabel="Modes retained", ylabel=r"Best-fit $\Delta m^2$ [eV$^2$]")
    ax = axes[1, 1]
    if ratio_rows:
        for point in ("null", "best"):
            subset = sorted((row for row in ratio_rows if row["point"] == point), key=lambda row: row["k"])
            ax.plot(ks, [row["profiled_residual"] for row in subset], ".-", label=point.capitalize())
        ax.set(xlabel="Modes retained", ylabel="First FHC ratio bin, data minus profiled MC")
    else:
        for field, label in (("ue4_sq", r"$|U_{e4}|^2$"), ("umu4_sq", r"$|U_{\mu4}|^2$"),
                             ("utau4_sq", r"$|U_{\tau4}|^2$")):
            ax.plot(ks, [row[field] for row in fits], ".-", label=label)
        ax.set(xlabel="Modes retained", ylabel="Best-fit mixing")
    ax.legend()
    fig.savefig(out / f"{case}_convergence.png", dpi=160)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--root", action="store_true", help="Legacy MC-only study")
    mode.add_argument("--current-root", action="store_true",
                      help="Current residual Flux covariance and shared-ratio statistics")
    mode.add_argument("--arrays", type=Path, help="Saved profile_arrays.npz")
    p.add_argument("--spectra", type=Path, help="Saved normalized spectra.csv for array mode")
    p.add_argument("--reference", type=Path, help="Same-tune normalized data_fits.csv for ROOT mode")
    p.add_argument("--out", type=Path, default=Path("flux_mode_convergence_v4p3p1"))
    p.add_argument("--cases", nargs="+", choices=(*TAGS, *MASKED_CASES), default=list(TAGS))
    p.add_argument("--ks", nargs="+", type=int, default=[20, 25, 30, 33, 40, 54])
    p.add_argument("--tolerance", type=float, default=.1, help="Absolute chi2 tolerance for certification")
    p.add_argument("--mass-steps", type=int, default=101)
    p.add_argument("--seed", type=int, default=431)
    p.add_argument("--plots", action="store_true",
                   help="Save convergence plots, or ROOT post-fit spectra for --sterile-k")
    p.add_argument("--ntuple-tag", default="MAD", help="Label/prefix for post-fit spectrum plots")
    p.add_argument("--pulls-from", type=Path,
                   help="Recover pull curves from a completed current ROOT scan without refitting")
    p.add_argument("--null-only", action="store_true",
                   help="Scan every K at zero oscillations, skipping all sterile refits")
    p.add_argument("--tag-suffix", default="",
                   help="Append to both stitched input tags, e.g. _updatedStat_newEL")
    p.add_argument("--sterile-k", type=int, nargs="+",
                   help="Fit each requested K. For ratios_dropped, exclude ratios only from the pull solve.")
    p.add_argument("--include-full-basis", action="store_true",
                   help="Also fit all available modes: 54 for all-bin cases, 30 for ratios_dropped")
    optimizer = p.add_mutually_exclusive_group()
    optimizer.add_argument("--optimizer", choices=("mode-study", "production"), default="mode-study",
                           help="Sterile fit search: mode-study includes the mass scan; production uses DoFit.")
    optimizer.add_argument("--production-optimizer", dest="optimizer", action="store_const", const="production",
                           help="Alias for --optimizer production.")
    args = p.parse_args()
    if not np.isfinite(args.tolerance) or args.tolerance <= 0 or args.mass_steps < 3 or min(args.ks) < 0 or args.seed < 0:
        p.error("Use positive finite tolerance, at least three mass steps, and nonnegative K and seed.")
    if args.root and args.reference is None:
        p.error("--root requires --reference.")
    if args.arrays and args.spectra is None:
        p.error("--arrays requires --spectra.")
    if not args.current_root and (args.null_only or args.pulls_from or args.tag_suffix or
                                  args.sterile_k is not None or
                                  set(args.cases).intersection(MASKED_CASES)):
        p.error("Masked cases and the current ROOT scan options require --current-root.")
    if args.null_only and args.pulls_from:
        p.error("Choose either --null-only or --pulls-from.")
    if args.sterile_k is not None and (min(args.sterile_k) < 0 or args.null_only or args.pulls_from):
        p.error("--sterile-k requires a nonnegative K and a fresh sterile fit.")
    if args.optimizer == "production" and args.sterile_k is None:
        p.error("--optimizer production requires --sterile-k.")
    if args.include_full_basis and args.sterile_k is None:
        p.error("--include-full-basis requires --sterile-k.")
    args.out.mkdir(parents=True, exist_ok=True)
    bounds, fixed = [], []
    if args.current_root:
        null_curves = {}
        sterile_fits = []
        for case, A, C, d, m, h, first_ratio_bin, keep, profile_keep in current_root_inputs(args):
            if args.sterile_k is not None:
                ks = list(dict.fromkeys(args.sterile_k +
                           ([len(profile_keep)] if args.include_full_basis else [])))
                for k in ks:
                    out = args.out if len(ks) == 1 else args.out / f"k{k}"
                    out.mkdir(parents=True, exist_ok=True)
                    sterile_fits.append(run_fixed_k(case, A, C, d, m, h, first_ratio_bin, keep,
                                                    profile_keep, args, k, out))
                    write_csv(args.out/"sterile_refits.csv", sterile_fits)
                continue
            if args.null_only:
                null_curves[case] = null_pull_curve(case, A, C, d, m, args)
                continue
            if args.pulls_from:
                replay_pulls(case, A, C, d, keep, args)
                continue
            a, b = root_study(case, A, C, d, m, {}, h, {}, args, first_ratio_bin, keep)
            bounds.extend(a)
            fixed.extend(b)
            write_csv(args.out/"mode_bounds.csv", bounds)
            write_csv(args.out/"fixed_points.csv", fixed)
        if args.null_only and args.plots:
            plot_combined_null_curves(null_curves, args.out)
    elif args.root:
        with args.reference.open(newline="") as f:
            refs = {r["case"]: r for r in csv.DictReader(f) if r["mode"] == "normalized"
                    and r["basis"] == "mean_centered" and int(r["nprof"]) == 1000
                    and float(r["lam"]) == 1000}
        if not set(args.cases) <= refs.keys():
            p.error("Missing full normalized P8 reference rows.")
        for case,A,C,d,m,templates,checks,(h,_,_) in root_inputs(args):
            if len(A) != 1000 or not any(r["comparison"] == "MC Flux band" and
                                        r["basis"] == "mean_centered" and r["passed"] for r in checks):
                raise RuntimeError(f"{case}: normalized P8 prior did not close to ROOT.")
            a,b = root_study(case,A,C,d,m,templates,h,refs,args)
            bounds.extend(a)
            fixed.extend(b)
            write_csv(args.out/"mode_bounds.csv",bounds)
            write_csv(args.out/"fixed_points.csv",fixed)
    else:
        with args.spectra.open(newline="") as f:
            spectra = list(csv.DictReader(f))
        with np.load(args.arrays) as data:
            for case in args.cases:
                A,C,d,m = [data[TAGS[case]+"__"+k] for k in ("A","C_nonflux","data","mc")]
                rows = sorted([r for r in spectra if r["case"] == case and r["point"] == "bf"
                               and r["basis"] == "mean_centered" and int(r["nprof"]) == len(A)
                               and float(r["lam"]) == len(A)],key=lambda r:int(r["bin"]))
                if [int(r["bin"]) for r in rows] != list(range(1,len(d)+1)):
                    p.error(f"{case}: expected exactly one full-P8 best-fit prediction per bin.")
                np.testing.assert_allclose([float(r["data"]) for r in rows],d)
                a,b,_,_ = fixed_checks(case,A,C,d,m,{"null":m,"saved_best":np.array(
                    [float(r["raw_mc"]) for r in rows])},args.tolerance)
                bounds.extend(a)
                fixed.extend(b)
        write_csv(args.out/"mode_bounds.csv",bounds)
        write_csv(args.out/"fixed_points.csv",fixed)
    print(f"Results: {args.out.resolve()}")


if __name__ == "__main__":
    main()

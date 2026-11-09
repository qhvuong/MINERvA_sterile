#!/usr/bin/env python3
"""v4.3.1 full-P8 mode convergence. Does not change the standard fitter."""
import argparse
import csv
from pathlib import Path

import numpy as np

from compare_normalized_flux_profile import TAGS, mean_center, profile, root_inputs, vec, write_csv
from test_electron_holdouts import fit_nested
from test_electron_mixing import MAX, PARAMS, parameter_row


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


def root_study(case, A, C, d, m, templates, h, refs, args):
    from scipy.linalg import cho_factor, cho_solve

    ref = refs[case]
    x_ref = np.array([float(ref[name]) for name in PARAMS])
    x_ref[0] /= 100.
    def prediction(x):
        h.OscillateHistogram(100*x[0], *x[1:], False, False)
        return vec(h.GetOscillatedHistogram())

    np.testing.assert_allclose(prediction(np.zeros(4)), m, rtol=1e-7, atol=1e-8)
    rows, points, F, certified = fixed_checks(
        case, A, C, d, m, {"null": m, "saved_best": prediction(x_ref), **templates}, args.tolerance)
    kmax = len(F)
    ks = [kmax] + sorted((set(args.ks) | {max(0, certified-1), certified}) - {kmax})
    if max(ks) > kmax:
        raise ValueError(f"{case}: cannot retain more than {kmax} modes.")
    masses = np.unique(np.r_[np.linspace(0, 100, args.mass_steps), 100*x_ref[0], 19.26, 96.21])
    fits, curves, spectra = [], [], []
    full_points = None
    for k in ks:
        factor = cho_factor(C+F[:k].T@F[:k])
        def score(x):
            r = d-prediction(x)
            return float(r@cho_solve(factor, r))

        def report(msg):
            print(f"{case} K={k}: {msg}", flush=True)
        null = score(np.zeros(4))
        if k == kmax:
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
                         delta_chi2_shift=(null-q)-(nullfull-qfull),
                         max_null_relative_grid_difference=delta_curve_diff,
                         max_best_relative_grid_difference=profile_diff,
                         max_optimizer_bound_violation=worst_violation,
                         **parameter_row(x)))
        curves += [dict(case=case,k=k,chi2=v,delta_vs_null=null-v,delta_from_best=v-q,
                        **parameter_row(theta)) for v,theta in curve]
        for name, theta in (("null", np.zeros(4)), ("best", x)):
            mc = prediction(theta)
            _, _, _, shift = profile(C, F[:k], 1., d-mc)
            spectra += [dict(case=case,k=k,point=name,bin=i+1,data=d[i],mc=mc[i],
                             flux_shift=shift[i],profiled_mc=mc[i]+shift[i]) for i in range(len(d))]
        report(f"delta chi2={null-q:.5f}; largest observed profile difference={profile_diff:.5f}")
        if worst_violation > 1e-5:
            report(f"WARNING: optimizer solutions violate the covariance bound by {worst_violation:.5g}.")
        for name, entries in (("refits",fits),("mass_profiles",curves),("profiled_spectra",spectra)):
            write_csv(args.out/f"{case}_{name}.csv",entries)
    return rows, points


def main():
    p = argparse.ArgumentParser(description=__doc__)
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--root", action="store_true")
    mode.add_argument("--arrays", type=Path, help="Saved profile_arrays.npz")
    p.add_argument("--spectra", type=Path, help="Saved normalized spectra.csv for array mode")
    p.add_argument("--reference", type=Path, help="Same-tune normalized data_fits.csv for ROOT mode")
    p.add_argument("--out", type=Path, default=Path("flux_mode_convergence_v4p3p1"))
    p.add_argument("--cases", nargs="+", choices=TAGS, default=list(TAGS))
    p.add_argument("--ks", nargs="+", type=int, default=[10,20,30,33,40])
    p.add_argument("--tolerance", type=float, default=.1, help="Absolute chi2 tolerance for certification")
    p.add_argument("--mass-steps", type=int, default=101)
    p.add_argument("--seed", type=int, default=431)
    args = p.parse_args()
    if not np.isfinite(args.tolerance) or args.tolerance <= 0 or args.mass_steps < 3 or min(args.ks) < 0 or args.seed < 0:
        p.error("Use positive finite tolerance, at least three mass steps, and nonnegative K and seed.")
    if args.root and args.reference is None:
        p.error("--root requires --reference.")
    if args.arrays and args.spectra is None:
        p.error("--arrays requires --spectra.")
    args.out.mkdir(parents=True, exist_ok=True)
    bounds, fixed = [], []
    if args.root:
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

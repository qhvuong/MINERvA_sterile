#!/usr/bin/env python3
"""Fit pure muon disappearance with both CC electron-sensitive blocks removed."""
import argparse
import csv
import json
from functools import lru_cache
from pathlib import Path

import numpy as np

from compare_normalized_flux_profile import (
    BASE_TUNE, TUNES, histogram_tag, mean_center, profile, root_inputs, vec, write_csv,
)
from test_electron_holdouts import make_holdout
from test_electron_mixing import MAX

CORE = ("fhc_elastic", "fhc_numu_selection", "rhc_numu_selection")
U_MAX = MAX[2]
BASE_MASS = 19.262691987596064


def amplitude(u):
    return 4*u*(1-u)


def fit_mass_grid(predict, null_mc, data, V, keep, grid, reference_mass, report, normalize=None):
    # SciPy must load after PyROOT on GPVM.
    from scipy.linalg import cho_factor, cho_solve
    from scipy.optimize import minimize_scalar

    factor = cho_factor(V[np.ix_(keep, keep)])
    null_pred = null_mc if normalize is None else normalize(null_mc)[0]
    residual = (data-null_pred)[keep]
    r_weighted = cho_solve(factor, residual)
    null_q = float(residual @ r_weighted)

    @lru_cache(maxsize=None)
    def at_mass(mass):
        if normalize is not None:
            # Cache the ROOT response once per mass; normalization is
            # profiled in NumPy at each trial oscillation amplitude.
            response = (predict(mass, .15)-null_mc) / amplitude(.15)

            def score(a):
                mc, _ = normalize(null_mc+a*response)
                r = (data-mc)[keep]
                return float(r @ cho_solve(factor, r))

            grid_a = np.linspace(0., amplitude(U_MAX), 17)
            scores = [score(a) for a in grid_a]
            options = list(zip(scores, grid_a))
            intervals = [(grid_a[0], grid_a[1]), (grid_a[-2], grid_a[-1])]
            for i in range(1, len(grid_a)-1):
                if scores[i] <= min(scores[i-1], scores[i+1]) and scores[i] < max(scores[i-1], scores[i+1]):
                    intervals.append((grid_a[i-1], grid_a[i+1]))
            for bounds in intervals:
                opt = minimize_scalar(score, method="bounded", bounds=bounds,
                                      options={"xatol": 1e-9})
                options.append((float(opt.fun), float(opt.x)))
            q, a = min(options)
            return q, a / (2*(1 + np.sqrt(1-a)))

        # For this 30-bin core with Ue^2=Utau^2=0 the prediction is linear
        # in the survival amplitude at fixed mass. Fit the amplitude exactly.
        probe_u = .15
        response = (predict(mass, probe_u)-null_mc)[keep] / amplitude(probe_u)
        denom = float(response @ cho_solve(factor, response))
        a = float(np.clip((response @ r_weighted)/denom, 0., amplitude(U_MAX))) if denom > 0 else 0.
        u = a / (2*(1 + np.sqrt(1-a)))
        r = residual-a*response
        q = float(r @ cho_solve(factor, r))
        return q, u

    points = []
    for i, mass in enumerate(grid, 1):
        q, u = at_mass(float(mass))
        points.append((q, float(mass), u))
        if i == 1 or i % 50 == 0 or i == len(grid):
            report(f"mass grid {i}/{len(grid)}")
    candidates = points.copy()
    for i in range(1, len(grid)-1):
        q = points[i][0]
        if q <= min(points[i-1][0], points[i+1][0]) and q < max(points[i-1][0], points[i+1][0]):
            opt = minimize_scalar(lambda mass: at_mass(mass)[0], method="bounded",
                                  bounds=(grid[i-1], grid[i+1]), options={"xatol": 1e-7})
            mass = float(opt.x)
            q, u = at_mass(mass)
            candidates.append((q, mass, u))
    q, u = at_mass(reference_mass)
    anchored = (q, reference_mass, u)
    best = min([*candidates, anchored, (null_q, 0., 0.)])
    return (null_q, 0., 0.), best, anchored, points


def rate_profiler(data, V, keep, blocks):
    """Profile one nonnegative multiplicative rate per retained muon block."""
    from scipy.linalg import solve_triangular
    from scipy.optimize import nnls

    lower = np.linalg.cholesky(V[np.ix_(keep, keep)])
    block_masks = {name: np.isin(keep, indices) for name, indices in blocks.items()}
    block_masks = {name: mask for name, mask in block_masks.items() if mask.any()}
    elastic = ~np.logical_or.reduce(list(block_masks.values()))
    dw = solve_triangular(lower, data[keep], lower=True)

    def normalize(mc):
        retained = mc[keep]
        X = np.column_stack([retained*mask for mask in block_masks.values()])
        Xw = solve_triangular(lower, X, lower=True)
        yw = dw - solve_triangular(lower, retained*elastic, lower=True)
        # Unconstrained GLS is exact whenever all scales are positive.
        scales = np.linalg.lstsq(Xw, yw, rcond=None)[0]
        if np.any(scales < 0.):
            scales = nnls(Xw, yw)[0]
        pred = mc.copy()
        for scale, mask in zip(scales, block_masks.values()):
            pred[keep[mask]] *= scale
        return pred, dict(zip(block_masks, scales))

    return normalize


def heldout_rate(data_held, mc_held, cond_shift, V, keep, held):
    """Give the omitted block its own scale for a conditional shape score."""
    from scipy.linalg import cho_factor, cho_solve

    cross = V[np.ix_(held, keep)]
    kept = cho_factor(V[np.ix_(keep, keep)])
    cond = cho_factor(V[np.ix_(held, held)] - cross @ cho_solve(kept, cross.T))
    direction = cho_solve(cond, mc_held)
    return max(0., float(direction @ (data_held-cond_shift) / (direction @ mc_held)))


def check_same_data(case, cfg, groups, data, h, core):
    """Check data counts and sample bin edges against the v4.3.1 file."""
    import ROOT

    tag = histogram_tag(case)
    osc = cfg.parent
    base_cfg = json.loads((osc / f"HIST_CONFIG_{tag}.json").read_text())
    for name in CORE:
        for key in ("start", "end"):
            if base_cfg[name][key] != groups[name][key]:
                raise RuntimeError(f"{name}: stitched bin indices differ from {BASE_TUNE}.")
    path = osc / "rootfiles" / f"NuE_stitched_hists_{tag}.root"
    f = ROOT.TFile.Open(str(path), "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open baseline data file {path}.")
    def edges(hist):
        axis = hist.GetXaxis()
        return [axis.GetBinLowEdge(i) for i in range(1, hist.GetNbinsX()+2)]

    try:
        hist = f.Get("sample_data")
        if not hist:
            raise RuntimeError(f"No sample_data in {path}.")
        np.testing.assert_allclose(data[core], vec(hist)[core], rtol=0., atol=1e-8,
                                   err_msg=f"Core observed data differ from {BASE_TUNE}")
        for name in CORE:
            hist = f.Get("data_" + name)
            other = h.data_hists.get(name)
            if not hist or not other:
                raise RuntimeError(f"{name}: cannot verify the sample bin edges in both files.")
            np.testing.assert_allclose(edges(other), edges(hist), rtol=0., atol=1e-10,
                                       err_msg=f"{name} bin edges differ from {BASE_TUNE}")
    finally:
        f.Close()
    print(f"{case}: 30 observed bins and sample bin edges match {BASE_TUNE}", flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", type=Path, help="Optional same-tune electron holdout holdout_fits.csv")
    p.add_argument("--tune", choices=TUNES, default=BASE_TUNE)
    p.add_argument("--out", type=Path)
    p.add_argument("--cases", nargs="+", choices=("direct_all", "ratio_all"), default=["direct_all"])
    p.add_argument("--grid-steps", type=int, default=401)
    p.add_argument("--float-muon-normalizations", action="store_true",
                   help="Remove absolute muon-rate information with one free scale per beam mode.")
    args = p.parse_args()
    if args.grid_steps < 3:
        p.error("Use at least three mass grid points.")
    refs = {}
    if args.reference:
        with args.reference.open(newline="") as f:
            refs = {(r["case"], r["model"]): r for r in csv.DictReader(f) if r["holdout"] == "both"}
        for case in args.cases:
            for model in ("null", "no_e"):
                if (case, model) not in refs:
                    p.error(f"Missing both-omitted reference row for {case}/{model}.")
                row = refs[case, model]
                if abs(float(row["utau4_sq"])) > 1e-8 or abs(float(row["ue4_sq"])) > 1e-8:
                    p.error(f"{case}/{model}: expected pure muon-disappearance reference.")
    if args.out is None:
        kind = "muon_shape_holdouts" if args.float_muon_normalizations else "muon_holdouts"
        args.out = Path(f"{kind}_{args.tune}")
    args.out.mkdir(parents=True, exist_ok=True)
    comparisons, fits, profiles, spectra = [], [], [], []
    grid = np.linspace(0, 100, args.grid_steps)
    for case, A, C, d, m, _, checks, backend in root_inputs(args):
        if len(A) != 1000 or not any(r["comparison"] == "MC Flux band"
                                    and r["basis"] == "mean_centered" and r["passed"] for r in checks):
            raise RuntimeError(f"{case}: normalized flux prior did not close to ROOT.")
        h, cfg, _ = backend
        groups = json.loads(cfg.read_text())
        indices = {name: np.arange(groups[name]["start"], groups[name]["end"]+1) for name in CORE}
        core = np.concatenate([indices[name] for name in CORE])
        if len(set(core)) != 30 or [len(indices[n]) for n in CORE] != [6, 12, 12]:
            raise RuntimeError(f"{case}: unexpected 30-bin core definition.")
        if args.tune != BASE_TUNE:
            check_same_data(case, cfg, groups, d, h, core)
        labels = [(name, int(i)+1, j+1) for name in CORE for j, i in enumerate(indices[name])]
        mask = np.ix_(core, core)
        B = mean_center(A)[:, core]
        C, data = C[mask], d[core]
        flux = B.T @ B / len(B)
        V = C + flux
        mu_blocks = {"fhc": np.arange(6, 18), "rhc": np.arange(18, 30)}

        @lru_cache(maxsize=None)
        def predict(mass, u):
            h.OscillateHistogram(float(mass), 0., float(u), 0., False, False)
            return vec(h.GetOscillatedHistogram())[core]

        null_mc = predict(0., 0.)
        np.testing.assert_allclose(null_mc, m[core], rtol=1e-7, atol=1e-8)
        for mass in (1., 19.26, 80.):
            response = (predict(mass, .15)-null_mc) / amplitude(.15)
            np.testing.assert_allclose(predict(mass, .04), null_mc+amplitude(.04)*response,
                                       rtol=1e-8, atol=1e-8)
        reference_mass = BASE_MASS
        if refs:
            for model in ("null", "no_e"):
                row = refs[case, model]
                r = data - predict(float(row["dm2"]), float(row["umu4_sq"]))
                np.testing.assert_allclose(r @ np.linalg.solve(V, r), float(row["chi2_retained"]),
                                           atol=1e-5, rtol=0.)
            if args.tune == BASE_TUNE:
                reference_mass = float(refs[case, "no_e"]["dm2"])

        for key, held in (("none", np.array([], dtype=int)),
                          ("fhc", np.arange(6, 18)), ("rhc", np.arange(18, 30))):
            keep = np.setdiff1d(np.arange(len(core)), held)
            has_holdout = held.size > 0
            shape = args.float_muon_normalizations
            normalize = rate_profiler(data, V, keep, mu_blocks) if shape else None
            meta = dict(case=case, tune=args.tune, holdout=key, retained_bins=len(keep), held_bins=len(held),
                        normalization="free_per_muon_block" if shape else "nominal")
            def report(message):
                print(f"{case} drop_{key}: {message}", flush=True)

            if has_holdout:
                _, _, evaluate, cond_sigma = make_holdout(C, flux, held)
            null, best, anchored, points = fit_mass_grid(
                predict, null_mc, data, V, keep, grid, reference_mass, report, normalize)
            if refs and not has_holdout and not shape:
                np.testing.assert_allclose(best[0], float(refs[case, "no_e"]["chi2_retained"]),
                                           atol=1e-5, rtol=0.)
            held_scores = {}
            for model, (q, mass, u) in (("null", null), ("best", best), ("fixed_core_mass", anchored)):
                raw_mc = predict(mass, u)
                mc, rates = normalize(raw_mc) if shape else (raw_mc.copy(), {})
                r = data-mc
                c_mask = np.ix_(keep, keep)
                total, residual, penalty, shift_keep = profile(C[c_mask], B[:, keep], len(B), r[keep])
                np.testing.assert_allclose(q, total, rtol=1e-8, atol=1e-7)
                held_q = None
                held_scale = None
                if has_holdout:
                    train, held_q, shift, cond_shift = evaluate(r)
                    if shape:
                        held_scale = heldout_rate(data[held], mc[held], cond_shift, V, keep, held)
                        mc[held] *= held_scale
                        r = data-mc
                        train, held_q, shift, cond_shift = evaluate(r)
                    np.testing.assert_allclose(train, q, rtol=1e-8, atol=1e-7)
                    np.testing.assert_allclose(shift[keep], shift_keep, rtol=1e-8, atol=1e-7)
                    np.testing.assert_allclose(train+held_q, r @ np.linalg.solve(V, r),
                                               rtol=1e-8, atol=1e-7)
                    for j, i in enumerate(held):
                        sample, bin_global, bin_local = labels[i]
                        spectra.append(dict(**meta, model=model, sample=sample, bin=bin_global,
                                            sample_bin=bin_local, data=data[i], raw_mc=raw_mc[i],
                                            rate_adjusted_mc=mc[i],
                                            flux_only_prediction=mc[i]+shift[i],
                                            conditional_prediction=mc[i]+cond_shift[j],
                                            conditional_sigma=cond_sigma[j]))
                held_scores[model] = held_q
                fits.append(dict(**meta, model=model, dm2=mass, umu4_sq=u, chi2_retained=q,
                                 chi2_residual=residual, flux_penalty=penalty, chi2_held_conditional=held_q,
                                 fitted_fhc_scale=rates.get("fhc"), fitted_rhc_scale=rates.get("rhc"),
                                 held_shape_scale=held_scale))
            def transfer(model):
                return held_scores["null"]-held_scores[model] if has_holdout else None

            comparisons.append(dict(**meta, best_dm2=best[1], best_umu4_sq=best[2],
                                    delta_chi2_retained=null[0]-best[0],
                                    delta_chi2_held_conditional=transfer("best"),
                                    fixed_core_mass_dm2=reference_mass,
                                    fixed_core_mass_delta_retained=null[0]-anchored[0],
                                    fixed_core_mass_delta_held_conditional=transfer("fixed_core_mass")))
            profiles += [dict(**meta, dm2=mass, umu4_sq=u, chi2=q, delta_chi2=null[0]-q)
                         for q, mass, u in points]
            report(f"best dm2={best[1]:.4f}; retained gain={null[0]-best[0]:.4f}; "
                   f"held gain={transfer('best')}")
            for name, rows in (("muon_holdout_comparisons", comparisons), ("muon_holdout_fits", fits),
                               ("mass_profiles", profiles), ("heldout_spectra", spectra)):
                write_csv(args.out / f"{name}.csv", rows)
    print("Results:", args.out.resolve())


if __name__ == "__main__":
    main()

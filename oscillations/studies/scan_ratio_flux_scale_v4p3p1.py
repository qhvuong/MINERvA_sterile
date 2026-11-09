#!/usr/bin/env python3
"""Compare v4.3.1 ratio mass profiles at s=1 and 5, then hold out FHC ratio bin 1."""
import argparse
import csv
from functools import lru_cache
import os
from pathlib import Path
import sys

import numpy as np

TAG = "predNueel_p8_v4p3p1"
MAX = np.array([.394, .489, .718])
STARTS = ([0, 0, 0], [.025, .044, 0], [.05, .01, 0],
          [.02, .02, .5], [.22, .26, 0], [.38, .48, .10], [.15, .12, .70])
FIELDS = ("dm2", "ue4_sq", "umu4_sq", "utau4_sq")


def feasible(x, upper):
    return (np.all(np.isfinite(x)) and np.all(x >= -1e-8)
            and np.all(x <= upper + 1e-8) and x[-3:].sum() <= 1 + 1e-8)


def sanitize(x, upper):
    x = np.clip(x, 0, upper)
    x[-3:] /= max(1., x[-3:].sum())
    return x


def fit_local(objective, starts, upper):
    from scipy.optimize import LinearConstraint, minimize

    cons = LinearConstraint([np.r_[np.zeros(len(upper) - 3), np.ones(3)]], -np.inf, 1.)
    tries = []
    for start in {tuple(x) for x in starts}:
        x = sanitize(np.asarray(start, dtype=float), upper)
        tries.append((objective(x), x, "seed"))
        result = minimize(objective, x, method="SLSQP", bounds=list(zip(np.zeros(len(x)), upper)),
                          constraints=cons, options={"ftol": 1e-9, "maxiter": 140})
        if feasible(result.x, upper):
            y = sanitize(result.x, upper)
            tries.append((objective(y), y, "converged" if result.success else "not_converged"))
    finite = [row for row in tries if np.isfinite(row[0])]
    if not finite:
        raise RuntimeError("No finite feasible oscillation candidate")
    return min(finite, key=lambda row: (row[0], row[2] == "not_converged"))


def scan(label, score, masses, seeds=(), previous=None):
    rows = []
    warm = np.zeros(3)
    for i, mass in enumerate(masses):
        if previous is None:
            starts = [*STARTS, warm, *(x[1:] for x in seeds)]
        else:
            starts = [previous[i][1][1:], warm, *(x[1:] for x in seeds)]
        value, warm, state = fit_local(lambda u: score(np.r_[mass, u]), starts, MAX)
        if previous is not None and previous[i][0] <= value:
            value, point, state = previous[i]
            warm = point[1:]
        rows.append((value, np.r_[mass, warm], state))
        if (i + 1) % 40 == 0 or i + 1 == len(masses):
            print(f"{label}: fixed mass {i+1}/{len(masses)}", flush=True)
    return rows


def polish(score, rows, seeds):
    valleys = [rows[0], rows[-1]]
    valleys.extend(row for left, row, right in zip(rows, rows[1:], rows[2:])
                   if row[0] <= min(left[0], right[0]) and row[0] < max(left[0], right[0]))
    starts = [*(x for _, x, _ in sorted(valleys, key=lambda row: row[0])[:12]), *seeds]
    candidates = []
    upper = np.r_[1., MAX]
    for point in starts:
        x = np.r_[point[0]/100., point[1:]]
        value, y, state = fit_local(lambda z: score(np.r_[100*z[0], z[1:]]), [x], upper)
        candidates.append((value, np.r_[100*y[0], y[1:]], state))
    return sorted(candidates, key=lambda row: row[0])


def make_score(data, cov, predict, keep):
    from scipy.linalg import solve_triangular

    lower = np.linalg.cholesky(cov[np.ix_(keep, keep)])

    def score(x):
        r = data[keep] - predict(x)[keep]
        if not np.all(np.isfinite(r)):
            return 1e30
        y = solve_triangular(lower, r, lower=True, check_finite=False)
        return float(y @ y)

    return score


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("ratio_scale_scan_v4p3p1"))
    parser.add_argument("--mass-steps", type=int, default=401, help="Grid points from 0 to 100 eV²")
    parser.add_argument("--seed", action="append", type=float, nargs=4,
                        metavar=("DM2", "UE2", "UMU2", "UTAU2"),
                        help="Optional saved best fit; repeat to include both prior scales")
    args = parser.parse_args()
    if args.mass_steps < 3:
        parser.error("--mass-steps must be at least 3")
    seeds = [np.array([19.259947, .025, .04433, 0.])]
    for seed in args.seed or []:
        if not feasible(np.array(seed), np.r_[100., MAX]):
            parser.error(f"Unphysical seed {seed}")
        seeds.append(sanitize(np.array(seed), np.r_[100., MAX]))

    osc = Path(__file__).resolve().parents[1]
    cfg = osc / f"HIST_CONFIG_{TAG}.json"
    rootfile = osc / "rootfiles" / f"NuE_stitched_hists_{TAG}.root"
    for path in (cfg, rootfile):
        if not path.is_file():
            parser.error(f"Missing {path}")
    sys.path[:0] = [str(osc), str(osc.parent)]
    os.environ.setdefault("CCNUEROOT", str(osc.parent))
    sys.argv = [sys.argv[0]]
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    from fitData_normalizedFluxPrior import NormalizedFluxPrior, vector
    from tools.Helper import GetSliceIndices, sin_average
    from tools.StitchedHistogram import StitchedHistogram
    from scipy.linalg import cho_factor, cho_solve

    h = StitchedHistogram("sample", include_shared_ratio_stats=True)
    h.SetHistConfig(str(cfg))
    h.Load(str(rootfile))
    dh = h.GetDataHistogram()
    d, m = vector(dh), vector(h.GetMCHistogram())
    C = h.GetCovarianceMatrix(True)
    root_full = h.GetCovarianceMatrix(False)
    band = dh.GetVertErrorBand("Flux")
    A = np.asarray(h.GetAMatrix(), dtype=float)
    if A.shape != (1000, len(d)) or not band or band.GetNHists() != len(A):
        raise RuntimeError("Expected aligned 1000-universe Data and MC Flux bands")
    shifts = np.array([vector(band.GetHist(i)) - d for i in range(len(A))]) - A
    prior = NormalizedFluxPrior(C, shifts)
    prior.assert_flux_covariance(root_full - C)
    F = prior.flux_cov
    ratio_mask = vector(h.ratio_id)
    elastic_mask = vector(h.elastic_id)
    if not (np.isin(ratio_mask, [0., 1.]).all() and np.isin(elastic_mask, [0., 1.]).all()):
        raise RuntimeError("Expected binary stitched sample masks")
    ratio = ratio_mask.astype(bool)
    electron, muon, swapped = map(vector, (h.nue_hist, h.numu_hist, h.swap_hist))

    @lru_cache(maxsize=2048)
    def oscillation_response(mass):
        return tuple(4 * cv * np.array([sin_average(i, mass, template)
                                       for i in range(1, len(d)+1)])
                     for cv, template in ((electron, h.nue_template), (muon, h.numu_template),
                                          (swapped, h.swap_template)))

    def predict(x):
        mass, e, u, t = x
        de, du, ds = oscillation_response(float(mass))
        nue = electron - e*(1-e)*de + e*u*ds
        numu = muon - u*(1-u)*du
        pred = nue + numu
        pred[ratio] = np.divide(numu[ratio], nue[ratio], out=np.zeros(ratio.sum()),
                                where=nue[ratio] != 0.)
        tau = elastic_mask * (t*(u*du + e*de))
        return pred + tau if tau.sum() > 0 else pred

    def check_prediction(x):
        h.OscillateHistogram(*map(float, x), False, False)
        np.testing.assert_allclose(predict(x), vector(h.GetOscillatedHistogram()),
                                   rtol=1e-7, atol=1e-8)

    zero = np.zeros(4)
    np.testing.assert_allclose(predict(zero), m, rtol=1e-7, atol=1e-8)
    for point in (zero, *seeds, [2, .04, .11, .5], [39.2, .015, .11, .22], [98, .2, .25, .3]):
        check_prediction(point)
    print("ROOT covariance and oscillation prediction checks passed.", flush=True)

    ratio_bins = GetSliceIndices(str(cfg), "", ["fhc_ratio"])
    if not ratio_bins:
        raise RuntimeError("No FHC ratio bins in the histogram configuration")
    held = ratio_bins[0]
    if not ratio[held]:
        raise RuntimeError("The requested holdout is not a ratio bin")
    all_bins = np.arange(len(d))
    keep = all_bins[all_bins != held]
    masses = np.unique(np.r_[np.linspace(0., 100., args.mass_steps), [x[0] for x in seeds]])
    covs = {s: C + s*s*F for s in (1, 5)}
    scores = {s: make_score(d, covs[s], predict, all_bins) for s in covs}
    preliminary = {s: scan(f"s={s} initial", scores[s], masses, seeds) for s in scores}
    basins = {s: polish(scores[s], preliminary[s], seeds) for s in scores}
    cross_seeds = [*seeds, *(basins[s][0][1] for s in scores)]

    profiles, fits, candidates = {}, {}, {}
    for s in scores:
        label = f"all_s{s}"
        profiles[label] = scan(label, scores[s], masses, cross_seeds, preliminary[s])
        candidates[label] = polish(scores[s], profiles[label], cross_seeds)
        fits[label] = (scores[s](zero), *candidates[label][0])
    label = "holdout_s1"
    kept_score = make_score(d, covs[1], predict, keep)
    held_seeds = [*cross_seeds, fits["all_s1"][2], fits["all_s5"][2]]
    profiles[label] = scan(label, kept_score, masses, held_seeds)
    candidates[label] = polish(kept_score, profiles[label], held_seeds)
    fits[label] = (kept_score(zero), *candidates[label][0])

    args.out.mkdir(parents=True, exist_ok=True)
    for label, (null, best, point, state) in fits.items():
        check_prediction(point)
        print(f"{label}: null={null:.6f}, best={best:.6f}, gain={null-best:.6f}; "
              f"dm2, e2, mu2, tau2 = {point.tolist()} ({state})", flush=True)
    write_csv(args.out / "fits.csv", [
        dict(case=label, chi2_null=null, chi2_best=best, delta_chi2=null-best,
             status=state, **dict(zip(FIELDS, point)))
        for label, (null, best, point, state) in fits.items()])
    for filename, collection in (("mass_profiles.csv", profiles), ("polished_candidates.csv", candidates)):
        write_csv(args.out / filename, [
            dict(case=label, chi2=value, delta_from_best=value-fits[label][1],
                 status=state, **dict(zip(FIELDS, point)))
            for label, rows in collection.items() for value, point, state in rows])

    cov = covs[1]
    factor = cho_factor(cov[np.ix_(keep, keep)])
    cond_weights = cho_solve(factor, cov[keep, held])
    cond_sigma = np.sqrt(cov[held, held] - cov[held, keep] @ cond_weights)
    rows = []
    for name, point in (("null", zero), ("heldout_refit", fits["holdout_s1"][2]),
                        ("all_bin_best_diagnostic", fits["all_s1"][2])):
        pred = predict(point)
        residual = d[keep] - pred[keep]
        conditional = pred[held] + cond_weights @ residual
        flux_shift = F[held, keep] @ cho_solve(factor, residual)
        pull = (d[held] - conditional)/cond_sigma
        retained = kept_score(point)
        np.testing.assert_allclose(retained + pull**2, scores[1](point), atol=1e-7, rtol=1e-9)
        rows.append(dict(model=name, root_bin=held+1, data=d[held], raw=pred[held],
                         flux_profiled_from_retained=pred[held]+flux_shift,
                         conditional_prediction=conditional, conditional_sigma=cond_sigma,
                         conditional_pull=pull, conditional_chi2=pull**2, retained_chi2=retained))
        print(f"held bin {held+1} {name}: conditional prediction={conditional:.4f}, "
              f"pull={pull:+.4f}; data={d[held]:.4f}", flush=True)
    print(f"Held-bin conditional gain (refit vs null): "
          f"{rows[0]['conditional_chi2']-rows[1]['conditional_chi2']:+.6f}", flush=True)
    write_csv(args.out / "heldout_bin.csv", rows)


if __name__ == "__main__":
    main()

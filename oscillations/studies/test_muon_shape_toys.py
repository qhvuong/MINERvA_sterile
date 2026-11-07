#!/usr/bin/env python3
"""Calibrate the combined 30-bin v4.3.1 shape-only muon scan under the null."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from compare_normalized_flux_profile import BASE_TUNE, TUNES, mean_center, root_inputs, vec, write_csv
from test_muon_holdouts import CORE, U_MAX, amplitude, rate_profiler


def prepare_directions(null_mc, response, amp, L):
    """Whiten each muon rate direction and the unscaled elastic response."""
    from scipy.linalg import solve_triangular

    mc = (null_mc[None, None, :] + response[:, None, :]*amp[None, :, None]).reshape(-1, 30)
    dirs = []
    for lo, hi in ((6, 18), (18, 30)):
        block = np.zeros_like(mc)
        block[:, lo:hi] = mc[:, lo:hi]
        X = solve_triangular(L, block.T, lower=True).T.copy()
        norm = np.linalg.norm(X, axis=1)
        if np.any(norm <= 0.):
            raise ValueError("Muon template has no rate direction.")
        dirs.append(X / norm[:, None])
    x, y = dirs
    delta_elastic = np.zeros_like(mc)
    delta_elastic[:, :6] = mc[:, :6]-null_mc[:6]
    e = solve_triangular(L, delta_elastic.T, lower=True).T.copy()
    rho = np.sum(x*y, axis=1)
    if np.any(1-rho*rho < 1e-10):
        raise ValueError("Muon rate directions are numerically collinear.")
    return x, y, rho, e, np.sum(x*e, axis=1), np.sum(y*e, axis=1), np.sum(e*e, axis=1)


def gain_nonnegative(a, b, rho):
    """Squared norm of the projection onto the cone of positive rates."""
    s = a-rho*b
    t = b-rho*a
    both = (a*s+b*t)/(1-rho*rho)
    edge = np.maximum(np.maximum(a, b), 0.)**2
    return np.where((s >= 0.) & (t >= 0.), both, edge)


def scan_batch(w, design, nmass, namp):
    """Scan identical fine and nested coarse grids for every spectrum."""
    x, y, rho, e, xe, ye, ee = design
    null = gain_nonnegative(w @ x[0], w @ y[0], rho[0])
    best = null.copy()
    coarse = null.copy()
    index = np.zeros(len(w), dtype=int)
    for start in range(0, nmass*namp, 1024):
        stop = min(start+1024, nmass*namp)
        s = slice(start, stop)
        a = x[s] @ w.T - xe[s, None]
        b = y[s] @ w.T - ye[s, None]
        scores = gain_nonnegative(a, b, rho[s, None]) + 2*(e[s] @ w.T) - ee[s, None]
        imax = scores.argmax(axis=0)
        value = scores[imax, np.arange(len(w))]
        improve = value > best
        index[improve] = start+imax[improve]
        best = np.maximum(best, value)
        j = np.arange(start, stop)
        coarse_mask = ((j//namp) % 2 == 0) & ((j % namp) % 2 == 0)
        if np.any(coarse_mask):
            coarse = np.maximum(coarse, scores[coarse_mask].max(axis=0))
    return best-null, coarse-null, index, np.sum(w*w, axis=1)-null


def tail(q, observed):
    from scipy.stats import binomtest

    count = int(np.count_nonzero(q >= observed))
    lo, hi = binomtest(count, len(q)).proportion_ci()
    return dict(observed_q=observed, n_toys=len(q), exceedances=count,
                p_mc=(count+1)/(len(q)+1), tail_ci95_low=lo, tail_ci95_high=hi)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", type=Path, required=True, help="Shape-only muon_holdout_fits.csv")
    p.add_argument("--tune", choices=TUNES, default=BASE_TUNE)
    p.add_argument("--out", type=Path)
    p.add_argument("--toys", type=int, default=10000)
    p.add_argument("--seed", type=int, default=431)
    p.add_argument("--mass-steps", type=int, default=401)
    p.add_argument("--amplitude-steps", type=int, default=129)
    p.add_argument("--batch-size", type=int, default=256)
    args = p.parse_args()
    if args.toys < 1 or args.batch_size < 1 or args.seed < 0:
        p.error("Use positive toy and batch counts and a nonnegative seed.")
    if any(n < 3 or n % 2 == 0 for n in (args.mass_steps, args.amplitude_steps)):
        p.error("Grid sizes must be odd and at least 3 for the paired coarse-grid check.")
    args.cases = ["direct_all"]
    if args.out is None:
        args.out = Path(f"muon_shape_toys_{args.tune}")
    with args.reference.open(newline="") as f:
        refs = {r["model"]: r for r in csv.DictReader(f) if r["case"] == "direct_all"
                and r["holdout"] == "none" and r["normalization"] == "free_per_muon_block"
                and r.get("tune", BASE_TUNE) == args.tune}
    if not {"null", "best"} <= refs.keys():
        p.error(f"Expected combined null and best rows from the {args.tune} shape-only holdout run.")

    _, A, C, d, m, _, checks, (h, cfg, _) = next(root_inputs(args))
    if len(A) != 1000 or not any(r["comparison"] == "MC Flux band"
                                and r["basis"] == "mean_centered" and r["passed"] for r in checks):
        raise RuntimeError("Normalized flux prior did not close to ROOT.")
    cfg = json.loads(cfg.read_text())
    blocks = [np.arange(cfg[k]["start"], cfg[k]["end"]+1) for k in CORE]
    core = np.concatenate(blocks)
    if len(set(core)) != 30 or list(map(len, blocks)) != [6, 12, 12]:
        raise RuntimeError("Unexpected 30-bin core.")
    B = mean_center(A)[:, core]
    V = C[np.ix_(core, core)] + B.T @ B/len(B)
    L = np.linalg.cholesky(V)
    data = d[core]
    from scipy.linalg import solve_triangular

    def predict(mass, u):
        h.OscillateHistogram(float(mass), 0., float(u), 0., False, False)
        return vec(h.GetOscillatedHistogram())[core]

    null_mc = predict(0., 0.)
    np.testing.assert_allclose(null_mc, m[core], rtol=1e-7, atol=1e-8)
    normalize = rate_profiler(data, V, np.arange(30),
                              {"fhc": np.arange(6, 18), "rhc": np.arange(18, 30)})
    fitted_null, scales = normalize(null_mc)
    for name in ("fhc", "rhc"):
        np.testing.assert_allclose(scales[name], float(refs["null"][f"fitted_{name}_scale"]), atol=1e-7)
    for row in (refs["null"], refs["best"]):
        mc = predict(float(row["dm2"]), float(row["umu4_sq"]))
        r = data-normalize(mc)[0]
        np.testing.assert_allclose(r @ np.linalg.solve(V, r), float(row["chi2_retained"]), atol=1e-5)

    masses = np.linspace(0., 100., args.mass_steps)
    amps = np.linspace(0., amplitude(U_MAX), args.amplitude_steps)
    response = np.empty((len(masses), 30))
    for i, mass in enumerate(masses):
        response[i] = (predict(mass, .15)-null_mc) / amplitude(.15)
        if i % 100 == 0 or i == len(masses)-1:
            print(f"ROOT templates: {i+1}/{len(masses)}", flush=True)
    for mass in (1., 19.26, 80.):
        slope = (predict(mass, .15)-null_mc)/amplitude(.15)
        np.testing.assert_allclose(predict(mass, .04), null_mc+amplitude(.04)*slope,
                                   rtol=1e-8, atol=1e-8)
    design = prepare_directions(null_mc, response, amps, L)
    elastic = null_mc.copy()
    elastic[6:] = 0.
    observed_w = solve_triangular(L, data-elastic, lower=True)[None, :]
    q, q_coarse, _, q0 = scan_batch(observed_w, design, len(masses), len(amps))
    np.testing.assert_allclose(q0[0], float(refs["null"]["chi2_retained"]), atol=1e-5)
    q_ref = float(refs["null"]["chi2_retained"])-float(refs["best"]["chi2_retained"])
    if abs(q[0]-q_ref) > .05:
        raise RuntimeError(f"Observed grid q={q[0]:.5f} does not reproduce refined q={q_ref:.5f}. "
                           "Increase --mass-steps and/or --amplitude-steps.")
    print(f"Observed fine q={q[0]:.6f}; coarse q={q_coarse[0]:.6f}; refined q={q_ref:.6f}",
          flush=True)

    # In whitened coordinates, correlated Gaussian fluctuations are iid
    # standard normals. Flux covariance is already in V; do not add it twice.
    mean = solve_triangular(L, fitted_null-elastic, lower=True)
    rng = np.random.default_rng(args.seed)
    draws = np.empty((args.toys, 2))
    best_index = np.empty(args.toys, dtype=int)
    for i in range(0, args.toys, args.batch_size):
        j = min(i+args.batch_size, args.toys)
        w = mean + rng.standard_normal((j-i, len(data)))
        draws[i:j, 0], draws[i:j, 1], best_index[i:j], _ = scan_batch(
            w, design, len(masses), len(amps))
        if i == 0 or j//1000 != i//1000 or j == args.toys:
            print(f"Toys: {j}/{args.toys}", flush=True)

    args.out.mkdir(parents=True, exist_ok=True)
    best_amp = amps[best_index % len(amps)]
    write_csv(args.out / "toy_statistics.csv", [
        dict(toy=i, q_fine=qf, q_coarse=qc, best_dm2=masses[k//len(amps)],
             best_umu4_sq=a/(2*(1+np.sqrt(1-a))))
        for i, ((qf, qc), k, a) in enumerate(zip(draws, best_index, best_amp))])
    summaries = [dict(tune=args.tune, grid=name, mass_steps=nm, amplitude_steps=na, seed=args.seed,
                      observed_refined_q=q_ref, **tail(draws[:, col], obs))
                 for name, nm, na, col, obs in
                 (("fine", len(masses), len(amps), 0, q[0]),
                  ("coarse", (len(masses)+1)//2, (len(amps)+1)//2, 1, q_coarse[0]))]
    write_csv(args.out / "toy_summary.csv", summaries)
    for row in summaries:
        print(f"{row['grid']}: {row['exceedances']}/{args.toys} exceed, "
              f"p={row['p_mc']:.6g}, 95% MC interval "
              f"[{row['tail_ci95_low']:.6g}, {row['tail_ci95_high']:.6g}]", flush=True)
    print(f"Results: {args.out.resolve()}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Nested sterile fits with the full, mean-centered v4.3.1 P8 MC flux prior."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from compare_normalized_flux_profile import (
    NULL, TAGS, mean_center, profile, reference_profile, root_inputs, vec, write_csv,
)

# Fit coordinates: dm2 / 100, then the three squared mixing elements.
MAX = np.array([1., .394, .489, .718])
PARAMS = ("dm2", "ue4_sq", "umu4_sq", "utau4_sq")
FREE = {"mu_only": [0, 2], "no_e": [0, 2, 3], "full": [0, 1, 2, 3]}


def minimize_nested(score, model, seeds, grid_steps, seed):
    from scipy.optimize import LinearConstraint, differential_evolution, minimize, minimize_scalar

    free = FREE[model]
    bounds = [(0., MAX[i]) for i in free]
    constraint = LinearConstraint([[float(i != 0) for i in free]], -np.inf, 1.)

    def expand(x):
        full = np.zeros(4)
        full[free] = x
        return full

    def objective(x):
        return score(expand(x))

    candidates, scan = [], []

    def consider(label, x, success=True):
        x = np.asarray(x, dtype=float)
        if (np.all(np.isfinite(x)) and np.all(x >= -1e-9)
                and np.all(x <= MAX + 1e-9) and np.sum(x[1:]) <= 1. + 1e-9):
            x = x.clip(0., MAX)
            if x[1:].sum() > 1.:
                x[1:] /= x[1:].sum()
            value = score(x)
            if np.isfinite(value):
                candidates.append((value, x, label, bool(success)))

    for i, x in enumerate(seeds):
        consider(f"seed_{i}", expand(np.asarray(x)[free]))

    if model != "full":
        # Trace the restricted profile across mass before floating mass in the
        # final local fits. Multiple starts cover the mu/tau tradeoff.
        for dm in np.linspace(0., 1., grid_steps):
            if model == "mu_only":
                def one(u):
                    return objective([dm, u])

                res = minimize_scalar(one, bounds=bounds[1], method="bounded",
                                      options={"xatol": 1e-7})
                u = min([0., MAX[2], res.x], key=one)
                x, success = expand([dm, u]), bool(res.success)
            else:
                def at_mass(mix):
                    return objective(np.r_[dm, mix])

                mix_cons = LinearConstraint([[1., 1.]], -np.inf, 1.)
                starts = [[0., 0.], [.04, 0.], [.12, .5], [.28, .72]]
                options = [(at_mass(z), np.asarray(z), True) for z in starts[:3]]
                for start in starts:
                    z = minimize(at_mass, np.minimum(start, MAX[2:]), method="SLSQP",
                                 bounds=bounds[1:], constraints=mix_cons,
                                 options={"ftol": 1e-9, "maxiter": 120})
                    if np.all(np.isfinite(z.x)) and z.x.sum() <= 1. + 1e-9:
                        options.append((float(z.fun), z.x, bool(z.success)))
                _, mix, success = min(options, key=lambda z: z[0])
                x = expand(np.r_[dm, mix])
            scan.append((score(x), x, success))

        # Polish distinct mass valleys, rather than the neighboring grid
        # points of just one valley. Keep the endpoints as explicit candidates.
        valleys = [scan[0], scan[-1]]
        valleys += [s for i, s in enumerate(scan[1:-1], 1)
                    if s[0] <= scan[i-1][0] and s[0] <= scan[i+1][0]]
        for i, (_, x, success) in enumerate(sorted(valleys, key=lambda s: s[0])[:12]):
            consider(f"grid_valley_{i}", x, success)

    de = differential_evolution(objective, bounds, constraints=(constraint,),
                                seed=seed, popsize=12, maxiter=160, tol=1e-8,
                                polish=False)
    consider("differential_evolution", expand(de.x), de.success)
    for _, x, label, _ in tuple(candidates):
        result = minimize(objective, x[free], bounds=bounds, constraints=constraint,
                          method="SLSQP", options={"ftol": 1e-11, "maxiter": 300})
        consider(label + "_polished", expand(result.x), result.success)
    if not candidates:
        raise RuntimeError(f"No feasible candidate for {model}")
    return min(candidates, key=lambda z: z[0]), scan, candidates


def parameter_row(x):
    return dict(zip(PARAMS, [100*x[0], *x[1:]]))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", type=Path, required=True,
                   help="data_fits.csv from the mean-centered normalized audit")
    p.add_argument("--out", type=Path, default=Path("electron_mixing_v4p3p1"))
    p.add_argument("--cases", nargs="+", choices=TAGS, default=list(TAGS))
    p.add_argument("--grid-steps", type=int, default=201)
    p.add_argument("--seed", type=int, default=431)
    args = p.parse_args()
    if args.grid_steps < 3 or args.seed < 0:
        p.error("Use at least three mass grid points and a nonnegative seed.")

    with args.reference.open(newline="") as f:
        refs = {r["case"]: r for r in csv.DictReader(f) if r["mode"] == "normalized"
                and r["basis"] == "mean_centered" and int(r["nprof"]) == 1000
                and float(r["lam"]) == 1000.}
    if not set(args.cases) <= refs.keys():
        p.error("Reference CSV must include a centered Nprof=1000, lambda=1000 row for each case.")

    args.out.mkdir(parents=True, exist_ok=True)
    fits, tests, scans, candidates, samples, spectra = [], [], [], [], [], []
    for case, A, C, d, m, _, checks, backend in root_inputs(args):
        # root_inputs has initialized PyROOT/PlotUtils by the first yield.
        # Load SciPy afterwards, as in the working normalized audit, so it
        # cannot load the system C++ runtime before the CVMFS ROOT runtime.
        from scipy.linalg import cho_factor, cho_solve

        if len(A) != 1000 or not any(r["comparison"] == "MC Flux band"
                                    and r["basis"] == "mean_centered" and r["passed"]
                                    for r in checks):
            raise RuntimeError(f"{case}: full centered MC flux prior did not close to ROOT.")
        h, cfg, _ = backend
        groups = json.loads(cfg.read_text())
        bins = [i for info in groups.values() for i in range(info["start"], info["end"]+1)]
        if sorted(bins) != list(range(len(d))):
            raise RuntimeError(f"{case}: sample bins do not partition the spectrum.")
        B = mean_center(A)
        cf = cho_factor(C)
        eff = cho_factor(C + B.T @ B / len(B))

        def prediction(x):
            h.OscillateHistogram(100*x[0], *x[1:], False, False)
            return vec(h.GetOscillatedHistogram())

        def score(x):
            r = d - prediction(x)
            return float(r @ cho_solve(eff, r))

        zero = np.array(NULL)
        np.testing.assert_allclose(prediction(zero), m, rtol=1e-7, atol=1e-8)
        null = score(zero)
        ref = refs[case]
        ref_x = np.array([float(ref[k]) for k in PARAMS])
        ref_x[0] /= 100.
        ref_score = score(ref_x)
        np.testing.assert_allclose([null, ref_score],
                                   [float(ref[k]) for k in ("chi2_null", "chi2_bf")],
                                   atol=1e-5, rtol=0)
        results = {"null": (null, zero, "null", True)}
        seeds = [zero, ref_x, np.array([.02, 0., .1, 0.])]
        for model in FREE:
            best, scan, options = minimize_nested(score, model, seeds, args.grid_steps, args.seed)
            results[model] = best
            seeds.append(best[1])
            scans += [dict(case=case, model=model, chi2=q, success=ok, **parameter_row(x))
                      for q, x, ok in scan]
            candidates += [dict(case=case, model=model, chi2=q, source=label, success=ok,
                                **parameter_row(x)) for q, x, label, ok in options]
            print(f"{case} {model}: chi2={best[0]:.6f} from {best[2]}; {parameter_row(best[1])}",
                  flush=True)

        mu, no_e, full = [results[k][0] for k in FREE]
        if max(mu-null, no_e-mu, full-no_e, full-ref_score) > 1e-6:
            raise RuntimeError(f"{case}: nested minimum ordering failed.")
        tests.append(dict(case=case, q_mu_only=null-mu, q_no_e=null-no_e, q_full=null-full,
                          gain_allowing_tau_at_zero_e=mu-no_e,
                          gain_allowing_e=no_e-full, full_minus_previous=full-ref_score))

        contributions = {}
        for model, (q, x, label, ok) in results.items():
            mc = prediction(x)
            r = d - mc
            total, residual, penalty, shift = profile(C, B, len(B), r)
            np.testing.assert_allclose([total, penalty], reference_profile(C, B, len(B), r),
                                       rtol=1e-8, atol=1e-7)
            np.testing.assert_allclose(q, total, rtol=1e-8, atol=1e-7)
            postfit_r = r - shift
            per_bin = postfit_r * cho_solve(cf, postfit_r)
            np.testing.assert_allclose(per_bin.sum(), residual, rtol=1e-8, atol=1e-7)
            contributions[model] = {name: per_bin[b["start"]:b["end"]+1].sum()
                                    for name, b in groups.items()}
            contributions[model]["flux_penalty"] = penalty
            fits.append(dict(case=case, model=model, chi2=total, residual=residual,
                             penalty=penalty, delta_vs_null=null-total, source=label,
                             success=ok, **parameter_row(x)))
            for name, b in groups.items():
                for i in range(b["start"], b["end"]+1):
                    spectra.append(dict(case=case, model=model, sample=name, bin=i+1,
                                        data=d[i], raw_mc=mc[i], flux_shift=shift[i],
                                        profiled_mc=mc[i]+shift[i]))
        for name in contributions["full"]:
            for left, right in (("null", "mu_only"), ("mu_only", "no_e"), ("no_e", "full")):
                before, after = contributions[left][name], contributions[right][name]
                samples.append(dict(case=case, comparison=f"{left}_to_{right}", component=name,
                                    before=before, after=after, chi2_improvement=before-after))

        for name, records in (("nested_fits", fits), ("nested_comparisons", tests),
                              ("mass_profiles", scans), ("optimizer_candidates", candidates),
                              ("sample_improvements", samples), ("spectra", spectra)):
            write_csv(args.out / f"{name}.csv", records)
    print("Results:", args.out.resolve())


if __name__ == "__main__":
    main()

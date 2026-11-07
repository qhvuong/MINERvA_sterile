#!/usr/bin/env python3
"""Refit v4.3.1 with FHC/RHC electron-sensitive samples held out."""
import argparse
import csv
import inspect
import json
from pathlib import Path

import numpy as np

from compare_normalized_flux_profile import TAGS, mean_center, profile, root_inputs, vec, write_csv
from test_electron_mixing import PARAMS, minimize_nested, parameter_row

BLOCKS = {
    "ratio_all": ("fhc_ratio", "rhc_ratio"),
    "direct_all": ("fhc_nue_selection", "rhc_nue_selection"),
}
HOLDOUTS = {"fhc": (0,), "rhc": (1,), "both": (0, 1)}


def fit_nested(score, model, seeds, grid_steps, seed, report):
    # Older copies of test_electron_mixing.py use the same optimizer without
    # a progress callback. Only pass it when the imported helper accepts it.
    kwargs = {}
    if "report" in inspect.signature(minimize_nested).parameters:
        kwargs["report"] = report
    return minimize_nested(score, model, seeds, grid_steps, seed, **kwargs)


def make_holdout(C, flux, held):
    # Imported lazily so the GPVM ROOT C++ runtime loads before SciPy.
    from scipy.linalg import cho_factor, cho_solve

    V = C + flux
    keep = np.setdiff1d(np.arange(len(V)), held)
    k, hh = np.ix_(keep, keep), np.ix_(held, held)
    hk = np.ix_(held, keep)
    factor = cho_factor(V[k])
    cross = V[hk]
    conditional_cov = V[hh] - cross @ cho_solve(factor, cross.T)
    conditional_factor = cho_factor(conditional_cov)

    def evaluate(residual):
        alpha = cho_solve(factor, residual[keep])
        flux_shift = flux[:, keep] @ alpha
        conditional_shift = cross @ alpha
        leftover = residual[held] - conditional_shift
        return (float(residual[keep] @ alpha),
                float(leftover @ cho_solve(conditional_factor, leftover)),
                flux_shift, conditional_shift)

    def training_score(residual):
        r = residual[keep]
        return float(r @ cho_solve(factor, r))

    return keep, training_score, evaluate, np.sqrt(np.diag(conditional_cov))


def read_references(path):
    refs = {}
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            x = np.array([float(row[k]) for k in PARAMS])
            x[0] /= 100.
            refs[row["case"], row["model"]] = (x, float(row["chi2"]))
    return refs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", type=Path, required=True, help="All-bin nested_fits.csv")
    p.add_argument("--out", type=Path, default=Path("electron_holdouts_v4p3p1"))
    p.add_argument("--cases", nargs="+", choices=TAGS, default=list(TAGS))
    p.add_argument("--holdouts", nargs="+", choices=HOLDOUTS, default=list(HOLDOUTS))
    p.add_argument("--grid-steps", type=int, default=201)
    p.add_argument("--seed", type=int, default=431)
    args = p.parse_args()
    if args.grid_steps < 3 or args.seed < 0:
        p.error("Use at least three mass grid points and a nonnegative seed.")
    refs = read_references(args.reference)
    required = {(case, model) for case in args.cases for model in ("null", "no_e", "full")}
    if required - refs.keys():
        p.error(f"Missing reference rows: {sorted(required - refs.keys())}")
    args.out.mkdir(parents=True, exist_ok=True)
    comparisons, fits, spectra, scans, candidates = [], [], [], [], []

    for case, A, C, d, m, _, checks, backend in root_inputs(args):
        print(f"{case}: checking ROOT flux covariance and all-bin references", flush=True)
        if len(A) != 1000 or not any(r["comparison"] == "MC Flux band"
                                    and r["basis"] == "mean_centered" and r["passed"]
                                    for r in checks):
            raise RuntimeError(f"{case}: full centered MC flux prior did not close to ROOT.")
        h, cfg, _ = backend
        groups = json.loads(cfg.read_text())
        labels = {i: (name, i-b["start"]+1) for name, b in groups.items()
                  for i in range(b["start"], b["end"]+1)}
        if set(labels) != set(range(len(d))):
            raise RuntimeError(f"{case}: sample bins do not cover the spectrum.")
        B = mean_center(A)
        flux = B.T @ B / len(B)
        V = C + flux

        def prediction(x):
            h.OscillateHistogram(100*x[0], *x[1:], False, False)
            return vec(h.GetOscillatedHistogram())

        ref = {k: refs[case, k] for k in ("null", "no_e", "full")}
        np.testing.assert_allclose(prediction(ref["null"][0]), m, rtol=1e-7, atol=1e-8)
        for x, expected in ref.values():
            r = d - prediction(x)
            np.testing.assert_allclose(r @ np.linalg.solve(V, r), expected, rtol=0., atol=1e-5)
        q_all = ref["no_e"][1] - ref["full"][1]

        for key in args.holdouts:
            names = tuple(BLOCKS[case][i] for i in HOLDOUTS[key])
            held = np.array(sorted(i for i, (name, _) in labels.items() if name in names))
            if any(name not in groups for name in names):
                raise RuntimeError(f"{case}: missing holdout sample in {names}")
            keep, training_score, evaluate, cond_sigma = make_holdout(C, flux, held)
            mask = np.ix_(keep, keep)
            meta = dict(case=case, holdout=key, removed="+".join(names))

            def score(x):
                return training_score(d - prediction(x))

            def report(model, message):
                print(f"{case} drop_{key} {model}: {message}", flush=True)

            seeds = [x for x, _ in ref.values()]
            zero = ref["null"][0]
            results = {"null": (score(zero), zero, "null", True)}
            for model in ("no_e", "full"):
                report(model, "starting fit")
                best, grid, options = fit_nested(
                    score, model, seeds, args.grid_steps, args.seed,
                    report=lambda msg: report(model, msg))
                results[model] = best
                seeds.append(best[1])
                scans += [dict(**meta, model=model, chi2=q, success=ok, **parameter_row(x))
                          for q, x, ok in grid]
                candidates += [dict(**meta, model=model, chi2=q, source=source, success=ok,
                                    **parameter_row(x)) for q, x, source, ok in options]
            null, no_e, full = (results[k][0] for k in ("null", "no_e", "full"))
            if max(no_e-null, full-no_e) > 1e-6:
                raise RuntimeError(f"{case} {key}: nested minimum ordering failed.")
            held_scores = {}
            for model, (q, x, source, ok) in results.items():
                mc = prediction(x)
                r = d - mc
                train, held_q, shift, cond_shift = evaluate(r)
                total, residual, penalty, shift_keep = profile(C[mask], B[:, keep], len(B), r[keep])
                np.testing.assert_allclose([train, total], [q, q], rtol=1e-8, atol=1e-7)
                np.testing.assert_allclose(shift[keep], shift_keep, rtol=1e-8, atol=1e-7)
                np.testing.assert_allclose(train+held_q, r @ np.linalg.solve(V, r),
                                           rtol=1e-8, atol=1e-7)
                fits.append(dict(**meta, model=model, chi2_retained=q, residual=residual,
                                 penalty=penalty, chi2_held_conditional=held_q,
                                 source=source, success=ok, **parameter_row(x)))
                held_scores[model] = held_q
                for j, i in enumerate(held):
                    name, local_bin = labels[i]
                    spectra.append(dict(**meta, model=model, sample=name, bin=i+1,
                                        sample_bin=local_bin, data=d[i], raw_mc=mc[i],
                                        flux_only_prediction=mc[i]+shift[i],
                                        conditional_prediction=mc[i]+cond_shift[j],
                                        conditional_sigma=cond_sigma[j]))
            comparisons.append(dict(**meta, retained_bins=len(keep), held_bins=len(held),
                                    gain_allowing_e_all_bins=q_all,
                                    gain_allowing_e_retained=no_e-full,
                                    gain_allowing_e_held_conditional=held_scores["no_e"]-held_scores["full"],
                                    q_no_e_retained=null-no_e, q_full_retained=null-full))
            print(f"{case} drop_{key}: q_e retained={no_e-full:.4f}; "
                  f"held conditional improvement={held_scores['no_e']-held_scores['full']:.4f}; "
                  f"all-bin q_e={q_all:.4f}", flush=True)
            for name, records in (("holdout_comparisons", comparisons), ("holdout_fits", fits),
                                  ("heldout_spectra", spectra), ("mass_profiles", scans),
                                  ("optimizer_candidates", candidates)):
                write_csv(args.out / f"{name}.csv", records)
    print("Results:", args.out.resolve())


if __name__ == "__main__":
    main()

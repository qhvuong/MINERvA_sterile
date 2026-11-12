#!/usr/bin/env python3
"""Profile the unconstrained flavor mixing using only retained direct-CC data."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from audit_direct_cc_holdouts import HOLDOUTS, SAMPLES
from test_electron_holdouts import make_holdout
from test_electron_mixing import MAX, PARAMS, parameter_row
from test_flux_mode_convergence import current_root_inputs, mode_basis, vec, write_csv


def fit_fixed_mixing(score, coordinate, value, seeds, seed):
    """Separate profile-study search, with fixed mixing and free mass."""
    from scipy.optimize import LinearConstraint, differential_evolution, minimize

    free = [i for i in range(4) if i != coordinate]
    bounds = [(0., MAX[i]) for i in free]
    constraint = LinearConstraint([[float(i != 0) for i in free]], -np.inf, 1.-value)

    def expand(z):
        point = np.zeros(4)
        point[coordinate] = value
        point[free] = z
        return point

    def objective(z):
        return score(expand(z))

    candidates = []
    def consider(z, success):
        point = expand(np.asarray(z))
        if (np.isfinite(point).all() and np.all(point >= -1e-9) and
                np.all(point <= MAX+1e-9) and point[1:].sum() <= 1.+1e-9):
            point = np.clip(point, 0., MAX)
            others = [i for i in (1, 2, 3) if i != coordinate]
            if point[others].sum() > 1.-value:
                point[others] *= (1.-value)/point[others].sum()
            q = score(point)
            if np.isfinite(q):
                candidates.append((q, point, bool(success)))

    starts = [np.asarray(x).copy() for x in seeds]
    starts.extend(np.array([mass, .02, .04, 0.]) for mass in (.01, .10, .19, .5, .96))
    for point in starts:
        point[coordinate] = value
        others = [i for i in (1, 2, 3) if i != coordinate]
        point[others] *= min(1., (1.-value)/max(point[others].sum(), 1e-12))
        consider(point[free], True)
    result = differential_evolution(objective, bounds, constraints=(constraint,), seed=seed,
                                    popsize=10, maxiter=120, tol=1e-8, polish=False)
    consider(result.x, result.success)
    for _, point, _ in tuple(candidates):
        result = minimize(objective, point[free], method="SLSQP", bounds=bounds,
                          constraints=(constraint,), options={"ftol": 1e-10, "maxiter": 250})
        consider(result.x, result.success)
    if not candidates:
        raise RuntimeError("No feasible fixed-mixing solution")
    return min(candidates, key=lambda row: row[0])


def plot_results(rows, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    for column, holdout in enumerate(("electron", "muon")):
        subset = [row for row in rows if row["holdout"] == holdout]
        x = [row["fixed_mixing"] for row in subset]
        axes[0, column].plot(x, [row["delta_train"] for row in subset], ".-")
        axes[0, column].set(title=f"{holdout.capitalize()} samples held out",
                            ylabel=r"Training $\chi^2$ above lowest found")
        axes[1, column].plot(x, [row["gain_held_conditional"] for row in subset], ".-")
        axes[1, column].axhline(0., color="black", linewidth=.8)
        axes[1, column].set(xlabel=r"$|U_{e4}|^2$" if column == 0 else r"$|U_{\mu4}|^2$",
                            ylabel=r"Conditional held-out gain over null")
    fig.savefig(out / "holdout_mixing_profiles.png", dpi=170)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path, help="Production sterile_refits.csv")
    parser.add_argument("--holdout-fits", required=True, type=Path, help="Completed audit holdout_fits.csv")
    parser.add_argument("--k", default=13, type=int)
    parser.add_argument("--tag-suffix", default="_updatedStat_newEL")
    parser.add_argument("--steps", type=int, default=31)
    parser.add_argument("--seed", type=int, default=431)
    parser.add_argument("--plots", action="store_true")
    parser.add_argument("--out", type=Path, default=Path("holdout_mixing_profiles_k13"))
    args = parser.parse_args()
    if not 1 <= args.k <= 54 or args.steps < 3 or args.seed < 0:
        parser.error("Use K=1..54, at least three grid steps and a nonnegative seed.")
    with args.reference.open(newline="") as stream:
        references = [r for r in csv.DictReader(stream) if r["case"] == "direct_all" and int(r["k"]) == args.k]
    if len(references) != 1:
        parser.error("Expected one direct_all reference for the requested K.")
    ref = references[0]
    def point(row):
        x = np.array([float(row[name]) for name in PARAMS])
        x[0] /= 100.
        return x
    xref = point(ref)
    with args.holdout_fits.open(newline="") as stream:
        fits = {(r["holdout"], r["model"]): r for r in csv.DictReader(stream)}
    args.cases, args.sterile_k = ["direct_all"], [args.k]
    _, A, C, data, mc, hist, _, _, _ = next(current_root_inputs(args))
    config = json.loads(Path(hist.hist_config).read_text())
    blocks = {name: np.arange(config[name]["start"], config[name]["end"]+1) for name in SAMPLES}
    np.testing.assert_array_equal(np.sort(np.concatenate(list(blocks.values()))), np.arange(len(data)))
    _, _, _, F = mode_basis(A, C)
    flux = F[:args.k].T @ F[:args.k]
    V = C+flux
    def prediction(x):
        hist.OscillateHistogram(100*x[0], *x[1:], False, False)
        return vec(hist.GetOscillatedHistogram())
    np.testing.assert_allclose(prediction(np.zeros(4)), mc, rtol=1e-7, atol=1e-8)
    for x, key in ((np.zeros(4), "chi2_null"), (xref, "chi2_bf")):
        r = data-prediction(x)
        np.testing.assert_allclose(r @ np.linalg.solve(V, r), float(ref[key]), atol=1e-5, rtol=0.)

    args.out.mkdir(parents=True, exist_ok=True)
    rows, benchmarks = [], []
    for holdout, coordinate in (("electron", 1), ("muon", 2)):
        held = np.sort(np.concatenate([blocks[name] for name in HOLDOUTS[holdout]]))
        _, train, evaluate, sigma = make_holdout(C, flux, held)
        best_row = fits[holdout, "training_best"]
        xbest = point(best_row)
        def score(x):
            return train(data-prediction(x))
        for model, x in (("null", np.zeros(4)), ("training_best", xbest), ("all_fit_reference", xref)):
            q, h, _, _ = evaluate(data-prediction(x))
            expected = fits[holdout, model]
            np.testing.assert_allclose([q, h], [float(expected["chi2_train"]), float(expected["chi2_held_conditional"])],
                                       atol=1e-5, rtol=0.)
            benchmarks.append(dict(holdout=holdout, model=model, chi2_train=q, chi2_held_conditional=h,
                                   **parameter_row(x)))
        null_held = float(fits[holdout, "null"]["chi2_held_conditional"])
        grid = np.unique(np.r_[np.linspace(0., MAX[coordinate], args.steps), xref[coordinate], xbest[coordinate]])
        solutions = {}
        warm = xbest
        # The reverse pass checks for a different mass basin missed by continuation.
        for pass_id, values in enumerate((grid, grid[::-1])):
            for index, value in enumerate(values):
                seeds = [xref, xbest, warm]
                if float(value) in solutions:
                    seeds.append(solutions[float(value)][1])
                result = fit_fixed_mixing(score, coordinate, float(value), seeds,
                                          args.seed+pass_id*1000+index)
                if float(value) not in solutions or result[0] < solutions[float(value)][0]:
                    solutions[float(value)] = result
                warm = solutions[float(value)][1]
                print(f"{holdout}: pass {pass_id+1}, point {index+1}/{len(grid)}, mixing={value:.6g}, "
                      f"training chi2={solutions[float(value)][0]:.6f}", flush=True)
            current, current_spectra = [], []
            baseline = min(float(best_row["chi2_train"]), *(q for q, _, _ in solutions.values()))
            for value, (q, x, success) in sorted(solutions.items()):
                pred = prediction(x)
                tq, hq, _, shift = evaluate(data-pred)
                np.testing.assert_allclose(tq, q, atol=1e-7, rtol=1e-9)
                current.append(dict(holdout=holdout, fixed_mixing=value, chi2_train=q, delta_train=q-baseline,
                                    chi2_held_conditional=hq, gain_held_conditional=null_held-hq,
                                    optimizer="profile-study", success=success, **parameter_row(x)))
                current_spectra.extend(dict(holdout=holdout, fixed_mixing=value, bin=int(i)+1, data=data[i],
                                            conditional_prediction=pred[i]+shift[j], conditional_sigma=sigma[j])
                                       for j, i in enumerate(held))
            write_csv(args.out / f"{holdout}_profile.csv", current)
            write_csv(args.out / f"{holdout}_spectra.csv", current_spectra)
        rows.extend(current)
    write_csv(args.out / "benchmarks.csv", benchmarks)
    write_csv(args.out / "mixing_profiles.csv", rows)
    if args.plots:
        plot_results(rows, args.out)


if __name__ == "__main__":
    main()

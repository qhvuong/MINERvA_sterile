#!/usr/bin/env python3
"""Direct-CC sample holdouts using the production optimizer and fixed P8 modes."""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from test_electron_holdouts import make_holdout
from test_electron_mixing import MAX, PARAMS, parameter_row
from test_flux_mode_convergence import current_root_inputs, mass_profile, mode_basis, vec, write_csv


SAMPLES = ("fhc_elastic", "fhc_nue_selection", "rhc_nue_selection",
           "fhc_numu_selection", "rhc_numu_selection")
HOLDOUTS = dict(
    electron=SAMPLES[1:3],
    muon=SAMPLES[3:5],
    elastic=SAMPLES[:1],
    fhc_electron=SAMPLES[1:2],
    rhc_electron=SAMPLES[2:3],
    fhc_muon=SAMPLES[3:4],
    rhc_muon=SAMPLES[4:5],
    fhc_cc=(SAMPLES[1], SAMPLES[3]),
    rhc_cc=(SAMPLES[2], SAMPLES[4]),
)


def sample_attribution(V, residuals, blocks):
    """Average each block's incremental gain over all sample orderings."""
    from scipy.linalg import cho_factor, cho_solve

    n = len(blocks)
    gains = {0: 0.}
    for mask in range(1, 1 << n):
        keep = np.concatenate([block for i, block in enumerate(blocks) if mask & (1 << i)])
        r = residuals[:, keep]
        factor = cho_factor(V[np.ix_(keep, keep)])
        q = np.sum(r * cho_solve(factor, r.T).T, axis=1)
        gains[mask] = float(q[0] - q[1])
    values = []
    for i in range(n):
        bit = 1 << i
        values.append(sum((gains[mask | bit] - gain) /
                          (n * math.comb(n-1, bin(mask).count("1")))
                          for mask, gain in gains.items() if not mask & bit))
    np.testing.assert_allclose(sum(values), gains[(1 << n)-1], atol=1e-8, rtol=1e-9)
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True,
                        help="sterile_refits.csv with the production direct_all K=13 row")
    parser.add_argument("--tag-suffix", default="_updatedStat_newEL")
    parser.add_argument("--k", type=int, default=13)
    parser.add_argument("--holdouts", nargs="+", choices=HOLDOUTS, default=list(HOLDOUTS))
    parser.add_argument("--out", type=Path, default=Path("direct_cc_holdouts_k13"))
    args = parser.parse_args()
    if not 0 <= args.k <= 54:
        parser.error("--k must be between 0 and 54.")
    with args.reference.open(newline="") as stream:
        rows = [row for row in csv.DictReader(stream)
                if row["case"] == "direct_all" and int(row["k"]) == args.k]
    if len(rows) != 1:
        parser.error("Expected one direct_all reference row for the requested K.")
    ref = rows[0]
    ref_point = np.array([float(ref[name]) for name in PARAMS])
    ref_point[0] /= 100.
    args.cases = ["direct_all"]
    args.sterile_k = [args.k]
    _, A, C, d, m, hist, _, _, _ = next(current_root_inputs(args))

    # Import after ROOT initialization on GPVM.
    from scipy.optimize import LinearConstraint, minimize
    from tools.Fitters import OscillationFitter

    cfg = json.loads(Path(hist.hist_config).read_text())
    blocks = {name: np.arange(cfg[name]["start"], cfg[name]["end"]+1) for name in SAMPLES}
    np.testing.assert_array_equal(np.sort(np.concatenate(list(blocks.values()))), np.arange(len(d)))
    _, _, _, F = mode_basis(A, C)
    flux = F[:args.k].T @ F[:args.k]
    V = C + flux

    def prediction(x):
        hist.OscillateHistogram(100*x[0], *x[1:], False, False)
        return vec(hist.GetOscillatedHistogram())

    zero = np.zeros(4)
    np.testing.assert_allclose(prediction(zero), m, rtol=1e-7, atol=1e-8)
    ref_residuals = np.array([d-m, d-prediction(ref_point)])
    q_ref = [float(r @ np.linalg.solve(V, r)) for r in ref_residuals]
    np.testing.assert_allclose(q_ref, [float(ref["chi2_null"]), float(ref["chi2_bf"])],
                               rtol=0., atol=1e-5)

    args.out.mkdir(parents=True, exist_ok=True)
    write_csv(args.out / "all_fit_sample_attribution.csv",
              [dict(sample=name, gain_averaged_over_sample_orderings=value)
               for name, value in zip(SAMPLES, sample_attribution(V, ref_residuals, list(blocks.values())))])
    fits, summary, spectra = [], [], []
    for name in dict.fromkeys(args.holdouts):
        held = np.sort(np.concatenate([blocks[s] for s in HOLDOUTS[name]]))
        keep, train, evaluate, sigma = make_holdout(C, flux, held)

        def score(x):
            return train(d-prediction(x))

        class HeldOutFitter(OscillationFitter):
            def CalChi2(self, x):
                return score(x)

        print(f"{name}: fitting {len(keep)} bins and holding out {len(held)}", flush=True)
        fitter = HeldOutFitter(hist, exclude="", marginalize_flux=False, hist_config=hist.hist_config)
        best, params = fitter.DoFit()
        x = np.array([params["m"]/100., params["ue4"], params["umu4"], params["utau4"]])
        np.testing.assert_allclose(best, score(x), rtol=1e-9, atol=1e-8)
        anchored, x_anchor = mass_profile(score, x, ref_point, [100*ref_point[0]])[0]
        # Polish the anchored solution with free mass on the same training objective.
        polish = minimize(score, x_anchor, method="SLSQP", bounds=list(zip(np.zeros(4), MAX)),
                          constraints=LinearConstraint([[0., 1., 1., 1.]], -np.inf, 1.),
                          options={"ftol": 1e-10, "maxiter": 300})
        if (np.all(polish.x >= 0) and np.all(polish.x <= MAX) and
                polish.x[1:].sum() <= 1. and polish.fun < anchored):
            if polish.fun < best:
                x, best = polish.x, float(polish.fun)
        if anchored < best:
            x, best = x_anchor, anchored

        scores = {}
        for model, point in (("null", zero), ("training_best", x),
                             ("training_fixed_reference_mass", x_anchor), ("all_fit_reference", ref_point)):
            mc = prediction(point)
            r = d-mc
            train_q, held_q, flux_shift, cond_shift = evaluate(r)
            np.testing.assert_allclose(train_q+held_q, r @ np.linalg.solve(V, r), rtol=1e-9, atol=1e-7)
            fits.append(dict(holdout=name, model=model, chi2_train=train_q,
                             chi2_held_conditional=held_q, **parameter_row(point)))
            scores[model] = (train_q, held_q)
            for j, i in enumerate(held):
                spectra.append(dict(holdout=name, model=model, bin=int(i)+1, data=d[i],
                                    cv=mc[i], flux_only_prediction=mc[i]+flux_shift[i],
                                    conditional_prediction=mc[i]+cond_shift[j], conditional_sigma=sigma[j]))
        q0, h0 = scores["null"]
        qb, hb = scores["training_best"]
        qa, ha = scores["training_fixed_reference_mass"]
        qr, hr = scores["all_fit_reference"]
        summary.append(dict(holdout=name, samples="+".join(HOLDOUTS[name]), n_train=len(keep),
                            n_held=len(held), gain_train=q0-qb, gain_held_conditional=h0-hb,
                            anchored_gain_train=q0-qa, anchored_gain_held_conditional=h0-ha,
                            all_fit_gain_train=q0-qr, all_fit_gain_held_conditional=h0-hr,
                            **parameter_row(x)))
        print(f"{name}: train gain={q0-qb:+.4f}; held gain={h0-hb:+.4f}; "
              f"held gain at fixed reference mass={h0-ha:+.4f}", flush=True)
        write_csv(args.out / "holdout_summary.csv", summary)
        write_csv(args.out / "holdout_fits.csv", fits)
        write_csv(args.out / "heldout_spectra.csv", spectra)


if __name__ == "__main__":
    main()

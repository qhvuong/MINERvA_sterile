#!/usr/bin/env python3
"""Cross-score the logged v4.3.1 candidates using the production oscillation and flux solve."""
import argparse
import csv
import json
import os
from pathlib import Path
import sys

import numpy as np

POINTS = {
    "null": (0., 0., 0., 0.),
    "ratio_all_BF": (0.2842873055830976, .394, .489, 0.),
    "ratio_excluded_BF": (21.265051982341323, 0.015353785709121017,
                          0.06169659950574969, 0.),
    "direct_BF": (18.8770015348253, 0.05364601633383929,
                  0.026969706695506964, .718),
}
LOGGED = {
    ("ratio_all", "null"): 47.179266,
    ("ratio_excluded", "null"): 108.324207,
    ("direct_all", "null"): 48.716424,
    ("ratio_all", "ratio_all_BF"): 32.15929107178868,
    ("ratio_excluded", "ratio_excluded_BF"): 43.79023714453416,
    ("direct_all", "direct_BF"): 40.86446802443094,
}


def vec(h):
    return np.array([h.GetBinContent(i) for i in range(1, h.GetNbinsX()+1)])


def partition(C, r, other, interest):
    """Exact Schur split of r.T @ inv(C) @ r."""
    B = C[np.ix_(other, other)]
    X = C[np.ix_(other, interest)]
    base = np.linalg.solve(B, r[other])
    cond = r[interest] - X.T @ base
    S = C[np.ix_(interest, interest)] - X.T @ np.linalg.solve(B, X)
    return float(r[other] @ base), float(cond @ np.linalg.solve(S, cond))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=Path("cross_score_v4p3p1.csv"))
    args = p.parse_args()
    if any(os.environ.get(k) == "1" for k in
           ("USE_RELATIVE_FLUX_A", "USE_HYBRID_EXCLUDED_FLUX_COV")):
        p.error("Unset USE_RELATIVE_FLUX_A and USE_HYBRID_EXCLUDED_FLUX_COV for fixed-A scoring.")
    osc = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(osc.parent))
    os.environ.setdefault("CCNUEROOT", str(osc.parent))
    sys.argv = [sys.argv[0]]
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    from tools.StitchedHistogram import StitchedHistogram
    from tools.Fitters import Statistics

    rows = []
    for direct in (False, True):
        tag = "predNueel_{}p8_v4p3p1".format("noRatio_" if direct else "")
        cfg_path = osc / ("HIST_CONFIG_" + tag + ".json")
        cfg = json.loads(cfg_path.read_text())
        h = StitchedHistogram("sample")
        h.SetHistConfig(str(cfg_path))
        h.Load(str(osc / "rootfiles" / ("NuE_stitched_hists_" + tag + ".root")))
        d = vec(h.GetDataHistogram())
        C = np.asarray(h.GetCovarianceMatrix(sansFlux=True))
        names = np.full(len(d), "", dtype=object)
        for name, lim in cfg.items():
            names[lim["start"]:lim["end"]+1] = name
        interest = np.flatnonzero(np.isin(names, ("fhc_nue_selection", "rhc_nue_selection")
                                            if direct else ("fhc_ratio", "rhc_ratio")))
        other = np.setdiff1d(np.arange(len(d)), interest)
        if len(interest) == 0 or np.any(names == ""):
            raise ValueError("Unexpected sample map: " + str(cfg_path))
        ratio_bins = np.flatnonzero(vec(h.ratio_id) > .5)
        if (direct and len(ratio_bins)) or (not direct and not np.array_equal(ratio_bins, interest)):
            raise ValueError("ROOT ratio mode disagrees with sample map")
        cases = [("direct_all", "")] if direct else [("ratio_all", ""), ("ratio_excluded", "ratio")]
        for label, theta in POINTS.items():
            use_osc = label != "null"
            if use_osc:
                h.OscillateHistogram(*theta, False, False)
            for case, exclude in cases:
                stat = Statistics(h, exclude=exclude, lam=1.,
                                  profile_n_universes=50, hist_config=str(cfg_path))
                total, penalty = stat.Chi2DataMC(marginalize=True, useOsc=use_osc)
                ff = stat.GetFluxFitter(useOsc=use_osc)
                m, _ = ff.MarginalizeFlux()
                core, cond = partition(C, d-np.asarray(m), other, interest)
                if not np.isclose(total, core+cond+penalty, rtol=1e-7, atol=1e-7):
                    raise ValueError("Schur split does not close for {} {}".format(case, label))
                logged = LOGGED.get((case, label))
                print("{:<15} {:<20} chi2={:9.5f} core={:8.4f} conditional={:8.4f} penalty={:8.4f}".format(
                    case, label, total, core, cond, penalty), flush=True)
                rows.append(dict(case=case, point=label, dm2=theta[0], ue4_sq=theta[1],
                                 umu4_sq=theta[2], utau4_sq=theta[3], chi2=total, core=core,
                                 conditional_ratio_or_direct_e=cond, penalty=penalty,
                                 logged_chi2=logged, difference_from_log=None if logged is None else total-logged))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    mismatches = [r for r in rows if r["difference_from_log"] is not None
                  and abs(r["difference_from_log"]) > 1e-4]
    print("Wrote", args.out, "| logged-score mismatches:", len(mismatches))
    return int(bool(mismatches))


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3

"""
Numerical closure test for JOINT and MENU1PI tune weights.

Baseline:
    v1_D2_LECOH_DiffCOH

Predictions:
    JOINT:
        baseline * JOINT diagnostic weight

    MENU1PI:
        baseline * MENU1PI diagnostic weight

Compare reconstructed predictions with actual tune output.

Important:
    This is histogram-level closure, not exact event-level closure.
    Weight values are approximated using the TH2 Y-bin centers.

Run:
    python3 studies/check_mnv_tune_closure.py --mode all
"""

import argparse
import csv
import math
from pathlib import Path

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)


# ============================================================
# INPUT CONFIGURATION
# ============================================================

SAMPLES = {
    "FHC": {
        "directory": "/exp/minerva/data/users/qvuong/nu_mu",
        "playlist": "leFHC",
        "sample": "CCnumu",
    },
    "RHC": {
        "directory": "/exp/minerva/data/users/qvuong/antinu_mu",
        "playlist": "le5",
        "sample": "CCnumubar",
    },
}

BASELINE = "v1_D2_LECOH_DiffCOH"

TESTS = {
    "JOINT": {
        "tag": "v1_D2_JOINT_LECOH_DiffCOH",
        "map": "JOINTWeight_trueEnu",
    },
    "MENU1PI": {
        "tag": "v4p3p0_D2_MENU1PI_LECOH_DiffCOH",
        "map": "MENU1PIWeight_trueEnu",
    },
}

ENERGY_HIST = "Flux"

PEAK_LOW = 2.5
PEAK_HIGH = 3.8


# ============================================================
# FILE ACCESS
# ============================================================

def input_path(mode, tag):

    c = SAMPLES[mode]

    return Path(c["directory"]) / (
        f"kin_dist_mc{c['playlist']}_"
        f"Validation_{c['sample']}_{tag}_MAD.root"
    )


def open_file(mode, tag):

    path = input_path(mode, tag)

    if not path.exists():
        raise FileNotFoundError(path)

    f = ROOT.TFile.Open(str(path), "READ")

    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {path}")

    return f


def get_hist(f, name, dimension):

    obj = f.Get(name)

    if not obj:
        raise RuntimeError(
            f"Histogram '{name}' not found in {f.GetName()}"
        )

    if hasattr(obj, "GetCVHistoWithStatError"):
        h = obj.GetCVHistoWithStatError()
    else:
        h = obj.Clone()

    if not h:
        raise RuntimeError(f"Cannot clone histogram {name}")

    actual_dim = 2 if h.InheritsFrom("TH2") else 1

    if actual_dim != dimension:
        raise RuntimeError(
            f"{name}: expected {dimension}D, "
            f"found {actual_dim}D"
        )

    h.SetDirectory(0)

    return h


def get_pot(f):

    tree = f.Get("Meta")

    if not tree:
        raise RuntimeError(
            f"No Meta tree in {f.GetName()}"
        )

    if not tree.GetBranch("POT_Used"):
        raise RuntimeError(
            f"No POT_Used branch in {f.GetName()}"
        )

    pot = sum(float(e.POT_Used) for e in tree)

    if pot <= 0:
        raise RuntimeError(
            f"Invalid POT: {pot}"
        )

    return pot


# ============================================================
# BINNING VALIDATION
# ============================================================

def same_x_binning(a, b):

    if a.GetNbinsX() != b.GetNbinsX():
        return False

    ax = a.GetXaxis()
    bx = b.GetXaxis()

    for i in range(1, a.GetNbinsX() + 2):

        xa = ax.GetBinLowEdge(i)
        xb = bx.GetBinLowEdge(i)

        if not math.isclose(
            xa, xb, rel_tol=1e-9, abs_tol=1e-7
        ):
            return False

    return True


# ============================================================
# REWEIGHT BASELINE DIAGNOSTIC MAP
# ============================================================

def predict_from_map(h2):

    """
    X axis: true neutrino energy
    Y axis: independently evaluated tune weight
    Z axis: baseline weighted event yield

    For every energy bin:

      N_baseline_from_map = sum_y N(x,y)

      N_predicted = sum_y N(x,y) * weight_y

      Delta_predicted = sum_y N(x,y) * (weight_y - 1)

    Uses bin centers as approximate weight values.
    """

    results = []

    yaxis = h2.GetYaxis()

    for ix in range(1, h2.GetNbinsX() + 1):

        low = h2.GetXaxis().GetBinLowEdge(ix)
        high = h2.GetXaxis().GetBinUpEdge(ix)

        baseline = 0.0
        predicted = 0.0

        suppressed = 0.0
        enhanced = 0.0
        unity = 0.0

        for iy in range(1, h2.GetNbinsY() + 1):

            content = h2.GetBinContent(ix, iy)

            weight = yaxis.GetBinCenter(iy)

            baseline += content
            predicted += content * weight

            if weight < 1.0 - 1e-8:
                suppressed += content

            elif weight > 1.0 + 1e-8:
                enhanced += content

            else:
                unity += content

        # Entries outside the diagnostic weight-axis range
        # cannot be assigned a trustworthy weight.
        flow_content = sum(
            h2.GetBinContent(ix, iy)
            for iy in (0, h2.GetNbinsY() + 1)
        )

        delta = predicted - baseline

        mean_weight = (
            predicted / baseline
            if baseline > 0
            else float("nan")
        )

        results.append({
            "energy_low": low,
            "energy_high": high,
            "map_baseline": baseline,
            "predicted": predicted,
            "predicted_delta": delta,
            "mean_weight": mean_weight,
            "suppressed": suppressed,
            "enhanced": enhanced,
            "unity": unity,
            "weight_flow": flow_content,
        })

    return results


# ============================================================
# COMPARISON
# ============================================================

def compare(mode, test_name, test, output_dir):

    print("\n" + "=" * 110)
    print(f"{mode} | {test_name} WEIGHT CLOSURE")
    print("=" * 110)

    print("Baseline:", BASELINE)
    print("Tune:    ", test["tag"])
    print("Map:     ", test["map"])

    f_base = open_file(mode, BASELINE)
    f_tune = open_file(mode, test["tag"])

    try:

        hbase = get_hist(f_base, ENERGY_HIST, 1)
        htune = get_hist(f_tune, ENERGY_HIST, 1)

        hmap = get_hist(f_base, test["map"], 2)

        pot_base = get_pot(f_base)
        pot_tune = get_pot(f_tune)

        htune.Scale(pot_base / pot_tune)

        print(f"\nBaseline POT: {pot_base:.8e}")
        print(f"Tune POT:     {pot_tune:.8e}")
        print(f"POT scale:    {pot_base / pot_tune:.8f}")

        if not same_x_binning(hbase, htune):
            raise RuntimeError(
                "Baseline and tuned Flux binning differ"
            )

        if not same_x_binning(hbase, hmap):
            raise RuntimeError(
                "Flux and diagnostic weight-map X binning differ.\n"
                "Need a common-binning projection or new histograms."
            )

        predicted_bins = predict_from_map(hmap)

        rows = []

        print("\nBIN-BY-BIN RESULTS")
        print("-" * 110)

        print(
            f"{'Energy':>14}"
            f"{'Base':>13}"
            f"{'Observed':>13}"
            f"{'Predicted':>13}"
            f"{'Obs dN':>13}"
            f"{'Pred dN':>13}"
            f"{'Residual':>13}"
            f"{'<w>':>9}"
        )

        for i, prediction in enumerate(predicted_bins, 1):

            low = prediction["energy_low"]
            high = prediction["energy_high"]

            base = hbase.GetBinContent(i)
            tuned = htune.GetBinContent(i)

            observed_delta = tuned - base
            predicted_delta = prediction["predicted_delta"]

            predicted_total = base + predicted_delta

            residual = observed_delta - predicted_delta

            map_baseline = prediction["map_baseline"]

            coverage = (
                map_baseline / base
                if base > 0
                else float("nan")
            )

            row = {
                "mode": mode,
                "test": test_name,
                "energy_low": low,
                "energy_high": high,

                "baseline": base,
                "observed": tuned,
                "predicted": predicted_total,

                "observed_delta": observed_delta,
                "predicted_delta": predicted_delta,

                "residual": residual,

                "observed_change_percent": (
                    100.0 * observed_delta / base
                    if base > 0 else float("nan")
                ),

                "predicted_change_percent": (
                    100.0 * predicted_delta / base
                    if base > 0 else float("nan")
                ),

                "mean_diagnostic_weight":
                    prediction["mean_weight"],

                "map_baseline": map_baseline,
                "map_coverage": coverage,

                "suppressed_yield":
                    prediction["suppressed"],

                "enhanced_yield":
                    prediction["enhanced"],

                "unity_yield":
                    prediction["unity"],

                "weight_flow":
                    prediction["weight_flow"],
            }

            rows.append(row)

            print(
                f"{low:5.2f}-{high:<7.2f}"
                f"{base:13.2f}"
                f"{tuned:13.2f}"
                f"{predicted_total:13.2f}"
                f"{observed_delta:13.2f}"
                f"{predicted_delta:13.2f}"
                f"{residual:13.2f}"
                f"{prediction['mean_weight']:9.4f}"
            )

        # ----------------------------------------------------
        # Integrated summaries
        # ----------------------------------------------------

        regions = [
            ("Full visible spectrum", None, None),
            ("Flux peak 2.5-3.8 GeV", 2.5, 3.8),
            ("2.0-5.0 GeV", 2.0, 5.0),
            ("Above 5.0 GeV", 5.0, None),
        ]

        print("\n" + "=" * 110)
        print("INTEGRATED REGIONS")
        print("=" * 110)

        for name, emin, emax in regions:

            selected = []

            for row in rows:

                low = row["energy_low"]
                high = row["energy_high"]

                if emin is not None and low < emin - 1e-7:
                    continue

                if emax is not None and high > emax + 1e-7:
                    continue

                selected.append(row)

            if not selected:
                continue

            base = sum(r["baseline"] for r in selected)
            observed = sum(r["observed"] for r in selected)

            observed_delta = sum(
                r["observed_delta"] for r in selected
            )

            predicted_delta = sum(
                r["predicted_delta"] for r in selected
            )

            residual = observed_delta - predicted_delta

            observed_pct = (
                100 * observed_delta / base
                if base else float("nan")
            )

            predicted_pct = (
                100 * predicted_delta / base
                if base else float("nan")
            )

            print(f"\n{name}")
            print(f"  Baseline yield:        {base:15.3f}")
            print(f"  Tuned yield:           {observed:15.3f}")

            print(f"  Observed delta:        {observed_delta:15.3f}")
            print(f"  Predicted delta:       {predicted_delta:15.3f}")

            print(f"  Observed change:       {observed_pct:14.4f}%")
            print(f"  Predicted change:      {predicted_pct:14.4f}%")

            print(f"  Residual:              {residual:15.3f}")

            if observed_delta != 0:
                print(
                    "  Residual / obs delta: "
                    f"{100 * residual / observed_delta:14.4f}%"
                )

            suppressed = sum(
                r["suppressed_yield"] for r in selected
            )
            enhanced = sum(
                r["enhanced_yield"] for r in selected
            )
            unity = sum(
                r["unity_yield"] for r in selected
            )

            total = suppressed + enhanced + unity

            print("\n  Baseline-weighted diagnostic populations:")

            if total > 0:
                print(
                    f"    Weight < 1:         "
                    f"{100 * suppressed / total:10.3f}%"
                )
                print(
                    f"    Weight ~ 1:         "
                    f"{100 * unity / total:10.3f}%"
                )
                print(
                    f"    Weight > 1:         "
                    f"{100 * enhanced / total:10.3f}%"
                )

            coverage_base = sum(
                r["map_baseline"] for r in selected
            )

            flow = sum(
                r["weight_flow"] for r in selected
            )

            coverage = (
                100 * coverage_base / base
                if base else float("nan")
            )

            print(
                f"\n  Diagnostic map coverage: {coverage:.4f}%"
            )
            print(
                f"  Diagnostic Y-flow yield: {flow:.3f}"
            )

            if abs(coverage - 100.0) > 0.5:
                print(
                    "  [WARNING] Map population does not fully "
                    "match the baseline Flux population."
                )

            if abs(flow) > 1e-6:
                print(
                    "  [WARNING] Diagnostic weight overflow/"
                    "underflow exists. Predicted closure "
                    "may be incomplete."
                )

        # ----------------------------------------------------
        # CSV output
        # ----------------------------------------------------

        output_dir.mkdir(parents=True, exist_ok=True)

        csv_file = output_dir / (
            f"{mode}_{test_name}_weight_closure.csv"
        )

        with open(csv_file, "w", newline="") as f:

            writer = csv.DictWriter(
                f,
                fieldnames=list(rows[0].keys()),
            )

            writer.writeheader()
            writer.writerows(rows)

        print("\nSaved:", csv_file)

    finally:
        f_base.Close()
        f_tune.Close()


# ============================================================
# MAIN
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--mode",
        choices=("FHC", "RHC", "all"),
        default="all",
    )

    parser.add_argument(
        "--output",
        default="mnv_tune_weight_closure",
    )

    args = parser.parse_args()

    modes = (
        list(SAMPLES.keys())
        if args.mode == "all"
        else [args.mode]
    )

    for mode in modes:

        for test_name, test in TESTS.items():

            try:
                compare(
                    mode,
                    test_name,
                    test,
                    Path(args.output) / mode,
                )

            except Exception as e:
                print(
                    f"\n[ERROR] {mode} {test_name}: {e}"
                )


if __name__ == "__main__":
    main()
#!/usr/bin/env python3

"""
Decompose the MC changes from JOINT and MENU1PI relative to V1.

Uses existing ROOT files only.

Studies:
  1. Interaction type (QE, RES, DIS, COH, 2p2h)
  2. Beam-sign charged pion multiplicity
  3. Tune eligibility
  4. True Q2
  5. Raw W
  6. Joint W vs Q2
  7. Joint true Enu vs Q2, including flux peak
  8. Numerical checks that different projections account
     for the same total yield change

This script does NOT use data.

Run:
  python3 studies/decompose_mnv_tune_changes.py --mode all \
      > tune_decomposition.log 2>&1
"""

import argparse
import math
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)


# ============================================================
# CONFIGURATION
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
        "eligibility": "JOINT_eligible",
        "weight": "JOINT_cv_weight",
    },
    "MENU1PI": {
        "tag": "v4p3p0_D2_MENU1PI_LECOH_DiffCOH",
        "eligibility": "MENU1PI_eligible",
        "weight": "MENU1PI_cv_weight",
    },
}

Q2_BINS = [
    (0.0, 0.1),
    (0.1, 0.2),
    (0.2, 0.3),
    (0.3, 0.5),
    (0.5, 0.7),
    (0.7, 1.0),
    (1.0, float("inf")),
]

W_BINS = [
    (0.0, 1.0),
    (1.0, 1.1),
    (1.1, 1.2),
    (1.2, 1.3),
    (1.3, 1.4),
    (1.4, 1.5),
    (1.5, 1.7),
    (1.7, 2.0),
    (2.0, float("inf")),
]

ENERGY_REGIONS = [
    ("Full spectrum", 0.0, float("inf")),
    ("2.0-2.5 GeV", 2.0, 2.5),
    ("Flux peak 2.5-3.8 GeV", 2.5, 3.8),
    ("3.8-5.0 GeV", 3.8, 5.0),
    ("5.0-10.0 GeV", 5.0, 10.0),
]

INT_TYPES = {
    1: "QE",
    2: "RES",
    3: "DIS",
    4: "COH",
    8: "2p2h",
}


# ============================================================
# ROOT UTILITIES
# ============================================================

def file_path(mode, tag):
    c = SAMPLES[mode]

    return (
        f"{c['directory']}/kin_dist_mc{c['playlist']}_"
        f"Validation_{c['sample']}_{tag}_MAD.root"
    )


def open_file(mode, tag):
    path = file_path(mode, tag)

    f = ROOT.TFile.Open(path, "READ")

    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {path}")

    return f


def get_hist(f, name):
    obj = f.Get(name)

    if not obj:
        raise RuntimeError(
            f"Missing '{name}' in {f.GetName()}"
        )

    if hasattr(obj, "GetCVHistoWithStatError"):
        h = obj.GetCVHistoWithStatError()
    else:
        h = obj.Clone()

    if not h:
        raise RuntimeError(f"Cannot clone {name}")

    h.SetDirectory(0)
    return h


def get_pot(f):
    tree = f.Get("Meta")

    if not tree or not tree.GetBranch("POT_Used"):
        raise RuntimeError("Missing Meta/POT_Used")

    pot = sum(float(e.POT_Used) for e in tree)

    if pot <= 0:
        raise RuntimeError("Invalid POT")

    return pot


def compatible(a, b):
    if a.GetNbinsX() != b.GetNbinsX():
        return False

    ax = a.GetXaxis()
    bx = b.GetXaxis()

    for i in range(1, a.GetNbinsX() + 2):
        if not math.isclose(
            ax.GetBinLowEdge(i),
            bx.GetBinLowEdge(i),
            rel_tol=0,
            abs_tol=1e-7,
        ):
            return False

    a2 = a.InheritsFrom("TH2")
    b2 = b.InheritsFrom("TH2")

    if a2 != b2:
        return False

    if a2:
        if a.GetNbinsY() != b.GetNbinsY():
            return False

        ay = a.GetYaxis()
        by = b.GetYaxis()

        for i in range(1, a.GetNbinsY() + 2):
            if not math.isclose(
                ay.GetBinLowEdge(i),
                by.GetBinLowEdge(i),
                rel_tol=0,
                abs_tol=1e-7,
            ):
                return False

    return True


def load_pair(fbase, ftune, name, scale):
    a = get_hist(fbase, name)
    b = get_hist(ftune, name)

    if not compatible(a, b):
        raise RuntimeError(f"Binning mismatch for {name}")

    b.Scale(scale)
    return a, b


# ============================================================
# INTEGRALS
# ============================================================

def integral_1d(h, low=-float("inf"), high=float("inf")):
    total = 0.0

    for i in range(1, h.GetNbinsX() + 1):
        x = h.GetXaxis().GetBinCenter(i)

        if low <= x < high:
            total += h.GetBinContent(i)

    return total


def integral_2d(
    h,
    xlow=-float("inf"),
    xhigh=float("inf"),
    ylow=-float("inf"),
    yhigh=float("inf"),
):
    total = 0.0

    for ix in range(1, h.GetNbinsX() + 1):
        x = h.GetXaxis().GetBinCenter(ix)

        if not xlow <= x < xhigh:
            continue

        for iy in range(1, h.GetNbinsY() + 1):
            y = h.GetYaxis().GetBinCenter(iy)

            if ylow <= y < yhigh:
                total += h.GetBinContent(ix, iy)

    return total


def pct(a, b):
    return 100.0 * a / b if b else float("nan")


def describe_change(name, base, tuned, total_delta):
    delta = tuned - base

    print(
        f"{name:32s}"
        f"{base:14.2f}"
        f"{tuned:14.2f}"
        f"{delta:14.2f}"
        f"{pct(delta, base):12.3f}%"
        f"{pct(delta, total_delta):13.2f}%"
    )

    return delta


def header():
    print(
        f"{'Region':32s}"
        f"{'Baseline':>14s}"
        f"{'Tuned':>14s}"
        f"{'Delta':>14s}"
        f"{'Local change':>13s}"
        f"{'Share of dN':>14s}"
    )


def check_closure(label, components, total):
    summed = sum(components)
    residual = total - summed

    print(
        f"\n  CLOSURE [{label}]"
        f"\n    Reference delta: {total:14.4f}"
        f"\n    Summed delta:    {summed:14.4f}"
        f"\n    Residual:        {residual:14.4f}"
    )

    if abs(residual) > max(1.0, 1e-4 * abs(total)):
        print(
            "    [NOTE] Nonclosure may come from uncovered "
            "bins, flow bins, or non-exhaustive categories."
        )


# ============================================================
# ONE-DIMENSIONAL DECOMPOSITION
# ============================================================

def audit_regions(
    fbase, ftune, hist_name, regions, scale, reference_delta
):
    a, b = load_pair(fbase, ftune, hist_name, scale)

    print("\n" + "=" * 110)
    print(hist_name)
    print("=" * 110)

    header()

    changes = []

    for name, low, high in regions:
        va = integral_1d(a, low, high)
        vb = integral_1d(b, low, high)

        changes.append(
            describe_change(name, va, vb, reference_delta)
        )

    check_closure(hist_name, changes, reference_delta)


def audit_interaction_types(
    fbase, ftune, scale, reference_delta
):
    a, b = load_pair(fbase, ftune, "mc_intType", scale)

    print("\n" + "=" * 110)
    print("GENIE INTERACTION TYPES")
    print("=" * 110)

    header()

    changes = []

    for i in range(1, a.GetNbinsX() + 1):
        x = a.GetXaxis().GetBinCenter(i)

        code = int(round(x))
        name = INT_TYPES.get(code, f"Type {code}")

        va = a.GetBinContent(i)
        vb = b.GetBinContent(i)

        changes.append(
            describe_change(
                f"{name} ({code})",
                va, vb, reference_delta
            )
        )

    check_closure(
        "interaction types",
        changes,
        reference_delta
    )


def audit_pion_counts(
    fbase, ftune, scale, reference_delta
):
    for name in [
        "true_n_beamsign_pion",
        "true_n_piplus",
        "true_n_piminus",
        "true_n_pi0",
        "true_n_pions",
    ]:
        regions = [
            ("0 pions", -0.5, 0.5),
            ("1 pion", 0.5, 1.5),
            ("2+ pions", 1.5, float("inf")),
        ]

        audit_regions(
            fbase, ftune, name, regions,
            scale, reference_delta
        )


def audit_eligibility(
    fbase, ftune, name, scale, reference_delta
):
    a, b = load_pair(fbase, ftune, name, scale)

    print("\n" + "=" * 110)
    print(f"ELIGIBILITY: {name}")
    print("=" * 110)

    header()

    n0a = integral_1d(a, -0.5, 0.5)
    n0b = integral_1d(b, -0.5, 0.5)

    n1a = integral_1d(a, 0.5, 1.5)
    n1b = integral_1d(b, 0.5, 1.5)

    d0 = describe_change(
        "Ineligible", n0a, n0b, reference_delta
    )
    d1 = describe_change(
        "Eligible", n1a, n1b, reference_delta
    )

    total_base = n0a + n1a
    eligible_fraction = n1a / total_base if total_base else 0.0
    mean_weight = n1b / n1a if n1a else float("nan")

    predicted_fractional_change = (
        eligible_fraction * (mean_weight - 1.0)
    )

    print("\nELIGIBILITY EXPLANATION")
    print(f"  Baseline total:             {total_base:.3f}")
    print(f"  Baseline eligible:          {n1a:.3f}")
    print(f"  Baseline ineligible:        {n0a:.3f}")
    print(f"  Eligible fraction:         {100 * eligible_fraction:.4f}%")
    print(f"  Mean eligible correction:  {mean_weight:.6f}")
    print(
        f"  Change from eligibility:   "
        f"{100 * predicted_fractional_change:.4f}%"
    )
    print(
        f"  Observed eligible delta:   {d1:.3f}"
    )
    print(
        f"  Observed ineligible delta: {d0:.3f}"
    )

    check_closure(
        "eligibility",
        [d0, d1],
        reference_delta
    )


# ============================================================
# TWO-DIMENSIONAL DECOMPOSITION
# ============================================================

def audit_w_q2(fbase, ftune, scale, reference_delta):
    a, b = load_pair(
        fbase, ftune, "trueW_trueQ2", scale
    )

    print("\n" + "=" * 110)
    print("JOINT W-Q2 DISTRIBUTION")
    print("X = RAW W, Y = TRUE Q2")
    print("=" * 110)

    header()

    all_changes = []

    for wlow, whigh in W_BINS:

        print(
            f"\nW = {wlow:.2f} to "
            f"{whigh if math.isfinite(whigh) else 'inf'} GeV"
        )

        for qlow, qhigh in Q2_BINS:

            va = integral_2d(
                a, wlow, whigh, qlow, qhigh
            )
            vb = integral_2d(
                b, wlow, whigh, qlow, qhigh
            )

            label = (
                f"Q2 {qlow:.2f}-"
                f"{qhigh if math.isfinite(qhigh) else 'inf'}"
            )

            delta = vb - va
            all_changes.append(delta)

            # Avoid printing regions with no meaningful change.
            if abs(delta) < 0.01:
                continue

            describe_change(
                label, va, vb, reference_delta
            )

    check_closure(
        "W-Q2 grid",
        all_changes,
        reference_delta
    )


def audit_energy_q2(
    fbase, ftune, scale, reference_delta
):
    a, b = load_pair(
        fbase, ftune, "trueEnu_trueQ2", scale
    )

    print("\n" + "=" * 110)
    print("TRUE ENERGY VS TRUE Q2")
    print("X = TRUE ENU, Y = TRUE Q2")
    print("=" * 110)

    header()

    for ename, elow, ehigh in ENERGY_REGIONS:

        base_energy = integral_2d(
            a, elow, ehigh
        )
        tuned_energy = integral_2d(
            b, elow, ehigh
        )

        energy_delta = tuned_energy - base_energy

        print("\n" + "-" * 110)
        print(ename)
        print(
            f"Baseline = {base_energy:.2f}, "
            f"Tuned = {tuned_energy:.2f}, "
            f"Delta = {energy_delta:.2f}, "
            f"Change = {pct(energy_delta, base_energy):.3f}%"
        )
        print("-" * 110)

        changes = []

        for qlow, qhigh in Q2_BINS:
            va = integral_2d(
                a, elow, ehigh, qlow, qhigh
            )
            vb = integral_2d(
                b, elow, ehigh, qlow, qhigh
            )

            label = (
                f"Q2 {qlow:.2f}-"
                f"{qhigh if math.isfinite(qhigh) else 'inf'}"
            )

            changes.append(
                describe_change(
                    label,
                    va, vb,
                    energy_delta
                )
            )

        check_closure(
            f"{ename} / Q2",
            changes,
            energy_delta
        )


# ============================================================
# MAIN
# ============================================================

def run(mode, test_name, test):
    fbase = open_file(mode, BASELINE)
    ftune = open_file(mode, test["tag"])

    try:
        pot_base = get_pot(fbase)
        pot_tune = get_pot(ftune)
        scale = pot_base / pot_tune

        # Use the interaction-type histogram as the
        # common full-population normalization reference.
        hbase, htune = load_pair(
            fbase, ftune, "mc_intType", scale
        )

        base_total = integral_1d(hbase)
        tuned_total = integral_1d(htune)

        total_delta = tuned_total - base_total

        print("\n" + "#" * 110)
        print(f"{mode} | {test_name} vs V1 PROXY")
        print("#" * 110)

        print("Baseline:", BASELINE)
        print("Tune:    ", test["tag"])

        print(f"Baseline POT: {pot_base:.8e}")
        print(f"Tune POT:     {pot_tune:.8e}")
        print(f"POT scale:    {scale:.8f}")

        print(f"\nBaseline selected yield: {base_total:.3f}")
        print(f"Tuned selected yield:    {tuned_total:.3f}")
        print(f"Absolute difference:     {total_delta:+.3f}")
        print(
            f"Relative difference:     "
            f"{pct(total_delta, base_total):+.4f}%"
        )

        audit_interaction_types(
            fbase, ftune, scale, total_delta
        )

        audit_eligibility(
            fbase, ftune,
            test["eligibility"],
            scale, total_delta
        )

        audit_pion_counts(
            fbase, ftune, scale, total_delta
        )

        audit_regions(
            fbase, ftune, "mc_w_raw",
            [
                (
                    f"W {low:.2f}-"
                    f"{high if math.isfinite(high) else 'inf'}",
                    low, high
                )
                for low, high in W_BINS
            ],
            scale, total_delta
        )

        audit_regions(
            fbase, ftune, "Q2",
            [
                (
                    f"Q2 {low:.2f}-"
                    f"{high if math.isfinite(high) else 'inf'}",
                    low, high
                )
                for low, high in Q2_BINS
            ],
            scale, total_delta
        )

        audit_w_q2(
            fbase, ftune, scale, total_delta
        )

        audit_energy_q2(
            fbase, ftune, scale, total_delta
        )

    finally:
        fbase.Close()
        ftune.Close()


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--mode",
        choices=["FHC", "RHC", "all"],
        default="all",
    )

    args = parser.parse_args()

    modes = (
        ["FHC", "RHC"]
        if args.mode == "all"
        else [args.mode]
    )

    for mode in modes:
        for name, test in TESTS.items():
            run(mode, name, test)


if __name__ == "__main__":
    main()
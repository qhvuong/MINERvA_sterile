#!/usr/bin/env python3

import argparse
import csv
import math
from pathlib import Path
from collections import OrderedDict

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.AddDirectory(False)


# ============================================================
# INPUT
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

TAGS = [
    "v1_D2_JOINT_LECOH_DiffCOH",
    "v1_D2_LECOH_DiffCOH",
    "v4_D2_MENU1PI_LECOH_DiffCOH",
    "v4_D2_MENU1PI_MECOH_DiffCOH",
    "v4p3p0_D2_MENU1PI_LECOH_DiffCOH",
    "v4p3p0_D2_MENU1PI_MECOH_DiffCOH",
]

REFERENCE = "v1_D2_LECOH_DiffCOH"

COMPARISONS = [
    (
        "JOINT vs V1 without JOINT",
        "v1_D2_JOINT_LECOH_DiffCOH",
        "v1_D2_LECOH_DiffCOH",
    ),
    (
        "V4p3p0 MENU1PI vs V1 without JOINT",
        "v4p3p0_D2_MENU1PI_LECOH_DiffCOH",
        "v1_D2_LECOH_DiffCOH",
    ),
    (
        "V4p3p0 MENU1PI vs V1 JOINT",
        "v4p3p0_D2_MENU1PI_LECOH_DiffCOH",
        "v1_D2_JOINT_LECOH_DiffCOH",
    ),
    (
        "V4p3p0 vs V4",
        "v4p3p0_D2_MENU1PI_LECOH_DiffCOH",
        "v4_D2_MENU1PI_LECOH_DiffCOH",
    ),
    (
        "V4p3p0 MECOH vs LECOH",
        "v4p3p0_D2_MENU1PI_MECOH_DiffCOH",
        "v4p3p0_D2_MENU1PI_LECOH_DiffCOH",
    ),
    (
        "V4 MECOH vs LECOH",
        "v4_D2_MENU1PI_MECOH_DiffCOH",
        "v4_D2_MENU1PI_LECOH_DiffCOH",
    ),
]


# ============================================================
# HISTOGRAMS
# ============================================================

HISTS_1D = OrderedDict([
    ("Flux", "True Neutrino Energy"),
    ("EN4", "Reconstructed Energy"),
    ("E_avail", "Available Energy"),

    ("Q2", "True Q2"),

    ("mc_w_raw", "Raw True W"),
    ("trueW", "True W"),
    ("trueWExp", "True W Exp"),

    ("true_n_piplus", "Pi Plus Multiplicity"),
    ("true_n_piminus", "Pi Minus Multiplicity"),
    ("true_n_pi0", "Pi0 Multiplicity"),
    ("true_n_pions", "Total Pion Multiplicity"),
    ("true_n_beamsign_pion", "Beam Sign Pion Multiplicity"),

    ("mc_intType", "GENIE Interaction Type"),

    ("MENU1PI_eligible", "MENU1PI Eligibility"),
    ("JOINT_eligible", "JOINT Eligibility"),

    ("MENU1PI_cv_weight", "MENU1PI CV Weight"),
    ("JOINT_cv_weight", "JOINT CV Weight"),
])

HISTS_2D = OrderedDict([
    ("trueEnu_trueQ2", "True Enu vs True Q2"),
    ("trueW_trueQ2", "True W vs True Q2"),

    ("MENU1PIWeight_trueQ2", "MENU1PI Weight vs True Q2"),
    ("JOINTWeight_trueQ2", "JOINT Weight vs True Q2"),

    ("MENU1PIWeight_trueEnu", "MENU1PI Weight vs True Enu"),
    ("JOINTWeight_trueEnu", "JOINT Weight vs True Enu"),
])


# ============================================================
# PHYSICS REGIONS
# ============================================================

# Ranges are interpreted using bin centers.
# Boundary bins are therefore approximations if the requested
# boundaries are not histogram-bin edges.

REGIONS_1D = {
    "Flux": [
        ("1.5-2.5 GeV", 1.5, 2.5),
        ("2.5-3.8 GeV", 2.5, 3.8),
        ("3.8-5.0 GeV", 3.8, 5.0),
        ("5.0-10.0 GeV", 5.0, 10.0),
    ],
    "EN4": [
        ("0.8-2.5 GeV", 0.8, 2.5),
        ("2.5-3.8 GeV", 2.5, 3.8),
        ("3.8-5.5 GeV", 3.8, 5.5),
        ("5.5-10.0 GeV", 5.5, 10.0),
    ],
    "Q2": [
        ("Q2 < 0.1", 0.0, 0.1),
        ("0.1 <= Q2 < 0.2", 0.1, 0.2),
        ("0.2 <= Q2 < 0.5", 0.2, 0.5),
        ("0.5 <= Q2 < 1.0", 0.5, 1.0),
        ("Q2 >= 1.0", 1.0, float("inf")),
    ],
    "mc_w_raw": [
        ("W < 1.4", 0.0, 1.4),
        ("1.4 <= W < 2.0", 1.4, 2.0),
        ("W >= 2.0", 2.0, float("inf")),
    ],
    "true_n_beamsign_pion": [
        ("0 pion", -0.5, 0.5),
        ("exactly 1 pion", 0.5, 1.5),
        ("2+ pions", 1.5, float("inf")),
    ],
}

WEIGHT_WINDOWS = [
    ("Suppression: w < 0.8", -float("inf"), 0.8),
    ("Moderate suppression: 0.8-0.95", 0.8, 0.95),
    ("Near unity: 0.95-1.05", 0.95, 1.05),
    ("Moderate enhancement: 1.05-1.2", 1.05, 1.2),
    ("Enhancement: w >= 1.2", 1.2, float("inf")),
]


# ============================================================
# UTILITIES
# ============================================================

def ratio(a, b):
    return a / b if b != 0 else float("nan")


def percent_change(a, b):
    return 100.0 * (ratio(a, b) - 1.0)


def fmt(x):
    if not math.isfinite(x):
        return "nan"
    return f"{x:.6f}"


def input_path(mode, tag):
    c = SAMPLES[mode]
    return (
        Path(c["directory"])
        / (
            f"kin_dist_mc{c['playlist']}_"
            f"Validation_{c['sample']}_{tag}_MAD.root"
        )
    )


def read_pot(f):
    tree = f.Get("Meta")

    if not tree or not tree.GetBranch("POT_Used"):
        return None

    pot = sum(float(e.POT_Used) for e in tree)

    return pot if pot > 0 else None


def get_hist(f, name, dimension=1, scale=1.0):
    obj = f.Get(name)

    if not obj:
        return None

    if hasattr(obj, "GetCVHistoWithStatError"):
        h = obj.GetCVHistoWithStatError()
    else:
        h = obj.Clone()

    if not h:
        return None

    actual_dimension = 2 if h.InheritsFrom("TH2") else 1

    if actual_dimension != dimension:
        return None

    h.SetDirectory(0)

    if scale != 1.0:
        h.Scale(scale)

    return h


def same_binning(a, b):
    if a.GetNbinsX() != b.GetNbinsX():
        return False

    ax = a.GetXaxis()
    bx = b.GetXaxis()

    for i in range(1, a.GetNbinsX() + 2):
        if abs(
            ax.GetBinLowEdge(i) - bx.GetBinLowEdge(i)
        ) > 1e-7:
            return False

    if a.InheritsFrom("TH2") != b.InheritsFrom("TH2"):
        return False

    if a.InheritsFrom("TH2"):
        if a.GetNbinsY() != b.GetNbinsY():
            return False

        ay = a.GetYaxis()
        by = b.GetYaxis()

        for i in range(1, a.GetNbinsY() + 2):
            if abs(
                ay.GetBinLowEdge(i) - by.GetBinLowEdge(i)
            ) > 1e-7:
                return False

    return True


def integral_1d(h, low=-float("inf"), high=float("inf")):
    result = 0.0

    for i in range(1, h.GetNbinsX() + 1):
        x = h.GetXaxis().GetBinCenter(i)

        if low <= x < high:
            result += h.GetBinContent(i)

    return result


def total_1d(h):
    return h.Integral(1, h.GetNbinsX())


def weighted_mean_x(h):
    sumw = 0.0
    sumwx = 0.0

    for i in range(1, h.GetNbinsX() + 1):
        w = h.GetBinContent(i)
        x = h.GetXaxis().GetBinCenter(i)

        sumw += w
        sumwx += w * x

    return sumwx / sumw if sumw else float("nan")


def write_csv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


# ============================================================
# 1D DIAGNOSTICS
# ============================================================

def audit_1d(mode, files, scales, output):

    rows = []
    region_rows = []
    missing = []

    for name, label in HISTS_1D.items():

        hs = OrderedDict()

        for tag, f in files.items():
            h = get_hist(f, name, 1, scales[tag])

            if h is None:
                missing.append((tag, name))
            else:
                hs[tag] = h

        if REFERENCE not in hs:
            print(f"\n[MISSING REFERENCE HIST] {name}")
            continue

        href = hs[REFERENCE]

        print("\n" + "=" * 110)
        print(f"[{mode}] {label}: {name}")
        print("=" * 110)

        print(
            f"{'Tag':60s}"
            f"{'Yield':>15s}"
            f"{'Ratio':>12s}"
            f"{'Change %':>13s}"
            f"{'Mean X':>13s}"
        )

        for tag, h in hs.items():

            if not same_binning(h, href):
                print("[BINNING MISMATCH]", tag, name)
                continue

            n = total_1d(h)
            nref = total_1d(href)

            r = ratio(n, nref)
            pct = percent_change(n, nref)
            mean = weighted_mean_x(h)

            print(
                f"{tag:60s}"
                f"{n:15.2f}"
                f"{r:12.5f}"
                f"{pct:13.3f}"
                f"{mean:13.5f}"
            )

            rows.append({
                "mode": mode,
                "histogram": name,
                "tag": tag,
                "yield": n,
                "reference_yield": nref,
                "ratio": r,
                "change_percent": pct,
                "mean_x": mean,
            })

        # ----------------------------------------------------
        # Bin-by-bin comparison
        # ----------------------------------------------------

        print("\nBIN-BY-BIN DIFFERENCES")

        for title, numerator, denominator in COMPARISONS:

            if numerator not in hs or denominator not in hs:
                continue

            hn = hs[numerator]
            hd = hs[denominator]

            if not same_binning(hn, hd):
                continue

            print(f"\n  {title}")
            print(
                f"{'Bin':>15s}"
                f"{'Old':>15s}"
                f"{'New':>15s}"
                f"{'Delta':>15s}"
                f"{'Change %':>14s}"
            )

            for i in range(1, hn.GetNbinsX() + 1):

                low = hn.GetXaxis().GetBinLowEdge(i)
                high = hn.GetXaxis().GetBinUpEdge(i)

                old = hd.GetBinContent(i)
                new = hn.GetBinContent(i)

                print(
                    f"{low:6.2f}-{high:<7.2f}"
                    f"{old:15.2f}"
                    f"{new:15.2f}"
                    f"{new - old:15.2f}"
                    f"{percent_change(new, old):14.2f}"
                )

        # ----------------------------------------------------
        # Targeted regions
        # ----------------------------------------------------

        regions = list(REGIONS_1D.get(name, []))

        if name in ("MENU1PI_cv_weight", "JOINT_cv_weight"):
            regions.extend(WEIGHT_WINDOWS)

        if name in ("MENU1PI_eligible", "JOINT_eligible"):
            regions.extend([
                ("Ineligible", -0.5, 0.5),
                ("Eligible", 0.5, 1.5),
            ])

        if name in (
            "true_n_piplus",
            "true_n_piminus",
            "true_n_pi0",
            "true_n_pions",
        ):
            regions.extend([
                ("0 pions", -0.5, 0.5),
                ("1 pion", 0.5, 1.5),
                ("2+ pions", 1.5, float("inf")),
            ])

        if name == "mc_intType":
            regions.extend([
                ("QE", 0.5, 1.5),
                ("RES", 1.5, 2.5),
                ("DIS", 2.5, 3.5),
                ("COH", 3.5, 4.5),
                ("2p2h", 7.5, 8.5),
            ])

        for region_name, low, high in regions:

            print(f"\nREGION: {region_name}")

            for title, numerator, denominator in COMPARISONS:

                if numerator not in hs or denominator not in hs:
                    continue

                hn = hs[numerator]
                hd = hs[denominator]

                if not same_binning(hn, hd):
                    continue

                new = integral_1d(hn, low, high)
                old = integral_1d(hd, low, high)

                delta = new - old
                pct = percent_change(new, old)

                print(
                    f"  {title:45s}"
                    f" old={old:13.2f}"
                    f" new={new:13.2f}"
                    f" delta={delta:13.2f}"
                    f" change={pct:9.3f}%"
                )

                region_rows.append({
                    "mode": mode,
                    "histogram": name,
                    "region": region_name,
                    "comparison": title,
                    "numerator_tag": numerator,
                    "denominator_tag": denominator,
                    "old_yield": old,
                    "new_yield": new,
                    "delta": delta,
                    "change_percent": pct,
                })

    write_csv(
        output / "one_dimensional_summary.csv",
        rows,
        [
            "mode", "histogram", "tag", "yield",
            "reference_yield", "ratio",
            "change_percent", "mean_x",
        ],
    )

    write_csv(
        output / "one_dimensional_regions.csv",
        region_rows,
        [
            "mode", "histogram", "region",
            "comparison", "numerator_tag",
            "denominator_tag", "old_yield",
            "new_yield", "delta", "change_percent",
        ],
    )

    return missing


# ============================================================
# 2D DIAGNOSTICS
# ============================================================

def integral_2d(h, xlow, xhigh, ylow, yhigh):

    result = 0.0

    for ix in range(1, h.GetNbinsX() + 1):

        x = h.GetXaxis().GetBinCenter(ix)

        if not xlow <= x < xhigh:
            continue

        for iy in range(1, h.GetNbinsY() + 1):

            y = h.GetYaxis().GetBinCenter(iy)

            if ylow <= y < yhigh:
                result += h.GetBinContent(ix, iy)

    return result


def get_2d_regions(name):

    if name == "trueEnu_trueQ2":
        # X = true Enu, Y = true Q2
        return [
            ("peak, Q2 < 0.1", 2.5, 3.8, 0.0, 0.1),
            ("peak, Q2 < 0.2", 2.5, 3.8, 0.0, 0.2),
            ("peak, Q2 < 0.5", 2.5, 3.8, 0.0, 0.5),
            ("peak, Q2 >= 0.5", 2.5, 3.8, 0.5, float("inf")),
            ("2-5 GeV, Q2 < 0.2", 2.0, 5.0, 0.0, 0.2),
        ]

    if name == "trueW_trueQ2":
        # X = W, Y = Q2
        return [
            ("W < 1.4, Q2 < 0.1", 0.0, 1.4, 0.0, 0.1),
            ("W < 1.4, Q2 < 0.2", 0.0, 1.4, 0.0, 0.2),
            ("W < 1.4, Q2 < 0.5", 0.0, 1.4, 0.0, 0.5),
            ("W < 1.4, Q2 >= 0.5", 0.0, 1.4, 0.5, float("inf")),
            ("W >= 1.4, Q2 < 0.2", 1.4, float("inf"), 0.0, 0.2),
            ("W >= 1.4, Q2 >= 0.2", 1.4, float("inf"), 0.2, float("inf")),
        ]

    if name.endswith("_trueEnu"):
        # X = true Enu, Y = evaluated diagnostic weight
        return [
            ("peak, weight < 0.8", 2.5, 3.8, -float("inf"), 0.8),
            ("peak, weight 0.8-1.0", 2.5, 3.8, 0.8, 1.0),
            ("peak, weight >= 1.0", 2.5, 3.8, 1.0, float("inf")),
            ("2-5 GeV, weight < 0.8", 2.0, 5.0, -float("inf"), 0.8),
            ("2-5 GeV, weight >= 1.0", 2.0, 5.0, 1.0, float("inf")),
        ]

    if name.endswith("_trueQ2"):
        # X = Q2, Y = evaluated diagnostic weight
        return [
            ("Q2 < 0.1, weight < 0.8", 0.0, 0.1, -float("inf"), 0.8),
            ("Q2 < 0.2, weight < 0.8", 0.0, 0.2, -float("inf"), 0.8),
            ("Q2 < 0.5, weight < 0.8", 0.0, 0.5, -float("inf"), 0.8),
            ("Q2 < 0.5, weight >= 1.0", 0.0, 0.5, 1.0, float("inf")),
        ]

    return []


def audit_2d(mode, files, scales, output):

    rows = []
    missing = []

    for name, title in HISTS_2D.items():

        hs = OrderedDict()

        for tag, f in files.items():

            h = get_hist(f, name, 2, scales[tag])

            if h is None:
                missing.append((tag, name))
            else:
                hs[tag] = h

        if REFERENCE not in hs:
            print(f"\n[MISSING 2D REFERENCE] {name}")
            continue

        print("\n" + "=" * 110)
        print(f"[{mode}] 2D: {title}")
        print("=" * 110)

        regions = get_2d_regions(name)

        for region_name, xlow, xhigh, ylow, yhigh in regions:

            print(f"\nREGION: {region_name}")

            for title, numerator, denominator in COMPARISONS:

                if numerator not in hs or denominator not in hs:
                    continue

                hn = hs[numerator]
                hd = hs[denominator]

                if not same_binning(hn, hd):
                    continue

                new = integral_2d(
                    hn, xlow, xhigh, ylow, yhigh
                )

                old = integral_2d(
                    hd, xlow, xhigh, ylow, yhigh
                )

                delta = new - old
                pct = percent_change(new, old)

                print(
                    f"  {title:45s}"
                    f" old={old:13.2f}"
                    f" new={new:13.2f}"
                    f" delta={delta:13.2f}"
                    f" change={pct:9.3f}%"
                )

                rows.append({
                    "mode": mode,
                    "histogram": name,
                    "region": region_name,
                    "comparison": title,
                    "numerator_tag": numerator,
                    "denominator_tag": denominator,
                    "old_yield": old,
                    "new_yield": new,
                    "delta": delta,
                    "change_percent": pct,
                })

    write_csv(
        output / "two_dimensional_regions.csv",
        rows,
        [
            "mode", "histogram", "region", "comparison",
            "numerator_tag", "denominator_tag",
            "old_yield", "new_yield", "delta",
            "change_percent",
        ],
    )

    return missing


# ============================================================
# POT AND FILE VALIDATION
# ============================================================

def open_samples(mode):

    files = OrderedDict()

    print(f"\n{'=' * 110}")
    print(f"OPENING {mode}")
    print("=" * 110)

    for tag in TAGS:

        path = input_path(mode, tag)

        if not path.exists():
            print("[MISSING FILE]", path)
            continue

        f = ROOT.TFile.Open(str(path), "READ")

        if not f or f.IsZombie():
            print("[INVALID ROOT FILE]", path)
            continue

        files[tag] = f
        print("[FOUND]", path.name)

    return files


def get_scales(files, apply_pot_scale=True):

    pots = {
        tag: read_pot(f)
        for tag, f in files.items()
    }

    scales = {
        tag: 1.0
        for tag in files
    }

    reference_pot = pots.get(REFERENCE)

    if apply_pot_scale and reference_pot:

        for tag in files:

            if pots[tag]:
                scales[tag] = reference_pot / pots[tag]
            else:
                print(
                    f"[WARNING] Missing POT for {tag}; "
                    "using scale 1"
                )

    elif apply_pot_scale:
        print(
            "[WARNING] Reference POT unavailable; "
            "using scale 1 for all files"
        )

    return pots, scales


# ============================================================
# MAIN
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--mode",
        choices=["FHC", "RHC", "all"],
        default="all",
    )

    parser.add_argument(
        "--output",
        default="mnv_tune_numeric_audit",
    )

    parser.add_argument(
        "--no-pot-scale",
        action="store_true",
        help="Disable POT scaling",
    )

    args = parser.parse_args()

    modes = (
        ["FHC", "RHC"]
        if args.mode == "all"
        else [args.mode]
    )

    for mode in modes:

        files = open_samples(mode)

        if REFERENCE not in files:
            print(f"[ERROR] Reference missing for {mode}")

            for f in files.values():
                f.Close()

            continue

        try:

            output = Path(args.output) / mode
            output.mkdir(parents=True, exist_ok=True)

            pots, scales = get_scales(
                files,
                apply_pot_scale=not args.no_pot_scale,
            )

            print("\nPOT SUMMARY")

            for tag in files:
                print(
                    f"{tag:60s}"
                    f"POT={pots[tag]} "
                    f"scale={scales[tag]:.8f}"
                )

            missing_1d = audit_1d(
                mode, files, scales, output
            )

            missing_2d = audit_2d(
                mode, files, scales, output
            )

            missing = missing_1d + missing_2d

            if missing:

                print("\nMISSING HISTOGRAMS")

                for tag, name in missing:
                    print(tag, name)

            print(f"\n[DONE] {mode}: {output}")

        finally:

            for f in files.values():
                f.Close()


if __name__ == "__main__":
    main()
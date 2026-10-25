#!/usr/bin/env python3

"""
Final validation of MINERvA JOINT / MENU1PI tune comparisons.

No selection rerun required.

Checks:
  1. POT consistency across all six MC configurations.
  2. Normalization differences between configurations.
  3. Whether JOINT/MENU1PI changes are confined to their
     expected eligibility, interaction-type, pion and W regions.
  4. Full histogram-integral closure, including underflow
     and overflow.
  5. Whether the stand-alone Q2 histogram's missing yield
     explains its previous nonclosure.
  6. Whether changes between V4 and V4p3p0 are negligible
     or nonzero.
  7. Source-code audit: find weight definitions, diagnostic
     calls, and CV-weight construction.

The source-code audit is STATIC. It prints evidence to inspect;
it does not assert that matching symbol names prove identical
runtime behavior.

Example:

python3 studies/final_mnv_tune_audit.py \
    --mode all \
    --source-root . \
    --output final_tune_audit \
    > final_tune_audit.log 2>&1
"""

import argparse
import csv
import math
import re
from pathlib import Path

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

TAGS = {
    "V1": "v1_D2_LECOH_DiffCOH",

    "JOINT": "v1_D2_JOINT_LECOH_DiffCOH",

    "V4_MENU": "v4_D2_MENU1PI_LECOH_DiffCOH",

    "V4_MENU_MECOH": "v4_D2_MENU1PI_MECOH_DiffCOH",

    "V4P3_MENU": "v4p3p0_D2_MENU1PI_LECOH_DiffCOH",

    "V4P3_MENU_MECOH":
        "v4p3p0_D2_MENU1PI_MECOH_DiffCOH",
}

COMPARISONS = [
    ("V1", "JOINT", "JOINT only"),
    ("V1", "V4_MENU", "V1 -> V4 MENU1PI"),
    ("V4_MENU", "V4P3_MENU", "V4 -> V4p3p0"),
    ("V4_MENU", "V4_MENU_MECOH", "V4 LECOH -> MECOH"),
    (
        "V4P3_MENU",
        "V4P3_MENU_MECOH",
        "V4p3p0 LECOH -> MECOH",
    ),
    ("V1", "V4P3_MENU", "V1 -> V4p3p0 MENU1PI"),
]

CLOSURE_COMPARISONS = [
    ("V1", "JOINT", "JOINT"),
    ("V1", "V4P3_MENU", "MENU1PI"),
]

# All of these histograms already appeared in our earlier
# numeric diagnostic outputs.
HISTOGRAMS = [
    "mc_intType",
    "JOINT_eligible",
    "MENU1PI_eligible",
    "true_n_beamsign_pion",
    "mc_w_raw",
    "Q2",
    "trueW_trueQ2",
    "trueEnu_trueQ2",
    "Flux",
]

SOURCE_PATTERNS = {
    "CV weight construction": [
        r"\bGetStandardWeight\s*\(",
        r"\bGetMyLowQ2PiWeight\s*\(",
        r"\bGetLowQ2PiWeight\s*\(",
    ],
    "Tune configuration and routing": [
        r"LowQ2PiWeightChannel",
        r"MENU1PI",
        r"\bJOINT\b",
        r"\bgetWeight\s*\(",
    ],
    "Eligibility and diagnostic implementation": [
        r"IsCCResDiagnostic",
        r"IsCCNucleonPionDiagnostic",
        r"CountBeamSignChargedPions",
        r"JOINT_eligible",
        r"MENU1PI_eligible",
        r"JOINT_cv_weight",
        r"MENU1PI_cv_weight",
        r"JOINTWeight_trueEnu",
        r"MENU1PIWeight_trueEnu",
    ],
}

SOURCE_SUFFIXES = {
    ".py", ".cpp", ".cxx", ".cc", ".c",
    ".h", ".hpp", ".hxx",
}

SKIP_DIRS = {
    ".git", "__pycache__", "build", "node_modules",
    ".venv", "venv", ".tox",
}

EPS = 1e-7


# ============================================================
# HELPERS
# ============================================================

def section(title):
    print()
    print("=" * 104)
    print(title)
    print("=" * 104)


def input_path(mode, tag):
    c = SAMPLES[mode]
    return (
        Path(c["directory"])
        / f"kin_dist_mc{c['playlist']}_"
          f"Validation_{c['sample']}_{tag}_MAD.root"
    )


def open_file(mode, key):
    path = input_path(mode, TAGS[key])

    if not path.is_file():
        raise FileNotFoundError(path)

    f = ROOT.TFile.Open(str(path), "READ")

    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {path}")

    return f


def get_pot(f):
    tree = f.Get("Meta")

    if not tree or not tree.GetBranch("POT_Used"):
        raise RuntimeError(f"Missing POT_Used in {f.GetName()}")

    pot = sum(float(e.POT_Used) for e in tree)

    if pot <= 0:
        raise RuntimeError(f"Invalid POT in {f.GetName()}")

    return pot


def get_hist(f, name):
    obj = f.Get(name)

    if not obj:
        raise RuntimeError(
            f"Missing histogram {name} in {f.GetName()}"
        )

    if hasattr(obj, "GetCVHistoWithStatError"):
        h = obj.GetCVHistoWithStatError()
    else:
        h = obj.Clone()

    if not h:
        raise RuntimeError(f"Could not clone {name}")

    h.SetDirectory(0)
    return h


def dimension(h):
    return 2 if h.InheritsFrom("TH2") else 1


def axes_match(a, b):
    if dimension(a) != dimension(b):
        return False

    axes = [(a.GetXaxis(), b.GetXaxis())]

    if dimension(a) == 2:
        axes.append((a.GetYaxis(), b.GetYaxis()))

    for aa, bb in axes:
        if aa.GetNbins() != bb.GetNbins():
            return False

        for i in range(1, aa.GetNbins() + 2):
            if not math.isclose(
                aa.GetBinLowEdge(i),
                bb.GetBinLowEdge(i),
                rel_tol=0,
                abs_tol=EPS,
            ):
                return False

    return True


def get_pair(files, pots, left, right, name):
    a = get_hist(files[left], name)
    b = get_hist(files[right], name)

    if not axes_match(a, b):
        raise RuntimeError(
            f"Binning differs: {name}, {left} -> {right}"
        )

    b.Scale(pots[left] / pots[right])

    return a, b


def visible_integral(h):
    if dimension(h) == 1:
        return sum(
            h.GetBinContent(ix)
            for ix in range(1, h.GetNbinsX() + 1)
        )

    return sum(
        h.GetBinContent(ix, iy)
        for ix in range(1, h.GetNbinsX() + 1)
        for iy in range(1, h.GetNbinsY() + 1)
    )


def full_integral(h):
    if dimension(h) == 1:
        return sum(
            h.GetBinContent(ix)
            for ix in range(0, h.GetNbinsX() + 2)
        )

    return sum(
        h.GetBinContent(ix, iy)
        for ix in range(0, h.GetNbinsX() + 2)
        for iy in range(0, h.GetNbinsY() + 2)
    )


def flow_integral(h):
    return full_integral(h) - visible_integral(h)


def integral_1d_bins(h, predicate):
    total = 0.0

    for i in range(1, h.GetNbinsX() + 1):
        low = h.GetXaxis().GetBinLowEdge(i)
        high = h.GetXaxis().GetBinUpEdge(i)
        center = h.GetXaxis().GetBinCenter(i)

        if predicate(low, high, center):
            total += h.GetBinContent(i)

    return total


def integral_1d_range(h, low, high):
    """
    Use only bins wholly contained in [low, high].

    This avoids silently including neighboring physical
    regions through a bin-center approximation.
    """
    return integral_1d_bins(
        h,
        lambda lo, hi, center:
            lo >= low - EPS and hi <= high + EPS,
    )


def abs_close(a, b, tolerance=1e-5):
    return abs(a - b) <= tolerance


def percent(numerator, denominator):
    if denominator == 0:
        return float("nan")

    return 100.0 * numerator / denominator


def print_delta(label, base, tuned):
    delta = tuned - base

    print(
        f"{label:35s}"
        f"{base:15.3f}"
        f"{tuned:15.3f}"
        f"{delta:15.3f}"
        f"{percent(delta, base):12.4f}%"
    )

    return delta


# ============================================================
# 1. POT AND CONFIGURATION COMPARISONS
# ============================================================

def compare_configurations(files, pots):
    section("1. POT AND CONFIGURATION COMPARISONS")

    for key in TAGS:
        print(
            f"{key:19s}"
            f"POT = {pots[key]:.10e}"
            f"  {TAGS[key]}"
        )

    section("SELECTED YIELD DIFFERENCES BY CONFIGURATION")

    print(
        f"{'Comparison':35s}"
        f"{'Baseline':>15s}"
        f"{'Tuned':>15s}"
        f"{'Delta':>15s}"
        f"{'Change':>13s}"
    )

    rows = []

    for left, right, label in COMPARISONS:
        a, b = get_pair(
            files, pots, left, right, "mc_intType"
        )

        na = visible_integral(a)
        nb = visible_integral(b)

        delta = print_delta(label, na, nb)

        rows.append({
            "comparison": label,
            "left": left,
            "right": right,
            "baseline": na,
            "tuned": nb,
            "delta": delta,
            "percent_change": percent(delta, na),
        })

    return rows


# ============================================================
# 2. HISTOGRAM FLOW AND POPULATION CONSISTENCY
# ============================================================

def audit_flows(files, pots, left, right):
    section(
        f"2. HISTOGRAM FLOW AUDIT: {left} -> {right}"
    )

    a_ref, b_ref = get_pair(
        files, pots, left, right, "mc_intType"
    )

    ref_base = visible_integral(a_ref)
    ref_tuned = visible_integral(b_ref)
    ref_delta = ref_tuned - ref_base

    print(f"Reference baseline: {ref_base:.6f}")
    print(f"Reference tuned:    {ref_tuned:.6f}")
    print(f"Reference delta:    {ref_delta:+.6f}")

    print()
    print(
        f"{'Histogram':25s}"
        f"{'Base visible':>15s}"
        f"{'Base flow':>14s}"
        f"{'Tune flow':>14s}"
        f"{'Visible dN':>15s}"
        f"{'Full dN':>15s}"
        f"{'Full dN - ref':>16s}"
    )

    for name in HISTOGRAMS:
        try:
            a, b = get_pair(
                files, pots, left, right, name
            )

            va = visible_integral(a)
            vb = visible_integral(b)

            fa = flow_integral(a)
            fb = flow_integral(b)

            visible_delta = vb - va
            full_delta = full_integral(b) - full_integral(a)

            print(
                f"{name:25s}"
                f"{va:15.3f}"
                f"{fa:14.3f}"
                f"{fb:14.3f}"
                f"{visible_delta:15.3f}"
                f"{full_delta:15.3f}"
                f"{full_delta - ref_delta:16.3f}"
            )

            if abs(full_delta - ref_delta) > 1.0:
                print(
                    f"    [NOTE] {name} does not reproduce "
                    "the reference total even with flow bins."
                )

        except Exception as exc:
            print(f"  [ERROR] {name}: {exc}")


# ============================================================
# 3. APPLICABILITY CHECKS
# ============================================================

def check_expected_zero(label, delta, reference):
    tolerance = max(1e-4, abs(reference) * 1e-7)

    status = (
        "PASS"
        if abs(delta) <= tolerance
        else "CHECK"
    )

    print(
        f"  [{status}] {label:44s}"
        f"delta = {delta:+.6f}"
    )


def audit_applicability(files, pots, left, right, tune):
    section(
        f"3. APPLICABILITY: {tune}, {left} -> {right}"
    )

    a_ref, b_ref = get_pair(
        files, pots, left, right, "mc_intType"
    )

    total_delta = (
        visible_integral(b_ref) - visible_integral(a_ref)
    )

    eligibility_name = f"{tune}_eligible"

    # --------------------------------------------------------
    # Eligibility
    # --------------------------------------------------------

    a, b = get_pair(
        files, pots, left, right, eligibility_name
    )

    ia = integral_1d_range(a, -0.5, 0.5)
    ib = integral_1d_range(b, -0.5, 0.5)

    ea = integral_1d_range(a, 0.5, 1.5)
    eb = integral_1d_range(b, 0.5, 1.5)

    d_ineligible = ib - ia
    d_eligible = eb - ea

    print("\nEligibility:")

    print(f"  Baseline eligible yield: {ea:.6f}")
    print(f"  Tuned eligible yield:    {eb:.6f}")
    print(f"  Eligible fraction:       {percent(ea, ea + ia):.5f}%")

    mean = eb / ea if ea else float("nan")

    print(f"  Mean eligible correction: {mean:.7f}")
    print(f"  Eligible delta:           {d_eligible:+.6f}")
    print(f"  Ineligible delta:         {d_ineligible:+.6f}")

    check_expected_zero(
        "Ineligible population unchanged",
        d_ineligible,
        total_delta,
    )

    print(
        "  Eligible closure residual: "
        f"{total_delta - d_eligible:+.6f}"
    )

    # --------------------------------------------------------
    # GENIE interaction type
    # --------------------------------------------------------

    print("\nInteraction types:")

    expected_types = (
        {2} if tune == "JOINT" else {2, 3}
    )

    for i in range(1, a_ref.GetNbinsX() + 1):
        code = round(a_ref.GetXaxis().GetBinCenter(i))

        va = a_ref.GetBinContent(i)
        vb = b_ref.GetBinContent(i)

        delta = vb - va

        print(
            f"  Type {code:2d}: "
            f"base={va:13.3f} "
            f"tuned={vb:13.3f} "
            f"delta={delta:+13.3f}"
        )

        if code not in expected_types:
            check_expected_zero(
                f"Type {code} unchanged",
                delta,
                total_delta,
            )

    # --------------------------------------------------------
    # Pion multiplicity
    # --------------------------------------------------------

    a, b = get_pair(
        files, pots, left, right,
        "true_n_beamsign_pion",
    )

    print("\nBeam-sign charged pion multiplicity:")

    for label, low, high in [
        ("0 pions", -0.5, 0.5),
        ("1 pion", 0.5, 1.5),
        ("2+ pions", 1.5, float("inf")),
    ]:
        va = integral_1d_range(a, low, high)
        vb = integral_1d_range(b, low, high)

        d = vb - va

        print(
            f"  {label:12s}"
            f"base={va:13.3f} "
            f"tuned={vb:13.3f} "
            f"delta={d:+13.3f}"
        )

        if tune == "MENU1PI" and label != "1 pion":
            check_expected_zero(
                f"{label} unchanged",
                d,
                total_delta,
            )

    # --------------------------------------------------------
    # Raw W
    # --------------------------------------------------------

    a, b = get_pair(
        files, pots, left, right, "mc_w_raw"
    )

    print("\nRaw W:")

    for label, low, high in [
        ("W < 1.0", -float("inf"), 1.0),
        ("1.0 <= W < 1.4", 1.0, 1.4),
        ("W >= 1.4", 1.4, float("inf")),
    ]:
        va = integral_1d_range(a, low, high)
        vb = integral_1d_range(b, low, high)

        d = vb - va

        print(
            f"  {label:18s}"
            f"base={va:13.3f} "
            f"tuned={vb:13.3f} "
            f"delta={d:+13.3f}"
        )

        if tune == "MENU1PI" and label == "W >= 1.4":
            check_expected_zero(
                "W >= 1.4 unchanged",
                d,
                total_delta,
            )

    # --------------------------------------------------------
    # JOINT Q2 boundary
    # --------------------------------------------------------

    if tune == "JOINT":
        a, b = get_pair(
            files, pots, left, right, "Q2"
        )

        va = integral_1d_range(
            a, 0.7, float("inf")
        )
        vb = integral_1d_range(
            b, 0.7, float("inf")
        )

        check_expected_zero(
            "Q2 >= 0.7 unchanged",
            vb - va,
            total_delta,
        )


# ============================================================
# 4. COMPARE CHANGES BETWEEN VERSIONS
# ============================================================

def audit_version_differences(files, pots):
    section("4. V4 / V4p3p0 / MECOH CROSS-CHECKS")

    checks = [
        (
            "V4_MENU",
            "V4P3_MENU",
            "V4 -> V4p3p0",
        ),
        (
            "V4_MENU",
            "V4_MENU_MECOH",
            "V4 LECOH -> MECOH",
        ),
        (
            "V4P3_MENU",
            "V4P3_MENU_MECOH",
            "V4p3p0 LECOH -> MECOH",
        ),
    ]

    for left, right, title in checks:
        print("\n" + "-" * 104)
        print(title)
        print("-" * 104)

        for name in [
            "mc_intType",
            "mc_w_raw",
            "Q2",
            "Flux",
            "trueW_trueQ2",
            "trueEnu_trueQ2",
        ]:
            try:
                a, b = get_pair(
                    files, pots, left, right, name
                )

                va = visible_integral(a)
                vb = visible_integral(b)

                d = vb - va

                print(
                    f"{name:22s}"
                    f" delta={d:+14.4f}"
                    f"  relative={percent(d, va):+10.5f}%"
                )

            except Exception as exc:
                print(f"[ERROR] {name}: {exc}")


# ============================================================
# 5. CHECK THE WEIGHT HISTOGRAM'S UNITY BIN
# ============================================================

def audit_unity_bin(files, pots):
    section("5. DIAGNOSTIC WEIGHT BINNING CHECK")

    for tune in ("JOINT", "MENU1PI"):
        name = f"{tune}Weight_trueEnu"

        try:
            h = get_hist(files["V1"], name)

            if dimension(h) != 2:
                print(f"[ERROR] {name} is not TH2")
                continue

            yaxis = h.GetYaxis()
            iy = yaxis.FindBin(1.0)

            low = yaxis.GetBinLowEdge(iy)
            high = yaxis.GetBinUpEdge(iy)
            center = yaxis.GetBinCenter(iy)

            population = sum(
                h.GetBinContent(ix, iy)
                for ix in range(1, h.GetNbinsX() + 1)
            )

            print(f"\n{tune}")
            print(f"  Weight histogram: {name}")
            print(f"  Unity bin:        {iy}")
            print(f"  Low edge:         {low:.8f}")
            print(f"  High edge:        {high:.8f}")
            print(f"  Bin center:       {center:.8f}")
            print(f"  Yield in bin:     {population:.4f}")

            print(
                "  Artificial shift if every event in this "
                "bin were assigned the center:"
            )
            print(
                f"    {population * (center - 1.0):+.4f}"
            )

            print(
                "  [NOTE] This only diagnoses histogram-level "
                "reconstruction bias. It does not change "
                "event-level weights."
            )

        except Exception as exc:
            print(f"[ERROR] {name}: {exc}")


# ============================================================
# 6. SOURCE-CODE AUDIT
# ============================================================

def source_files(root):
    for path in root.rglob("*"):
        if not path.is_file():
            continue

        if path.suffix.lower() not in SOURCE_SUFFIXES:
            continue

        if any(part in SKIP_DIRS for part in path.parts):
            continue

        yield path


def source_audit(root, output_dir):
    section("6. STATIC SOURCE-CODE AUDIT")

    root = Path(root).resolve()

    if not root.is_dir():
        print(f"[ERROR] Source directory not found: {root}")
        return

    print("Source root:", root)
    print()
    print(
        "This audit finds candidate call sites and surrounding "
        "source lines. Review the full control flow manually."
    )

    matches = []
    files_examined = 0

    compiled = {
        group: [re.compile(p) for p in patterns]
        for group, patterns in SOURCE_PATTERNS.items()
    }

    for path in source_files(root):
        try:
            if path.stat().st_size > 3_000_000:
                continue

            lines = path.read_text(
                encoding="utf-8", errors="replace"
            ).splitlines()

        except OSError:
            continue

        files_examined += 1

        for number, line in enumerate(lines, 1):
            for group, patterns in compiled.items():
                if not any(p.search(line) for p in patterns):
                    continue

                matches.append({
                    "group": group,
                    "path": str(path.relative_to(root)),
                    "line": number,
                    "content": line.strip(),
                })

    print(f"\nFiles examined: {files_examined}")
    print(f"Matching source lines: {len(matches)}")

    audit_txt = output_dir / "source_weight_audit.txt"

    # Print context around every match to a separate file,
    # rather than flooding the main log.
    with audit_txt.open("w") as out:
        for group in SOURCE_PATTERNS:
            out.write("\n" + "=" * 100 + "\n")
            out.write(group + "\n")
            out.write("=" * 100 + "\n")

            group_matches = [
                m for m in matches
                if m["group"] == group
            ]

            print(
                f"\n{group}: {len(group_matches)} matching lines"
            )

            paths = sorted(set(
                m["path"] for m in group_matches
            ))

            for path in paths[:25]:
                print(f"  {path}")

            if len(paths) > 25:
                print(
                    f"  ... and {len(paths) - 25} more files"
                )

            for m in group_matches:
                out.write(
                    f"\n{m['path']}:{m['line']} "
                    f"[{m['group']}]\n"
                )

                path = root / m["path"]

                try:
                    lines = path.read_text(
                        encoding="utf-8",
                        errors="replace",
                    ).splitlines()
                except OSError:
                    continue

                first = max(1, m["line"] - 5)
                last = min(len(lines), m["line"] + 5)

                for line_number in range(first, last + 1):
                    out.write(
                        f"{line_number:5d} | "
                        f"{lines[line_number - 1]}\n"
                    )

    print(f"\nSaved source audit: {audit_txt}")

    print("\nIMPORTANT MANUAL QUESTIONS:")
    print(
        "  1. Is GetMyLowQ2PiWeight called exactly once "
        "in the production CV weight?"
    )
    print(
        "  2. Are other tune-weight functions also applied "
        "along that same call chain?"
    )
    print(
        "  3. Does the diagnostic use the same channel, Q2, "
        "target and eligibility conditions?"
    )
    print(
        "  4. Are V1 and V4p3p0 otherwise identical in "
        "the actual configured weight chain?"
    )
    print(
        "  5. Do the evaluated JOINT and MENU1PI curves "
        "agree with their independent definitions?"
    )

    print(
        "\n[NOTE] A search match is not proof that a "
        "function is executed. Conditional branches, "
        "inheritance, and runtime configuration must "
        "still be inspected."
    )


# ============================================================
# MAIN
# ============================================================

def run_mode(mode):
    section(f"FINAL MINERVA TUNE AUDIT: {mode}")

    files = {}
    pots = {}

    try:
        for key in TAGS:
            files[key] = open_file(mode, key)
            pots[key] = get_pot(files[key])

        rows = compare_configurations(files, pots)

        for left, right, tune in CLOSURE_COMPARISONS:
            audit_flows(files, pots, left, right)

            audit_applicability(
                files, pots, left, right, tune
            )

        audit_version_differences(files, pots)

        audit_unity_bin(files, pots)

        return rows

    finally:
        for f in files.values():
            f.Close()


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--mode",
        choices=["FHC", "RHC", "all"],
        default="all",
    )

    parser.add_argument(
        "--source-root",
        default=None,
        help=(
            "Directory containing analysis and weight code. "
            "Can be omitted to run ROOT checks only."
        ),
    )

    parser.add_argument(
        "--output",
        default="final_tune_audit",
    )

    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    modes = (
        ["FHC", "RHC"]
        if args.mode == "all"
        else [args.mode]
    )

    rows = []

    for mode in modes:
        try:
            for row in run_mode(mode):
                rows.append({
                    "mode": mode,
                    **row,
                })
        except Exception as exc:
            print(f"\n[ERROR] {mode}: {exc}")
            raise

    if rows:
        csv_path = output_dir / "configuration_comparisons.csv"

        with csv_path.open("w", newline="") as out:
            writer = csv.DictWriter(
                out,
                fieldnames=list(rows[0].keys()),
            )

            writer.writeheader()
            writer.writerows(rows)

        print(f"\nSaved comparison CSV: {csv_path}")

    if args.source_root:
        source_audit(
            args.source_root,
            output_dir,
        )
    else:
        section("SOURCE AUDIT NOT RUN")
        print(
            "Pass --source-root /path/to/source "
            "to examine the weight implementation."
        )

    section("AUDIT COMPLETE")
    print(
        "Interpret PASS/CHECK messages together with "
        "the histogram flow report and source audit."
    )


if __name__ == "__main__":
    main()
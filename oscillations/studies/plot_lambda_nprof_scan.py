#!/usr/bin/env python3

import os
import re
import glob
import csv
import math
import argparse
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# CASE_PATTERNS = {
#     "excludeRatio_maskRatio": "prodNueel exclude ratio, mask ratio",
#     "includeRatio": "prodNueel include ratio",
#     "profileOnlyRatio_maskNonRatio": "prodNueel only ratio",
#     "noRatio": "prodNueel_noRatio",
# }
# CASE_PATTERNS = {
#     "excludeRatio_maskRatio":
#         r"$\nu+e$ + CC ratios + CC$\nu_\mu$, ratios excluded from profiling",

#     "noRatio_p8_profileAll":
#         r"$\nu+e$ + CC$\nu_e$ + CC$\nu_\mu$, all profiled",

#     "profileOnlyRatio":
#         r"Profile only CC ratios",

#     "_p8_profileAll":
#         r"$\nu+e$ + CC ratios + CC$\nu_\mu$, all profiled",
# }
MASKED_RATIO_CASE = r"$\nu+e$ + CC$\nu_\mu$, ratios omitted from profiling and $\chi^2$"

CASE_PATTERNS = {
    "excludeRatio_maskRatio": MASKED_RATIO_CASE,
    "excludeRatio":
        r"$\nu+e$ + CC ratios + CC$\nu_\mu$, ratios excluded from profiling",

    "directCC_profileAll":
        r"$\nu+e$ + CC$\nu_e$ + CC$\nu_\mu$, all profiled",

    "ratioOnly":
        r"Profile only CC ratios",

    "ratio_profileAll":
        r"$\nu+e$ + CC ratios + CC$\nu_\mu$, all profiled",
}
CASE_SAFE_NAMES = {
    MASKED_RATIO_CASE: "excludeRatio_maskRatio",
    r"$\nu+e$ + CC ratios + CC$\nu_\mu$, ratios excluded from profiling":
        "excludeRatio",

    r"$\nu+e$ + CC ratios + CC$\nu_\mu$, all profiled":
        "includeAll_ratioConfig",

    r"Profile only CC ratios":
        "profileOnlyRatio",

    r"$\nu+e$ + CC$\nu_e$ + CC$\nu_\mu$, all profiled":
        "directCC_includeAll",
}
CASE_ORDER = [
    MASKED_RATIO_CASE,
    r"$\nu+e$ + CC$\nu_e$ + CC$\nu_\mu$, all profiled",
    r"$\nu+e$ + CC ratios + CC$\nu_\mu$, ratios excluded from profiling",
    r"$\nu+e$ + CC ratios + CC$\nu_\mu$, all profiled",
    r"Profile only CC ratios",
]

def safe_case_name(case):
    if case not in CASE_SAFE_NAMES:
        raise RuntimeError(
            "No safe filename mapping defined for case: {}".format(case)
        )

    return CASE_SAFE_NAMES[case]


def finite(x):
    try:
        return math.isfinite(float(x))
    except Exception:
        return False


def read_csv(path):
    rows = []
    with open(path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            out = {}
            for k, v in row.items():
                try:
                    out[k] = float(v)
                except Exception:
                    out[k] = v
            rows.append(out)

    rows.sort(key=lambda r: r["lambda"])
    return rows


def infer_case_and_nprof(path):
    base = os.path.basename(path)
    stem = os.path.splitext(base)[0]

    m = re.search(r"Nprof(\d+)", stem)
    if not m:
        raise RuntimeError("Could not infer Nprof from filename: {}".format(base))
    nprof = int(m.group(1))

    # Current scan filenames include tune tags between p8 and profileAll.
    # Match their suffix explicitly before the older filename conventions.
    mode = re.search(r"_(profileAll|profileOnlyRatio|excludeRatio_maskRatio)_Nprof\d+(?:_|$)", stem)
    if mode:
        name = mode.group(1)
        if name == "excludeRatio_maskRatio":
            return MASKED_RATIO_CASE, nprof
        if name == "profileOnlyRatio":
            return CASE_PATTERNS["ratioOnly"], nprof
        key = "directCC_profileAll" if "noRatio" in stem else "ratio_profileAll"
        return CASE_PATTERNS[key], nprof

    case = None
    for key, label in CASE_PATTERNS.items():
        if key in stem:
            case = label
            break

    if case is None:
        raise RuntimeError("Could not infer case from filename: {}".format(base))

    return case, nprof


# def get_lambda1_row(rows):
#     return min(rows, key=lambda r: abs(float(r["lambda"]) - 1.0))
def get_lambda1_row(rows):
    matches = [
        r for r in rows
        if abs(float(r["lambda"]) - 1.0) < 1e-10
    ]

    if len(matches) != 1:
        raise RuntimeError(
            "Expected exactly one lambda=1 row, found {}".format(
                len(matches)
            )
        )

    return matches[0]


LCURVE_YLABELS = {
    "norm2": r"Unweighted flux pull $\|a\|^2$",
    "penalty": r"Flux penalty $\lambda a^T a$",
}


def lcurve_points(rows, y_mode):
    key = "norm_a_null" if y_mode == "norm2" else "penalty_null"
    points = []
    for r in rows:
        if all(finite(r.get(k)) for k in ("resid_null", key, "lambda")):
            y = float(r[key])
            points.append((float(r["resid_null"]),
                           y ** 2 if y_mode == "norm2" else y,
                           float(r["lambda"])))
    return points


def draw_lcurve(rows, y_mode, label, linewidth):
    points = lcurve_points(rows, y_mode)
    if not points:
        raise ValueError("No finite L-curve points for {}".format(label))
    x, y, lam = zip(*points)
    line, = plt.plot(x, y, marker="o", linewidth=linewidth, markersize=4, label=label)
    for xi, yi, li in points:
        if math.isclose(li, 1., rel_tol=0., abs_tol=1e-10):
            plt.scatter([xi], [yi], marker="D", s=30, color=line.get_color(),
                        edgecolors="black", linewidths=.5, zorder=3)
        if any(math.isclose(li, v, rel_tol=1e-9)
               for v in (0.01, 0.1, 1., 2., 5., 10., 100.)):
            plt.annotate(r"$\lambda={:g}$".format(li), (xi, yi), fontsize=6,
                         xytext=(3, 3), textcoords="offset points")


def finish_lcurve(title, out, y_mode):
    plt.xlabel(r"Residual $\chi^2$")
    plt.ylabel(LCURVE_YLABELS[y_mode])
    plt.title(title)
    plt.grid(True, alpha=0.35)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    print("Wrote", out)


def plot_lcurves_by_case(grouped, outdir, y_mode="norm2"):
    """One plot per configuration, with a curve for each Nprof."""
    suffix = "_norm2" if y_mode == "norm2" else ""
    for case in CASE_ORDER:
        if case not in grouped:
            continue
        plt.figure(figsize=(8, 6))
        for nprof, rows in sorted(grouped[case].items()):
            draw_lcurve(rows, y_mode, "Nprof={}".format(nprof), 1.5)
        out = os.path.join(outdir, "lcurve_by_nprof_{}{}.png".format(
            safe_case_name(case), suffix))
        finish_lcurve("Null L-curve: {}".format(case), out, y_mode)


def plot_lcurves_by_nprof(grouped, outdir, y_mode="norm2"):
    """One plot per Nprof, with a curve for each configuration."""
    suffix = "_norm2" if y_mode == "norm2" else ""
    for nprof in sorted({n for by_nprof in grouped.values() for n in by_nprof}):
        plt.figure(figsize=(8, 6))
        for case in CASE_ORDER:
            if nprof in grouped.get(case, {}):
                draw_lcurve(grouped[case][nprof], y_mode, case, 1.8)
        out = os.path.join(outdir, "lcurve_compare_cases_Nprof{}{}.png".format(nprof, suffix))
        finish_lcurve("Null L-curve comparison (Nprof = {})".format(nprof), out, y_mode)


def plot_lambda1_vs_nprof(grouped, outdir):
    """
    For each case, show lambda=1 diagnostics vs Nprof.
    """
    metrics = [
        ("chi2_null", r"Total profiled null $\chi^2$", "lambda1_chi2_null_vs_nprof.png"),
        ("resid_null", r"Residual null $\chi^2$", "lambda1_resid_null_vs_nprof.png"),
        ("penalty_null", r"Flux penalty", "lambda1_penalty_null_vs_nprof.png"),
        ("norm_a_null", r"$||a||$", "lambda1_norm_a_vs_nprof.png"),
        ("max_abs_a_null", r"$\max |a_i|$", "lambda1_max_abs_a_vs_nprof.png"),
    ]

    for key, ylabel, outname in metrics:
        plt.figure(figsize=(8, 5.5))

        # for case, by_nprof in grouped.items():
        for case in CASE_ORDER:
            if case not in grouped:
                continue
            by_nprof = grouped[case]
            xs = []
            ys = []

            for nprof in sorted(by_nprof):
                r = get_lambda1_row(by_nprof[nprof])
                if finite(r.get(key)):
                    xs.append(nprof)
                    ys.append(float(r[key]))

            plt.plot(xs, ys, marker="o", linewidth=2, label=case)

        plt.xlabel("Number of profiled flux universes")
        plt.ylabel(ylabel)
        plt.title("{} at lambda = 1".format(ylabel))
        plt.grid(True, alpha=0.35)
        plt.legend(fontsize=8)
        plt.tight_layout()

        out = os.path.join(outdir, outname)
        plt.savefig(out, dpi=200)
        plt.close()
        print("Wrote", out)


def plot_vs_lambda_grid(grouped, outdir):
    """
    For each case, plot residual/penalty/norm as function of lambda,
    with different curves for Nprof.
    """
    metrics = [
        ("chi2_null", r"Total profiled null $\chi^2$", "chi2_null"),
        ("resid_null", r"Residual null $\chi^2$", "resid_null"),
        ("penalty_null", r"Flux penalty", "penalty_null"),
        ("norm_a_null", r"$||a||$", "norm_a_null"),
        ("max_abs_a_null", r"$\max |a_i|$", "max_abs_a_null"),
    ]

    # for case, by_nprof in grouped.items():
    for case in CASE_ORDER:
        if case not in grouped:
            continue
        by_nprof = grouped[case]
        # safe = case.replace(" ", "_").replace(",", "").replace("/", "-")
        safe = safe_case_name(case)

        for key, ylabel, tag in metrics:
            plt.figure(figsize=(8, 5.5))

            for nprof in sorted(by_nprof):
                rows = by_nprof[nprof]
                x = []
                y = []

                for r in rows:
                    if finite(r.get("lambda")) and finite(r.get(key)):
                        x.append(float(r["lambda"]))
                        y.append(float(r[key]))

                plt.plot(x, y, marker="o", linewidth=1.5, markersize=4, label="Nprof={}".format(nprof))

            plt.xscale("log")
            plt.xlabel(r"$\lambda$")
            plt.ylabel(ylabel)
            plt.title("{}: {}".format(case, ylabel))
            plt.grid(True, alpha=0.35)
            plt.legend(fontsize=8)
            plt.tight_layout()

            out = os.path.join(outdir, "{}_vs_lambda_{}.png".format(tag, safe))
            plt.savefig(out, dpi=200)
            plt.close()
            print("Wrote", out)


def write_lambda1_summary(grouped, outpath):
    fields = ["case", "nprof", "lambda", "chi2_null", "resid_null",
              "penalty_null", "norm_a_null", "max_abs_a_null"]
    with open(outpath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for case in CASE_ORDER:
            if case not in grouped:
                continue
            for nprof in sorted(grouped[case]):
                r = get_lambda1_row(grouped[case][nprof])
                writer.writerow(dict(case=case, nprof=nprof,
                                     **{key: r[key] for key in fields[2:]}))
    print("Wrote", outpath)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--glob",
        # default="/exp/minerva/data/users/qvuong/surfaces/lambda_scans/nprof/*Nprof*_fixedFlux.csv",
        # default="/exp/minerva/data/users/qvuong/surfaces/csvs/p8_lambda_scans/*Nprof*.csv",
        default="/exp/minerva/data/users/qvuong/surfaces/lambda_scans/p8_fixedFluxCV/nprof/*Nprof*.csv",
    )
    parser.add_argument(
        "--outdir",
        # default="/exp/minerva/data/users/qvuong/surfaces/lambda_scans/plots_nprof_lcurve_fixedFlux",
        default="/exp/minerva/data/users/qvuong/surfaces/plots/p8_lambda_scans",
    )
    parser.add_argument(
        "--lcurve-y", choices=("norm2", "penalty", "both"), default="norm2",
        help="L-curve y-axis: unweighted norm squared (default), lambda-weighted penalty, or both.",
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    files = sorted(glob.glob(args.glob))
    print("Found {} files".format(len(files)))

    grouped = defaultdict(dict)

    for path in files:
        case, nprof = infer_case_and_nprof(path)
        rows = read_csv(path)
        if nprof in grouped[case]:
            raise RuntimeError("Duplicate input for {} at Nprof={}: {}".format(case, nprof, path))
        grouped[case][nprof] = rows
        print("Loaded case='{}' Nprof={} rows={} path={}".format(case, nprof, len(rows), path))

    if len(grouped) == 0:
        raise RuntimeError("No datasets found")

    modes = ("norm2", "penalty") if args.lcurve_y == "both" else (args.lcurve_y,)
    for mode in modes:
        plot_lcurves_by_case(grouped, args.outdir, mode)
        plot_lcurves_by_nprof(grouped, args.outdir, mode)
    plot_lambda1_vs_nprof(grouped, args.outdir)
    plot_vs_lambda_grid(grouped, args.outdir)
    write_lambda1_summary(grouped, os.path.join(args.outdir, "lambda1_nprof_summary.csv"))

    print("\nDone. Output:", args.outdir)


if __name__ == "__main__":
    main()
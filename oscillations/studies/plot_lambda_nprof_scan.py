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
CASE_PATTERNS = {
    "excludeRatio_maskRatio":
        r"$\nu+e$ + CC ratios + CC$\nu_\mu$, ratios excluded from profiling",

    "prodNueel_noRatio_p8_profileAll":
        r"$\nu+e$ + CC$\nu_e$ + CC$\nu_\mu$, all profiled",

    "profileOnlyRatio":
        r"Profile only CC ratios",

    "prodNueel_p8_profileAll":
        r"$\nu+e$ + CC ratios + CC$\nu_\mu$, all profiled",
}
CASE_SAFE_NAMES = {
    r"$\nu+e$ + CC ratios + CC$\nu_\mu$, ratios excluded from profiling":
        "excludeRatio",

    r"$\nu+e$ + CC ratios + CC$\nu_\mu$, all profiled":
        "includeAll_ratioConfig",

    r"Profile only CC ratios":
        "profileOnlyRatio",

    r"$\nu+e$ + CC$\nu_e$ + CC$\nu_\mu$, all profiled":
        "directCC_includeAll",
}


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


def plot_lcurves_by_case(grouped, outdir):
    """
    One L-curve plot per case. Curves are different Nprof values.
    """
    for case, by_nprof in grouped.items():
        plt.figure(figsize=(8, 6))

        for nprof in sorted(by_nprof):
            rows = by_nprof[nprof]
            x = []
            y = []
            lam = []

            for r in rows:
                if finite(r.get("resid_null")) and finite(r.get("penalty_null")):
                    x.append(float(r["resid_null"]))
                    y.append(float(r["penalty_null"]))
                    lam.append(float(r["lambda"]))

            plt.plot(x, y, marker="o", linewidth=1.5, markersize=4, label="Nprof={}".format(nprof))

            for xi, yi, li in zip(x, y, lam):
                if li in [0.1, 1, 2, 5, 10, 100]:
                    plt.text(xi, yi, "{:g}".format(li), fontsize=6)

        plt.xlabel(r"Residual $\chi^2$")
        plt.ylabel(r"Flux penalty $\lambda a^T a$")
        plt.title("Null L-curve: {}".format(case))
        plt.grid(True, alpha=0.35)
        plt.legend(fontsize=8)
        plt.tight_layout()

        # safe = case.replace(" ", "_").replace(",", "").replace("/", "-")
        safe = safe_case_name(case)
        out = os.path.join(outdir, "lcurve_by_nprof_{}.png".format(safe))
        plt.savefig(out, dpi=200)
        plt.close()
        print("Wrote", out)


def plot_lcurves_by_nprof(grouped, outdir):
    """
    One L-curve plot per Nprof.
    Curves are the different analysis cases.

    grouped structure:
        grouped[case][nprof] = rows
    """
    all_nprof = sorted(
        set(
            nprof
            for case, by_nprof in grouped.items()
            for nprof in by_nprof.keys()
        )
    )

    for nprof in all_nprof:
        plt.figure(figsize=(8, 6))

        plotted_any = False

        for case, by_nprof in grouped.items():
            if nprof not in by_nprof:
                continue

            rows = by_nprof[nprof]

            x = []
            y = []
            lam = []

            for r in rows:
                if finite(r.get("resid_null")) and finite(r.get("penalty_null")) and finite(r.get("lambda")):
                    x.append(float(r["resid_null"]))
                    y.append(float(r["penalty_null"]))
                    lam.append(float(r["lambda"]))

            if len(x) == 0:
                continue

            plotted_any = True
            plt.plot(x, y, marker="o", linewidth=1.8, markersize=4, label=case)

            # Optional: label a few lambda points
            for xi, yi, li in zip(x, y, lam):
                if li in [0.1, 1, 2, 5, 10, 100]:
                    plt.text(xi, yi, "{:g}".format(li), fontsize=6)

        if not plotted_any:
            plt.close()
            continue

        plt.xlabel(r"Residual $\chi^2$")
        plt.ylabel(r"Flux penalty $\lambda a^T a$")
        plt.title("Null L-curve comparison across cases (Nprof = {})".format(nprof))
        plt.grid(True, alpha=0.35)
        plt.legend(fontsize=8)
        plt.tight_layout()

        out = os.path.join(outdir, "lcurve_compare_cases_Nprof{}.png".format(nprof))
        plt.savefig(out, dpi=200)
        plt.close()
        print("Wrote", out)


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

        for case, by_nprof in grouped.items():
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

    for case, by_nprof in grouped.items():
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
    with open(outpath, "w") as f:
        f.write("case,nprof,lambda,chi2_null,resid_null,penalty_null,norm_a_null,max_abs_a_null\n")

        for case, by_nprof in grouped.items():
            for nprof in sorted(by_nprof):
                r = get_lambda1_row(by_nprof[nprof])
                f.write(
                    "{},{},{:.6g},{:.8g},{:.8g},{:.8g},{:.8g},{:.8g}\n".format(
                        case,
                        nprof,
                        r["lambda"],
                        r["chi2_null"],
                        r["resid_null"],
                        r["penalty_null"],
                        r["norm_a_null"],
                        r["max_abs_a_null"],
                    )
                )

    print("Wrote", outpath)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--glob",
        # default="/exp/minerva/data/users/qvuong/surfaces/lambda_scans/nprof/*Nprof*_fixedFlux.csv",
        default="/exp/minerva/data/users/qvuong/surfaces/csvs/p8_lambda_scans/*Nprof*.csv",
    )
    parser.add_argument(
        "--outdir",
        # default="/exp/minerva/data/users/qvuong/surfaces/lambda_scans/plots_nprof_lcurve_fixedFlux",
        default="/exp/minerva/data/users/qvuong/surfaces/plots/p8_lambda_scans",
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    files = sorted(glob.glob(args.glob))
    print("Found {} files".format(len(files)))

    grouped = defaultdict(dict)

    for path in files:
        case, nprof = infer_case_and_nprof(path)
        rows = read_csv(path)
        grouped[case][nprof] = rows
        print("Loaded case='{}' Nprof={} rows={} path={}".format(case, nprof, len(rows), path))

    if len(grouped) == 0:
        raise RuntimeError("No datasets found")

    plot_lcurves_by_case(grouped, args.outdir)
    plot_lcurves_by_nprof(grouped, args.outdir)
    plot_lambda1_vs_nprof(grouped, args.outdir)
    plot_vs_lambda_grid(grouped, args.outdir)
    write_lambda1_summary(grouped, os.path.join(args.outdir, "lambda1_nprof_summary.csv"))

    print("\nDone. Output:", args.outdir)


if __name__ == "__main__":
    main()
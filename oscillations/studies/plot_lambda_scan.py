#!/usr/bin/env python3

import os
import csv
import glob
import argparse
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_FILES = [
    "/exp/minerva/data/users/qvuong/surfaces/csvs/data_lambda_scan_prodNueel_excludeRatio_maskRatio_nullOnly.csv",
    "/exp/minerva/data/users/qvuong/surfaces/csvs/data_lambda_scan_prodNueel_includeAll_maskNone_nullOnly.csv",
    "/exp/minerva/data/users/qvuong/surfaces/csvs/data_lambda_scan_prodNueel_profileOnlyRatio_maskNonRatio_nullOnly.csv",
    "/exp/minerva/data/users/qvuong/surfaces/csvs/data_lambda_scan_prodNueel_noRatio_includeAll_maskNone_nullOnly.csv",
]


LABEL_MAP = {
    "data_lambda_scan_prodNueel_excludeRatio_maskRatio_nullOnly": "prodNueel exclude ratio",
    "data_lambda_scan_prodNueel_includeAll_maskNone_nullOnly": "prodNueel include ratio",
    "data_lambda_scan_prodNueel_profileOnlyRatio_maskNonRatio_nullOnly": "prodNueel only ratio",
    "data_lambda_scan_prodNueel_noRatio_includeAll_maskNone_nullOnly": "prodNueel_noRatio",
    "data_lcurve_prodNueel_profileOnlyRatio_maskNonRatio_fixedFlux": "prodNueel only ratio",
    "data_lcurve_prodNueel_excludeRatio_maskRatio_fixedFlux": "prodNueel exclude ratio, mask ratio",
    "data_lcurve_prodNueel_includeRatio_fixedFlux": "prodNueel include ratio",
    "data_lcurve_prodNueel_noRatio_fixedFlux": "prodNueel_noRatio",
}


def nice_label(path):
    base = os.path.basename(path)
    base = os.path.splitext(base)[0]
    return LABEL_MAP.get(base, base.replace("lambda_scan_", "").replace("_", " "))


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


def get_col(rows, key):
    return [r[key] for r in rows]

def finite(x):
    try:
        return math.isfinite(float(x))
    except Exception:
        return False


def get_xy(rows, xkey, ykey):
    x = []
    y = []

    for r in rows:
        if xkey not in r or ykey not in r:
            continue
        if not finite(r[xkey]) or not finite(r[ykey]):
            continue

        x.append(float(r[xkey]))
        y.append(float(r[ykey]))

    return x, y

def plot_vs_lambda(datasets, ykey, ylabel, outpath, logx=True):
    plt.figure(figsize=(8, 5.5))

    for label, rows in datasets:
        x, y = get_xy(rows, "lambda", ykey)

        if len(x) == 0:
            print("WARNING: no valid values for", ykey, "in", label)
            continue

        plt.plot(x, y, marker="o", linewidth=2, label=label)

    if logx:
        plt.xscale("log")

    plt.xlabel(r"$\lambda$")
    plt.ylabel(ylabel)
    plt.grid(True, alpha=0.35)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def plot_l_curve(datasets, outpath, which="bf"):
    """
    L-curve: penalty vs residual chi2.
    which = "bf" or "null"

    x-axis: residual chi2
    y-axis: flux penalty
    """
    if which == "bf":
        xkey = "resid_bf"
        ykey = "penalty_bf"
        title = r"Data BF L-curve"
    elif which == "null":
        xkey = "resid_null"
        ykey = "penalty_null"
        title = r"Data null L-curve"
    else:
        raise ValueError("which must be bf or null")

    plt.figure(figsize=(7, 6))

    for label, rows in datasets:
        x = []
        y = []
        lam = []

        for r in rows:
            if not finite(r.get(xkey, None)) or not finite(r.get(ykey, None)) or not finite(r.get("lambda", None)):
                continue

            x.append(float(r[xkey]))
            y.append(float(r[ykey]))
            lam.append(float(r["lambda"]))

        if len(x) == 0:
            print("WARNING: no valid L-curve values for", label, which)
            continue

        plt.plot(x, y, marker="o", linewidth=2, label=label)

        for xi, yi, li in zip(x, y, lam):
            plt.text(xi, yi, "{:g}".format(li), fontsize=7)

    plt.xlabel(r"Residual $\chi^2$")
    plt.ylabel("Flux penalty")
    plt.title(title + r" ($\lambda$ labels shown)")
    plt.grid(True, alpha=0.35)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

def plot_lambda_vs_residual(datasets, outpath, which="null"):
    """
    Swapped-axis residual scan.

    x-axis: residual chi2
    y-axis: lambda
    """
    if which == "null":
        xkey = "resid_null"
        title = r"Data null residual $\chi^2$ scan"
    elif which == "bf":
        xkey = "resid_bf"
        title = r"Data BF residual $\chi^2$ scan"
    else:
        raise ValueError("which must be null or bf")

    plt.figure(figsize=(7, 6))

    for label, rows in datasets:
        x = []
        y = []

        for r in rows:
            if not finite(r.get(xkey, None)) or not finite(r.get("lambda", None)):
                continue

            x.append(float(r[xkey]))
            y.append(float(r["lambda"]))

        if len(x) == 0:
            print("WARNING: no valid residual scan values for", label, which)
            continue

        plt.plot(x, y, marker="o", linewidth=2, label=label)

    plt.yscale("log")
    plt.xlabel(r"Residual $\chi^2$")
    plt.ylabel(r"$\lambda$")
    plt.title(title)
    plt.grid(True, alpha=0.35)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

def plot_chi2_decomposition(datasets, outdir):
    """
    For each dataset, make one plot showing BF residual, BF penalty, BF total vs lambda.
    """
    for label, rows in datasets:
        safe = (
            label.replace(" ", "_")
            .replace(",", "")
            .replace("/", "-")
            .replace("(", "")
            .replace(")", "")
        )

        x = get_col(rows, "lambda")
        resid = get_col(rows, "resid_bf")
        pen = get_col(rows, "penalty_bf")
        total = get_col(rows, "chi2_bf")

        plt.figure(figsize=(8, 5.5))
        plt.plot(x, resid, marker="o", linewidth=2, label="BF residual")
        plt.plot(x, pen, marker="o", linewidth=2, label="BF penalty")
        plt.plot(x, total, marker="o", linewidth=2, label="BF total")

        plt.xscale("log")
        plt.xlabel(r"$\lambda$")
        plt.ylabel(r"$\chi^2$")
        plt.title(label)
        plt.grid(True, alpha=0.35)
        plt.legend(fontsize=9)
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "chi2_decomp_{}.png".format(safe)), dpi=200)
        plt.close()


def write_summary_table(datasets, outpath):
    """
    Write a compact summary at lambda=1 if available.
    """
    with open(outpath, "w") as f:
        f.write(
            "label,lambda,chi2_null,chi2_bf,delta_chi2,"
            "resid_null,penalty_null,resid_bf,penalty_bf,"
            "norm_a_null,norm_a_bf,dm2,ue4,umu4,utau4\n"
        )

        for label, rows in datasets:
            # Prefer lambda=1.0; otherwise use closest to 1.
            best = min(rows, key=lambda r: abs(r["lambda"] - 1.0))

            f.write(
                "{},{:.6g},{:.6g},{:.6g},{:.6g},"
                "{:.6g},{:.6g},{:.6g},{:.6g},"
                "{:.6g},{:.6g},{:.6g},{:.6g},{:.6g},{:.6g}\n".format(
                    label,
                    best["lambda"],
                    best["chi2_null"],
                    best["chi2_bf"],
                    best["delta_chi2"],
                    best["resid_null"],
                    best["penalty_null"],
                    best["resid_bf"],
                    best["penalty_bf"],
                    best["norm_a_null"],
                    best["norm_a_bf"],
                    best["dm2"],
                    best["ue4"],
                    best["umu4"],
                    best["utau4"],
                )
            )


def main():
    parser = argparse.ArgumentParser(
        description="Plot lambda-scan / L-curve diagnostics from CSV files."
    )
    parser.add_argument(
        "--files",
        nargs="+",
        default=DEFAULT_FILES,
        help="CSV files to plot. Default: known lambda_scan files in csvs/",
    )
    parser.add_argument(
        "--glob",
        default=None,
        help='Optional glob pattern, e.g. "csvs/lambda_scan_*.csv". Overrides --files.',
    )
    parser.add_argument(
        "--outdir",
        default="plots/lambda_scan",
        help="Output directory for plots.",
    )

    args = parser.parse_args()

    if args.glob:
        files = sorted(glob.glob(args.glob))
    else:
        files = args.files

    os.makedirs(args.outdir, exist_ok=True)

    datasets = []
    for path in files:
        if not os.path.exists(path):
            print("WARNING: missing file:", path)
            continue

        rows = read_csv(path)
        if len(rows) == 0:
            print("WARNING: no rows in:", path)
            continue

        datasets.append((nice_label(path), rows))
        print("Loaded", path, "with", len(rows), "lambda points")

    if len(datasets) == 0:
        raise RuntimeError("No valid datasets found.")

    # Main comparison plots.
    plot_vs_lambda(
        datasets,
        "delta_chi2",
        r"$\Delta\chi^2 = \chi^2_{\rm null} - \chi^2_{\rm BF}$",
        os.path.join(args.outdir, "delta_chi2_vs_lambda.png"),
    )

    plot_vs_lambda(
        datasets,
        "chi2_null",
        r"$\chi^2_{\rm null}$",
        os.path.join(args.outdir, "chi2_null_vs_lambda.png"),
    )

    plot_vs_lambda(
        datasets,
        "chi2_bf",
        r"$\chi^2_{\rm BF}$",
        os.path.join(args.outdir, "chi2_bf_vs_lambda.png"),
    )

    plot_vs_lambda(
        datasets,
        "norm_a_null",
        r"$||a_{\rm null}||$",
        os.path.join(args.outdir, "norm_a_null_vs_lambda.png"),
    )

    plot_vs_lambda(
        datasets,
        "norm_a_bf",
        r"$||a_{\rm BF}||$",
        os.path.join(args.outdir, "norm_a_bf_vs_lambda.png"),
    )

    plot_vs_lambda(
        datasets,
        "max_abs_a_null",
        r"$\max |a_{\rm null}|$",
        os.path.join(args.outdir, "max_abs_a_null_vs_lambda.png"),
    )

    plot_vs_lambda(
        datasets,
        "max_abs_a_bf",
        r"$\max |a_{\rm BF}|$",
        os.path.join(args.outdir, "max_abs_a_bf_vs_lambda.png"),
    )

    plot_vs_lambda(
        datasets,
        "dm2",
        r"Best-fit $\Delta m^2$ [eV$^2$]",
        os.path.join(args.outdir, "dm2_vs_lambda.png"),
    )

    plot_vs_lambda(
        datasets,
        "ue4",
        r"Best-fit $|U_{e4}|^2$",
        os.path.join(args.outdir, "ue4_vs_lambda.png"),
    )

    plot_vs_lambda(
        datasets,
        "umu4",
        r"Best-fit $|U_{\mu4}|^2$",
        os.path.join(args.outdir, "umu4_vs_lambda.png"),
    )

    # L-curves.
    plot_l_curve(
        datasets,
        os.path.join(args.outdir, "l_curve_bf_penalty_vs_residual.png"),
        which="bf",
    )

    plot_l_curve(
        datasets,
        os.path.join(args.outdir, "l_curve_null_penalty_vs_residual.png"),
        which="null",
    )

    # Swapped-axis residual scans.
    plot_lambda_vs_residual(
        datasets,
        os.path.join(args.outdir, "lambda_vs_resid_null.png"),
        which="null",
    )

    plot_lambda_vs_residual(
        datasets,
        os.path.join(args.outdir, "lambda_vs_resid_bf.png"),
        which="bf",
    )

    # Per-configuration decomposition.
    plot_chi2_decomposition(datasets, args.outdir)

    # Summary CSV.
    write_summary_table(
        datasets,
        os.path.join(args.outdir, "lambda1_summary.csv"),
    )

    print("\nWrote plots to:", args.outdir)
    print("Summary table:", os.path.join(args.outdir, "lambda1_summary.csv"))


if __name__ == "__main__":
    main()
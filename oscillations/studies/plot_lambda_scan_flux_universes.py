#!/usr/bin/env python3

import os
import csv
import glob
import argparse
import math
from collections import defaultdict

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LABEL_MAP = {
    "lambda_scan_fluxUniverseThrows_prodNueel_excludeRatio_maskRatio_allChunks_nullOnly": "prodNueel exclude ratio",
    "lambda_scan_fluxUniverseThrows_prodNueel_includeAll_maskNone_allChunks_nullOnly": "prodNueel include ratio",
    "lambda_scan_fluxUniverseThrows_prodNueel_profileOnlyRatio_maskNonRatio_allChunks_nullOnly": "prodNueel only ratio",
    "lambda_scan_fluxUniverseThrows_prodNueel_noRatio_includeAll_maskNone_allChunks_nullOnly": "prodNueel_noRatio",

    # optional: also map chunk files if you ever plot them directly
    "lambda_scan_fluxUniverseThrows_prodNueel_excludeRatio_maskRatio": "prodNueel exclude ratio",
    "lambda_scan_fluxUniverseThrows_prodNueel_includeAll_maskNone": "prodNueel include ratio",
    "lambda_scan_fluxUniverseThrows_prodNueel_profileOnlyRatio_maskNonRatio": "prodNueel only ratio",
    "lambda_scan_fluxUniverseThrows_prodNueel_noRatio_includeAll_maskNone": "prodNueel_noRatio",
}

def read_csv(path):
    rows = []

    with open(path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            out = {}
            for k, v in row.items():
                if v is None or v == "":
                    out[k] = np.nan
                    continue

                try:
                    out[k] = float(v)
                except Exception:
                    out[k] = v

            rows.append(out)

    return rows


def nice_label(path):
    base = os.path.basename(path)
    base = os.path.splitext(base)[0]

    # exact match first
    if base in LABEL_MAP:
        return LABEL_MAP[base]

    # fallback substring matching
    if "prodNueel_noRatio" in base:
        return "prodNueel_noRatio"
    if "prodNueel" in base and "profileOnlyRatio" in base:
        return "prodNueel only ratio"
    if "prodNueel" in base and "excludeRatio" in base:
        return "prodNueel exclude ratio"
    if "prodNueel" in base and "includeAll" in base:
        return "prodNueel include ratio"

    return base


def finite(x):
    try:
        return np.isfinite(float(x))
    except Exception:
        return False


def group_by_lambda(rows, ykey):
    """
    Return sorted lambdas and per-lambda arrays of y values over thrown universes.
    """
    groups = defaultdict(list)

    for r in rows:
        if "lambda" not in r or ykey not in r:
            continue

        lam = r["lambda"]
        y = r[ykey]

        if finite(lam) and finite(y):
            groups[float(lam)].append(float(y))

    lambdas = sorted(groups.keys())
    values = [np.asarray(groups[lam], dtype=float) for lam in lambdas]

    return np.asarray(lambdas, dtype=float), values


def summarize_by_lambda(rows, ykey):
    lambdas, values = group_by_lambda(rows, ykey)

    med = []
    p16 = []
    p84 = []
    p05 = []
    p95 = []
    ymin = []
    ymax = []
    n = []

    for arr in values:
        arr = arr[np.isfinite(arr)]
        n.append(len(arr))

        if len(arr) == 0:
            med.append(np.nan)
            p16.append(np.nan)
            p84.append(np.nan)
            p05.append(np.nan)
            p95.append(np.nan)
            ymin.append(np.nan)
            ymax.append(np.nan)
        else:
            med.append(np.percentile(arr, 50))
            p16.append(np.percentile(arr, 16))
            p84.append(np.percentile(arr, 84))
            p05.append(np.percentile(arr, 5))
            p95.append(np.percentile(arr, 95))
            ymin.append(np.min(arr))
            ymax.append(np.max(arr))

    return {
        "lambda": lambdas,
        "median": np.asarray(med),
        "p16": np.asarray(p16),
        "p84": np.asarray(p84),
        "p05": np.asarray(p05),
        "p95": np.asarray(p95),
        "min": np.asarray(ymin),
        "max": np.asarray(ymax),
        "n": np.asarray(n),
    }


def get_universe_curves(rows, xkey, ykey):
    """
    Return dict: thrown_universe -> sorted arrays x, y.
    """
    curves = defaultdict(list)

    for r in rows:
        if "thrown_universe" not in r:
            continue
        if xkey not in r or ykey not in r:
            continue
        if not finite(r[xkey]) or not finite(r[ykey]):
            continue

        u = int(r["thrown_universe"])
        curves[u].append((float(r[xkey]), float(r[ykey])))

    out = {}
    for u, pts in curves.items():
        pts = sorted(pts, key=lambda p: p[0])
        x = np.asarray([p[0] for p in pts])
        y = np.asarray([p[1] for p in pts])
        out[u] = (x, y)

    return out


def plot_summary_vs_lambda(datasets, ykey, ylabel, outpath, logx=True, show_individual=False):
    plt.figure(figsize=(8, 5.5))

    for label, rows in datasets:
        summary = summarize_by_lambda(rows, ykey)

        x = summary["lambda"]
        med = summary["median"]
        p16 = summary["p16"]
        p84 = summary["p84"]

        if len(x) == 0:
            print("WARNING: no valid values for", ykey, "in", label)
            continue

        if show_individual:
            curves = get_universe_curves(rows, "lambda", ykey)
            for u, (xu, yu) in curves.items():
                plt.plot(xu, yu, linewidth=0.7, alpha=0.18)

        plt.fill_between(x, p16, p84, alpha=0.20)
        plt.plot(x, med, marker="o", linewidth=2, label=label)

    if logx:
        plt.xscale("log")

    plt.xlabel(r"$\lambda$")
    plt.ylabel(ylabel)
    plt.grid(True, alpha=0.35)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

def plot_lambda_vs_summary(datasets, ykey, xlabel, outpath, show_individual=False):
    """
    Swapped-axis version of summary-vs-lambda plot.

    x-axis: metric, e.g. residual chi2
    y-axis: lambda

    The band is the 16--84% spread in the metric direction.
    """
    plt.figure(figsize=(7, 6))

    for label, rows in datasets:
        summary = summarize_by_lambda(rows, ykey)

        lam = summary["lambda"]
        med = summary["median"]
        p16 = summary["p16"]
        p84 = summary["p84"]

        if len(lam) == 0:
            print("WARNING: no valid values for", ykey, "in", label)
            continue

        if show_individual:
            curves = get_universe_curves(rows, "lambda", ykey)
            for u, (lam_u, y_u) in curves.items():
                plt.plot(y_u, lam_u, linewidth=0.7, alpha=0.18)

        line = plt.plot(
            med,
            lam,
            marker="o",
            linewidth=2,
            label=label,
        )[0]

        color = line.get_color()

        plt.fill_betweenx(
            lam,
            p16,
            p84,
            color=color,
            alpha=0.20,
            linewidth=0,
        )

    plt.yscale("log")
    plt.xlabel(xlabel)
    plt.ylabel(r"$\lambda$")
    plt.title(r"Stored flux-universe Asimov null residual $\chi^2$ scan (median with 16--84% band)")
    plt.grid(True, alpha=0.35)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

def plot_lcurve_summary(datasets, outpath, which="null", show_individual=False, show_band=False):
    """
    L-curve using median residual and penalty over thrown flux universes.

    x-axis: residual chi2
    y-axis: flux penalty

    If show_band=True, draw a shaded 16--84% penalty band around the median curve.
    """
    if which == "null":
        resid_key = "resid_null"
        penalty_key = "penalty_null"
        title = "Stored flux-universe Asimov null L-curve"
    elif which == "bf":
        resid_key = "resid_bf"
        penalty_key = "penalty_bf"
        title = "Stored flux-universe Asimov BF L-curve"
    else:
        raise ValueError("which must be null or bf")

    plt.figure(figsize=(7, 6))

    for label, rows in datasets:
        s_resid = summarize_by_lambda(rows, resid_key)
        s_pen = summarize_by_lambda(rows, penalty_key)

        lambdas = s_resid["lambda"]
        if len(lambdas) == 0:
            print("WARNING: no valid L-curve values for", label, which)
            continue

        # Median L-curve
        x = s_resid["median"]
        y = s_pen["median"]

        # Penalty band at each lambda.
        # This keeps the L-curve clean and gives a chi2_null-style shaded band.
        y_low = s_pen["p16"]
        y_high = s_pen["p84"]

        if show_individual:
            resid_curves = get_universe_curves(rows, "lambda", resid_key)
            pen_curves = get_universe_curves(rows, "lambda", penalty_key)

            common = sorted(set(resid_curves.keys()).intersection(set(pen_curves.keys())))
            for u in common:
                lam_r, res_u = resid_curves[u]
                lam_p, pen_u = pen_curves[u]

                if len(lam_r) != len(lam_p) or np.max(np.abs(lam_r - lam_p)) > 1e-12:
                    continue

                plt.plot(res_u, pen_u, linewidth=0.8, alpha=0.15)
                plt.scatter(res_u, pen_u, s=8, alpha=0.15)

        line = plt.plot(
            x,
            y,
            marker="o",
            linewidth=2,
            label=label,
        )[0]

        color = line.get_color()

        if show_band:
            plt.fill_between(
                x,
                y_low,
                y_high,
                color=color,
                alpha=0.20,
                linewidth=0,
            )

        for xi, yi, li in zip(x, y, lambdas):
            plt.text(xi, yi, "{:g}".format(li), fontsize=7)

    plt.xlabel("Residual " + r"$\chi^2$")
    plt.ylabel("Flux penalty")

    if show_band:
        plt.title(title + r"  (median with 16--84% band)")
    else:
        plt.title(title + r"  (median over flux universes)")

    plt.grid(True, alpha=0.35)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def plot_delta_chi2_distribution_at_lambda(datasets, target_lambda, outpath):
    plt.figure(figsize=(8, 5.5))

    any_valid = False

    for label, rows in datasets:
        vals = []

        for r in rows:
            if "lambda" not in r or "delta_chi2" not in r:
                continue
            if not finite(r["lambda"]) or not finite(r["delta_chi2"]):
                continue

            if abs(float(r["lambda"]) - target_lambda) < 1e-9:
                vals.append(float(r["delta_chi2"]))

        if len(vals) == 0:
            continue

        any_valid = True
        plt.hist(vals, bins=25, histtype="step", linewidth=2, label=label)

    if not any_valid:
        plt.close()
        print("WARNING: no valid delta_chi2 values at lambda =", target_lambda)
        return

    plt.xlabel(r"$\Delta\chi^2 = \chi^2_{\rm null} - \chi^2_{\rm BF}$")
    plt.ylabel("Number of flux universes")
    plt.title(r"Stored flux-universe Asimovs at $\lambda={}$".format(target_lambda))
    plt.grid(True, alpha=0.35)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def write_lambda_summary(datasets, outpath):
    keys = [
        "chi2_null",
        "resid_null",
        "penalty_null",
        "norm_a_null",
        "max_abs_a_null",
        "chi2_bf",
        "resid_bf",
        "penalty_bf",
        "norm_a_bf",
        "max_abs_a_bf",
        "delta_chi2",
        "dm2",
        "ue4",
        "umu4",
        "utau4",
    ]

    with open(outpath, "w") as f:
        f.write("label,lambda,n_universes,key,median,p16,p84,p05,p95,min,max\n")

        for label, rows in datasets:
            lambdas = sorted(
                set(
                    float(r["lambda"])
                    for r in rows
                    if "lambda" in r and finite(r["lambda"])
                )
            )

            for key in keys:
                if key not in rows[0]:
                    continue

                for lam in lambdas:
                    vals = []
                    for r in rows:
                        if "lambda" not in r or key not in r:
                            continue
                        if not finite(r["lambda"]) or not finite(r[key]):
                            continue
                        if abs(float(r["lambda"]) - lam) < 1e-9:
                            vals.append(float(r[key]))

                    if len(vals) == 0:
                        continue

                    arr = np.asarray(vals, dtype=float)
                    f.write(
                        "{},{:.8g},{},{},{:.8g},{:.8g},{:.8g},{:.8g},{:.8g},{:.8g},{:.8g}\n".format(
                            label,
                            lam,
                            len(arr),
                            key,
                            np.percentile(arr, 50),
                            np.percentile(arr, 16),
                            np.percentile(arr, 84),
                            np.percentile(arr, 5),
                            np.percentile(arr, 95),
                            np.min(arr),
                            np.max(arr),
                        )
                    )


def main():
    parser = argparse.ArgumentParser(
        description="Plot stored-flux-universe Asimov lambda-scan diagnostics."
    )

    parser.add_argument(
        "--files",
        nargs="+",
        default=None,
        help="CSV files to plot.",
    )

    parser.add_argument(
        "--glob",
        default=None,
        help='Optional glob pattern, e.g. "csvs/lambda_scan_fluxUniverses_*.csv". Overrides --files.',
    )

    parser.add_argument(
        "--outdir",
        default="plots/lambda_scan_flux_universes",
        help="Output directory for plots.",
    )

    parser.add_argument(
        "--show-individual",
        action="store_true",
        help="Overlay individual thrown-universe curves faintly behind median curves.",
    )

    parser.add_argument(
        "--target-lambda",
        type=float,
        default=1.0,
        help="Lambda value for delta-chi2 distribution plot.",
    )

    args = parser.parse_args()

    if args.glob:
        files = sorted(glob.glob(args.glob))
    elif args.files:
        files = args.files
    else:
        files = [
            "csvs/lambda_scan_fluxUniverses_prodNueel_excludeRatio_maskRatio_all_nullOnly.csv"
        ]

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

        universes = sorted(
            set(
                int(r["thrown_universe"])
                for r in rows
                if "thrown_universe" in r and finite(r["thrown_universe"])
            )
        )

        lambdas = sorted(
            set(
                float(r["lambda"])
                for r in rows
                if "lambda" in r and finite(r["lambda"])
            )
        )

        print(
            "Loaded",
            path,
            "rows =",
            len(rows),
            "universes =",
            len(universes),
            "lambdas =",
            lambdas,
        )

    if len(datasets) == 0:
        raise RuntimeError("No valid datasets found.")

    plot_summary_vs_lambda(
        datasets,
        "chi2_null",
        r"Null total $\chi^2$",
        os.path.join(args.outdir, "chi2_null_vs_lambda_median_band.png"),
        show_individual=args.show_individual,
    )

    plot_summary_vs_lambda(
        datasets,
        "resid_null",
        r"Null residual $\chi^2$",
        os.path.join(args.outdir, "resid_null_vs_lambda_median_band.png"),
        show_individual=args.show_individual,
    )

    plot_summary_vs_lambda(
        datasets,
        "penalty_null",
        "Null flux penalty",
        os.path.join(args.outdir, "penalty_null_vs_lambda_median_band.png"),
        show_individual=args.show_individual,
    )

    plot_summary_vs_lambda(
        datasets,
        "norm_a_null",
        r"$||a_{\rm null}||$",
        os.path.join(args.outdir, "norm_a_null_vs_lambda_median_band.png"),
        show_individual=args.show_individual,
    )

    plot_summary_vs_lambda(
        datasets,
        "max_abs_a_null",
        r"$\max |a_{\rm null}|$",
        os.path.join(args.outdir, "max_abs_a_null_vs_lambda_median_band.png"),
        show_individual=args.show_individual,
    )

    # BF plots only work if the CSV has finite BF values.
    plot_summary_vs_lambda(
        datasets,
        "delta_chi2",
        r"$\Delta\chi^2 = \chi^2_{\rm null} - \chi^2_{\rm BF}$",
        os.path.join(args.outdir, "delta_chi2_vs_lambda_median_band.png"),
        show_individual=args.show_individual,
    )

    plot_summary_vs_lambda(
        datasets,
        "chi2_bf",
        r"BF total $\chi^2$",
        os.path.join(args.outdir, "chi2_bf_vs_lambda_median_band.png"),
        show_individual=args.show_individual,
    )

    plot_summary_vs_lambda(
        datasets,
        "dm2",
        r"BF $\Delta m^2$ [eV$^2$]",
        os.path.join(args.outdir, "dm2_vs_lambda_median_band.png"),
        show_individual=args.show_individual,
    )

    plot_summary_vs_lambda(
        datasets,
        "ue4",
        r"BF $|U_{e4}|^2$",
        os.path.join(args.outdir, "ue4_vs_lambda_median_band.png"),
        show_individual=args.show_individual,
    )

    plot_summary_vs_lambda(
        datasets,
        "umu4",
        r"BF $|U_{\mu4}|^2$",
        os.path.join(args.outdir, "umu4_vs_lambda_median_band.png"),
        show_individual=args.show_individual,
    )
 
    plot_lambda_vs_summary(
        datasets,
        "resid_null",
        r"Residual $\chi^2$",
        os.path.join(args.outdir, "lambda_vs_resid_null_median_band.png"),
        show_individual=args.show_individual,
    )

    # Median-only null L-curve
    plot_lcurve_summary(
        datasets,
        os.path.join(args.outdir, "l_curve_null_median_penalty_vs_residual.png"),
        which="null",
        show_individual=args.show_individual,
        show_band=False,
    )

    # Median + shaded-band null L-curve
    plot_lcurve_summary(
        datasets,
        os.path.join(args.outdir, "l_curve_null_median_band_penalty_vs_residual.png"),
        which="null",
        show_individual=args.show_individual,
        show_band=True,
    )

    # Median-only BF L-curve
    plot_lcurve_summary(
        datasets,
        os.path.join(args.outdir, "l_curve_bf_median_penalty_vs_residual.png"),
        which="bf",
        show_individual=args.show_individual,
        show_band=False,
    )

    # Median + shaded-band BF L-curve
    plot_lcurve_summary(
        datasets,
        os.path.join(args.outdir, "l_curve_bf_median_band_penalty_vs_residual.png"),
        which="bf",
        show_individual=args.show_individual,
        show_band=True,
    )

    plot_delta_chi2_distribution_at_lambda(
        datasets,
        args.target_lambda,
        os.path.join(args.outdir, "delta_chi2_distribution_lambda{}.png".format(args.target_lambda)),
    )

    write_lambda_summary(
        datasets,
        os.path.join(args.outdir, "lambda_summary_by_universe_distribution.csv"),
    )

    print("\nWrote plots to:", args.outdir)
    print(
        "Summary table:",
        os.path.join(args.outdir, "lambda_summary_by_universe_distribution.csv"),
    )


if __name__ == "__main__":
    main()
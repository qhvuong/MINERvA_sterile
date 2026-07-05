#!/usr/bin/env python3

import os
import glob
import csv
import argparse
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_rows(paths):
    rows = []

    for path in paths:
        with open(path, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                out = {}
                for k, v in row.items():
                    try:
                        out[k] = float(v)
                    except Exception:
                        out[k] = v

                out["source_file"] = os.path.basename(path)
                rows.append(out)

    return rows


def summarize_by_nprof(rows, key):
    nvals = sorted(set(int(r["profile_n_universes"]) for r in rows))

    out = {
        "nprof": [],
        "median": [],
        "p16": [],
        "p84": [],
        "mean": [],
        "std": [],
    }

    for n in nvals:
        vals = np.array(
            [
                float(r[key])
                for r in rows
                if int(r["profile_n_universes"]) == n and np.isfinite(float(r[key]))
            ],
            dtype=float,
        )

        if len(vals) == 0:
            continue

        out["nprof"].append(n)
        out["median"].append(np.median(vals))
        out["p16"].append(np.percentile(vals, 16))
        out["p84"].append(np.percentile(vals, 84))
        out["mean"].append(np.mean(vals))
        out["std"].append(np.std(vals))

    for k in out:
        out[k] = np.array(out[k])

    return out


def write_summary(rows, out_csv):
    keys = [
        "chi2_null",
        "resid_null",
        "penalty_null",
        "norm_a_null",
        "max_abs_a_null",
    ]

    nvals = sorted(set(int(r["profile_n_universes"]) for r in rows))

    with open(out_csv, "w") as f:
        writer = csv.writer(f)

        header = ["profile_n_universes", "n_throws"]
        for key in keys:
            header += [
                key + "_median",
                key + "_p16",
                key + "_p84",
                key + "_mean",
                key + "_std",
            ]

        writer.writerow(header)

        for n in nvals:
            subset = [r for r in rows if int(r["profile_n_universes"]) == n]
            row_out = [n, len(subset)]

            for key in keys:
                vals = np.array(
                    [
                        float(r[key])
                        for r in subset
                        if key in r and np.isfinite(float(r[key]))
                    ],
                    dtype=float,
                )

                if len(vals) == 0:
                    row_out += [np.nan, np.nan, np.nan, np.nan, np.nan]
                else:
                    row_out += [
                        np.median(vals),
                        np.percentile(vals, 16),
                        np.percentile(vals, 84),
                        np.mean(vals),
                        np.std(vals),
                    ]

            writer.writerow(row_out)


def plot_summary(rows, key, ylabel, outpath, title):
    s = summarize_by_nprof(rows, key)

    if len(s["nprof"]) == 0:
        print("WARNING: no valid values for", key)
        return

    x = s["nprof"]
    med = s["median"]
    p16 = s["p16"]
    p84 = s["p84"]

    plt.figure(figsize=(7, 6))

    line = plt.plot(
        x,
        med,
        marker="o",
        linewidth=2,
        label="Flux-universe Asimov median",
    )[0]

    color = line.get_color()

    plt.fill_between(
        x,
        p16,
        p84,
        color=color,
        alpha=0.20,
        linewidth=0,
        label="16--84% spread",
    )

    plt.xlabel("Number of flux universes in profiling basis")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, alpha=0.35)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

def plot_vs_residual(rows, ykey, ylabel, outpath, title):
    """
    L-curve-like profile-N summary.

    x-axis: median residual chi2
    y-axis: median quantity, e.g. penalty or pull norm

    Band is drawn in y direction using 16--84% spread.
    """
    s_x = summarize_by_nprof(rows, "resid_null")
    s_y = summarize_by_nprof(rows, ykey)

    if len(s_x["nprof"]) == 0 or len(s_y["nprof"]) == 0:
        print("WARNING: no valid values for residual-vs-{}".format(ykey))
        return

    # Assumes same Nprof grid.
    x = s_x["median"]
    y = s_y["median"]
    y_low = s_y["p16"]
    y_high = s_y["p84"]
    nprof = s_x["nprof"]

    plt.figure(figsize=(7, 6))

    line = plt.plot(
        x,
        y,
        marker="o",
        linewidth=2,
        label="Flux-universe Asimov median",
    )[0]

    color = line.get_color()

    plt.fill_between(
        x,
        y_low,
        y_high,
        color=color,
        alpha=0.20,
        linewidth=0,
        label="16--84% spread",
    )

    for xi, yi, ni in zip(x, y, nprof):
        plt.text(xi, yi, "{}".format(int(ni)), fontsize=7)

    plt.xlabel(r"Residual $\chi^2$")
    plt.ylabel(ylabel)
    plt.title(title + r" ($N_{\rm prof}$ labels shown)")
    plt.grid(True, alpha=0.35)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--pattern",
        required=True,
        help="Glob pattern for profileN CSV files.",
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for summary and plots.",
    )

    parser.add_argument(
        "--label",
        default="prodNueel exclude ratio",
        help="Label used in plot titles.",
    )

    args = parser.parse_args()

    paths = sorted(glob.glob(args.pattern))
    if len(paths) == 0:
        raise RuntimeError("No files matched pattern: {}".format(args.pattern))

    os.makedirs(args.outdir, exist_ok=True)

    print("Matched files:")
    for p in paths:
        print(" ", p)

    rows = read_rows(paths)
    print("Loaded rows =", len(rows))

    summary_csv = os.path.join(args.outdir, "profileN_summary.csv")
    write_summary(rows, summary_csv)
    print("Wrote", summary_csv)

    title_prefix = args.label + r", $\lambda=1$ null-only"

    plot_summary(
        rows,
        "chi2_null",
        r"Null $\chi^2$",
        os.path.join(args.outdir, "profileN_chi2_null.png"),
        title_prefix + r": total null $\chi^2$",
    )

    plot_summary(
        rows,
        "resid_null",
        r"Residual $\chi^2$",
        os.path.join(args.outdir, "profileN_resid_null.png"),
        title_prefix + r": residual $\chi^2$",
    )

    plot_summary(
        rows,
        "penalty_null",
        "Flux penalty",
        os.path.join(args.outdir, "profileN_penalty_null.png"),
        title_prefix + ": flux penalty",
    )

    plot_summary(
        rows,
        "norm_a_null",
        r"$||a_{\rm null}||$",
        os.path.join(args.outdir, "profileN_norm_a_null.png"),
        title_prefix + r": flux-pull norm",
    )

    plot_summary(
        rows,
        "max_abs_a_null",
        r"$\max |a_{\rm null}|$",
        os.path.join(args.outdir, "profileN_max_abs_a_null.png"),
        title_prefix + r": max flux pull",
    )

    # L-curve-like plots with residual chi2 on x-axis.
    plot_vs_residual(
        rows,
        "chi2_null",
        r"Null $\chi^2$",
        os.path.join(args.outdir, "profileN_chi2_null_vs_resid_null.png"),
        title_prefix + r": total null $\chi^2$ vs residual $\chi^2$",
    )

    plot_vs_residual(
        rows,
        "penalty_null",
        "Flux penalty",
        os.path.join(args.outdir, "profileN_penalty_null_vs_resid_null.png"),
        title_prefix + r": flux penalty vs residual $\chi^2$",
    )

    plot_vs_residual(
        rows,
        "norm_a_null",
        r"$||a_{\rm null}||$",
        os.path.join(args.outdir, "profileN_norm_a_null_vs_resid_null.png"),
        title_prefix + r": flux-pull norm vs residual $\chi^2$",
    )

    plot_vs_residual(
        rows,
        "max_abs_a_null",
        r"$\max |a_{\rm null}|$",
        os.path.join(args.outdir, "profileN_max_abs_a_null_vs_resid_null.png"),
        title_prefix + r": max flux pull vs residual $\chi^2$",
    )

    print("Wrote plots to", args.outdir)


if __name__ == "__main__":
    main()
#!/usr/bin/env python3

import os
import csv
import argparse
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


SAMPLES = [
    {
        "label": "prodNueel exclude ratio",
        "summary": "/exp/minerva/data/users/qvuong/surfaces/plots/profileN_scan_prodNueel_excludeRatio_maskRatio_lambda1_nullOnly/profileN_summary.csv",
    },
    {
        "label": "prodNueel include ratio",
        "summary": "/exp/minerva/data/users/qvuong/surfaces/plots/profileN_scan_prodNueel_includeRatio_maskNone_lambda1_nullOnly/profileN_summary.csv",
    },
    {
        "label": "prodNueel_noRatio",
        "summary": "/exp/minerva/data/users/qvuong/surfaces/plots/profileN_scan_prodNueel_noRatio_includeAll_maskNone_lambda1_nullOnly/profileN_summary.csv",
    },
]


def read_summary(path):
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

    rows = sorted(rows, key=lambda r: int(r["profile_n_universes"]))
    return rows


def get_array(rows, key):
    return np.array([float(r[key]) for r in rows], dtype=float)


def plot_quantity(samples, key, ylabel, title, outpath):
    plt.figure(figsize=(7, 6))

    for sample in samples:
        rows = read_summary(sample["summary"])

        x = get_array(rows, "profile_n_universes")
        med = get_array(rows, key + "_median")
        p16 = get_array(rows, key + "_p16")
        p84 = get_array(rows, key + "_p84")

        line = plt.plot(
            x,
            med,
            marker="o",
            linewidth=2,
            label=sample["label"],
        )[0]

        color = line.get_color()

        plt.fill_between(
            x,
            p16,
            p84,
            color=color,
            alpha=0.18,
            linewidth=0,
        )

    plt.xlabel("Number of flux universes in profiling basis")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, alpha=0.35)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def plot_vs_residual(samples, ykey, ylabel, title, outpath):
    plt.figure(figsize=(7, 6))

    for sample in samples:
        rows = read_summary(sample["summary"])

        x = get_array(rows, "resid_null_median")
        y = get_array(rows, ykey + "_median")
        y16 = get_array(rows, ykey + "_p16")
        y84 = get_array(rows, ykey + "_p84")
        nprof = get_array(rows, "profile_n_universes")

        line = plt.plot(
            x,
            y,
            marker="o",
            linewidth=2,
            label=sample["label"],
        )[0]

        color = line.get_color()

        plt.fill_between(
            x,
            y16,
            y84,
            color=color,
            alpha=0.18,
            linewidth=0,
        )

        for xi, yi, ni in zip(x, y, nprof):
            plt.text(xi, yi, str(int(ni)), fontsize=6)

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
        "--outdir",
        default="/exp/minerva/data/users/qvuong/surfaces/plots/profileN_scan_allsamples_lambda1_nullOnly",
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    for s in SAMPLES:
        if not os.path.exists(s["summary"]):
            raise RuntimeError("Missing summary file: {}".format(s["summary"]))

    plot_quantity(
        SAMPLES,
        "resid_null",
        r"Residual $\chi^2$",
        r"Null-only profile-$N$ scan: residual $\chi^2$",
        os.path.join(args.outdir, "profileN_resid_null_allsamples.png"),
    )

    plot_quantity(
        SAMPLES,
        "chi2_null",
        r"Null $\chi^2$",
        r"Null-only profile-$N$ scan: total null $\chi^2$",
        os.path.join(args.outdir, "profileN_chi2_null_allsamples.png"),
    )

    plot_quantity(
        SAMPLES,
        "penalty_null",
        "Flux penalty",
        r"Null-only profile-$N$ scan: flux penalty",
        os.path.join(args.outdir, "profileN_penalty_null_allsamples.png"),
    )

    plot_quantity(
        SAMPLES,
        "norm_a_null",
        r"$||a_{\rm null}||$",
        r"Null-only profile-$N$ scan: flux-pull norm",
        os.path.join(args.outdir, "profileN_norm_a_null_allsamples.png"),
    )

    plot_quantity(
        SAMPLES,
        "max_abs_a_null",
        r"$\max |a_{\rm null}|$",
        r"Null-only profile-$N$ scan: max flux pull",
        os.path.join(args.outdir, "profileN_max_abs_a_null_allsamples.png"),
    )

    plot_vs_residual(
        SAMPLES,
        "penalty_null",
        "Flux penalty",
        r"Null-only profile-$N$ scan: penalty vs residual $\chi^2$",
        os.path.join(args.outdir, "profileN_penalty_null_vs_resid_null_allsamples.png"),
    )

    plot_vs_residual(
        SAMPLES,
        "norm_a_null",
        r"$||a_{\rm null}||$",
        r"Null-only profile-$N$ scan: flux-pull norm vs residual $\chi^2$",
        os.path.join(args.outdir, "profileN_norm_a_null_vs_resid_null_allsamples.png"),
    )

    plot_vs_residual(
        SAMPLES,
        "max_abs_a_null",
        r"$\max |a_{\rm null}|$",
        r"Null-only profile-$N$ scan: max flux pull vs residual $\chi^2$",
        os.path.join(args.outdir, "profileN_max_abs_a_null_vs_resid_null_allsamples.png"),
    )

    print("Wrote all-sample profile-N plots to", args.outdir)


if __name__ == "__main__":
    main()
#!/usr/bin/env python3

import os
import sys
import argparse
import csv
import json
import shutil
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args_strip_project_args():
    parser = argparse.ArgumentParser(
        description="Check nominal vs unoscillated-template projection vs raw BF oscillated prediction."
    )

    parser.add_argument("--hist-config-tag", default="prodNueel")
    parser.add_argument("--outdir", default="/exp/minerva/data/users/qvuong/surfaces/plots/check_unoscillated_projection")

    parser.add_argument("--dm2", type=float, default=13.8977)
    parser.add_argument("--ue4", type=float, default=0.05914)
    parser.add_argument("--umu4", type=float, default=0.03146)
    parser.add_argument("--utau4", type=float, default=0.0)

    parser.add_argument(
        "--zero-dm2",
        type=float,
        default=13.8977,
        help="dm2 value used for zero-mixing projection. Mixing is set to zero, so dm2 should not matter unless code has a template-side effect.",
    )

    parser.add_argument(
        "--samples",
        default="all",
        help='Comma-separated sample names to plot, or "all".',
    )

    args, remaining = parser.parse_known_args()

    # Hide script-only args from AnalysisConfig/project parsers.
    sys.argv = [sys.argv[0]] + remaining

    return args


args = parse_args_strip_project_args()

import ROOT
ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

try:
    import PlotUtils
except Exception:
    PlotUtils = None

from config.AnalysisConfig import AnalysisConfig
from tools.StitchedHistogram import StitchedHistogram


def mkdir(path):
    if path:
        os.makedirs(path, exist_ok=True)


def safe_div(num, den):
    num = np.asarray(num, dtype=float)
    den = np.asarray(den, dtype=float)
    out = np.full_like(num, np.nan, dtype=float)
    good = np.abs(den) > 0
    out[good] = num[good] / den[good]
    return out


def rms_finite(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return np.nan
    return float(np.sqrt(np.mean(x * x)))


def max_abs_finite(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return np.nan
    return float(np.max(np.abs(x)))


def get_nominal_vector(histogram):
    return np.array(histogram.GetMCHistogram())[1:-1].astype(float)


def get_oscillated_vector(histogram):
    return np.array(histogram.GetOscillatedHistogram())[1:-1].astype(float)


def load_prediction(file_path, osc=None):
    """
    Load stitched prediction vector.

    osc=None:
        nominal unoscillated MC from GetMCHistogram()

    osc=(dm2, ue4, umu4, utau4):
        calls OscillateHistogram, then reads the actual oscillated prediction
        from GetOscillatedHistogram().
    """
    h = StitchedHistogram("sample")
    h.Load(file_path)

    if osc is None:
        return get_nominal_vector(h)

    dm2, ue4, umu4, utau4 = osc
    h.OscillateHistogram(dm2, ue4, umu4, utau4)

    return get_oscillated_vector(h)


def load_hist_config(hist_config_path):
    with open(hist_config_path, "r") as f:
        cfg = json.load(f)

    # Sort by global start bin.
    items = []
    for sample, info in cfg.items():
        if "start" not in info or "end" not in info:
            continue
        items.append((sample, int(info["start"]), int(info["end"])))

    items = sorted(items, key=lambda x: x[1])
    return items


def selected_samples(all_samples, sample_spec):
    if sample_spec in [None, "", "all", "All", "ALL"]:
        return set([s[0] for s in all_samples])

    keep = set([x.strip() for x in sample_spec.split(",") if x.strip()])
    return keep


def write_csv(rows, path):
    if len(rows) == 0:
        return

    mkdir(os.path.dirname(path))

    with open(path, "w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def plot_sample(sample, x, nominal, zero_proj, bf_raw, outdir):
    mkdir(outdir)

    ratio_zero = safe_div(zero_proj, nominal)
    ratio_bf = safe_div(bf_raw, nominal)

    frac_zero = ratio_zero - 1.0
    frac_bf = ratio_bf - 1.0

    # Absolute prediction overlay.
    plt.figure(figsize=(8, 5.5))
    plt.plot(x, nominal, marker="o", linewidth=2, label="nominal 1D MC")
    plt.plot(x, zero_proj, marker="s", linewidth=1.7, linestyle="--", label="zero-mixing projection")
    plt.plot(x, bf_raw, marker="^", linewidth=1.7, linestyle=":", label="raw BF oscillated")
    plt.xlabel("Local bin")
    plt.ylabel("Prediction")
    plt.title(sample + ": prediction comparison")
    plt.grid(True, alpha=0.35)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, sample + "_prediction_overlay.png"), dpi=200)
    plt.close()

    # Ratio to nominal.
    plt.figure(figsize=(8, 5.5))
    plt.axhline(1.0, color="black", linewidth=1)
    plt.plot(x, ratio_zero, marker="s", linewidth=1.7, linestyle="--", label="zero-mixing / nominal")
    plt.plot(x, ratio_bf, marker="^", linewidth=1.7, linestyle=":", label="raw BF / nominal")
    plt.xlabel("Local bin")
    plt.ylabel("Ratio to nominal")
    plt.title(sample + ": ratio to nominal")
    plt.grid(True, alpha=0.35)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, sample + "_ratio_to_nominal.png"), dpi=200)
    plt.close()

    # Fractional difference.
    plt.figure(figsize=(8, 5.5))
    plt.axhline(0.0, color="black", linewidth=1)
    plt.plot(x, frac_zero, marker="s", linewidth=1.7, linestyle="--", label="zero-mixing projection")
    plt.plot(x, frac_bf, marker="^", linewidth=1.7, linestyle=":", label="raw BF oscillated")
    plt.xlabel("Local bin")
    plt.ylabel("(prediction - nominal) / nominal")
    plt.title(sample + ": fractional difference from nominal")
    plt.grid(True, alpha=0.35)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, sample + "_fractional_difference.png"), dpi=200)
    plt.close()


def main():
    ccnueroot = os.environ.get("CCNUEROOT")
    if ccnueroot is None:
        raise RuntimeError("CCNUEROOT is not set")

    plot_tag = args.hist_config_tag
    if plot_tag in [None, "", "none"]:
        plot_tag = getattr(AnalysisConfig, "hist_config_tag", "default")

    root_file = os.path.join(
        ccnueroot,
        "oscillations",
        "rootfiles",
        "NuE_stitched_hists_{}.root".format(plot_tag),
    )

    hist_config = "HIST_CONFIG_{}.json".format(plot_tag)
    if not os.path.exists(hist_config):
        raise RuntimeError("Missing hist config: {}".format(hist_config))

    # Existing fitter code expects HIST_CONFIG.json in local directory.
    shutil.copyfile(hist_config, "HIST_CONFIG.json")

    samples = load_hist_config(hist_config)
    keep_samples = selected_samples(samples, args.samples)

    mkdir(args.outdir)

    print("\n===== check_unoscillated_projection setup =====")
    print("plot_tag      =", plot_tag)
    print("root_file     =", root_file)
    print("hist_config   =", hist_config)
    print("outdir        =", args.outdir)
    print("BF dm2        =", args.dm2)
    print("BF ue4        =", args.ue4)
    print("BF umu4       =", args.umu4)
    print("BF utau4      =", args.utau4)
    print("zero dm2      =", args.zero_dm2)
    print("samples       =", args.samples)
    print("available samples:")
    for s, start, end in samples:
        print("  {:25s} {:4d} {:4d}".format(s, start, end))

    nominal = load_prediction(root_file, osc=None)

    # This is the key diagnostic:
    # If this differs from nominal, the zero-mixing call is changing the histogram
    # through the oscillation/template projection machinery.
    zero_proj = load_prediction(
        root_file,
        osc=(args.zero_dm2, 0.0, 0.0, 0.0),
    )

    bf_raw = load_prediction(
        root_file,
        osc=(args.dm2, args.ue4, args.umu4, args.utau4),
    )

    print("max |zero - nominal| =", np.max(np.abs(zero_proj - nominal)))
    print("max |bf - nominal|   =", np.max(np.abs(bf_raw - nominal)))

    bin_rows = []
    summary_rows = []

    for sample, start, end in samples:
        if sample not in keep_samples:
            continue

        nom = nominal[start:end + 1]
        zer = zero_proj[start:end + 1]
        bf = bf_raw[start:end + 1]

        x = np.arange(1, len(nom) + 1)

        ratio_zero = safe_div(zer, nom)
        ratio_bf = safe_div(bf, nom)

        frac_zero = ratio_zero - 1.0
        frac_bf = ratio_bf - 1.0

        sample_outdir = os.path.join(args.outdir, sample)
        plot_sample(sample, x, nom, zer, bf, sample_outdir)

        for ilocal in range(len(nom)):
            bin_rows.append({
                "sample": sample,
                "local_bin": ilocal + 1,
                "global_bin_zero_based": start + ilocal,
                "global_bin_one_based": start + ilocal + 1,
                "nominal": float(nom[ilocal]),
                "zero_mixing_projection": float(zer[ilocal]),
                "raw_bf_oscillated": float(bf[ilocal]),
                "zero_over_nominal": float(ratio_zero[ilocal]) if np.isfinite(ratio_zero[ilocal]) else np.nan,
                "bf_over_nominal": float(ratio_bf[ilocal]) if np.isfinite(ratio_bf[ilocal]) else np.nan,
                "frac_zero_minus_nominal": float(frac_zero[ilocal]) if np.isfinite(frac_zero[ilocal]) else np.nan,
                "frac_bf_minus_nominal": float(frac_bf[ilocal]) if np.isfinite(frac_bf[ilocal]) else np.nan,
            })

        summary_rows.append({
            "sample": sample,
            "n_bins": len(nom),
            "nominal_sum": float(np.sum(nom)),
            "zero_projection_sum": float(np.sum(zer)),
            "raw_bf_sum": float(np.sum(bf)),
            "zero_projection_sum_over_nominal": float(np.sum(zer) / np.sum(nom)) if np.sum(nom) != 0 else np.nan,
            "raw_bf_sum_over_nominal": float(np.sum(bf) / np.sum(nom)) if np.sum(nom) != 0 else np.nan,
            "zero_max_abs_frac_diff": max_abs_finite(frac_zero),
            "zero_rms_frac_diff": rms_finite(frac_zero),
            "bf_max_abs_frac_diff": max_abs_finite(frac_bf),
            "bf_rms_frac_diff": rms_finite(frac_bf),
        })

    write_csv(
        bin_rows,
        os.path.join(args.outdir, "per_bin_comparison.csv"),
    )

    write_csv(
        summary_rows,
        os.path.join(args.outdir, "sample_summary.csv"),
    )

    print("\n===== sample summary =====")
    for r in summary_rows:
        print(
            "{sample:25s}  zero max|frac|={zero_max_abs_frac_diff:.4g}  "
            "zero RMS={zero_rms_frac_diff:.4g}  "
            "BF max|frac|={bf_max_abs_frac_diff:.4g}  "
            "BF RMS={bf_rms_frac_diff:.4g}".format(**r)
        )

    print("\nWrote:")
    print("  ", os.path.join(args.outdir, "per_bin_comparison.csv"))
    print("  ", os.path.join(args.outdir, "sample_summary.csv"))
    print("  sample plots under:", args.outdir)

    print("\nInterpretation guide:")
    print("  If zero-mixing projection differs from nominal, investigate template projection / normalization / bin mapping.")
    print("  If zero-mixing projection matches nominal but raw BF wiggles, investigate oscillation weighting or ratio amplification.")
    print("  For ratio samples, large BF wiggles can be amplified by a small/noisy CCnue denominator.")


if __name__ == "__main__":
    main()

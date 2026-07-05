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
        description="Check whether CC ratio BF wiggles come from numerator or denominator."
    )

    parser.add_argument("--hist-config-tag", default="prodNueel_noRatio")
    parser.add_argument(
        "--outdir",
        default="/exp/minerva/data/users/qvuong/surfaces/plots/check_ratio_components",
    )

    parser.add_argument("--dm2", type=float, default=13.8977)
    parser.add_argument("--ue4", type=float, default=0.05914)
    parser.add_argument("--umu4", type=float, default=0.03146)
    parser.add_argument("--utau4", type=float, default=0.0)

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

    out = {}
    for sample, info in cfg.items():
        if "start" not in info or "end" not in info:
            continue
        out[sample] = (int(info["start"]), int(info["end"]))

    return out


def slice_sample(vec, cfg, sample):
    if sample not in cfg:
        raise RuntimeError("Sample {} not found in HIST_CONFIG".format(sample))

    start, end = cfg[sample]
    return vec[start:end + 1]


def write_csv(rows, path):
    if len(rows) == 0:
        return

    mkdir(os.path.dirname(path))

    with open(path, "w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def plot_ratio_component(mode, x, num_nom, den_nom, ratio_nom, num_bf, den_bf, ratio_bf, outdir):
    mkdir(outdir)

    num_bf_over_nom = safe_div(num_bf, num_nom)
    den_bf_over_nom = safe_div(den_bf, den_nom)
    ratio_bf_over_nom = safe_div(ratio_bf, ratio_nom)

    # 1. Absolute numerator/denominator.
    plt.figure(figsize=(8, 5.5))
    plt.plot(x, num_nom, marker="o", linewidth=2, label="CCnumu nominal")
    plt.plot(x, num_bf, marker="o", linewidth=1.7, linestyle=":", label="CCnumu raw BF")
    plt.plot(x, den_nom, marker="s", linewidth=2, label="CCnue nominal")
    plt.plot(x, den_bf, marker="s", linewidth=1.7, linestyle=":", label="CCnue raw BF")
    plt.xlabel("Local bin")
    plt.ylabel("Prediction")
    plt.title(mode + ": numerator and denominator")
    plt.grid(True, alpha=0.35)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, mode + "_num_den_prediction_overlay.png"), dpi=200)
    plt.close()

    # 2. BF / nominal for numerator, denominator, and ratio.
    plt.figure(figsize=(8, 5.5))
    plt.axhline(1.0, color="black", linewidth=1)
    plt.plot(x, num_bf_over_nom, marker="o", linewidth=2, label="CCnumu BF / nominal")
    plt.plot(x, den_bf_over_nom, marker="s", linewidth=2, label="CCnue BF / nominal")
    plt.plot(x, ratio_bf_over_nom, marker="^", linewidth=2, label="ratio BF / nominal")
    plt.xlabel("Local bin")
    plt.ylabel("BF / nominal")
    plt.title(mode + ": fractional source of ratio distortion")
    plt.grid(True, alpha=0.35)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, mode + "_bf_over_nominal_components.png"), dpi=200)
    plt.close()

    # 3. Fractional change.
    plt.figure(figsize=(8, 5.5))
    plt.axhline(0.0, color="black", linewidth=1)
    plt.plot(x, num_bf_over_nom - 1.0, marker="o", linewidth=2, label="CCnumu fractional change")
    plt.plot(x, den_bf_over_nom - 1.0, marker="s", linewidth=2, label="CCnue fractional change")
    plt.plot(x, ratio_bf_over_nom - 1.0, marker="^", linewidth=2, label="ratio fractional change")
    plt.xlabel("Local bin")
    plt.ylabel("(BF - nominal) / nominal")
    plt.title(mode + ": numerator, denominator, and ratio fractional changes")
    plt.grid(True, alpha=0.35)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, mode + "_fractional_changes.png"), dpi=200)
    plt.close()

    # 4. Ratio itself.
    plt.figure(figsize=(8, 5.5))
    plt.plot(x, ratio_nom, marker="o", linewidth=2, label="nominal CCnumu / CCnue")
    plt.plot(x, ratio_bf, marker="^", linewidth=1.7, linestyle=":", label="raw BF CCnumu / CCnue")
    plt.xlabel("Local bin")
    plt.ylabel("CCnumu / CCnue")
    plt.title(mode + ": constructed ratio")
    plt.grid(True, alpha=0.35)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, mode + "_constructed_ratio.png"), dpi=200)
    plt.close()


def add_rows(rows, mode, x, num_nom, den_nom, ratio_nom, num_bf, den_bf, ratio_bf):
    num_bf_over_nom = safe_div(num_bf, num_nom)
    den_bf_over_nom = safe_div(den_bf, den_nom)
    ratio_bf_over_nom = safe_div(ratio_bf, ratio_nom)

    for i in range(len(x)):
        rows.append({
            "mode": mode,
            "local_bin": int(x[i]),
            "numu_nominal": float(num_nom[i]),
            "nue_nominal": float(den_nom[i]),
            "ratio_nominal": float(ratio_nom[i]),
            "numu_bf": float(num_bf[i]),
            "nue_bf": float(den_bf[i]),
            "ratio_bf": float(ratio_bf[i]),
            "numu_bf_over_nominal": float(num_bf_over_nom[i]) if np.isfinite(num_bf_over_nom[i]) else np.nan,
            "nue_bf_over_nominal": float(den_bf_over_nom[i]) if np.isfinite(den_bf_over_nom[i]) else np.nan,
            "ratio_bf_over_nominal": float(ratio_bf_over_nom[i]) if np.isfinite(ratio_bf_over_nom[i]) else np.nan,
            "numu_frac_change": float(num_bf_over_nom[i] - 1.0) if np.isfinite(num_bf_over_nom[i]) else np.nan,
            "nue_frac_change": float(den_bf_over_nom[i] - 1.0) if np.isfinite(den_bf_over_nom[i]) else np.nan,
            "ratio_frac_change": float(ratio_bf_over_nom[i] - 1.0) if np.isfinite(ratio_bf_over_nom[i]) else np.nan,
        })


def add_summary(summary_rows, mode, num_nom, den_nom, ratio_nom, num_bf, den_bf, ratio_bf):
    num_bf_over_nom = safe_div(num_bf, num_nom)
    den_bf_over_nom = safe_div(den_bf, den_nom)
    ratio_bf_over_nom = safe_div(ratio_bf, ratio_nom)

    summary_rows.append({
        "mode": mode,
        "n_bins": len(num_nom),

        "numu_sum_nominal": float(np.sum(num_nom)),
        "numu_sum_bf": float(np.sum(num_bf)),
        "numu_sum_bf_over_nominal": float(np.sum(num_bf) / np.sum(num_nom)) if np.sum(num_nom) != 0 else np.nan,
        "numu_max_abs_frac_change": max_abs_finite(num_bf_over_nom - 1.0),
        "numu_rms_frac_change": rms_finite(num_bf_over_nom - 1.0),

        "nue_sum_nominal": float(np.sum(den_nom)),
        "nue_sum_bf": float(np.sum(den_bf)),
        "nue_sum_bf_over_nominal": float(np.sum(den_bf) / np.sum(den_nom)) if np.sum(den_nom) != 0 else np.nan,
        "nue_max_abs_frac_change": max_abs_finite(den_bf_over_nom - 1.0),
        "nue_rms_frac_change": rms_finite(den_bf_over_nom - 1.0),

        "ratio_sum_nominal": float(np.sum(ratio_nom)),
        "ratio_sum_bf": float(np.sum(ratio_bf)),
        "ratio_sum_bf_over_nominal": float(np.sum(ratio_bf) / np.sum(ratio_nom)) if np.sum(ratio_nom) != 0 else np.nan,
        "ratio_max_abs_frac_change": max_abs_finite(ratio_bf_over_nom - 1.0),
        "ratio_rms_frac_change": rms_finite(ratio_bf_over_nom - 1.0),
    })


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

    shutil.copyfile(hist_config, "HIST_CONFIG.json")

    cfg = load_hist_config(hist_config)

    required = [
        "fhc_numu_selection",
        "fhc_nue_selection",
        "rhc_numu_selection",
        "rhc_nue_selection",
    ]

    for sample in required:
        if sample not in cfg:
            raise RuntimeError(
                "Required sample {} not found in {}. "
                "Use --hist-config-tag prodNueel_noRatio or another config with absolute nue and numu samples.".format(
                    sample,
                    hist_config,
                )
            )

    mkdir(args.outdir)

    print("\n===== check_ratio_components setup =====")
    print("plot_tag    =", plot_tag)
    print("root_file   =", root_file)
    print("hist_config =", hist_config)
    print("outdir      =", args.outdir)
    print("BF dm2      =", args.dm2)
    print("BF ue4      =", args.ue4)
    print("BF umu4     =", args.umu4)
    print("BF utau4    =", args.utau4)

    nominal = load_prediction(root_file, osc=None)
    bf_raw = load_prediction(root_file, osc=(args.dm2, args.ue4, args.umu4, args.utau4))

    rows = []
    summary_rows = []

    modes = [
        ("fhc", "fhc_numu_selection", "fhc_nue_selection"),
        ("rhc", "rhc_numu_selection", "rhc_nue_selection"),
    ]

    for mode, num_sample, den_sample in modes:
        num_nom = slice_sample(nominal, cfg, num_sample)
        den_nom = slice_sample(nominal, cfg, den_sample)

        num_bf = slice_sample(bf_raw, cfg, num_sample)
        den_bf = slice_sample(bf_raw, cfg, den_sample)

        ratio_nom = safe_div(num_nom, den_nom)
        ratio_bf = safe_div(num_bf, den_bf)

        x = np.arange(1, len(num_nom) + 1)

        outdir_mode = os.path.join(args.outdir, mode)
        plot_ratio_component(
            mode,
            x,
            num_nom,
            den_nom,
            ratio_nom,
            num_bf,
            den_bf,
            ratio_bf,
            outdir_mode,
        )

        add_rows(rows, mode, x, num_nom, den_nom, ratio_nom, num_bf, den_bf, ratio_bf)
        add_summary(summary_rows, mode, num_nom, den_nom, ratio_nom, num_bf, den_bf, ratio_bf)

    write_csv(rows, os.path.join(args.outdir, "per_bin_ratio_component_comparison.csv"))
    write_csv(summary_rows, os.path.join(args.outdir, "ratio_component_summary.csv"))

    print("\n===== ratio component summary =====")
    for r in summary_rows:
        print(
            "{mode:4s}  "
            "numu sum BF/nom={numu_sum_bf_over_nominal:.4g}, max|frac|={numu_max_abs_frac_change:.4g}  "
            "nue sum BF/nom={nue_sum_bf_over_nominal:.4g}, max|frac|={nue_max_abs_frac_change:.4g}  "
            "ratio sum BF/nom={ratio_sum_bf_over_nominal:.4g}, max|frac|={ratio_max_abs_frac_change:.4g}".format(**r)
        )

    print("\nWrote:")
    print("  ", os.path.join(args.outdir, "per_bin_ratio_component_comparison.csv"))
    print("  ", os.path.join(args.outdir, "ratio_component_summary.csv"))
    print("  plots under:", args.outdir)

    print("\nInterpretation:")
    print("  If CCnumu and CCnue fractional changes are both smooth but their ratio wiggles, the division is amplifying differences.")
    print("  If CCnue has much larger low-E fractional changes, the denominator is driving the ratio wiggles.")
    print("  If CCnumu has larger low-E fractional changes, the numerator is driving the ratio wiggles.")


if __name__ == "__main__":
    main()

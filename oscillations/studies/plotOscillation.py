#!/usr/bin/env python3

import os
import sys
import shutil
import argparse

import ROOT
import PlotUtils

from config.AnalysisConfig import AnalysisConfig
from tools.StitchedHistogram import *
from tools.PlotHistogram import *


# ============================================================
# ARGUMENTS
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot stitched distributions for a fixed oscillation hypothesis."
    )

    parser.add_argument(
        "--hist-config-tag",
        required=True,
        help="Histogram configuration tag, e.g. predNueel_p8",
    )

    parser.add_argument(
        "--dm2",
        required=True,
        type=float,
        help="Delta m^2_41 [eV^2]",
    )

    parser.add_argument(
        "--ue4",
        required=True,
        type=float,
        help="|U_e4|^2",
    )

    parser.add_argument(
        "--umu4",
        required=True,
        type=float,
        help="|U_mu4|^2",
    )

    parser.add_argument(
        "--utau4",
        default=0.0,
        type=float,
        help="|U_tau4|^2 (default: 0)",
    )

    parser.add_argument(
        "--exclude",
        default=None,
        help="Samples to exclude",
    )

    parser.add_argument(
        "--lam",
        default=None,
        type=float,
        help="Flux profiling regularization lambda",
    )

    parser.add_argument(
        "--profile-flux",
        action="store_true",
        help="Use flux profiling in plotting",
    )

    parser.add_argument(
        "--no-profile-flux",
        action="store_true",
        help="Disable flux profiling in plotting",
    )

    parser.add_argument(
        "--profile-only",
        default=None,
    )

    parser.add_argument(
        "--profile-n-universes",
        default=None,
        type=int,
    )

    parser.add_argument(
        "--output-tag",
        default=None,
        help="Optional output tag",
    )

    return parser.parse_args()


# ============================================================
# ROOT / PLOT SETTINGS
# ============================================================

ROOT.gROOT.SetBatch(True)

MNVPLOTTER = PlotUtils.MnvPlotter()
MNVPLOTTER.draw_normalized_to_bin_width = False
MNVPLOTTER.legend_text_size = 0.03
MNVPLOTTER.mc_bkgd_color = 46
MNVPLOTTER.mc_bkgd_line_color = 46
MNVPLOTTER.legend_n_columns = 3
MNVPLOTTER.data_bkgd_color = 12
MNVPLOTTER.data_bkgd_style = 24
MNVPLOTTER.axis_maximum = 500
MNVPLOTTER.height_nspaces_per_hist = 1.2
MNVPLOTTER.width_xspace_per_letter = 0.4
MNVPLOTTER.legend_offset_x = 0.15

ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    args = parse_args()

    ccnueroot = os.environ.get("CCNUEROOT")
    if ccnueroot is None:
        raise RuntimeError("CCNUEROOT is not set.")

    # --------------------------------------------------------
    # Histogram configuration
    # --------------------------------------------------------

    plot_tag = args.hist_config_tag

    filename = "rootfiles/NuE_stitched_hists_{}.root".format(plot_tag)
    file_path = "{}/oscillations/{}".format(ccnueroot, filename)

    hist_config = "HIST_CONFIG_{}.json".format(plot_tag)

    if not os.path.exists(hist_config):
        raise RuntimeError(
            "Missing requested hist config file: {}".format(hist_config)
        )

    shutil.copyfile(hist_config, "HIST_CONFIG.json")

    # --------------------------------------------------------
    # Plot/profiling configuration
    # --------------------------------------------------------

    lambda_value = args.lam
    if lambda_value is None:
        lambda_value = getattr(AnalysisConfig, "lambdaValue", 1)

    exclude_samples = args.exclude
    if exclude_samples is None:
        exclude_samples = getattr(AnalysisConfig, "exclude", [])

    if exclude_samples is None:
        exclude_samples = ""
    elif isinstance(exclude_samples, str):
        if exclude_samples.strip().lower() in ["none", ""]:
            exclude_samples = ""

    profileFlux = getattr(AnalysisConfig, "profileFlux", False)

    if args.profile_flux:
        profileFlux = True

    if args.no_profile_flux:
        profileFlux = False

    profile_only = args.profile_only
    profile_n_universes = args.profile_n_universes

    mask_spec = None

    # --------------------------------------------------------
    # Oscillation hypothesis
    # --------------------------------------------------------

    res = {
        "m": args.dm2,
        "ue4": args.ue4,
        "umu4": args.umu4,
        "utau4": args.utau4,
    }

    print("\n===== fixed oscillation hypothesis =====")
    print("file         =", file_path)
    print("hist config  =", hist_config)
    print("dm2          =", res["m"], "eV^2")
    print("|Ue4|^2      =", res["ue4"])
    print("|Umu4|^2     =", res["umu4"])
    print("|Utau4|^2    =", res["utau4"])
    print("profileFlux  =", profileFlux)
    print("lambda       =", lambda_value)
    print("exclude      =", exclude_samples)
    print("profile_only =", profile_only)
    print("N universes  =", profile_n_universes)
    print("========================================\n")

    # --------------------------------------------------------
    # Load stitched histogram
    # --------------------------------------------------------

    sample_histogram = StitchedHistogram("sample")
    sample_histogram.SetHistConfig(hist_config)
    sample_histogram.Load(file_path)

    CheckElasticPredictionInput(
        sample_histogram,
        hist_config=hist_config,
        verbose=True,
    )

    # --------------------------------------------------------
    # Apply requested oscillation hypothesis
    # --------------------------------------------------------

    sample_histogram.OscillateHistogram(
        res["m"],
        res["ue4"],
        res["umu4"],
        res["utau4"],
    )

    sample_histogram.SetPlottingStyle()

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------

    if args.output_tag is None:
        output_tag = "{}_dm2_{:.3f}_ue4_{:.4f}_umu4_{:.4f}_utau4_{:.4f}".format(
            plot_tag,
            res["m"],
            res["ue4"],
            res["umu4"],
            res["utau4"],
        )
    else:
        output_tag = args.output_tag

    plotter = PlottingContainer(
        "oscillation_hypothesis_{}".format(output_tag),
        sample_histogram,
        hist_config=hist_config,
    )

    plotter.SetExclude(exclude_samples)
    plotter.SetLambda(lambda_value)
    plotter.SetMaskSpec(mask_spec)
    plotter.SetProfileOnly(profile_only)
    plotter.SetProfileNUniverses(profile_n_universes)

    plotter.PlotOscillationEffects(
        res,
        AnalysisConfig.ntuple_tag,
        useMarg=profileFlux,
        plotSamples=False,
        plot_tag=output_tag,
    )

    print("\nFinished plotting oscillation hypothesis.")
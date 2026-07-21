import os
import logging, sys
import copy
import ROOT
import PlotUtils
import numpy as np
import shutil
import argparse

def parse_fitdata_args():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--hist-config-tag", default=None)
    parser.add_argument("--exclude", default=None)
    parser.add_argument("--lam", default=None, type=float)
    parser.add_argument("--profile-flux", action="store_true")
    parser.add_argument("--no-profile-flux", action="store_true")
    parser.add_argument("--profile-only", default=None)
    parser.add_argument("--profile-n-universes", default=None, type=int)

    args, remaining = parser.parse_known_args()
    sys.argv = [sys.argv[0]] + remaining
    return args

fit_args = parse_fitdata_args()

np.set_printoptions(precision=4)
np.set_printoptions(linewidth=1520)
np.set_printoptions(threshold=sys.maxsize)

from scipy import optimize, integrate
from config.AnalysisConfig import AnalysisConfig
import shutil
import argparse
ccnueroot = os.environ.get('CCNUEROOT')

import math
# import psutil
import time
from array import array

from tools.PlotLibrary import HistHolder
from tools.StitchedHistogram import *
from tools.PlotHistogram import *
from tools.Fitters import *

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

MNVPLOTTER = PlotUtils.MnvPlotter()
MNVPLOTTER.draw_normalized_to_bin_width=False
MNVPLOTTER.legend_text_size = 0.04
MNVPLOTTER.mc_bkgd_color = 46
MNVPLOTTER.mc_bkgd_line_color = 46
MNVPLOTTER.legend_n_columns = 3
MNVPLOTTER.data_bkgd_color = 12  # gray
MNVPLOTTER.data_bkgd_style = 24  # circle
MNVPLOTTER.axis_maximum = 500

MNVPLOTTER.height_nspaces_per_hist = 1.2
MNVPLOTTER.width_xspace_per_letter = .4
MNVPLOTTER.legend_text_size        = .03
MNVPLOTTER.legend_offset_x         = .15

ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)


if __name__ == "__main__":

    plot_tag = fit_args.hist_config_tag
    if plot_tag is None:
        plot_tag = getattr(AnalysisConfig, "hist_config_tag", "default")

    if plot_tag in [None, "", "none"]:
        plot_tag = "default"

    filename = "rootfiles/NuE_stitched_hists_{}.root".format(plot_tag)
    file_path = "{}/oscillations/{}".format(ccnueroot, filename)

    hist_config = "HIST_CONFIG_{}.json".format(plot_tag)

    if not os.path.exists(hist_config):
        raise RuntimeError("Missing requested hist config file: {}".format(hist_config))

    shutil.copyfile(hist_config, "HIST_CONFIG.json")

    print("plot_tag    =", plot_tag)
    print("file        =", file_path)
    print("hist_config =", hist_config)
    print("Copied {} -> HIST_CONFIG.json".format(hist_config))

    lambda_value = fit_args.lam
    if lambda_value is None:
        lambda_value = getattr(AnalysisConfig, "lambdaValue", 1)

    exclude_samples = fit_args.exclude
    if exclude_samples is None:
        exclude_samples = getattr(AnalysisConfig, "exclude", [])

    profileFlux = getattr(AnalysisConfig, "profileFlux", False)
    if fit_args.profile_flux:
        profileFlux = True
    if fit_args.no_profile_flux:
        profileFlux = False

    profile_only = fit_args.profile_only
    profile_n_universes = fit_args.profile_n_universes

    if exclude_samples is None:
        exclude_samples = ""
    elif isinstance(exclude_samples, str):
        if exclude_samples.strip().lower() in ["none", ""]:
            exclude_samples = ""

    print("\n===== fitData setup =====")
    print("file        =", file_path)
    print("profileFlux =", profileFlux)
    print("lambda      =", lambda_value)
    print("exclude     =", exclude_samples)
    print("profile_only =", profile_only)
    print("profile_n_universes =", profile_n_universes)


    sample_histogram = StitchedHistogram("sample")
    sample_histogram.Load(file_path)

    # mc = np.array(sample_histogram.GetMCHistogram())[1:-1]

    # sample_histogram.OscillateHistogram(
    #     40.0,   # arbitrary dm2
    #     0.0,
    #     0.0,
    #     0.0,
    # )

    # osc_null = np.array(
    #     sample_histogram.GetOscillatedHistogram()
    # )[1:-1]

    # diff = osc_null - mc
    # tol = 1e-8

    # if np.max(np.abs(diff)) > tol:
    #     raise RuntimeError(
    #         "Zero-mixing closure failed: max abs diff = {}".format(
    #             np.max(np.abs(diff))
    #         )
    #     )

    # print("\n===== ZERO-MIXING CLOSURE =====")
    # print("max abs diff =", np.max(np.abs(diff)))
    # print("norm diff    =", np.linalg.norm(diff))
    # print("sum diff     =", np.sum(diff))

    # for i, d in enumerate(diff):
    #     if abs(d) > 1e-8:
    #         print(
    #             "bin {:3d}: mc={:14.8f} osc0={:14.8f} diff={:14.8e}".format(
    #                 i,
    #                 mc[i],
    #                 osc_null[i],
    #                 d,
    #             )
    #         )
    # print("===== END ZERO-MIXING CLOSURE =====\n")

    # # Reload fresh histogram so the real fit starts from a clean nominal state.
    # sample_histogram = StitchedHistogram("sample")
    # sample_histogram.Load(file_path)


    # mc = np.array(sample_histogram.GetMCHistogram())[1:-1]

    # sample_histogram.OscillateHistogram(
    #     0.0,
    #     0.15,
    #     0.02,
    #     0.20,
    # )

    # osc_dm2zero = np.array(
    #     sample_histogram.GetOscillatedHistogram()
    # )[1:-1]

    # diff_dm2zero = osc_dm2zero - mc
    # if np.max(np.abs(diff_dm2zero)) > tol:
    #     raise RuntimeError(
    #         "dm2=0 closure failed: max abs diff = {}".format(
    #             np.max(np.abs(diff_dm2zero))
    #         )
    #     )

    # print("\n===== DM2-ZERO CLOSURE =====")
    # print("max abs diff =", np.max(np.abs(diff_dm2zero)))
    # print("norm diff    =", np.linalg.norm(diff_dm2zero))
    # print("sum diff     =", np.sum(diff_dm2zero))

    # for i, d in enumerate(diff_dm2zero):
    #     if abs(d) > 1e-8:
    #         print(
    #             "bin {:3d}: mc={:14.8f} osc_dm2zero={:14.8f} diff={:14.8e}".format(
    #                 i,
    #                 mc[i],
    #                 osc_dm2zero[i],
    #                 d,
    #             )
    #         )

    # print("===== END DM2-ZERO CLOSURE =====\n")

    # # Reload fresh histogram so the real fit starts from a clean nominal state.
    # sample_histogram = StitchedHistogram("sample")
    # sample_histogram.Load(file_path)

    # # # Direct CCnue configuration
    # # dm2_test = 40.749
    # # umu4_test = 0.017403
    # # utau4_test = 0.232719

    # # # CC ratios, included in profiling
    # # dm2_test = 40.895
    # # umu4_test = 0.021220
    # # utau4_test = 0.023624

    # # CC ratios, excluded from profiling
    # dm2_test = 41.422
    # umu4_test = 0.017750
    # utau4_test = 0.660

    # print("\n===== UE4 BOUNDARY SCAN =====")
    # print("dm2 fixed  =", dm2_test)
    # print("umu4 fixed =", umu4_test)
    # print("utau4 fixed=", utau4_test)

    # scan_results = []

    # for ue4_test in np.linspace(0.0, 0.30, 31):
    #     sample_histogram = StitchedHistogram("sample")
    #     sample_histogram.Load(file_path)

    #     sample_histogram.OscillateHistogram(
    #         dm2_test,
    #         ue4_test,
    #         umu4_test,
    #         utau4_test,
    #     )

    #     stat_test = Statistics(
    #         sample_histogram,
    #         exclude=exclude_samples,
    #         lam=lambda_value,
    #         mask_spec=None,
    #         profile_only=profile_only,
    #         profile_n_universes=profile_n_universes,
    #     )

    #     chi2_test, penalty_test = stat_test.Chi2DataMC(
    #         marginalize=profileFlux,
    #         useOsc=True,
    #     )

    #     residual_test = chi2_test - penalty_test
    #     scan_results.append(
    #         (ue4_test, chi2_test, residual_test, penalty_test)
    #     )

    #     print(
    #         "ue4={:.3f}  chi2={:.8f}  residual={:.8f}  penalty={:.8f}".format(
    #             ue4_test,
    #             chi2_test,
    #             residual_test,
    #             penalty_test,
    #         )
    #     )

    # best_scan_point = min(scan_results, key=lambda x: x[1])

    # print("\n===== UE4 SCAN MINIMUM =====")
    # print("ue4     =", best_scan_point[0])
    # print("chi2    =", best_scan_point[1])
    # print("residual=", best_scan_point[2])
    # print("penalty =", best_scan_point[3])


    # print("Loaded sample_mc nbins =", sample_histogram.GetMCHistogram().GetNbinsX())
    # print("Loaded sample_data nbins =", sample_histogram.GetDataHistogram().GetNbinsX())
    # print("Loaded sample_mc integral =", sample_histogram.GetMCHistogram().Integral())
    # print("Loaded sample_data integral =", sample_histogram.GetDataHistogram().Integral())

    # print("\n===== covariance check =====")
    # print("external covariances =", sample_histogram.external_covariances.keys())
    # print("covariance shape     =", sample_histogram.GetCovarianceMatrix(False).shape)
    # print("sansFlux cov shape   =", sample_histogram.GetCovarianceMatrix(True).shape)

    # if "fhc_elastic" in sample_histogram.external_covariances:
    #     print("First 6x6 full covariance block:")
    #     print(sample_histogram.GetCovarianceMatrix(False)[:6, :6])

    #     print("First 6x6 sans-flux covariance block:")
    #     print(sample_histogram.GetCovarianceMatrix(True)[:6, :6])
    # else:
    #     print("No fhc_elastic external covariance loaded. This is expected if elastic is excluded.")

    print("\n===== covariance decomposition check =====")

    hmc = sample_histogram.GetMCHistogram()
    nbins = hmc.GetNbinsX()

    cov_full = sample_histogram.GetCovarianceMatrix(False)
    cov_sans = sample_histogram.GetCovarianceMatrix(True)

    print("nbins                  =", nbins)
    print("cov_full shape         =", cov_full.shape)
    print("cov_sansFlux shape     =", cov_sans.shape)

    flux_cov = cov_full - cov_sans

    print("\nTrace checks:")
    print("trace full             =", np.trace(cov_full))
    print("trace sansFlux         =", np.trace(cov_sans))
    print("trace flux difference  =", np.trace(flux_cov))

    # print("\nDiagonal checks by bin:")
    # for i in range(nbins):
    #     print(
    #         "bin {:3d}: full={:12.5e}  sans={:12.5e}  fluxdiff={:12.5e}  mc={:12.5e}".format(
    #             i+1,
    #             cov_full[i, i],
    #             cov_sans[i, i],
    #             flux_cov[i, i],
    #             hmc.GetBinContent(i+1),
    #         )
    #     )


    mask_spec = None
    # mask_spec = {
    #     "fhc_ratio": [1, 9, 10],
    #     # "rhc_ratio": [1, 9, 10],
    # }

    # stat = Statistics(sample_histogram, lam=lambda_value, exclude=exclude_samples)
    stat = Statistics(
        sample_histogram,
        exclude=exclude_samples,
        lam=lambda_value,
        mask_spec=mask_spec,
        profile_only=profile_only,
        profile_n_universes=profile_n_universes,
    )


    chi2_null, penalty = stat.Chi2DataMC(marginalize=profileFlux)
    print("null chi2: {:.3f}  penalty: {:.3f}".format(chi2_null, penalty))

    fitter = OscillationFitter(
        sample_histogram,
        marginalize_flux=profileFlux,
        exclude=exclude_samples,
        lam=lambda_value,
        mask_spec=mask_spec,
        profile_only=profile_only,
        profile_n_universes=profile_n_universes,
    )


    chi2_fit, res = fitter.DoFit()

    print("Data fit: delta chi2 = {:.3f} = {:.3f} - {:.3f}".format(
        chi2_null - chi2_fit, chi2_null, chi2_fit
    ))
    print("Best fit params:")
    print("   delta m^2 = {:.3f} eV^2 +- {:.4f}".format(res["m"], 0))
    print("   U_e4^2    = {:.3f}      +- {:.4f}".format(res["ue4"], 0))
    print("   U_mu4^2   = {:.5f}    +- {:.4f}".format(res["umu4"], 0))
    print("   U_tau4^2  = {:.3f}      +- {:.4f}".format(res["utau4"], 0))

    sample_histogram.OscillateHistogram(res["m"], res["ue4"], res["umu4"], res["utau4"])
    sample_histogram.SetPlottingStyle()

    # plot_tag = os.path.basename(filename)
    # plot_tag = plot_tag.replace("NuE_stitched_hists_", "")
    # plot_tag = plot_tag.replace(".root", "")

    # print("plot_tag =", plot_tag)

    plotter = PlottingContainer("fitted_histogram_{}".format(plot_tag), sample_histogram)
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
        plot_tag=plot_tag
    )

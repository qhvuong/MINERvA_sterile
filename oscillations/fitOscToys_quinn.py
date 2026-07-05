import os
import logging, sys
import ROOT
import PlotUtils
import numpy as np
from scipy import optimize, integrate
from scipy.stats import multivariate_normal
import argparse
import shutil
import math
from array import array

def parse_fitasimov_args():
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

fit_args = parse_fitasimov_args()

ccnueroot = os.environ.get("CCNUEROOT")
process = os.environ.get("PROCESS")

from tools.PlotLibrary import HistHolder
from tools.Fitters import *
from tools.StitchedHistogram import *
from config.AnalysisConfig import AnalysisConfig

seed = int(os.environ.get("RANDOM_SEED", "0"))
if seed != 0:
    np.random.seed(seed)
    print("Using RANDOM_SEED =", seed)

# Reuse helpers from the validated Asimov/null-toy script.
from fitAsimovs_quinn import (
    print_cov_diagnostics,
    throw_multivariate_psd,
    ThrowPoissons,
)

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

ccnueroot = os.environ.get("CCNUEROOT")
process = os.environ.get("PROCESS")

def make_safe_tag(x):
    if x is None:
        return "none"
    x = str(x).strip()
    if x == "":
        return "none"
    return x.replace(",", "-").replace("/", "-")

def get_true_osc_params():
    true_dm2 = float(os.environ.get(
        "TRUE_DM2",
        getattr(AnalysisConfig, "single_dm2", 0.0) or 0.0
    ))

    true_ue4 = float(os.environ.get(
        "TRUE_UE4",
        getattr(AnalysisConfig, "single_ue4", 0.0)
    ))

    true_umu4 = float(os.environ.get(
        "TRUE_UMU4",
        getattr(AnalysisConfig, "single_umu4", 0.0)
    ))

    true_utau4 = float(os.environ.get(
        "TRUE_UTAU4",
        getattr(AnalysisConfig, "single_utau4", 0.0)
    ))

    return true_dm2, true_ue4, true_umu4, true_utau4


def GetOscillatedMeanVector(histogram, dm2, ue4, umu4, utau4):
    """
    Build the true expected mean vector for toys by oscillating the nominal MC.

    The returned vector excludes underflow/overflow and has length equal to
    the stitched global-bin count.
    """

    print("\n===== Building oscillated toy mean =====")
    print("true dm2   =", dm2)
    print("true ue4   =", ue4)
    print("true umu4  =", umu4)
    print("true utau4 =", utau4)

    histogram.OscillateHistogram(
        dm2,
        ue4,
        umu4,
        utau4,
        False,
        False,
    )

    osc_hist = histogram.GetOscillatedHistogram()
    osc_mean = np.array(osc_hist)[1:-1]

    cv_mean = np.array(histogram.GetMCHistogram())[1:-1]

    print("mean vector bins      =", len(osc_mean))
    print("nominal MC sum        =", np.sum(cv_mean))
    print("oscillated MC sum     =", np.sum(osc_mean))
    print("osc/nominal sum ratio =", np.sum(osc_mean) / np.sum(cv_mean) if np.sum(cv_mean) != 0 else np.nan)
    print("max abs osc-CV shift  =", np.max(np.abs(osc_mean - cv_mean)))

    return osc_mean


def ThrowSystematicsAroundMean(
    histogram,
    pred_vals,
    throwFlux=False,
    n_samples=50,
    doDiagnostics=False,
):
    """
    Same strategy as fitAsimovs_quinn.ThrowSystematics(), except the central
    mean is supplied explicitly.

    For throwFlux=True:
      oscillated mean
      + independent flux systematic throw
      + independent non-flux systematic throw
      + Poisson later

    For throwFlux=False:
      oscillated mean
      + full systematics-only Gaussian throw
      + Poisson later
    """

    pred_vals = np.asarray(pred_vals, dtype=float)

    # Full covariance currently includes stat + all systematics.
    V_full = histogram.GetCovarianceMatrix(False)
    V_sansFlux = histogram.GetCovarianceMatrix(sansFlux=True)
    V_flux = V_full - V_sansFlux

    # Same residual-hist stat convention used by StitchedHistogram.SetCovarianceMatrices().
    h_test = histogram.GetDataHistogram().Clone()
    h_test.Add(histogram.GetMCHistogram(), -1)

    V_stat = TMatrix_to_Numpy(h_test.GetStatErrorMatrix())[1:-1, 1:-1]

    # Remove stat from Gaussian throw covariance.
    # The final Poisson throw handles counting statistics.
    V_full_systOnly = V_full - V_stat
    V_sansFlux_systOnly = V_sansFlux - V_stat

    # Force symmetry.
    V_flux = 0.5 * (V_flux + V_flux.T)
    V_full_systOnly = 0.5 * (V_full_systOnly + V_full_systOnly.T)
    V_sansFlux_systOnly = 0.5 * (V_sansFlux_systOnly + V_sansFlux_systOnly.T)

    if doDiagnostics:
        print_cov_diagnostics("V_full = stat + all syst", V_full)
        print_cov_diagnostics("V_sansFlux = stat + nonflux syst", V_sansFlux)
        print_cov_diagnostics("V_stat", V_stat)
        print_cov_diagnostics("V_flux", V_flux)
        print_cov_diagnostics("V_full_systOnly", V_full_systOnly)
        print_cov_diagnostics("V_sansFlux_systOnly", V_sansFlux_systOnly)

    if throwFlux:
        # Independent flux throw per toy.
        flux_throws = throw_multivariate_psd(
            pred_vals,
            V_flux,
            size=n_samples,
        )

        flux_throws = np.asarray(flux_throws)
        if flux_throws.ndim == 1:
            flux_throws = flux_throws.reshape(1, -1)

        # Independent non-flux systematic throw around each flux-shifted toy mean.
        sys_throws = []
        for flux_mean in flux_throws:
            toy_mean = throw_multivariate_psd(
                flux_mean,
                V_sansFlux_systOnly,
                size=1,
            )
            toy_mean = np.asarray(toy_mean).reshape(-1)
            sys_throws.append(toy_mean)

        sys_throws = np.asarray(sys_throws)

    else:
        sys_throws = throw_multivariate_psd(
            pred_vals,
            V_full_systOnly,
            size=n_samples,
        )

    sys_throws = np.asarray(sys_throws)
    if sys_throws.ndim == 1:
        sys_throws = sys_throws.reshape(1, -1)

    print("ThrowSystematicsAroundMean shape:", sys_throws.shape)
    return sys_throws


def RunNominalFit(
    histogram,
    profileFlux=True,
    lam=1,
    exclude="",
    mask_spec=None,
    profile_only=None,
    profile_n_universes=None,
):
    print("\n===== Nominal no-throw fit =====")
    print("profileFlux =", profileFlux)
    print("lambda      =", lam)
    print("exclude     =", exclude)
    print("profile_only =", profile_only)
    print("profile_n_universes =", profile_n_universes)

    stat = Statistics(
        histogram,
        lam=lam,
        exclude=exclude,
        mask_spec=mask_spec,
        profile_only=profile_only,
        profile_n_universes=profile_n_universes,
    )

    chi2_null, null_penalty = stat.Chi2DataMC(
        marginalize=profileFlux
    )

    fitter = OscillationFitter(
        histogram,
        lam=lam,
        exclude=exclude,
        marginalize_flux=profileFlux,
        mask_spec=mask_spec,
        profile_only=profile_only,
        profile_n_universes=profile_n_universes,
    )

    chi2_fit, res = fitter.DoFit()

    dchi2 = chi2_null - chi2_fit

    print("\nNominal fit result:")
    print("  chi2_null = {:.6f}".format(chi2_null))
    print("  null penalty = {:.6f}".format(null_penalty))
    print("  chi2_fit  = {:.6f}".format(chi2_fit))
    print("  delta chi2 = {:.6f}".format(dchi2))
    print("  best-fit dm2   = {:.6f}".format(res["m"]))
    print("  best-fit ue4   = {:.6f}".format(res["ue4"]))
    print("  best-fit umu4  = {:.6f}".format(res["umu4"]))
    print("  best-fit utau4 = {:.6f}".format(res["utau4"]))

    return dchi2, chi2_null, chi2_fit, res

def FitToyExperimentsAtTruth(
    histogram,
    experiments,
    true_dm2,
    true_ue4,
    true_umu4,
    true_utau4,
    profileFlux=True,
    lam=1,
    exclude="",
    mask_spec=None,
    profile_only=None,
    profile_n_universes=None,
):
    stitched_data = histogram.GetDataHistogram()

    experiments = np.asarray(experiments)
    if experiments.ndim == 1:
        experiments = experiments.reshape(1, -1)

    print("FitToyExperimentsAtTruth experiments shape:", experiments.shape)

    dchi2s = []

    for itoy, toy in enumerate(experiments):
        print("\n===== Fitting oscillated toy {} / {} =====".format(itoy + 1, len(experiments)))

        # Build pseudo-data histogram exactly like fitAsimovs_quinn.FitToyExperiments().
        weights = stitched_data.Clone().GetCVHistoWithStatError()

        for i in range(1, weights.GetNbinsX() + 1):
            weight = stitched_data.GetBinContent(i) / toy[i - 1] if toy[i - 1] != 0 else stitched_data.GetBinContent(i)
            weights.SetBinContent(i, weight)
            weights.SetBinError(i, 0)

        data_histogram = stitched_data.Clone()
        data_histogram.DivideSingle(data_histogram, weights)
        histogram.SetDataHistogram(data_histogram)

        stat = Statistics(
            histogram,
            lam=lam,
            exclude=exclude,
            mask_spec=mask_spec,
            profile_only=profile_only,
            profile_n_universes=profile_n_universes,
        )

        # Chi2 at the truth/BF point.
        histogram.OscillateHistogram(true_dm2, true_ue4, true_umu4, true_utau4)

        chi2_truth, penalty_truth = stat.Chi2DataMC(
            marginalize=profileFlux,
            useOsc=True,
            usePseudo=False,
        )

        # Free fit of the same toy.
        fitter = OscillationFitter(
            histogram,
            lam=lam,
            exclude=exclude,
            marginalize_flux=profileFlux,
            mask_spec=mask_spec,
            profile_only=profile_only,
            profile_n_universes=profile_n_universes,
        )

        chi2_best, res = fitter.DoFit()

        dchi2 = chi2_truth - chi2_best
        dchi2s.append(dchi2)

        print("truth chi2 = {:.6f}, penalty = {:.6f}".format(chi2_truth, penalty_truth))
        print("best chi2  = {:.6f}".format(chi2_best))
        print("toy dchi2  = {:.6f}".format(dchi2))
        print(
            "toy BF     = dm2 {:.6f}, ue4 {:.6f}, umu4 {:.6f}, utau4 {:.6f}".format(
                res["m"], res["ue4"], res["umu4"], res["utau4"]
            )
        )

        # Reset original data before next toy.
        histogram.SetDataHistogram(stitched_data)

    histogram.SetDataHistogram(stitched_data)
    return np.asarray(dchi2s)


if __name__ == "__main__":
    hist_config_tag = fit_args.hist_config_tag
    if hist_config_tag is None:
        hist_config_tag = getattr(AnalysisConfig, "hist_config_tag", "default")

    if hist_config_tag in [None, "", "none"]:
        hist_config_tag = "default"


    filename = "rootfiles/NuE_stitched_hists_{}.root".format(hist_config_tag)
    file_path = "{}/oscillations/{}".format(ccnueroot, filename)

    hist_config = "HIST_CONFIG_{}.json".format(hist_config_tag)
    if not os.path.exists(hist_config):
        raise RuntimeError("Missing requested hist config file: {}".format(hist_config))

    shutil.copyfile(hist_config, "HIST_CONFIG.json")

    print("hist_config_tag =", hist_config_tag)
    print("Using hist config:", hist_config)
    print("Copied {} -> HIST_CONFIG.json".format(hist_config))

    sample_histogram = StitchedHistogram("sample")
    sample_histogram.Load(file_path)

    if fit_args.profile_flux:
        profileFlux = True
    elif fit_args.no_profile_flux:
        profileFlux = False
    else:
        profileFlux = getattr(AnalysisConfig, "profileFlux", True)

    lam = fit_args.lam
    if lam is None:
        lam = getattr(AnalysisConfig, "lambdaValue", 1)

    exclude = fit_args.exclude
    if exclude is None:
        exclude = getattr(AnalysisConfig, "exclude", "")

    if exclude is None:
        exclude = ""
    elif isinstance(exclude, str) and exclude.strip().lower() in ["none", ""]:
        exclude = ""

    profile_only = fit_args.profile_only
    if profile_only is not None and profile_only.strip().lower() in ["", "none"]:
        profile_only = None

    profile_n_universes = fit_args.profile_n_universes

    mask_spec = None

    profile_tag = "profiledFlux" if profileFlux else "noFluxProfile"

    if exclude in [None, "", "none"]:
        exclude_tag = "includeAll"
    else:
        exclude_tag = "exclude{}".format(make_safe_tag(exclude))

    if profile_n_universes is None:
        nprof_tag = "NprofAll"
    else:
        nprof_tag = "Nprof{}".format(profile_n_universes)

    profile_only_tag = ""
    if profile_only not in [None, "", "none"]:
        profile_only_tag = "_profileOnly{}".format(make_safe_tag(profile_only))

    analysis_tag = "{}_{}_{}{}".format(
        profile_tag,
        exclude_tag,
        nprof_tag,
        profile_only_tag,
    )

    true_dm2, true_ue4, true_umu4, true_utau4 = get_true_osc_params()

    mode = "profiledFlux" if profileFlux else "noFluxProfile"

    nominal_fit_only = os.environ.get("NOMINAL_FIT_ONLY", "0") == "1"
    print("nominal fit only =", nominal_fit_only)
    n_samples = int(os.environ.get("N_TOYS", "50"))
    cov_diagnostic_only = os.environ.get("COV_DIAGNOSTIC_ONLY", "0") == "1"
    do_cov_diagnostics = os.environ.get("DO_COV_DIAGNOSTICS", "0") == "1"


    process_tag = process if process is not None else "local"

    osc_tag = "dm2_{:.6g}_ue4_{:.6g}_umu4_{:.6g}_utau4_{:.6g}".format(
        true_dm2,
        true_ue4,
        true_umu4,
        true_utau4,
    )
    osc_tag = (
        osc_tag
        .replace(".", "p")
        .replace("-", "m")
        .replace("+", "p")
    )

    print("\n===== FitOscToys setup =====")
    print("file        =", file_path)
    print("profileFlux =", profileFlux)
    print("mode        =", mode)
    print("lambda      =", lam)
    print("exclude     =", exclude)
    print("profile_only =", profile_only)
    print("profile_n_universes =", profile_n_universes)
    print("n_samples   =", n_samples)
    print("process     =", process_tag)
    print("output_dir  =", AnalysisConfig.output_dir)
    print("true dm2    =", true_dm2)
    print("true ue4    =", true_ue4)
    print("true umu4   =", true_umu4)
    print("true utau4  =", true_utau4)
    print("cov diag only      =", cov_diagnostic_only)
    print("do cov diagnostics =", do_cov_diagnostics)

    print("\n===== Loaded stitched covariance check =====")
    print("external covariances =", sample_histogram.external_covariances.keys())
    print("full covariance shape =", sample_histogram.GetCovarianceMatrix(False).shape)
    print("sans-flux covariance shape =", sample_histogram.GetCovarianceMatrix(True).shape)

    if nominal_fit_only:
        RunNominalFit(
            sample_histogram,
            profileFlux=profileFlux,
            lam=lam,
            exclude=exclude,
            mask_spec=mask_spec,
            profile_only=profile_only,
            profile_n_universes=profile_n_universes,
        )
        print("\nDone nominal no-throw fit. Exiting before oscillated toy throws.")
        sys.exit(0)


    osc_mean = GetOscillatedMeanVector(
        sample_histogram,
        true_dm2,
        true_ue4,
        true_umu4,
        true_utau4,
    )

    sys_throws = ThrowSystematicsAroundMean(
        sample_histogram,
        pred_vals=osc_mean,
        throwFlux=profileFlux,
        n_samples=n_samples,
        doDiagnostics=do_cov_diagnostics,
    )

    if cov_diagnostic_only:
        print("\nDone oscillated covariance/throw diagnostic. Exiting before Poisson + fits.")
        sys.exit(0)

    experiments = ThrowPoissons(
        sys_throws,
        sample_histogram,
        throwFlux=profileFlux,
    )

    dchi2s = FitToyExperimentsAtTruth(
        sample_histogram,
        experiments,
        true_dm2=true_dm2,
        true_ue4=true_ue4,
        true_umu4=true_umu4,
        true_utau4=true_utau4,
        profileFlux=profileFlux,
        lam=lam,
        exclude=exclude,
        mask_spec=mask_spec,
        profile_only=profile_only,
        profile_n_universes=profile_n_universes,
    )

    outname = "{}/sample_dchi2s_{}_{}_truthOsc_{}_{}.csv".format(
        AnalysisConfig.output_dir,
        hist_config_tag,
        analysis_tag,
        osc_tag,
        process_tag,
    )

    np.savetxt(outname, dchi2s, delimiter=",")

    print("Saved:", outname)
    print("dchi2 min/max/mean =", np.nanmin(dchi2s), np.nanmax(dchi2s), np.nanmean(dchi2s))
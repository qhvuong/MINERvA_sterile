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

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

# Kevin test
# throw first from the flux covariance using cholesky decomp
# throw rest of the systematics around that result
# throw poisson from that result
# do rest as normal

def print_cov_diagnostics(name, V, max_print=8):
    V = np.asarray(V)

    print("\n===== {} diagnostics =====".format(name))
    print("shape =", V.shape)
    print("min   =", np.min(V))
    print("max   =", np.max(V))
    print("diag min/max =", np.min(np.diag(V)), np.max(np.diag(V)))
    print("symmetric max abs diff =", np.max(np.abs(V - V.T)))

    eigvals = np.linalg.eigvalsh(0.5 * (V + V.T))
    print("eig min/max =", np.min(eigvals), np.max(eigvals))
    print("num negative eigvals =", np.sum(eigvals < -1e-10))
    print("num tiny/negative eigvals =", np.sum(eigvals < 1e-12))

    print("\n{} top-left {}x{} block:".format(name, max_print, max_print))
    print(V[:max_print, :max_print])

    print("\n{} diagonal first {} bins:".format(name, max_print))
    print(np.diag(V)[:max_print])

def throw_multivariate_psd(mean, cov, size=1, tol=1e-10):
    """
    Draw from N(mean, cov) for a symmetric positive-semidefinite covariance.
    Allows exact zero eigenvalues, unlike Cholesky.

    Returns:
      size == 1 : shape (nbins,)
      size > 1  : shape (size, nbins)
    """
    mean = np.asarray(mean, dtype=float)
    cov = np.asarray(cov, dtype=float)

    # Force exact symmetry.
    cov = 0.5 * (cov + cov.T)

    eigvals, eigvecs = np.linalg.eigh(cov)

    min_eig = np.min(eigvals)
    if min_eig < -tol:
        print("WARNING: covariance has negative eigenvalue:", min_eig)

    # Clip tiny negative numerical noise.
    eigvals = np.clip(eigvals, 0.0, None)

    z = np.random.normal(size=(size, len(mean)))
    sqrt_cov = eigvecs @ np.diag(np.sqrt(eigvals))

    throws = mean + z @ sqrt_cov.T

    if size == 1:
        return throws[0]

    return throws


def ThrowSystematics(
    histogram,
    throwFlux=False,
    useDataSubMCCov=True,
    n_samples=50,
    doDiagnostics=False,
):
    pred_vals = np.array(histogram.GetMCHistogram())[1:-1]

    # Full covariance currently includes stat + all systematics.
    V_full = histogram.GetCovarianceMatrix(False)
    V_sansFlux = histogram.GetCovarianceMatrix(sansFlux=True)
    V_flux = V_full - V_sansFlux

    # Get the stat covariance from the same residual histogram convention
    # used in StitchedHistogram.SetCovarianceMatrices().
    h_test = histogram.GetDataHistogram().Clone()
    h_test.Add(histogram.GetMCHistogram(), -1)

    V_stat = TMatrix_to_Numpy(h_test.GetStatErrorMatrix())[1:-1, 1:-1]

    # Remove stat from the Gaussian throw covariance.
    # Poisson will handle statistical fluctuations later.
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
        # Same strategy as Ryan's script:
        # one independent flux throw per toy.
        flux_throws = throw_multivariate_psd(
            pred_vals,
            V_flux,
            size=n_samples,
        )

        flux_throws = np.asarray(flux_throws)
        if flux_throws.ndim == 1:
            flux_throws = flux_throws.reshape(1, -1)

        # Then one independent non-flux systematic throw per toy,
        # centered on that toy's flux-shifted mean.
        sys_throws = []
        for flux_mean in flux_throws:
            toy_mean = throw_multivariate_psd(
                flux_mean,
                V_sansFlux_systOnly,
                size=1,
            )
            sys_throws.append(toy_mean)

        sys_throws = np.asarray(sys_throws)

    else:
        # No explicit flux profiling:
        # throw all systematics together, but still no stat in the Gaussian covariance.
        sys_throws = throw_multivariate_psd(
            pred_vals,
            V_full_systOnly,
            size=n_samples,
        )

    sys_throws = np.asarray(sys_throws)
    if sys_throws.ndim == 1:
        sys_throws = sys_throws.reshape(1, -1)

    print("ThrowSystematics shape:", sys_throws.shape)
    return sys_throws
# def ThrowSystematics(
#     histogram,
#     throwFlux=False,
#     useDataSubMCCov=True,
#     n_samples=50,
#     doDiagnostics=False,
# ):
#     includeStatError = False
#     errorAsFraction  = True
#     useOnlyShapeErrors = False

#     pred_vals = np.array(histogram.GetMCHistogram())[1:-1]

#     V_full = histogram.GetCovarianceMatrix()
#     V_sansFlux = histogram.GetCovarianceMatrix(sansFlux=True)
#     V_flux = V_full - V_sansFlux

#     if doDiagnostics:
#         print_cov_diagnostics("V_full", V_full)
#         print_cov_diagnostics("V_sansFlux", V_sansFlux)
#         print_cov_diagnostics("V_flux = V_full - V_sansFlux", V_flux)
#         print_cov_diagnostics("V_flux @ V_flux.T", V_flux @ V_flux.T)

#     if useDataSubMCCov:
#         covMatrix = V_full
#     else:
#         covMatrixTmp  = histogram.GetMCHistogram().GetTotalErrorMatrix(
#             includeStatError, errorAsFraction, useOnlyShapeErrors
#         )
#         covMatrixTmp += histogram.GetPseudoHistogram().GetTotalErrorMatrix(
#             includeStatError, errorAsFraction, useOnlyShapeErrors
#         )

#         mc_hist = histogram.GetMCHistogram()
#         for i in range(mc_hist.GetNbinsX()+1):
#             for j in range(i, mc_hist.GetNbinsX()+1):
#                 covMatrixTmp[i][j] = (
#                     covMatrixTmp[i][j]
#                     * mc_hist.GetBinContent(i)
#                     * mc_hist.GetBinContent(j)
#                 )
#                 covMatrixTmp[j][i] = covMatrixTmp[i][j]

#         covMatrix = np.asarray(matrix(covMatrixTmp))[1:-1, 1:-1]

#     if throwFlux:
#         # V_flux is already the flux covariance.
#         # Do NOT do V_flux @ V_flux.T.
#         fluxCov = 0.5 * (V_flux + V_flux.T)

#         # First throw flux around the CV prediction.
#         flux_throws = throw_multivariate_psd(
#             pred_vals,
#             fluxCov,
#             size=1,
#         )

#         # Then throw non-flux systematics around the flux-shifted prediction.
#         sansFluxCov = 0.5 * (V_sansFlux + V_sansFlux.T)

#         sys_throws = throw_multivariate_psd(
#             flux_throws,
#             sansFluxCov,
#             size=n_samples,
#         )

#     else:
#         covMatrix = 0.5 * (covMatrix + covMatrix.T)

#         sys_throws = throw_multivariate_psd(
#             pred_vals,
#             covMatrix,
#             size=n_samples,
#         )

#     sys_throws = np.asarray(sys_throws)
#     if sys_throws.ndim == 1:
#         sys_throws = sys_throws.reshape(1, -1)

#     print("ThrowSystematics shape:", sys_throws.shape)
#     return sys_throws

def ThrowPoissons(lambdas, histogram, throwFlux=False):
    lambdas = np.asarray(lambdas)

    # Force shape = (n_toys, n_bins)
    if lambdas.ndim == 1:
        lambdas = lambdas.reshape(1, -1)

    print("ThrowPoissons input shape:", lambdas.shape)

    throws = []
    for lam in lambdas:
        while lam[lam < 0].any():
            logging.error("Negative value in sys throws, rerunning this throw...")
            lam = ThrowSystematics(
                histogram,
                throwFlux=throwFlux,
                n_samples=1,
                doDiagnostics=False,
            )[0]

        throw = np.random.poisson(lam)
        throws.append(throw)

    throws = np.asarray(throws)

    if throws.ndim == 1:
        throws = throws.reshape(1, -1)

    print("ThrowPoissons output shape:", throws.shape)

    return throws

def FitToyExperiments(
    histogram,
    experiments,
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

    print("FitToyExperiments experiments shape:", experiments.shape)

    results = []

    for itoy, toy in enumerate(experiments):
        weights = stitched_data.Clone().GetCVHistoWithStatError()

        for i in range(1, weights.GetNbinsX()+1):
            weight = stitched_data.GetBinContent(i) / toy[i-1] if toy[i-1] != 0 else stitched_data.GetBinContent(i)
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
        results.append(dchi2)

        histogram.SetDataHistogram(stitched_data)

        logging.info(
            "Toy {}: chi2_null = {:.6f}, chi2_fit = {:.6f}, dchi2 = {:.6f}, penalty = {:.6f}".format(
                itoy, chi2_null, chi2_fit, dchi2, null_penalty
            )
        )

    histogram.SetDataHistogram(stitched_data)
    return np.array(results)

def SplitList(alist,wanted_parts=1):
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] for i in range(wanted_parts) ]


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

    profile_only = fit_args.profile_only
    if profile_only is not None and profile_only.strip().lower() in ["", "none"]:
        profile_only = None

    profile_n_universes = fit_args.profile_n_universes
    if exclude is None:
        exclude = ""
    elif isinstance(exclude, str):
        if exclude.strip().lower() in ["none", ""]:
            exclude = ""

    mask_spec = None


    mode = "profiledFlux" if profileFlux else "noFluxProfile"

    # Grid-safe defaults:
    #   old behavior: 50 toys/job, no covariance-only exit
    # Local override examples:
    #   N_TOYS=2 COV_DIAGNOSTIC_ONLY=1 DO_COV_DIAGNOSTICS=1 python fitAsimovs_quinn.py ...
    n_samples = int(os.environ.get("N_TOYS", "50"))
    cov_diagnostic_only = os.environ.get("COV_DIAGNOSTIC_ONLY", "0") == "1"
    do_cov_diagnostics = os.environ.get("DO_COV_DIAGNOSTICS", "0") == "1"
    nominal_fit_only = os.environ.get("NOMINAL_FIT_ONLY", "0") == "1"

    process_tag = process if process is not None else "local"

    print("\n===== FitAsimovs setup =====")
    print("file        =", file_path)
    print("profileFlux =", profileFlux)
    print("mode        =", mode)
    print("lambda      =", lam)
    print("profile_only =", profile_only)
    print("profile_n_universes =", profile_n_universes)
    print("exclude     =", exclude)
    print("n_samples   =", n_samples)
    print("process     =", process_tag)
    print("output_dir  =", AnalysisConfig.output_dir)
    print("cov diag only      =", cov_diagnostic_only)
    print("do cov diagnostics =", do_cov_diagnostics)
    print("nominal fit only =", nominal_fit_only)

    print("\n===== Loaded stitched covariance check =====")
    print("external covariances =", sample_histogram.external_covariances.keys())
    print("full covariance shape =", sample_histogram.GetCovarianceMatrix(False).shape)
    print("sans-flux covariance shape =", sample_histogram.GetCovarianceMatrix(True).shape)

    if "fhc_elastic" in sample_histogram.external_covariances:
        print("first 6x6 full covariance block:")
        print(sample_histogram.GetCovarianceMatrix(False)[:6, :6])

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
        print("\nDone nominal no-throw fit. Exiting before toys.")
        sys.exit(0)

    sys_throws = ThrowSystematics(
        sample_histogram,
        throwFlux=profileFlux,
        n_samples=n_samples,
        doDiagnostics=do_cov_diagnostics,
    )

    if cov_diagnostic_only:
        print("\nDone covariance/throw diagnostic. Exiting before Poisson + fits.")
        sys.exit(0)

    experiments = ThrowPoissons(
        sys_throws,
        sample_histogram,
        throwFlux=profileFlux,
    )

    dchi2s = FitToyExperiments(
        sample_histogram,
        experiments,
        profileFlux=profileFlux,
        lam=lam,
        exclude=exclude,
        mask_spec=mask_spec,
        profile_only=profile_only,
        profile_n_universes=profile_n_universes,
    )

    exclude_tag = "includeAll" if exclude in [None, ""] else "exclude" + str(exclude).replace(",", "-")

    profile_only_tag = ""
    if profile_only not in [None, "", "none"]:
        profile_only_tag = "_profileOnly" + str(profile_only).replace(",", "-")

    nprof_tag = "NprofAll" if profile_n_universes is None else "Nprof{}".format(profile_n_universes)

    outname = "{}/sample_dchi2s_{}_{}_{}_{}{}_{}.csv".format(
        AnalysisConfig.output_dir,
        hist_config_tag,
        mode,
        exclude_tag,
        nprof_tag,
        profile_only_tag,
        process_tag,
    )

    np.savetxt(outname, dchi2s, delimiter=",")

    print("Saved:", outname)
    print("dchi2 min/max/mean =", np.nanmin(dchi2s), np.nanmax(dchi2s), np.nanmean(dchi2s))

# if __name__ == "__main__":
#     filename = "rootfiles/NuE_stitched_hists_newNuE.root"
#     file_path = "{}/oscillations/{}".format(ccnueroot, filename)

#     sample_histogram = StitchedHistogram("sample")
#     sample_histogram.Load(file_path)

#     profileFlux = getattr(AnalysisConfig, "profileFlux", True)
#     lam = getattr(AnalysisConfig, "lambdaValue", 1)
#     exclude = getattr(AnalysisConfig, "exclude", "none")
#     mode = "profiledFlux" if profileFlux else "noFluxProfile"

#     if profileFlux:
#         if exclude in [None, "", "none"]:
#             exclude = "ratio,elastic"
#         else:
#             parts = [x.strip() for x in exclude.split(",") if x.strip()]
#             for needed in ["ratio", "elastic"]:
#                 if needed not in parts:
#                     parts.append(needed)
#             exclude = ",".join(parts)

#     # Number of toys per grid job.
#     # Keep 50 as old default unless you add a config option later.
#     n_samples = 2

#     # PROCESS can be None for local tests.
#     process_tag = process if process is not None else "local"

#     mode = "profiledFlux" if profileFlux else "noFluxProfile"

#     print("\n===== FitAsimovs setup =====")
#     print("profileFlux =", profileFlux)
#     print("mode        =", mode)
#     print("lambda      =", lam)
#     print("exclude     =", exclude)
#     print("n_samples   =", n_samples)
#     print("process     =", process_tag)
#     print("output_dir  =", AnalysisConfig.output_dir)

#     # For no flux profiling:
#     #   throw from the full covariance and fit with full covariance.
#     #
#     # For flux profiling:
#     #   the original script used throwFlux=True, separating flux and sans-flux.
#     #   Keep that behavior for profiled mode.
#     sys_throws = ThrowSystematics(
#         sample_histogram,
#         throwFlux=profileFlux,
#         n_samples=n_samples
#     )

#     experiments = ThrowPoissons(sys_throws, sample_histogram)

#     dchi2s = FitToyExperiments(
#         sample_histogram,
#         experiments,
#         profileFlux=profileFlux,
#         lam=lam,
#         exclude=exclude
#     )

#     process_tag = process if process is not None else "local"

#     print("AnalysisConfig.output_dir", AnalysisConfig.output_dir)
#     outname = "{}/sample_dchi2s_{}_{}.csv".format(
#         AnalysisConfig.output_dir,
#         mode,
#         process_tag
#     )

#     np.savetxt(outname, dchi2s, delimiter=',')

#     print("Saved:", outname)
#     print("dchi2 min/max/mean =", np.nanmin(dchi2s), np.nanmax(dchi2s), np.nanmean(dchi2s))

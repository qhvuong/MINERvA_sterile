import os
import logging, sys
import ROOT
import PlotUtils
import numpy as np
np.set_printoptions(precision=1)
np.set_printoptions(linewidth=1520)
np.set_printoptions(threshold=sys.maxsize)
from scipy import optimize, integrate
from scipy.stats import multivariate_normal
import argparse
ccnueroot = os.environ.get('CCNUEROOT')
process = os.environ.get('PROCESS')

import math
from array import array

#insert path for modules of this package.
from tools.PlotLibrary import HistHolder
from tools.Fitters import *
from tools.StitchedHistogram import *

from config.AnalysisConfig import AnalysisConfig

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

def ThrowSystematics(histogram,throwFlux=False,useDataSubMCCov=True,n_samples=50):
    includeStatError = False
    errorAsFraction  = True
    useOnlyShapeErrors = False

    pred_vals = np.array(histogram.GetMCHistogram())[1:-1] # store MC bin contents excluding over/underflow bins

    if useDataSubMCCov:
        covMatrix = histogram.GetCovarianceMatrix()
    else:
        covMatrixTmp  = histogram.GetMCHistogram().GetTotalErrorMatrix(includeStatError, errorAsFraction, useOnlyShapeErrors)
        covMatrixTmp += histogram.GetPseudoHistogram().GetTotalErrorMatrix(includeStatError, errorAsFraction, useOnlyShapeErrors)

        mc_hist = histogram.GetMCHistogram()
        for i in range(mc_hist.GetNbinsX()+1):
            for j in range(i,mc_hist.GetNbinsX()+1):
                covMatrixTmp[i][j] = covMatrixTmp[i][j] * mc_hist.GetBinContent(i) * mc_hist.GetBinContent(j)
                covMatrixTmp[j][i] = covMatrixTmp[i][j]

        covMatrix = np.asarray(matrix(covMatrixTmp))[1:-1,1:-1] # convert to numpy array, exclude over/underflow bins

    if throwFlux:
        fluxCov = histogram.GetCovarianceMatrix() - histogram.GetCovarianceMatrix(sansFlux=True)
        fluxCov = fluxCov @ np.conj(fluxCov.T)
        flux_throws = multivariate_normal.rvs(mean=pred_vals,cov=fluxCov) # randomly sample covariance matrix around null hypothesis

        sansFluxCov = histogram.GetCovarianceMatrix(sansFlux=True)
        sys_throws  = multivariate_normal.rvs(mean=flux_throws,cov=sansFluxCov,size=n_samples) # randomly sample covariance matrix around null hypothesis
    else:
        sys_throws = multivariate_normal.rvs(mean=pred_vals,cov=covMatrix,size=n_samples) # randomly sample covariance matrix around null hypothesis

    sys_throws = np.asarray(sys_throws)

    # scipy can return shape (nbins,) instead of (1, nbins) for one throw.
    # Force always shape = (n_toys, n_bins).
    if sys_throws.ndim == 1:
        sys_throws = sys_throws.reshape(1, -1)

    print("ThrowSystematics shape:", sys_throws.shape)

    return sys_throws

def ThrowPoissons(lambdas, histogram):
    lambdas = np.asarray(lambdas)

    # Force shape = (n_toys, n_bins)
    if lambdas.ndim == 1:
        lambdas = lambdas.reshape(1, -1)

    print("ThrowPoissons input shape:", lambdas.shape)

    throws = []
    for lam in lambdas:
        while lam[lam < 0].any():
            logging.error("Negative value in sys throws, rerunning this throw...")
            lam = ThrowSystematics(histogram, n_samples=1)[0]

        throw = np.random.poisson(lam)
        throws.append(throw)

    throws = np.asarray(throws)

    if throws.ndim == 1:
        throws = throws.reshape(1, -1)

    print("ThrowPoissons output shape:", throws.shape)

    return throws

def FitToyExperiments(histogram, experiments, profileFlux=True, lam=1, exclude="none"):
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

        stat = Statistics(histogram, lam=lam, exclude=exclude)

        chi2_null, null_penalty = stat.Chi2DataMC(
            marginalize=profileFlux
        )

        fitter = OscillationFitter(
            histogram,
            lam=lam,
            exclude=exclude,
            marginalize_flux=profileFlux
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

if __name__ == "__main__":
    filename = "rootfiles/NuE_stitched_hists.root"
    file_path = "{}/oscillations/{}".format(ccnueroot, filename)

    sample_histogram = StitchedHistogram("sample")
    sample_histogram.Load(file_path)

    profileFlux = getattr(AnalysisConfig, "profileFlux", True)
    lam = getattr(AnalysisConfig, "lambdaValue", 1)
    exclude = getattr(AnalysisConfig, "exclude", "none")
    mode = "profiledFlux" if profileFlux else "noFluxProfile"

    # Number of toys per grid job.
    # Keep 50 as old default unless you add a config option later.
    n_samples = 50

    # PROCESS can be None for local tests.
    process_tag = process if process is not None else "local"

    mode = "profiledFlux" if profileFlux else "noFluxProfile"

    print("\n===== FitAsimovs setup =====")
    print("profileFlux =", profileFlux)
    print("mode        =", mode)
    print("lambda      =", lam)
    print("exclude     =", exclude)
    print("n_samples   =", n_samples)
    print("process     =", process_tag)
    print("output_dir  =", AnalysisConfig.output_dir)

    # For no flux profiling:
    #   throw from the full covariance and fit with full covariance.
    #
    # For flux profiling:
    #   the original script used throwFlux=True, separating flux and sans-flux.
    #   Keep that behavior for profiled mode.
    sys_throws = ThrowSystematics(
        sample_histogram,
        throwFlux=profileFlux,
        n_samples=1
    )

    experiments = ThrowPoissons(sys_throws, sample_histogram)

    dchi2s = FitToyExperiments(
        sample_histogram,
        experiments,
        profileFlux=profileFlux,
        lam=lam,
        exclude=exclude
    )

    process_tag = process if process is not None else "local"

    print("AnalysisConfig.output_dir", AnalysisConfig.output_dir)
    outname = "{}/sample_dchi2s_{}_{}.csv".format(
        AnalysisConfig.output_dir,
        mode,
        process_tag
    )

    np.savetxt(outname, dchi2s, delimiter=',')

    print("Saved:", outname)
    print("dchi2 min/max/mean =", np.nanmin(dchi2s), np.nanmax(dchi2s), np.nanmean(dchi2s))

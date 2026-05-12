import os
import logging, sys
import copy
import ROOT
import PlotUtils
import numpy as np
np.set_printoptions(precision=4)
np.set_printoptions(linewidth=1520)
np.set_printoptions(threshold=sys.maxsize)
from scipy import optimize, integrate
from config.AnalysisConfig import AnalysisConfig

import argparse
ccnueroot = os.environ.get('CCNUEROOT')

import math
import psutil
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

    filename = "rootfiles/NuE_stitched_hists.root"
    file_path = "{}/oscillations/{}".format(ccnueroot, filename)

    lambda_value = getattr(AnalysisConfig, "lambdaValue", 1)
    exclude_samples = getattr(AnalysisConfig, "exclude", [])

    sample_histogram = StitchedHistogram("sample")
    sample_histogram.Load(file_path)

    stat = Statistics(sample_histogram, lam=lambda_value, exclude=exclude_samples)
    chi2_null, penalty = stat.Chi2DataMC(marginalize=False)
    print("null chi2: {:.3f}".format(chi2_null))

    # fitter = OscillationFitter(sample_histogram, lam=lambda_value, exclude=exclude_samples)
    fitter = OscillationFitter(
        sample_histogram,
        lam=lambda_value,
        exclude=exclude_samples,
        marginalize_flux=False,
    )
    chi2_fit, res = fitter.DoFit()

    print("Data fit: delta chi2 = {:.3f} = {:.3f} - {:.3f}".format(chi2_null - chi2_fit, chi2_null, chi2_fit))
    print("Best fit params:")
    print("   delta m^2 = {:.3f} eV^2 +- {:.4f}".format(res["m"], 0))
    print("   U_e4^2    = {:.3f}      +- {:.4f}".format(res["ue4"], 0))
    print("   U_mu4^2   = {:.5f}    +- {:.4f}".format(res["umu4"], 0))
    print("   U_tau4^2  = {:.3f}      +- {:.4f}".format(res["utau4"], 0))

    sample_histogram.OscillateHistogram(res["m"], res["ue4"], res["umu4"], res["utau4"])
    sample_histogram.SetPlottingStyle()

    plotter = PlottingContainer("fitted_histogram", sample_histogram)
    plotter.SetExclude(exclude_samples)
    plotter.SetLambda(lambda_value)

    # plotter.PlotOscillationEffects(res, AnalysisConfig.ntuple_tag, plotSamples=False)
    plotter.PlotOscillationEffects(
        res,
        AnalysisConfig.ntuple_tag,
        useMarg=False,
        plotSamples=False
    )
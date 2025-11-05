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
#import multiprocessing
import time
#import threading
#nthreads = 4
from array import array

#insert path for modules of this package.
from tools.PlotLibrary import HistHolder
from Tools.Histogram import *
from Tools.PlotHistogram import *

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

MNVPLOTTER = PlotUtils.MnvPlotter()
#config MNVPLOTTER:
#MNVPLOTTER.chi2_use_overflow_err = True
MNVPLOTTER.draw_normalized_to_bin_width=False
MNVPLOTTER.legend_text_size = 0.04
#MNVPLOTTER.extra_top_margin = -.035# go slightly closer to top of pad
MNVPLOTTER.mc_bkgd_color = 46 
MNVPLOTTER.mc_bkgd_line_color = 46
MNVPLOTTER.legend_n_columns = 3
#MNVPLOTTER.mc_line_width = 0
MNVPLOTTER.data_bkgd_color = 12 #gray
MNVPLOTTER.data_bkgd_style = 24 #circle  
MNVPLOTTER.axis_maximum = 500 #circle  

#legend entries are closer
MNVPLOTTER.height_nspaces_per_hist = 1.2
MNVPLOTTER.width_xspace_per_letter = .4
MNVPLOTTER.legend_text_size        = .03
MNVPLOTTER.legend_offset_x           = .15

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

if __name__ == "__main__":
    filename = "rootfiles/NuE_stitched_hists.root"
    
    file_path = "{}/oscillations/{}".format(ccnueroot,filename)

    sample_histogram = StitchedHistogram("sample")
    sample_histogram.Load(file_path)

    invCov = sample_histogram.GetInverseCovarianceMatrix(sansFlux=True)

    chi2_null,penalty = Chi2DataMC(sample_histogram,invCov=invCov,marginalize=True,lam=AnalysisConfig.lambdaValue,exclude=AnalysisConfig.exclude)
    print("null chi2: {:.3f}".format(chi2_null))

    fitter = Fitter(sample_histogram,invCov=invCov,lam=AnalysisConfig.lambdaValue,exclude=AnalysisConfig.exclude)
    chi2_fit,res = fitter.DoFit()

    print("Data fit: delta chi2 = {:.3f} = {:.3f} - {:.3f}".format(chi2_null-chi2_fit,chi2_null,chi2_fit))
    print("Best fit params:")
    print("   delta m^2 = {:.3f} eV^2 +- {:.4f}".format(res['m'],0))
    print("   U_e4^2    = {:.3f}      +- {:.4f}".format(res['ue4'],0))
    print("   U_mu4^2   = {:.5f}    +- {:.4f}".format(res['umu4'],0))
    print("   U_tau4^2  = {:.3f}      +- {:.4f}".format(res['utau4'],0))
    sample_histogram.SetPlottingStyle()

    plotter = PlottingContainer("fitted_histogram",sample_histogram)
    plotter.SetExclude(AnalysisConfig.exclude)
    plotter.SetInverseCovariance(invCov)
    plotter.SetLambda(AnalysisConfig.lambdaValue)

    plotter.PlotOscillationEffects(res,AnalysisConfig.ntuple_tag,plotSamples=True)
    #plotter.PlotFluxMarginalizationEffects(res,"bestfit")

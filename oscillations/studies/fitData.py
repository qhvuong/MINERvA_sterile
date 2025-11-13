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
from Tools.PlotTools import *

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
    BLUEARC = "/exp/minerva/data/users/{}/surfaces".format(os.environ["USER"])
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--grid",
                        action="store_true",
                        default = False,
                        help = "Run macro on grid, Input/Output path must be updated to avoid direct access to BlueArc"
    )
    parser.add_argument("-o", "--output",
                        dest = "output_dir",
                        help="Use alternate location for output file.",
                        default=BLUEARC
    )
    parser.add_argument("-m", "--delta_m",
                        dest = "delta_m",
                        help="Delta m^2 value to probe.",
                        type=float,
                        default=0
    )
    parser.add_argument("-Ue4", "--U_e4",
                        dest = "U_e4",
                        help="U_e4 parameter to probe.",
                        type=float,
                        default=0
    )
    parser.add_argument("-Umu4", "--U_mu4",
                        dest = "U_mu4",
                        help="U_mu4 parameter to probe.",
                        type=float,
                        default=0
    )
    parser.add_argument("-U", "--U_tau4",
                        dest = "U_tau4",
                        help="U_tau4 parameter to probe.",
                        type=float,
                        default=0
    )
    parser.add_argument("--pseudodata",
                        dest = "pseudodata",
                        default = False,
                        action="store_true",
    )

    args = parser.parse_args()
    outdir_surface = args.output_dir
    delta_m = args.delta_m
    runongrid = args.grid
    U_e4 = args.U_e4
    U_mu4 = args.U_mu4
    U_tau4 = args.U_tau4

    filename = ""
    filename = "NuE_stitched_hists.root"
    
    file_path = "{}/oscillations/{}".format(ccnueroot,filename)

    sample_histogram = StitchedHistogram("sample")
    sample_histogram.Load(file_path)

    invCovariance = sample_histogram.GetInverseCovarianceMatrix()

    chi2_null,penalty = Chi2DataMC(sample_histogram,invCov=invCovariance,marginalize=True)
    print("null chi2: {:.3f}".format(chi2_null))

    fitter = Fitter(sample_histogram,invCov=invCovariance)
    chi2_fit,res = fitter.DoFit()

    print("Data fit: delta chi2 = {:.3f} = {:.3f} - {:.3f}".format(chi2_null-chi2_fit,chi2_null,chi2_fit))
    print("Best fit params:")
    print("   delta m^2 = {:.3f} eV^2 +- {:.4f}".format(res['m'],0))
    print("   U_e4^2    = {:.3f}      +- {:.4f}".format(res['ue4'],0))
    print("   U_mu4^2   = {:.5f}    +- {:.4f}".format(res['umu4'],0))
    print("   U_tau4^2  = {:.3f}      +- {:.4f}".format(res['utau4'],0))

    PlotOscillationEffects(sample_histogram,res,"bestfit",plotSamples=False)
    PlotOscillationRatios(sample_histogram,res,"bestfit")
    PlotFluxMarginalizationEffects(sample_histogram,res,"bestfit")

import os
import logging, sys
import copy
import ROOT
import PlotUtils
import numpy as np
from scipy import optimize, integrate

import argparse
ccnueroot = os.environ.get('CCNUEROOT')

import math
from array import array

#insert path for modules of this package.
from tools.PlotLibrary import HistHolder
from Tools.Histogram import *
from Tools.PlotHistogram import *
from config.AnalysisConfig import AnalysisConfig

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

np.set_printoptions(precision=4)
np.set_printoptions(linewidth=1520)
np.set_printoptions(threshold=sys.maxsize)

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

def CompareResults(histogram):
    dataHist = histogram.GetDataHistogram()
    mcHist = histogram.GetMCHistogram()

    data = np.array(dataHist)[1:-1]
    mc = np.array(mcHist)[1:-1]

    universes = histogram.GetFluxUniverses()
    invCov = histogram.GetInverseCovarianceMatrix(sansFlux=True) 

    A = histogram.GetAMatrix()
    V = invCov    
    C = data - mc
    I = np.identity(len(universes))

    L = 2 * A @ V @ C
    Q = A @ V @ A.T + I
    solution = np.linalg.inv(Q) @ L/2

if __name__ == "__main__":
    filename = "NuE_stitched_hists.root"
    file_path = "{}/oscillations/rootfiles/{}".format(ccnueroot,filename)

    sample_histogram = StitchedHistogram("sample")
    sample_histogram.Load(file_path)
    sample_histogram.SetPlottingStyle()
    invCov=sample_histogram.GetInverseCovarianceMatrix(sansFlux=True)

    n4 = {
            'm':7,
            'ue4':.14,
            'umu4':.01,
            'utau4':0.
        }
    nulltest = {
            'm':0,
            'ue4':0,
            'umu4':0,
            'utau4':0
        }
    bestfit = {
            'm':0.18,
            'ue4':0.15,
            'umu4':0.41,
            'utau4':0.0
        }

    plotter = PlottingContainer("test",sample_histogram)
    plotter.SetInverseCovariance(invCov)
    plotter.SetLambda(AnalysisConfig.lambdaValue)
    plotter.SetExclude(AnalysisConfig.exclude)
    nullSolution,nullPen = FluxSolution(sample_histogram,invCov=invCov)
    sample_histogram.PlotSamples(fluxSolution=nullSolution,plotName="AllSamples")
    sample_histogram.PlotStitchedHistogram(fluxSolution=nullSolution,plotName="StitchedHistogram")

    plotter.PlotOscillationEffects(bestfit,"BestFit",plotSamples=True)
    plotter.PlotOscillationEffects(n4,"n4",plotSamples=True)
    plotter.PlotOscillationEffects(nulltest,"nulltest",plotSamples=True)
    U_mu4s = 0.41*np.logspace(-5,0,100)
    U_mu4s[0] = 0
    for i,umu4 in enumerate(U_mu4s):
        n4['umu4'] = umu4
        plotter.PlotOscillationEffects(n4,"Neutrino4_{:02d}".format(i),plotSamples=True)
    exit()


    #for i,tau4 in enumerate(np.linspace(0,0.66,20)):
    #    nulltest['m'] = 10
    #    nulltest['ue4'] = 0.1
    #    nulltest['umu4'] = 0.01
    #    nulltest['utau4'] = tau4
    #    plotter.PlotOscillationEffects(nulltest,"Utau4_{:02d}".format(i),plotSamples=True)

    for i,dm2 in enumerate(np.logspace(-1,2,100)):
        n4['m'] = dm2
        n4['ue4'] = 0.1
        n4['umu4'] = 0.01
        plotter.PlotOscillationEffects(n4,"Deltam2_{:02d}".format(i),plotSamples=True)
    

    #invCov = np.loadtxt("data_minus_mc_COV.csv",delimiter=',')
    #invCov = np.linalg.inv(invCov)

    #OscillateHistogram(sample_histogram, n4['m'], n4['ue4'], n4['umu4'], n4['utau4'])
    #plotter.SetInverseCovariance(invCov)
    #plotter.PlotScatteringIntegrals()
    #plotter.PlotFluxReweight()
    #plotter.PlotProfileEffects()
    #plotter.PlotFluxProfilingEffects("profile_effect")
    #plotter.PlotOscillationEffects(n4,"Neutrino4",plotSamples=True)
    

    #OscillateHistogram(sample_histogram, 0, 0, 0, 0)
    #plotter.PlotOscillationEffects(nulltest,"nulltest_excludeRatio",plotSamples=False)
    #nullSolution,nullPen = FluxSolution(sample_histogram,invCov=invCov)

    #sample_histogram.PlotSamples(fluxSolution=nullSolution,plotName="AllSamples")
    #PlotOscillationEffects(sample_histogram,n4,"Neutrino4",plotSamples=True)
    #PlotOscillationRatios(sample_histogram,n4,"Neutrino4")
    #PlotFluxMarginalizationEffects(sample_histogram,n4,"Neutrino4")
    #PlotSampleMarginalizationEffects(sample_histogram)

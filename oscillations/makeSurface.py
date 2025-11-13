import os
import logging, sys
import ROOT
import PlotUtils
import numpy as np
from scipy import optimize, integrate
np.set_printoptions(precision=5)
np.set_printoptions(linewidth=1520)
np.set_printoptions(threshold=sys.maxsize)
ccnueroot = os.environ.get('CCNUEROOT')

import math
from array import array

#insert path for modules of this package.
from tools.PlotLibrary import HistHolder
from Tools.FitTools import *
from Tools.PlotTools import *
from config.AnalysisConfig import AnalysisConfig

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

def MakeSurface(histogram,outdir,deltam=1,U_tau4=0,makePlot=False,exclude="",lam=1):
    U_mu4s = 0.41*np.logspace(-5,0,100)
    U_mu4s[0] = 0
    U_e4s = 0.15*np.logspace(-4,0,100)
    U_e4s[0] = 0

    arrShape = (np.shape(U_mu4s)[0],np.shape(U_e4s)[0])

    asimov_surface    = np.zeros(arrShape,dtype='f')
    data_surface      = np.zeros(arrShape,dtype='f')
    data_penalties    = np.zeros(arrShape,dtype='f')
    asimov_penalties  = np.zeros(arrShape,dtype='f')
    fits              = np.zeros(arrShape,dtype='f')
    count = 0

    for i in range(U_mu4s.shape[0]):
        count+=1
        for j in range(U_e4s.shape[0]):
            U_mu4 = U_mu4s[i]
            U_e4  = U_e4s[j]

            OscillateHistogram(histogram, deltam, U_e4, U_mu4, U_tau4)

            chi2_data,data_penalty = Chi2DataMC(histogram,invCov=histogram.GetInverseCovarianceMatrix(sansFlux=True),marginalize=True,useOsc=True,exclude=exclude,lam=lam)
            chi2_asimov,asimov_penalty = Chi2DataMC(histogram,invCov=histogram.GetInverseCovarianceMatrix(sansFlux=True),marginalize=True,useOsc=True,usePseudo=True,exclude=exclude,lam=lam)

            if makePlot:
                res={"m":deltam,"ue4":U_e4,"umu4":U_mu4,"utau4":0}
                PlotOscillationEffects(histogram,res,"asimov",False,useMarg=False,usePseudo=True,exclue=exclude,lam=lam)
                PlotOscillationEffects(histogram,res,"asimov_marg",True,useMarg=True,usePseudo=True,exclude=exclude,lam=lam)

            data_surface[i,j] = chi2_data
            data_penalties[i,j] = data_penalty
            asimov_surface[i,j] = chi2_asimov
            asimov_penalties[i,j] = asimov_penalty
        logging.info("{:.2f}% done with chi2s. Current data, asimov chi2s = {:.4f}, {:.4f}".format(100*count/(U_mu4s.shape[0]),data_surface[i,-1],asimov_surface[i,-1]))

    np.save('{}/chi2_surface_data_m_{}_Ue4_{}.dat'.format(outdir,deltam,U_e4),data_surface)
    np.save('{}/chi2_surface_pseudodata_m_{}_Ue4_{}.dat'.format(outdir,deltam,U_e4),asimov_surface)
    np.save('{}/chi2_penalty_data_m_{}_Ue4_{}.dat'.format(outdir,deltam,U_e4),data_penalties)
    np.save('{}/chi2_penalty_pseudodata_m_{}_Ue4_{}.dat'.format(outdir,deltam,U_e4),asimov_penalties)

if __name__ == "__main__":
    filename = "NuE_stitched_hists.root"
    file_path = "{}/oscillations/rootfiles/{}".format(ccnueroot,filename)

    sample_histogram = StitchedHistogram("sample")
    sample_histogram.Load(file_path)

    delta_ms = np.logspace(-1,2,100)
    delta_m = delta_ms[AnalysisConfig.delta_m]

    cat_to_exclude = AnalysisConfig.exclude_systematic
    sample_histogram.RemoveSystematics(cat_to_exclude)

    if not AnalysisConfig.grid: # surface plot
        m_toloop = np.logspace(-1,2,100)
        for m in m_toloop:
            print("running over delta_m^2 = {}".format(m))
            MakeSurface(sample_histogram,AnalysisConfig.output_dir,m,AnalysisConfig.U_tau4,lam=AnalysisConfig.lambdaValue,exclude=AnalysisConfig.exclude)
    else:
        MakeSurface(sample_histogram,AnalysisConfig.output_dir,delta_m,AnalysisConfig.U_tau4,lam=AnalysisConfig.lambdaValue,exclude=AnalysisConfig.exclude)

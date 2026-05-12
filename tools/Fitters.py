import os
import time
import logging, sys
import argparse
import math
try:
    import psutil
except ImportError:
    psutil = None

from array import array

import numpy as np
from scipy import optimize

ccnueroot = os.environ.get('CCNUEROOT')

import ROOT
import ctypes
import PlotUtils
#insert path for modules of this package.
from tools.PlotLibrary import HistHolder
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS 

from tools.StitchedHistogram import *
from tools.Helper import *

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

MNVPLOTTER = PlotUtils.MnvPlotter()
MNVPLOTTER.error_summary_group_map.clear()
for k,v in CONSOLIDATED_ERROR_GROUPS.items():
    vec = ROOT.vector("std::string")()
    for vs in v :
        vec.push_back(vs)
    MNVPLOTTER.error_summary_group_map[k]= vec

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

class OscillationFitter():
    def __init__(self, histogram, exclude="ratio", lam=1, marginalize_flux=True):
        self.hist = histogram
        self.tol = 1e-12
        self.exclude = exclude
        self.lam = lam
        self.marginalize_flux = marginalize_flux
        self.statistic = Statistics(histogram, exclude, lam)

    def DoFit(self):
        x0 = [0.0,0.0,0.0,0.0]
        bounds = np.array([[0.0,1.0],[0.0,0.15],[0.0,0.41],[0,0.66]], dtype = float)
        cons = optimize.LinearConstraint([[0,1,1,1]],-np.inf,1)

        null = optimize.minimize(fun=self.CalChi2,x0=x0,tol=self.tol,options={"maxiter":40},method="SLSQP",bounds=bounds,constraints=cons)
        null_chi2 = float(null.fun)
        print("fit near null: {}".format(null.fun))

        res = optimize.differential_evolution(func=self.CalChi2,tol=self.tol,bounds=bounds,polish=False,x0=x0,maxiter=100,disp=True,constraints=cons)
        new_x0 = res.x
        print("best fit: {}".format(res.fun))

        res = optimize.minimize(fun=self.CalChi2,x0=new_x0,tol=self.tol,options={"maxiter":20},method="SLSQP",bounds=bounds,constraints=cons)
        chi2 = float(res.fun)
        print("polish fit: {}".format(res.fun))

        if null_chi2 < chi2:
            chi2 = null_chi2
            res = null

        return(chi2,{"m":res.x[0]*100,"ue4":res.x[1],"umu4":res.x[2],"utau4":res.x[3]})

    def CalChi2(self,x):
        ms = x[0]*100
        U_e4 = x[1]
        U_mu4 = x[2]
        U_tau4 = x[3]
        self.hist.OscillateHistogram(ms,U_e4,U_mu4,U_tau4,False,False)

        chi2, penalty = self.statistic.Chi2DataMC(
            marginalize=self.marginalize_flux,
            useOsc=True
        )

        return chi2

class FluxFitter():
    def __init__(self,histogram,exclude,lam,usePseudo=False,useOsc=False):
        self.hist = histogram
        self.exclude = exclude
        self.lam = lam

        if usePseudo:
            self.dataHist = histogram.GetPseudoHistogram()
        else:
            self.dataHist = histogram.GetDataHistogram()

        if useOsc:
            self.mcHist = histogram.GetOscillatedHistogram()
        else:
            self.mcHist = histogram.GetMCHistogram()

        self.SolveFluxSolution()

    def SetHistogram(self,hist):
        self.hist = histogram
        if usePseudo:
            self.dataHist= histogram.GetPseudoHistogram()
        else:
            self.dataHist = histogram.GetDataHistogram()

        if useOsc:
            self.mcHist = histogram.GetOscillatedHistogram()
        else:
            self.mcHist = histogram.GetMCHistogram()

        self.SolveFluxSolution()

    def SolveFluxSolution(self):
        data = np.array(self.dataHist)[1:-1]
        mc = np.array(self.mcHist)[1:-1]
        universes = self.hist.GetFluxUniverses()
        invCov = self.hist.GetInverseCovarianceMatrix(sansFlux=True)
        A = self.hist.GetAMatrix()

        sliceInds = GetSliceIndices("HIST_CONFIG.json",self.exclude,self.hist.keys)

        data = slicer(data,sliceInds)
        mc   = slicer(mc,sliceInds)
        A    = slicer(A,sliceInds,axis=1)
        V    = slicer(invCov,sliceInds)
        C = data - mc
        I = np.identity(len(universes))

        L = 2 * A @ V @ C
        Q = A @ V @ A.T + I * self.lam
        self.fluxSolution = np.linalg.inv(Q) @ L/2

    def SetFluxSolution(self,solution):
        self.fluxSolution = solution

    def GetFluxSolution(self,):
        return(self.fluxSolution)

    def MarginalizeFlux(self):
        A    = self.hist.GetAMatrix()
        solution = self.fluxSolution
        
        penalty = solution @ solution * self.lam
        new_cv = np.array(self.mcHist)[1:-1] + solution @ A
        return(new_cv,penalty)

    def ReweightToFluxSolution(self,histogram):
        mc = np.array(histogram)[1:-1]
        band = histogram.GetVertErrorBand("Flux")
        nhists = self.mcHist.GetVertErrorBand("Flux").GetNHists()
        universes = np.array([np.array(band.GetHist(l))[1:-1] for l in range(nhists)])
        cv_table = np.array([mc for l in range(len(universes))])
        A = universes - cv_table

        weights = histogram.GetCVHistoWithStatError()
        new_cv = mc + self.fluxSolution @ A

        for j in range(1,weights.GetNbinsX()+1):
            weight = weights.GetBinContent(j) / new_cv[j-1] if new_cv[j-1] != 0 else weights.GetBinContent(j)
            weights.SetBinContent(j,weight)
            weights.SetBinError(j,0)

        histogram.DivideSingle(histogram,weights)

class Statistics():
    def __init__(self,histogram,exclude="ratio",lam=1):
        self.hist = histogram
        self.exclude = exclude
        self.lam = lam
        self.nulFluxFitter = None
        self.oscFluxFitter = None

    def GetFluxFitter(self,useOsc=False):
        if useOsc:
            return(self.oscFluxFitter)
        else:
            return(self.nulFluxFitter)

    def SetHistogram(self,histogram):
        self.hist = histogram

    def Chi2DataMC(self,marginalize=True,usePseudo=False,useOsc=False): 
        ##### Get self.hists to calculate chi2 between #####
        if useOsc:
            mcHist = self.hist.GetOscillatedHistogram()
        else:
            mcHist = self.hist.GetMCHistogram()

        if usePseudo:
            dataHist = self.hist.GetPseudoHistogram()
        else:
            dataHist = self.hist.GetDataHistogram()

        #We get the number of bins and make sure it's compatible with the NxN matrix given
        if dataHist.GetNbinsX() != mcHist.GetNbinsX():
            logging.error("breaking error in Chi2DataMC")
            logging.error("The number of bins from Data ({}) and MC ({}) histograms differ. Returning -1.".format(dataHist.GetNbinsX(),mcHist.GetNbinsX()))
            return(-1)

        mc = np.array(mcHist)[1:-1] # store MC bin contents excluding over/underflow bins
        data = np.array(dataHist)[1:-1] # store data bin contents excluding over/underflow bins 
       
        invCov = self.hist.GetInverseCovarianceMatrix(sansFlux=False)
        penalty = 0

        # Do we want to marginalize over the flux systematic before calculating chi2
        if marginalize:
            fluxFitter = FluxFitter(self.hist,self.exclude,self.lam,usePseudo,useOsc)
            oldMc = mc
            mc,penalty = fluxFitter.MarginalizeFlux()
            invCov = self.hist.GetInverseCovarianceMatrix(sansFlux=True)

            if useOsc:
                self.oscFluxFitter = fluxFitter
            else:
                self.nulFluxFitter = fluxFitter

        # ----- Calculate chi2 value ----= #
        diff = data - mc
        chi2 = diff.T @ invCov @ diff + penalty # @ is numpy efficient matrix multiplication

        if abs(chi2) > 1e30:
            logging.error("chi2 has invalid value: {}".format(chi2))
            print("chi2 has invalid value: {}".format(chi2))
            return(-1)

        return(chi2,penalty)

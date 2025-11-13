import os
import time
import logging, sys
import argparse
import math
import psutil

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

from Tools.Histogram import *
from Tools.Helper import *

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

MNVPLOTTER = PlotUtils.MnvPlotter()
MNVPLOTTER.error_summary_group_map.clear();
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

class Fitter():
    def __init__(self,histogram,invCov=None,remakeCov=False,useNewUniverses=False,exclude="",lam=1):
        self.hist = histogram
        self.remakeCov = remakeCov
        self.useNewUniverses = useNewUniverses
        self.invCov = invCov
        self.tol = 1e-12
        self.exclude = exclude
        self.lam = lam

    def DoFit(self):
        x0 = [0.0,0.0,0.0,0.0]
        bounds = np.array([[0.0,1.0],[0.0,0.15],[0.0,0.41],[0,0.66]], dtype = float)
        cons = optimize.LinearConstraint([[0,1,1,1]],-np.inf,1)

        null = optimize.minimize(fun=self.CalChi2,x0=x0,tol=self.tol,options={"maxiter":20},method="SLSQP",bounds=bounds,constraints=cons)
        null_chi2 = float(null.fun)
        print("fit near null: {}".format(null.fun))

        res = optimize.differential_evolution(func=self.CalChi2,tol=self.tol,bounds=bounds,polish=False,x0=x0,maxiter=50,disp=True,constraints=cons)
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
        OscillateHistogram(self.hist,ms,U_e4,U_mu4,U_tau4,False,False)

        # Do we want to recalculate the same invCovariance matrix for each hypothesis?
        if self.remakeCov:
            self.invCov = None
        else:
            self.invCov = self.hist.GetInverseCovarianceMatrix(sansFlux=True)

        chi2,penalty = Chi2DataMC(self.hist,invCov=self.invCov,marginalize=True,useOsc=True,remakeCov=self.remakeCov,useNewUniverses=self.useNewUniverses,exclude=self.exclude,lam=self.lam)
        return(chi2)


def FluxSolution(histogram,invCov=None,useOsc=False,usePseudo=False,exclude="ratio",lam=1):
    if usePseudo:
        dataHist = histogram.GetPseudoHistogram()
    else:
        dataHist = histogram.GetDataHistogram()

    if useOsc:
        mcHist = histogram.GetOscillatedHistogram()
    else:
        mcHist = histogram.GetMCHistogram()

    data = np.array(dataHist)[1:-1]
    mc = np.array(mcHist)[1:-1]
    universes = histogram.GetFluxUniverses()

    if invCov is None:
        invCov = histogram.GetInverseCovarianceMatrix(sansFlux=True)

    A = histogram.GetAMatrix()

    sliceInds = GetSliceIndices("HIST_CONFIG.json",exclude,histogram.keys)

    data = slicer(data,sliceInds)
    mc   = slicer(mc,sliceInds)
    A    = slicer(A,sliceInds,axis=1)
    V    = slicer(invCov,sliceInds)

    C = data - mc
    I = np.identity(len(universes))

    L = 2 * A @ V @ C
    Q = A @ V @ A.T + I * lam
    solution = np.linalg.inv(Q) @ L/2
    penalty = lam * solution @ solution
    return(solution,penalty)

def ReweightCV(histogram,fluxSolution,cv=None,mc=None):
    if cv is None:
        cv = np.array(histogram)[1:-1]
    
    band = histogram.GetVertErrorBand("Flux")
    nhists = band.GetNHists()
    nhists = 100
    universes = np.array([np.array(band.GetHist(l))[1:-1] for l in range(nhists)])
    cv_table = np.array([cv for l in range(len(universes))])
    A = universes - cv_table

    if mc is None:
        mc = cv

    new_cv = mc + fluxSolution @ A

    weights = histogram.GetCVHistoWithStatError()
    for j in range(1,weights.GetNbinsX()+1):
        weight = weights.GetBinContent(j) / new_cv[j-1] if new_cv[j-1] != 0 else weights.GetBinContent(j)
        weights.SetBinContent(j,weight)
        weights.SetBinError(j,0)

    histogram.DivideSingle(histogram,weights)
    return(weights)

def MarginalizeFlux(histogram,invCov=None,fluxSolution=None,useOsc=False,usePseudo=False,setHists=False,remakeCov=False,useNewUniverses=False,exclude=None,lam=1):
    if usePseudo:
        dataHist = histogram.GetPseudoHistogram()
    else:
        dataHist = histogram.GetDataHistogram()

    if useOsc:
        mcHist = histogram.GetOscillatedHistogram()
    else:
        mcHist = histogram.GetMCHistogram()

    data = np.array(dataHist)[1:-1]
    mc = np.array(mcHist)[1:-1]
    universes = histogram.GetFluxUniverses()

    if invCov is None:
        invCov = histogram.GetInverseCovarianceMatrix(sansFlux=True)

    if useNewUniverses:
        if "Flux" in mcHist.GetVertErrorBandNames():
            band = mcHist.GetVertErrorBand("Flux")
            nhists = band.GetNHists()
            universes = np.array([np.array(band.GetHist(i))[1:-1] for i in range(nhists)])
            cv_table = np.array([mc for i in range(len(universes))])
            A = universes - cv_table
        else:
            raise ValueError("No flux universes in mc histogram")
    else:
        A = histogram.GetAMatrix()

    sliceInds = GetSliceIndices("HIST_CONFIG.json",exclude,histogram.keys)
    
    data = slicer(data,sliceInds)
    mc   = slicer(mc,sliceInds)
    A    = slicer(A,sliceInds,axis=1)
    V    = slicer(invCov,sliceInds)

    if fluxSolution is None:
        C = data - mc
        I = np.identity(len(universes))

        L = 2 * A @ V @ C
        Q = A @ V @ A.T + I * lam
        solution = np.linalg.inv(Q) @ L/2
    else:
        solution = fluxSolution

    penalty = solution @ solution * lam
    new_cv = np.array(mcHist)[1:-1] + solution @ histogram.GetAMatrix()

    if remakeCov:
        ##### grab new invCovariance matrix to reflect reweighted MC values #####
        weights = mcHist.GetCVHistoWithStatError()
        for i in range(1,weights.GetNbinsX()+1):
            weight = weights.GetBinContent(i) / new_cv[i-1] if new_cv[i-1] != 0 else weights.GetBinContent(i)
            weights.SetBinContent(i,weight)
            weights.SetBinError(i,0)

        #get the invCovariance matrix
        useOnlyShapeErrors = False
        includeStatError   = True
        errorAsFraction    = False

        dataHist.PopVertErrorBand("Flux")
        mcHist.PopVertErrorBand("Flux")
     
        new_mc = mcHist.Clone()
        new_mc.DivideSingle(new_mc,weights)
        dataHist.Add(new_mc,-1)

        new_invCov = TMatrix_to_Numpy(dataHist.GetTotalErrorMatrix(includeStatError,errorAsFraction,useOnlyShapeErrors))[1:-1,1:-1]
        new_invCov = new_invCov - TMatrix_to_Numpy(dataHist.GetSysErrorMatrix("Flux"))[1:-1,1:-1]
        new_invCov = np.linalg.inv(new_invCov)
    else:
        new_invCov = invCov

    if setHists:
        weights = mcHist.GetCVHistoWithStatError()
        for i in range(1,weights.GetNbinsX()+1):
            weight = weights.GetBinContent(i) / new_cv[i-1] if new_cv[i-1] != 0 else weights.GetBinContent(i)
            weights.SetBinContent(i,weight)
            weights.SetBinError(i,0)

        new_mc = mcHist.Clone()
        new_mc.DivideSingle(new_mc,weights)

        if useOsc:
            histogram.SetOscHistogram(new_mc)
        else:
            histogram.SetMCHistogram(new_mc)

    return(new_cv,new_invCov,penalty)

def Chi2DataMC(histogram,marginalize=False,fluxSolution=None,useOsc=False,usePseudo=False,invCov=None,setHists=False,remakeCov=False,useNewUniverses=False,exclude="",lam=1):
    ##### Get histograms to calculate chi2 between #####
    if useOsc:
        mcHist = histogram.GetOscillatedHistogram()
    else:
        mcHist = histogram.GetMCHistogram()

    if usePseudo:
        dataHist = histogram.GetPseudoHistogram()
    else:
        dataHist = histogram.GetDataHistogram()

    #We get the number of bins and make sure it's compatible with the NxN matrix given
    if dataHist.GetNbinsX() != mcHist.GetNbinsX():
        logging.error("breaking error in Chi2DataMC")
        logging.error("The number of bins from Data ({}) and MC ({}) histograms differ. Returning -1.".format(dataHist.GetNbinsX(),mcHist.GetNbinsX()))
        return(-1)

    mc = np.array(mcHist)[1:-1] # store MC bin contents excluding over/underflow bins
    data = np.array(dataHist)[1:-1] # store data bin contents excluding over/underflow bins 
    
    # Do we want to marginalize over the flux systematic before calculating chi2
    if marginalize:
        mc,invCov,penalty = MarginalizeFlux(histogram,usePseudo=usePseudo,fluxSolution=fluxSolution,invCov=invCov,useOsc=useOsc,setHists=setHists,remakeCov=remakeCov,exclude=exclude,lam=lam)
    else:
        penalty = 0

    # ----- Get invCovariance matrix for chi2 calculation ----- #
    if invCov is None:
        useOnlyShapeErrors = False
        includeStatError   = True
        errorAsFraction    = False

        h_test = dataHist.Clone()
        h_test.Add(mcHist,-1)
        toInvCov = TMatrix_to_Numpy(h_test.GetTotalErrorMatrix(includeStatError,errorAsFraction,useOnlyShapeErrors))[1:-1,1:-1]
        invCov = np.linalg.inv(toInvCov)

    # ----- Calculate chi2 value ----= #
    diff = data - mc
    chi2 = diff.T @ invCov @ diff + penalty # @ is numpy efficient matrix multiplication

    if abs(chi2) > 1e30:
        logging.error("chi2 has invalid value: {}".format(chi2))
        print("chi2 has invalid value: {}".format(chi2))
        return(-1)

    return(chi2,penalty)

def InvertID(hist):
    for i in range(hist.GetNbinsX()+1):
        if hist.GetBinContent(i) == 0:
            hist.SetBinContent(i,1.0)
            hist.SetBinError(i,0.0)
        elif hist.GetBinContent(i) == 1:
            hist.SetBinContent(i,0.0)
            hist.SetBinError(i,0.0)

def OscillateSubHistogram(histogram,name,m,U_e4,U_mu4,U_tau4):
    # elastic_id is 1 only for the nueel samples
    # ratio_id is 1 only for the ratio samples
    # if 0 is False, if 1 is True
    hist_nueTemp = histogram.nue_templates[name]
    hist_numuTemp = histogram.numu_templates[name]
    hist_swapTemp = histogram.swap_templates[name]

    hist_nue_energy = histogram.nue_hists[name].Clone()
    hist_numu_energy = histogram.numu_hists[name].Clone()
    hist_swap_energy = histogram.swap_hists[name].Clone()

    hist = histogram.nue_hists[name].Clone()

    nue_weights = hist.GetCVHistoWithStatError().Clone()
    numu_weights = hist.GetCVHistoWithStatError().Clone()
    nuenutau_weights = hist.GetCVHistoWithStatError().Clone()
    numunue_weights = hist.GetCVHistoWithStatError().Clone()
    numunutau_weights = hist.GetCVHistoWithStatError().Clone()

    for i in range(0,hist.GetNbinsX() + 1):
        nue_sin = sin_average(i,m,hist_nueTemp,False)
        numu_sin = sin_average(i,m,hist_numuTemp,False)
        swap_sin = sin_average(i,m,hist_swapTemp,False)

        P_ee = float(1 - 4*U_e4*(1-U_e4)*nue_sin)
        
        if ('elastic' in name):
            scale = 1
            P_mue = float(4*(U_e4)*(U_mu4)*swap_sin) * scale
        else:
            P_mue = float(4*(U_e4)*(U_mu4)*swap_sin)

        P_mumu = float(1 - 4*U_mu4*(1-U_mu4)*numu_sin)

        P_mutau = float(4*U_tau4*U_mu4*numu_sin)

        P_etau = float(4*U_e4*U_tau4*nue_sin)

        nue_weights.SetBinContent(i,P_ee)
        numu_weights.SetBinContent(i,P_mumu)
        nuenutau_weights.SetBinContent(i,P_etau)
        numunue_weights.SetBinContent(i,P_mue)
        numunutau_weights.SetBinContent(i,P_mutau)

        nue_weights.SetBinError(i,0)
        numu_weights.SetBinError(i,0)
        nuenutau_weights.SetBinError(i,0)
        numunue_weights.SetBinError(i,0)
        numunutau_weights.SetBinError(i,0)

    numu = hist_numu_energy.Clone()
    numu.MultiplySingle(numu,numu_weights)

    nue = hist_nue_energy.Clone()
    nue.MultiplySingle(nue,nue_weights)

    numunue = hist_swap_energy.Clone()
    numunue.MultiplySingle(numunue,numunue_weights)

    numunutau = hist_numu_energy.Clone()
    numunutau.MultiplySingle(numunutau,numunutau_weights)

    nuenutau  = hist_nue_energy.Clone()
    nuenutau.MultiplySingle(nuenutau,nuenutau_weights)
    nutau = numunutau.Clone()
    nutau.Add(nuenutau)

    total = nue.Clone()
    total.Reset()

    nue.SetFillColor(ROOT.kRed)
    numu.SetFillColor(ROOT.kBlue)
    numunue.SetFillColor(ROOT.kBlue)
    nutau.SetFillColor(ROOT.kGray)

    nue.SetLineColor(ROOT.kRed)
    numu.SetLineColor(ROOT.kBlue)
    numunue.SetLineColor(ROOT.kBlue)
    nutau.SetLineColor(ROOT.kGray)

    nue.SetLineWidth(0)
    numu.SetLineWidth(0)
    numunue.SetLineWidth(0)
    nutau.SetLineWidth(0)
    histogram.data_hists[name].SetLineWidth(1)

    nue.SetFillStyle(3244)
    numu.SetFillStyle(3744)
    numunue.SetFillStyle(3244)
    nutau.SetFillStyle(3409)

    oscHists = []
    nue.SetTitle("#nu_{e}")
    numu.SetTitle("#nu_{#mu}")
    numunue.SetTitle("#nu_{#mu}#rightarrow #nu_{e}")
    nutau.SetTitle("#nu_{#mu}#rightarrow #nu_{#tau}")
    histogram.data_hists[name].SetTitle("Oscillated {}".format(name))

    if 'elastic' in name:
        oscHists.append(nue)
        oscHists.append(numunue)
        oscHists.append(numu)
        oscHists.append(nutau)
        total.Add(numu)
        total.Add(numunue)
        total.Add(nutau)
        total.Add(nue)
    elif 'imd' in name or 'numu' in name:
        oscHists.append(numu)
        total.Add(numu)
    elif 'ratio' in name:
        oscHists.append(nue)
        oscHists.append(numunue)
        oscHists.append(numu)
        total.Add(nue)
        total.Add(numunue)
        total.Add(numu)
    elif 'nue' in name:
        oscHists.append(nue)
        oscHists.append(numunue)
        total.Add(nue)
        total.Add(numunue)
    else:
        print("No histogram added for {}".format(name))

    return(oscHists,total)

def OscillateHistogram(histogram, m, U_e4, U_mu4, U_tau4,fitPseudodata=False,fitFluxUniverses=False):
    # elastic_id is 1 only for the nueel samples
    # ratio_id is 1 only for the ratio samples
    # if 0 is False, if 1 is True

    hist_nueTemp = histogram.nue_template
    hist_numuTemp = histogram.numu_template
    hist_swapTemp = histogram.swap_template

    hist_nue_energy = histogram.nue_hist.Clone()
    hist_numu_energy = histogram.numu_hist.Clone()
    hist_swap_energy = histogram.swap_hist.Clone()
    hist_nutau_id = histogram.elastic_id.Clone()
    hist_ratio_id = histogram.ratio_id.Clone()
    hist_inv_ratio_id = hist_ratio_id.Clone()

    InvertID(hist_inv_ratio_id)

    if fitPseudodata:
        hist = histogram.GetPseudoHistogram()
    else:
        hist = histogram.GetMCHistogram()

    nue_weights = hist.GetCVHistoWithStatError().Clone()
    numu_weights = hist.GetCVHistoWithStatError().Clone()
    nuenutau_weights = hist.GetCVHistoWithStatError().Clone()
    numunue_weights = hist.GetCVHistoWithStatError().Clone()
    numunutau_weights = hist.GetCVHistoWithStatError().Clone()

    for i in range(0,hist.GetNbinsX() + 1):
        nue_sin = sin_average(i,m,hist_nueTemp)
        numu_sin = sin_average(i,m,hist_numuTemp)
        swap_sin = sin_average(i,m,hist_swapTemp)

        P_ee = float(1 - 4*U_e4*(1-U_e4)*nue_sin)

        P_mue = float(4*(U_e4)*(U_mu4)*swap_sin)

        P_mumu = float(1 - 4*U_mu4*(1-U_mu4)*numu_sin)

        P_mutau = float(4*U_tau4*U_mu4*numu_sin)

        P_etau = float(4*U_e4*U_tau4*nue_sin)

        nue_weights.SetBinContent(i,P_ee)
        numu_weights.SetBinContent(i,P_mumu)
        nuenutau_weights.SetBinContent(i,P_etau)
        numunue_weights.SetBinContent(i,P_mue)
        numunutau_weights.SetBinContent(i,P_mutau)

        nue_weights.SetBinError(i,0)
        numu_weights.SetBinError(i,0)
        nuenutau_weights.SetBinError(i,0)
        numunue_weights.SetBinError(i,0)
        numunutau_weights.SetBinError(i,0)
        hist_ratio_id.SetBinError(i,0)
        hist_nutau_id.SetBinError(i,0)

    numu = hist_numu_energy.Clone()
    numu.MultiplySingle(numu,numu_weights)

    nue = hist_nue_energy.Clone()
    nue.MultiplySingle(nue,nue_weights)

    numunue = hist_swap_energy.Clone()
    numunue.MultiplySingle(numunue,numunue_weights)
    nue.Add(numunue)

    numunutau = hist_numu_energy.Clone()
    numunutau.MultiplySingle(numunutau,numunutau_weights)

    nuenutau  = hist_nue_energy.Clone()
    nuenutau.MultiplySingle(nuenutau,nuenutau_weights)
    nutau = numunutau.Clone()
    nutau.Add(nuenutau)
    nutau.MultiplySingle(nutau,hist_nutau_id)

    ratio = numu.Clone()
    ratio.Divide(ratio,nue)
    ratio.MultiplySingle(ratio,hist_ratio_id)

    nue.MultiplySingle(nue,hist_inv_ratio_id)
    numu.MultiplySingle(numu,hist_inv_ratio_id)

    osc = numu.Clone()
    osc.Add(nue)
    osc.Add(ratio)
    if nutau.Integral() > 0:
        osc.Add(nutau)

    histogram.SetOscHistogram(osc)

def sin_average(q=0,dm2=0,template=None,yaxis=True):
    avgsin = 0
    total_N = 0

    if yaxis:
        length = template.GetNbinsY()+1
        axis = template.GetYaxis()
    else:
        length = template.GetNbinsX()+1
        axis = template.GetXaxis()

    for b in range(length):
        lowEdge = axis.GetBinLowEdge(b)
        upEdge = axis.GetBinUpEdge(b)
        bin_width = upEdge-lowEdge
        bin_center = (upEdge+lowEdge)/2
        N_bin = template.GetBinContent(q,b)
        total_N+=N_bin
        if N_bin == 0:
            continue
        nue_sin = np.sin(1.27*dm2*bin_center)**2
        avgsin += nue_sin * N_bin
    if total_N != 0:
        avgsin = avgsin / total_N
    return(avgsin)

import os
import time
import logging, sys
import argparse
import math
import psutil
from array import array

#os.environ["OPENBLAS_NUM_THREADS"] = "1"
#os.environ["MKL_NUM_THREADS"] = "1"
#os.environ["NUMEXPR_NUM_THREADS"] = "1"
#os.environ["OMP_NUM_THREADS"] = "1"
#os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
import numpy as np
from scipy import optimize,linalg

ccnueroot = os.environ.get('CCNUEROOT')

import ROOT
from root_numpy import matrix
import PlotUtils
#insert path for modules of this package.
from tools.PlotLibrary import HistHolder
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS 

from Tools.OscHistogram import *
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

def constraint(t):
    U_e4  = t[1]
    U_mu4 = t[2]
    U_tau4= t[3]
    return 1-(U_e4 + U_mu4 + U_tau4)

class MinimizeStopper(object):
    def __init__(self):
        self.start = time.time()
    def __call__(self, xk=None):
        elapsed = time.time() - self.start
        print("Elapsed: %.3f sec" % elapsed)

def DoFit(histogram):
    x0 = [0.0,0.0,0.0,0.0]
    bounds = np.array([[0.0,1.0],[0.0,0.15],[0.0,0.41],[0,0.66]], dtype = float)
    cons = [{"type":"ineq","fun":constraint}]

    null = optimize.minimize(fun=CalChi2,x0=x0,tol=1e-4,options={"maxiter":20},args=(histogram,),method="SLSQP",bounds=bounds,constraints=optimize.LinearConstraint([[0,1,1,1]],-np.inf,1))
    null_chi2 = float(null.fun)
    print("fit near null: {}".format(null.fun))

    res = optimize.differential_evolution(func=CalChi2,bounds=bounds,polish=False,x0=x0,args=(histogram,),maxiter=50,disp=True,constraints=optimize.LinearConstraint([[0,1,1,1]],-np.inf,1))
    new_x0 = res.x
    print("best fit: {}".format(res.fun))

    res = optimize.minimize(fun=CalChi2,x0=new_x0,tol=1e-4,options={"maxiter":20},args=(histogram,),method="SLSQP",bounds=bounds,constraints=optimize.LinearConstraint([[0,1,1,1]],-np.inf,1))
    chi2 = float(res.fun)
    print("polish fit: {}".format(res.fun))

    if null_chi2 < chi2:
        chi2 = null_chi2
        res = null
    return(chi2,{"m":res.x[0]*100,"ue4":res.x[1],"umu4":res.x[2],"utau4":res.x[3]})

def FitFluxUniverses(histogram,plot=False,useOsc=False):
    if useOsc:
        mc_hist = histogram.GetOscillatedHistogram()
    else:
        mc_hist = histogram.GetMCHistogram()
        
    mc = np.array(mc_hist)[1:-1]
    data = np.array(histogram.GetDataHistogram())[1:-1]

    band = mc_hist.GetVertErrorBand("Flux")
    I = np.identity(band.GetNHists())
    A = np.array([np.array(band.GetHist(i))[1:-1] for i in range(band.GetNHists())])
    V = histogram.GetMCCov() + histogram.GetDataCov()
    C = data - mc
    
    L = 2 * np.dot(np.dot(A,np.linalg.inv(V)),C)
    Q = np.dot(np.dot(A,np.linalg.inv(V)),A.T)+I
    solution = np.linalg.solve(Q,L/2)
    new_cv = mc + np.dot(solution,A-np.array([mc for i in range(band.GetNHists())]))    

    weights = mc_hist.GetCVHistoWithStatError()
    for i in range(1,weights.GetNbinsX()+1):
        weight = weights.GetBinContent(i) / new_cv[i-1] if new_cv[i-1] != 0 else weights.GetBinContent(i)
        weights.SetBinContent(i,weight)
        weights.SetBinError(i,0)

    new_mc = mc_hist.Clone()
    new_mc.DivideSingle(new_mc,weights)

    #new_nue = histogram.nue_hist.Clone()
    #new_nue.DivideSingle(new_nue,weights)
    #histogram.nue_hist = new_nue

    #new_numu = histogram.numu_hist.Clone()
    #new_numu.DivideSingle(new_numu,weights)
    #histogram.numu_hist = new_numu

    #new_swap = histogram.swap_hist.Clone()
    #new_swap.DivideSingle(new_swap,weights)
    #histogram.swap_hist = new_swap

    if useOsc:
        histogram.SetOscHistogram(new_mc)
    else:
        histogram.SetMCHistogram(new_mc)
        histogram.RemoveFluxSystematic()

    if plot:
        cv_ratio = data/new_cv
        mc_ratio = data/mc

        x = np.arange(len(mc),dtype=float)

        c1 = ROOT.TCanvas()

        mg = ROOT.TMultiGraph()

        g1 = ROOT.TGraph(len(x), x, cv_ratio)
        g3 = ROOT.TGraph(len(x), x, mc_ratio)

        g1.SetLineColor(ROOT.kRed)
        g3.SetLineColor(ROOT.kBlue)

        g1.SetLineWidth(3)
        g3.SetLineWidth(3)

        g1.SetTitle("After Flux Calc.")
        g3.SetTitle("Before Flux Calc.")

        mg.Add(g1)
        mg.Add(g3)

        mg.GetYaxis().SetTitle("Data Ratio")

        mg.Draw("AL")
        c1.BuildLegend()

        c1.Print("flux_fit_results.png")

    return(solution)

def fitNorm(hist_data, templates, hist_mc):
    x0 = [1,1,1,1]
    mc = hist_mc.Clone()
    data = hist_data.Clone()

    res = optimize.minimize(NormCalc, x0,tol=0.,method="Nelder-Mead",args=(mc,templates,data),options={"disp":True})
    chi2 = float(res.fun)
    bestfit = MuonNorm(hist_mc,templates,res.x[0],res.x[1],res.x[2],res.x[3])
    return(chi2,bestfit,{"fhc numu":res.x[0],"rhc numu":res.x[1],"fhc nue":res.x[2],"rhc nue":res.x[3]})

def NormCalc(x,hist_mc,templates,hist_data):
    fhc_norm = x[0]
    rhc_norm = x[1]
    fhc_nue_norm = x[2]
    rhc_nue_norm = x[3]
    normalized = MuonNorm(hist_mc,templates,fhc_norm,rhc_norm,fhc_nue_norm,rhc_nue_norm)
    PlotNorms(hist_mc,hist_data,normalized,templates,[fhc_norm,rhc_norm,fhc_nue_norm,rhc_nue_norm])
    chi2 = Chi2DataMC(hist_data,normalized)
    return(chi2)

def PlotNorms(hist_mc,hist_data,fitHist,templates,params=[]):
    chi2_model = Chi2DataMC(hist_data,fitHist)
    chi2_null = Chi2DataMC(hist_data,hist_mc)

    c1 = ROOT.TCanvas()
    margin = .12
    bottomFraction = .2
    overall = ROOT.TCanvas("Data/MC")
    top = ROOT.TPad("DATAMC", "DATAMC", 0, bottomFraction, 1, 1)
    bottom = ROOT.TPad("Ratio", "Ratio", 0, 0, 1, bottomFraction+margin)

    top.Draw()
    bottom.Draw()

    top.cd()
    top.SetLogy()
    
    fitHist.GetXaxis().SetTitle("Bin number")
    fitHist.GetYaxis().SetTitle("Entries")

    nullRatio =  hist_data.Clone()
    oscRatio =  fitHist.Clone()

    fitHist.SetLineColor(ROOT.kRed)
    fitHist.SetLineWidth(3)
    hist_mc.SetLineColor(ROOT.kBlue)
    hist_mc.SetLineWidth(3)
    hist_mc.GetYaxis().SetTitle("Nevents")
    hist_mc.Draw("hist")
    fitHist.Draw("hist same")
    hist_data.Draw("same")

    osc = fitHist.GetCVHistoWithError()
    osc.SetLineColor(ROOT.kRed)
    osc.SetLineWidth(3)
    osc.SetMarkerStyle(0)
    osc.SetFillColorAlpha(ROOT.kPink + 1, 0.3)
    osc.Draw("E2 SAME")

    null = hist_mc.GetCVHistoWithError()
    null.SetLineColor(ROOT.kBlue)
    null.SetLineWidth(3)
    null.SetMarkerStyle(0)
    null.SetFillColorAlpha(ROOT.kBlue + 1, 0.3)
    null.Draw("E2 SAME")

    leg = ROOT.TLegend()
    if len(params) > 0:
        leg.SetHeader("fhc muon:{:.2f} rhc muon:{:.2f} fhc elec:{:.2f} rhc elec:{:.2f}".format(params[0],params[1],params[2],params[3]))
    leg.AddEntry(hist_data,"Data","p")
    leg.AddEntry(fitHist,"Best Fit #chi^{2}="+"{:.1f}".format(chi2_model),"l")
    #leg.AddEntry(fitHist,"Oscillation #chi^{2}="+"{:.1f}".format(chi2_model),"l")
    #leg.AddEntry(fitHist,"RAA #chi^{2}="+"{:.1f}".format(chi2_model),"l")
    leg.AddEntry(hist_mc,"Null Hypothesis #chi^{2}="+"{:.1f}".format(chi2_null),"l")
    leg.Draw()

    oscRatio.Divide(oscRatio, hist_mc)
    nullRatio.Divide(nullRatio,hist_mc)

    bottom.cd()
    bottom.SetTopMargin(0)
    bottom.SetBottomMargin(0.3)

    nullErrors = hist_mc.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    for whichBin in range(0, nullErrors.GetXaxis().GetNbins()+1): 
        nullErrors.SetBinError(whichBin, max(nullErrors.GetBinContent(whichBin), 1e-9))
        nullErrors.SetBinContent(whichBin, 1)

    oscRatio.SetTitle("")
    oscRatio.SetLineColor(ROOT.kRed)
    nullRatio.SetLineColor(ROOT.kBlue)
    oscRatio.SetLineWidth(3)
    nullRatio.SetLineWidth(3)
    oscRatio.SetTitleSize(0)

    #Error envelope for the MC
    nullErrors.SetLineWidth(0)
    nullErrors.SetMarkerStyle(0)
    nullErrors.SetFillColorAlpha(ROOT.kBlue + 1, 0.4)
    nullErrors.Draw("E2")

    nullErrors.GetYaxis().SetTitle("#splitline{Ratio to}{Null Hypothesis}")
    nullErrors.GetYaxis().SetLabelSize(.13)
    nullErrors.GetYaxis().SetTitleSize(0.1)
    nullErrors.GetYaxis().SetTitleOffset(0.6)
    nullErrors.GetYaxis().SetNdivisions(505) #5 minor divisions between 5 major divisions.  I'm trying to match a specific paper here.
    nullErrors.GetXaxis().SetTitleSize(0.16)
    nullErrors.GetXaxis().SetTitleOffset(0.9)
    nullErrors.GetXaxis().SetLabelSize(.15)
    nullErrors.GetXaxis().SetTitle("Bin Number")
    nullErrors.SetMinimum(0.5)
    nullErrors.SetMaximum(1.5)

    
    #Draw the data ratios
    oscRatio.SetMinimum(0.5)
    oscRatio.SetMaximum(1.5)
    nullRatio.SetMinimum(0.5)
    nullRatio.SetMaximum(1.5)
    nullRatio.SetLineColorAlpha(ROOT.kBlue+1,0.6)
    oscRatio.Draw("same")
    nullRatio.Draw("same")

    #Draw a flat line at 1 for oscRatio of MC to itself
    straightLine = nullErrors.Clone()
    straightLine.SetLineColor(ROOT.kBlue)
    straightLine.SetLineWidth(3)
    straightLine.SetFillColor(0)
    straightLine.Draw("HIST L SAME")

    leg1 = ROOT.TLegend(.5,.5,.9,.9)
    leg1.AddEntry(nullRatio,"Data/Null Hypothesis","p")
    leg1.AddEntry(oscRatio,"Model/Null Hypothesis","l")
    leg1.AddEntry(straightLine,"Null/Null Hypothesis","l")
    #leg1.Draw()

    top.cd()
    if len(params) > 0:
        overall.Print("fit_norm_{:.2f}.png".format(chi2_model))
    else:
        overall.Print("fit_muonnorm.png")

def MuonNorm(hist_mc,templates,fhc_norm,rhc_norm,fhc_nue_norm,rhc_nue_norm):
    hist_nueTemp = templates["nue_template"]
    hist_numuTemp = templates["numu_template"]
    hist_swapTemp = templates["swap_template"]

    hist_nue_energy = templates["nue"].Clone()
    hist_numu_energy = templates["numu"].Clone()
    hist_nutau_id = templates["elastic_id"].Clone()
    hist_ratio_id = templates["ratio_id"].Clone()
    hist_fhc_id = templates["beam_id"].Clone()
    hist_rhc_id = hist_fhc_id.Clone()
    hist_inv_ratio_id = hist_ratio_id.Clone()
    hist_swap_energy = templates["swap"].Clone()

    InvertID(hist_inv_ratio_id)
    InvertID(hist_rhc_id)

    fhc_numu_energy = hist_numu_energy.Clone()
    rhc_numu_energy = hist_numu_energy.Clone()
    fhc_nue_energy = hist_nue_energy.Clone()
    rhc_nue_energy = hist_nue_energy.Clone()

    fhc_numu_energy.MultiplySingle(fhc_numu_energy,hist_fhc_id)
    rhc_numu_energy.MultiplySingle(rhc_numu_energy,hist_rhc_id)
    fhc_nue_energy.MultiplySingle(fhc_nue_energy,hist_fhc_id)
    rhc_nue_energy.MultiplySingle(rhc_nue_energy,hist_rhc_id)

    fhc_numu_energy.Scale(1/fhc_norm)
    rhc_numu_energy.Scale(1/rhc_norm)
    fhc_nue_energy.Scale(1/fhc_nue_norm)
    rhc_nue_energy.Scale(1/rhc_nue_norm)

    numu = fhc_numu_energy.Clone()
    numu.Add(rhc_numu_energy)
    nue = fhc_nue_energy.Clone()
    nue.Add(rhc_nue_energy)

    ratio = numu.Clone()
    ratio.Divide(ratio,nue)
    ratio.MultiplySingle(ratio,hist_ratio_id)

    nue.MultiplySingle(nue,hist_inv_ratio_id)
    numu.MultiplySingle(numu,hist_inv_ratio_id)

    hist = numu.Clone()
    hist.Add(nue)
    hist.Add(ratio)
    return(hist)

def Chi2DataMC(dataHist,mcHist,cov=None):
    #We get the number of bins and make sure it's compatible with the NxN matrix given
    if dataHist.GetNbinsX() != mcHist.GetNbinsX():
        logging.error("breaking error in Chi2DataMC")
        logging.error("The number of bins from Data ({}) and MC ({}) histograms differ. Returning -1.".format(dataHist.GetNbinsX(),mcHist.GetNbinsX()))
        return(-1)

    #get the covariance matrix
    useOnlyShapeErrors = False
    includeStatError   = True
    errorAsFraction    = False

    # ----- Get covariance matrix for chi2 calculation ----- #
    if cov is None:
        h_test = dataHist.Clone()
        h_test.Add(mcHist,-1)
        cov = np.asarray(matrix(h_test.GetTotalErrorMatrix(includeStatError,errorAsFraction,useOnlyShapeErrors)))[1:-1,1:-1]

    try:
        errorMatrix = np.linalg.inv(cov)
    except:
        logging.error("Data-MC covariance matrix couldn't be inverted. Returning -1")
        return(-1)

    mc = np.array(mcHist)[1:-1] # store MC bin contents excluding over/underflow bins
    data = np.array(dataHist)[1:-1] # store data bin contents excluding over/underflow bins 

    # ----- Calculate chi2 value ----= #
    diff = mc - data
    chi2 = diff.T @ errorMatrix @ diff # @ is numpy efficient matrix multiplication

    if abs(chi2) > 1e30:
        logging.error("chi2 has invalid value: {}".format(chi2))
        print("chi2 has invalid value: {}".format(chi2))
        return(-1)
    

    return(chi2)

def CalChi2(x,histogram):
    ms = x[0]*100
    U_e4 = x[1]
    U_mu4 = x[2]
    U_tau4 = x[3]
    OscillateHistogram(histogram,ms,U_e4,U_mu4,U_tau4,False,False)
    chi2 = Chi2DataMC(histogram.GetDataHistogram(),histogram.GetOscillatedHistogram())
    return(chi2)

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

    hist_nueTemp = histogram.samples_nue_templates[name]
    hist_numuTemp = histogram.samples_numu_templates[name]
    hist_swapTemp = histogram.samples_swap_templates[name]

    hist_nue_energy = histogram.samples_nue[name].Clone()
    hist_numu_energy = histogram.samples_numu[name].Clone()
    hist_swap_energy = histogram.samples_swap[name].Clone()

    hist = histogram.samples_nue[name].Clone()

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
    histogram.data_samples[name].SetLineWidth(1)

    nue.SetFillStyle(3244)
    numu.SetFillStyle(3744)
    numunue.SetFillStyle(3244)
    nutau.SetFillStyle(3409)

    TArray = []
    nue.SetTitle("#nu_{e}")
    numu.SetTitle("#nu_{#mu}")
    numunue.SetTitle("#nu_{#mu}->#nu_{e}")
    nutau.SetTitle("#nu_{#tau}")
    histogram.data_samples[name].SetTitle("Oscillated {}".format(name))

    if "elastic" in name:
        nue.Scale(2,'width')
        numu.Scale(2,'width')
        numunue.Scale(2,'width')
        nutau.Scale(2,'width')
        histogram.data_samples[name].Scale(2,'width')
        histogram.data_samples[name].GetXaxis().SetTitle("Electron Energy [ / 2 GeV ]")
    else:
        nue.Scale(1,'width')
        numu.Scale(1,'width')
        numunue.Scale(1,'width')
        histogram.data_samples[name].Scale(1,'width')
        histogram.data_samples[name].GetXaxis().SetTitle("Neutrino Energy Estimator [ GeV ]")

    if nue.Integral() > 0:
        TArray.append(nue)
    if numu.Integral() > 0:
        TArray.append(numu)
    if numunue.Integral() > 0 and ("elastic" in name or "nue" in name):
        TArray.append(numunue)
    if nutau.Integral() > 0 and "elastic" in name:
        TArray.append(nutau)

    return(TArray)

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

    # If we've already fit the flux universes, just reweight the MC to maintain flux fits 
    if "Flux" not in hist.GetVertErrorBandNames():
        nue = hist_nue_energy.Clone()
        numu = hist_numu_energy.Clone()
        ratio = hist_numu_energy.Clone()
        ratio.Divide(ratio,nue)

        nue.MultiplySingle(nue,hist_inv_ratio_id)
        numu.MultiplySingle(numu,hist_inv_ratio_id)
        ratio.MultiplySingle(ratio,hist_ratio_id)

        nue.Add(numu)
        nue.Add(ratio)

        weights = nue.Clone()
        weights.DivideSingle(weights,osc)

        osc = histogram.GetMCHistogram()
        osc.DivideSingle(osc,weights)


    histogram.SetOscHistogram(osc)
    if fitFluxUniverses:
        sol = FitFluxUniverses(histogram,False,True)
        histogram.RemoveFluxfromOsc()

def makeChi2Surface(histogram,outdir,dodeltachi2=False,deltam=1,U_e4=0,U_tau4=0):
    U_mu4s = 0.41*np.logspace(-3,0,60)
    U_mu4s[0] = 0
    asimov_surface = np.zeros(np.shape(U_mu4s)[0],dtype='f')
    data_surface   = np.zeros(np.shape(U_mu4s)[0],dtype='f')
    fits           = np.zeros(np.shape(U_mu4s)[0],dtype='f')
    count = 0
    for i in range(U_mu4s.shape[0]):
        count+=1
        U_mu4 = U_mu4s[i]
        if dodeltachi2:
            OscillateHistogram(histogram, deltam, U_e4, U_mu4, U_tau4)
            test_hist = histogram.GetOscillatedHistogram()
            
            minchi2,res = DoFit(histogram)
            chi2_null = Chi2DataMC(histogram.GetDataHistogram(),histogram.GetMCHistogram())

            asimov_surface[i] = minchi2 - chi2_null
            fits[i] = minchi2
        else:
            OscillateHistogram(histogram, deltam, U_e4, U_mu4, U_tau4)
            test_hist = histogram.GetOscillatedHistogram()

            chi2_data       = Chi2DataMC(histogram.GetDataHistogram(),test_hist)
            data_surface[i] = chi2_data
            chi2_asimov       = Chi2DataMC(histogram.GetPseudoHistogram(),test_hist)
            asimov_surface[i] = chi2_asimov
        
        logging.info("{:.2f}% done with chi2s. Current chi2 = {:.4f}".format(100*count/(U_mu4s.shape[0]),data_surface[i]))

    if dodeltachi2:
        np.save('{}/deltachi2_surface_m_{}_Ue4_{}.dat'.format(outdir,deltam,U_e4),surface)
        np.save('{}/chi2_fit_chi2_m_{}_Ue4_{}.dat'.format(outdir,deltam,U_e4),fits)
    else:
        np.save('{}/chi2_surface_data_m_{}_Ue4_{}.dat'.format(outdir,deltam,U_e4),data_surface)
        np.save('{}/chi2_surface_pseudodata_m_{}_Ue4_{}.dat'.format(outdir,deltam,U_e4),asimov_surface)

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

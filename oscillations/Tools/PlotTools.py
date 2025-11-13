import os
import time
import logging, sys
import ROOT
import PlotUtils
import numpy as np
from scipy import optimize, integrate
import argparse
ccnueroot = os.environ.get('CCNUEROOT')
plotutils = os.environ.get('PLOTUTILSROOT')

import math
import psutil
import multiprocessing
import threading
nthreads = 4
from array import array
from Tools.FitTools import *
from Tools.Histogram import *
from Tools.Helper import *

#insert path for modules of this package.
from tools.PlotLibrary import HistHolder

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

MNVPLOTTER = PlotUtils.MnvPlotter()
MNVPLOTTER.draw_normalized_to_bin_width=False

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)
legend_text_size = .025

def UndoBinWidthNorm(histogram):
    for i in range(0,histogram.GetNbinsX()+1):
        width = histogram.GetBinWidth(i)
        cont = histogram.GetBinContent(i)
        new_cont = cont*width
        if new_cont != 0:
            histogram.SetBinContent(i,new_cont)

        for name in histogram.GetVertErrorBandNames():
            band = histogram.GetVertErrorBand(name)
            cont = band.GetBinContent(i)
            new_cont = cont*width
            if new_cont != 0:
                band.SetBinContent(i,new_cont)

            for j in range(band.GetNHists()):
                hist = band.GetHist(j)
                cont = hist.GetBinContent(i)
                new_cont = cont*width
                if new_cont != 0:
                    hist.SetBinContent(i,new_cont)


def PlotSampleMarginalizationEffects(sample_histogram):
    histogram = copy.deepcopy(sample_histogram)
    samplesToFit = [["elastic","imd","numu","ratio"],["elastic"],["elastic","imd"],["numu"],["elastic","imd","numu"],["elastic","imd","ratio"]]
    names = ["All","elastic","elastic_imd","numu","elastic_imd_numu","elastic_imd_ratio"]
    titles = ["#nu+e, IMD, CC #nu_{#mu}, CC #nu_{#mu}/#nu_{e}","#nu+e","#nu+e, IMD","CC #nu_{#mu}","#nu+e, IMD, CC #nu_{#mu}","#nu+e, IMD, CC #nu_{#mu}/#nu_{e}"]
    colors = [ROOT.kRed,ROOT.kBlue,ROOT.kViolet,ROOT.kOrange,ROOT.kTeal,ROOT.kGreen]
    hists = []
    flux_hists = []
    chi2s = []
    pens = []

    subSample = histogram.mc_hists["fhc_elastic"].Clone()
    band = subSample.GetVertErrorBand("Flux")
    nhists = band.GetNHists()
    universes = np.array([np.array(band.GetHist(l))[1:-1] for l in range(nhists)])
    fhc_integrals = universes.sum(axis=1)
    old_fhc = subSample.Integral()

    subSample = histogram.mc_hists["rhc_elastic"].Clone()
    band = subSample.GetVertErrorBand("Flux")
    nhists = band.GetNHists()
    universes = np.array([np.array(band.GetHist(l))[1:-1] for l in range(nhists)])
    rhc_integrals = universes.sum(axis=1)
    old_rhc = subSample.Integral()

    c0 = ROOT.TCanvas()

    mg = ROOT.TMultiGraph()
    mg.SetTitle(";#nu_{#mu}-mode Number of #nue events;#bar{#nu}_{#mu}-mode Number of #nue events")
    mg.GetHistogram().GetXaxis().SetLimits(1000,1800)
    mg.GetHistogram().GetYaxis().SetRangeUser(650,1200)

    g1 = ROOT.TGraph(len(fhc_integrals),fhc_integrals,rhc_integrals)
    ROOT.SetOwnership(g1, False)
    g1.SetTitle("Flux Universes")
    g1.SetMarkerStyle(4)
    g1.SetMarkerColorAlpha(ROOT.kBlack,1)
    g1.SetLineWidth(0)
    g2 = ROOT.TGraph()
    ROOT.SetOwnership(g2, False)
    g2.SetPoint(0,old_fhc,old_rhc)
    g2.SetTitle("Central Value")
    g2.SetMarkerStyle(41)
    g2.SetMarkerSize(3)
    g2.SetMarkerColor(ROOT.kRed)
    g2.SetLineWidth(0)
    g3 = ROOT.TGraphErrors(1,np.array([1338.0]),np.array([838.2]),np.array([99.7]),np.array([63.1]))
    ROOT.SetOwnership(g3, False)
    g3.SetTitle("Mean/RMS Before Constraint Result")
    g3.SetMarkerStyle(20)
    g3.SetMarkerColor(ROOT.kBlack)
    g4 = ROOT.TGraphErrors(1,np.array([1239.0]),np.array([778.4]),np.array([41.9]),np.array([31.2]))
    ROOT.SetOwnership(g4, False)
    g4.SetTitle("Mean/RMS After Constraint Result")
    g4.SetMarkerStyle(20)
    g4.SetMarkerColor(ROOT.kBlack)

    mg.Add(g1)
    mg.Add(g2)
    mg.Add(g3)
    mg.Add(g4)

    f_numu = ROOT.TFile.Open(plotutils+'/data/flux/flux-g4numiv6-pdg14-minervame1D1M1NWeightedAve.root')
    fhc_numu = f_numu.Get("flux_E_unweighted")
    f_numu.Close()
    f_nue = ROOT.TFile.Open(plotutils+'/data/flux/flux-g4numiv6-pdg12-minervame1D1M1NWeightedAve.root')
    fhc_nue = f_nue.Get("flux_E_unweighted")
    f_nue.Close()
    fhc_numu_univ = fhc_numu.GetVertErrorBand("Flux")
    fhc_nue_univ = fhc_nue.GetVertErrorBand("Flux")

    for i,samples in enumerate(samplesToFit):
        new_histogram = StitchedHistogram(names[i],True)

        for key in histogram.keys + ["fhc_ratio","rhc_ratio"]:
            for sample in samples:
                if sample in key:
                    if 'ratio' in sample:
                        for beam in ['fhc','rhc']:
                            h_mc_numu = histogram.mc_hists[beam+"_numu_selection"].Clone()
                            h_data_numu = histogram.data_hists[beam+"_numu_selection"].Clone()
                            h_mc_nue = histogram.mc_hists[beam+"_nue_selection"].Clone()
                            h_data_nue = histogram.data_hists[beam+"_nue_selection"].Clone()
                            h_mc_numu.Divide(h_mc_numu,h_mc_nue)
                            h_data_numu.Divide(h_data_numu,h_data_nue)
                            new_histogram.AddHistograms(beam+"_ratio",h_mc_numu.Clone(),h_data_numu.Clone())
                    else:
                        new_histogram.AddHistograms(key,histogram.mc_hists[key].Clone(),histogram.data_hists[key].Clone())

        new_histogram.Stitch()
        invCov=new_histogram.GetInverseCovarianceMatrix(sansFlux=True)

        fluxSolution,nullPen = FluxSolution(new_histogram,invCov=invCov)
        invCov = new_histogram.GetInverseCovarianceMatrix(sansFlux=True)
        plot_histogram = copy.deepcopy(histogram)
        chi2,penalty = Chi2DataMC(plot_histogram,fluxSolution=fluxSolution,invCov=invCov,setHists=True,marginalize=True)

        hist = plot_histogram.GetMCHistogram()
        hists.append(hist)
        chi2s.append(chi2)
        pens.append(penalty)

        subSample = histogram.mc_hists["fhc_elastic"].Clone()
        weights = ReweightCV(subSample,fluxSolution=fluxSolution)
        new_fhc = subSample.Integral()

        subSample = histogram.mc_hists["rhc_elastic"].Clone()
        weights = ReweightCV(subSample,fluxSolution=fluxSolution)
        new_rhc = subSample.Integral()

        new_fhc_numu = fhc_numu.Clone()
        new_fhc_nue = fhc_nue.Clone()
        weights = ReweightCV(new_fhc_numu,fluxSolution=fluxSolution)
        weights = ReweightCV(new_fhc_nue,fluxSolution=fluxSolution)

        flux_hists.append(new_fhc_nue)

        g_ = ROOT.TGraph()
        ROOT.SetOwnership(g_, False)
        g_.SetPoint(0,new_fhc,new_rhc)
        g_.SetTitle(titles[i])
        g_.SetMarkerStyle(29)
        g_.SetMarkerSize(3)
        g_.SetLineWidth(0)
        g_.SetMarkerColorAlpha(colors[i],0.4)
        mg.Add(g_)

    mg.Draw("AP")
    pad = ROOT.gPad
    pad.BuildLegend()
    c0.Print("plots/integrated_elastic_events.png")

    c1 = ROOT.TCanvas("C2", "canvas2", 1024, 640)
    c1.cd(1)
    h_null =  histogram.GetMCHistogram()
    h_data = histogram.GetDataHistogram()
    invCov = histogram.GetInverseCovarianceMatrix(sansFlux=True)
    chi2,pen = Chi2DataMC(histogram,invCov=invCov)

    nullRatio = h_data.Clone()
    nullRatio.Divide(nullRatio,h_null)

    nullErrors = h_null.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    for whichBin in range(0, nullErrors.GetXaxis().GetNbins()+1): 
        nullErrors.SetBinError(whichBin, max(nullErrors.GetBinContent(whichBin), 1e-9))
        nullErrors.SetBinContent(whichBin, 1)
    nullErrors.SetFillColorAlpha(ROOT.kPink + 1, 0.4)

    straightLine = nullErrors.Clone()
    straightLine.SetFillStyle(0)

    RatioAxis(nullErrors,MNVPLOTTER)
    nullErrors.SetMinimum(0)
    nullErrors.SetMaximum(2)
    nullErrors.GetYaxis().SetTitle("#splitline{Ratio to Null}{Hypothesis}")
    nullErrors.GetXaxis().SetTitleOffset(1.5)
    nullErrors.Draw("E2")
    nullRatio.Draw("SAME")

    c1.SetTopMargin(0.35)
    c1.SetRightMargin(0.05)
    leg = ROOT.TLegend(0.05, 0.675, 0.95, .975)
    leg.SetTextSize(legend_text_size);
    leg.SetNColumns(len(titles)//2) # // median N number of rows per column
    top = "Data"
    bottom = "#chi^{2}="+"{:.2f}".format(chi2)
    legend = "#splitline{%s}{%s}" % (top,bottom)
    leg.AddEntry(nullRatio,legend,"p")

    for i,hist in enumerate(hists):
        hist.Divide(hist,h_null)
        hist.SetLineColor(colors[i])
        hist.SetFillStyle(0)
        if i == 0:
            hist.SetLineStyle(2)
        hist.Draw("same hist l")

        top = titles[i]
        bottom = "#chi^{2}="+"{:.2f}+".format(chi2s[i]-pens[i])+"{:.2f} pen.".format(pens[i])
        legend = "#splitline{%s}{%s}" % (top,bottom)
        leg.AddEntry(hist,legend,"l")

    straightLine.Draw("same hist") 

    leg.Draw()

    c1.Print("plots/stitched_flux_marg_effects.png")

    new_bins = array('d',list(range(0,21)))
    UndoBinWidthNorm(fhc_nue)
    UndoBinWidthNorm(fhc_numu)

    fhc_numu = fhc_nue.Rebin(20,"hnew",new_bins)
    fhc_numu.Scale(1,"width")
    fhc_nue = fhc_nue.Rebin(20,"hnew",new_bins)
    fhc_nue.Scale(1,"width")

    for i,hist in enumerate(flux_hists):
        UndoBinWidthNorm(flux_hists[i])
        flux_hists[i] = hist.Rebin(20,str(i),new_bins)
        flux_hists[i].Scale(1,'width')

    fhc_numu.GetXaxis().SetRangeUser(0,20)
    fhc_numu.GetXaxis().SetTitle("Neutrino Energy")
    fhc_numu.SetTitle("FHC #nu_{#mu} Flux Prediction")

    fhc_nue.GetXaxis().SetRangeUser(0,20)
    fhc_nue.GetXaxis().SetTitle("Neutrino Energy")
    fhc_nue.SetTitle("FHC #nu_{e} Flux Prediction")

    PlotWithRatio(MNVPLOTTER,"plots/FluxReweighting.png",fhc_nue,hists=flux_hists,titles=titles,colors=colors)

def PlotFluxMarginalizationEffects(sample_histogram,parameters,name="",plotSamples=False,usePseudo=False):
    histogram = copy.deepcopy(sample_histogram)
    histogram.SetPlottingStyle()
    
    invCov=histogram.GetInverseCovarianceMatrix(sansFlux=True)

    nullSolution,nullPen = FluxSolution(histogram,invCov=invCov,usePseudo=usePseudo)
    oscSolution,oscPen   = FluxSolution(histogram,invCov=invCov,useOsc=True,usePseudo=usePseudo)

    nullSols = ROOT.TH1D("null_solutions","Null Solutions",20,-1,1)
    oscSols = ROOT.TH1D("osc_solutions","Best Fit Solutions",20,-1,1)
    for i in range(len(nullSolution)):
        nullSols.Fill(nullSolution[i])
        oscSols.Fill(oscSolution[i])

    nullSolution = nullSolution/np.max(nullSolution)
    oscSolution = oscSolution/np.max(oscSolution)

    chi2_null_marg,null_pen = Chi2DataMC(histogram,invCov=invCov,marginalize=True,usePseudo=usePseudo)
    chi2_model_marg,model_pen = Chi2DataMC(histogram,invCov=invCov,marginalize=True,useOsc=True,usePseudo=usePseudo)
    chi2_null,_ = Chi2DataMC(histogram,invCov=invCov,marginalize=False,usePseudo=usePseudo)
    chi2_model,_ = Chi2DataMC(histogram,invCov=invCov,marginalize=False,useOsc=True,usePseudo=usePseudo)
    h_data = histogram.GetPseudoHistogram() if usePseudo else histogram.GetDataHistogram()

    h_osc = histogram.GetOscillatedHistogram()
    _,_ = Chi2DataMC(histogram,invCov=invCov,useOsc=True,setHists=True,marginalize=True,usePseudo=usePseudo)
    h_osc_marg = histogram.GetOscillatedHistogram()
    h_null = histogram.GetMCHistogram()
    _,_ = Chi2DataMC(histogram,invCov=invCov,setHists=True,marginalize=True,usePseudo=usePseudo)
    h_null_marg = histogram.GetMCHistogram()

    c1 = ROOT.TCanvas("CMarg", "canvas", 1024, 640)
    c1.SetFillStyle(4000)

    lMargin = 0.12
    rMargin = 0.05
    bMargin = 0.15
    tMargin = 0.05

    c1.Divide(1,3)

    # ------------------------------------ First Pad -----------------------------------------
    c1.cd(1)
    pad = ROOT.gPad
    pad.SetMargin(lMargin, rMargin, bMargin, tMargin)

    nullSols.GetYaxis().SetTitle("Frequency")
    nullSols.GetXaxis().SetTitle("Marginalization Solution Pulls")
    RatioAxis(nullSols,MNVPLOTTER)

    nullSols.SetLineWidth(2)
    nullSols.SetLineColor(ROOT.kRed)
    oscSols.SetLineWidth(2)
    oscSols.SetLineColor(ROOT.kBlue)

    nullSols.Draw("hist")
    oscSols.Draw("hist same")

    leg1 = ROOT.TLegend(.3,.2)
    leg1.AddEntry(nullSols,"Null Hypothesis Pulls","l")
    leg1.AddEntry(oscSols,"Osc. Model Pulls","l")
    leg1.Draw()

    # ------------------------------------ Second Pad -----------------------------------------
    # plot mc - flux universe / mc

    c1.cd(2)
    pad = ROOT.gPad
    pad.SetMargin(lMargin, rMargin, bMargin, tMargin)

    straightLine = h_null.Clone()
    for whichBin in range(0, straightLine.GetXaxis().GetNbins()+1): 
        straightLine.SetBinContent(whichBin, 0)

    RatioAxis(straightLine,MNVPLOTTER)
    straightLine.SetMinimum(-.5)
    straightLine.SetMaximum(.5)
    straightLine.GetYaxis().SetTitle("#frac{Universe - CV}{CV}")

    straightLine.Draw("hist l")

    cv = h_null.GetCVHistoWithError()
    h_test = h_null.Clone()
    errband = h_test.GetVertErrorBand("Flux")
    for i in range(errband.GetNHists()):
        universe = errband.GetHist(i)
        universe.Add(cv,-1)
        universe.Divide(universe,cv)
        universe.SetMarkerSize(0)
        universe.DrawClone("same hist l")

    straightLine.Draw("hist l same")

    leg2 = ROOT.TLegend(.3,.2)
    leg2.AddEntry(straightLine,"Central Value","l")
    #leg2.Draw()

    # ------------------------------------ Third Pad -----------------------------------------
    # plot solution *  [mc - flux universe / mc ]

    c1.cd(3)
    pad = ROOT.gPad
    pad.SetMargin(lMargin, rMargin, bMargin, tMargin)

    line = straightLine.Clone()
    line.GetYaxis().SetTitle("#frac{solution#times(Universe - CV)}{CV}")
    line.Draw("hist l")

    cv = h_null.GetCVHistoWithError()
    errband = h_null.GetVertErrorBand("Flux")
    for i in range(errband.GetNHists()):
        universe = errband.GetHist(i)
        universe.Add(cv,-1)
        universe.Scale(nullSolution[i])
        universe.Divide(universe,cv)
        universe.SetMarkerSize(0)
        universe.DrawClone("same hist l")

    line.Draw("hist l same")

    leg3 = ROOT.TLegend(.3,.2)
    leg3.AddEntry(straightLine,"Central Value","l")
    leg3.Draw()

    c1.Print("plots/{}_marg_effects_{:.1f}_{:.3f}_{:.4f}.png".format(name,parameters["m"],parameters['ue4'],parameters['umu4']))

def PlotOscillationRatios(sample_histogram,parameters,name="",plotSamples=False,usePseudo=False):
    histogram = copy.deepcopy(sample_histogram)
    histogram.SetPlottingStyle()
    
    invCov=histogram.GetInverseCovarianceMatrix(sansFlux=True)

    nullSolution,nullPen = FluxSolution(histogram,invCov=invCov,usePseudo=usePseudo)
    oscSolution,oscPen   = FluxSolution(histogram,invCov=invCov,useOsc=True,usePseudo=usePseudo)

    nullSols = ROOT.TH1D("null_solutions","Null Solutions",20,-.5,.5)
    oscSols = ROOT.TH1D("osc_solutions","Best Fit Solutions",20,-.5,.5)
    for i in range(len(nullSolution)):
        nullSols.Fill(nullSolution[i])
        oscSols.Fill(oscSolution[i])

    h_null = histogram.GetMCHistogram()
    h_osc = histogram.GetOscillatedHistogram()
    h_data = histogram.GetPseudoHistogram() if usePseudo else histogram.GetDataHistogram()

    chi2_null,_ = Chi2DataMC(histogram,invCov=invCov,marginalize=False,usePseudo=usePseudo)
    chi2_model,_ = Chi2DataMC(histogram,invCov=invCov,marginalize=False,useOsc=True,usePseudo=usePseudo)

    chi2_model_marg,model_pen = Chi2DataMC(histogram,invCov=invCov,marginalize=True,useOsc=True,usePseudo=usePseudo,setHists=True)
    chi2_null_marg,null_pen = Chi2DataMC(histogram,invCov=invCov,marginalize=True,usePseudo=usePseudo,setHists=True)

    h_null_marg = histogram.GetMCHistogram()
    h_osc_marg = histogram.GetOscillatedHistogram()

    c1 = ROOT.TCanvas("C", "canvas", 1024, 640)
    c1.SetFillStyle(4000)

    Nx = 1
    Ny = 3

    lMargin = 0.12
    rMargin = 0.05
    bMargin = 0.15
    tMargin = 0.05

    CanvasPartition(c1, Nx, Ny, lMargin, rMargin, bMargin, tMargin)

    # ------------------------------------ First Pad -----------------------------------------
    pad1 = c1.FindObject("pad_{}_{}".format(0,2))
    pad1.cd()

    nullRatio = h_data.Clone()
    oscRatio  = h_osc.Clone()
    margRatio = h_null_marg.Clone()
    flux_null = h_null.Clone()

    nullRatio.Divide(nullRatio,h_null)
    oscRatio.Divide(oscRatio, h_null)
    margRatio.Divide(margRatio, h_null)

    for err in flux_null.GetVertErrorBandNames():
        if err != "Flux":
            flux_null.PopVertErrorBand(err)

    fluxErrors = flux_null.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    for whichBin in range(0, fluxErrors.GetXaxis().GetNbins()+1): 
        fluxErrors.SetBinError(whichBin, max(fluxErrors.GetBinContent(whichBin), 1e-9))
        fluxErrors.SetBinContent(whichBin, 1)
    fluxErrors.SetFillColorAlpha(ROOT.kPink + 1, 0.4)

    fluxCV = fluxErrors.Clone()
    fluxCV.SetLineColor(ROOT.kRed)
    fluxCV.SetLineWidth(2)
    fluxCV.SetFillColor(0)

    fluxErrors.SetMinimum(.5)
    fluxErrors.SetMaximum(1.5)
    fluxErrors.SetLineWidth(0)
    fluxErrors.GetYaxis().SetTitle("#splitline{Ratio to}{Null Flux CV}")

    margRatio.SetLineStyle(2)
    oscRatio.SetFillColor(0)

    RatioAxis(fluxErrors,MNVPLOTTER)
    fluxErrors.Draw("E2")
    oscRatio.Draw("same hist l")
    margRatio.Draw("same l")
    fluxCV.Draw("same hist")

    x1, y1, x2, y2 = XtoPad(0.05),YtoPad(0.65),XtoPad(.4),YtoPad(.95)

    leg0 = ROOT.TLegend(x1,y1,x2,y2)
    leg0.AddEntry(fluxErrors,"Null Hypothesis Flux CV","f")
    leg0.AddEntry(oscRatio,"Osc. Model","l")
    leg0.AddEntry(margRatio,"Null Marg'd","l")
    leg0.Draw()

    # ------------------------------------ Second Pad -----------------------------------------
    pad1 = c1.FindObject("pad_{}_{}".format(0,1))
    pad1.cd()

    nullErrors = h_null.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    for whichBin in range(0, nullErrors.GetXaxis().GetNbins()+1): 
        nullErrors.SetBinError(whichBin, max(nullErrors.GetBinContent(whichBin), 1e-9))
        nullErrors.SetBinContent(whichBin, 1)
    nullErrors.SetFillColorAlpha(ROOT.kPink + 1, 0.4)

    straightLine = nullErrors.Clone()
    straightLine.SetLineColor(ROOT.kRed)
    straightLine.SetLineWidth(2)
    straightLine.SetFillColor(0)

    nullErrors.SetMinimum(0)
    nullErrors.SetMaximum(2)
    nullErrors.GetYaxis().SetTitle("#splitline{Ratio to Null}{Hypothesis}")

    margRatio.SetLineStyle(2)
    oscRatio.SetFillColor(0)

    RatioAxis(nullErrors,MNVPLOTTER)
    nullErrors.Draw("E2")
    nullRatio.Draw("same")
    oscRatio.Draw("same hist l")
    margRatio.Draw("same l")
    straightLine.Draw("same hist")

    x1, y1, x2, y2 = XtoPad(0.05),YtoPad(0.65),XtoPad(.4),YtoPad(.95)

    leg1 = ROOT.TLegend(x1,y1,x2,y2)
    leg1.AddEntry(h_data,"Data","p")
    leg1.AddEntry(oscRatio,"Osc. Model #chi^{2}="+"{:.2f}".format(chi2_model),"l")
    leg1.AddEntry(straightLine,"Null #chi^{2}="+"{:.2f}".format(chi2_null),"l")
    leg1.AddEntry(margRatio,"Null Marg'd #chi^{2}="+"{:.2f}+{:.2f}".format(chi2_null_marg-null_pen,null_pen),"l")
    leg1.Draw()

    # ------------------------------------ Third Pad -----------------------------------------
    pad2 = c1.FindObject("pad_{}_{}".format(0,0))
    pad2.cd()

    nullMargRatio =  h_data.Clone()
    oscMargRatio =  h_osc_marg.Clone()

    nullMargRatio.Divide(nullMargRatio,h_null_marg)
    oscMargRatio.Divide(oscMargRatio, h_null_marg)

    h_null_marg.PopVertErrorBand("Flux")
    h_null_marg.AddMissingErrorBandsAndFillWithCV(h_data)
    nullMargErrors = h_null_marg.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    for whichBin in range(0, nullMargErrors.GetXaxis().GetNbins()+1): 
        nullMargErrors.SetBinError(whichBin, max(nullMargErrors.GetBinContent(whichBin), 1e-9))
        nullMargErrors.SetBinContent(whichBin, 1)
    nullMargErrors.SetFillColorAlpha(ROOT.kPink + 1, 0.4)

    RatioAxis(nullMargErrors,MNVPLOTTER)
    nullMargErrors.SetMinimum(0)
    nullMargErrors.SetMaximum(2)
    nullMargErrors.GetYaxis().SetTitle("#splitline{Ratio to Marg'd}{Null Hypothesis}")
    nullMargRatio.SetFillColor(0)
    nullMargRatio.SetLineWidth(2)
    oscMargRatio.SetFillColor(0)

    nullMargErrors.Draw("E2")
    nullMargRatio.Draw("same")
    oscMargRatio.Draw("same hist l")
    straightLine.Draw("same hist")

    x1, y1, x2, y2 = XtoPad(0.05),YtoPad(0.65),XtoPad(.4),YtoPad(.95)

    leg2 = ROOT.TLegend(x1,y1,x2,y2)
    leg2.AddEntry(h_data,"Data","p")
    leg2.AddEntry(oscMargRatio,"Osc. Model #chi^{2}="+"{:.2f}+{:.2f}".format(chi2_model_marg-model_pen,model_pen),"l")
    leg2.AddEntry(straightLine,"Null #chi^{2}="+"{:.2f}+{:.2f}".format(chi2_null_marg-null_pen,null_pen),"l")
    leg2.Draw()

    c1.Print("plots/{}_stitched_ratios_{:.1f}_{:.3f}_{:.4f}.png".format(name,parameters["m"],parameters['ue4'],parameters['umu4']))

def PlotOscillationGrid(sample_histogram,parameters,name="",exclude="",lam=1):
    plots = histogram.keys
    nullSolution,nullPen = FluxSolution(histogram,invCov=invCov,usePseudo=usePseudo,exclude=exclude,lam=lam)
    plots.append("fhc_ratio")
    plots.append("rhc_ratio")
    histogram.titles["fhc_ratio"] = "FHC CC #nu_{#mu}/#nu_{e} Ratio"
    histogram.titles["rhc_ratio"] = "RHC CC anti #nu_{#mu}/#nu_{e} Ratio"

    mc_hists = []
    osc_hists = []
    data_hists = []
    titles = []

    for i,plot in enumerate(plots):
        title = histogram.titles[plot]
        titles.append(title)

        if 'ratio' in plot:
            numu_hists,numu_total_hist = OscillateSubHistogram(histogram,plot[:3]+'_numu_selection',parameters["m"],parameters['ue4'],parameters['umu4'],parameters['utau4'])
            nue_hists,nue_total_hist = OscillateSubHistogram(histogram,plot[:3]+'_nue_selection',parameters["m"],parameters['ue4'],parameters['umu4'],parameters['utau4'])

            total_hist = numu_total_hist.Clone()
            total_hist.Divide(total_hist,nue_total_hist)
            swap_fraction = total_hist.Clone()
            norm_fraction = total_hist.Clone()

            hists = []
            for j,nue in enumerate(nue_hists):
                if "#mu" in nue.GetTitle():
                    swap_fraction.Divide(nue,swap_fraction)
                    norm_fraction.Add(swap_fraction,-1)

                    norm_fraction.SetTitle("#nu_{#mu}/#nu_{e}")
                    norm_fraction.SetFillStyle(3244)
                    norm_fraction.SetLineColor(ROOT.kRed)
                    norm_fraction.SetFillColor(ROOT.kRed)

                    swap_fraction.SetTitle("#nu_{#mu}/"+nue.GetTitle())
                    swap_fraction.SetFillStyle(3244)
                    swap_fraction.SetLineColor(ROOT.kBlue)
                    swap_fraction.SetFillColor(ROOT.kBlue)

                    hists.append(norm_fraction)
                    hists.append(swap_fraction)

            h_data = histogram.data_hists[plot[:3]+'_numu_selection'].Clone()
            h_data.Divide(h_data,histogram.data_hists[plot[:3]+'_nue_selection'])
            subSample = histogram.mc_hists[plot[:3]+'_numu_selection'].Clone()
            subSample.Divide(subSample,histogram.mc_hists[plot[:3]+'_nue_selection'])
        else:
            hists,total_hist = OscillateSubHistogram(histogram,plot,parameters["m"],parameters['ue4'],parameters['umu4'],parameters['utau4'])
            h_data = histogram.data_hists[plot].Clone()
            subSample = histogram.mc_hists[plot].Clone()

        cv = np.array(subSample)[1:-1]
        mc = np.array(total_hist)[1:-1]

        weights = ReweightCV(total_hist,fluxSolution=nullSolution,cv=cv,mc=mc)
        total_hist.PopVertErrorBand("Flux")
        total_hist.AddMissingErrorBandsAndFillWithCV(h_data)

        TArray = ROOT.TObjArray()
        for hist in hists:
            hist.DivideSingle(hist,weights)
            if 'elastic' in plot:
                hist.Scale(2,'width')
            elif 'ratio' not in plot:
                hist.Scale(1,'width')
            TArray.Add(hist)

        if "elastic" in plot:
            h_data.Scale(2,'width')
            total_hist.Scale(2,'width')
        elif 'ratio' not in plot:
            h_data.Scale(1,'width')
            total_hist.Scale(1,'width')

        mc_hists.append(TArray)
        osc_hists.append(total_hist)
        data_hists.append(h_data)

    overall = plot_osc_side_by_side(mc_hists, osc_hists, data_hists, titles,MNVPLOTTER)
    overall.Print("plots/grid_oscillated_{:.1f}_{:.3f}_{:.4f}.png".format(parameters["m"],parameters['ue4'],parameters['umu4']))

def PlotOscillationEffects(sample_histogram,parameters,name="",plotSamples=False,usePseudo=False,exclude="",lam=1):
    histogram = copy.deepcopy(sample_histogram)
    invCov=histogram.GetInverseCovarianceMatrix(sansFlux=True)

    h_null = histogram.GetMCHistogram()
    h_osc = histogram.GetOscillatedHistogram()
    h_data = histogram.GetPseudoHistogram() if usePseudo else histogram.GetDataHistogram()

    chi2_null,_ = Chi2DataMC(histogram,invCov=invCov,marginalize=False,usePseudo=usePseudo,exclude=exclude,lam=lam)
    chi2_model,_ = Chi2DataMC(histogram,invCov=invCov,marginalize=False,useOsc=True,usePseudo=usePseudo,exclude=exclude,lam=lam)

    chi2_model_marg,model_pen = Chi2DataMC(histogram,invCov=invCov,marginalize=True,useOsc=True,usePseudo=usePseudo,setHists=True,exclude=exclude)
    chi2_null_marg,null_pen = Chi2DataMC(histogram,invCov=invCov,marginalize=True,usePseudo=usePseudo,setHists=True,lam=lam,exclude=exclude)

    h_null_marg = histogram.GetMCHistogram()
    h_osc_marg = histogram.GetOscillatedHistogram()

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
    
    h_osc_marg.SetLineStyle(2)
    h_null_marg.SetLineStyle(2)

    h_null.SetTitle(name)
    
    h_null.Draw("hist")
    h_null_marg.Draw("hist same")
    h_osc.Draw("hist same")
    h_osc_marg.Draw("hist_ same")
    h_data.Draw("same")

    null = h_null.GetCVHistoWithError()
    null.SetLineColor(ROOT.kRed)
    null.SetLineWidth(2)
    null.SetMarkerStyle(0)
    null.SetFillColorAlpha(ROOT.kPink + 1, 0.3)
    null.Draw("E2 SAME")
    null_marg = h_null_marg.GetCVHistoWithError()
    null_marg.SetLineColor(ROOT.kPink)
    null_marg.SetLineWidth(2)
    null_marg.SetLineStyle(2)
    null_marg.SetMarkerStyle(0)
    null_marg.SetFillColorAlpha(ROOT.kPink + 1, 0.3)
    null_marg.Draw("E2 SAME")

    osc = h_osc.GetCVHistoWithError()
    osc.SetLineColor(ROOT.kBlue)
    osc.SetLineWidth(2)
    osc.SetMarkerStyle(0)
    osc.SetFillColorAlpha(ROOT.kBlue + 1, 0.3)
    osc.Draw("E2 SAME")
    osc_marg = h_osc_marg.GetCVHistoWithError()
    osc_marg.SetLineColor(ROOT.kBlue)
    osc_marg.SetLineWidth(2)
    osc_marg.SetLineStyle(2)
    osc_marg.SetMarkerStyle(0)
    osc_marg.SetFillColorAlpha(ROOT.kBlue + 1, 0.3)
    osc_marg.Draw("E2 SAME")

    top_text = "#Delta m^{2}="+"{:.1f}".format(parameters["m"])+" |U_{e4}|^{2}="+"{:.4f}".format(parameters["ue4"])
    bot_text = "|U_{#mu4}|^{2}="+"{:.5f}".format(parameters["umu4"])+" |U_{#tau4}|^{2}="+"{:.1f}".format(parameters["utau4"])
    header = "#splitline{%s}{%s}" % (top_text,bot_text)
    MNVPLOTTER.AddPlotLabel(header,.35,.25)

    leg = ROOT.TLegend(.28,.4)

    #leg.SetBorderSize(0)
    leg.AddEntry(h_data,"Data","p")
    top_text = "Null Hypothesis"
    bot_text = "#chi^{2}="+"{:.2f}".format(chi2_null)
    leg_text = "#splitline{%s}{%s}" % (top_text,bot_text)
    leg.AddEntry(h_null,leg_text,"l")
    top_text = "Null Hypothesis Marg."
    bot_text = "#chi^{2}="+"{:.2f} + {:.2f} penalty".format(chi2_null_marg-null_pen,null_pen)
    leg_text = "#splitline{%s}{%s}" % (top_text,bot_text)
    leg.AddEntry(h_null_marg,leg_text,"l")
    top_text = "Best Osc. Fit"
    bot_text = "#chi^{2}="+"{:.2f}".format(chi2_model)
    leg_text = "#splitline{%s}{%s}" % (top_text,bot_text)
    leg.AddEntry(h_osc,leg_text,"l")
    top_text = "Best Osc. Fit Marg."
    bot_text = "#chi^{2}="+"{:.2f} + {:.2f} penalty".format(chi2_model_marg-model_pen,model_pen)
    leg_text = "#splitline{%s}{%s}" % (top_text,bot_text)
    leg.AddEntry(h_osc_marg,leg_text,"l")
    leg.Draw()

    nullRatio =  h_data.Clone()
    nullMargRatio =  h_null_marg.Clone()
    oscRatio =  h_osc.Clone()
    oscMargRatio =  h_osc_marg.Clone()

    nullRatio.Divide(nullRatio,h_null)
    nullMargRatio.Divide(nullMargRatio,h_null)
    oscRatio.Divide(oscRatio, h_null)
    oscMargRatio.Divide(oscMargRatio, h_null)

    bottom.cd()
    bottom.SetTopMargin(0)
    bottom.SetBottomMargin(0.3)

    nullErrors = h_null.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    for whichBin in range(0, nullErrors.GetXaxis().GetNbins()+1): 
        nullErrors.SetBinError(whichBin, max(nullErrors.GetBinContent(whichBin), 1e-9))
        nullErrors.SetBinContent(whichBin, 1)

    nullRatio.SetLineColor(ROOT.kBlack)
    nullRatio.SetLineWidth(3)

    #Error envelope for the MC
    nullErrors.SetLineWidth(0)
    nullErrors.SetMarkerStyle(0)
    nullErrors.SetFillColorAlpha(ROOT.kBlue + 1, 0.4)
    nullErrors.GetYaxis().SetTitle("#splitline{Ratio to Null}{Hypothesis}")
    RatioAxis(nullErrors,MNVPLOTTER)
    nullErrors.GetXaxis().SetTitle("Bin Number")
    nullErrors.SetMinimum(0)
    nullErrors.SetMaximum(2)
    nullErrors.Draw("E2")

    #Draw the data ratios
    nullRatio.Draw("same")
    nullMargRatio.Draw("same hist l")
    oscRatio.Draw('same hist l')
    oscMargRatio.Draw('same hist l')

    #Draw a flat line at 1 for oscRatio of MC to itself
    straightLine = nullErrors.Clone()
    straightLine.SetLineColor(ROOT.kRed)
    straightLine.SetLineWidth(2)
    straightLine.SetFillColor(0)
    straightLine.Draw("HIST SAME")

    top.cd()

    overall.Print("plots/{}_stitched_{:.1f}_{:.3f}_{:.4f}.png".format(name,parameters["m"],parameters['ue4'],parameters['umu4']))

    if not plotSamples:
        return

    plots = histogram.keys
    nullSolution,nullPen = FluxSolution(histogram,invCov=invCov,usePseudo=usePseudo,exclude=exclude,lam=lam)
    plots.append("fhc_ratio")
    plots.append("rhc_ratio")
    histogram.titles["fhc_ratio"] = "FHC CC #nu_{#mu}/#nu_{e} Ratio"
    histogram.titles["rhc_ratio"] = "RHC CC anti #nu_{#mu}/#nu_{e} Ratio"

    for i,plot in enumerate(plots):
        title = histogram.titles[plot]
        margin = .12
        bottomFraction = .2
        overall = ROOT.TCanvas(plot)
        top = ROOT.TPad("DATAMC", "DATAMC", 0, bottomFraction, 1, 1)
        bottom = ROOT.TPad("Ratio", "Ratio", 0, 0, 1, bottomFraction+margin)

        top.Draw()
        bottom.Draw()

        top.cd()

        if 'ratio' in plot:
            numu_hists,numu_total_hist = OscillateSubHistogram(histogram,plot[:3]+'_numu_selection',parameters["m"],parameters['ue4'],parameters['umu4'],parameters['utau4'])
            nue_hists,nue_total_hist = OscillateSubHistogram(histogram,plot[:3]+'_nue_selection',parameters["m"],parameters['ue4'],parameters['umu4'],parameters['utau4'])

            total_hist = numu_total_hist.Clone()
            total_hist.Divide(total_hist,nue_total_hist)
            swap_fraction = total_hist.Clone()
            norm_fraction = total_hist.Clone()

            hists = []
            for j,nue in enumerate(nue_hists):
                if "#mu" in nue.GetTitle():
                    swap_fraction.Divide(nue,swap_fraction)
                    norm_fraction.Add(swap_fraction,-1)

                    norm_fraction.SetTitle("#nu_{#mu}/#nu_{e}")
                    norm_fraction.SetFillStyle(3244)
                    norm_fraction.SetLineColor(ROOT.kRed)
                    norm_fraction.SetFillColor(ROOT.kRed)

                    swap_fraction.SetTitle("#nu_{#mu}/"+nue.GetTitle())
                    swap_fraction.SetFillStyle(3244)
                    swap_fraction.SetLineColor(ROOT.kBlue)
                    swap_fraction.SetFillColor(ROOT.kBlue)

                    hists.append(norm_fraction)
                    hists.append(swap_fraction)

            h_data = histogram.data_hists[plot[:3]+'_numu_selection'].Clone()
            h_data.Divide(h_data,histogram.data_hists[plot[:3]+'_nue_selection'])
            subSample = histogram.mc_hists[plot[:3]+'_numu_selection'].Clone()
            subSample.Divide(subSample,histogram.mc_hists[plot[:3]+'_nue_selection'])
        else:
            hists,total_hist = OscillateSubHistogram(histogram,plot,parameters["m"],parameters['ue4'],parameters['umu4'],parameters['utau4'])
            h_data = histogram.data_hists[plot].Clone()
            subSample = histogram.mc_hists[plot].Clone()

        cv = np.array(subSample)[1:-1]
        mc = np.array(total_hist)[1:-1]

        weights = ReweightCV(total_hist,fluxSolution=nullSolution,cv=cv,mc=mc)
        total_hist.PopVertErrorBand("Flux")
        total_hist.AddMissingErrorBandsAndFillWithCV(h_data)

        TArray = ROOT.TObjArray()
        for hist in hists:
            hist.DivideSingle(hist,weights)
            if 'elastic' in plot:
                hist.Scale(2,'width')
            elif 'ratio' not in plot:
                hist.Scale(1,'width')
            TArray.Add(hist)

        if "elastic" in plot:
            h_data.Scale(2,'width')
            total_hist.Scale(2,'width')
        elif 'ratio' not in plot:
            h_data.Scale(1,'width')
            total_hist.Scale(1,'width')

        MNVPLOTTER.DrawDataStackedMC(h_data,TArray,1.,"TR","Data",-1,-1,-1) 
        for obj in top.GetListOfPrimitives():
            if type(obj) == ROOT.TH1D:
                obj.SetTitle(title)

        bottom.cd()
        bottom.SetTopMargin(0)
        bottom.SetBottomMargin(0.3)

        ratio = h_data.Clone()
        ratio.Divide(ratio, total_hist)

        #Now fill mcRatio with 1 for bin content and fractional error
        mcRatio = total_hist.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
        for whichBin in range(1, mcRatio.GetXaxis().GetNbins()+1): 
            mcRatio.SetBinError(whichBin, max(mcRatio.GetBinContent(whichBin), 1e-9))
            mcRatio.SetBinContent(whichBin, 1)

        #Error envelope for the MC
        mcRatio.SetLineColor(ROOT.kRed)
        mcRatio.SetLineWidth(2)
        mcRatio.SetMarkerStyle(0)
        mcRatio.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
        mcRatio.SetFillStyle(1001)
        mcRatio.GetYaxis().SetTitle("Data/Model")
        mcRatio.GetXaxis().SetTitle(ratio.GetXaxis().GetTitle())
        mcRatio.SetMinimum(0)
        mcRatio.SetMaximum(2)
        RatioAxis(mcRatio,MNVPLOTTER)
        mcRatio.Draw("E2")

        ratio.Draw('E1 X0 SAME')

        straightLine = mcRatio.Clone()
        straightLine.SetFillStyle(0)
        straightLine.Draw("HIST SAME")
        ROOT.gStyle.SetOptTitle(1)

        overall.Print("plots/{}_oscillated_{:.1f}_{:.3f}_{:.4f}.png".format(plot,parameters["m"],parameters['ue4'],parameters['umu4']))

def DataMCPlot(mnv_data,mnv_mc,cv_data,cv_mc):
    #mnv_data.Scale(1,"width")
    #mnv_mc.Scale(1,"width")
    #cv_data.Scale(1,"width")
    #cv_mc.Scale(1,"width")

    c1 = ROOT.TCanvas()
    margin = .12
    bottomFraction = .2
    overall = ROOT.TCanvas("Data/MC")
    top = ROOT.TPad("DATAMC", "DATAMC", 0, bottomFraction, 1, 1)
    bottom = ROOT.TPad("Ratio", "Ratio", 0, 0, 1, bottomFraction+margin)

    top.Draw()
    bottom.Draw()

    top.cd()
    #top.SetTopMargin(0.1)

    if mnv_data.GetBinContent(mnv_data.GetMaximumBin()) > 1e4:
        top.SetLogy()

    MNVPLOTTER.axis_minimum = 0.1
    cv_data.GetXaxis().SetTitle("Bin number")
    cv_data.GetYaxis().SetTitle("Entries")

    ratio = cv_data.Clone()
    data = cv_data.GetCVHistoWithError()
    data.SetMarkerColor(ROOT.kBlack)
    data.SetTitle(mnv_mc.GetName())
    #data.SetTitleSize(0.15,'t')
    #data.SetTitleOffset(0.1)
    #data.SetLabelSize(0.15)
    mnv_mc.SetLineColor(ROOT.kBlue)
    mnv_mc.SetLineWidth(3)
    mnv_mc.SetMarkerStyle(0)

    data.Draw("P")
    mnv_mc.Draw("HIST SAME")
    cv_mc.SetLineColor(ROOT.kRed)
    cv_mc.SetLineWidth(3)
    cv_mc.SetMarkerStyle(0)
    #cv_data.Draw("P SAME")
    cv_mc.Draw("HIST SAME")

    mc = mnv_mc.GetCVHistoWithError()
    mc.SetLineColor(ROOT.kBlue)
    mc.SetLineWidth(3)
    mc.SetMarkerStyle(0)
    mc.SetFillColorAlpha(ROOT.kBlue + 1, 0.4)
    mc.Draw("E2 SAME")
    cv = cv_mc.GetCVHistoWithError()
    cv.SetLineColor(ROOT.kRed)
    cv.SetLineWidth(3)
    cv.SetMarkerStyle(0)
    cv.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
    cv.Draw("E2 SAME")

    leg = ROOT.TLegend(.7,.7,1,1)
    chi2_model = Chi2DataMC(mnv_data,mnv_mc)
    leg.SetHeader("#chi^{2}_{model}="+"{:.2f}".format(chi2_model),"C")
    leg.AddEntry(cv_data,"CV data",'p')
    leg.AddEntry(cv_mc,"CV MC",'l')
    leg.AddEntry(mnv_mc,"Universe MC","l")
    leg.Draw()

    ratio.Divide(ratio, cv_mc)

    bottom.cd()
    bottom.SetTopMargin(0)
    bottom.SetBottomMargin(0.3)

    #Now fill mcRatio with 1 for bin content and fractional error
    mcRatio = mnv_mc.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    for whichBin in range(1, mcRatio.GetXaxis().GetNbins()+1): 
        mcRatio.SetBinError(whichBin, max(mcRatio.GetBinContent(whichBin), 1e-9))
        mcRatio.SetBinContent(whichBin, 1)

    ratio.SetTitle("")
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetLineWidth(3)
    ratio.SetTitleSize(0)

    ratio.GetYaxis().SetTitle("#splitline{Ratio to Null}{Hypothesis}")
    ratio.GetYaxis().SetLabelSize(.13)
    ratio.GetYaxis().SetTitleSize(0.1)
    ratio.GetYaxis().SetTitleOffset(0.6)
    ratio.GetYaxis().SetNdivisions(505) #5 minor divisions between 5 major divisions.  I'm trying to match a specific paper here.
    ratio.GetXaxis().SetTitleSize(0.16)
    ratio.GetXaxis().SetTitleOffset(0.9)
    ratio.GetXaxis().SetLabelSize(.15)
    ratio.GetXaxis().SetTitle("Bin Number")
    ratio.SetMinimum(0)
    ratio.SetMaximum(2)
    ratio.Draw()

    univratio = mnv_mc.Clone()
    univratio.Divide(univratio,cv_mc)
    univratio.SetLineColor(ROOT.kBlue)
    univratio.SetLineWidth(3)
    univratio.Draw("HIST SAME")

    #Error envelope for the MC
    mcRatio.SetLineColor(ROOT.kRed)
    mcRatio.SetLineWidth(3)
    mcRatio.SetMarkerStyle(0)
    mcRatio.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
    mcRatio.Draw("E2 SAME")

    #Draw a flat line at 1 for ratio of MC to itself
    straightLine = mcRatio.Clone()
    straightLine.SetFillStyle(0)
    straightLine.Draw("HIST SAME")
    ROOT.gStyle.SetOptTitle(1)
    c1.Update()
    #MNVPLOTTER.WritePreliminary(0.4, 0.05, 5e-2, True)
    overall.Print("plots/{}.png".format(mnv_mc.GetName()))

def DataMCCVPlot(cv_data,cv_mc,title="mc_stitched.png"):
    #mnv_data.Scale(1,"width")
    #mnv_mc.Scale(1,"width")
    #cv_data.Scale(1,"width")
    #cv_mc.Scale(1,"width")

    c1 = ROOT.TCanvas()
    margin = .12
    bottomFraction = .2
    overall = ROOT.TCanvas("Data/MC")
    top = ROOT.TPad("DATAMC", "DATAMC", 0, bottomFraction, 1, 1)
    bottom = ROOT.TPad("Ratio", "Ratio", 0, 0, 1, bottomFraction+margin)

    top.Draw()
    bottom.Draw()

    top.cd()
    if title != "mc_stitched.png":
        cv_data.GetXaxis().SetTitle("Energy Estimator (GeV)")
        cv_data.GetYaxis().SetTitle("Muon / Electron Flavor Ratio")
    else:
        cv_data.GetXaxis().SetTitle("Bin Number")
        cv_data.GetYaxis().SetTitle("Entries")
    cv_data.GetXaxis().SetTitle("Bin Number")
    cv_data.GetYaxis().SetTitle("Entries")
    top.SetLogy()
    MNVPLOTTER.axis_minimum = 0.1

    ratio = cv_data.GetCVHistoWithError().Clone()

    cv_data.SetMarkerColor(ROOT.kBlack)
    cv_data.SetLineColor(ROOT.kBlack)
    cv_data.SetMarkerSize(1)
    cv_data.SetMarkerStyle(20)
    cv_data.SetLineWidth(3)
    cv_data.SetMinimum(0.1)

    cv_data.GetCVHistoWithError().Draw("P")
    cv_mc.SetLineColor(ROOT.kRed)
    cv_mc.SetLineWidth(3)
    cv_mc.SetMarkerStyle(0)
    cv_data.Draw("P SAME")
    cv_mc.Draw("HIST SAME")
    cv = cv_mc.GetCVHistoWithError()
    cv.SetLineColor(ROOT.kRed)
    cv.SetLineWidth(3)
    cv.SetMarkerStyle(0)
    cv.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
    cv.Draw("E2 SAME")

    leg = ROOT.TLegend(.65,.35,.8,.55)
    leg.AddEntry(cv_data,"Data",'p')
    leg.AddEntry(cv_mc,"Null Hypothesis",'l')
    leg.Draw()

    ratio.Divide(ratio, cv_mc)

    bottom.cd()
    bottom.SetTopMargin(0)
    bottom.SetBottomMargin(0.3)

    #Now fill mcRatio with 1 for bin content and fractional error
    mcRatio = cv_mc.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    for whichBin in range(1, mcRatio.GetXaxis().GetNbins()+1): 
        mcRatio.SetBinError(whichBin, max(mcRatio.GetBinContent(whichBin), 1e-9))
        mcRatio.SetBinContent(whichBin, 1)

    ratio.SetTitle("")
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetLineWidth(3)
    ratio.SetTitleSize(0)

    ratio.GetYaxis().SetTitle("#splitline{Ratio to}{Null Hypothesis}")
    ratio.GetYaxis().SetLabelSize(.13)
    ratio.GetYaxis().SetTitleSize(0.1)
    ratio.GetYaxis().SetTitleOffset(0.6)
    ratio.GetYaxis().SetNdivisions(505) #5 minor divisions between 5 major divisions.  I'm trying to match a specific paper here.
    ratio.GetXaxis().SetTitleSize(0.16)
    ratio.GetXaxis().SetTitleOffset(0.9)
    ratio.GetXaxis().SetLabelSize(.15)
    if title != "mc_stitched.png":
        ratio.GetXaxis().SetTitle("Energy Estimator (GeV)")
    else:
        ratio.GetXaxis().SetTitle("Bin Number")
    ratio.SetMinimum(0)
    ratio.SetMaximum(2)
    ratio.Draw()

    #Error envelope for the MC
    mcRatio.SetLineColor(ROOT.kRed)
    mcRatio.SetLineWidth(3)
    mcRatio.SetMarkerStyle(0)
    mcRatio.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
    mcRatio.Draw("E2 SAME")

    #Draw a flat line at 1 for ratio of MC to itself
    straightLine = mcRatio.Clone()
    straightLine.SetFillStyle(0)
    straightLine.Draw("HIST SAME")
    top.cd()
    
    MNVPLOTTER.WritePreliminary(0.4, 0.05, 5e-2, True)
    overall.Print("plots/{}".format(title))

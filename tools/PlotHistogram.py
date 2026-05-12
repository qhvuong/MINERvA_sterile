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
try:
    import psutil
except ImportError:
    psutil = None
import multiprocessing
import threading
nthreads = 4
from array import array
from tools.Fitters import *
from tools.StitchedHistogram import *
from tools.Helper import *

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
legend_text_size = .035

class PlottingContainer:
    def __init__(self,tag,histogram):
        self.tag = tag
        self.histogram = histogram
        self.exclude = ""
        self.lam = 1
        self.flux_solution = None
        self.osc_parameters = {"m":0,"ue4":0,"umu4":0,"utau4":0}
        self.invCov = histogram.GetInverseCovarianceMatrix(sansFlux=True)

        self.exclude_samples = ["", "", "ratio"]
        self.titles = ["#nu+e, IMD, CC #nu_{#mu}, CC #nu_{#mu}/#nu_{e} #lambda = 1","#nu+e, IMD, CC #nu_{#mu}, CC #nu_{#mu}/#nu_{e} #lambda = 12", "#nu+e, IMD, CC #nu_{#mu} #lambda = 1"]
        self.colors = [ROOT.kRed,ROOT.kBlue,ROOT.kGreen]
        self.lams = [1,12,1]

    def SetExclude(self,exclude):
        self.exclude = exclude

    def SetHistogram(self,histogram):
        self.histogram = histogram

    def SetFluxSolution(self,solution):
        self.flux_solution = solution

    def SetOscParameters(self,oscDict):
        self.osc_parameters = oscDict
        self.histogram.OscillateHistogram(self.histogram, oscDict['m'], oscDict['ue4'], oscDict['umu4'], oscDict['utau4'])

    def SetLambda(self,lam):
        self.lam = lam

    def SetInverseCovariance(self,inv):
        self.invCov = inv

    def PlotScatteringIntegrals(self):
        subSample = self.histogram.mc_hists["fhc_elastic"].Clone()
        band = subSample.GetVertErrorBand("Flux")
        nhists = band.GetNHists()
        universes = np.array([np.array(band.GetHist(l))[1:-1] for l in range(nhists)])
        fhc_integrals = universes.sum(axis=1)
        old_fhc = subSample.Integral()

        subSample = self.histogram.mc_hists["rhc_elastic"].Clone()
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

        histogram = copy.deepcopy(self.histogram)
        statistic = Statistics(histogram,exclude=self.exclude,lam=self.lam)
        statistic.Chi2DataMC(marginalize=True)

        #for i,exclude in enumerate(self.exclude_samples):
        #    subSample = self.histogram.mc_hists["fhc_elastic"].Clone()
        #    statistic.GetFluxFitter().ReweightToFluxSolution(subSample)
        #    new_fhc = subSample.Integral()

        #    subSample = self.histogram.mc_hists["rhc_elastic"].Clone()
        #    statistic.GetFluxFitter().ReweightToFluxSolution(subSample)
        #    new_rhc = subSample.Integral()

        #    g_ = ROOT.TGraph()
        #    ROOT.SetOwnership(g_, False)
        #    g_.SetPoint(0,new_fhc,new_rhc)
        #    g_.SetTitle(self.titles[i])
        #    g_.SetMarkerStyle(29)
        #    g_.SetMarkerSize(3)
        #    g_.SetLineWidth(0)
        #    g_.SetMarkerColorAlpha(self.colors[i],0.4)
        #    mg.Add(g_)
        subSample = self.histogram.mc_hists["fhc_elastic"].Clone()
        statistic.GetFluxFitter().ReweightToFluxSolution(subSample)
        new_fhc = subSample.Integral()

        subSample = self.histogram.mc_hists["rhc_elastic"].Clone()
        statistic.GetFluxFitter().ReweightToFluxSolution(subSample)
        new_rhc = subSample.Integral()

        g_ = ROOT.TGraph()
        ROOT.SetOwnership(g_, False)
        g_.SetPoint(0,new_fhc,new_rhc)
        g_.SetTitle("Profiled Flux Solution - New")
        g_.SetMarkerStyle(29)
        g_.SetMarkerSize(3)
        g_.SetLineWidth(0)
        g_.SetMarkerColorAlpha(ROOT.kRed,0.4)
        mg.Add(g_)

        # TODO old sample with just 100 universes
        gOld = ROOT.TGraph()
        ROOT.SetOwnership(gOld, False)
        gOld.SetPoint(0,1187,780)
        gOld.SetTitle("Profiled Flux Solution - Old")
        gOld.SetMarkerStyle(29)
        gOld.SetMarkerSize(3)
        gOld.SetLineWidth(0)
        gOld.SetMarkerColorAlpha(ROOT.kBlue,0.4)
        mg.Add(gOld)

        mg.Draw("AP")
        pad = ROOT.gPad
        pad.BuildLegend()
        c0.Print("plots/integrated_elastic_events.png")

    def PlotFluxReweight(self,beam):
        beam = beam.upper()
        if beam == "FHC":
            fileExt = "minervame1D1M1NWeightedAve.root"
            pdgPref = ""
            oldNueRatio = [1.16,1.02,.914,.893,.925,.910,.855,.855,.878,.961,.995,.989,1.01,1.08,1.11,1.12,1.08,1.09,1.13,1.18]
            oldNumuRatio = [0.9069767961765031, 0.995348802548998, 0.8953489622339715, 0.8813954763376145, 0.9116278871708561, 0.90232570518215, 0.8488372006372495, 0.851162852591075, 0.893023310280146, 0.972093134663935, 1.0162790313935335, 1.009302288445355, 1.0302325172898905, 1.0906977647829694, 1.130232570518215, 1.1325582224720405, 1.097674507731148, 1.065116232030783, 1.139534965420219, 1.2]
        elif beam == "RHC":
            fileExt = "minervame6A.root"
            pdgPref = "-"
            oldNumuRatio = [1.10,1.02,.999,.918,.947,.954,.907,.941,.981,1.02,1.03,.992,1.01,1.03,1.02,1.03,.954,1.02,1.01,.951]
            oldNueRatio = [1.18,.893,1.16,1.05,.993,.984,1.03,1.02,1.06,1.08,1.11,1.08,1.06,1.29,1.31,1.31,1.3,1.3,1.3,1.3]
            oldNumuRatio = [0.8302326237465395, 0.981395316652641, 0.972093134663935, 0.9186046301190345, 0.948837253865574, 0.9604650878081056, 0.9232559340266855, 0.9511629058193996, 1.000000106456649, 1.055814050042077, 1.0674418839846085, 1.0395349121918944, 1.046511655140073, 1.1046512506793265, 1.1279069185643895, 1.0906977647829694, 1.032558169243716, 1.069767535938434, 1.1232558275700364, 1.088372112829144]

        f_numu = ROOT.TFile.Open(plotutils+'/data/flux/flux-g4numiv6-pdg'+pdgPref+"14-"+fileExt)
        numu = f_numu.Get("flux_E_unweighted")
        f_numu.Close()
        f_nue = ROOT.TFile.Open(plotutils+'/data/flux/flux-g4numiv6-pdg'+pdgPref+"12-"+fileExt)
        nue = f_nue.Get("flux_E_unweighted")
        f_nue.Close()

        numu_fluxes = []
        nue_fluxes = []

        histogram = copy.deepcopy(self.histogram)
        statistic = Statistics(histogram,exclude=self.exclude,lam=self.lam)
        statistic.Chi2DataMC(marginalize=True)

        #for i,exclude in enumerate(self.exclude_samples):
        #    new_numu = numu.Clone()
        #    new_nue = nue.Clone()
        #    statistic.GetFluxFitter().ReweightToFluxSolution(new_numu)
        #    statistic.GetFluxFitter().ReweightToFluxSolution(new_nue)

        #    nue_fluxes.append(new_nue)
        #    numu_fluxes.append(new_numu)
        new_numu = numu.Clone()
        new_nue = nue.Clone()
        statistic.GetFluxFitter().ReweightToFluxSolution(new_numu)
        statistic.GetFluxFitter().ReweightToFluxSolution(new_nue)

        nue_fluxes.append(new_nue)
        numu_fluxes.append(new_numu)
        titles = ["Profiled Flux Solution - New"]

        new_bins = array('d',list(range(0,21)))
        UndoBinWidthNorm(nue)
        UndoBinWidthNorm(numu)

        numu = numu.Rebin(20,"hnew",new_bins)
        numu.Scale(1,"width")
        nue = nue.Rebin(20,"hnew",new_bins)
        nue.Scale(1,"width")

        colors = [ROOT.kBlue,ROOT.kRed]

        for i in range(len(numu_fluxes)):
            UndoBinWidthNorm(numu_fluxes[i])
            numu_fluxes[i] = numu_fluxes[i].Rebin(20,str(i),new_bins)
            numu_fluxes[i].Scale(1,'width')

            UndoBinWidthNorm(nue_fluxes[i])
            nue_fluxes[i] = nue_fluxes[i].Rebin(20,str(i),new_bins)
            nue_fluxes[i].Scale(1,'width')

        oldNumuFlux = numu_fluxes[0].Clone()
        for i in range(oldNumuFlux.GetNbinsX()+1):
            if i < len(oldNumuRatio):
                oldNumuFlux.SetBinContent(i+1,numu.GetBinContent(i+1)*oldNumuRatio[i])
        oldNueFlux = nue_fluxes[0].Clone()
        for i in range(oldNueFlux.GetNbinsX()+1):
            if i < len(oldNueRatio):
                oldNueFlux.SetBinContent(i+1,nue.GetBinContent(i+1)*oldNueRatio[i])

        numu_fluxes.append(oldNumuFlux)
        nue_fluxes.append(oldNueFlux)
        titles.append("Profiled Flux Solution - Old")
        
        for i in range(len(numu_fluxes)):
            print(titles[i], np.array(numu_fluxes[i]) / np.array(numu))
            #print(np.array(nue_fluxes[i]) / np.array(nue))

        numu.GetXaxis().SetRangeUser(0,20)
        numu.GetXaxis().SetTitle("Neutrino Energy")
        numu.SetTitle(beam+" #nu_{#mu} Flux Prediction")

        nue.GetXaxis().SetRangeUser(0,20)
        nue.GetXaxis().SetTitle("Neutrino Energy")
        nue.SetTitle(beam+" #nu_{e} Flux Prediction")

        PlotWithRatio(MNVPLOTTER,"plots/"+beam+"_NuMuFlux_Reweight.png",numu,hists=numu_fluxes,titles=titles,colors=colors)
        PlotWithRatio(MNVPLOTTER,"plots/"+beam+"_NuEFlux_Reweight.png",nue,hists=nue_fluxes,titles=titles,colors=colors)

    def PlotProfileEffects(self):
        c1 = ROOT.TCanvas("C2", "canvas2", 1024, 640)
        c1.SetTopMargin(0.35)
        c1.SetRightMargin(0.05)

        histogram = copy.deepcopy(self.histogram)

        h_null =  histogram.GetMCHistogram()
        h_data = histogram.GetDataHistogram()
        invCov = self.invCov
        chi2,pen = Chi2DataMC(histogram,invCov=histogram.GetInverseCovarianceMatrix(sansFlux=False))

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
        nullErrors.SetMinimum(.7)
        nullErrors.SetMaximum(1.3)
        nullErrors.GetYaxis().SetTitle("#splitline{Ratio to Null}{Hypothesis}")
        nullErrors.GetXaxis().SetTitleOffset(1.5)
        nullErrors.Draw("E2")
        nullRatio.Draw("SAME")

        leg = ROOT.TLegend(0.15, 0.675, 0.95, .975)
        leg.SetTextSize(legend_text_size);
        #leg.SetNColumns(len(self.titles)//2) # // median N number of rows per column
        leg.SetNColumns(2) # // median N number of rows per column

        top = "Data w/o Profiling"
        bottom = "#chi^{2}="+"{:.2f}".format(chi2)
        legend = "#splitline{%s}{%s}" % (top,bottom)

        leg.AddEntry(nullRatio,legend,"p")

        hists = []
        for i,exclude in enumerate(self.exclude_samples):
            fluxSolution,nullPen = FluxSolution(histogram,invCov=self.invCov,exclude=exclude,lam=self.lams[i])
            chi2,penalty = Chi2DataMC(histogram,fluxSolution=fluxSolution,invCov=self.invCov,exclude=exclude,lam=self.lams[i],marginalize=True)

            hist = histogram.GetMCHistogram()
            weights = ReweightCV(hist,fluxSolution=fluxSolution)

            hist.SetName(exclude+"_{}".format(i))
            hist.SetTitle(self.titles[i])

            hist.Divide(hist,h_null)
            hist.SetFillStyle(0)
            if i == 0:
                hist.SetLineStyle(2)

            top = self.titles[i]
            bottom = "#chi^{2}="+"{:.2f}+".format(chi2-penalty)+"{:.2f} pen.".format(penalty)
            legend = "#splitline{%s}{%s}" % (top,bottom)
            hist.SetLineColor(self.colors[i])
            hist.SetTitle(legend)
            hists.append(hist)

        for hist in hists:
            hist.Draw("hist same l")
            leg.AddEntry(hist,hist.GetTitle(),"l")

        straightLine.Draw("same hist") 
        leg.Draw()

        c1.Print("plots/stitched_flux_marg_effects.png")

    def PlotOscillationEffects(self, parameters, name="", useMarg=True, plotSamples=False, usePseudo=False):
        histogram = copy.deepcopy(self.histogram)
        exclude = self.exclude
        lam = self.lam

        statistic = Statistics(histogram, exclude=exclude, lam=lam)

        chi2_model, model_pen = statistic.Chi2DataMC(
            marginalize=useMarg,
            useOsc=True,
            usePseudo=usePseudo
        )
        chi2_null, null_pen = statistic.Chi2DataMC(
            marginalize=useMarg,
            usePseudo=usePseudo
        )

        if useMarg:
            null_chi2_text = "Null Hyp. #chi^{{2}} = {:.2f} + {:.2f}".format(
                chi2_null - null_pen, null_pen
            )
            osc_chi2_text = "Osc. Model #chi^{{2}} = {:.2f} + {:.2f}".format(
                chi2_model - model_pen, model_pen
            )
        else:
            null_chi2_text = "Null Hyp. #chi^{{2}} = {:.2f}".format(chi2_null)
            osc_chi2_text = "Osc. Model #chi^{{2}} = {:.2f}".format(chi2_model)

        plot_texts = [
            null_chi2_text,
            osc_chi2_text,
            "#Delta m^{2} = " + "{}".format(parameters["m"]) + " eV^{2}",
            "|U_{e4}|^{2} = " + "{}".format(parameters["ue4"]),
            "|U_{#mu4}|^{2} = " + "{:.6f}".format(parameters["umu4"]),
            "|U_{#tau4}|^{2} = " + "{}".format(parameters["utau4"]),
        ]

        h_null = histogram.GetMCHistogram()
        h_osc = histogram.GetOscillatedHistogram()
        h_data = histogram.GetPseudoHistogram() if usePseudo else histogram.GetDataHistogram()

        h_prof = None
        if useMarg:
            h_prof = h_null.Clone()
            h_prof.SetLineColor(ROOT.kGreen)

            statistic.GetFluxFitter(useOsc=False).ReweightToFluxSolution(h_prof)
            statistic.GetFluxFitter(useOsc=True).ReweightToFluxSolution(h_osc)

        c1 = ROOT.TCanvas()
        margin = .12
        bottomFraction = .2
        overall = ROOT.TCanvas("Data/MC")
        top = ROOT.TPad("DATAMC", "DATAMC", 0, bottomFraction, 1, 1)
        bottom = ROOT.TPad("Ratio", "Ratio", 0, 0, 1, bottomFraction + margin)

        top.Draw()
        bottom.Draw()

        top.cd()
        top.SetLogy()

        h_null.SetTitle(name)

        h_null.Draw("hist")
        if useMarg and h_prof is not None:
            h_prof.Draw("hist same")
        h_osc.Draw("hist same")
        h_data.Draw("same")

        if useMarg and h_null.HasVertErrorBand("Flux"):
            h_null.PopVertErrorBand("Flux")
        null = h_null.GetCVHistoWithError()
        null.SetLineColor(ROOT.kRed)
        null.SetLineWidth(2)
        null.SetMarkerStyle(0)
        null.SetFillColorAlpha(ROOT.kPink + 1, 0.3)
        null.Draw("E2 SAME")

        if useMarg and h_osc.HasVertErrorBand("Flux"):
            h_osc.PopVertErrorBand("Flux")
        osc = h_osc.GetCVHistoWithError()
        osc.SetLineColor(ROOT.kBlue)
        osc.SetLineWidth(2)
        osc.SetMarkerStyle(0)
        osc.SetFillColorAlpha(ROOT.kBlue + 1, 0.3)
        osc.Draw("E2 SAME")

        top_text = "#Delta m^{2}=" + "{:.1f}".format(parameters["m"]) + " |U_{e4}|^{2}=" + "{:.4f}".format(parameters["ue4"])
        bot_text = "|U_{#mu4}|^{2}=" + "{:.5f}".format(parameters["umu4"]) + " |U_{#tau4}|^{2}=" + "{:.1f}".format(parameters["utau4"])
        header = "#splitline{%s}{%s}" % (top_text, bot_text)
        MNVPLOTTER.AddPlotLabel(header, .35, .25)

        leg = ROOT.TLegend(.28, .35)
        leg.AddEntry(h_data, "Data", "p")

        top_text = "Null Hypothesis"
        if useMarg:
            bot_text = "#chi^{{2}}={:.2f} + {:.2f} penalty".format(
                chi2_null - null_pen, null_pen
            )
        else:
            bot_text = "#chi^{{2}}={:.2f}".format(chi2_null)
        leg_text = "#splitline{%s}{%s}" % (top_text, bot_text)
        leg.AddEntry(h_null, leg_text, "l")

        if useMarg and h_prof is not None:
            top_text = "Null Hypothesis Profiled"
            bot_text = "#chi^{{2}}={:.2f} + {:.2f} penalty".format(
                chi2_null - null_pen, null_pen
            )
            leg_text = "#splitline{%s}{%s}" % (top_text, bot_text)
            leg.AddEntry(h_prof, leg_text, "l")

        top_text = "Best Osc. Fit"
        if useMarg:
            bot_text = "#chi^{{2}}={:.2f} + {:.2f} penalty".format(
                chi2_model - model_pen, model_pen
            )
        else:
            bot_text = "#chi^{{2}}={:.2f}".format(chi2_model)

        leg_text = "#splitline{%s}{%s}" % (top_text, bot_text)
        leg.AddEntry(h_osc, leg_text, "l")

        leg.Draw()

        nullRatio = h_data.Clone()
        oscRatio = h_osc.Clone()

        nullRatio.Divide(nullRatio, h_null)
        oscRatio.Divide(oscRatio, h_null)

        profRatio = None
        if useMarg and h_prof is not None:
            profRatio = h_prof.Clone()
            profRatio.Divide(profRatio, h_null)

        bottom.cd()
        bottom.SetTopMargin(0)
        bottom.SetBottomMargin(0.3)

        nullErrors = h_null.GetTotalError(False, True, False)
        for whichBin in range(0, nullErrors.GetXaxis().GetNbins() + 1):
            nullErrors.SetBinError(whichBin, max(nullErrors.GetBinContent(whichBin), 1e-9))
            nullErrors.SetBinContent(whichBin, 1)

        nullRatio.SetLineColor(ROOT.kBlack)
        nullRatio.SetLineWidth(3)

        if profRatio is not None:
            profRatio.SetLineColor(ROOT.kGreen)
            profRatio.SetLineWidth(3)

        oscRatio.SetLineColor(ROOT.kBlue)
        oscRatio.SetLineWidth(3)

        nullErrors.SetLineWidth(0)
        nullErrors.SetMarkerStyle(0)
        nullErrors.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
        nullErrors.GetYaxis().SetTitle("#splitline{Ratio to Null}{Hypothesis}")
        RatioAxis(nullErrors, MNVPLOTTER)
        nullErrors.GetXaxis().SetTitle("Bin Number")
        nullErrors.SetMinimum(.7)
        nullErrors.SetMaximum(1.3)
        nullErrors.Draw("E2")

        nullRatio.Draw("same")
        oscRatio.Draw("same hist l")
        if profRatio is not None:
            profRatio.Draw("same hist l")

        straightLine = nullErrors.Clone()
        straightLine.SetLineColor(ROOT.kRed)
        straightLine.SetLineWidth(2)
        straightLine.SetFillColor(0)
        straightLine.Draw("HIST SAME")

        top.cd()

        overall.Print(
            "plots/{}_stitched_{:.1f}_{:.3f}_{:.4f}.png".format(
                name, parameters["m"], parameters["ue4"], parameters["umu4"]
            )
        )

        if not plotSamples:
            return

        if not useMarg:
            print("WARNING: plotSamples=True currently expects flux fitters. Skipping sample plots because useMarg=False.")
            return

        plots = histogram.keys
        histogram = copy.deepcopy(self.histogram)

        nullSolution = statistic.GetFluxFitter(useOsc=False).GetFluxSolution()
        oscSolution = statistic.GetFluxFitter(useOsc=True).GetFluxSolution()

        mc_hists = []
        null_hists = []
        osc_hists = []
        data_hists = []
        titles = []

        plots = [item for item in plots if "nue" not in item]
        plots.insert(3, plots.pop(6))

        for i, plot in enumerate(plots):
            title = histogram.titles[plot]
            titles.append(title)

            hists = []

            if 'ratio' in plot:
                numu_hists, numu_total_hist = histogram.OscillateSubHistogram(
                    plot[:3] + '_numu_selection',
                    parameters["m"],
                    parameters['ue4'],
                    parameters['umu4'],
                    parameters['utau4']
                )
                nue_hists, nue_total_hist = histogram.OscillateSubHistogram(
                    plot[:3] + '_nue_selection',
                    parameters["m"],
                    parameters['ue4'],
                    parameters['umu4'],
                    parameters['utau4']
                )

                total_hist = numu_total_hist.Clone()
                total_hist.Divide(total_hist, nue_total_hist)

                numu_rat_hist = numu_total_hist.Clone()
                nue_rat_hist = nue_total_hist.Clone()

                for m in range(numu_total_hist.GetNbinsX() + 1):
                    frac = total_hist.GetBinContent(m) / numu_rat_hist.GetBinContent(m) if numu_rat_hist.GetBinContent(m) != 0 else 0
                    numu_rat_hist.SetBinContent(m, frac * total_hist.GetBinContent(m))
                    nue_rat_hist.SetBinContent(m, (1 - frac) * total_hist.GetBinContent(m))

                    if nue_rat_hist.HasVertErrorBand("Flux"):
                        for u in range(nue_rat_hist.GetVertErrorBand("Flux").GetNHists()):
                            numu_rat_hist.GetVertErrorBand("Flux").GetHist(u).SetBinContent(
                                m, numu_rat_hist.GetBinContent(m) * frac
                            )
                            nue_rat_hist.GetVertErrorBand("Flux").GetHist(u).SetBinContent(
                                m, nue_rat_hist.GetBinContent(m) * (1 - frac)
                            )

                numu_rat_hist.SetTitle("numu")
                nue_rat_hist.SetTitle("nue")
                hists.append(numu_rat_hist)
                hists.append(nue_rat_hist)

                h_data = histogram.data_hists[plot[:3] + '_numu_selection'].Clone()
                h_data.Divide(h_data, histogram.data_hists[plot[:3] + '_nue_selection'])

                subSample = histogram.mc_hists[plot[:3] + '_numu_selection'].Clone()
                subSample.Divide(subSample, histogram.mc_hists[plot[:3] + '_nue_selection'])

            else:
                hists, total_hist = histogram.OscillateSubHistogram(
                    plot,
                    parameters["m"],
                    parameters['ue4'],
                    parameters['umu4'],
                    parameters['utau4']
                )
                h_data = histogram.data_hists[plot].Clone()
                subSample = histogram.mc_hists[plot].Clone()

            null_hist = histogram.mc_hists[plot].Clone()

            statistic.GetFluxFitter(useOsc=False).ReweightToFluxSolution(null_hist)
            if null_hist.HasVertErrorBand("Flux"):
                null_hist.PopVertErrorBand("Flux")

            statistic.GetFluxFitter(useOsc=True).ReweightToFluxSolution(total_hist)
            if total_hist.HasVertErrorBand("Flux"):
                total_hist.PopVertErrorBand("Flux")
            total_hist.AddMissingErrorBandsAndFillWithCV(h_data)

            TArray = []
            for hist in hists:
                if hist.Integral() <= 0:
                    continue

                statistic.GetFluxFitter(useOsc=True).ReweightToFluxSolution(hist)

                if 'elastic' in plot:
                    hist.Scale(2, 'width')
                elif 'ratio' not in plot:
                    hist.Scale(1, 'width')

                if "tau" in hist.GetTitle():
                    hist.SetFillColorAlpha(ROOT.kGray, 0.6)
                elif "ratio" in plot:
                    hist.SetFillColorAlpha(ROOT.kViolet, 0.6)
                elif "mu" in hist.GetTitle():
                    hist.SetFillColorAlpha(ROOT.kBlue, 0.6)
                else:
                    hist.SetFillColorAlpha(ROOT.kRed, 0.6)

                hist.SetLineWidth(0)
                TArray.append(hist)

            if "elastic" in plot:
                h_data.Scale(2, 'width')
                total_hist.Scale(2, 'width')
                null_hist.Scale(2, 'width')
            elif 'ratio' not in plot:
                h_data.Scale(1, 'width')
                total_hist.Scale(1, 'width')
                null_hist.Scale(1, 'width')

            if "elastic" in plot:
                total_hist.SetTitle("#nu+e")
            elif "imd" in plot:
                total_hist.SetTitle("#nu_{#mu}+e^{-}#rightarrow #mu^{-}+#nu_{e}")
            elif "mu" in plot:
                total_hist.SetTitle("CC #nu_{#mu}")
            elif "ratio" in plot:
                total_hist.SetTitle("CC #nu_{#mu}/#nu_{e}")
            else:
                total_hist.SetTitle("CC #nu_{e}")

            if "fhc" in plot:
                total_hist.GetYaxis().SetTitle("#nu-mode    Events/GeV")
            else:
                total_hist.GetYaxis().SetTitle("#bar{#nu}-mode    Events/GeV")

            mc_hists.append(TArray)
            osc_hists.append(total_hist)
            null_hists.append(null_hist)
            data_hists.append(h_data)

        overall = plot_osc_side_by_side(
            mc_hists,
            null_hists,
            osc_hists,
            data_hists,
            titles,
            plot_texts,
            MNVPLOTTER,
            narrow_pads=[1, 4]
        )
        overall.Print(
            "plots/{}_oscillated_{:.1f}_{:.3f}_{:.4f}.png".format(
                name, parameters["m"], parameters['ue4'], parameters['umu4']
            )
        )
    # def PlotFluxMarginalizationEffects(self,parameters,name="",plotSamples=False,usePseudo=False):
    #     histogram = copy.deepcopy(self.histogram)
    #     histogram.SetPlottingStyle()
        
    #     invCov=histogram.GetInverseCovarianceMatrix(sansFlux=True)

    #     nullSolution,nullPen = FluxSolution(histogram,invCov=invCov,usePseudo=usePseudo)
    #     oscSolution,oscPen   = FluxSolution(histogram,invCov=invCov,useOsc=True,usePseudo=usePseudo)

    #     nullSols = ROOT.TH1D("null_solutions","Null Solutions",20,-1,1)
    #     oscSols = ROOT.TH1D("osc_solutions","Best Fit Solutions",20,-1,1)
    #     for i in range(len(nullSolution)):
    #         nullSols.Fill(nullSolution[i])
    #         oscSols.Fill(oscSolution[i])

    #     nullSolution = nullSolution/np.max(nullSolution)
    #     oscSolution = oscSolution/np.max(oscSolution)

    #     chi2_null_marg,null_pen = Chi2DataMC(histogram,invCov=invCov,marginalize=True,usePseudo=usePseudo)
    #     chi2_model_marg,model_pen = Chi2DataMC(histogram,invCov=invCov,marginalize=True,useOsc=True,usePseudo=usePseudo)
    #     chi2_null,_ = Chi2DataMC(histogram,invCov=invCov,marginalize=False,usePseudo=usePseudo)
    #     chi2_model,_ = Chi2DataMC(histogram,invCov=invCov,marginalize=False,useOsc=True,usePseudo=usePseudo)
    #     h_data = histogram.GetPseudoHistogram() if usePseudo else histogram.GetDataHistogram()

    #     h_osc = histogram.GetOscillatedHistogram()
    #     _,_ = Chi2DataMC(histogram,invCov=invCov,useOsc=True,setHists=True,marginalize=True,usePseudo=usePseudo)
    #     h_osc_marg = histogram.GetOscillatedHistogram()
    #     h_null = histogram.GetMCHistogram()
    #     _,_ = Chi2DataMC(histogram,invCov=invCov,setHists=True,marginalize=True,usePseudo=usePseudo)
    #     h_null_marg = histogram.GetMCHistogram()

    #     c1 = ROOT.TCanvas("CMarg", "canvas", 1024, 640)
    #     c1.SetFillStyle(4000)

    #     lMargin = 0.12
    #     rMargin = 0.05
    #     bMargin = 0.15
    #     tMargin = 0.05

    #     c1.Divide(1,3)

    #     # ------------------------------------ First Pad -----------------------------------------
    #     c1.cd(1)
    #     pad = ROOT.gPad
    #     pad.SetMargin(lMargin, rMargin, bMargin, tMargin)

    #     nullSols.GetYaxis().SetTitle("Frequency")
    #     nullSols.GetXaxis().SetTitle("Marginalization Solution Pulls")
    #     RatioAxis(nullSols,MNVPLOTTER)

    #     nullSols.SetLineWidth(2)
    #     nullSols.SetLineColor(ROOT.kRed)
    #     oscSols.SetLineWidth(2)
    #     oscSols.SetLineColor(ROOT.kBlue)

    #     nullSols.Draw("hist")
    #     oscSols.Draw("hist same")

    #     leg1 = ROOT.TLegend(.3,.2)
    #     leg1.AddEntry(nullSols,"Null Hypothesis Pulls","l")
    #     leg1.AddEntry(oscSols,"Osc. Model Pulls","l")
    #     leg1.Draw()

    #     # ------------------------------------ Second Pad -----------------------------------------
    #     # plot mc - flux universe / mc

    #     c1.cd(2)
    #     pad = ROOT.gPad
    #     pad.SetMargin(lMargin, rMargin, bMargin, tMargin)

    #     straightLine = h_null.Clone()
    #     for whichBin in range(0, straightLine.GetXaxis().GetNbins()+1): 
    #         straightLine.SetBinContent(whichBin, 0)

    #     RatioAxis(straightLine,MNVPLOTTER)
    #     straightLine.SetMinimum(-.3)
    #     straightLine.SetMaximum(.3)
    #     straightLine.GetYaxis().SetTitle("#frac{Universe - CV}{CV}")

    #     straightLine.Draw("hist l")

    #     cv = h_null.GetCVHistoWithError()
    #     h_test = h_null.Clone()
    #     errband = h_test.GetVertErrorBand("Flux")
    #     for i in range(errband.GetNHists()):
    #         universe = errband.GetHist(i)
    #         universe.Add(cv,-1)
    #         universe.Divide(universe,cv)
    #         universe.SetMarkerSize(0)
    #         universe.DrawClone("same hist l")

    #     straightLine.Draw("hist l same")

    #     leg2 = ROOT.TLegend(.3,.2)
    #     leg2.AddEntry(straightLine,"Central Value","l")
    #     #leg2.Draw()

    #     # ------------------------------------ Third Pad -----------------------------------------
    #     # plot solution *  [mc - flux universe / mc ]

    #     c1.cd(3)
    #     pad = ROOT.gPad
    #     pad.SetMargin(lMargin, rMargin, bMargin, tMargin)

    #     line = straightLine.Clone()
    #     line.GetYaxis().SetTitle("#frac{solution#times(Universe - CV)}{CV}")
    #     line.Draw("hist l")

    #     cv = h_null.GetCVHistoWithError()
    #     errband = h_null.GetVertErrorBand("Flux")
    #     for i in range(errband.GetNHists()):
    #         universe = errband.GetHist(i)
    #         universe.Add(cv,-1)
    #         universe.Scale(nullSolution[i])
    #         universe.Divide(universe,cv)
    #         universe.SetMarkerSize(0)
    #         universe.DrawClone("same hist l")

    #     line.Draw("hist l same")

    #     leg3 = ROOT.TLegend(.3,.2)
    #     leg3.AddEntry(straightLine,"Central Value","l")
    #     leg3.Draw()

    #     c1.Print("plots/{}_marg_effects_{:.1f}_{:.3f}_{:.4f}.png".format(name,parameters["m"],parameters['ue4'],parameters['umu4']))

    # def PlotFluxProfilingEffects(self,name=""):
    #     histogram = copy.deepcopy(self.histogram)
    #     exclude = self.exclude
    #     lam = self.lam

    #     statistic = Statistics(histogram,exclude=exclude,lam=lam)
    #     chi2_model,model_pen = statistic.Chi2DataMC(marginalize=True)
    #     chi2_null,null_pen = statistic.Chi2DataMC(marginalize=False)

    #     h_null = histogram.GetMCHistogram()
    #     h_data = histogram.GetDataHistogram()

    #     h_prof = histogram.GetMCHistogram()
    #     statistic.GetFluxFitter(useOsc=False).ReweightToFluxSolution(h_prof)

    #     c1 = ROOT.TCanvas()
    #     margin = .12
    #     bottomFraction = .2
    #     overall = ROOT.TCanvas("Data/MC")
    #     top = ROOT.TPad("DATAMC", "DATAMC", 0, bottomFraction, 1, 1)
    #     bottom = ROOT.TPad("Ratio", "Ratio", 0, 0, 1, bottomFraction+margin)

    #     top.Draw()
    #     bottom.Draw()

    #     top.cd()
    #     top.SetLogy()
        
    #     h_null.SetTitle(name)
    #     h_prof.SetLineColor(ROOT.kBlue)
        
    #     h_null.Draw("hist")
    #     h_prof.Draw("hist same")
    #     h_data.Draw("same")

    #     null = h_null.GetCVHistoWithError()
    #     null.SetLineColor(ROOT.kRed)
    #     null.SetLineWidth(2)
    #     null.SetMarkerStyle(0)
    #     null.SetFillColorAlpha(ROOT.kPink + 1, 0.3)
    #     null.Draw("E2 SAME")

    #     h_prof.PopVertErrorBand("Flux")
    #     prof = h_prof.GetCVHistoWithError()
    #     prof.SetLineColor(ROOT.kBlue)
    #     prof.SetLineWidth(2)
    #     prof.SetMarkerStyle(0)
    #     prof.SetFillColorAlpha(ROOT.kBlue + 1, 0.3)
    #     prof.Draw("E2 SAME")

    #     leg = ROOT.TLegend(.28,.35)
    #     leg.AddEntry(h_data,"Data","p")
    #     top_text = "Null Hypothesis"
    #     bot_text = "#chi^{2}="+"{:.2f}".format(chi2_null)
    #     leg_text = "#splitline{%s}{%s}" % (top_text,bot_text)
    #     leg.AddEntry(h_null,leg_text,"l")
    #     top_text = "Profiled Flux"
    #     bot_text = "#chi^{2}="+"{:.2f} + {:.2f} penalty".format(chi2_model-model_pen,model_pen)
    #     leg_text = "#splitline{%s}{%s}" % (top_text,bot_text)
    #     leg.AddEntry(h_prof,leg_text,"l")
    #     leg.Draw()

    #     nullRatio =  h_data.Clone()
    #     profRatio =  h_prof.Clone()

    #     nullRatio.Divide(nullRatio,h_null)
    #     profRatio.Divide(profRatio, h_null)

    #     bottom.cd()
    #     bottom.SetTopMargin(0)
    #     bottom.SetBottomMargin(0.3)

    #     nullErrors = h_null.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    #     for whichBin in range(0, nullErrors.GetXaxis().GetNbins()+1): 
    #         nullErrors.SetBinError(whichBin, max(nullErrors.GetBinContent(whichBin), 1e-9))
    #         nullErrors.SetBinContent(whichBin, 1)
    #     profErrors = h_prof.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    #     for whichBin in range(0, profErrors.GetXaxis().GetNbins()+1): 
    #         profErrors.SetBinError(whichBin, max(profErrors.GetBinContent(whichBin), 1e-9))
    #         profErrors.SetBinContent(whichBin,profRatio.GetBinContent(whichBin))

    #     nullRatio.SetLineColor(ROOT.kBlack)
    #     nullRatio.SetLineWidth(3)

    #     #Error envelope for the MC
    #     nullErrors.SetLineWidth(0)
    #     nullErrors.SetMarkerStyle(0)
    #     nullErrors.SetFillColorAlpha(ROOT.kPink + 1, 0.3)
    #     nullErrors.GetYaxis().SetTitle("#splitline{Ratio to Null}{Hypothesis}")
    #     RatioAxis(nullErrors,MNVPLOTTER)
    #     nullErrors.GetXaxis().SetTitle("Bin Number")
    #     nullErrors.SetMinimum(.7)
    #     nullErrors.SetMaximum(1.3)
    #     profErrors.SetLineWidth(0)
    #     profErrors.SetMarkerStyle(0)
    #     profErrors.SetFillColorAlpha(ROOT.kBlue + 1, 0.3)
    #     nullErrors.Draw("E2")
    #     profErrors.Draw("E2 same")

    #     #Draw the data ratios
    #     nullRatio.Draw("same")
    #     profRatio.Draw('same hist l')

    #     #Draw a flat line at 1 for profRatio of MC to itself
    #     straightLine = nullErrors.Clone()
    #     straightLine.SetLineColor(ROOT.kRed)
    #     straightLine.SetLineWidth(2)
    #     straightLine.SetFillColor(0)
    #     straightLine.Draw("HIST SAME")

    #     top.cd()

    #     overall.Print("plots/{}_stitched.png".format(name))

    # def PlotOscillationEffects(self,parameters,name="",useMarg=True,plotSamples=False,usePseudo=False):
    #     histogram = copy.deepcopy(self.histogram)
    #     exclude = self.exclude
    #     lam = self.lam

    #     statistic = Statistics(histogram,exclude=exclude,lam=lam)
    #     chi2_model,model_pen = statistic.Chi2DataMC(marginalize=useMarg,useOsc=True,usePseudo=usePseudo)
    #     chi2_null,null_pen = statistic.Chi2DataMC(marginalize=useMarg,usePseudo=usePseudo)

    #     plot_texts = ["Null Hyp. #chi^{2} = "+"{:.2f} + {:.2f}".format(chi2_null-null_pen,null_pen),
    #             "Osc. Model #chi^{2} = "+"{:.2f} + {:.2f}".format(chi2_model-model_pen,model_pen),
    #             "#Delta m^{2} = "+"{}".format(parameters["m"])+" eV^{2}",
    #             "|U_{e4}|^{2} = "+"{}".format(parameters["ue4"]),
    #             "|U_{#mu4}|^{2} = "+"{:.6f}".format(parameters["umu4"]),
    #             "|U_{#tau4}|^{2} = "+"{}".format(parameters["utau4"])
    #             ]

    #     h_null = histogram.GetMCHistogram()
    #     h_osc = histogram.GetOscillatedHistogram()
    #     h_data = histogram.GetPseudoHistogram() if usePseudo else histogram.GetDataHistogram()
    #     h_prof = h_null.Clone()
    #     h_prof.SetLineColor(ROOT.kGreen)
    #     statistic.GetFluxFitter(useOsc=False).ReweightToFluxSolution(h_prof)
    #     statistic.GetFluxFitter(useOsc=True).ReweightToFluxSolution(h_osc)

    #     c1 = ROOT.TCanvas()
    #     margin = .12
    #     bottomFraction = .2
    #     overall = ROOT.TCanvas("Data/MC")
    #     top = ROOT.TPad("DATAMC", "DATAMC", 0, bottomFraction, 1, 1)
    #     bottom = ROOT.TPad("Ratio", "Ratio", 0, 0, 1, bottomFraction+margin)

    #     top.Draw()
    #     bottom.Draw()

    #     top.cd()
    #     top.SetLogy()
        
    #     h_null.SetTitle(name)
        
    #     h_null.Draw("hist")
    #     h_prof.Draw("hist")
    #     h_osc.Draw("hist same")
    #     h_data.Draw("same")

    #     h_null.PopVertErrorBand("Flux")
    #     null = h_null.GetCVHistoWithError()
    #     null.SetLineColor(ROOT.kRed)
    #     null.SetLineWidth(2)
    #     null.SetMarkerStyle(0)
    #     null.SetFillColorAlpha(ROOT.kPink + 1, 0.3)
    #     null.Draw("E2 SAME")

    #     h_osc.PopVertErrorBand("Flux")
    #     osc = h_osc.GetCVHistoWithError()
    #     osc.SetLineColor(ROOT.kBlue)
    #     osc.SetLineWidth(2)
    #     osc.SetMarkerStyle(0)
    #     osc.SetFillColorAlpha(ROOT.kBlue + 1, 0.3)
    #     osc.Draw("E2 SAME")

    #     top_text = "#Delta m^{2}="+"{:.1f}".format(parameters["m"])+" |U_{e4}|^{2}="+"{:.4f}".format(parameters["ue4"])
    #     bot_text = "|U_{#mu4}|^{2}="+"{:.5f}".format(parameters["umu4"])+" |U_{#tau4}|^{2}="+"{:.1f}".format(parameters["utau4"])
    #     header = "#splitline{%s}{%s}" % (top_text,bot_text)
    #     MNVPLOTTER.AddPlotLabel(header,.35,.25)

    #     leg = ROOT.TLegend(.28,.35)

    #     leg.AddEntry(h_data,"Data","p")
    #     top_text = "Null Hypothesis"
    #     bot_text = "#chi^{2}="+"{:.2f} + {:.2f} penalty".format(chi2_null-null_pen,null_pen)
    #     leg_text = "#splitline{%s}{%s}" % (top_text,bot_text)
    #     leg.AddEntry(h_null,leg_text,"l")
    #     top_text = "Best Osc. Fit"
    #     bot_text = "#chi^{2}="+"{:.2f} + {:.2f} penalty".format(chi2_model-model_pen,model_pen)
    #     leg_text = "#splitline{%s}{%s}" % (top_text,bot_text)
    #     leg.AddEntry(h_osc,leg_text,"l")
    #     leg.Draw()

    #     nullRatio =  h_data.Clone()
    #     oscRatio =  h_osc.Clone()
    #     profRatio = h_prof.Clone()

    #     nullRatio.Divide(nullRatio,h_null)
    #     oscRatio.Divide(oscRatio, h_null)
    #     profRatio.Divide(profRatio, h_null)

    #     bottom.cd()
    #     bottom.SetTopMargin(0)
    #     bottom.SetBottomMargin(0.3)

    #     nullErrors = h_null.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    #     for whichBin in range(0, nullErrors.GetXaxis().GetNbins()+1): 
    #         nullErrors.SetBinError(whichBin, max(nullErrors.GetBinContent(whichBin), 1e-9))
    #         nullErrors.SetBinContent(whichBin, 1)

    #     nullRatio.SetLineColor(ROOT.kBlack)
    #     nullRatio.SetLineWidth(3)

    #     profRatio.SetLineWidth(3)

    #     #Error envelope for the MC
    #     nullErrors.SetLineWidth(0)
    #     nullErrors.SetMarkerStyle(0)
    #     nullErrors.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
    #     nullErrors.GetYaxis().SetTitle("#splitline{Ratio to Null}{Hypothesis}")
    #     RatioAxis(nullErrors,MNVPLOTTER)
    #     nullErrors.GetXaxis().SetTitle("Bin Number")
    #     nullErrors.SetMinimum(.7)
    #     nullErrors.SetMaximum(1.3)
    #     nullErrors.Draw("E2")

    #     #Draw the data ratios
    #     nullRatio.Draw("same")
    #     oscRatio.Draw('same hist l')
    #     profRatio.Draw("same hist l")

    #     #Draw a flat line at 1 for oscRatio of MC to itself
    #     straightLine = nullErrors.Clone()
    #     straightLine.SetLineColor(ROOT.kRed)
    #     straightLine.SetLineWidth(2)
    #     straightLine.SetFillColor(0)
    #     straightLine.Draw("HIST SAME")

    #     top.cd()

    #     overall.Print("plots/{}_stitched_{:.1f}_{:.3f}_{:.4f}.png".format(name,parameters["m"],parameters['ue4'],parameters['umu4']))

    #     if not plotSamples:
    #         return

    #     plots = histogram.keys
    #     histogram = copy.deepcopy(self.histogram)
    #     nullSolution = statistic.GetFluxFitter(useOsc=False).GetFluxSolution()
    #     oscSolution = statistic.GetFluxFitter(useOsc=True).GetFluxSolution()
    #     #plots.append("fhc_ratio")
    #     #plots.append("rhc_ratio")
    #     #histogram.titles["fhc_ratio"] = "FHC CC #nu_{#mu}/#nu_{e} Ratio"
    #     #histogram.titles["rhc_ratio"] = "RHC CC anti #nu_{#mu}/#nu_{e} Ratio"

    #     mc_hists = []
    #     null_hists = []
    #     osc_hists = []
    #     data_hists = []
    #     titles = []

    #     plots = [item for item in plots if "nue" not in item]
    #     plots.insert(3, plots.pop(6))

    #     for i,plot in enumerate(plots):
    #         title = histogram.titles[plot]
    #         titles.append(title)
    #         margin = .12
    #         bottomFraction = .2
    #         overall = ROOT.TCanvas(plot)
    #         top = ROOT.TPad("DATAMC", "DATAMC", 0, bottomFraction, 1, 1)
    #         bottom = ROOT.TPad("Ratio", "Ratio", 0, 0, 1, bottomFraction+margin)

    #         top.Draw()
    #         bottom.Draw()

    #         top.cd()

    #         hists = []
    #         if 'ratio' in plot:
    #             numu_hists,numu_total_hist = histogram.OscillateSubHistogram(plot[:3]+'_numu_selection',parameters["m"],parameters['ue4'],parameters['umu4'],parameters['utau4'])
    #             nue_hists,nue_total_hist = histogram.OscillateSubHistogram(plot[:3]+'_nue_selection',parameters["m"],parameters['ue4'],parameters['umu4'],parameters['utau4'])

    #             total_hist = numu_total_hist.Clone()
    #             total_hist.Divide(total_hist,nue_total_hist)
    #             swap_fraction = total_hist.Clone()
    #             norm_fraction = total_hist.Clone()
    #             numu_rat_hist = numu_total_hist.Clone()
    #             nue_rat_hist = nue_total_hist.Clone()

    #             for m in range(numu_total_hist.GetNbinsX()+1):
    #                 frac = total_hist.GetBinContent(m)/numu_rat_hist.GetBinContent(m) if numu_rat_hist.GetBinContent(m) != 0 else 0
    #                 numu_rat_hist.SetBinContent(m,frac*total_hist.GetBinContent(m))
    #                 nue_rat_hist.SetBinContent(m,(1-frac)*total_hist.GetBinContent(m))
    #                 for u in range(nue_rat_hist.GetVertErrorBand("Flux").GetNHists()):
    #                     numu_rat_hist.GetVertErrorBand("Flux").GetHist(u).SetBinContent(m,numu_rat_hist.GetBinContent(m)*frac)
    #                     nue_rat_hist.GetVertErrorBand("Flux").GetHist(u).SetBinContent(m,nue_rat_hist.GetBinContent(m)*(1-frac))

    #             numu_rat_hist.SetTitle("numu")
    #             nue_rat_hist.SetTitle("nue")
    #             hists.append(numu_rat_hist)
    #             hists.append(nue_rat_hist)

    #             """
    #             for j,nue in enumerate(nue_hists):
    #                 if "#mu" in nue.GetTitle():
    #                     swap_fraction.Divide(nue,swap_fraction)
    #                     norm_fraction.Add(swap_fraction,-1)

    #                     norm_fraction.SetTitle("#nu_{#mu}/#nu_{e}")
    #                     norm_fraction.SetFillStyle(3244)
    #                     norm_fraction.SetLineColor(ROOT.kRed)
    #                     norm_fraction.SetFillColor(ROOT.kRed)

    #                     swap_fraction.SetTitle("#nu_{#mu}/"+nue.GetTitle())
    #                     swap_fraction.SetFillStyle(3244)
    #                     swap_fraction.SetLineColor(ROOT.kBlue)
    #                     swap_fraction.SetFillColor(ROOT.kBlue)

    #                     hists.append(norm_fraction)
    #                     hists.append(swap_fraction)
    #             """
    #             h_data = histogram.data_hists[plot[:3]+'_numu_selection'].Clone()
    #             h_data.Divide(h_data,histogram.data_hists[plot[:3]+'_nue_selection'])
    #             subSample = histogram.mc_hists[plot[:3]+'_numu_selection'].Clone()
    #             subSample.Divide(subSample,histogram.mc_hists[plot[:3]+'_nue_selection'])
    #         else:
    #             hists,total_hist = histogram.OscillateSubHistogram(plot,parameters["m"],parameters['ue4'],parameters['umu4'],parameters['utau4'])
    #             h_data = histogram.data_hists[plot].Clone()
    #             subSample = histogram.mc_hists[plot].Clone()

    #         cv = np.array(subSample)[1:-1]
    #         mc = np.array(total_hist)[1:-1]

    #         null_hist = histogram.mc_hists[plot].Clone()

    #         statistic.GetFluxFitter(useOsc=False).ReweightToFluxSolution(null_hist)
    #         null_hist.PopVertErrorBand("Flux")

    #         statistic.GetFluxFitter(useOsc=True).ReweightToFluxSolution(total_hist)
    #         total_hist.PopVertErrorBand("Flux")
    #         total_hist.AddMissingErrorBandsAndFillWithCV(h_data)

    #         TArray = []
    #         for hist in hists:
    #             if hist.Integral() <= 0:
    #                 continue
    #             statistic.GetFluxFitter(useOsc=True).ReweightToFluxSolution(hist)
    #             if 'elastic' in plot:
    #                 hist.Scale(2,'width')
    #             elif 'ratio' not in plot:
    #                 hist.Scale(1,'width')

    #             if "tau" in hist.GetTitle():
    #                 hist.SetFillColorAlpha(ROOT.kGray, 0.6)
    #             elif "ratio" in plot:
    #                 hist.SetFillColorAlpha(ROOT.kViolet,0.6)
    #             elif "mu" in hist.GetTitle():
    #                 hist.SetFillColorAlpha(ROOT.kBlue, 0.6)
    #             else:
    #                 hist.SetFillColorAlpha(ROOT.kRed, 0.6)
    #             hist.SetLineWidth(0)
    #             TArray.append(hist)

    #         if "elastic" in plot:
    #             h_data.Scale(2,'width')
    #             total_hist.Scale(2,'width')
    #             null_hist.Scale(2,'width')
    #         elif 'ratio' not in plot:
    #             h_data.Scale(1,'width')
    #             total_hist.Scale(1,'width')
    #             null_hist.Scale(1,'width')

    #         if "elastic" in plot:
    #             total_hist.SetTitle("#nu+e")
    #         elif "imd" in plot:
    #             total_hist.SetTitle("#nu_{#mu}+e^{-}#rightarrow #mu^{-}+#nu_{e}")
    #         elif "mu" in plot:
    #             total_hist.SetTitle("CC #nu_{#mu}")
    #         elif "ratio" in plot:
    #             total_hist.SetTitle("CC #nu_{#mu}/#nu_{e}")
    #         else:
    #             total_hist.SetTitle("CC #nu_{e}")

    #         if "fhc" in plot:
    #             total_hist.GetYaxis().SetTitle("#nu-mode    Events/GeV")
    #         else:
    #             total_hist.GetYaxis().SetTitle("#bar{#nu}-mode    Events/GeV")

    #         mc_hists.append(TArray)
    #         osc_hists.append(total_hist)
    #         null_hists.append(null_hist)
    #         data_hists.append(h_data)
        
    #     overall = plot_osc_side_by_side(mc_hists,null_hists,osc_hists,data_hists,titles,plot_texts,MNVPLOTTER,narrow_pads=[1,4])
    #     overall.Print("plots/{}_oscillated_{:.1f}_{:.3f}_{:.4f}.png".format(name,parameters["m"],parameters['ue4'],parameters['umu4']))

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
    def __init__(self, tag, histogram, hist_config="HIST_CONFIG.json"):
        self.tag = tag
        self.histogram = histogram
        self.hist_config = hist_config
        self.exclude = ""
        self.lam = 1
        self.flux_solution = None
        self.osc_parameters = {"m":0,"ue4":0,"umu4":0,"utau4":0}
        self.invCov = histogram.GetInverseCovarianceMatrix(sansFlux=True)

        self.exclude_samples = ["", "", "ratio"]
        self.titles = ["#nu+e, IMD, CC #nu_{#mu}, CC #nu_{#mu}/#nu_{e} #lambda = 1","#nu+e, IMD, CC #nu_{#mu}, CC #nu_{#mu}/#nu_{e} #lambda = 12", "#nu+e, IMD, CC #nu_{#mu} #lambda = 1"]
        self.colors = [ROOT.kRed,ROOT.kBlue,ROOT.kGreen]
        self.lams = [1,12,1]
        self.mask_spec = None
        self.mask_bins = []
        self.profile_only = None
        self.profile_n_universes = None


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

    def SetMaskSpec(self, mask_spec):
        self.mask_spec = mask_spec
        self.mask_bins = self.GetMaskedBinIndicesForPlot(mask_spec)

    def SetProfileOnly(self, profile_only):
        self.profile_only = profile_only

    def SetProfileNUniverses(self, profile_n_universes):
        self.profile_n_universes = profile_n_universes

    def GetMaskedBinIndicesForPlot(self, mask_spec):
        if mask_spec is None:
            return []

        import json

        with open(self.hist_config, "r") as f:
            cfg = json.load(f)

        masked = []

        for sample, local_bins in mask_spec.items():
            if sample not in cfg:
                print(
                    "WARNING: sample {} not found in {}".format(
                        sample,
                        self.hist_config,
                    )
                )
                continue

            start = cfg[sample]["start"]  # zero-based

            for local_bin in local_bins:
                masked.append(start + local_bin - 1)

        masked = sorted(set(masked))

        print("\n===== Plot mask bins =====")
        print("mask_spec =", mask_spec)
        print("masked zero-based bins =", masked)
        print("masked ROOT bins       =", [x + 1 for x in masked])

        return masked

    def ApplyExternalCovarianceErrorsForPlot(
        self,
        hist,
        reference_hist=None,
        ratio_mode=False,
        set_content_to_one=False,
    ):
        """
        Plotting-only helper.

        For top panel:
            ratio_mode=False
            error = sqrt(cov_ii)

        For ratio panel:
            ratio_mode=True
            error = sqrt(cov_ii) / reference_hist bin content

        This only uses the diagonal of the covariance for drawing.
        The chi2 still uses the full covariance matrix.
        """
        h_out = hist.Clone(hist.GetName() + "_with_external_cov_errors")
        h_out.SetDirectory(0)

        if not hasattr(self.histogram, "external_covariances"):
            return h_out

        for sample_name, root_cov in self.histogram.external_covariances.items():
            if root_cov is None:
                continue

            inds = self.histogram.GetExternalCovarianceBinIndices(
                sample_name,
                hist_config=self.hist_config,
            )
            if len(inds) == 0:
                continue

            cov = TMatrix_to_Numpy(root_cov)

            if cov.shape[0] != len(inds) or cov.shape[1] != len(inds):
                raise RuntimeError(
                    "External covariance shape {} does not match {} bins for {}".format(
                        cov.shape, len(inds), sample_name
                    )
                )

            for a, idx0 in enumerate(inds):
                ibin = idx0 + 1
                sigma = math.sqrt(max(cov[a, a], 0.0))

                if ratio_mode:
                    if reference_hist is None:
                        raise RuntimeError("ratio_mode=True requires reference_hist")

                    denom = reference_hist.GetBinContent(ibin)
                    err = sigma / denom if denom != 0 else 0.0

                    h_out.SetBinError(ibin, err)

                    if set_content_to_one:
                        h_out.SetBinContent(ibin, 1.0)

                else:
                    h_out.SetBinError(ibin, sigma)

        return h_out

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
        statistic = Statistics(
            histogram,
            exclude=self.exclude,
            lam=self.lam,
            mask_spec=self.mask_spec,
            profile_only=self.profile_only,
            profile_n_universes=self.profile_n_universes,
            hist_config=self.hist_config,
        )      
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
        statistic = Statistics(
            histogram,
            exclude=self.exclude,
            lam=self.lam,
            mask_spec=self.mask_spec,
            profile_only=self.profile_only,
            profile_n_universes=self.profile_n_universes,
            hist_config=self.hist_config,
        )    
        statistic.Chi2DataMC(marginalize=True)

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

    def PlotOscillationEffects(
        self,
        res,
        ntuple_tag,
        useMarg=False,
        plotSamples=False,
        plot_tag="",
        usePseudo=False,
    ):
        histogram = copy.deepcopy(self.histogram)
        exclude = self.exclude
        lam = self.lam


        def make_mask_tag(mask_spec):
            if mask_spec is None:
                return ""

            parts = []
            for sample, bins in sorted(mask_spec.items()):
                bin_str = "-".join([str(b) for b in bins])
                parts.append("{}_bins{}".format(sample, bin_str))

            return "_mask_" + "__".join(parts)

        def make_profile_tag(exclude):
            if exclude is None or exclude == "" or exclude == "none":
                return "_profiledFlux_includeAll"

            safe = str(exclude).replace(",", "-").replace(" ", "")
            return "_profiledFlux_exclude{}".format(safe)


        import json

        def print_sample_ratio(label, h_ref, h_test):
            print("\n===== {} =====".format(label))

            with open(self.hist_config, "r") as f:
                cfg = json.load(f)

            for sample, info in cfg.items():
                start = info["start"] + 1  # ROOT bin index
                end   = info["end"] + 1

                ref_sum = 0.0
                test_sum = 0.0

                for ibin in range(start, end + 1):
                    ref_sum  += h_ref.GetBinContent(ibin)
                    test_sum += h_test.GetBinContent(ibin)

                ratio = test_sum / ref_sum if ref_sum != 0 else 0.0

                print("{:25s} ref={:14.6g} test={:14.6g} test/ref={:12.6f}".format(
                    sample, ref_sum, test_sum, ratio
                ))

        def print_flux_shift_by_sample(label, h_before, h_after):
            print("\n===== {} =====".format(label))

            with open(self.hist_config, "r") as f:
                cfg = json.load(f)

            for sample, info in cfg.items():
                start = info["start"] + 1  # ROOT bin index
                end   = info["end"] + 1

                before_sum = 0.0
                after_sum = 0.0

                for ibin in range(start, end + 1):
                    before_sum += h_before.GetBinContent(ibin)
                    after_sum  += h_after.GetBinContent(ibin)

                shift = after_sum - before_sum
                frac = shift / before_sum if before_sum != 0 else 0.0

                print("{:25s} before={:14.6g} after={:14.6g} shift={:14.6g} frac={:12.6f}".format(
                    sample, before_sum, after_sum, shift, frac
                ))

        def draw_masked_bin_lines(mask_bins, ymin=0.5, ymax=1.5):
            lines = []
            for idx0 in mask_bins:
                x = idx0 + 1
                line = ROOT.TLine(x, ymin, x, ymax)
                line.SetLineColor(ROOT.kGray + 2)
                line.SetLineStyle(2)
                line.SetLineWidth(2)
                line.Draw("SAME")
                lines.append(line)
            return lines

        def draw_masked_bin_boxes(mask_bins, axis_hist, ymin=0.5, ymax=1.5):
            boxes = []

            for idx0 in mask_bins:
                ibin = idx0 + 1  # ROOT bin index

                x1 = axis_hist.GetXaxis().GetBinLowEdge(ibin)
                x2 = axis_hist.GetXaxis().GetBinUpEdge(ibin)

                box = ROOT.TBox(x1, ymin, x2, ymax)
                box.SetFillColorAlpha(ROOT.kGray + 3, 0.25)
                box.SetLineColor(ROOT.kGray + 4)
                box.SetLineStyle(2)
                box.SetLineWidth(1)
                box.Draw("SAME")

                boxes.append(box)

            return boxes



        parameters = res
        name = ntuple_tag


        h_stored_before_reset = histogram.GetOscillatedHistogram().Clone(
            "h_stored_before_reset"
        )

        histogram.OscillateHistogram(
            parameters["m"],
            parameters["ue4"],
            parameters["umu4"],
            parameters["utau4"],
            False,
            False,
        )

        h_recomputed_bf = histogram.GetOscillatedHistogram()

        stored = np.array(h_stored_before_reset)[1:-1]
        recomputed = np.array(h_recomputed_bf)[1:-1]

        ratio = np.divide(
            stored,
            recomputed,
            out=np.ones_like(stored),
            where=np.abs(recomputed) > 1e-12,
        )

        print("\n===== stored oscillation vs explicitly recomputed BF =====")
        print("max |stored/recomputed - 1| =", np.max(np.abs(ratio - 1.0)))
        print("mean |stored/recomputed - 1| =", np.mean(np.abs(ratio - 1.0)))


        statistic = Statistics(
            histogram,
            exclude=self.exclude,
            lam=self.lam,
            mask_spec=self.mask_spec,
            profile_only=self.profile_only,
            profile_n_universes=self.profile_n_universes,
            hist_config=self.hist_config,
        )


        # Raw/unprofiled chi2 values.
        # These correspond to the unprofiled red null curve and raw oscillated model.
        chi2_null_raw, null_pen_raw = statistic.Chi2DataMC(
            marginalize=False,
            usePseudo=usePseudo
        )
        chi2_model_raw, model_pen_raw = statistic.Chi2DataMC(
            marginalize=False,
            useOsc=True,
            usePseudo=usePseudo
        )

        # Profiled chi2 values.
        # These correspond to the green profiled null curve and profiled oscillated model.
        if useMarg:
            chi2_null, null_pen = statistic.Chi2DataMC(
                marginalize=True,
                usePseudo=usePseudo
            )
            chi2_model, model_pen = statistic.Chi2DataMC(
                marginalize=True,
                useOsc=True,
                usePseudo=usePseudo
            )
        else:
            chi2_null, null_pen = chi2_null_raw, null_pen_raw
            chi2_model, model_pen = chi2_model_raw, model_pen_raw

        if useMarg:
            null_chi2_text = "Raw Null #chi^{{2}} = {:.2f}".format(chi2_null_raw)
            prof_null_chi2_text = "Profiled Null #chi^{{2}} = {:.2f} + {:.2f}".format(
                chi2_null - null_pen, null_pen
            )
            osc_chi2_text = "Profiled Osc. #chi^{{2}} = {:.2f} + {:.2f}".format(
                chi2_model - model_pen, model_pen
            )
            plot_texts = [
                null_chi2_text,
                prof_null_chi2_text,
                osc_chi2_text,
                "#Delta m^{2} = " + "{}".format(parameters["m"]) + " eV^{2}",
                "|U_{e4}|^{2} = " + "{}".format(parameters["ue4"]),
                "|U_{#mu4}|^{2} = " + "{:.6f}".format(parameters["umu4"]),
                "|U_{#tau4}|^{2} = " + "{}".format(parameters["utau4"]),
            ]
        else:
            null_chi2_text = "Null Hyp. #chi^{{2}} = {:.2f}".format(chi2_null_raw)
            osc_chi2_text = "Osc. Model #chi^{{2}} = {:.2f}".format(chi2_model_raw)
            plot_texts = [
                null_chi2_text,
                osc_chi2_text,
                "#Delta m^{2} = " + "{}".format(parameters["m"]) + " eV^{2}",
                "|U_{e4}|^{2} = " + "{}".format(parameters["ue4"]),
                "|U_{#mu4}|^{2} = " + "{:.6f}".format(parameters["umu4"]),
                "|U_{#tau4}|^{2} = " + "{}".format(parameters["utau4"]),
            ]

        if self.mask_spec is not None:
            mask_label = "Masked: " + ", ".join(
                ["{} {}".format(k, v) for k, v in self.mask_spec.items()]
            )
            plot_texts.append(mask_label)

        h_null = histogram.GetMCHistogram()
        h_osc = histogram.GetOscillatedHistogram()
        h_data = histogram.GetPseudoHistogram() if usePseudo else histogram.GetDataHistogram()


        h_osc_raw = h_osc.Clone("h_osc_raw_before_flux_profile")
        h_osc_raw.SetLineColor(ROOT.kOrange + 7)
        h_osc_raw.SetLineStyle(2)
        h_osc_raw.SetLineWidth(3)


        def build_A_from_hist_flux_band(hist, band_name="Flux"):
            if not hist.HasVertErrorBand(band_name):
                raise RuntimeError("hist has no {} band: {}".format(band_name, hist.GetName()))

            band = hist.GetVertErrorBand(band_name)
            nhists = band.GetNHists()

            cv = np.array(hist)[1:-1]
            universes = np.array([
                np.array(band.GetHist(i))[1:-1]
                for i in range(nhists)
            ])

            return universes - cv[None, :]

        def print_direct_Aa_bins(label, Aa_null, Aa_bf, sliceInds, bins_to_print):
            print("\n===== {} =====".format(label))
            print("{:>20s} {:>6s} {:>6s} {:>12s} {:>12s} {:>12s}".format(
                "sample", "gbin", "lbin", "Aa_null", "Aa_BF", "diff"
            ))

            import json
            with open(self.hist_config, "r") as f:
                cfg = json.load(f)

            # Map global zero-based bin -> position inside sliced A
            pos = {idx0: j for j, idx0 in enumerate(sliceInds)}

            for sample, local_bins in bins_to_print.items():
                if sample not in cfg:
                    continue

                start = cfg[sample]["start"]  # zero-based global start

                for lbin in local_bins:
                    idx0 = start + lbin - 1

                    if idx0 not in pos:
                        print("{:>20s} {:6d} {:6d} {:>12s} {:>12s} {:>12s}".format(
                            sample, idx0 + 1, lbin, "excluded", "excluded", "excluded"
                        ))
                        continue

                    j = pos[idx0]

                    print("{:>20s} {:6d} {:6d} {:12.6g} {:12.6g} {:12.6g}".format(
                        sample,
                        idx0 + 1,
                        lbin,
                        Aa_null[j],
                        Aa_bf[j],
                        Aa_bf[j] - Aa_null[j]
                    ))


        h_prof = None
        if useMarg:
            h_prof = h_null.Clone()
            h_prof.SetLineColor(ROOT.kGreen)

            # statistic.GetFluxFitter(useOsc=False).ReweightToFluxSolution(h_prof)
            # statistic.GetFluxFitter(useOsc=True).ReweightToFluxSolution(h_osc)
            # Use the same stored A matrix as the chi2 calculation.
            statistic.GetFluxFitter(useOsc=False).ReweightToFluxSolutionStoredA(h_prof)
            statistic.GetFluxFitter(useOsc=True).ReweightToFluxSolutionStoredA(h_osc)


            null_fitter = statistic.GetFluxFitter(useOsc=False)
            bf_fitter   = statistic.GetFluxFitter(useOsc=True)

            a_null = null_fitter.GetFluxSolution()
            a_bf   = bf_fitter.GetFluxSolution()


            # Need same bin selection used inside the flux solve.
            # If your FluxFitter has a getter for sliceInds, use that.
            # Otherwise reconstruct it exactly the same way as FluxFitter.
            sliceInds = GetFluxSolveBinIndices(
                histogram.keys,
                exclude,
                profile_only=self.profile_only,
                mask_bins=self.mask_bins,
                hist_config=self.hist_config,
            )

            A_stored = histogram.GetAMatrix()

            if self.profile_n_universes is not None:
                A_stored = A_stored[:self.profile_n_universes, :]

            A_slice = A_stored[:, sliceInds]

            Aa_null = A_slice.T @ a_null
            Aa_bf   = A_slice.T @ a_bf
            dAa     = Aa_bf - Aa_null

            print("\n===== A source comparison =====")
            print("Using stored A matrix for both null and BF plotting.")
            print("A_stored shape =", A_stored.shape)
            print("A_stored norm  =", np.linalg.norm(A_stored[:, sliceInds]))

            print("\n===== direct flux-vector check =====")
            print("lambda =", lam)
            print("|a_null|        =", np.linalg.norm(a_null))
            print("|a_BF|          =", np.linalg.norm(a_bf))
            print("|a_BF-a_null|   =", np.linalg.norm(a_bf - a_null))
            print("max |a_BF-a_null| =", np.max(np.abs(a_bf - a_null)))
            print("penalty_null =", lam * np.dot(a_null, a_null))
            print("penalty_BF   =", lam * np.dot(a_bf, a_bf))

            print("\nA@a comparison on fitted bins:")
            print("|A a_null|      =", np.linalg.norm(Aa_null))
            print("|A a_BF|        =", np.linalg.norm(Aa_bf))
            print("|A(a_BF-a_null)|=", np.linalg.norm(dAa))
            print("max |A(a_BF-a_null)| =", np.max(np.abs(dAa)))

            raw_null_arr  = np.array(h_null)[1:-1]
            prof_null_arr = np.array(h_prof)[1:-1]

            raw_bf_arr  = np.array(h_osc_raw)[1:-1]
            prof_bf_arr = np.array(h_osc)[1:-1]

            hist_Aa_null = prof_null_arr[sliceInds] - raw_null_arr[sliceInds]
            hist_Aa_bf   = prof_bf_arr[sliceInds] - raw_bf_arr[sliceInds]

            print("\n===== A@a vs histogram shift sanity check =====")
            print("max diff null =", np.max(np.abs(Aa_null - hist_Aa_null)))
            print("max diff BF   =", np.max(np.abs(Aa_bf - hist_Aa_bf)))

            print_direct_Aa_bins(
                "direct A@a bin shifts, independent of plotted histograms",
                Aa_null,
                Aa_bf,
                sliceInds,
                {
                    "fhc_numu_selection": [9, 10],
                    "rhc_numu_selection": [9, 10],
                    "fhc_ratio": [9, 10],
                    "rhc_ratio": [9, 10],
                }
            )

        def print_bf_pull_bins(label, h_data, h_null_model, h_bf_model, invCov, max_bins=20):
            data = np.array(h_data)[1:-1]
            mc_null = np.array(h_null_model)[1:-1]
            mc_bf = np.array(h_bf_model)[1:-1]

            diff_null = data - mc_null
            diff_bf = data - mc_bf

            q_null = diff_null * (invCov @ diff_null)
            q_bf = diff_bf * (invCov @ diff_bf)
            dq = q_null - q_bf

            import json
            sample_ranges = []
            try:
                with open(self.hist_config, "r") as f:
                    cfg = json.load(f)
                for sample, info in cfg.items():
                    sample_ranges.append((sample, info["start"], info["end"]))
            except Exception:
                sample_ranges = []

            def sample_name(idx0):
                for sample, start, end in sample_ranges:
                    if start <= idx0 <= end:
                        return sample
                return "unknown"

            print("\n===== BF pull bin diagnostics: {} =====".format(label))
            print("sum q_null =", np.sum(q_null))
            print("sum q_bf   =", np.sum(q_bf))
            print("sum dq     =", np.sum(dq))

            order = np.argsort(dq)[::-1]

            print("\nTop bins where BF improves over null:")
            print("{:>5s} {:>25s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s}".format(
                "bin", "sample", "data", "null", "BF", "dq", "q_null", "q_BF"
            ))

            for idx0 in order[:max_bins]:
                print("{:5d} {:>25s} {:12.6g} {:12.6g} {:12.6g} {:12.6g} {:12.6g} {:12.6g}".format(
                    idx0 + 1,
                    sample_name(idx0),
                    data[idx0],
                    mc_null[idx0],
                    mc_bf[idx0],
                    dq[idx0],
                    q_null[idx0],
                    q_bf[idx0],
                ))

            print("\nTop bins where BF worsens relative to null:")
            order_bad = np.argsort(dq)

            for idx0 in order_bad[:max_bins]:
                print("{:5d} {:>25s} {:12.6g} {:12.6g} {:12.6g} {:12.6g} {:12.6g} {:12.6g}".format(
                    idx0 + 1,
                    sample_name(idx0),
                    data[idx0],
                    mc_null[idx0],
                    mc_bf[idx0],
                    dq[idx0],
                    q_null[idx0],
                    q_bf[idx0],
                ))

        def print_ratio_bins_table(label, h_data, h_raw_null, h_prof_null, h_raw_bf, h_prof_bf):
            print("\n===== {} =====".format(label))

            with open(self.hist_config, "r") as f:
                cfg = json.load(f)

            ratio_samples = [s for s in cfg.keys() if "ratio" in s]

            print(
                "{:>12s} {:>6s} {:>6s} "
                "{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} "
                "{:>12s} {:>12s} {:>12s}".format(
                    "sample", "gbin", "lbin",
                    "data", "rawNull", "green", "rawBF", "blue",
                    "data/green", "blue/green", "rawBF/rawNull"
                )
            )

            for sample in ratio_samples:
                start = cfg[sample]["start"] + 1  # ROOT bin index
                end   = cfg[sample]["end"] + 1

                for ibin in range(start, end + 1):
                    local_bin = ibin - start + 1

                    data = h_data.GetBinContent(ibin)
                    raw_null = h_raw_null.GetBinContent(ibin)
                    green = h_prof_null.GetBinContent(ibin)
                    raw_bf = h_raw_bf.GetBinContent(ibin)
                    blue = h_prof_bf.GetBinContent(ibin)

                    data_over_green = data / green if green != 0 else 0.0
                    blue_over_green = blue / green if green != 0 else 0.0
                    rawbf_over_rawnull = raw_bf / raw_null if raw_null != 0 else 0.0

                    print(
                        "{:>12s} {:6d} {:6d} "
                        "{:12.6g} {:12.6g} {:12.6g} {:12.6g} {:12.6g} "
                        "{:12.6f} {:12.6f} {:12.6f}".format(
                            sample,
                            ibin,
                            local_bin,
                            data,
                            raw_null,
                            green,
                            raw_bf,
                            blue,
                            data_over_green,
                            blue_over_green,
                            rawbf_over_rawnull,
                        )
                    )

        def print_ratio_num_den_bins(label, histogram, parameters, ratio_bins_to_check):
            print("\n===== {} =====".format(label))

            for sample, local_bins in ratio_bins_to_check.items():
                beam = sample[:3]
                numu_name = beam + "_numu_selection"
                nue_name  = beam + "_nue_selection"

                print("\n--- {} from {} / {} ---".format(sample, numu_name, nue_name))

                # Raw null numerator/denominator
                raw_numu_null = histogram.mc_hists[numu_name].Clone()
                raw_nue_null  = histogram.mc_hists[nue_name].Clone()

                # Data numerator/denominator
                data_numu = histogram.data_hists[numu_name].Clone()
                data_nue  = histogram.data_hists[nue_name].Clone()

                # Raw BF numerator/denominator from subhist oscillation
                _, raw_numu_bf = histogram.OscillateSubHistogram(
                    numu_name,
                    parameters["m"],
                    parameters["ue4"],
                    parameters["umu4"],
                    parameters["utau4"]
                )

                _, raw_nue_bf = histogram.OscillateSubHistogram(
                    nue_name,
                    parameters["m"],
                    parameters["ue4"],
                    parameters["umu4"],
                    parameters["utau4"]
                )

                print(
                    "{:>6s} "
                    "{:>12s} {:>12s} {:>12s} "
                    "{:>12s} {:>12s} {:>12s} "
                    "{:>12s} {:>12s} {:>12s}".format(
                        "lbin",
                        "dataNmu", "dataNe", "dataR",
                        "nullNmu", "nullNe", "nullR",
                        "bfNmu", "bfNe", "bfR"
                    )
                )

                for ibin in local_bins:
                    d_numu = data_numu.GetBinContent(ibin)
                    d_nue  = data_nue.GetBinContent(ibin)
                    d_ratio = d_numu / d_nue if d_nue != 0 else 0.0

                    n_numu = raw_numu_null.GetBinContent(ibin)
                    n_nue  = raw_nue_null.GetBinContent(ibin)
                    n_ratio = n_numu / n_nue if n_nue != 0 else 0.0

                    b_numu = raw_numu_bf.GetBinContent(ibin)
                    b_nue  = raw_nue_bf.GetBinContent(ibin)
                    b_ratio = b_numu / b_nue if b_nue != 0 else 0.0

                    print(
                        "{:6d} "
                        "{:12.6g} {:12.6g} {:12.6g} "
                        "{:12.6g} {:12.6g} {:12.6g} "
                        "{:12.6g} {:12.6g} {:12.6g}".format(
                            ibin,
                            d_numu, d_nue, d_ratio,
                            n_numu, n_nue, n_ratio,
                            b_numu, b_nue, b_ratio,
                        )
                    )

        def print_ratio_bin_errors(label, h_data, h_green, h_blue, invCov):
            print("\n===== {} =====".format(label))

            cov = np.linalg.pinv(invCov)
            sigma = np.sqrt(np.maximum(np.diag(cov), 0.0))

            with open(self.hist_config, "r") as f:
                cfg = json.load(f)

            for sample in ["fhc_ratio", "rhc_ratio"]:
                if sample not in cfg:
                    continue

                start = cfg[sample]["start"] + 1
                end = cfg[sample]["end"] + 1

                print("\n--- {} ---".format(sample))
                print("{:>6s} {:>6s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s}".format(
                    "gbin", "lbin", "data", "green", "blue", "sigma", "pull_g", "pull_b"
                ))

                for ibin in range(start, end + 1):
                    idx0 = ibin - 1
                    local = ibin - start + 1

                    data = h_data.GetBinContent(ibin)
                    green = h_green.GetBinContent(ibin)
                    blue = h_blue.GetBinContent(ibin)
                    sig = sigma[idx0]

                    pull_g = (data - green) / sig if sig != 0 else 0.0
                    pull_b = (data - blue) / sig if sig != 0 else 0.0

                    print("{:6d} {:6d} {:12.6g} {:12.6g} {:12.6g} {:12.6g} {:12.6f} {:12.6f}".format(
                        ibin, local, data, green, blue, sig, pull_g, pull_b
                    ))

        def print_direct_sample_bins_table(label, h_data, h_raw_null, h_prof_null, h_raw_bf, h_prof_bf, invCov, samples_to_print):
            print("\n===== {} =====".format(label))

            import json
            import numpy as np

            with open(self.hist_config, "r") as f:
                cfg = json.load(f)

            # Use diagonal errors from the covariance being used for the plotted/profiled comparison.
            cov = np.linalg.pinv(invCov)
            sigma = np.sqrt(np.maximum(np.diag(cov), 0.0))

            print(
                "{:>20s} {:>6s} {:>6s} "
                "{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} "
                "{:>12s} {:>12s} {:>12s} "
                "{:>12s} {:>12s}".format(
                    "sample", "gbin", "lbin",
                    "data", "rawNull", "green", "rawBF", "blue",
                    "data/green", "blue/green", "rawBF/rawNull",
                    "pull_g", "pull_b"
                )
            )

            for sample in samples_to_print:
                if sample not in cfg:
                    print(
                        "WARNING: sample {} not found in {}".format(
                            sample,
                            self.hist_config,
                        )
                    )
                    continue

                start = cfg[sample]["start"] + 1  # ROOT bin index
                end   = cfg[sample]["end"] + 1

                for ibin in range(start, end + 1):
                    idx0 = ibin - 1
                    local_bin = ibin - start + 1

                    data = h_data.GetBinContent(ibin)
                    raw_null = h_raw_null.GetBinContent(ibin)
                    green = h_prof_null.GetBinContent(ibin)
                    raw_bf = h_raw_bf.GetBinContent(ibin)
                    blue = h_prof_bf.GetBinContent(ibin)

                    data_over_green = data / green if green != 0 else 0.0
                    blue_over_green = blue / green if green != 0 else 0.0
                    rawbf_over_rawnull = raw_bf / raw_null if raw_null != 0 else 0.0

                    sig = sigma[idx0]
                    pull_g = (data - green) / sig if sig != 0 else 0.0
                    pull_b = (data - blue) / sig if sig != 0 else 0.0

                    print(
                        "{:>20s} {:6d} {:6d} "
                        "{:12.6g} {:12.6g} {:12.6g} {:12.6g} {:12.6g} "
                        "{:12.6f} {:12.6f} {:12.6f} "
                        "{:12.6f} {:12.6f}".format(
                            sample,
                            ibin,
                            local_bin,
                            data,
                            raw_null,
                            green,
                            raw_bf,
                            blue,
                            data_over_green,
                            blue_over_green,
                            rawbf_over_rawnull,
                            pull_g,
                            pull_b,
                        )
                    )

        if useMarg and h_prof is not None:
            print_ratio_bins_table(
                "all ratio bins: data, raw null, green profiled null, raw BF, blue profiled BF",
                h_data,
                h_null,
                h_prof,
                h_osc_raw,
                h_osc
            )
            print_ratio_num_den_bins(
                "ratio numerator/denominator check for high-pull bins",
                histogram,
                parameters,
                {
                    "fhc_ratio": [1, 5, 6, 7, 8, 9, 10],
                    "rhc_ratio": [1, 4, 11, 12],
                }
            )
            print_ratio_bin_errors(
                "ratio bin diagonal pulls using sansFlux covariance",
                h_data,
                h_prof,
                h_osc,
                histogram.GetInverseCovarianceMatrix(sansFlux=True)
            )
            print_direct_sample_bins_table(
                "direct CCnue samples: data, raw null, green profiled null, raw BF, blue profiled BF",
                h_data,
                h_null,
                h_prof,
                h_osc_raw,
                h_osc,
                histogram.GetInverseCovarianceMatrix(sansFlux=True),
                ["fhc_nue_selection", "rhc_nue_selection"]
            )



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
        # Raw BF before flux profiling, for visualization only.
        h_osc_raw.Draw("hist same")
        h_osc.Draw("hist same")
        h_data.Draw("same")

        if useMarg and h_null.HasVertErrorBand("Flux"):
            h_null.PopVertErrorBand("Flux")
        null = h_null.GetCVHistoWithError()
        null = self.ApplyExternalCovarianceErrorsForPlot(
            null,
            ratio_mode=False,
        )
        null.SetLineColor(ROOT.kRed)
        null.SetLineWidth(2)
        null.SetMarkerStyle(0)
        null.SetFillColorAlpha(ROOT.kPink + 1, 0.3)
        null.Draw("E2 SAME")

        if useMarg and h_osc.HasVertErrorBand("Flux"):
            h_osc.PopVertErrorBand("Flux")
        osc = h_osc.GetCVHistoWithError()
        osc = self.ApplyExternalCovarianceErrorsForPlot(
            osc,
            ratio_mode=False,
        )
        osc.SetLineColor(ROOT.kBlue)
        osc.SetLineWidth(2)
        osc.SetMarkerStyle(0)
        osc.SetFillColorAlpha(ROOT.kBlue + 1, 0.3)
        osc.Draw("E2 SAME")
        # if useMarg and h_null.HasVertErrorBand("Flux"):
        #     h_null.PopVertErrorBand("Flux")
        # null = h_null.GetCVHistoWithError()
        # null.SetLineColor(ROOT.kRed)
        # null.SetLineWidth(2)
        # null.SetMarkerStyle(0)
        # null.SetFillColorAlpha(ROOT.kPink + 1, 0.3)
        # null.Draw("E2 SAME")

        # if useMarg and h_osc.HasVertErrorBand("Flux"):
        #     h_osc.PopVertErrorBand("Flux")
        # osc = h_osc.GetCVHistoWithError()
        # osc.SetLineColor(ROOT.kBlue)
        # osc.SetLineWidth(2)
        # osc.SetMarkerStyle(0)
        # osc.SetFillColorAlpha(ROOT.kBlue + 1, 0.3)
        # osc.Draw("E2 SAME")

        # Redraw line curves on top of uncertainty bands
        h_null.Draw("hist same")
        if useMarg and h_prof is not None:
            h_prof.Draw("hist same")
        h_osc_raw.Draw("hist same")
        h_osc.Draw("hist same")
        h_data.Draw("same")

        bf_label = (
            "#Delta m^{{2}} = {:.1f} eV^{{2}}, "
            "|U_{{e4}}|^{{2}} = {:.3f}, "
            "|U_{{#mu4}}|^{{2}} = {:.3f}, "
            "|U_{{#tau4}}|^{{2}} = {:.3f}"
        ).format(
            parameters["m"],
            parameters["ue4"],
            parameters["umu4"],
            parameters["utau4"],
        )

        MNVPLOTTER.AddPlotLabel(bf_label, .53, .18, 0.04)



        # # Legend in top-right corner of the top pad.
        # # Coordinates are x1, y1, x2, y2 in normalized pad coordinates.
        # leg = ROOT.TLegend(0.55, 0.55, 0.90, 0.90)
        # leg.SetBorderSize(0)
        # leg.SetFillStyle(0)
        # leg.SetTextSize(0.022)
        leg = ROOT.TLegend(.28, .35)

        leg.AddEntry(h_data, "Data", "p")

        # Red curve: raw/unprofiled null hypothesis.
        top_text = "Null Hypothesis"
        bot_text = "#chi^{{2}}={:.2f}".format(chi2_null_raw)
        leg_text = "#splitline{%s}{%s}" % (top_text, bot_text)
        leg.AddEntry(h_null, leg_text, "l")

        # Green curve: profiled null hypothesis.
        if useMarg and h_prof is not None:
            top_text = "Null Hypothesis Profiled"
            bot_text = "#chi^{{2}}={:.2f} + {:.2f} penalty".format(
                chi2_null - null_pen, null_pen
            )
            leg_text = "#splitline{%s}{%s}" % (top_text, bot_text)
            leg.AddEntry(h_prof, leg_text, "l")

        top_text = "Best Osc. Fit Raw"
        bot_text = "#chi^{{2}}={:.2f}".format(chi2_model_raw)
        leg_text = "#splitline{%s}{%s}" % (top_text, bot_text)
        leg.AddEntry(h_osc_raw, leg_text, "l")

        # Blue curve: oscillated model.
        top_text = "Best Osc. Fit"
        if useMarg:
            bot_text = "#chi^{{2}}={:.2f} + {:.2f} penalty".format(
                chi2_model - model_pen, model_pen
            )
        else:
            bot_text = "#chi^{{2}}={:.2f}".format(chi2_model_raw)

        leg_text = "#splitline{%s}{%s}" % (top_text, bot_text)
        leg.AddEntry(h_osc, leg_text, "l")

        leg.Draw()

        nullRatio = h_data.Clone()
        oscRatio = h_osc.Clone()

        nullRatio.Divide(nullRatio, h_null)
        oscRatio.Divide(oscRatio, h_null)

        oscRawRatio = h_osc_raw.Clone()
        oscRawRatio.Divide(oscRawRatio, h_null)
        oscRawRatio.SetLineColor(ROOT.kOrange + 7)
        oscRawRatio.SetLineStyle(2)
        oscRawRatio.SetLineWidth(3)


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

        nullErrors = self.ApplyExternalCovarianceErrorsForPlot(
            nullErrors,
            reference_hist=h_null,
            ratio_mode=True,
            set_content_to_one=True,
        )

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
        nullErrors.SetMinimum(.5)
        nullErrors.SetMaximum(1.5)
        nullErrors.Draw("E2")

        mask_boxes = []
        if len(self.mask_bins) > 0:
            mask_boxes = draw_masked_bin_boxes(self.mask_bins, nullErrors, 0.5, 1.5)

        nullRatio.Draw("same")
        if profRatio is not None:
            profRatio.Draw("same hist l")
        oscRawRatio.Draw("same hist l")
        oscRatio.Draw("same hist l")

        straightLine = nullErrors.Clone()
        straightLine.SetLineColor(ROOT.kRed)
        straightLine.SetLineWidth(2)
        straightLine.SetFillColor(0)
        straightLine.Draw("HIST SAME")

        top.cd()
        tag = "_{}".format(plot_tag) if plot_tag else ""
        tag += make_profile_tag(exclude)

        # Only tag lambda when using a non-nominal value.
        if abs(float(lam) - 1.0) > 1e-12:
            lam_str = "{:g}".format(float(lam))
            lam_tag = lam_str.replace(".", "p").replace("-", "m")
            tag += "_lam{}".format(lam_tag)

        if self.profile_only not in [None, "", "none"]:
            tag += "_profileOnly{}".format(str(self.profile_only).replace(",", "-"))

        if self.profile_n_universes is not None:
            tag += "_Nprof{}".format(self.profile_n_universes)

        tag += make_mask_tag(self.mask_spec)

        overall.Print(
            "plots/{}_stitched{}_{:.1f}_{:.3f}_{:.4f}.png".format(
                name, tag, parameters["m"], parameters["ue4"], parameters["umu4"]
            )
        )

        if not plotSamples:
            return

        if not useMarg:
            print("WARNING: plotSamples=True currently expects flux fitters. Skipping sample plots because useMarg=False.")
            return


        histogram = copy.deepcopy(self.histogram)

        nullSolution = statistic.GetFluxFitter(useOsc=False).GetFluxSolution()
        oscSolution = statistic.GetFluxFitter(useOsc=True).GetFluxSolution()

        mc_hists = []
        null_hists = []
        osc_hists = []
        data_hists = []
        titles = []

        if "fhc_ratio" in histogram.keys:
            desired_plots = [
                "fhc_elastic",
                "fhc_ratio",
                "rhc_ratio",
                "fhc_numu_selection",
                "rhc_numu_selection",
            ]
        else:
            desired_plots = [
                "fhc_elastic",
                "fhc_nue_selection",
                "rhc_nue_selection",
                "fhc_numu_selection",
                "rhc_numu_selection",
            ]

        plots = [
            plot for plot in desired_plots
            if plot in histogram.keys
        ]

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

            # statistic.GetFluxFitter(useOsc=False).ReweightToFluxSolution(null_hist)
            statistic.GetFluxFitter(useOsc=False).ReweightToFluxSolutionStoredA(null_hist)
            if null_hist.HasVertErrorBand("Flux"):
                null_hist.PopVertErrorBand("Flux")

            # statistic.GetFluxFitter(useOsc=True).ReweightToFluxSolution(total_hist)
            statistic.GetFluxFitter(useOsc=True).ReweightToFluxSolutionStoredA(total_hist)
            if total_hist.HasVertErrorBand("Flux"):
                total_hist.PopVertErrorBand("Flux")
            total_hist.AddMissingErrorBandsAndFillWithCV(h_data)

            TArray = []
            for hist in hists:
                if hist.Integral() <= 0:
                    continue

                # statistic.GetFluxFitter(useOsc=True).ReweightToFluxSolution(hist)
                statistic.GetFluxFitter(useOsc=True).ReweightToFluxSolutionStoredA(hist)

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
            "plots/{}_oscillated{}_{:.1f}_{:.3f}_{:.4f}.png".format(
                name, tag, parameters["m"], parameters["ue4"], parameters["umu4"]
            )
        )

# Modified from the original 2-region backgroundFit script.
# Keeps the original RunMinimizer/error-band-CV procedure.
# Adds CV/universe diagnostic plotting helpers and reachable background-subtracted plots.

import os
import sys
import ROOT
import PlotUtils
import math
import copy
from array import array
from collections import OrderedDict

from tools.PlotLibrary import HistHolder
from config.AnalysisConfig import AnalysisConfig
from config import BackgroundFitConfig
from tools import Utilities,PlotTools
from config.UnfoldingConfig import HISTOGRAMS_TO_UNFOLD
from config.DrawingConfig import SignalOnly,Default_Plot_Type,Default_Scale,DefaultPlotters,DefaultSlicer,PLOTS_TO_MAKE,SignalChargedBackground
from config.SignalDef import SIGNAL_DEFINITION
mnvplotter = PlotUtils.MnvPlotter()

from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS 
mnvplotter.error_summary_group_map.clear()
for k,v in CONSOLIDATED_ERROR_GROUPS.items():
    vec = ROOT.vector("std::string")()
    for vs in v :
        vec.push_back(vs)
    mnvplotter.error_summary_group_map[k]= vec
# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)


def SetFractionalUncertaintyYAxis(mnvplotter, ymin=0.0, ymax=1.0):
    """Force fractional-uncertainty/error-summary plots to a fixed y range."""
    mnvplotter.axis_minimum = ymin
    mnvplotter.axis_maximum = ymax


def ResetPlotterYAxis(mnvplotter):
    """Return MnvPlotter y-axis controls to automatic scaling."""
    mnvplotter.axis_minimum = -1111
    mnvplotter.axis_maximum = -1111


def DrawBandUniversesFromHist(h, bandname, out_tag):
    """
    Draw:
      - main/parent CV of h
      - band CV of the requested vertical error band
      - all universes in that band
    """

    if h is None:
        print("Histogram is None")
        return

    if not h.HasVertErrorBand(bandname):
        print(f"{h.GetName()} has no vertical band {bandname}")
        return

    band = h.GetVertErrorBand(bandname)

    c = ROOT.TCanvas(f"c_{h.GetName()}_{bandname}_{out_tag}", "", 1200, 900)

    # Main CV
    h_main = ROOT.TH1D(h)
    h_main.SetDirectory(0)
    h_main.SetName(f"{h.GetName()}_{bandname}_mainCV")
    h_main.SetLineColor(ROOT.kRed)
    h_main.SetLineWidth(4)
    h_main.SetStats(0)

    # Band CV
    h_bandcv = ROOT.TH1D(band)
    h_bandcv.SetDirectory(0)
    h_bandcv.SetName(f"{h.GetName()}_{bandname}_bandCV")
    h_bandcv.SetLineColor(ROOT.kBlue)
    h_bandcv.SetLineWidth(3)
    h_bandcv.SetLineStyle(2)
    h_bandcv.SetStats(0)

    # Universe clones
    univ_hists = []
    ymax = max(h_main.GetMaximum(), h_bandcv.GetMaximum())

    for i in range(band.GetNHists()):
        hu = ROOT.TH1D(band.GetHist(i))
        hu.SetDirectory(0)
        hu.SetName(f"{h.GetName()}_{bandname}_univ_{i}")
        hu.SetLineColor(ROOT.kGray + 1)
        hu.SetLineWidth(1)
        hu.SetStats(0)
        univ_hists.append(hu)
        ymax = max(ymax, hu.GetMaximum())

    h_main.SetMaximum(1.25 * ymax if ymax > 0 else 1.0)
    h_main.Draw("HIST")

    for hu in univ_hists:
        hu.Draw("HIST SAME")

    h_bandcv.Draw("HIST SAME")
    h_main.Draw("HIST SAME")

    leg = ROOT.TLegend(0.55, 0.70, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(h_main, "Main CV", "l")
    leg.AddEntry(h_bandcv, "Band CV", "l")
    if len(univ_hists) > 0:
        leg.AddEntry(univ_hists[0], "Universes", "l")
    leg.Draw()

    c.Modified()
    c.Update()

    # keep references alive in PyROOT
    c._h_main = h_main
    c._h_bandcv = h_bandcv
    c._univ_hists = univ_hists
    c._leg = leg

    outname = AnalysisConfig.PlotPath(h.GetName(), "BandDebug", out_tag) + ".png"
    print("Saving to:", outname)
    c.SaveAs(outname)

def DrawAllBandsFromHist(h, out_prefix="BandDebug"):
    if h is None:
        return

    for bandname in h.GetVertErrorBandNames():
        DrawBandUniversesFromHist(h, str(bandname), f"{out_prefix}_{h.GetName()}_{bandname}")

def CheckBandCVsBinByBin(hist, label="", bands=None):
    if hist is None:
        return

    if bands is None:
        bands = [str(x) for x in hist.GetVertErrorBandNames()]

    print("\n===== Band CV bin-by-bin check:", label, "=====")

    for bandname in bands:
        band = hist.GetVertErrorBand(bandname)
        if not band:
            continue

        ndiff = 0
        maxdiff = 0.0

        for b in range(0, hist.GetNbinsX() + 2):
            main = hist.GetBinContent(b)
            bcv = band.GetBinContent(b)
            diff = bcv - main
            maxdiff = max(maxdiff, abs(diff))

            if abs(diff) > 1e-10:
                ndiff += 1
                print(
                    "  {} bin {:2d}: main={:.12g} bandCV={:.12g} diff={:.12g} ratio={:.12g}".format(
                        bandname,
                        b,
                        main,
                        bcv,
                        diff,
                        bcv / main if main != 0 else 0.0,
                    )
                )

        if ndiff == 0:
            print("  {}: OK, band CV matches main CV bin-by-bin".format(bandname))
        else:
            print("  {}: {} bins differ, maxdiff={}".format(bandname, ndiff, maxdiff))

def DrawBandUniversesFromHistExplicit(h, bandname, out_tag):
    if h is None:
        return

    if not h.HasVertErrorBand(bandname):
        print("No band", bandname, "on", h.GetName())
        return

    band = h.GetVertErrorBand(bandname)

    h_main = ROOT.TH1D(h)
    h_main.SetDirectory(0)
    h_main.SetName(f"{h.GetName()}_{bandname}_mainCV")
    h_main.SetStats(0)
    h_main.SetLineColor(ROOT.kBlack)
    h_main.SetLineWidth(4)
    h_main.SetLineStyle(2)   # dashed

    h_bandcv = ROOT.TH1D(band)
    h_bandcv.SetDirectory(0)
    h_bandcv.SetName(f"{h.GetName()}_{bandname}_bandCV")
    h_bandcv.SetStats(0)
    h_bandcv.SetLineColor(ROOT.kBlue)
    h_bandcv.SetLineWidth(2)
    h_bandcv.SetMarkerColor(ROOT.kBlue)
    h_bandcv.SetMarkerStyle(20)
    h_bandcv.SetMarkerSize(1.1)

    univs = []
    ymax = max(h_main.GetMaximum(), h_bandcv.GetMaximum())

    for i in range(band.GetNHists()):
        hu = ROOT.TH1D(band.GetHist(i))
        hu.SetDirectory(0)
        hu.SetName(f"{h.GetName()}_{bandname}_univ_{i}")
        hu.SetStats(0)
        hu.SetLineColor(ROOT.kGray + 1)
        hu.SetLineWidth(1)
        univs.append(hu)
        ymax = max(ymax, hu.GetMaximum())

    c = ROOT.TCanvas(f"c_{h.GetName()}_{bandname}_{out_tag}", "", 1200, 900)

    h_main.SetMaximum(1.25 * ymax if ymax > 0 else 1.0)
    h_main.Draw("HIST")

    for hu in univs:
        hu.Draw("HIST SAME")

    # draw band CV as line+markers
    h_bandcv.Draw("HIST SAME")
    h_bandcv.Draw("P SAME")

    # redraw main CV on top so the dashed line stays visible
    h_main.Draw("HIST SAME")

    leg = ROOT.TLegend(0.55, 0.68, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(h_main, "Main CV (dashed)", "l")
    leg.AddEntry(h_bandcv, "Band CV (markers)", "lp")
    if univs:
        leg.AddEntry(univs[0], "Universes", "l")
    leg.Draw()

    c._h_main = h_main
    c._h_bandcv = h_bandcv
    c._univs = univs
    c._leg = leg

    outname = AnalysisConfig.PlotPath(h.GetName(), "BandDebug", out_tag) + ".png"
    print("Saving to:", outname)
    c.SaveAs(outname)

def DrawAllBandsFromHistExplicit(h, out_prefix):
    for bandname in h.GetVertErrorBandNames():
        DrawBandUniversesFromHistExplicit(h, str(bandname), "{}_{}".format(out_prefix, bandname))

def CloneHistHolderShallow(holder):
    clone = copy.copy(holder)
    clone.hists = {k: v.Clone(f"{v.GetName()}_clone") for k, v in holder.hists.items() if v is not None}
    return clone

def CheckTwoUniverseSymmetry(h, bandname, label="", bins=None):
    """
    For a 2-universe vertical error band, check whether the CV is centered
    between universe 0 and universe 1 bin-by-bin.
    """
    if h is None:
        return

    if not h.HasVertErrorBand(bandname):
        print(f"[CheckTwoUniverseSymmetry] {h.GetName()} has no band {bandname}")
        return

    band = h.GetVertErrorBand(bandname)

    if band.GetNHists() != 2:
        print(
            f"[CheckTwoUniverseSymmetry] {h.GetName()} {bandname} "
            f"has {band.GetNHists()} universes, not 2. Skipping."
        )
        return

    if bins is None:
        bins = range(1, h.GetNbinsX() + 1)

    up = band.GetHist(0)
    dn = band.GetHist(1)

    print("\n===== Two-universe symmetry check:", label, h.GetName(), bandname, "=====")

    for b in bins:
        cv = h.GetBinContent(b)
        u0 = up.GetBinContent(b)
        u1 = dn.GetBinContent(b)

        midpoint = 0.5 * (u0 + u1)
        halfspread = 0.5 * abs(u0 - u1)
        offset = midpoint - cv

        frac_offset = offset / cv if cv != 0 else 0.0
        frac_halfspread = halfspread / cv if cv != 0 else 0.0

        same_side = ((u0 - cv) * (u1 - cv)) > 0

        print(
            "bin {:2d}  x=[{:.3f},{:.3f}]  "
            "CV={:.6g}  u0={:.6g}  u1={:.6g}  "
            "mid={:.6g}  mid-CV={:.6g}  "
            "halfspread={:.6g}  frac_offset={:.6g}  frac_halfspread={:.6g}  same_side={}".format(
                b,
                h.GetXaxis().GetBinLowEdge(b),
                h.GetXaxis().GetBinUpEdge(b),
                cv,
                u0,
                u1,
                midpoint,
                offset,
                halfspread,
                frac_offset,
                frac_halfspread,
                same_side,
            )
        )

def CheckAllTwoUniverseBands(h, label="", bins=None):
    if h is None:
        return

    for bandname in h.GetVertErrorBandNames():
        bandname = str(bandname)
        band = h.GetVertErrorBand(bandname)
        if band and band.GetNHists() == 2:
            CheckTwoUniverseSymmetry(h, bandname, label=label, bins=bins)




def Get1DScaleFactor(variable_hists,scale_hists):
    scale_dict = {}
    comparable_scale = MakeComparableMnvHXD(variable_hists.GetHist(), scale_hists, False)
    for cate in variable_hists.hists:
        
        if variable_hists.hists[cate] is None:
            continue
        scaled =variable_hists.hists[cate].Clone()
        try:
            scale = comparable_scale[BackgroundFitConfig.CATEGORY_FACTORS[cate]]
            scaled.Multiply(scaled,scale)  
            scale_dict[BackgroundFitConfig.CATEGORY_FACTORS[cate]] = scaled.Integral()/variable_hists.hists[cate].Integral() if variable_hists.hists[cate].Integral() != 0 else 0
        except KeyError:
            pass
        del scaled
    return scale_dict

def MakeComparableMnvHXD(hist, scale_hist, y_axis=False):
    new_scale = {} 
    for cate in scale_hist: 
        new_scale[cate] = hist.Clone()
    xbins = hist.GetNbinsX()+2 #including under/overflows.
    for i in range(0,hist.GetSize()):
        nx = i%xbins
        ny = i//xbins
        scale_bin_entry = hist.GetYaxis().GetBinCenter(ny) if y_axis else hist.GetXaxis().GetBinCenter(nx)
        for cate in scale_hist:
            k = scale_hist[cate].FindBin(scale_bin_entry)
            new_scale[cate].SetBinContent(i,scale_hist[cate].GetBinContent(k))
            new_scale[cate].SetBinError(i,scale_hist[cate].GetBinError(k))
            
            for bandname in new_scale[cate].GetErrorBandNames():
                errorband = new_scale[cate].GetVertErrorBand(bandname)
                errorband.SetBinContent(i,scale_hists[cate].GetCVHistoWithStatError().GetBinContent(i))
                for ith in range(errorband.GetNHists()):
                    errorband.GetHist(ith).SetBinContent(i,scale_hist[cate].GetVertErrorBand(bandname).GetHist(ith).GetBinContent(k))
                    errorband.GetHist(ith).SetBinError(i,scale_hist[cate].GetVertErrorBand(bandname).GetHist(ith).GetBinError(k))

    return new_scale

def WriteScaleToMnvH1D(hist, scale, scale_err = None,  errorband=None,i=None):
    for group in scale:	
        if errorband is None:
            universe_hist = hist[group]
        elif i is None:
            universe_hist = hist[group].GetVertErrorBand(errorband)
        else:
            universe_hist= hist[group].GetVertErrorBand(errorband).GetHist(i)
        
        for q in range(0,universe_hist.GetNbinsX()+1):
            universe_hist.SetBinContent(q,scale[group].GetBinContent(q))
            universe_hist.SetBinError(q,scale[group].GetBinError(q))

def RunUniverseMinimizer(datasideband_histholders, datasignal_histholders, mcsideband_histholders, mcsignal_histholders, error_band = None, i = None):
    index = 0
    data_sideband = datasideband_histholders[index].GetHist()
    data_signal = datasignal_histholders[index].GetHist()
    mc_sidebandBKG = mcsideband_histholders[index].GetHist().Clone()
    mc_sidebandBKG.Reset()
    mc_signalBKG = mcsignal_histholders[index].GetHist().Clone()
    mc_signalBKG.Reset()
    mc_sidebandSIG = mcsideband_histholders[index].GetHist().Clone()
    mc_sidebandSIG.Reset()
    mc_signalSIG = mcsignal_histholders[index].GetHist().Clone()
    mc_signalSIG.Reset()

    for cate in mcsignal_histholders[index].hists:
        if cate == "Total":
            continue
        if cate not in SIGNAL_DEFINITION and cate != "NuEElastic":
            mc_sidebandBKG.Add(mcsideband_histholders[index].hists[cate])
            mc_signalBKG.Add(mcsignal_histholders[index].hists[cate])
        if cate in SIGNAL_DEFINITION:
            mc_sidebandSIG.Add(mcsideband_histholders[index].hists[cate])
            mc_signalSIG.Add(mcsignal_histholders[index].hists[cate])

    mc_signalNUEEL = mcsignal_histholders[index].hists["NuEElastic"].Clone()
    mc_sidebandNUEEL = mcsideband_histholders[index].hists["NuEElastic"].Clone()

    ### want a MC-like pseudodata signal region to avoid preliminary unblinding
    if AnalysisConfig.pseudodata:
        for q in range(0,data_signal.GetNbinsX()+1):
            data_signal.SetBinContent(q,mcsignal_histholders[index].hists["Total"].GetBinContent(q))
        for q in range(0,data_sideband.GetNbinsX()+1):
            data_sideband.SetBinContent(q,mcsideband_histholders[index].hists["Total"].GetBinContent(q))

    if error_band is not None and i is not None:
        mc_sidebandBKG = mc_sidebandBKG.GetVertErrorBand(error_band).GetHist(i).Clone()
        mc_signalBKG = mc_signalBKG.GetVertErrorBand(error_band).GetHist(i).Clone()
        mc_sidebandSIG = mc_sidebandSIG.GetVertErrorBand(error_band).GetHist(i).Clone()
        mc_signalSIG = mc_signalSIG.GetVertErrorBand(error_band).GetHist(i).Clone()
        mc_sidebandNUEEL = mc_sidebandNUEEL.GetVertErrorBand(error_band).GetHist(i).Clone()
        mc_signalNUEEL= mc_signalNUEEL.GetVertErrorBand(error_band).GetHist(i).Clone()
    elif error_band is not None:
        mc_sidebandBKG = mc_sidebandBKG.GetVertErrorBand(error_band).Clone()
        mc_signalBKG = mc_signalBKG.GetVertErrorBand(error_band).Clone()
        mc_sidebandSIG = mc_sidebandSIG.GetVertErrorBand(error_band).Clone()
        mc_signalSIG = mc_signalSIG.GetVertErrorBand(error_band).Clone()
        mc_sidebandNUEEL = mc_sidebandNUEEL.GetVertErrorBand(error_band).Clone()
        mc_signalNUEEL= mc_signalNUEEL.GetVertErrorBand(error_band).Clone()

    bkgscale = (mc_sidebandSIG * (mc_signalNUEEL - data_signal) + mc_signalSIG * (data_sideband - mc_sidebandNUEEL))/(mc_sidebandBKG * mc_signalSIG - mc_sidebandSIG * mc_signalBKG)
    sigscale = (mc_sidebandBKG * (data_signal - mc_signalNUEEL) + mc_signalBKG * (mc_sidebandNUEEL - data_sideband)) / (mc_sidebandBKG * mc_signalSIG - mc_sidebandSIG * mc_signalBKG)
    predscale = (data_sideband - mc_sidebandSIG - mc_sidebandNUEEL) / mc_sidebandBKG
    scales = {"signal":sigscale,"background":bkgscale,"prediction":predscale}

    return scales

def RunMinimizer(datasideband_histholders,datasignal_histholders, mcsideband_histholders, mcsignal_histholders,scale_hists):
    scales = RunUniverseMinimizer(datasideband_histholders,datasignal_histholders,mcsideband_histholders,mcsignal_histholders) 
    WriteScaleToMnvH1D(scale_hists,scales,None)
    hists = scale_hists
    print("Done with CV scale")

    #errorbands:
    for error_band in (mcsideband_histholders[0].GetHist().GetErrorBandNames()):
        #do errorband hist
        scales = RunUniverseMinimizer(datasideband_histholders,datasignal_histholders,mcsideband_histholders,mcsignal_histholders,error_band) 
        WriteScaleToMnvH1D(hists,scales,None,error_band)
        print("Done with error band histogram {}".format(error_band))

        for i in range(mcsideband_histholders[0].GetHist().GetVertErrorBand(error_band).GetNHists()):
            #do errorband universes 
            scales = RunUniverseMinimizer(datasideband_histholders,datasignal_histholders,mcsideband_histholders,mcsignal_histholders,error_band,i) 
            WriteScaleToMnvH1D(scale_hists,scales,None,error_band,i)

def TuneMC(hist_holder, scale_hists, x_axis=False, y_axis=False, prediction=False):
    if (x_axis and y_axis):
        return None # shouldnt happend
    elif not (x_axis or y_axis):
        try:
            ScaleCategories1D(hist_holder,scale_hists,prediction) #scale_dict is a global variable
        except AttributeError:
            return False
    else:
        comparable_scale = MakeComparableMnvHXD(hist_holder.GetHist(),scale_hists,y_axis)
        try:
            #ScaleCategories(hist_holder,comparable_scale)
            ScaleCategories(hist_holder,comparable_scale,prediction)
        except AttributeError:
            return False
    hist_holder.ResumTotal()
    return True

def ScaleCategories(hist_holder,scale_hists,prediction=False):
    for cate in hist_holder.hists:
        if cate == "Total":
            continue
        try:
            if not prediction:
                if cate not in SIGNAL_DEFINITION and cate != "NuEElastic":
                    scale = scale_hists["background"]
                    hist_holder.hists[cate].Multiply(hist_holder.hists[cate],scale)
                if cate in SIGNAL_DEFINITION:
                    scale = scale_hists["signal"]
                    hist_holder.hists[cate].Multiply(hist_holder.hists[cate],scale)
            else:
                if cate not in SIGNAL_DEFINITION and cate != "NuEElastic":
                    scale = scale_hists["background"]
                    hist_holder.hists[cate].Multiply(hist_holder.hists[cate],scale)

        except KeyError:
            print("KeyError with {} in {}".format(cate,hist_holder.sideband))
            continue

def ScaleCategories1D(hist_holder,scale_dict):
    for cate in hist_holder.hists:
        try:
            scale = scale_dict[BackgroundFitConfig.CATEGORY_FACTORS[cate]]
            hist_holder.hists[cate].Scale(scale,bin_width_normalize=True) 
        except KeyError:
            continue

def BackgroundSubtraction(data_hists, mc_hists, pred_hists, errs=None):
    if not data_hists.POT_scaled:
        data_hists.POTScale(False)
    if not mc_hists.POT_scaled:
        mc_hists.POTScale(False)
    if not pred_hists.POT_scaled:
        pred_hists.POTScale(False)
    out_data = data_hists.GetHist().Clone()
    out_mc = pred_hists.hists["Total"].Clone()
    out_data.AddMissingErrorBandsAndFillWithCV(out_mc)

    for group in mc_hists.hists:
        if group == "Total":
                continue
        elif group not in SIGNAL_DEFINITION:
            SubtractPoissonHistograms(out_data,mc_hists.hists[group]) #data tuned signal
            SubtractPoissonHistograms(out_mc,pred_hists.hists[group]) #no oscillation predicted signal

    return out_data,out_mc

def GetBackground(mc_hists):
    out_bkg = mc_hists.hists["Total"].Clone("bkgTotal")
    out_bkg.Reset()

    for group in mc_hists.hists:
        if group == "Total":
                continue
        elif group not in SIGNAL_DEFINITION:
            out_bkg.Add(mc_hists.hists[group])
    return out_bkg

def SubtractPoissonHistograms(h,h1):
    errors = []
    for i in range(h.GetSize()):
        errors.append(math.sqrt(h.GetBinError(i)**2 + h1.GetBinError(i)**2))
    h.Add(h1,-1)
    for i in range(h.GetSize()):
        h.SetBinError(i,errors[i])
    return h

def GetScaledDataMC(hist,datafile,mcfile,region):
    data_hist = HistHolder(hist,datafile,region,False,pot_scale)
    mc_hist = HistHolder(hist,mcfile,region,True,pot_scale)
    pred_hist = HistHolder(hist,mcfile,region,True,pot_scale)
    fit_on_axis = scaled_hist_name.upper() in data_hist.plot_name.upper() or "estimator" in data_hist.plot_name.lower()
    if fit_on_axis: # fit_on_axis = True
        fit_on_yaxis = ("_"+scaled_hist_name).upper() in data_hist.plot_name.upper() or "_estimator" in data_hist.plot_name.lower()
        print(("fit {} on {} axis".format(data_hist.plot_name, "y" if fit_on_yaxis else "x")))
        TuneMC(mc_hist, scale_hists, not fit_on_yaxis , fit_on_yaxis )
        TuneMC(pred_hist, scale_hists, not fit_on_yaxis , fit_on_yaxis, True)
    else:
        print(("not fitting {} on any axis".format(data_hist.plot_name)))
        variable_hist = HistHolder(BackgroundFitConfig.HIST_TO_FIT,mcfile,region,True,pot_scale) 
        scale_dict = Get1DScaleFactor(variable_hist,scale_hists)
        TuneMC(mc_hist, scale_dict, False , False )
        TuneMC(pred_hist, scale_hists, False, False, True)
    return data_hist,mc_hist,pred_hist

def MakeRatio(signalHist,sidebandHist,normsignalHist,normsidebandHist,config):
    #scale
    if "scale" in config:
        config["scale"](signalHist)
        config["scale"](sidebandHist)
        config["scale"](normsignalHist)
        config["scale"](normsidebandHist)
    else: 
        Default_Scale(signalHist)
        Default_Scale(sidebandHist)
        Default_Scale(normsignalHist)
        Default_Scale(normsidebandHist)
    sig_bkg = signalHist.hists["Total"].Clone("bkgTotal")
    sig_bkg.Reset()
    sid_bkg = sidebandHist.hists["Total"].Clone("sigTotal")
    sid_bkg.Reset()
    normsig_bkg = normsignalHist.hists["Total"].Clone("bkgTotal")
    normsig_bkg.Reset()
    normsid_bkg = normsidebandHist.hists["Total"].Clone("sigTotal")
    normsid_bkg.Reset()

    for group in signalHist.hists:
        if group == "Total":
                continue
        elif group not in SIGNAL_DEFINITION:
            sig_bkg.Add(signalHist.hists[group])
            sid_bkg.Add(sidebandHist.hists[group])
            normsig_bkg.Add(normsignalHist.hists[group])
            normsid_bkg.Add(normsidebandHist.hists[group])

    c1 = ROOT.TCanvas()
    sig_bkg.GetVertErrorBand("Flux").DrawAll("hist",True)
    c1.Print("{}_post_tuneSigBkgFlux.png".format(signalHist.plot_name))
    c1 = ROOT.TCanvas()
    sid_bkg.GetVertErrorBand("Flux").DrawAll("hist",True)
    c1.Print("{}_post_tuneSidBkgFlux.png".format(signalHist.plot_name))

    sig_errorband = sig_bkg.GetVertErrorBand("Flux")
    sid_errorband = sid_bkg.GetVertErrorBand("Flux")
    normsig_errorband = normsig_bkg.GetVertErrorBand("Flux")
    normsid_errorband = normsid_bkg.GetVertErrorBand("Flux")

    c1 = ROOT.TCanvas()
    fluxerr = sig_bkg.GetVertErrorBand("Flux").GetErrorBand(True,False).Clone()
    normfluxerr = normsig_bkg.GetVertErrorBand("Flux").GetErrorBand(True,False).Clone()
    normfluxerr.Divide(normfluxerr,fluxerr)
    normfluxerr.Draw()
    c1.Print("{}_sig_bkgPostTune_Fluxerrband.png".format(signalHist.plot_name))
    fluxerr = sid_bkg.GetVertErrorBand("Flux").GetErrorBand(True,False).Clone()
    normfluxerr = normsid_bkg.GetVertErrorBand("Flux").GetErrorBand(True,False).Clone()
    normfluxerr.Divide(normfluxerr,fluxerr)
    normfluxerr.Draw()
    c1.Print("{}_sid_bkgPostTune_Fluxerrband.png".format(signalHist.plot_name))

    sig_bkg.Scale(1/sig_bkg.Integral())
    sid_bkg.Scale(1/sid_bkg.Integral())
    sig_bkg.Divide(sig_bkg,sid_bkg)
    sig_bkg.GetXaxis().SetTitle("E_{available} + E_{lepton}")
    sig_bkg.GetYaxis().SetTitle("Ratio")
    sig_bkg.SetTitle("Signal/Sideband Background Ratio")
    #total.GetYaxis().SetTitle("Background Fraction per Bin")
    mnvplotter.DrawMCWithErrorBand(sig_bkg)

    PlotTools.Print(AnalysisConfig.PlotPath("EN4_ratio","Combined","bkgratioN4_tune"))

    c1 = ROOT.TCanvas()
    mnvplotter.DrawErrorSummary(sig_bkg,"TR",True,True,0)
    c1.Print("{}_post_tuneSigBkgErrSummary.png".format(signalHist.plot_name))

def MakePlot1(data_hists,mc_hists,config):
    if not (data_hists.valid and mc_hists.valid):
        return False
    if "scale" in config: 
        config["scale"](data_hists)
        config["scale"](mc_hists)
    else: 
        Default_Scale(data_hists)
        Default_Scale(mc_hists)

    c1 = ROOT.TCanvas()
    mnvplotter.DrawDataMCWithErrorBand(data_hists.GetHist(),mc_hists.GetHist(),1,"TR")
    c1.Print("{}_test_datamc.png".format(data_hists.sideband))
    mc_list,color,title = mc_hists.GetCateList(SignalChargedBackground)
    TArray = ROOT.TObjArray()
    for i in range(len(mc_list)):
        if color:
            mc_list[i].SetFillColor(color[i])
        if title:
            mc_list[i].SetTitle(title[i])

        TArray.Add(mc_list[i])
    c1.Clear()
    if data_hists.GetHist() and TArray:
        mnvplotter.DrawDataStackedMC(data_hists.GetHist(),TArray,1,"TR","Data",0,0,1001)
        c1.Print("{}_teststacked.png".format(data_hists.sideband))

    CanvasConfig = config.setdefault("canvasconfig",lambda x:True)
    PlotType = config.setdefault("plot_type",Default_Plot_Type)
    slicer = config.setdefault("slicer", DefaultSlicer(data_hists))
    draw_seperate_legend = config.setdefault("draw_seperate_legend",data_hists.dimension!=1 and PlotType != "migration")
    try:
        custom_tag = config["tag"]+PlotType if "tag" in config else PlotType+AnalysisConfig.bkgTune_tag
        if PlotType == "custom":
            plotfunction,hists=config["getplotters"](data_hists,mc_hists)
        else:
            if "args" in config:
                args = config["args"]
            elif "args" in DefaultPlotters[PlotType]:
                args = DefaultPlotters[PlotType]["args"]
            else:
                args = None
            if args is None:
                plotfunction,hists = DefaultPlotters[PlotType]["func"](data_hists,mc_hists)
            else:
                plotfunction,hists = DefaultPlotters[PlotType]["func"](data_hists,mc_hists,*args)
            PlotTools.MakeGridPlot(slicer,plotfunction,hists,CanvasConfig,draw_seperate_legend)
            PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name,sideband,custom_tag))
            print("plot {} made.".format(data_hists.plot_name))
    except KeyError as e:
        print("plot {} not made.".format(data_hists.plot_name))
        print(e)
        return False
    return True

def MakePlot(data_hists,mc_hists,config):
    if not (data_hists.valid and mc_hists.valid):
        return False
    #scale
    #mc_hists.ResumTotal()
    # if not mc_hists.POT_scaled and not data_hists.POT_scaled: 
    #     if "scale" in config: 
    #         config["scale"](data_hists)
    #         config["scale"](mc_hists)
    #     else: 
    #         Default_Scale(data_hists)
    #         Default_Scale(mc_hists)
    if "scale" in config:
        if not data_hists.POT_scaled:
            config["scale"](data_hists)
        if not mc_hists.POT_scaled:
            config["scale"](mc_hists)
    else:
        if not data_hists.POT_scaled:
            Default_Scale(data_hists)
        if not mc_hists.POT_scaled:
            Default_Scale(mc_hists)
    CanvasConfig = config.setdefault("canvasconfig",lambda x:True)
    PlotType = config.setdefault("plot_type",Default_Plot_Type)
    typeBool = PlotType!="migration" and PlotType!="category_hist" and PlotType!="hist2d"
    slicer = config.setdefault("slicer", DefaultSlicer(data_hists)) if typeBool else PlotTools.IdentitySlicer
    #slicer = config.setdefault("slicer", DefaultSlicer(data_hists))
    #draw_seperate_legend = config.setdefault("draw_seperate_legend",data_hists.dimension!=1 and PlotType != "migration")
    draw_seperate_legend = config.setdefault("draw_seperate_legend",data_hists.dimension!=1 and (PlotType != "migration" or PlotType != "category_hist" or PlotType != "hist2d"))
    try:
        custom_tag = config["tag"]+PlotType if "tag" in config else PlotType
        if PlotType == "custom":
            plotfunction,hists=config["getplotters"](data_hists,mc_hists)
        elif PlotType == "category_hist":
            if "args" in config:
                args = config["args"]
            elif "args" in DefaultPlotters[PlotType]:
                args = DefaultPlotters[PlotType]["args"]
            else:
                args = None
            categories = args[0]
            for category in categories:
                plotfunction,hists = DefaultPlotters[PlotType]["func"](data_hists,mc_hists,categories[category])
                PlotTools.MakeGridPlot(PlotTools.IdentitySlicer,plotfunction,hists,draw_seperate_legend=False,title=category)
                PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name,sideband,category))
                print("plot {} made for category {}.".format(data_hists.plot_name,category))
        else:
            if "args" in config:
                args = config["args"]
            elif "args" in DefaultPlotters[PlotType]:
                args = DefaultPlotters[PlotType]["args"]
            else:
                args = None
            if args is None:
                plotfunction,hists = DefaultPlotters[PlotType]["func"](data_hists,mc_hists)
            else:
                plotfunction,hists = DefaultPlotters[PlotType]["func"](data_hists,mc_hists,*args)

            if PlotType == "2Dstacked":
                PlotTools.SumGridPlots(slicer,plotfunction,hists,draw_seperate_legend=False)
            else:
                PlotTools.MakeGridPlot(slicer,plotfunction,hists,draw_seperate_legend=False)
            PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name,sideband,custom_tag))
            print("plot {} made.".format(data_hists.plot_name))
    except KeyError as e:
        print("plot {} not made.".format(data_hists.plot_name))
        print(e)
        return False
    return True

if __name__ == "__main__":
    #input knobs
    playlist=AnalysisConfig.playlist
    type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
    datafile,mcfile,pot_scale = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag)
    #output knobs:
    background_fit_tag = AnalysisConfig.bkgTune_tag
    scalefile=ROOT.TFile.Open(AnalysisConfig.BackgroundFitPath(playlist,background_fit_tag),"RECREATE")
    BackgroundFitConfig.SetGlobalParameter(background_fit_tag)

    datasideband_histholders = []
    mcsideband_histholders = []
    datasignal_histholders = []
    mcsignal_histholders = []
    scaled_hist_name = None

    #fit scale histograms
    scale_hists = {"signal":None,"background":None}
    sel_histholder = HistHolder(BackgroundFitConfig.HIST_TO_FIT,mcfile,"Signal",True,pot_scale)
    sid_histholder = HistHolder(BackgroundFitConfig.HIST_TO_FIT,mcfile,"dEdX",True,pot_scale)
    scale_hists["signal"] = sel_histholder.GetHist().Clone()
    scale_hists["signal"].Reset()
    scale_hists["signal"].GetYaxis().SetTitle("Scale Factor")
    scale_hists["signal"].SetTitle("Signal Scale Factor")
    scale_hists["background"] = sid_histholder.GetHist().Clone()
    scale_hists["background"].Reset()
    scale_hists["background"].GetYaxis().SetTitle("Scale Factor")
    scale_hists["background"].SetTitle("Background Scale Factor")
    scale_hists["prediction"] = sid_histholder.GetHist().Clone()
    scale_hists["prediction"].Reset()
    scale_hists["prediction"].GetYaxis().SetTitle("Scale Factor")
    scale_hists["prediction"].SetTitle("Background Scale Factor")
    if scaled_hist_name is None:
        scaled_hist_name = sel_histholder.plot_name

    for region in AnalysisConfig.sidebands:
        datasideband_histholders.append(HistHolder(BackgroundFitConfig.HIST_OBSERVABLE,datafile,region,False))
        mcsideband_histholders.append(HistHolder(BackgroundFitConfig.HIST_OBSERVABLE,mcfile,region,True,pot_scale))
        mcsideband_histholders[-1].POTScale(False)

    datasignal_histholders.append(HistHolder(BackgroundFitConfig.HIST_OBSERVABLE,datafile,"Signal",False))
    mcsignal_histholders.append(HistHolder(BackgroundFitConfig.HIST_OBSERVABLE,mcfile,"Signal",True,pot_scale))
    mcsignal_histholders[-1].POTScale(False)
    signalHist = HistHolder(BackgroundFitConfig.HIST_OBSERVABLE,mcfile,"Signal",True,pot_scale)
    #mcsignal_histholders[-1].POTScale(False)
    signalHist.POTScale(False)
    mc_prediction = signalHist.GetHist().Clone()
    mc_prediction.Reset()
    for cate in signalHist.hists:
        if cate in SIGNAL_DEFINITION:
            mc_prediction.Add(signalHist.hists[cate])

    RunMinimizer(datasideband_histholders,datasignal_histholders,mcsideband_histholders,mcsignal_histholders,scale_hists)

    # Optional diagnostic plots: draw the main CV, the band CV, and all universes.
    # Uncomment these when debugging a specific input histogram or the fitted scales.
    # h_signal_debug = HistHolder(BackgroundFitConfig.HIST_TO_FIT, mcfile, "Signal", True, pot_scale)
    # h_signal_debug.POTScale(False)
    # DrawAllBandsFromHistExplicit(h_signal_debug.GetHist(), "RawMC_Signal")
    # DrawAllBandsFromHistExplicit(scale_hists["signal"], "ScaleHist_signal")
    # DrawAllBandsFromHistExplicit(scale_hists["background"], "ScaleHist_background")
    # DrawAllBandsFromHistExplicit(scale_hists["prediction"], "ScaleHist_prediction")


    region = "Signal"
    for factor in scale_hists: 
        hist = scale_hists[factor]
        hist.SetXTitle("E_{estimator}")
        c1 = ROOT.TCanvas()
        #mnvplotter.axis_maximum = 3.0
        mnvplotter.DrawMCWithErrorBand(hist)
        #c1.Print("{}_scales.png".format(factor))
        PlotTools.Print(AnalysisConfig.PlotPath("EN4_scales",factor,"N4_tune"),mnvplotter,c1)
        c1 = ROOT.TCanvas()
        SetFractionalUncertaintyYAxis(mnvplotter, 0.0, 1.0)
        mnvplotter.DrawErrorSummary(hist,"TR",True,True,0)
        PlotTools.Print(AnalysisConfig.PlotPath("EN4_scale_errors",factor,"N4_tune"),mnvplotter,c1)
        ResetPlotterYAxis(mnvplotter)
        hist.Write("{}_Scale_Factor".format(factor))

    for hist in HISTOGRAMS_TO_UNFOLD:
        data_hist,mc_hist,pred_hist = GetScaledDataMC(hist,datafile,mcfile,"dEdX")
        for _h in mc_hist.hists:
            if mc_hist.hists[_h]:
                htemp = mc_hist.hists[_h]
                htemp.Write("EN4_dEdX_{}".format(_h))
 
        data_hist,mc_hist,pred_hist = GetScaledDataMC(hist,datafile,mcfile,region)
        mc_hist.GetHist().Write(data_hist.plot_name)
        subbedData, subbedMC = BackgroundSubtraction(data_hist,mc_hist,pred_hist)
        subbedData.Write(data_hist.plot_name+"_data_bkgSubbed") #added here
        mc_prediction.Write(data_hist.plot_name+"_predicted_Signal") #added here

    #change playlist
    type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
    datafile,mcfile,pot_scale = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag)

    for config in PLOTS_TO_MAKE:
        postfit_config = config.copy()
        postfit_config["tag"] = postfit_config.get("tag", "") + "postfit_"

        data_sighist, signalHist, pred_hist_sig = GetScaledDataMC(
            config["name"] if "name" in config else config,
            datafile,
            mcfile,
            "Signal"
        )
        data_sidehist, sidebandHist, pred_hist_sid = GetScaledDataMC(
            config["name"] if "name" in config else config,
            datafile,
            mcfile,
            "dEdX"
        )

        normsignalHist = HistHolder(config["name"] if "name" in config else config, mcfile, "Signal", True, pot_scale)
        normsidebandHist = HistHolder(config["name"] if "name" in config else config, mcfile, "dEdX", True, pot_scale)

        sideband_group = config.setdefault("sideband_group", ["Signal"] + AnalysisConfig.sidebands)

        if "Front dEdX" in config['name']:
            sideband = "Scaled"
            normsignalHist.Add(normsidebandHist)
            data_sighist.Add(data_sidehist)
            signalHist.Add(sidebandHist)
            pred_hist_sid.Add(pred_hist_sig)
            MakePlot(data_sighist, signalHist, postfit_config)

            if False:
                for cate in list(signalHist.hists.keys()):
                    if cate in SIGNAL_DEFINITION:
                        signalHist.hists[cate].Reset()
                    elif cate != "Total":
                        normsignalHist.hists[cate].Reset()
                signalHist.Add(normsignalHist)
                signalHist.ResumTotal()
                MakePlot(data_sighist, signalHist, postfit_config)

        if isinstance(sideband_group, list):
            for sideband in sideband_group:
                data_hist, mc_hist, pred_hist = GetScaledDataMC(
                    config["name"] if "name" in config else config,
                    datafile,
                    mcfile,
                    sideband
                )

                # 1) Fully tuned MC: signal + background scaled.
                mc_postfit_config = postfit_config.copy()
                mc_postfit_config["tag"] = postfit_config.get("tag", "") + "mcHist_"

                # 2) Prediction MC: background scaled, signal left nominal.
                pred_postfit_config = postfit_config.copy()
                pred_postfit_config["tag"] = postfit_config.get("tag", "") + "predHist_"

                # 3) Background-subtracted data vs predicted signal-like MC.
                sub_postfit_config = postfit_config.copy()
                sub_postfit_config["tag"] = postfit_config.get("tag", "") + "bkgSub_"

                if sideband == "Signal" and AnalysisConfig.pseudodata:
                    MakePlot(datasignal_histholders[0], mc_hist, mc_postfit_config)
                    MakePlot(datasignal_histholders[0], pred_hist, pred_postfit_config)
                else:
                    MakePlot(data_hist, mc_hist, mc_postfit_config)
                    MakePlot(data_hist, pred_hist, pred_postfit_config)

                if sideband == "Signal":
                    subbedData, subbedMC = BackgroundSubtraction(data_hist, mc_hist, pred_hist)

                    c = ROOT.TCanvas("c_sub", "c_sub", 1200, 1000)
                    mnvplotter.DrawDataMCWithErrorBand(subbedData, subbedMC, 1.0, "TR")
                    PlotTools.Print(
                        AnalysisConfig.PlotPath(data_hist.plot_name, sideband, sub_postfit_config["tag"]),
                        mnvplotter,
                        c
                    )

                    c2 = ROOT.TCanvas("c_sub_err", "c_sub_err", 1200, 1000)
                    SetFractionalUncertaintyYAxis(mnvplotter, 0.0, 1.0)
                    mnvplotter.DrawErrorSummary(subbedData, "TR", True, True, 0)
                    PlotTools.Print(
                        AnalysisConfig.PlotPath(data_hist.plot_name + "_err", sideband, sub_postfit_config["tag"]),
                        mnvplotter,
                        c2
                    )
                    ResetPlotterYAxis(mnvplotter)

                    # Optional stacked signal-category comparison using the nominal signal category histograms.
                    mc_list, color, title = normsignalHist.GetCateList(SignalOnly)
                    c3 = ROOT.TCanvas(f"c_signalCats_{sideband}_{data_hist.plot_name}", "c_signalCats", 1200, 1000)
                    c3.Divide(*PlotTools.CalMXN(1))
                    c3.cd(1)
                    pad = c3.GetPad(1)
                    pad.SetRightMargin(0.15)
                    pad.SetLeftMargin(.15)
                    pad.SetTopMargin(0.08)
                    pad.SetBottomMargin(0.2)

                    TArray = ROOT.TObjArray()
                    for i in range(len(mc_list)):
                        h = mc_list[i]
                        if not h:
                            continue

                        nonzero = False
                        for b in range(0, h.GetNbinsX() + 2):
                            if h.GetBinContent(b) != 0:
                                nonzero = True
                                break

                        if not nonzero:
                            print(f"Skipping empty stacked component: {h.GetName()}")
                            continue

                        if color:
                            h.SetFillColor(color[i])
                        if title:
                            h.SetTitle(title[i])

                        TArray.Add(h)

                    subbedData.GetXaxis().SetTitle("Energy_{estimator}")
                    mnvplotter.DrawDataStackedMC(subbedData, TArray, pot_scale, "TR", "Data", 0, 0, 1001)
                    PlotTools.Print(AnalysisConfig.PlotPath("data_signalCats", sideband, "N4_tune"), mnvplotter, c3)

                    c4 = ROOT.TCanvas(f"c_subbedErr_{sideband}_{data_hist.plot_name}", "c_subbedErr", 1200, 1000)
                    subbedData.SetTitle("Background Subtracted Data")
                    SetFractionalUncertaintyYAxis(mnvplotter, 0.0, 1.0)
                    mnvplotter.DrawErrorSummary(subbedData, "TR", True, True, 0)
                    PlotTools.Print(AnalysisConfig.PlotPath("data_subbedErr", sideband, "N4_tune"), mnvplotter, c4)
                    ResetPlotterYAxis(mnvplotter)

        else:
            # assuming sideband_group is a tuple of name, and list of sidebands
            sideband = sideband_group[0]
            sidebands = sideband_group[1]
            data_hist, mc_hist, pred_hist = GetScaledDataMC(
                config["name"] if "name" in config else config,
                datafile,
                mcfile,
                sidebands[0]
            )
            for _ in range(1, len(sidebands)):
                data_hist_tmp, mc_hist_tmp, pred_hist_tmp = GetScaledDataMC(
                    config["name"] if "name" in config else config,
                    datafile,
                    mcfile,
                    sidebands[_]
                )
                data_hist.Add(data_hist_tmp)
                mc_hist.Add(mc_hist_tmp)
            MakePlot(data_hist, mc_hist, postfit_config)

    #make bkg subtracted data histogram

    datafile.Close()
    mcfile.Close()
    scalefile.Close()
 

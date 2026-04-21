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

def PrintBinSummary(hist, name, bins=None):
    nb = hist.GetNbinsX()
    if bins is None:
        bins = range(1, nb + 1)

    print(f"\n--- {name} ---")
    for b in bins:
        c = hist.GetBinContent(b)
        e = hist.GetBinError(b)
        frac = e / c if c != 0 else float("inf")
        print(
            f"bin {b:2d}  "
            f"x=[{hist.GetXaxis().GetBinLowEdge(b):.3f}, {hist.GetXaxis().GetBinUpEdge(b):.3f}]  "
            f"content={c:.6g}  error={e:.6g}  frac={frac:.6g}"
        )

def UseCCnueMatrixOldStyle():
    return AnalysisConfig.bkgTune_tag in BackgroundFitConfig.CCNUE_MATRIX_OLDSTYLE_TAGS

def IsBackgroundCate(cate):
    if cate == "Total":
        return False
    if IsFixedCate(cate):
        return False
    return not IsSignalCate(cate)

def IsFixedCate(cate):
    if UseCCnueMatrixOldStyle():
        return cate in BackgroundFitConfig.CCNUE_FIXED_KEYS
    return cate == "NuEElastic"

def IsSignalCate(cate):
    if UseCCnueMatrixOldStyle():
        return cate in BackgroundFitConfig.CCNUE_SIGNAL_KEYS
    return cate in SIGNAL_DEFINITION

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
                # errorband.SetBinContent(i,scale_hists[cate].GetCVHistoWithStatError().GetBinContent(i))
                errorband.SetBinContent(i, scale_hist[cate].GetCVHistoWithStatError().GetBinContent(i))
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

        if IsBackgroundCate(cate):
            mc_sidebandBKG.Add(mcsideband_histholders[index].hists[cate])
            mc_signalBKG.Add(mcsignal_histholders[index].hists[cate])

        elif IsSignalCate(cate):
            mc_sidebandSIG.Add(mcsideband_histholders[index].hists[cate])
            mc_signalSIG.Add(mcsignal_histholders[index].hists[cate])

    mc_signalNUEEL = mcsignal_histholders[index].GetHist().Clone("mc_signal_fixed")
    mc_signalNUEEL.Reset()
    mc_sidebandNUEEL = mcsideband_histholders[index].GetHist().Clone("mc_sideband_fixed")
    mc_sidebandNUEEL.Reset()

    for cate in mcsignal_histholders[index].hists:
        if cate == "Total":
            continue
        if IsFixedCate(cate):
            mc_signalNUEEL.Add(mcsignal_histholders[index].hists[cate])
            mc_sidebandNUEEL.Add(mcsideband_histholders[index].hists[cate])

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


def RunUniverseMinimizer_RegularizedCV(
    datasideband_histholders,
    datasignal_histholders,
    mcsideband_histholders,
    mcsignal_histholders,
    lam_bkg=0.0,
    lam_sig=0.0,
):
    """
    CV-only global fit across bins with first-difference (slope) regularization:
        chi2 = chi2_data
             + lam_bkg * sum_b (B_{b+1} - B_b)^2
             + lam_sig * sum_b (S_{b+1} - S_b)^2

    Keeps the old 2-region structure:
      - sideband equation
      - signal equation
    and solves all visible bins together.
    """
    index = 0

    data_sideband = datasideband_histholders[index].GetHist()
    data_signal   = datasignal_histholders[index].GetHist()

    # Build CV component sums exactly like the old fit
    mc_sidebandBKG = mcsideband_histholders[index].GetHist().Clone("mc_sidebandBKG_cv")
    mc_sidebandBKG.Reset()
    mc_signalBKG = mcsignal_histholders[index].GetHist().Clone("mc_signalBKG_cv")
    mc_signalBKG.Reset()

    mc_sidebandSIG = mcsideband_histholders[index].GetHist().Clone("mc_sidebandSIG_cv")
    mc_sidebandSIG.Reset()
    mc_signalSIG = mcsignal_histholders[index].GetHist().Clone("mc_signalSIG_cv")
    mc_signalSIG.Reset()

    for cate in mcsignal_histholders[index].hists:
        if cate == "Total":
            continue
        if IsBackgroundCate(cate):
            mc_sidebandBKG.Add(mcsideband_histholders[index].hists[cate])
            mc_signalBKG.Add(mcsignal_histholders[index].hists[cate])
        elif IsSignalCate(cate):
            mc_sidebandSIG.Add(mcsideband_histholders[index].hists[cate])
            mc_signalSIG.Add(mcsignal_histholders[index].hists[cate])

    mc_signalFIX = mcsignal_histholders[index].GetHist().Clone("mc_signalFIX_cv")
    mc_signalFIX.Reset()
    mc_sidebandFIX = mcsideband_histholders[index].GetHist().Clone("mc_sidebandFIX_cv")
    mc_sidebandFIX.Reset()

    for cate in mcsignal_histholders[index].hists:
        if cate == "Total":
            continue
        if IsFixedCate(cate):
            mc_signalFIX.Add(mcsignal_histholders[index].hists[cate])
            mc_sidebandFIX.Add(mcsideband_histholders[index].hists[cate])

    # Pseudodata option
    if AnalysisConfig.pseudodata:
        for q in range(0, data_signal.GetNbinsX() + 1):
            data_signal.SetBinContent(q, mcsignal_histholders[index].hists["Total"].GetBinContent(q))
        for q in range(0, data_sideband.GetNbinsX() + 1):
            data_sideband.SetBinContent(q, mcsideband_histholders[index].hists["Total"].GetBinContent(q))

    nb = data_signal.GetNbinsX()   # visible bins only
    npar = 2 * nb                  # [B1..BN, S1..SN]

    H = ROOT.TMatrixD(npar, npar)
    y = ROOT.TVectorD(npar)

    def iB(b):
        return b - 1           # visible bins are 1..nb

    def iS(b):
        return nb + (b - 1)

    # -------------------------
    # Data term: two equations per bin
    # -------------------------
    for b in range(1, nb + 1):
        A_sb_bkg = float(mc_sidebandBKG.GetBinContent(b))
        A_sb_sig = float(mc_sidebandSIG.GetBinContent(b))
        rhs_sb   = float(data_sideband.GetBinContent(b) - mc_sidebandFIX.GetBinContent(b))

        A_sig_bkg = float(mc_signalBKG.GetBinContent(b))
        A_sig_sig = float(mc_signalSIG.GetBinContent(b))
        rhs_sig   = float(data_signal.GetBinContent(b) - mc_signalFIX.GetBinContent(b))

        jb = iB(b)
        js = iS(b)

        # sideband row contributes [A_sb_bkg, A_sb_sig]
        H[jb][jb] += A_sb_bkg * A_sb_bkg
        H[jb][js] += A_sb_bkg * A_sb_sig
        H[js][jb] += A_sb_bkg * A_sb_sig
        H[js][js] += A_sb_sig * A_sb_sig
        y[jb]     += A_sb_bkg * rhs_sb
        y[js]     += A_sb_sig * rhs_sb

        # signal row contributes [A_sig_bkg, A_sig_sig]
        H[jb][jb] += A_sig_bkg * A_sig_bkg
        H[jb][js] += A_sig_bkg * A_sig_sig
        H[js][jb] += A_sig_bkg * A_sig_sig
        H[js][js] += A_sig_sig * A_sig_sig
        y[jb]     += A_sig_bkg * rhs_sig
        y[js]     += A_sig_sig * rhs_sig

    # -------------------------
    # Slope penalty: (x_{b+1} - x_b)^2
    # -------------------------
    def add_slope_penalty(i1, i2, lam):
        H[i1][i1] += lam
        H[i1][i2] += -lam
        H[i2][i1] += -lam
        H[i2][i2] += lam

    for b in range(1, nb):
        if lam_bkg > 0.0:
            add_slope_penalty(iB(b), iB(b + 1), lam_bkg)
        if lam_sig > 0.0:
            add_slope_penalty(iS(b), iS(b + 1), lam_sig)

    # Solve
    svd = ROOT.TDecompSVD(H)
    xvec = ROOT.TVectorD(y)
    ok = svd.Solve(xvec)

    # Output histograms
    bkgscale = mc_sidebandBKG.Clone("bkgscale_regcv")
    sigscale = mc_signalSIG.Clone("sigscale_regcv")
    predscale = mc_sidebandBKG.Clone("predscale_regcv")

    bkgscale.Reset()
    sigscale.Reset()
    predscale.Reset()

    # Keep under/overflow at 1
    for q in range(0, nb + 2):
        bkgscale.SetBinContent(q, 1.0)
        sigscale.SetBinContent(q, 1.0)
        predscale.SetBinContent(q, 1.0)
        bkgscale.SetBinError(q, 0.0)
        sigscale.SetBinError(q, 0.0)
        predscale.SetBinError(q, 0.0)

    if ok:
        for b in range(1, nb + 1):
            bv = max(0.0, float(xvec[iB(b)]))
            sv = max(0.0, float(xvec[iS(b)]))
            bkgscale.SetBinContent(b, bv)
            sigscale.SetBinContent(b, sv)
    else:
        print("[REG] SVD solve failed in RunUniverseMinimizer_RegularizedCV, falling back to 1.0")

    # Keep prediction in the old style for now
    for b in range(1, nb + 1):
        denom = mc_sidebandBKG.GetBinContent(b)
        if denom != 0.0:
            pv = (data_sideband.GetBinContent(b)
                  - mc_sidebandSIG.GetBinContent(b)
                  - mc_sidebandFIX.GetBinContent(b)) / denom
            predscale.SetBinContent(b, pv)
        else:
            predscale.SetBinContent(b, 1.0)

    print(f"[REG-FIT] lam_bkg={lam_bkg} lam_sig={lam_sig}")
    print("[REG-FIT] signal    ", [round(sigscale.GetBinContent(b), 3) for b in range(1, nb + 1)])
    print("[REG-FIT] background", [round(bkgscale.GetBinContent(b), 3) for b in range(1, nb + 1)])

    return {"signal": sigscale, "background": bkgscale, "prediction": predscale}



# def RunMinimizer(datasideband_histholders,datasignal_histholders, mcsideband_histholders, mcsignal_histholders,scale_hists):
#     scales = RunUniverseMinimizer(datasideband_histholders,datasignal_histholders,mcsideband_histholders,mcsignal_histholders) 
#     WriteScaleToMnvH1D(scale_hists,scales,None)
#     hists = scale_hists
#     print("Done with CV scale")

#     #errorbands:
#     for error_band in (mcsideband_histholders[0].GetHist().GetErrorBandNames()):
#         #do errorband hist
#         scales = RunUniverseMinimizer(datasideband_histholders,datasignal_histholders,mcsideband_histholders,mcsignal_histholders,error_band) 
#         WriteScaleToMnvH1D(hists,scales,None,error_band)
#         print("Done with error band histogram {}".format(error_band))

#         for i in range(mcsideband_histholders[0].GetHist().GetVertErrorBand(error_band).GetNHists()):
#             #do errorband universes 
#             scales = RunUniverseMinimizer(datasideband_histholders,datasignal_histholders,mcsideband_histholders,mcsignal_histholders,error_band,i) 
#             WriteScaleToMnvH1D(scale_hists,scales,None,error_band,i)
def RunMinimizer(datasideband_histholders, datasignal_histholders,
                 mcsideband_histholders, mcsignal_histholders, scale_hists):

    reg = getattr(BackgroundFitConfig, "REGULATION_PARAMETER", 0.0)

    # CV
    if reg > 0.0:
        # regularized contents
        scales = RunUniverseMinimizer_RegularizedCV(
            datasideband_histholders,
            datasignal_histholders,
            mcsideband_histholders,
            mcsignal_histholders,
            lam_bkg=reg,
            lam_sig=reg,
        )

        # unregularized CV once, only to steal central bin errors
        raw_scales = RunUniverseMinimizer(
            datasideband_histholders,
            datasignal_histholders,
            mcsideband_histholders,
            mcsignal_histholders
        )

        print("\n===== RAW scales before copying errors =====")
        PrintBinSummary(raw_scales["background"], "RAW background scale", bins=[1,2,3,4,5])
        PrintBinSummary(raw_scales["signal"],     "RAW signal scale",     bins=[1,2,3,4,5])
        PrintBinSummary(raw_scales["prediction"], "RAW prediction scale", bins=[1,2,3,4,5])

        print("\n===== REG scales before copying errors =====")
        PrintBinSummary(scales["background"], "REG background scale", bins=[1,2,3,4,5])
        PrintBinSummary(scales["signal"],     "REG signal scale",     bins=[1,2,3,4,5])
        PrintBinSummary(scales["prediction"], "REG prediction scale", bins=[1,2,3,4,5])

        for name in ("background", "signal", "prediction"):
            h_reg = scales[name]
            h_raw = raw_scales[name]
            for b in range(0, h_reg.GetNbinsX() + 2):
                h_reg.SetBinError(b, h_raw.GetBinError(b))

        print("\n===== REG scales after copying RAW errors =====")
        PrintBinSummary(scales["background"], "REG+RAWERR background scale", bins=[1,2,3,4,5])
        PrintBinSummary(scales["signal"],     "REG+RAWERR signal scale",     bins=[1,2,3,4,5])
        PrintBinSummary(scales["prediction"], "REG+RAWERR prediction scale", bins=[1,2,3,4,5])

    else:
        scales = RunUniverseMinimizer(
            datasideband_histholders,
            datasignal_histholders,
            mcsideband_histholders,
            mcsignal_histholders
        )

        print("\n===== UNREGULARIZED CV scales =====")
        PrintBinSummary(scales["background"], "background scale", bins=[1,2,3,4,5])
        PrintBinSummary(scales["signal"],     "signal scale",     bins=[1,2,3,4,5])
        PrintBinSummary(scales["prediction"], "prediction scale", bins=[1,2,3,4,5])

    WriteScaleToMnvH1D(scale_hists, scales, None)
    print("Done with CV scale")

    print("\n===== scale_hists after WriteScaleToMnvH1D =====")
    PrintBinSummary(scale_hists["background"], "stored background scale", bins=[1,2,3,4,5])
    PrintBinSummary(scale_hists["signal"],     "stored signal scale",     bins=[1,2,3,4,5])
    PrintBinSummary(scale_hists["prediction"], "stored prediction scale", bins=[1,2,3,4,5])

    # error bands: keep original unsmoothed behavior
    for error_band in mcsideband_histholders[0].GetHist().GetErrorBandNames():
        scales = RunUniverseMinimizer(
            datasideband_histholders,
            datasignal_histholders,
            mcsideband_histholders,
            mcsignal_histholders,
            error_band
        )
        WriteScaleToMnvH1D(scale_hists, scales, None, error_band)
        print("Done with error band histogram {}".format(error_band))

        for i in range(mcsideband_histholders[0].GetHist().GetVertErrorBand(error_band).GetNHists()):
            scales = RunUniverseMinimizer(
                datasideband_histholders,
                datasignal_histholders,
                mcsideband_histholders,
                mcsignal_histholders,
                error_band,
                i
            )
            WriteScaleToMnvH1D(scale_hists, scales, None, error_band, i)    


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
                if IsBackgroundCate(cate):
                    scale = scale_hists["background"]
                    hist_holder.hists[cate].Multiply(hist_holder.hists[cate], scale)
                elif IsSignalCate(cate):
                    scale = scale_hists["signal"]
                    hist_holder.hists[cate].Multiply(hist_holder.hists[cate], scale)
            else:
                if IsBackgroundCate(cate):
                    scale = scale_hists["background"]
                    hist_holder.hists[cate].Multiply(hist_holder.hists[cate], scale)

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

def BackgroundSubtraction(data_hists, mc_hists, pred_hists, errs = None):
    data_hists.POTScale(False)
    mc_hists.POTScale(False)
    pred_hists.POTScale(False)
    out_data = data_hists.GetHist().Clone()
    out_mc = pred_hists.hists["Total"].Clone()
    out_data.AddMissingErrorBandsAndFillWithCV(out_mc)

    for group in mc_hists.hists:
        if group == "Total":
            continue
        if IsBackgroundCate(group) or IsFixedCate(group):
            SubtractPoissonHistograms(out_data, mc_hists.hists[group])
            SubtractPoissonHistograms(out_mc, pred_hists.hists[group])

    return out_data,out_mc

def GetBackground(mc_hists):
    out_bkg = mc_hists.hists["Total"].Clone("bkgTotal")
    out_bkg.Reset()

    for group in mc_hists.hists:
        if group == "Total":
            continue
        if IsBackgroundCate(group) or IsFixedCate(group):
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
    if not mc_hists.POT_scaled and not data_hists.POT_scaled: 
        if "scale" in config: 
            config["scale"](data_hists)
            config["scale"](mc_hists)
        else: 
            Default_Scale(data_hists)
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
    # input knobs
    playlist = AnalysisConfig.playlist
    type_path_map = {
        t: AnalysisConfig.SelectionHistoPath(playlist, t == "data", False)
        for t in AnalysisConfig.data_types
    }
    datafile, mcfile, pot_scale = Utilities.getFilesAndPOTScale(
        playlist, type_path_map, AnalysisConfig.ntuple_tag
    )

    # output knobs
    background_fit_tag = AnalysisConfig.bkgTune_tag
    scalefile = ROOT.TFile.Open(
        AnalysisConfig.BackgroundFitPath(playlist, background_fit_tag), "RECREATE"
    )
    BackgroundFitConfig.SetGlobalParameter(background_fit_tag)

    reg = getattr(BackgroundFitConfig, "REGULATION_PARAMETER", 0.0)
    reg_tag = f"{background_fit_tag}_lam{reg:.3f}".replace(".", "p")

    datasideband_histholders = []
    mcsideband_histholders = []
    datasignal_histholders = []
    mcsignal_histholders = []
    scaled_hist_name = None

    # fit scale histograms
    scale_hists = {"signal": None, "background": None}
    sel_histholder = HistHolder(BackgroundFitConfig.HIST_TO_FIT, mcfile, "Signal", True, pot_scale)
    sid_histholder = HistHolder(BackgroundFitConfig.HIST_TO_FIT, mcfile, "dEdX", True, pot_scale)

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
        datasideband_histholders.append(
            HistHolder(BackgroundFitConfig.HIST_OBSERVABLE, datafile, region, False)
        )
        mcsideband_histholders.append(
            HistHolder(BackgroundFitConfig.HIST_OBSERVABLE, mcfile, region, True, pot_scale)
        )
        mcsideband_histholders[-1].POTScale(False)

    datasignal_histholders.append(
        HistHolder(BackgroundFitConfig.HIST_OBSERVABLE, datafile, "Signal", False)
    )
    mcsignal_histholders.append(
        HistHolder(BackgroundFitConfig.HIST_OBSERVABLE, mcfile, "Signal", True, pot_scale)
    )
    mcsignal_histholders[-1].POTScale(False)

    signalHist = HistHolder(BackgroundFitConfig.HIST_OBSERVABLE, mcfile, "Signal", True, pot_scale)
    signalHist.POTScale(False)

    mc_prediction = signalHist.GetHist().Clone()
    mc_prediction.Reset()
    for cate in signalHist.hists:
        if IsSignalCate(cate):
            mc_prediction.Add(signalHist.hists[cate])

    RunMinimizer(
        datasideband_histholders,
        datasignal_histholders,
        mcsideband_histholders,
        mcsignal_histholders,
        scale_hists,
    )

    region = "Signal"

    # write / plot scale histograms
    for factor in scale_hists:
        hist = scale_hists[factor]
        hist.SetXTitle("E_{estimator}")

        c1 = ROOT.TCanvas()
        mnvplotter.DrawMCWithErrorBand(hist)
        PlotTools.Print(
            AnalysisConfig.PlotPath("EN4_scales", factor, reg_tag),
            mnvplotter,
            c1,
        )

        c1 = ROOT.TCanvas()
        mnvplotter.DrawErrorSummary(hist, "TR", True, True, 0)
        PlotTools.Print(
            AnalysisConfig.PlotPath("EN4_scale_errors", factor, reg_tag),
            mnvplotter,
            c1,
        )

        hist.Write("{}_Scale_Factor".format(factor))

    for hist in HISTOGRAMS_TO_UNFOLD:
        data_hist, mc_hist, pred_hist = GetScaledDataMC(hist, datafile, mcfile, "dEdX")
        for _h in mc_hist.hists:
            if mc_hist.hists[_h]:
                htemp = mc_hist.hists[_h]
                htemp.Write("EN4_dEdX_{}".format(_h))

        data_hist, mc_hist, pred_hist = GetScaledDataMC(hist, datafile, mcfile, region)
        mc_hist.GetHist().Write(data_hist.plot_name)
        subbedData, subbedMC = BackgroundSubtraction(data_hist, mc_hist, pred_hist)
        subbedData.Write(data_hist.plot_name + "_data_bkgSubbed")
        mc_prediction.Write(data_hist.plot_name + "_predicted_Signal")

    # reopen files if needed for plotting
    type_path_map = {
        t: AnalysisConfig.SelectionHistoPath(playlist, t == "data", False)
        for t in AnalysisConfig.data_types
    }
    datafile, mcfile, pot_scale = Utilities.getFilesAndPOTScale(
        playlist, type_path_map, AnalysisConfig.ntuple_tag
    )

    for config in PLOTS_TO_MAKE:
        postfit_config = config.copy()
        postfit_config["tag"] = postfit_config.get("tag", "") + "postfit_"

        data_sighist, signalHist, pred_hist_sig = GetScaledDataMC(
            config["name"] if "name" in config else config, datafile, mcfile, "Signal"
        )
        data_sidehist, sidebandHist, pred_hist_sid = GetScaledDataMC(
            config["name"] if "name" in config else config, datafile, mcfile, "dEdX"
        )
        normsignalHist = HistHolder(
            config["name"] if "name" in config else config, mcfile, "Signal", True, pot_scale
        )
        normsidebandHist = HistHolder(
            config["name"] if "name" in config else config, mcfile, "dEdX", True, pot_scale
        )

        sideband_group = config.setdefault("sideband_group", ["Signal"] + AnalysisConfig.sidebands)

        if "Front dEdX" in config["name"]:
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

        elif isinstance(sideband_group, list):
            for sideband in sideband_group:
                data_hist, mc_hist, pred_hist = GetScaledDataMC(
                    config["name"] if "name" in config else config,
                    datafile,
                    mcfile,
                    sideband,
                )
                if sideband == "Signal" and AnalysisConfig.pseudodata:
                    MakePlot(datasignal_histholders[0], pred_hist, postfit_config)
                else:
                    MakePlot(data_hist, pred_hist, postfit_config)
                    continue

                    subbedData, subbedMC = BackgroundSubtraction(data_hist, mc_hist, pred_hist)
                    mc_list, color, title = normsignalHist.GetCateList(SignalOnly)
                    c = ROOT.TCanvas("c2", "c2", 1200, 1000)
                    c.Divide(*PlotTools.CalMXN(1))
                    c.cd(1)
                    pad = c.GetPad(1)
                    pad.SetRightMargin(0.15)
                    pad.SetLeftMargin(0.15)
                    pad.SetTopMargin(0.08)
                    pad.SetBottomMargin(0.2)

                    TArray = ROOT.TObjArray()
                    for i in range(len(mc_list)):
                        if color:
                            mc_list[i].SetFillColor(color[i])
                        if title:
                            mc_list[i].SetTitle(title[i])
                        if mc_list[i]:
                            TArray.Add(mc_list[i])

                    subbedData.GetXaxis().SetTitle("Energy_{estimator}")
                    mnvplotter.DrawDataStackedMC(subbedData, TArray, pot_scale, "TR", "Data", 0, 0, 1001)
                    PlotTools.Print(
                        AnalysisConfig.PlotPath("data_signalCats", sideband, reg_tag),
                        mnvplotter,
                        c,
                    )

                    subbedData.SetTitle("Backgrounded Subtracted Data")
                    mnvplotter.DrawErrorSummary(subbedData, "TR", True, True, 0)
                    PlotTools.Print(
                        AnalysisConfig.PlotPath("data_subbedErr", sideband, reg_tag),
                        mnvplotter,
                        c,
                    )

        else:
            # assuming sideband_group is a tuple of name, and list of sidebands
            sideband = sideband_group[0]
            sidebands = sideband_group[1]

            data_hist, mc_hist = GetScaledDataMC(
                config["name"] if "name" in config else config,
                datafile,
                mcfile,
                sidebands[0],
            )
            for _ in range(1, len(sidebands)):
                data_hist_tmp, mc_hist_tmp = GetScaledDataMC(
                    config["name"] if "name" in config else config,
                    datafile,
                    mcfile,
                    sidebands[_],
                )
                data_hist.Add(data_hist_tmp)
                mc_hist.Add(mc_hist_tmp)

            MakePlot(data_hist, mc_hist, postfit_config)

    datafile.Close()
    mcfile.Close()
    scalefile.Close()
 

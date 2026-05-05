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
mnvplotter.error_summary_group_map.clear();
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

def PrintLastBinCheck(mnvhist, name):
    nb = mnvhist.GetNbinsX()

    cv = mnvhist.GetBinContent(nb)

    total_abs_hist = mnvhist.GetTotalError(False, False, False)
    total_frac_hist = mnvhist.GetTotalError(False, True, False)

    print(f"\n===== Last visible bin check: {name} =====")
    print("nbinsX =", nb)
    print("x range =", mnvhist.GetXaxis().GetBinLowEdge(nb), "to", mnvhist.GetXaxis().GetBinUpEdge(nb))
    print("CV =", cv)
    print("TOTAL abs =", total_abs_hist.GetBinContent(nb))
    print("TOTAL frac (GetTotalError) =", total_frac_hist.GetBinContent(nb))
    print("TOTAL frac (abs/CV) =", total_abs_hist.GetBinContent(nb) / cv if cv != 0 else 0.0)

    # also print overflow just in case
    of = nb + 1
    print("\nOverflow bin:")
    print("CV overflow =", mnvhist.GetBinContent(of))
    print("TOTAL abs overflow =", total_abs_hist.GetBinContent(of))
    print("TOTAL frac overflow =", total_frac_hist.GetBinContent(of))

def PrintFitInputs(data_sideband, data_signal,
                   mc_sidebandBKG, mc_sidebandSIG, mc_sidebandNUEEL,
                   mc_signalBKG, mc_signalSIG, mc_signalNUEEL,
                   bins=None):
    nb = data_signal.GetNbinsX()
    if bins is None:
        bins = range(1, nb + 1)

    print("\n--- Fit inputs ---")
    for b in bins:
        denom = (
            mc_sidebandBKG.GetBinContent(b) * mc_signalSIG.GetBinContent(b)
            - mc_sidebandSIG.GetBinContent(b) * mc_signalBKG.GetBinContent(b)
        )
        print(
            f"bin {b:2d}  "
            f"x=[{data_signal.GetXaxis().GetBinLowEdge(b):.3f}, {data_signal.GetXaxis().GetBinUpEdge(b):.3f}]  "
            f"Dsb={data_sideband.GetBinContent(b):.6g}  "
            f"Dsig={data_signal.GetBinContent(b):.6g}  "
            f"SB_BKG={mc_sidebandBKG.GetBinContent(b):.6g}  "
            f"SB_SIG={mc_sidebandSIG.GetBinContent(b):.6g}  "
            f"SB_NuEEl={mc_sidebandNUEEL.GetBinContent(b):.6g}  "
            f"SIG_BKG={mc_signalBKG.GetBinContent(b):.6g}  "
            f"SIG_SIG={mc_signalSIG.GetBinContent(b):.6g}  "
            f"SIG_NuEEl={mc_signalNUEEL.GetBinContent(b):.6g}  "
            f"denom={denom:.6g}"
        )

def PrintScaledContributions(mc_sidebandBKG, mc_sidebandSIG, mc_sidebandNUEEL,
                             mc_signalBKG, mc_signalSIG, mc_signalNUEEL,
                             bkgscale, sigscale, bins=None):
    nb = bkgscale.GetNbinsX()
    if bins is None:
        bins = range(1, nb + 1)

    print("\n--- Postfit contributions ---")
    for b in bins:
        sb_bkg = mc_sidebandBKG.GetBinContent(b) * bkgscale.GetBinContent(b)
        sb_sig = mc_sidebandSIG.GetBinContent(b) * sigscale.GetBinContent(b)
        sb_fix = mc_sidebandNUEEL.GetBinContent(b)

        sig_bkg = mc_signalBKG.GetBinContent(b) * bkgscale.GetBinContent(b)
        sig_sig = mc_signalSIG.GetBinContent(b) * sigscale.GetBinContent(b)
        sig_fix = mc_signalNUEEL.GetBinContent(b)

        print(
            f"bin {b:2d}  "
            f"SB: bkg={sb_bkg:.3f} sig={sb_sig:.3f} fix={sb_fix:.3f} total={sb_bkg+sb_sig+sb_fix:.3f}  "
            f"SIG: bkg={sig_bkg:.3f} sig={sig_sig:.3f} fix={sig_fix:.3f} total={sig_bkg+sig_sig+sig_fix:.3f}"
        )

def PrintErrorBreakdownConsistent(mnvhist, name, bins=None):
    import math

    nb = mnvhist.GetNbinsX()
    if bins is None:
        bins = range(1, nb + 1)

    total_cov = mnvhist.GetTotalErrorMatrix(False, False, False)
    stat_cov  = mnvhist.GetStatErrorMatrix()

    print(f"\n===== Consistent error breakdown: {name} =====")
    for b in bins:
        xlo = mnvhist.GetXaxis().GetBinLowEdge(b)
        xhi = mnvhist.GetXaxis().GetBinUpEdge(b)
        cv  = mnvhist.GetBinContent(b)

        # use ROOT bin index directly
        ib = b

        stat2  = stat_cov[ib][ib]
        total2 = total_cov[ib][ib]

        stat  = math.sqrt(stat2) if stat2 > 0 else 0.0
        total = math.sqrt(total2) if total2 > 0 else 0.0

        stat_frac  = stat / cv if cv != 0 else float("inf")
        total_frac = total / cv if cv != 0 else float("inf")

        print(
            f"\nbin {b:2d}  x=[{xlo:.3f}, {xhi:.3f}]  "
            f"CV={cv:.6g}  STAT={stat:.6g} ({stat_frac:.6g})  "
            f"TOTAL_SYST={total:.6g} ({total_frac:.6g})"
        )

        quad_sum = 0.0
        for bandname in mnvhist.GetErrorBandNames():
            label = str(bandname)
            cov = mnvhist.GetVertErrorBand(bandname).CalcCovMx(False, False)
            err2 = cov[ib][ib]
            err = math.sqrt(err2) if err2 > 0 else 0.0
            frac = err / cv if cv != 0 else float("inf")
            quad_sum += err2
            print(f"   {label:25} : {err:.6g} ({frac:.6g})")

        quad = math.sqrt(quad_sum) if quad_sum > 0 else 0.0
        quad_frac = quad / cv if cv != 0 else float("inf")
        print(f"   {'QUAD_SUM_OF_BANDS':25} : {quad:.6g} ({quad_frac:.6g})")

def PrintGroupedFracForPlot(mnvhist, mnvplotter, name, bins=None):
    import math

    nb = mnvhist.GetNbinsX()
    if bins is None:
        bins = range(1, nb + 1)

    available = set(str(x) for x in mnvhist.GetErrorBandNames())

    print(f"\n===== Grouped fractional uncertainties for plot: {name} =====")
    print("Available bands on hist:")
    print(sorted(available))

    for b in bins:
        cv = mnvhist.GetBinContent(b)
        xlo = mnvhist.GetXaxis().GetBinLowEdge(b)
        xhi = mnvhist.GetXaxis().GetBinUpEdge(b)

        print(f"\nbin {b:2d} x=[{xlo:.3f}, {xhi:.3f}] CV={cv:.6g}")

        for item in mnvplotter.error_summary_group_map:
            groupname = str(item.first)
            names = item.second

            err2 = 0.0
            used = []
            skipped = []

            for i in range(names.size()):
                bandname = str(names[i])

                if bandname not in available:
                    skipped.append(bandname)
                    continue

                band = mnvhist.GetVertErrorBand(bandname)
                if not band:
                    skipped.append(bandname)
                    continue

                cov = band.CalcCovMx(False, False)
                err2 += cov[b][b]
                used.append(bandname)

            err = math.sqrt(err2) if err2 > 0 else 0.0
            frac = err / cv if cv != 0 else float("inf")

            print(f"   {groupname:25} : abs={err:.6g}  frac={frac:.6g}")
            if skipped:
                print(f"      skipped: {skipped}")

def BuildGroupedFracHist(mnvhist, group_map, group_name):
    import math

    h = mnvhist.GetCVHistoWithStatError().Clone(f"{mnvhist.GetName()}_{group_name}_frac")
    h.Reset()

    nb = mnvhist.GetNbinsX()
    for b in range(1, nb + 1):
        cv = mnvhist.GetBinContent(b)
        err2 = 0.0

        for bandname in group_map[group_name]:
            if bandname not in [str(x) for x in mnvhist.GetErrorBandNames()]:
                continue
            band = mnvhist.GetVertErrorBand(bandname)
            if not band:
                continue
            cov = band.CalcCovMx(False, False)
            err2 += cov[b][b]

        err = math.sqrt(err2) if err2 > 0 else 0.0
        frac = err / cv if cv != 0 else 0.0
        h.SetBinContent(b, frac)
        h.SetBinError(b, 0.0)

    h.GetYaxis().SetTitle("Fractional uncertainty")
    return h

group_map = {
    "Flux": ["Flux"],
    "Electron Reconstruction": [
        "eltheta", "elE_ECAL", "elE_HCAL", "elE_Tracker", "electron_scale"
    ],
    "MnvTunes": [
        "RPA_HighQ2", "RPA_LowQ2", "Low_Recoil_2p2h_Tune",
        "LowQ2Pi", "fsi_weight", "SuSA_Valencia_Weight", "MK_model"
    ],
    "Interaction model": [
        "GENIE_AGKYxF1pi", "GENIE_AhtBY", "GENIE_BhtBY",
        "GENIE_CCQEPauliSupViaKF", "GENIE_CV1uBY", "GENIE_CV2uBY",
        "GENIE_D2_MaRES", "GENIE_D2_NormCCRES", "GENIE_EP_MvRES",
        "GENIE_EtaNCEL", "GENIE_FrAbs_N", "GENIE_FrAbs_pi",
        "GENIE_FrCEx_N", "GENIE_FrCEx_pi", "GENIE_FrElas_N",
        "GENIE_FrElas_pi", "GENIE_FrInel_N", "GENIE_FrPiProd_N",
        "GENIE_FrPiProd_pi", "GENIE_MFP_N", "GENIE_MFP_pi",
        "GENIE_MaNCEL", "GENIE_MaZExpCCQE", "GENIE_NormDISCC",
        "GENIE_NormNCRES", "GENIE_RDecBR1gamma", "GENIE_Rvn1pi",
        "GENIE_Rvn2pi", "GENIE_Rvp1pi", "GENIE_Rvp2pi",
        "GENIE_Theta_Delta2Npi", "GENIE_VecFFCCQEshape"
    ],
    "Detector model": [
        "beam_angle", "Leakage_Uncertainty", "Target_Mass_CH",
        "response_p", "response_meson", "response_em",
        "response_other", "response_xtalk"
    ],
}

def PrintCustomGroupedFrac(mnvhist, group_map, bins=None):
    import math

    if bins is None:
        bins = range(1, mnvhist.GetNbinsX() + 1)

    available = set(str(x) for x in mnvhist.GetErrorBandNames())

    print("\n===== Custom grouped fractional uncertainties =====")
    for b in bins:
        cv = mnvhist.GetBinContent(b)
        print(f"\nbin {b}  CV={cv:.6g}")
        for g, bands in group_map.items():
            err2 = 0.0
            used = []
            for bandname in bands:
                if bandname not in available:
                    continue
                band = mnvhist.GetVertErrorBand(bandname)
                if not band:
                    continue
                cov = band.CalcCovMx(False, False)
                err2 += cov[b][b]
                used.append(bandname)
            err = math.sqrt(err2) if err2 > 0 else 0.0
            frac = err / cv if cv != 0 else 0.0
            print(f"   {g:25} : abs={err:.6g}  frac={frac:.6g}")

def PrintAllBandErrorsForBin(mnvhist, name, b=None):
    if b is None:
        b = mnvhist.GetNbinsX()

    cv = mnvhist.GetBinContent(b)
    xlo = mnvhist.GetXaxis().GetBinLowEdge(b)
    xhi = mnvhist.GetXaxis().GetBinUpEdge(b)

    total_abs_hist = mnvhist.GetTotalError(False, False, False)
    total_frac_hist = mnvhist.GetTotalError(False, True, False)
    stat_abs_hist = mnvhist.GetStatError(False)
    stat_frac_hist = mnvhist.GetStatError(True)

    print(f"\n===== Full raw-band error dump: {name} =====")
    print(f"bin {b} x=[{xlo:.3f}, {xhi:.3f}]")
    print(f"CV = {cv:.12g}")
    print(f"STAT  abs = {stat_abs_hist.GetBinContent(b):.12g}")
    print(f"STAT  frac = {stat_frac_hist.GetBinContent(b):.12g}")
    print(f"TOTAL abs = {total_abs_hist.GetBinContent(b):.12g}")
    print(f"TOTAL frac = {total_frac_hist.GetBinContent(b):.12g}")
    print(f"TOTAL abs/CV = {(total_abs_hist.GetBinContent(b)/cv if cv != 0 else 0.0):.12g}")

    rows = []
    for bandname in mnvhist.GetErrorBandNames():
        label = str(bandname)
        band = mnvhist.GetVertErrorBand(bandname)
        if not band:
            continue

        abs_hist = band.GetErrorBand(False, False)
        frac_hist = band.GetErrorBand(True, False)

        abs_val = abs_hist.GetBinContent(b)
        frac_val = frac_hist.GetBinContent(b)
        frac_from_abs = abs_val / cv if cv != 0 else 0.0

        rows.append((max(frac_val, frac_from_abs), label, abs_val, frac_val, frac_from_abs))

    rows.sort(reverse=True)

    print("\n--- raw bands sorted by max(frac_plot, frac_abs_over_cv) ---")
    for _, label, abs_val, frac_val, frac_from_abs in rows:
        print(
            f"{label:25}  "
            f"abs={abs_val:.12g}  "
            f"frac_plot={frac_val:.12g}  "
            f"frac_abs_over_cv={frac_from_abs:.12g}"
        )

def PrintAllGroupedErrorsForBin(mnvhist, group_map, name, b=None):
    import math

    if b is None:
        b = mnvhist.GetNbinsX()

    cv = mnvhist.GetBinContent(b)
    xlo = mnvhist.GetXaxis().GetBinLowEdge(b)
    xhi = mnvhist.GetXaxis().GetBinUpEdge(b)

    available = set(str(x) for x in mnvhist.GetErrorBandNames())

    print(f"\n===== Full grouped error dump: {name} =====")
    print(f"bin {b} x=[{xlo:.3f}, {xhi:.3f}]")
    print(f"CV = {cv:.12g}")

    quad_abs = 0.0
    for g, bands in group_map.items():
        err2 = 0.0
        used = []
        for bandname in bands:
            if bandname not in available:
                continue
            band = mnvhist.GetVertErrorBand(bandname)
            if not band:
                continue
            abs_hist = band.GetErrorBand(False, False)
            err = abs_hist.GetBinContent(b)
            err2 += err * err
            used.append(bandname)

        abs_val = math.sqrt(err2) if err2 > 0 else 0.0
        frac_val = abs_val / cv if cv != 0 else 0.0
        quad_abs += err2

        print(
            f"{g:25}  abs={abs_val:.12g}  frac={frac_val:.12g}  used={used}"
        )

    quad_abs = math.sqrt(quad_abs) if quad_abs > 0 else 0.0
    quad_frac = quad_abs / cv if cv != 0 else 0.0
    print(f"\nGROUP quad sum abs  = {quad_abs:.12g}")
    print(f"GROUP quad sum frac = {quad_frac:.12g}")

def PrintSignalCategoryYieldChanges(prefit_holder, postfit_holder, bins=None):
    if bins is None:
        bins = range(1, prefit_holder.GetHist().GetNbinsX() + 1)

    print("\n===== Signal-category yield changes =====")
    for cate in prefit_holder.hists:
        if cate == "Total":
            continue
        if cate not in SIGNAL_DEFINITION:
            continue

        print(f"\n--- {cate} ---")
        for b in bins:
            pre = prefit_holder.hists[cate].GetBinContent(b)
            post = postfit_holder.hists[cate].GetBinContent(b)
            ratio = (post / pre) if pre != 0 else 0.0
            print(f"bin {b:2d}: pre={pre:.6g}  post={post:.6g}  post/pre={ratio:.6g}")

def PrintSummedSignalYieldChanges(prefit_holder, postfit_holder, bins=None):
    if bins is None:
        bins = range(1, prefit_holder.GetHist().GetNbinsX() + 1)

    print("\n===== Summed signal-component changes =====")
    for b in bins:
        pre = 0.0
        post = 0.0
        for cate in prefit_holder.hists:
            if cate == "Total":
                continue
            if cate in SIGNAL_DEFINITION:
                pre += prefit_holder.hists[cate].GetBinContent(b)
                post += postfit_holder.hists[cate].GetBinContent(b)

        ratio = (post / pre) if pre != 0 else 0.0
        print(f"bin {b:2d}: pre={pre:.6g}  post={post:.6g}  post/pre={ratio:.6g}")

def CloneHistHolderShallow(holder):
    clone = copy.copy(holder)
    clone.hists = {k: v.Clone(f"{v.GetName()}_clone") for k, v in holder.hists.items() if v is not None}
    return clone

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
                errorband.SetBinContent(i, scale_hist[cate].GetCVHistoWithStatError().GetBinContent(k))
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

    # if error_band is None and i is None:
    #     PrintFitInputs(
    #         data_sideband, data_signal,
    #         mc_sidebandBKG, mc_sidebandSIG, mc_sidebandNUEEL,
    #         mc_signalBKG, mc_signalSIG, mc_signalNUEEL,
    #         bins=[1,2,3,4,5]
    #     )

    bkgscale = (mc_sidebandSIG * (mc_signalNUEEL - data_signal) + mc_signalSIG * (data_sideband - mc_sidebandNUEEL))/(mc_sidebandBKG * mc_signalSIG - mc_sidebandSIG * mc_signalBKG)
    sigscale = (mc_sidebandBKG * (data_signal - mc_signalNUEEL) + mc_signalBKG * (mc_sidebandNUEEL - data_sideband)) / (mc_sidebandBKG * mc_signalSIG - mc_sidebandSIG * mc_signalBKG)
    predscale = (data_sideband - mc_sidebandSIG - mc_sidebandNUEEL) / mc_sidebandBKG
    scales = {"signal":sigscale,"background":bkgscale,"prediction":predscale}

    # if error_band is None and i is None:
    #     PrintScaledContributions(
    #         mc_sidebandBKG, mc_sidebandSIG, mc_sidebandNUEEL,
    #         mc_signalBKG, mc_signalSIG, mc_signalNUEEL,
    #         bkgscale, sigscale,
    #         bins=[1,2,3,4,5]
    #     )
    #     print("\n===== UNREGULARIZED CV scales =====")
    #     PrintBinSummary(scales["background"], "background scale", bins=[1,2,3,4,5])
    #     PrintBinSummary(scales["signal"],     "signal scale",     bins=[1,2,3,4,5])
    #     PrintBinSummary(scales["prediction"], "prediction scale", bins=[1,2,3,4,5])

    return scales

def RunMinimizer(datasideband_histholders,datasignal_histholders, mcsideband_histholders, mcsignal_histholders,scale_hists):
    scales = RunUniverseMinimizer(
        datasideband_histholders,
        datasignal_histholders,
        mcsideband_histholders,
        mcsignal_histholders
    )
    WriteScaleToMnvH1D(scale_hists,scales,None)
    hists = scale_hists
    print("Done with CV scale")

    # print("\n===== scale_hists after WriteScaleToMnvH1D =====")
    # PrintBinSummary(scale_hists["background"], "stored background scale", bins=[1,2,3,4,5])
    # PrintBinSummary(scale_hists["signal"],     "stored signal scale",     bins=[1,2,3,4,5])
    # PrintBinSummary(scale_hists["prediction"], "stored prediction scale", bins=[1,2,3,4,5])

    # errorbands:
    for error_band in mcsideband_histholders[0].GetHist().GetErrorBandNames():
        scales = RunUniverseMinimizer(
            datasideband_histholders,
            datasignal_histholders,
            mcsideband_histholders,
            mcsignal_histholders,
            error_band
        )
        WriteScaleToMnvH1D(hists,scales,None,error_band)
        print("Done with error band histogram {}".format(error_band))

        for i in range(mcsideband_histholders[0].GetHist().GetVertErrorBand(error_band).GetNHists()):
            scales = RunUniverseMinimizer(
                datasideband_histholders,
                datasignal_histholders,
                mcsideband_histholders,
                mcsignal_histholders,
                error_band,i
            )
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

    # save a prefit copy of MC before tuning
    mc_hist_prefit = CloneHistHolderShallow(mc_hist)

    fit_on_axis = scaled_hist_name.upper() in data_hist.plot_name.upper() or "estimator" in data_hist.plot_name.lower()
    if fit_on_axis:
        fit_on_yaxis = ("_"+scaled_hist_name).upper() in data_hist.plot_name.upper() or "_estimator" in data_hist.plot_name.lower()
        TuneMC(mc_hist, scale_hists, not fit_on_yaxis , fit_on_yaxis )
        TuneMC(pred_hist, scale_hists, not fit_on_yaxis , fit_on_yaxis, True)
    else:
        variable_hist = HistHolder(BackgroundFitConfig.HIST_TO_FIT,mcfile,region,True,pot_scale)
        scale_dict = Get1DScaleFactor(variable_hist,scale_hists)
        TuneMC(mc_hist, scale_dict, False , False )
        TuneMC(pred_hist, scale_hists, False, False, True)
        del scale_dict
        del variable_hist

    # # compare before/after tuning here
    # if region == "Signal":
    #     PrintSummedSignalYieldChanges(mc_hist_prefit, mc_hist)
    #     PrintSignalCategoryYieldChanges(mc_hist_prefit, mc_hist)

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

    # PrintErrorBreakdownConsistent(scale_hists["signal"], "signal scale", bins=[1,2,3,4,5])
    # PrintErrorBreakdownConsistent(scale_hists["background"], "background scale", bins=[1,2,3,4,5])
    # PrintErrorBreakdownConsistent(scale_hists["prediction"], "prediction scale", bins=[1,2,3,4,5])

    # PrintGroupedFracForPlot(scale_hists["signal"], mnvplotter, "signal scale", bins=[1,2,3,4,5])

    # PrintLastBinCheck(scale_hists["background"], "background scale")
    # PrintAllBandErrorsForBin(scale_hists["background"], "background scale", b=12)
    # PrintAllGroupedErrorsForBin(scale_hists["background"], group_map, "background scale", b=12)


    # flux_band = scale_hists["signal"].GetVertErrorBand("Flux")
    # h_flux_frac = flux_band.GetErrorBand(True, False)
    # h_flux_abs  = flux_band.GetErrorBand(False, False)

    # print("Flux GetErrorBand frac bin1 =", h_flux_frac.GetBinContent(1))
    # print("Flux GetErrorBand abs  bin1 =", h_flux_abs.GetBinContent(1))

    # tot_frac = scale_hists["signal"].GetTotalError(False, True, False)
    # tot_abs  = scale_hists["signal"].GetTotalError(False, False, False)

    # print("Total GetTotalError frac bin1 =", tot_frac.GetBinContent(1))
    # print("Total GetTotalError abs  bin1 =", tot_abs.GetBinContent(1))

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
        #mnvplotter.axis_maximum = 1.0

        # mnvplotter.error_summary_group_map.clear()

        # def vec(lst):
        #     v = ROOT.vector("std::string")()
        #     for x in lst:
        #         v.push_back(x)
        #     return v

        # mnvplotter.error_summary_group_map["Flux"] = vec(["Flux"])
        # mnvplotter.error_summary_group_map["Electron Reconstruction"] = vec([
        #     "eltheta", "elE_ECAL", "elE_HCAL", "elE_Tracker", "electron_scale"
        # ])
        # mnvplotter.error_summary_group_map["MnvTunes"] = vec([
        #     "RPA_HighQ2", "RPA_LowQ2", "Low_Recoil_2p2h_Tune",
        #     "LowQ2Pi", "fsi_weight", "SuSA_Valencia_Weight", "MK_model"
        # ])
        # mnvplotter.error_summary_group_map["Interaction model"] = vec([
        #     "GENIE_AGKYxF1pi", "GENIE_AhtBY", "GENIE_BhtBY",
        #     "GENIE_CCQEPauliSupViaKF", "GENIE_CV1uBY", "GENIE_CV2uBY",
        #     "GENIE_D2_MaRES", "GENIE_D2_NormCCRES", "GENIE_EP_MvRES",
        #     "GENIE_EtaNCEL", "GENIE_FrAbs_N", "GENIE_FrAbs_pi",
        #     "GENIE_FrCEx_N", "GENIE_FrCEx_pi", "GENIE_FrElas_N",
        #     "GENIE_FrElas_pi", "GENIE_FrInel_N", "GENIE_FrPiProd_N",
        #     "GENIE_FrPiProd_pi", "GENIE_MFP_N", "GENIE_MFP_pi",
        #     "GENIE_MaNCEL", "GENIE_MaZExpCCQE", "GENIE_NormDISCC",
        #     "GENIE_NormNCRES", "GENIE_RDecBR1gamma", "GENIE_Rvn1pi",
        #     "GENIE_Rvn2pi", "GENIE_Rvp1pi", "GENIE_Rvp2pi",
        #     "GENIE_Theta_Delta2Npi", "GENIE_VecFFCCQEshape"
        # ])
        # mnvplotter.error_summary_group_map["Detector model"] = vec([
        #     "beam_angle", "Leakage_Uncertainty", "Target_Mass_CH",
        #     "response_p", "response_meson", "response_em",
        #     "response_other", "response_xtalk"
        # ])

        # mnvplotter.DrawErrorSummary(hist,"TR",True,True,0)
        mnvplotter.DrawErrorSummaryDerived(hist, "TR", True, True, 0, False, "", True, "", False, "HIST")

        PlotTools.Print(AnalysisConfig.PlotPath("EN4_scale_errors",factor,"N4_tune"),mnvplotter,c1)
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

        b = 1  # or whichever bin looks suspicious

        h = subbedData
        cv = h.GetBinContent(b)

        flux_band = h.GetVertErrorBand("Flux")
        h_flux_abs = flux_band.GetErrorBand(False, False)
        h_flux_frac = flux_band.GetErrorBand(True, False)

        print("CV =", cv)
        print("Flux abs =", h_flux_abs.GetBinContent(b))
        print("Flux frac plot =", h_flux_frac.GetBinContent(b))
        print("Flux abs/CV =", h_flux_abs.GetBinContent(b) / cv if cv != 0 else 0.0)

        subbedData.Write(data_hist.plot_name+"_data_bkgSubbed") #added here
        mc_prediction.Write(data_hist.plot_name+"_predicted_Signal") #added here

    #change playlist
    type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
    datafile,mcfile,pot_scale = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag)

    for config in PLOTS_TO_MAKE:
        postfit_config = config.copy()
        postfit_config["tag"] = postfit_config.get("tag", "") + "postfit_"
        data_sighist,signalHist,pred_hist_sig = GetScaledDataMC(config["name"] if "name" in config else config,datafile,mcfile,"Signal")
        data_sidehist,sidebandHist,pred_hist_sid = GetScaledDataMC(config["name"] if "name" in config else config,datafile,mcfile,"dEdX")
        normsignalHist = HistHolder(config["name"] if "name" in config else config,mcfile,"Signal",True,pot_scale)
        normsidebandHist = HistHolder(config["name"] if "name" in config else config,mcfile,"dEdX",True,pot_scale)
        #MakeRatio(signalHist,sidebandHist,normsignalHist,normsidebandHist,config)
        sideband_group =  config.setdefault("sideband_group",["Signal"]+AnalysisConfig.sidebands)
        if "Front dEdX" in config['name']:
            sideband = "Scaled"
            normsignalHist.Add(normsidebandHist)
            data_sighist.Add(data_sidehist)
            signalHist.Add(sidebandHist)
            pred_hist_sid.Add(pred_hist_sig)
            MakePlot(data_sighist,signalHist,postfit_config)

            if False:
                for cate in list(signalHist.hists.keys()):
                    if cate in SIGNAL_DEFINITION:
                        signalHist.hists[cate].Reset()
                    elif cate != "Total":
                        normsignalHist.hists[cate].Reset()
                signalHist.Add(normsignalHist)
                signalHist.ResumTotal()
                MakePlot(data_sighist,signalHist,postfit_config)

        # elif isinstance(sideband_group,list):
        #     for sideband in sideband_group:
        #         data_hist,mc_hist,pred_hist = GetScaledDataMC(config["name"] if "name" in config else config,datafile,mcfile,sideband)
        #         if sideband == "Signal" and AnalysisConfig.pseudodata:
        #             MakePlot(datasignal_histholders[0],pred_hist,postfit_config)
        #         else:
        #             MakePlot(data_hist,pred_hist,postfit_config)
        #             #MakePlot(data_hist,normsignalHist,config) 
        #             continue

        if isinstance(sideband_group,list):
            for sideband in sideband_group:
                data_hist, mc_hist, pred_hist = GetScaledDataMC(
                    config["name"] if "name" in config else config,
                    datafile, mcfile, sideband
                )

                # 1) Fully tuned MC: signal + background scaled
                mc_postfit_config = postfit_config.copy()
                mc_postfit_config["tag"] = postfit_config.get("tag", "") + "mcHist_"

                # 2) Prediction MC: background only scaled
                pred_postfit_config = postfit_config.copy()
                pred_postfit_config["tag"] = postfit_config.get("tag", "") + "predHist_"

                # 3) Background-subtracted data vs predicted signal-like MC
                sub_postfit_config = postfit_config.copy()
                sub_postfit_config["tag"] = postfit_config.get("tag", "") + "bkgSub_"

                if sideband == "Signal" and AnalysisConfig.pseudodata:
                    MakePlot(datasignal_histholders[0], mc_hist, mc_postfit_config)
                    MakePlot(datasignal_histholders[0], pred_hist, pred_postfit_config)
                else:
                    MakePlot(data_hist, mc_hist, mc_postfit_config)
                    MakePlot(data_hist, pred_hist, pred_postfit_config)

                # Background-subtracted data and predicted signal
                if sideband == "Signal":
                    subbedData, subbedMC = BackgroundSubtraction(data_hist, mc_hist, pred_hist)

                    # Wrap these in HistHolder-like plotting only if your MakePlot expects HistHolder.
                    # Otherwise do a direct plotting call here.
                    c = ROOT.TCanvas("c_sub", "c_sub", 1200, 1000)
                    mnvplotter.DrawDataMCWithErrorBand(subbedData, subbedMC, 1.0, "TR")
                    PlotTools.Print(
                        AnalysisConfig.PlotPath(data_hist.plot_name, sideband, sub_postfit_config["tag"]),
                        mnvplotter, c
                    )

                    c2 = ROOT.TCanvas("c_sub_err", "c_sub_err", 1200, 1000)
                    # mnvplotter.DrawErrorSummary(subbedData, "TR", True, True, 0)
                    mnvplotter.DrawErrorSummaryDerived(subbedData, "TR", True, True, 0, False, "", True, "", False, "HIST")
                    PlotTools.Print(
                        AnalysisConfig.PlotPath(data_hist.plot_name + "_err", sideband, sub_postfit_config["tag"]),
                        mnvplotter, c2
                    )

                    subbedData, subbedMC = BackgroundSubtraction(data_hist,mc_hist,pred_hist)
                    mc_list,color,title = normsignalHist.GetCateList(SignalOnly)
                    # c = ROOT.TCanvas("c2","c2",1200,1000)
                    c = ROOT.TCanvas(f"c2_{sideband}_{data_hist.plot_name}", "c2", 1200, 1000)
                    c.Divide(*PlotTools.CalMXN(1))
                    c.cd(1)
                    pad = c.GetPad(1)
                    pad.SetRightMargin(0.15)
                    pad.SetLeftMargin(.15)
                    pad.SetTopMargin(0.08)
                    pad.SetBottomMargin(0.2)
                    TArray = ROOT.TObjArray()
                    # for i in range(len(mc_list)):
                    #     if color:
                    #         mc_list[i].SetFillColor(color[i])
                    #     if title:
                    #         mc_list[i].SetTitle(title[i])

                    #     if mc_list[i]:
                    #         TArray.Add(mc_list[i])
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
                    #mnvplotter.axis_maximum = 0.4
                    mnvplotter.DrawDataStackedMC(subbedData,TArray,pot_scale,"TR","Data",0,0,1001)
                    PlotTools.Print(AnalysisConfig.PlotPath("data_signalCats",sideband,"N4_tune"),mnvplotter,c)
                    subbedData.SetTitle("Backgrounded Subtracted Data")
                    # mnvplotter.DrawErrorSummary(subbedData,"TR",True,True,0)
                    mnvplotter.DrawErrorSummaryDerived(subbedData, "TR", True, True, 0, False, "", True, "", False, "HIST")
                    PlotTools.Print(AnalysisConfig.PlotPath("data_subbedErr",sideband,"N4_tune"),mnvplotter,c)
                    

        else:
            #assuing sideband_group is a tuple of name, and list of sidebands
            sideband = sideband_group[0]
            sidebands = sideband_group[1]
            data_hist,mc_hist = GetScaledDataMC(config["name"] if "name" in config else config,datafile,mcfile,sidebands[0])
            for _ in range(1,len(sidebands)):
                data_hist_tmp,mc_hist_tmp = GetScaledDataMC(config["name"] if "name" in config else config,datafile,mcfile,sidebands[_])
                data_hist.Add(data_hist_tmp)
                mc_hist.Add(mc_hist_tmp)
            MakePlot(data_hist,mc_hist,postfit_config)

    #make bkg subtracted data histogram

    datafile.Close()
    mcfile.Close()
    scalefile.Close()
 

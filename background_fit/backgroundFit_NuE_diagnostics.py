# Modified from the original 2-region backgroundFit script.
# Keeps the original RunMinimizer/error-band-CV procedure.
# Replaces the 2x2 analytic solve with a recipe-driven 4-region matrix solve.

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
from tools import Utilities, PlotTools
from config.UnfoldingConfig import HISTOGRAMS_TO_UNFOLD
from config.DrawingConfig import (
    SignalOnly,
    Default_Plot_Type,
    Default_Scale,
    DefaultPlotters,
    DefaultSlicer,
    PLOTS_TO_MAKE,
    SignalChargedBackground,
)

mnvplotter = PlotUtils.MnvPlotter()

from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS
mnvplotter.error_summary_group_map.clear()
for k, v in CONSOLIDATED_ERROR_GROUPS.items():
    vec = ROOT.vector("std::string")()
    for vs in v:
        vec.push_back(vs)
    mnvplotter.error_summary_group_map[k] = vec

# Stop ROOT from owning histograms.
ROOT.TH1.AddDirectory(False)


def SetFractionalUncertaintyYAxis(mnvplotter, ymin=0.0, ymax=1.0):
    mnvplotter.axis_minimum = ymin
    mnvplotter.axis_maximum = ymax


def ResetPlotterYAxis(mnvplotter):
    mnvplotter.axis_minimum = -1111
    mnvplotter.axis_maximum = -1111


def CheckBandCVsBinByBin(hist, label="", bands=None):
    if hist is None:
        return

    if bands is None:
        bands = [str(x) for x in hist.GetVertErrorBandNames()]

    print("\n===== Band CV bin-by-bin check:", label, "=====")

    for bandname in bands:
        if not hist.HasVertErrorBand(bandname):
            print(f"  {bandname}: missing")
            continue

        band = hist.GetVertErrorBand(bandname)
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
            print(f"  {bandname}: OK, band CV matches main CV bin-by-bin")
        else:
            print(f"  {bandname}: {ndiff} bins differ, maxdiff={maxdiff}")


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
    h_main.SetLineStyle(2)

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

    h_bandcv.Draw("HIST SAME")
    h_bandcv.Draw("P SAME")
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


def DrawAllBandsFromHistExplicit(h, out_prefix, bands=None):
    if h is None:
        return

    if bands is None:
        bands = [str(x) for x in h.GetVertErrorBandNames()]

    for bandname in bands:
        DrawBandUniversesFromHistExplicit(
            h,
            str(bandname),
            f"{out_prefix}_{h.GetName()}_{bandname}",
        )


def FillVertBandCVsFromParent(h):
    """
    Copy parent CV contents into each vertical error-band CV.
    This avoids empty/stale band CVs after derived histogram operations.
    """
    if h is None:
        return

    for bandname in h.GetErrorBandNames():
        band = h.GetVertErrorBand(bandname)
        if not band:
            continue

        for i in range(h.GetSize()):
            band.SetBinContent(i, h.GetBinContent(i))
            band.SetBinError(i, h.GetBinError(i))


def PrintScaleBandCVSyncSummary(scale_hists, label=""):
    print(f"\n===== Scale hist band-CV sync summary: {label} =====")
    for comp, hist in scale_hists.items():
        if hist is None:
            continue

        main_int = hist.Integral()
        print("\nScale hist:", comp, "main integral =", main_int)

        for bandname in hist.GetVertErrorBandNames():
            band = hist.GetVertErrorBand(bandname)
            ratio = band.Integral() / main_int if main_int else "NA"
            print("  {} bandCV/main = {}".format(bandname, ratio))


def MakeComparableMnvHXD(target_hist, scale_hists, y_axis=False):
    """
    Map 1D scale histograms onto the binning of target_hist, using x or y bin centers.
    Propagates vertical error-band universes.
    """
    new_scale = {comp: target_hist.Clone() for comp in scale_hists}
    xbins = target_hist.GetNbinsX() + 2

    for comp, sh in scale_hists.items():
        out = new_scale[comp]
        out.Reset()

        for i in range(target_hist.GetSize()):
            nx = i % xbins
            ny = i // xbins

            coord = (
                target_hist.GetYaxis().GetBinCenter(ny)
                if y_axis
                else target_hist.GetXaxis().GetBinCenter(nx)
            )
            k = sh.FindBin(coord)

            out.SetBinContent(i, sh.GetBinContent(k))
            out.SetBinError(i, sh.GetBinError(k))

            for bandname in sh.GetErrorBandNames():
                if bandname not in out.GetErrorBandNames():
                    continue

                out_band = out.GetVertErrorBand(bandname)
                src_band = sh.GetVertErrorBand(bandname)
                if not out_band or not src_band:
                    continue

                out_band.SetBinContent(i, sh.GetBinContent(k))
                out_band.SetBinError(i, sh.GetBinError(k))

                nh = min(out_band.GetNHists(), src_band.GetNHists())
                for u in range(nh):
                    out_u = out_band.GetHist(u)
                    src_u = src_band.GetHist(u)
                    out_u.SetBinContent(i, src_u.GetBinContent(k))
                    out_u.SetBinError(i, src_u.GetBinError(k))

        FillVertBandCVsFromParent(out)

    return new_scale


def WriteScaleToMnvH1D(out_hists, in_hists, scale_err=None, errorband=None, i=None):
    """
    Same flow as the original 2-region writer:
      - errorband=None writes main CV scale hist
      - errorband set, i=None writes the error-band CV/container
      - errorband set, i set writes universe i
    """
    for group in in_hists:
        if group not in out_hists:
            print("[WriteScaleToMnvH1D] missing output key", group)
            continue

        if errorband is None:
            universe_hist = out_hists[group]
        elif i is None:
            band = out_hists[group].GetVertErrorBand(errorband)
            if not band:
                print(f"[WriteScaleToMnvH1D] missing band {errorband} on {group}")
                continue
            universe_hist = band
        else:
            band = out_hists[group].GetVertErrorBand(errorband)
            if not band:
                print(f"[WriteScaleToMnvH1D] missing band {errorband} on {group}")
                continue
            universe_hist = band.GetHist(i)

        nb = min(universe_hist.GetNbinsX(), in_hists[group].GetNbinsX())
        for q in range(0, nb + 2):
            universe_hist.SetBinContent(q, in_hists[group].GetBinContent(q))
            universe_hist.SetBinError(q, in_hists[group].GetBinError(q))


def _hist_view(h, error_band=None, i=None):
    """
    Histogram view for the fit.

    Important for the original band-CV flow:
    when error_band is set and i is None, return a TH1D copy of the band's internal CV,
    not the full MnvVertErrorBand container. This avoids adding bands with different
    universe counts while preserving the 2-region CV/band/universe procedure.
    """
    if h is None:
        return None

    if error_band is None:
        return h

    if not h.HasVertErrorBand(error_band):
        return None

    band = h.GetVertErrorBand(error_band)

    if i is None:
        out = ROOT.TH1D(band)
        out.SetDirectory(0)
        return out

    return band.GetHist(i)


def _clone_empty_like(h, name):
    out = h.Clone(name)
    out.Reset()
    return out


def BuildComponentSums(region_holder, recipe, error_band=None, i=None):
    """
    Build one summed MC histogram per fitted component for one region.
    Components and category mapping are supplied by the matrix recipe.
    """
    href = _hist_view(region_holder.GetHist(), error_band, i)
    if href is None:
        raise RuntimeError(f"Could not get reference hist for region {region_holder.sideband}, band={error_band}, i={i}")

    out = {comp: _clone_empty_like(href, f"sum_{region_holder.sideband}_{comp}") for comp in recipe.components}

    for cate, hcate in region_holder.hists.items():
        if hcate is None or cate == "Total":
            continue
        if cate in getattr(recipe, "fixed_cates", []):
            continue

        comp = recipe.cate_to_comp(cate)
        if comp is None or comp not in out:
            continue

        huse = _hist_view(hcate, error_band, i)
        if huse is None:
            continue

        out[comp].Add(huse)

    return out


def BuildFixedSum(region_holder, recipe, error_band=None, i=None):
    """Sum fixed categories for one region. These are moved to the RHS: data - fixed."""
    href = _hist_view(region_holder.GetHist(), error_band, i)
    if href is None:
        raise RuntimeError(f"Could not get fixed reference hist for region {region_holder.sideband}, band={error_band}, i={i}")

    fixed = _clone_empty_like(href, f"fixed_{region_holder.sideband}")

    for cate in getattr(recipe, "fixed_cates", []):
        hcate = region_holder.hists.get(cate, None)
        if hcate is None:
            continue
        huse = _hist_view(hcate, error_band, i)
        if huse is None:
            continue
        fixed.Add(huse)

    return fixed


def SolveScalesMatrix(A_cols, b_vec, prior=None, clamp_nonneg=True):
    """
    Solve A x = b using SVD.
    No regularization is applied: all singular modes are kept, matching the current setup.
    """
    K = len(A_cols)
    R = len(b_vec)

    if prior is None:
        prior = [1.0] * K

    total_mc = sum(float(A_cols[k][r]) for k in range(K) for r in range(R))
    if total_mc == 0.0:
        return prior

    A = ROOT.TMatrixD(R, K)
    for r in range(R):
        for k in range(K):
            A[r][k] = float(A_cols[k][r])

    b = ROOT.TVectorD(R)
    for r in range(R):
        b[r] = float(b_vec[r])

    svd = ROOT.TDecompSVD(A)
    sig = svd.GetSig()
    nsv = sig.GetNrows()

    if nsv == 0:
        return prior

    svals = [float(sig[j]) for j in range(nsv)]
    if max(svals) <= 0.0:
        return prior

    U = svd.GetU()
    V = svd.GetV()

    # keep all singular modes with positive singular values
    x = [0.0] * K
    for j, sj in enumerate(svals):
        if sj <= 0.0:
            continue

        yj = 0.0
        for r in range(R):
            yj += float(U[r][j]) * float(b[r])

        coef = yj / sj
        for k in range(K):
            x[k] += float(V[k][j]) * coef

    if clamp_nonneg:
        x = [max(0.0, xi) for xi in x]

    return x

def _matrix_inverse_tsvd(A, kmax=None, tol=None):
    """
    Return pseudo-inverse A^+ for A with shape R x K.
    For square 4x4 full-rank case, this acts like an inverse.
    """
    R = A.GetNrows()
    K = A.GetNcols()

    svd = ROOT.TDecompSVD(A)
    sig = svd.GetSig()
    U = svd.GetU()
    V = svd.GetV()

    nsv = sig.GetNrows()
    svals = [float(sig[i]) for i in range(nsv)]
    smax = max(svals) if svals else 0.0

    if smax <= 0.0:
        return None

    if kmax is not None:
        keep = list(range(max(0, min(int(kmax), nsv))))
    elif tol is not None:
        keep = [i for i in range(nsv) if svals[i] >= float(tol) * smax]
    else:
        keep = list(range(nsv))

    if not keep:
        return None

    # A^+ has shape K x R
    Aplus = ROOT.TMatrixD(K, R)
    Aplus.Zero()

    # A^+ = V S^+ U^T
    for i in keep:
        si = svals[i]
        if si <= 0:
            continue

        inv_si = 1.0 / si

        for k in range(K):
            for r in range(R):
                Aplus[k][r] += float(V[k][i]) * inv_si * float(U[r][i])

    return Aplus

def _scale_stat_errors_from_data(A_cols, data_err_vec, kmax=None, tol=None):
    """
    Propagate data statistical errors through x = A^+ b.

    A_cols: list over components, each containing region yields.
            Shape: K x R in Python-list form.
    data_err_vec: data/fixed-RHS errors per region. Shape: R.

    Returns list of stat errors for each component.
    """
    K = len(A_cols)
    R = len(data_err_vec)

    A = ROOT.TMatrixD(R, K)
    for r in range(R):
        for k in range(K):
            A[r][k] = float(A_cols[k][r])

    Aplus = _matrix_inverse_tsvd(A, kmax=kmax, tol=tol)
    if Aplus is None:
        return [0.0] * K

    errs = []

    for k in range(K):
        var = 0.0
        for r in range(R):
            var += (float(Aplus[k][r]) * float(data_err_vec[r])) ** 2
        errs.append(math.sqrt(max(0.0, var)))

    return errs

def RunUniverseMinimizer4Region(data_holders, mc_holders, recipe, error_band=None, i=None):
    """
    Matrix version of the original RunUniverseMinimizer.

    Original 2-region solve:
      unknowns = signal, background
      equations = Signal region, dEdX region

    This version:
      unknowns = recipe.components, e.g. BkgPhoton, BkgCC, BkgOther, Signal
      equations = recipe.regions, e.g. Signal, dEdX, Eavail, Etheta
    """
    ref_region = recipe.regions[0]
    href = data_holders[ref_region].GetHist()

    tag = "CV" if error_band is None else (f"{error_band}_bandCV" if i is None else f"{error_band}_{i}")

    scales = {}
    for comp in recipe.components:
        scales[comp] = href.Clone(f"scale_{comp}_{tag}")
        scales[comp].Reset()

    # Build component sums, fixed sums, and data views by region.
    comp_sums = {}
    fixed_sums = {}
    data_views = {}

    for region in recipe.regions:
        comp_sums[region] = BuildComponentSums(mc_holders[region], recipe, error_band, i)
        fixed_sums[region] = BuildFixedSum(mc_holders[region], recipe, error_band, i)
        data_views[region] = data_holders[region].GetHist()

    # Pseudodata follows the original approach: replace data with MC total.
    if AnalysisConfig.pseudodata:
        for region in recipe.regions:
            mctot = mc_holders[region].hists.get("Total", None)
            if mctot is None:
                continue

            mct = _hist_view(mctot, error_band, i)
            if mct is None:
                continue

            for b in range(0, data_views[region].GetNbinsX() + 2):
                data_views[region].SetBinContent(b, mct.GetBinContent(b))
                data_views[region].SetBinError(b, mct.GetBinError(b))

    min_total_mc = getattr(recipe, "min_total_mc", None)
    if min_total_mc is None:
        min_total_mc = getattr(BackgroundFitConfig, "MIN_TOTAL_FLOATED_MC", 1.0)




    kmax_used = getattr(recipe, "kreg", None)
    tol_used = getattr(recipe, "tol", None)

    for b in range(0, href.GetNbinsX() + 2):
        b_vec = []
        b_err_vec = []

        for region in recipe.regions:
            data_val = data_views[region].GetBinContent(b)
            fixed_val = fixed_sums[region].GetBinContent(b)

            rhs = data_val - fixed_val
            b_vec.append(float(rhs))

            # Statistical uncertainty on RHS = data - fixed.
            # This includes data stat and fixed-category MC stat.
            data_err = data_views[region].GetBinError(b)
            fixed_err = fixed_sums[region].GetBinError(b)
            rhs_err = math.sqrt(data_err * data_err + fixed_err * fixed_err)

            b_err_vec.append(float(rhs_err))

        A_cols = []
        total_mc = 0.0

        for comp in recipe.components:
            col = []

            for region in recipe.regions:
                val = comp_sums[region][comp].GetBinContent(b)
                col.append(float(val))
                total_mc += float(val)

            A_cols.append(col)

        if total_mc < min_total_mc:
            x = [1.0] * len(recipe.components)
            xerr = [0.0] * len(recipe.components)
        else:
            x = SolveScalesMatrix(
                A_cols,
                b_vec,
                prior=[1.0] * len(recipe.components),
                clamp_nonneg=True,
            )

            xerr = _scale_stat_errors_from_data(
                A_cols,
                b_err_vec,
                kmax=kmax_used,
                tol=tol_used,
            )

        for ic, comp in enumerate(recipe.components):
            scales[comp].SetBinContent(b, x[ic])
            scales[comp].SetBinError(b, xerr[ic])
    # for b in range(0, href.GetNbinsX() + 2):
    #     b_vec = []
    #     for region in recipe.regions:
    #         rhs = data_views[region].GetBinContent(b) - fixed_sums[region].GetBinContent(b)
    #         b_vec.append(float(rhs))

    #     A_cols = []
    #     total_mc = 0.0
    #     for comp in recipe.components:
    #         col = []
    #         for region in recipe.regions:
    #             val = comp_sums[region][comp].GetBinContent(b)
    #             col.append(float(val))
    #             total_mc += float(val)
    #         A_cols.append(col)

    #     if total_mc < min_total_mc:
    #         x = [1.0] * len(recipe.components)
    #     else:
    #         x = SolveScalesMatrix(
    #             A_cols,
    #             b_vec,
    #             prior=[1.0] * len(recipe.components),
    #             clamp_nonneg=True,
    #         )

    #     for ic, comp in enumerate(recipe.components):
    #         scales[comp].SetBinContent(b, x[ic])
    #         scales[comp].SetBinError(b, 0.0)

    return scales


def RunMinimizer4Region(data_holders, mc_holders, recipe, scale_hists):
    """
    Same CV -> error-band CV -> universes flow as the original 2-region RunMinimizer.
    """
    scales = RunUniverseMinimizer4Region(data_holders, mc_holders, recipe)
    WriteScaleToMnvH1D(scale_hists, scales, None)
    print("Done with CV scale")

    any_hist = mc_holders[recipe.regions[0]].GetHist()

    for error_band in any_hist.GetErrorBandNames():
        scales = RunUniverseMinimizer4Region(data_holders, mc_holders, recipe, error_band)
        WriteScaleToMnvH1D(scale_hists, scales, None, error_band)
        print("Done with error band histogram {}".format(error_band))

        band = any_hist.GetVertErrorBand(error_band)
        if not band:
            continue

        for i in range(band.GetNHists()):
            scales = RunUniverseMinimizer4Region(data_holders, mc_holders, recipe, error_band, i)
            WriteScaleToMnvH1D(scale_hists, scales, None, error_band, i)


def ScaleCategories(hist_holder, scale_hists, recipe, prediction=False, signal_comp_name="Signal"):
    for cate, h in hist_holder.hists.items():
        if h is None or cate == "Total":
            continue

        if cate in getattr(recipe, "fixed_cates", []):
            FillVertBandCVsFromParent(h)
            continue

        comp = recipe.cate_to_comp(cate)
        if comp is None or comp not in scale_hists:
            FillVertBandCVsFromParent(h)
            continue

        # Prediction plot: tune backgrounds only, leave signal component nominal.
        if prediction and comp == signal_comp_name:
            FillVertBandCVsFromParent(h)
            continue

        try:
            h.Multiply(h, scale_hists[comp])
            FillVertBandCVsFromParent(h)
        except Exception:
            print("KeyError/multiply failure with {} in {}".format(cate, hist_holder.sideband))
            FillVertBandCVsFromParent(h)
            continue


def TuneMC(hist_holder, scale_hists, recipe, x_axis=False, y_axis=False, prediction=False, signal_comp_name="Signal"):
    if x_axis and y_axis:
        return None

    try:
        if x_axis or y_axis:
            comparable_scale = MakeComparableMnvHXD(hist_holder.GetHist(), scale_hists, y_axis)
            ScaleCategories(hist_holder, comparable_scale, recipe, prediction, signal_comp_name)
        else:
            ScaleCategories(hist_holder, scale_hists, recipe, prediction, signal_comp_name)
    except AttributeError:
        return False

    hist_holder.ResumTotal()
    return True


def BackgroundSubtraction(data_hists, mc_hists, pred_hists, recipe, signal_comp_name="Signal"):
    if not data_hists.POT_scaled:
        data_hists.POTScale(False)
    if not mc_hists.POT_scaled:
        mc_hists.POTScale(False)
    if not pred_hists.POT_scaled:
        pred_hists.POTScale(False)

    out_data = data_hists.GetHist().Clone()
    out_mc = pred_hists.hists["Total"].Clone()
    out_data.AddMissingErrorBandsAndFillWithCV(out_mc)

    for cate, h in mc_hists.hists.items():
        if h is None or cate == "Total":
            continue
        if cate in getattr(recipe, "fixed_cates", []):
            continue

        comp = recipe.cate_to_comp(cate)
        if comp is None:
            continue

        if comp != signal_comp_name:
            SubtractPoissonHistograms(out_data, mc_hists.hists[cate])
            SubtractPoissonHistograms(out_mc, pred_hists.hists[cate])

    return out_data, out_mc


def SubtractPoissonHistograms(h, h1):
    errors = []
    for i in range(h.GetSize()):
        errors.append(math.sqrt(h.GetBinError(i) ** 2 + h1.GetBinError(i) ** 2))
    h.Add(h1, -1)
    for i in range(h.GetSize()):
        h.SetBinError(i, errors[i])
    return h


def GetScaledDataMC(hist, datafile, mcfile, region, recipe, scale_hists, pot_scale, scaled_hist_name, signal_comp_name):
    data_hist = HistHolder(hist, datafile, region, False, pot_scale)
    mc_hist = HistHolder(hist, mcfile, region, True, pot_scale)
    pred_hist = HistHolder(hist, mcfile, region, True, pot_scale)

    fit_on_axis = scaled_hist_name.upper() in data_hist.plot_name.upper() or "estimator" in data_hist.plot_name.lower()

    if fit_on_axis:
        fit_on_yaxis = ("_" + scaled_hist_name).upper() in data_hist.plot_name.upper() or "_estimator" in data_hist.plot_name.lower()
        print("fit {} on {} axis".format(data_hist.plot_name, "y" if fit_on_yaxis else "x"))
        TuneMC(mc_hist, scale_hists, recipe, not fit_on_yaxis, fit_on_yaxis, False, signal_comp_name)
        TuneMC(pred_hist, scale_hists, recipe, not fit_on_yaxis, fit_on_yaxis, True, signal_comp_name)
    else:
        # For this 4-region matrix fit, scales are functions of the fitted estimator.
        # Apply them directly only to histograms on the fitted axis. For other variables,
        # this mirrors the modern matrix-fit behavior and avoids deriving a scalar average.
        print("not fitting {} on any axis; applying estimator scales directly".format(data_hist.plot_name))
        TuneMC(mc_hist, scale_hists, recipe, False, False, False, signal_comp_name)
        TuneMC(pred_hist, scale_hists, recipe, False, False, True, signal_comp_name)

    return data_hist, mc_hist, pred_hist


def MakePlot(data_hists, mc_hists, config):
    if not (data_hists.valid and mc_hists.valid):
        return False

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

    CanvasConfig = config.setdefault("canvasconfig", lambda x: True)
    PlotType = config.setdefault("plot_type", Default_Plot_Type)
    typeBool = PlotType != "migration" and PlotType != "category_hist" and PlotType != "hist2d"
    slicer = config.setdefault("slicer", DefaultSlicer(data_hists)) if typeBool else PlotTools.IdentitySlicer
    draw_seperate_legend = config.setdefault(
        "draw_seperate_legend",
        data_hists.dimension != 1 and (PlotType != "migration" or PlotType != "category_hist" or PlotType != "hist2d"),
    )

    try:
        custom_tag = config["tag"] + PlotType if "tag" in config else PlotType

        if PlotType == "custom":
            plotfunction, hists = config["getplotters"](data_hists, mc_hists)
        elif PlotType == "category_hist":
            if "args" in config:
                args = config["args"]
            elif "args" in DefaultPlotters[PlotType]:
                args = DefaultPlotters[PlotType]["args"]
            else:
                args = None
            categories = args[0]
            for category in categories:
                plotfunction, hists = DefaultPlotters[PlotType]["func"](data_hists, mc_hists, categories[category])
                PlotTools.MakeGridPlot(PlotTools.IdentitySlicer, plotfunction, hists, draw_seperate_legend=False, title=category)
                PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name, data_hists.sideband, category))
                print("plot {} made for category {}.".format(data_hists.plot_name, category))
            return True
        else:
            if "args" in config:
                args = config["args"]
            elif "args" in DefaultPlotters[PlotType]:
                args = DefaultPlotters[PlotType]["args"]
            else:
                args = None

            if args is None:
                plotfunction, hists = DefaultPlotters[PlotType]["func"](data_hists, mc_hists)
            else:
                plotfunction, hists = DefaultPlotters[PlotType]["func"](data_hists, mc_hists, *args)

            if PlotType == "2Dstacked":
                PlotTools.SumGridPlots(slicer, plotfunction, hists, draw_seperate_legend=False)
            else:
                PlotTools.MakeGridPlot(slicer, plotfunction, hists, draw_seperate_legend=False)

            PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name, data_hists.sideband, custom_tag))
            print("plot {} made.".format(data_hists.plot_name))

    except KeyError as e:
        print("plot {} not made.".format(data_hists.plot_name))
        print(e)
        return False

    return True


if __name__ == "__main__":
    playlist = AnalysisConfig.playlist
    type_path_map = {t: AnalysisConfig.SelectionHistoPath(playlist, t == "data", False) for t in AnalysisConfig.data_types}
    datafile, mcfile, pot_scale = Utilities.getFilesAndPOTScale(playlist, type_path_map, AnalysisConfig.ntuple_tag)

    background_fit_tag = AnalysisConfig.bkgTune_tag
    scalefile = ROOT.TFile.Open(AnalysisConfig.BackgroundFitPath(playlist, background_fit_tag), "RECREATE")
    BackgroundFitConfig.SetGlobalParameter(background_fit_tag)

    recipe = BackgroundFitConfig.GetMatrixRecipe(background_fit_tag)
    if recipe is None:
        raise RuntimeError(
            f"No matrix-fit recipe found for tag={background_fit_tag}. "
            "This 4-region script requires a recipe with regions, components, cate_to_comp, and fixed_cates."
        )

    print(f"[MatrixFit] Using recipe: {recipe.name} regions={recipe.regions} comps={recipe.components} kreg={getattr(recipe, 'kreg', None)}")

    fit_axis = recipe.hist_to_fit if getattr(recipe, "hist_to_fit", None) is not None else BackgroundFitConfig.HIST_TO_FIT
    fit_obs = recipe.hist_observable if getattr(recipe, "hist_observable", None) is not None else BackgroundFitConfig.HIST_OBSERVABLE

    ref_holder = HistHolder(fit_axis, mcfile, recipe.regions[0], True, pot_scale)
    ref_holder.POTScale(False)
    ref_for_binning = ref_holder.GetHist()
    scaled_hist_name = ref_holder.plot_name

    if "Signal" in recipe.components:
        signal_comp_name = "Signal"
    elif "SignalLike" in recipe.components:
        signal_comp_name = "SignalLike"
    else:
        raise RuntimeError(f"Cannot determine signal component from recipe.components={recipe.components}")

    data_holders = {}
    mc_holders = {}

    print("fit_obs =", fit_obs)
    for region in recipe.regions:
        data_holders[region] = HistHolder(fit_obs, datafile, region, False)
        mc_holders[region] = HistHolder(fit_obs, mcfile, region, True, pot_scale)
        mc_holders[region].POTScale(False)

    print("\n[DEBUG] category -> component mapping coverage")
    for reg in recipe.regions:
        keys = [k for k in mc_holders[reg].hists.keys() if k != "Total" and mc_holders[reg].hists[k] is not None]
        unmapped = [k for k in keys if recipe.cate_to_comp(k) is None and k not in getattr(recipe, "fixed_cates", [])]
        print(f"  region={reg} Nkeys={len(keys)} example={keys[:10]}")
        print(f"    unmapped N={len(unmapped)} example={unmapped[:15]}")

    # Scale histograms: one per fitted component.
    scale_hists = {}
    for comp in recipe.components:
        h = ref_for_binning.Clone(f"{comp}_Scale_Factor")
        h.Reset()
        h.GetYaxis().SetTitle("Scale Factor")
        h.SetTitle(f"{comp} Scale Factor")
        scale_hists[comp] = h

    RunMinimizer4Region(data_holders, mc_holders, recipe, scale_hists)
    PrintScaleBandCVSyncSummary(scale_hists, "after RunMinimizer4Region")

    scalefile.cd()

    axis_tag = ref_holder.plot_name
    for comp, hist in scale_hists.items():
        hist.SetXTitle(hist.GetXaxis().GetTitle())

        c1 = ROOT.TCanvas()
        mnvplotter.DrawMCWithErrorBand(hist)
        PlotTools.Print(AnalysisConfig.PlotPath(f"{axis_tag}_scales", comp, background_fit_tag), mnvplotter, c1)

        c1 = ROOT.TCanvas()
        SetFractionalUncertaintyYAxis(mnvplotter, 0.0, 1.0)
        mnvplotter.DrawErrorSummary(hist, "TR", True, True, 0)
        PlotTools.Print(AnalysisConfig.PlotPath(f"{axis_tag}_scale_errors", comp, background_fit_tag), mnvplotter, c1)
        ResetPlotterYAxis(mnvplotter)

        hist.Write(f"{comp}_Scale_Factor")

    # Predicted nominal signal shape in Signal region.
    signal_region = "Signal" if "Signal" in recipe.regions else recipe.regions[0]
    signalHist = HistHolder(fit_obs, mcfile, signal_region, True, pot_scale)
    signalHist.POTScale(False)

    mc_prediction = signalHist.GetHist().Clone("predicted_signal_total")
    mc_prediction.Reset()
    for cate, h in signalHist.hists.items():
        if h is None or cate == "Total":
            continue
        if cate in getattr(recipe, "fixed_cates", []):
            continue
        if recipe.cate_to_comp(cate) == signal_comp_name:
            mc_prediction.Add(h)

    # Write background-subtracted unfolding inputs.
    sideband_for_debug = "dEdX" if "dEdX" in recipe.regions else (recipe.regions[1] if len(recipe.regions) > 1 else None)

    for hist in HISTOGRAMS_TO_UNFOLD:
        if sideband_for_debug is not None:
            data_hist_sb, mc_hist_sb, pred_hist_sb = GetScaledDataMC(
                hist, datafile, mcfile, sideband_for_debug, recipe, scale_hists, pot_scale, scaled_hist_name, signal_comp_name
            )
            for _h in mc_hist_sb.hists:
                if mc_hist_sb.hists[_h]:
                    mc_hist_sb.hists[_h].Write(f"POSTFIT_{sideband_for_debug}_{_h}")

        data_hist, mc_hist, pred_hist = GetScaledDataMC(
            hist, datafile, mcfile, signal_region, recipe, scale_hists, pot_scale, scaled_hist_name, signal_comp_name
        )
        mc_hist.GetHist().Write(data_hist.plot_name)

        subbedData, subbedMC = BackgroundSubtraction(data_hist, mc_hist, pred_hist, recipe, signal_comp_name)
        subbedData.Write(data_hist.plot_name + "_data_bkgSubbed")
        subbedMC.Write(data_hist.plot_name + "_mc_bkgSubbed")
        mc_prediction.Write(data_hist.plot_name + "_predicted_Signal")

    # Reopen input files for plotting, as in the original script.
    type_path_map = {t: AnalysisConfig.SelectionHistoPath(playlist, t == "data", False) for t in AnalysisConfig.data_types}
    datafile, mcfile, pot_scale = Utilities.getFilesAndPOTScale(playlist, type_path_map, AnalysisConfig.ntuple_tag)

    for config in PLOTS_TO_MAKE:
        postfit_config = config.copy()
        postfit_config["tag"] = postfit_config.get("tag", "") + "postfit_"

        sideband_group = config.setdefault("sideband_group", recipe.regions)

        if isinstance(sideband_group, list):
            for sideband in sideband_group:
                data_hist, mc_hist, pred_hist = GetScaledDataMC(
                    config["name"] if "name" in config else config,
                    datafile,
                    mcfile,
                    sideband,
                    recipe,
                    scale_hists,
                    pot_scale,
                    scaled_hist_name,
                    signal_comp_name,
                )

                mc_postfit_config = postfit_config.copy()
                mc_postfit_config["tag"] = postfit_config.get("tag", "") + "mcHist_"

                pred_postfit_config = postfit_config.copy()
                pred_postfit_config["tag"] = postfit_config.get("tag", "") + "predHist_"

                sub_postfit_config = postfit_config.copy()
                sub_postfit_config["tag"] = postfit_config.get("tag", "") + "bkgSub_"

                if sideband == signal_region and AnalysisConfig.pseudodata:
                    MakePlot(data_holders[signal_region], mc_hist, mc_postfit_config)
                    MakePlot(data_holders[signal_region], pred_hist, pred_postfit_config)
                else:
                    MakePlot(data_hist, mc_hist, mc_postfit_config)
                    MakePlot(data_hist, pred_hist, pred_postfit_config)

                if sideband == signal_region:
                    subbedData, subbedMC = BackgroundSubtraction(data_hist, mc_hist, pred_hist, recipe, signal_comp_name)

                    c = ROOT.TCanvas("c_sub", "c_sub", 1200, 1000)
                    mnvplotter.DrawDataMCWithErrorBand(subbedData, subbedMC, 1.0, "TR")
                    PlotTools.Print(
                        AnalysisConfig.PlotPath(data_hist.plot_name, sideband, sub_postfit_config["tag"]),
                        mnvplotter,
                        c,
                    )

                    c2 = ROOT.TCanvas("c_sub_err", "c_sub_err", 1200, 1000)
                    SetFractionalUncertaintyYAxis(mnvplotter, 0.0, 1.0)
                    mnvplotter.DrawErrorSummary(subbedData, "TR", True, True, 0)
                    PlotTools.Print(
                        AnalysisConfig.PlotPath(data_hist.plot_name + "_err", sideband, sub_postfit_config["tag"]),
                        mnvplotter,
                        c2,
                    )
                    ResetPlotterYAxis(mnvplotter)

        else:
            sideband = sideband_group[0]
            sidebands = sideband_group[1]

            data_hist, mc_hist, pred_hist = GetScaledDataMC(
                config["name"] if "name" in config else config,
                datafile,
                mcfile,
                sidebands[0],
                recipe,
                scale_hists,
                pot_scale,
                scaled_hist_name,
                signal_comp_name,
            )

            for idx in range(1, len(sidebands)):
                data_hist_tmp, mc_hist_tmp, pred_hist_tmp = GetScaledDataMC(
                    config["name"] if "name" in config else config,
                    datafile,
                    mcfile,
                    sidebands[idx],
                    recipe,
                    scale_hists,
                    pot_scale,
                    scaled_hist_name,
                    signal_comp_name,
                )
                data_hist.Add(data_hist_tmp)
                mc_hist.Add(mc_hist_tmp)
                pred_hist.Add(pred_hist_tmp)

            mc_postfit_config = postfit_config.copy()
            mc_postfit_config["tag"] = postfit_config.get("tag", "") + "mcHist_"
            pred_postfit_config = postfit_config.copy()
            pred_postfit_config["tag"] = postfit_config.get("tag", "") + "predHist_"

            MakePlot(data_hist, mc_hist, mc_postfit_config)
            MakePlot(data_hist, pred_hist, pred_postfit_config)

    datafile.Close()
    mcfile.Close()
    scalefile.Close()

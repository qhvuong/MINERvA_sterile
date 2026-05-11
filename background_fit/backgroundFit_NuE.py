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

# NOTE:
# For recipe-driven matrix-fit, you should NOT rely on SIGNAL_DEFINITION inside fit/subtraction logic.
# Keep the import only if other plotting utilities still depend on it.
# from config.SignalDef import SIGNAL_DEFINITION

mnvplotter = PlotUtils.MnvPlotter()

# from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS
# mnvplotter.error_summary_group_map.clear()
# for k, v in CONSOLIDATED_ERROR_GROUPS.items():
#     vec = ROOT.vector("std::string")()
#     for vs in v:
#         vec.push_back(vs)
#     mnvplotter.error_summary_group_map[k] = vec

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
mnvplotter.error_summary_group_map.clear()
for k, v in group_map.items():
    vec = ROOT.vector("std::string")()
    for vs in v:
        vec.push_back(vs)
    mnvplotter.error_summary_group_map[k] = vec


# Stop ROOT from owning histograms (avoid segfaults / double-deletes)
ROOT.TH1.AddDirectory(False)


def FillScaleBandCVsFromCV(scale_hists):
    """
    After the main CV scale factors are filled, copy that CV into the
    CV histogram stored inside each vertical error band.

    This is needed because DrawErrorSummary/CalcCovMx expects the error
    band's internal CV histo to be non-empty.
    """
    for comp, h in scale_hists.items():
        for bandname in h.GetErrorBandNames():
            band = h.GetVertErrorBand(bandname)
            if not band:
                continue

            # MnvVertErrorBand acts histogram-like for its CV contents
            for b in range(0, h.GetNbinsX() + 2):
                band.SetBinContent(b, h.GetBinContent(b))
                band.SetBinError(b, h.GetBinError(b))


def FillVertBandCVsFromParent(h):
    """
    For a derived MnvH1D/MnvH2D, copy the parent CV contents into the
    internal CV histogram of every vertical error band.

    This avoids ErrBandEmptyCVError in PlotUtils covariance/error calculations.
    Universe histograms are not changed.
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


def PrintGroupedErrorsForBin(h, group_map, bins_to_check=None, label=""):
    if bins_to_check is None:
        bins_to_check = [2, 3, 4, 5, 6]

    available = set(str(x) for x in h.GetErrorBandNames())

    print(f"\n[GroupedErrors] {label} {h.GetName()}")
    for b in bins_to_check:
        cv = h.GetBinContent(b)
        print(f"  bin {b}: CV={cv:.6g}")

        for group, bands in group_map.items():
            abs2 = 0.0
            used = []

            for band in bands:
                if band not in available:
                    continue

                eb = h.GetVertErrorBand(band)
                if not eb:
                    continue

                herr = eb.GetErrorBand(False, False)
                err = herr.GetBinContent(b)
                abs2 += err * err
                used.append(band)

            abs_err = math.sqrt(abs2)
            frac = abs_err / abs(cv) if cv != 0 else 0.0

            print(
                f"    {group:25s} abs={abs_err:.6g} "
                f"frac={frac:.6g} used={len(used)}"
            )


# -----------------------------------------------------------------------------
# Helper: make scale histograms "comparable" to a target hist (useful for 2D axes)
# -----------------------------------------------------------------------------
def MakeComparableMnvHXD(target_hist, scale_hists, y_axis=False):
    """
    Map 1D scale histograms onto the binning of `target_hist` (1D or 2D),
    using x (or y) bin centers. Propagates vertical error-band universes.
    """
    new_scale = {comp: target_hist.Clone() for comp in scale_hists}

    xbins = target_hist.GetNbinsX() + 2  # include under/overflow

    for comp, sh in scale_hists.items():
        out = new_scale[comp]
        out.Reset()

        for i in range(target_hist.GetSize()):
            nx = i % xbins
            ny = i // xbins

            coord = (
                target_hist.GetYaxis().GetBinCenter(ny) if y_axis
                else target_hist.GetXaxis().GetBinCenter(nx)
            )
            k = sh.FindBin(coord)

            # Parent CV
            out.SetBinContent(i, sh.GetBinContent(k))
            out.SetBinError(i, sh.GetBinError(k))

            # Vertical error-band universes
            for bandname in sh.GetErrorBandNames():
                if bandname not in out.GetErrorBandNames():
                    continue

                out_band = out.GetVertErrorBand(bandname)
                src_band = sh.GetVertErrorBand(bandname)

                # Fill the error band's internal CV too
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


def WriteScaleToMnvH1D(out_hists, in_hists, errorband=None, uni_idx=None):
    for key, hin in in_hists.items():
        if key not in out_hists:
            continue

        hout = out_hists[key]

        if errorband is None:
            target = hout
        else:
            if uni_idx is None:
                raise RuntimeError(
                    f"Do not write fitted scales directly to error-band container for {errorband}. "
                    "Use FillScaleBandCVsFromCV for band CV, or pass a universe index."
                )

            band = hout.GetVertErrorBand(errorband)
            if not band:
                print(f"[WriteScaleToMnvH1D] Missing band {errorband} on {key}")
                continue

            target = band.GetHist(uni_idx)

        for b in range(0, target.GetNbinsX() + 2):
            target.SetBinContent(b, hin.GetBinContent(b))
            target.SetBinError(b, hin.GetBinError(b))


def _get_region_hist(hist_holder, error_band=None, uni_idx=None, use_universe=True):
    h = hist_holder.GetHist()
    if (not use_universe) or error_band is None:
        return h
    band = h.GetVertErrorBand(error_band)
    if not band:
        return None
    return band.GetHist(uni_idx) if uni_idx is not None else band


def _build_component_sums(recipe, mc_holder, error_band=None, uni_idx=None):
    """
    Build summed MC histograms per floated component for ONE region.
    Returns: dict[comp] -> MnvH1D (cloned & filled)
    """
    href = _get_region_hist(mc_holder, error_band, uni_idx)
    out = {comp: href.Clone(f"sum_{mc_holder.sideband}_{comp}") for comp in recipe.components}
    for comp in out:
        out[comp].Reset()

    for cate, hcate in mc_holder.hists.items():
        if hcate is None or cate == "Total":
            continue
        if cate in recipe.fixed_cates:
            continue

        comp = recipe.cate_to_comp(cate)
        if comp is None:
            continue
        if comp not in out:
            continue

        # get the correct universe view of this category histogram
        huse = hcate
        if error_band is not None:
            if uni_idx is None:
                huse = hcate.GetVertErrorBand(error_band)
            else:
                huse = hcate.GetVertErrorBand(error_band).GetHist(uni_idx)

        out[comp].Add(huse)

    return out


def _build_fixed_sum(recipe, mc_holder, error_band=None, uni_idx=None):
    """
    Sum of fixed categories for ONE region (moved to RHS: data - fixed).
    Returns: MnvH1D (cloned & filled)
    """
    href = _get_region_hist(mc_holder, error_band, uni_idx)
    fixed = href.Clone(f"fixed_{mc_holder.sideband}")
    fixed.Reset()

    if not recipe.fixed_cates:
        return fixed

    for cate in recipe.fixed_cates:
        if cate not in mc_holder.hists or mc_holder.hists[cate] is None:
            continue
        hcate = mc_holder.hists[cate]
        huse = hcate
        if error_band is not None:
            if uni_idx is None:
                huse = hcate.GetVertErrorBand(error_band)
            else:
                huse = hcate.GetVertErrorBand(error_band).GetHist(uni_idx)
        fixed.Add(huse)

    return fixed

def _solve_scales_per_bin_tsvd(
    A_cols,
    b_vec,
    prior=None,
    kmax=None,          # e.g. recipe.kreg == 3 means "keep 3 modes"
    tol=None,           # optional alternative: relative cutoff; if set, used to define keep set
    clamp_nonneg=True,
    min_total_mc=None,  # low-stat fallback threshold (total floated MC across all regions & comps)
    want_diag=False,
):
    """
    Solve A x ~= b by TSVD on A (R x K).

    A_cols: list length K, each is list length R: A_cols[k][r]
    b_vec : list length R

    Regularization:
      - If kmax is not None: keep top kmax singular values
      - Else if tol is not None: keep s_i >= tol*smax
      - Else: keep all

    Low-stat fallback:
      - If min_total_mc is not None and sum(|A|) (or sum(A) if nonnegative) < threshold: return prior

    Diagnostics (if want_diag=True):
      returns (x, diag_dict)
    """
    K = len(A_cols)
    R = len(b_vec)

    if prior is None:
        prior = [1.0] * K

    # Early exit: no MC content at all
    total_mc = 0.0
    for k in range(K):
        for r in range(R):
            total_mc += float(A_cols[k][r])

    if min_total_mc is not None and total_mc < float(min_total_mc):
        diag = {
            "used_fallback": True,
            "total_mc": total_mc,
            "nsv": min(R, K),
            "n_kept": 0,
            "svals": [],
            "smax": 0.0,
            "smin_kept": 0.0,
        }
        return (prior, diag) if want_diag else prior

    if total_mc == 0.0:
        diag = {
            "used_fallback": True,
            "total_mc": 0.0,
            "nsv": min(R, K),
            "n_kept": 0,
            "svals": [],
            "smax": 0.0,
            "smin_kept": 0.0,
        }
        return (prior, diag) if want_diag else prior

    # Build A (R x K) and b (R)
    A = ROOT.TMatrixD(R, K)
    for r in range(R):
        for k in range(K):
            A[r][k] = float(A_cols[k][r])

    b = ROOT.TVectorD(R)
    for r in range(R):
        b[r] = float(b_vec[r])

    svd = ROOT.TDecompSVD(A)

    sig = svd.GetSig()  # length min(R,K)
    nsv = sig.GetNrows()
    if nsv == 0:
        diag = {"used_fallback": True, "total_mc": total_mc, "nsv": 0, "n_kept": 0, "svals": [], "smax": 0.0, "smin_kept": 0.0}
        return (prior, diag) if want_diag else prior

    svals = [float(sig[i]) for i in range(nsv)]
    smax = max(svals) if svals else 0.0
    if smax <= 0.0:
        diag = {"used_fallback": True, "total_mc": total_mc, "nsv": nsv, "n_kept": 0, "svals": svals, "smax": smax, "smin_kept": 0.0}
        return (prior, diag) if want_diag else prior

    # Decide which modes to keep
    if kmax is not None:
        kk = max(0, min(int(kmax), nsv))
        keep = list(range(kk))  # singular values are ordered descending in ROOT
    elif tol is not None:
        keep = [i for i in range(nsv) if svals[i] >= float(tol) * smax]
    else:
        keep = list(range(nsv))

    if len(keep) == 0:
        diag = {"used_fallback": True, "total_mc": total_mc, "nsv": nsv, "n_kept": 0, "svals": svals, "smax": smax, "smin_kept": 0.0}
        return (prior, diag) if want_diag else prior

    U = svd.GetU()  # (R x nsv)
    V = svd.GetV()  # (K x nsv)

    # y = U^T b
    y = [0.0] * nsv
    for i in range(nsv):
        acc = 0.0
        for r in range(R):
            acc += float(U[r][i]) * float(b[r])
        y[i] = acc

    # x = sum_{i in keep} V[:,i] * (y_i / s_i)
    x = [0.0] * K
    for i in keep:
        si = svals[i]
        if si == 0.0:
            continue
        coef = y[i] / si
        for k in range(K):
            x[k] += float(V[k][i]) * coef

    if clamp_nonneg:
        x = [max(0.0, xi) for xi in x]

    smin_kept = min(svals[i] for i in keep) if keep else 0.0
    diag = {
        "used_fallback": False,
        "total_mc": total_mc,
        "nsv": nsv,
        "n_kept": len(keep),
        "svals": svals,
        "smax": smax,
        "smin_kept": smin_kept,
    }
    return (x, diag) if want_diag else x


def RunUniverseMinimizer_MatrixFit(
    recipe,
    data_holders,
    mc_holders,
    error_band=None,
    uni_idx=None,
):
    """
    Solve for scale factors per bin for all floated components in the NuE
    multi-region matrix-fit recipe.

    data_holders: dict region -> HistHolder for data
    mc_holders  : dict region -> HistHolder for MC

    Returns
    -------
    scales: dict[component] -> MnvH1D scale histogram
    """
    # Use CV data histogram for binning. Data should not have universe variations.
    href = _get_region_hist(data_holders[recipe.regions[0]], use_universe=False)

    tag = "CV" if error_band is None else (
        f"{error_band}[bandCV]" if uni_idx is None else f"{error_band}[{uni_idx}]"
    )

    # Output scale histograms, one per floated component
    scales = {}
    for comp in recipe.components:
        scales[comp] = href.Clone(f"scale_{comp}_{tag}")
        scales[comp].Reset()

    # Diagnostics bookkeeping
    diag = {
        "n_bins_total": 0,
        "n_bins_fallback": 0,
        "n_bins_dropped": 0,
        "drop_by_nkept": {},
        "fallback_bins": [],
        "dropped_bins": [],
        "worst_cond_bins": [],
    }

    min_total_mc = getattr(recipe, "min_total_mc", None)
    if min_total_mc is None:
        min_total_mc = getattr(BackgroundFitConfig, "MIN_TOTAL_FLOATED_MC", 1.0)

    # TSVD settings
    # kmax_used = getattr(recipe, "kmax", None)
    # tol = getattr(recipe, "tol", 0.02)
    kmax_used = None
    tol = 0.0

    # Build component sums, fixed sums, and data views for each region
    comp_sums_by_region = {}
    fixed_by_region = {}
    data_by_region = {}

    for reg in recipe.regions:
        comp_sums_by_region[reg] = _build_component_sums(
            recipe,
            mc_holders[reg],
            error_band=error_band,
            uni_idx=uni_idx,
        )

        fixed_by_region[reg] = _build_fixed_sum(
            recipe,
            mc_holders[reg],
            error_band=error_band,
            uni_idx=uni_idx,
        )

        # Data always uses CV. Only MC gets shifted by universes.
        data_by_region[reg] = _get_region_hist(
            data_holders[reg],
            use_universe=False,
        )

    # Pseudodata option: replace data with MC total for the same universe view
    if AnalysisConfig.pseudodata:
        for reg in recipe.regions:
            mctot = mc_holders[reg].hists.get("Total", None)
            if mctot is None:
                continue

            if error_band is None:
                mct = mctot
            else:
                band = mctot.GetVertErrorBand(error_band)
                if not band:
                    continue
                mct = band if uni_idx is None else band.GetHist(uni_idx)

            for b in range(0, data_by_region[reg].GetNbinsX() + 2):
                data_by_region[reg].SetBinContent(b, mct.GetBinContent(b))
                data_by_region[reg].SetBinError(b, mct.GetBinError(b))

    nb = href.GetNbinsX()

    for b in range(0, nb + 2):  # include underflow and overflow
        # RHS vector: data - fixed components, one entry per region
        b_vec = []
        for reg in recipe.regions:
            bval = (
                data_by_region[reg].GetBinContent(b)
                - fixed_by_region[reg].GetBinContent(b)
            )
            b_vec.append(float(bval))

        # A matrix stored by columns:
        # A_cols[component_index][region_index] = floated MC component yield
        A_cols = []
        total_mc = 0.0

        for comp in recipe.components:
            col = []
            for reg in recipe.regions:
                val = comp_sums_by_region[reg][comp].GetBinContent(b)
                col.append(float(val))
                total_mc += float(val)
            A_cols.append(col)

        if total_mc < min_total_mc:
            x = [1.0] * len(recipe.components)
            d = {
                "used_fallback": True,
                "total_mc": total_mc,
                "n_kept": 0,
                "smax": 0.0,
                "smin_kept": 0.0,
                "svals": [],
            }
        else:
            x, d = _solve_scales_per_bin_tsvd(
                A_cols,
                b_vec,
                prior=[1.0] * len(recipe.components),
                kmax=kmax_used,
                tol=tol,
                clamp_nonneg=True,
                min_total_mc=None,
                want_diag=True,
            )

        diag["n_bins_total"] += 1

        xcenter = (
            href.GetXaxis().GetBinCenter(b)
            if 1 <= b <= href.GetNbinsX()
            else float("nan")
        )

        if d["used_fallback"]:
            diag["n_bins_fallback"] += 1
            diag["fallback_bins"].append((b, xcenter, d["total_mc"]))
        else:
            full_rank = min(len(recipe.regions), len(recipe.components))

            if d["n_kept"] < full_rank:
                diag["n_bins_dropped"] += 1
                diag["drop_by_nkept"][d["n_kept"]] = (
                    diag["drop_by_nkept"].get(d["n_kept"], 0) + 1
                )
                diag["dropped_bins"].append(
                    (b, xcenter, d["n_kept"], d["svals"], d["total_mc"])
                )

            if 1 <= b <= href.GetNbinsX():
                ratio = (
                    d["smin_kept"] / d["smax"]
                    if d["smax"] > 0 and d["smin_kept"] > 0
                    else 0.0
                )
                diag["worst_cond_bins"].append((b, xcenter, ratio, d["n_kept"]))

        # Store fitted scale factors
        for ic, comp in enumerate(recipe.components):
            scales[comp].SetBinContent(b, x[ic])
            scales[comp].SetBinError(b, 0.0)

    # # Diagnostics summary
    # reg_mode = []
    # if kmax_used is not None:
    #     reg_mode.append(f"kmax={kmax_used}")
    # if tol is not None:
    #     reg_mode.append(f"tol={tol}")
    # reg_mode_str = " ".join(reg_mode) if reg_mode else "no_tsvd_cut"

    # print(
    #     f"[TSVD-DIAG {tag}] bins={diag['n_bins_total']} "
    #     f"fallback={diag['n_bins_fallback']} "
    #     f"dropped={diag['n_bins_dropped']} "
    #     f"min_total_mc={min_total_mc} {reg_mode_str}"
    # )

    # if diag["drop_by_nkept"]:
    #     items = ", ".join(
    #         f"kept{nk}:{ct}" for nk, ct in sorted(diag["drop_by_nkept"].items())
    #     )
    #     print(f"[TSVD-DIAG {tag}] dropped breakdown: {items}")

    # diag["worst_cond_bins"].sort(key=lambda t: t[2])
    # for bb, xc, ratio, nk in diag["worst_cond_bins"][:5]:
    #     if math.isnan(xc):
    #         print(
    #             f"[TSVD-DIAG {tag}] worst bin b={bb} "
    #             f"smin/smax={ratio:.3g} n_kept={nk}"
    #         )
    #     else:
    #         print(
    #             f"[TSVD-DIAG {tag}] worst bin b={bb} x={xc:.3g} "
    #             f"smin/smax={ratio:.3g} n_kept={nk}"
    #         )

    # if diag["n_bins_fallback"] > 0:
    #     for bb, xc, tmc in diag["fallback_bins"][:10]:
    #         if math.isnan(xc):
    #             print(f"[TSVD-DIAG {tag}] fallback bin b={bb} totalMC={tmc:.3g}")
    #         else:
    #             print(
    #                 f"[TSVD-DIAG {tag}] fallback bin b={bb} "
    #                 f"x={xc:.3g} totalMC={tmc:.3g}"
    #             )

    return scales

def RunMinimizer_MatrixFit(recipe, data_holders, mc_holders, scale_hists_out):
    # --------------------
    # CV fit
    # --------------------
    cv_scales = RunUniverseMinimizer_MatrixFit(
        recipe,
        data_holders,
        mc_holders,
        error_band=None,
        uni_idx=None,
    )

    WriteScaleToMnvH1D(
        scale_hists_out,
        cv_scales,
        errorband=None,
        uni_idx=None,
    )
    FillScaleBandCVsFromCV(scale_hists_out)
    print("[MatrixFit] Done with CV scales")

    any_hist = mc_holders[recipe.regions[0]].GetHist()

    # --------------------
    # Universe fits
    # --------------------
    for error_band in any_hist.GetErrorBandNames():
        band = any_hist.GetVertErrorBand(error_band)
        if not band:
            print(f"[MatrixFit] Skip missing band {error_band}")
            continue

        nuni = band.GetNHists()

        for i in range(nuni):
            uni_scales = RunUniverseMinimizer_MatrixFit(
                recipe,
                data_holders,
                mc_holders,
                error_band=error_band,
                uni_idx=i,
            )

            WriteScaleToMnvH1D(
                scale_hists_out,
                uni_scales,
                errorband=error_band,
                uni_idx=i,
            )

        print(f"[MatrixFit] Done with universes for band {error_band} (N={nuni})")

    # Safe final refresh of each band's internal CV
    FillScaleBandCVsFromCV(scale_hists_out)

# -------------------------------------------------------------------------
# NEW: apply scales to MC holders using recipe mapping
# -------------------------------------------------------------------------
def ScaleCategories_Recipe(recipe, hist_holder, scale_hists_mapped, prediction=False, signal_comp_name="Signal"):
    """
    Multiply each category histogram by the appropriate scale histogram,
    based on recipe.cate_to_comp.

    prediction=False: scale all floated comps
    prediction=True : scale only background comps, leave signal comp untouched
    """
    for cate, h in hist_holder.hists.items():
        if h is None or cate == "Total":
            continue

        if cate in recipe.fixed_cates:
            # Fixed categories are still used in plots/subtraction, so make sure
            # their error-band CVs are valid.
            FillVertBandCVsFromParent(h)
            continue

        comp = recipe.cate_to_comp(cate)
        if comp is None:
            FillVertBandCVsFromParent(h)
            continue

        if comp not in scale_hists_mapped:
            FillVertBandCVsFromParent(h)
            continue

        if prediction and comp == signal_comp_name:
            FillVertBandCVsFromParent(h)
            continue

        try:
            h.Multiply(h, scale_hists_mapped[comp])
            FillVertBandCVsFromParent(h)
        except Exception:
            print(f"[MatrixFit] Multiply failed for cate={cate} region={hist_holder.sideband} comp={comp}")
            FillVertBandCVsFromParent(h)
            continue


def TuneMC_Recipe(recipe, hist_holder, scale_hists, x_axis=False, y_axis=False, prediction=False, signal_comp_name="Signal"):
    if x_axis and y_axis:
        return False

    if x_axis or y_axis:
        mapped = MakeComparableMnvHXD(hist_holder.GetHist(), scale_hists, y_axis=y_axis)
        ScaleCategories_Recipe(
            recipe, hist_holder, mapped,
            prediction=prediction,
            signal_comp_name=signal_comp_name
        )
    else:
        ScaleCategories_Recipe(
            recipe, hist_holder, scale_hists,
            prediction=prediction,
            signal_comp_name=signal_comp_name
        )

    hist_holder.ResumTotal()

    for cate, h in hist_holder.hists.items():
        if h is not None:
            FillVertBandCVsFromParent(h)

    return True


# -------------------------------------------------------------------------
# NEW: recipe-driven background subtraction
# -------------------------------------------------------------------------
def BackgroundSubtraction_Recipe(recipe, data_hists, tuned_mc_hists, pred_hists, signal_comp_name="Signal"):
    """
    Subtract everything that is NOT mapped to the signal component.
    tuned_mc_hists: tuned MC used to subtract from data
    pred_hists: prediction histogram holder (typically tuned only bkg or tuned all depending on your convention)
    """
    data_hists.POTScale(False)
    tuned_mc_hists.POTScale(False)
    pred_hists.POTScale(False)

    out_data = data_hists.GetHist().Clone()
    out_mc = pred_hists.hists["Total"].Clone()
    out_data.AddMissingErrorBandsAndFillWithCV(out_mc)

    for cate, h in tuned_mc_hists.hists.items():
        if h is None or cate == "Total":
            continue
        if cate in recipe.fixed_cates:
            continue

        comp = recipe.cate_to_comp(cate)
        if comp is None:
            continue

        if comp != signal_comp_name:
            SubtractPoissonHistograms(out_data, tuned_mc_hists.hists[cate])
            SubtractPoissonHistograms(out_mc,   pred_hists.hists[cate])

    return out_data, out_mc


# -------------------------------------------------------------------------
# NEW: scaled data/mc getter (recipe-driven)
# -------------------------------------------------------------------------
def GetScaledDataMC_Recipe(hist, datafile, mcfile, region, recipe, scale_hists, pot_scale, scaled_hist_name=None):
    data_hist = HistHolder(hist, datafile, region, False, pot_scale)
    mc_hist   = HistHolder(hist, mcfile,   region, True,  pot_scale)
    pred_hist = HistHolder(hist, mcfile,   region, True,  pot_scale)

    # Determine if this plot is along the fitted axis
    fit_on_axis = False
    fit_on_yaxis = False
    if scaled_hist_name is not None:
        fit_on_axis = (scaled_hist_name.upper() in data_hist.plot_name.upper()) or ("estimator" in data_hist.plot_name.lower())
        if fit_on_axis:
            fit_on_yaxis = ("_" + scaled_hist_name).upper() in data_hist.plot_name.upper() or "_estimator" in data_hist.plot_name.lower()

    if fit_on_axis:
        print(f"[MatrixFit] fit {data_hist.plot_name} on {'y' if fit_on_yaxis else 'x'} axis")
        TuneMC_Recipe(recipe, mc_hist,  scale_hists, x_axis=not fit_on_yaxis, y_axis=fit_on_yaxis, prediction=False)
        TuneMC_Recipe(recipe, pred_hist, scale_hists, x_axis=not fit_on_yaxis, y_axis=fit_on_yaxis, prediction=True)
    else:
        # In matrix-fit world, you usually still apply the same per-bin scales even if variable != fit variable.
        # That is: scale factors are functions of estimator, not of this variable. So no axis mapping.
        print(f"[MatrixFit] not fitting {data_hist.plot_name} on axis; applying 1D estimator scales without mapping")
        TuneMC_Recipe(recipe, mc_hist,  scale_hists, x_axis=False, y_axis=False, prediction=False)
        TuneMC_Recipe(recipe, pred_hist, scale_hists, x_axis=False, y_axis=False, prediction=True)

    return data_hist, mc_hist, pred_hist

def GetBackground_Recipe(recipe, mc_hists, signal_comp_name="Signal"):
    out_bkg = mc_hists.hists["Total"].Clone("bkgTotal")
    out_bkg.Reset()

    for cate, h in mc_hists.hists.items():
        if h is None or cate == "Total":
            continue
        if cate in recipe.fixed_cates:
            # fixed categories are still "background" for plotting purposes unless you want otherwise
            out_bkg.Add(h)
            continue

        comp = recipe.cate_to_comp(cate)
        if comp is None:
            continue

        if comp != signal_comp_name:
            out_bkg.Add(h)

    return out_bkg

def SubtractPoissonHistograms(h, h1):
    # save combined errors first, then subtract, then apply errors
    errs = []
    for i in range(h.GetSize()):
        errs.append(math.sqrt(h.GetBinError(i)**2 + h1.GetBinError(i)**2))
    h.Add(h1, -1)
    for i in range(h.GetSize()):
        h.SetBinError(i, errs[i])
    return h

def GetSignalCompName(recipe):
    if "Signal" in recipe.components:
        return "Signal"
    if "SignalLike" in recipe.components:
        return "SignalLike"
    raise RuntimeError(f"Cannot determine signal component name from {recipe.components}")


# def MakeRatio_Recipe(
#     recipe,
#     signalHist, sidebandHist,
#     normsignalHist, normsidebandHist,
#     config,
#     signal_comp_name="Signal",
#     include_fixed_in_bkg=True,
#     flux_band_name="Flux",
# ):
#     # scale
#     if "scale" in config:
#         config["scale"](signalHist)
#         config["scale"](sidebandHist)
#         config["scale"](normsignalHist)
#         config["scale"](normsidebandHist)
#     else:
#         Default_Scale(signalHist)
#         Default_Scale(sidebandHist)
#         Default_Scale(normsignalHist)
#         Default_Scale(normsidebandHist)

#     def _sum_background(holder, name):
#         out = holder.hists["Total"].Clone(name)
#         out.Reset()

#         for cate, h in holder.hists.items():
#             if h is None or cate == "Total":
#                 continue

#             if cate in recipe.fixed_cates:
#                 if include_fixed_in_bkg:
#                     out.Add(h)
#                 continue

#             comp = recipe.cate_to_comp(cate)
#             if comp is None:
#                 continue

#             if comp != signal_comp_name:
#                 out.Add(h)

#         return out

#     sig_bkg     = _sum_background(signalHist,     "sig_bkg")
#     sid_bkg     = _sum_background(sidebandHist,   "sid_bkg")
#     normsig_bkg = _sum_background(normsignalHist, "normsig_bkg")
#     normsid_bkg = _sum_background(normsidebandHist,"normsid_bkg")

#     # Optional debug: Flux band checks
#     def _has_band(h, band):
#         try:
#             return band in list(h.GetErrorBandNames())
#         except Exception:
#             return False

#     if _has_band(sig_bkg, flux_band_name) and _has_band(sid_bkg, flux_band_name):
#         c1 = ROOT.TCanvas()
#         sig_bkg.GetVertErrorBand(flux_band_name).DrawAll("hist", True)
#         c1.Print(f"{signalHist.plot_name}_post_tuneSigBkg{flux_band_name}.png")

#         c1 = ROOT.TCanvas()
#         sid_bkg.GetVertErrorBand(flux_band_name).DrawAll("hist", True)
#         c1.Print(f"{signalHist.plot_name}_post_tuneSidBkg{flux_band_name}.png")

#         # Compare flux error bands pre/post tune (normalized ratio)
#         c1 = ROOT.TCanvas()
#         fluxerr = sig_bkg.GetVertErrorBand(flux_band_name).GetErrorBand(True, False).Clone()
#         normfluxerr = normsig_bkg.GetVertErrorBand(flux_band_name).GetErrorBand(True, False).Clone()
#         normfluxerr.Divide(normfluxerr, fluxerr)
#         normfluxerr.Draw()
#         c1.Print(f"{signalHist.plot_name}_sig_bkgPostTune_{flux_band_name}errband.png")

#         c1 = ROOT.TCanvas()
#         fluxerr = sid_bkg.GetVertErrorBand(flux_band_name).GetErrorBand(True, False).Clone()
#         normfluxerr = normsid_bkg.GetVertErrorBand(flux_band_name).GetErrorBand(True, False).Clone()
#         normfluxerr.Divide(normfluxerr, fluxerr)
#         normfluxerr.Draw()
#         c1.Print(f"{signalHist.plot_name}_sid_bkgPostTune_{flux_band_name}errband.png")

#     # Ratio of normalized background shapes: Signal-region-bkg / Sideband-region-bkg
#     isig = sig_bkg.Integral()
#     isd  = sid_bkg.Integral()
#     if isig <= 0 or isd <= 0:
#         print(f"[MakeRatio_Recipe] Skip ratio: zero integral (sig={isig}, sid={isd})")
#         return

#     sig_bkg.Scale(1.0 / isig)
#     sid_bkg.Scale(1.0 / isd)

#     sig_bkg.Divide(sig_bkg, sid_bkg)
#     sig_bkg.GetXaxis().SetTitle("E_{available} + E_{lepton}")
#     sig_bkg.GetYaxis().SetTitle("Ratio")
#     sig_bkg.SetTitle("Signal/Sideband Background Ratio")

#     mnvplotter.DrawMCWithErrorBand(sig_bkg)
#     PlotTools.Print(AnalysisConfig.PlotPath("EN4_ratio", "Combined", "bkgratio_matrixfit"))

#     c1 = ROOT.TCanvas()
#     mnvplotter.DrawErrorSummary(sig_bkg, "TR", True, True, 0)
#     c1.Print(f"{signalHist.plot_name}_post_tuneSigBkgErrSummary.png")


def MakePlot(data_hists, mc_hists, config):
    if not (data_hists.valid and mc_hists.valid):
        return False

    # Apply scaling only if neither is POT-scaled yet
    if not mc_hists.POT_scaled and not data_hists.POT_scaled:
        if "scale" in config:
            config["scale"](data_hists)
            config["scale"](mc_hists)
        else:
            Default_Scale(data_hists)
            Default_Scale(mc_hists)

    CanvasConfig = config.setdefault("canvasconfig", lambda x: True)
    PlotType = config.setdefault("plot_type", Default_Plot_Type)

    # Certain plot types expect identity slicing
    typeBool = PlotType not in ("migration", "category_hist", "hist2d")
    slicer = config.setdefault("slicer", DefaultSlicer(data_hists)) if typeBool else PlotTools.IdentitySlicer

    # Separate legend only for multi-panel plots where it helps
    draw_seperate_legend = config.setdefault(
        "draw_seperate_legend",
        data_hists.dimension != 1 and (PlotType not in ("migration", "category_hist", "hist2d"))
    )

    try:
        custom_tag = (config.get("tag", "") + PlotType) if "tag" in config else PlotType

        if PlotType == "custom":
            plotfunction, hists = config["getplotters"](data_hists, mc_hists)

        elif PlotType == "category_hist":
            # category_hist is special: loop categories and print each separately
            if "args" in config:
                args = config["args"]
            elif "args" in DefaultPlotters[PlotType]:
                args = DefaultPlotters[PlotType]["args"]
            else:
                args = None

            categories = args[0]
            for category in categories:
                plotfunction, hists = DefaultPlotters[PlotType]["func"](data_hists, mc_hists, categories[category])
                PlotTools.MakeGridPlot(
                    PlotTools.IdentitySlicer,
                    plotfunction,
                    hists,
                    draw_seperate_legend=False,
                    title=category
                )
                PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name, data_hists.sideband, category))
                print(f"plot {data_hists.plot_name} made for category {category}.")

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
        print(f"plot {data_hists.plot_name} made.")

    except KeyError as e:
        print(f"plot {data_hists.plot_name} not made.")
        print(e)
        return False

    return True


if __name__ == "__main__":
    # -------------------- INPUT KNOBS --------------------
    playlist = AnalysisConfig.playlist
    type_path_map = {
        t: AnalysisConfig.SelectionHistoPath(playlist, t == "data", False)
        for t in AnalysisConfig.data_types
    }
    datafile, mcfile, pot_scale = Utilities.getFilesAndPOTScale(
        playlist, type_path_map, AnalysisConfig.ntuple_tag
    )

    # -------------------- OUTPUT KNOBS --------------------
    background_fit_tag = AnalysisConfig.bkgTune_tag
    scalefile = ROOT.TFile.Open(
        AnalysisConfig.BackgroundFitPath(playlist, background_fit_tag),
        "RECREATE"
    )

    # Apply legacy global parameters (HIST_TO_FIT, etc.)
    BackgroundFitConfig.SetGlobalParameter(background_fit_tag)

    # Pick matrix-fit recipe (CCnue vs nueElastic)
    recipe = BackgroundFitConfig.GetMatrixRecipe(background_fit_tag)
    if recipe is None:
        raise RuntimeError(
            f"No matrix-fit recipe for tag={background_fit_tag}. "
            f"Available: {list(getattr(BackgroundFitConfig,'MATRIX_RECIPE_MAP',{}).keys())}"
        )
    print(f"[MatrixFit] Using recipe: {recipe.name} regions={recipe.regions} comps={recipe.components} kreg={recipe.kreg}")

    # ------------------------------------------------------------
    # Determine fit axis (EN4 for CCnue, Eel for nu+e)
    # ------------------------------------------------------------
    fit_axis = (
        recipe.hist_to_fit
        if hasattr(recipe, "hist_to_fit") and recipe.hist_to_fit is not None
        else BackgroundFitConfig.HIST_TO_FIT
    )

    # This holder defines the binning of the scale-factor histograms
    ref_holder = HistHolder(fit_axis, mcfile, recipe.regions[0], True, pot_scale)
    ref_for_binning = ref_holder.GetHist()
    scaled_hist_name = ref_holder.plot_name

    # Determine which component name is "signal" in this recipe
    if "Signal" in recipe.components:
        signal_comp_name = "Signal"
    elif "SignalLike" in recipe.components:
        signal_comp_name = "SignalLike"
    else:
        raise RuntimeError(f"Cannot determine signal component name from recipe.components={recipe.components}")

    # -------------------- BUILD REGION HIST HOLDERS (DATA/MC) --------------------
    # Dict keyed by region name, because recipes differ between CCnue and nu+e
    data_holders = {}
    mc_holders = {}
    fit_obs = recipe.hist_observable if recipe.hist_observable is not None else BackgroundFitConfig.HIST_OBSERVABLE
    print("fit_obs =", fit_obs)

    for region in recipe.regions:
        data_holders[region] = HistHolder(fit_obs, datafile, region, False)
        mc_holders[region]   = HistHolder(fit_obs, mcfile,   region, True, pot_scale)
        mc_holders[region].POTScale(False)

    print("\n[DEBUG] category -> component mapping coverage")
    for reg in recipe.regions:
        keys = [k for k in mc_holders[reg].hists.keys() if k not in ("Total",) and mc_holders[reg].hists[k] is not None]
        print(f"  region={reg} Nkeys={len(keys)} example={keys[:10]}")
        unmapped = []
        for k in keys:
            c = recipe.cate_to_comp(k)
            if c is None:
                unmapped.append(k)
        print(f"    unmapped N={len(unmapped)} example={unmapped[:15]}")


    # Keep a reference "Signal" MC HistHolder for later
    signal_mc_holder = mc_holders["Signal"]

    # -------------------- INIT SCALE HISTOGRAMS (ONE PER COMPONENT) --------------------
    # Use the fit_axis-defined ref_for_binning from above (EN4 for CCnue, Eel for nu+e)
    scale_hists = {}
    for comp in recipe.components:
        h = ref_for_binning.Clone(f"{comp}_Scale_Factor")
        h.Reset()
        h.GetYaxis().SetTitle("Scale Factor")
        h.SetTitle(f"{comp} Scale Factor")
        scale_hists[comp] = h

    # -------------------- RUN MATRIX FIT (CV + UNIVERSES) --------------------
    RunMinimizer_MatrixFit(recipe, data_holders, mc_holders, scale_hists)
    print("\n[DEBUG] scale hist CV ranges:")
    for comp, h in scale_hists.items():
        mn =  1e9
        mx = -1e9
        for b in range(0, h.GetNbinsX() + 2):
            v = h.GetBinContent(b)
            if v < mn: mn = v
            if v > mx: mx = v
        print(f"  {comp:>10s}: min={mn:.3g} max={mx:.3g}")

    print("\n[DEBUG] scale hist CV ranges:")
    for comp, h in scale_hists.items():
        mn =  1e9
        mx = -1e9
        for b in range(0, h.GetNbinsX() + 2):
            v = h.GetBinContent(b)
            if v < mn: mn = v
            if v > mx: mx = v
        print(f"  {comp:>10s}: min={mn:.3g} max={mx:.3g}")


    print("\n[DEBUG] grouped fractional uncertainties for scale factors")
    for comp, h in scale_hists.items():
        PrintGroupedErrorsForBin(
            h,
            CONSOLIDATED_ERROR_GROUPS,
            bins_to_check=[2, 3, 4, 5, 6],
            label=f"scale {comp}",
        )


    # -------------------- WRITE / PLOT SCALE HISTOGRAMS --------------------
    scalefile.cd()

    axis_tag = "Eel"

    for comp, hist in scale_hists.items():
        # Optional: keep a real physics axis label
        hist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())  # no-op; keeps whatever it has
        # or set explicitly:
        # hist.GetXaxis().SetTitle("E_{estimator} (GeV)" if axis_tag == "EN4" else "E_{e} (GeV)")

        c1 = ROOT.TCanvas()
        mnvplotter.DrawMCWithErrorBand(hist)
        PlotTools.Print(AnalysisConfig.PlotPath(f"{axis_tag}_scales", comp, background_fit_tag), mnvplotter, c1)

        c1 = ROOT.TCanvas()
        # mnvplotter.DrawErrorSummary(hist, "TR", True, True, 0)
        mnvplotter.DrawErrorSummaryDerived(hist, "TR", True, True, 0, False, "", True, "", True, "HIST")
        PlotTools.Print(AnalysisConfig.PlotPath(f"{axis_tag}_scale_errors", comp, background_fit_tag), mnvplotter, c1)

        hist.Write(f"{comp}_Scale_Factor")

    # -------------------- BUILD "PREDICTED SIGNAL" SHAPE (OPTIONAL) --------------------
    # Sum categories that map to the recipe's signal component (in the Signal region).
    mc_prediction = signal_mc_holder.GetHist().Clone("predicted_signal_total")
    mc_prediction.Reset()
    for cate, h in signal_mc_holder.hists.items():
        if h is None or cate == "Total":
            continue
        if cate in recipe.fixed_cates:
            continue
        if recipe.cate_to_comp(cate) == signal_comp_name:
            mc_prediction.Add(h)

    # -------------------- APPLY TUNE + BACKGROUND SUBTRACTION FOR UNFOLD INPUTS --------------------
    region_signal = "Signal"
    region_sideband_for_debug = "dEdX" if "dEdX" in recipe.regions else (recipe.regions[1] if len(recipe.regions) > 1 else None)

    for hist in HISTOGRAMS_TO_UNFOLD:
        # Debug: write post-fit breakdown in a representative sideband
        if region_sideband_for_debug is not None:
            data_hist_sb, mc_hist_sb, pred_hist_sb = GetScaledDataMC_Recipe(
                hist, datafile, mcfile, region_sideband_for_debug,
                recipe=recipe,
                scale_hists=scale_hists,
                pot_scale=pot_scale,
                scaled_hist_name=scaled_hist_name
            )
            for _h, hh in mc_hist_sb.hists.items():
                if hh:
                    hh.Write(f"POSTFIT_{region_sideband_for_debug}_{_h}")

        # Signal region: do bkg subtraction
        data_hist, mc_hist, pred_hist = GetScaledDataMC_Recipe(
            hist, datafile, mcfile, region_signal,
            recipe=recipe,
            scale_hists=scale_hists,
            pot_scale=pot_scale,
            scaled_hist_name=scaled_hist_name
        )

        mc_hist.GetHist().Write(data_hist.plot_name)

        subbedData, subbedMC = BackgroundSubtraction_Recipe(
            recipe, data_hist, mc_hist, pred_hist, signal_comp_name=signal_comp_name
        )
        subbedData.Write(data_hist.plot_name + "_data_bkgSubbed")
        mc_prediction.Write(data_hist.plot_name + "_predicted_Signal")

        # -------------------- PLOT BACKGROUND-SUBTRACTED RESULT --------------------
        c_sub = ROOT.TCanvas(f"c_sub_{data_hist.plot_name}", f"c_sub_{data_hist.plot_name}", 800, 600)

        subbedData.SetTitle("Background-subtracted data vs predicted signal")
        subbedData.GetYaxis().SetTitle("Events")
        subbedData.GetXaxis().SetTitle(data_hist.GetHist().GetXaxis().GetTitle())

        mnvplotter.DrawDataMCWithErrorBand(subbedData, subbedMC, 1.0, "TR")
        PlotTools.Print(
            AnalysisConfig.PlotPath(data_hist.plot_name + "_bkgSubtracted", region_signal, background_fit_tag),
            mnvplotter,
            c_sub,
        )

        # # Optional: data-only background-subtracted plot
        # c_sub_data = ROOT.TCanvas(f"c_sub_data_{data_hist.plot_name}", f"c_sub_data_{data_hist.plot_name}", 800, 600)

        # subbedData.SetTitle("Background-subtracted data")
        # subbedData.Draw("E1")
        # PlotTools.Print(
        #     AnalysisConfig.PlotPath(data_hist.plot_name + "_bkgSubtractedDataOnly", region_signal, background_fit_tag),
        #     mnvplotter,
        #     c_sub_data,
        # )



    # -------------------- (RE)OPEN FILES IF NEEDED FOR PLOTTING --------------------
    type_path_map = {
        t: AnalysisConfig.SelectionHistoPath(playlist, t == "data", False)
        for t in AnalysisConfig.data_types
    }
    datafile, mcfile, pot_scale = Utilities.getFilesAndPOTScale(
        playlist, type_path_map, AnalysisConfig.ntuple_tag
    )

    # -------------------- MAKE POSTFIT PLOTS --------------------
    for config in PLOTS_TO_MAKE:
        postfit_config = config.copy()
        postfit_config["tag"] = postfit_config.get("tag", "") + "postfit_"

        sideband_group = postfit_config.setdefault("sideband_group", recipe.regions)

        if isinstance(sideband_group, list):
            for sideband in sideband_group:
                data_hist, mc_hist, pred_hist = GetScaledDataMC_Recipe(
                    config["name"] if "name" in config else config,
                    datafile, mcfile, sideband,
                    recipe=recipe,
                    scale_hists=scale_hists,
                    pot_scale=pot_scale,
                    scaled_hist_name=scaled_hist_name
                )

                if sideband == "Signal" and AnalysisConfig.pseudodata:
                    # Use MC prediction as "data" for plotting (avoid unblinding)
                    MakePlot(data_holders["Signal"], pred_hist, postfit_config)
                else:
                    MakePlot(data_hist, pred_hist, postfit_config)

        else:
            # tuple mode: (name, [sidebands...])
            sideband = sideband_group[0]
            sidebands = sideband_group[1]

            data_hist, mc_hist, pred_hist = GetScaledDataMC_Recipe(
                config["name"] if "name" in config else config,
                datafile, mcfile, sidebands[0],
                recipe=recipe,
                scale_hists=scale_hists,
                pot_scale=pot_scale,
                scaled_hist_name=scaled_hist_name
            )

            for _ in range(1, len(sidebands)):
                dtmp, mtmp, ptmp = GetScaledDataMC_Recipe(
                    config["name"] if "name" in config else config,
                    datafile, mcfile, sidebands[_],
                    recipe=recipe,
                    scale_hists=scale_hists,
                    pot_scale=pot_scale,
                    scaled_hist_name=scaled_hist_name
                )
                data_hist.Add(dtmp)
                mc_hist.Add(mtmp)
                pred_hist.Add(ptmp)

            MakePlot(data_hist, pred_hist, postfit_config)

    # -------------------- CLOSE FILES --------------------
    datafile.Close()
    mcfile.Close()
    scalefile.Close()

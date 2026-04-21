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
# For recipe-driven matrix-fit, you should NOT rely on SIGNAL_DEFINATION inside fit/subtraction logic.
# Keep the import only if other plotting utilities still depend on it.
# from config.SignalDef import SIGNAL_DEFINATION

mnvplotter = PlotUtils.MnvPlotter()

from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS
mnvplotter.error_summary_group_map.clear()
for k, v in CONSOLIDATED_ERROR_GROUPS.items():
    vec = ROOT.vector("std::string")()
    for vs in v:
        vec.push_back(vs)
    mnvplotter.error_summary_group_map[k] = vec

# Stop ROOT from owning histograms (avoid segfaults / double-deletes)
ROOT.TH1.AddDirectory(False)



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

            # CV
            out.SetBinContent(i, sh.GetBinContent(k))
            out.SetBinError(i,   sh.GetBinError(k))

            # vertical error bands (universes)
            for bandname in sh.GetErrorBandNames():
                if bandname not in out.GetErrorBandNames():
                    # if clone didn't carry it for some reason, skip safely
                    continue

                out_band = out.GetVertErrorBand(bandname)
                src_band = sh.GetVertErrorBand(bandname)

                nh = min(out_band.GetNHists(), src_band.GetNHists())
                for u in range(nh):
                    out_u = out_band.GetHist(u)
                    src_u = src_band.GetHist(u)
                    out_u.SetBinContent(i, src_u.GetBinContent(k))
                    out_u.SetBinError(i,   src_u.GetBinError(k))

    return new_scale



# -----------------------------------------------------------------------------
# Helper: write scale histograms into output MnvH1D objects (CV / band / universe)
# -----------------------------------------------------------------------------
# def WriteScaleToMnvH1D(out_hists, in_hists, errorband=None, uni_idx=None):
#     """
#     Write CV if errorband=None.
#     Write one universe if errorband!=None and uni_idx is not None.
#     """
#     for key, hin in in_hists.items():
#         if key not in out_hists:
#             continue

#         hout = out_hists[key]

#         if errorband is None:
#             target = hout
#         else:
#             if uni_idx is None:
#                 # Don't allow writing into the band container
#                 raise RuntimeError(
#                     f"WriteScaleToMnvH1D called with errorband={errorband} but uni_idx=None. "
#                     "Write universes only: pass uni_idx."
#                 )
#             target = hout.GetVertErrorBand(errorband).GetHist(uni_idx)

#         for b in range(0, target.GetNbinsX() + 2):
#             target.SetBinContent(b, hin.GetBinContent(b))
#             target.SetBinError(b,   hin.GetBinError(b))
def WriteScaleToMnvH1D(out_hists, in_hists, errorband=None, uni_idx=None):
    for key, hin in in_hists.items():
        if key not in out_hists:
            continue

        hout = out_hists[key]

        if errorband is None:
            target = hout
        else:
            band = hout.GetVertErrorBand(errorband)
            if not band:
                continue
            if uni_idx is None:
                target = band
            else:
                target = band.GetHist(uni_idx)

        for b in range(0, target.GetNbinsX() + 2):
            target.SetBinContent(b, hin.GetBinContent(b))
            target.SetBinError(b, hin.GetBinError(b))


# -----------------------------------------------------------------------------
# LEGACY: 1D integral-based scale factor extraction (NOT used for matrix-fit)
# -----------------------------------------------------------------------------
# def Get1DScaleFactor(variable_hists, scale_hists):
#     """
#     Legacy behavior: compute an overall integral scale per group by applying
#     the scale hist bin-by-bin and taking the ratio of integrals.
#     This conflicts conceptually with matrix-fit where you already solved per-bin scales.
#     Keep only if you must support old 'CATEGORY_FACTORS' workflows.
#     """
#     scale_dict = {}
#     comparable_scale = MakeComparableMnvHXD(variable_hists.GetHist(), scale_hists, y_axis=False)
#     for cate in variable_hists.hists:
#         if variable_hists.hists[cate] is None:
#             continue
#         scaled = variable_hists.hists[cate].Clone()
#         try:
#             group = BackgroundFitConfig.CATEGORY_FACTORS[cate]
#             scale = comparable_scale[group]
#             scaled.Multiply(scaled, scale)
#             denom = variable_hists.hists[cate].Integral()
#             scale_dict[group] = (scaled.Integral() / denom) if denom != 0 else 0
#         except KeyError:
#             pass
#         del scaled
#     return scale_dict


# -------------------------------------------------------------------------
# NEW: matrix fit per bin (works for any #regions and #components)
# -------------------------------------------------------------------------

# def _get_region_hist(hist_holder, error_band=None, uni_idx=None):
#     """
#     Return the appropriate MnvH1D to use for a region given (CV / errorband / universe).
#     """
#     h = hist_holder.GetHist()
#     if error_band is None:
#         return h
#     if uni_idx is None:
#         return h.GetVertErrorBand(error_band)
#     return h.GetVertErrorBand(error_band).GetHist(uni_idx)
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


# def _solve_scales_per_bin(A_cols, b_vec, kreg=0.0, prior=None):
#     K = len(A_cols)
#     R = len(b_vec)

#     if prior is None:
#         prior = [1.0] * K

#     # Early exit: no MC content in any floated component in this bin
#     if all(sum(abs(A_cols[k][r]) for r in range(R)) == 0.0 for k in range(K)):
#         return prior

#     # Build normal equations
#     m = ROOT.TMatrixD(K, K)
#     y = ROOT.TVectorD(K)

#     for i in range(K):
#         # y_i = sum_r A_ir * b_r   (note: A_cols stores columns, so A_ir = A_cols[i][r])
#         s = 0.0
#         for r in range(R):
#             s += A_cols[i][r] * b_vec[r]
#         y[i] = s

#         for j in range(K):
#             sij = 0.0
#             for r in range(R):
#                 sij += A_cols[i][r] * A_cols[j][r]
#             if i == j:
#                 sij += float(kreg)
#             m[i][j] = sij

#     svd = ROOT.TDecompSVD(m)

#     # IMPORTANT: Solve(b) overwrites b with x and returns bool in PyROOT
#     xvec = ROOT.TVectorD(y)      # make a copy; will be overwritten into the solution
#     ok = svd.Solve(xvec)         # <- ok is bool; xvec becomes the solution if ok==True

#     if not ok:
#         return prior

#     x = [float(xvec[i]) for i in range(K)]
#     x = [max(0.0, xi) for xi in x]   # clamp physical

#     return x

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


def _solve_scales_all_bins_smooth_2comp(
    A_cols_by_bin,   # list over bins: [ [col_comp0, col_comp1], ... ], each col is list over regions
    b_vec_by_bin,    # list over bins: [ [b_r...], ... ]
    prior=(1.0, 1.0),
    lam_bkg=1.0,
    lam_sig=1.0,
    curvature=False,
    clamp_nonneg=True,
    fallback_det_tol=1e-12,
    use_original_binwise=False,
    use_prior=False,
    prior_bkg=0.2,
    prior_sig=0.5,
):
    """
    Solve for 2 component scale factors across bins.

    Modes
    -----
    1) use_original_binwise=True
       Solve each bin independently (the "original" fit).
       If use_prior=True, add a per-bin prior penalty toward `prior`.

    2) use_original_binwise=False
       Solve one global penalized least-squares problem across all bins.
       If curvature=True, smooth with second-difference penalty.
       If curvature=False, smooth with first-difference penalty.
       If use_prior=True, add a per-bin prior penalty toward `prior`.

    Inputs
    ------
    A_cols_by_bin[b] = [col_bkg, col_sig]
      where each column is a list over regions.
      So for bin b:
        col_bkg[r] = MC_bkg(region r, bin b)
        col_sig[r] = MC_sig(region r, bin b)

    b_vec_by_bin[b][r] = data(region r, bin b) - fixed(region r, bin b)

    Returns
    -------
    xb, xs
      lists of background and signal scale factors for each bin.
    """
    nbins = len(A_cols_by_bin)

    # ------------------------------------------------------------------
    # MODE 1: original per-bin analytic / least-squares solution
    # ------------------------------------------------------------------
    if use_original_binwise:
        xb = []
        xs = []

        for b in range(nbins):
            A_cols = A_cols_by_bin[b]
            d = b_vec_by_bin[b]

            # columns across regions
            a0 = A_cols[0]   # background component
            a1 = A_cols[1]   # signal component

            # Build normal equations for this bin:
            #   M x = y
            # where M = A^T A, y = A^T d
            m00 = 0.0
            m01 = 0.0
            m11 = 0.0
            y0  = 0.0
            y1  = 0.0

            for r in range(len(d)):
                m00 += a0[r] * a0[r]
                m01 += a0[r] * a1[r]
                m11 += a1[r] * a1[r]
                y0  += a0[r] * d[r]
                y1  += a1[r] * d[r]

            # Optional prior: add penalty
            #   prior_bkg * (sb - prior[0])^2 + prior_sig * (ss - prior[1])^2
            if use_prior:
                m00 += prior_bkg
                y0  += prior_bkg * prior[0]

                m11 += prior_sig
                y1  += prior_sig * prior[1]

            det = m00 * m11 - m01 * m01

            # Fallback if singular / nearly singular
            if abs(det) < fallback_det_tol:
                sb = prior[0]
                ss = prior[1]
            else:
                sb = ( y0 * m11 - y1 * m01) / det
                ss = (-y0 * m01 + y1 * m00) / det

            if clamp_nonneg:
                sb = max(0.0, sb)
                ss = max(0.0, ss)

            xb.append(sb)
            xs.append(ss)

        return xb, xs

    # ------------------------------------------------------------------
    # MODE 2: global smooth fit across all bins
    # ------------------------------------------------------------------
    npar = 2 * nbins  # [bkg bin0..N-1, sig bin0..N-1]

    H = ROOT.TMatrixD(npar, npar)  # normal matrix
    y = ROOT.TVectorD(npar)        # RHS

    def ibkg(b):
        return b

    def isig(b):
        return nbins + b

    # -----------------------------
    # data term: sum_b ||A_b x_b - d_b||^2
    # -----------------------------
    for b in range(nbins):
        A_cols = A_cols_by_bin[b]
        d = b_vec_by_bin[b]

        # A is R x 2, stored as columns
        a0 = A_cols[0]  # bkg
        a1 = A_cols[1]  # sig

        ata00 = 0.0
        ata01 = 0.0
        ata11 = 0.0
        atd0  = 0.0
        atd1  = 0.0

        for r in range(len(d)):
            ata00 += a0[r] * a0[r]
            ata01 += a0[r] * a1[r]
            ata11 += a1[r] * a1[r]
            atd0  += a0[r] * d[r]
            atd1  += a1[r] * d[r]

        i0 = ibkg(b)
        i1 = isig(b)

        H[i0][i0] += ata00
        H[i0][i1] += ata01
        H[i1][i0] += ata01
        H[i1][i1] += ata11

        y[i0] += atd0
        y[i1] += atd1

    # -----------------------------
    # smoothness term
    # -----------------------------
    def add_penalty(indices, coeffs, lam):
        # adds lam * (sum_k coeffs[k] * x[idx_k])^2
        for i in range(len(indices)):
            ii = indices[i]
            ci = coeffs[i]
            for j in range(len(indices)):
                jj = indices[j]
                cj = coeffs[j]
                H[ii][jj] += lam * ci * cj

    if curvature:
        # second differences: x_{b+1} - 2x_b + x_{b-1}
        for b in range(1, nbins - 1):
            add_penalty([ibkg(b-1), ibkg(b), ibkg(b+1)], [1.0, -2.0, 1.0], lam_bkg)
            add_penalty([isig(b-1), isig(b), isig(b+1)], [1.0, -2.0, 1.0], lam_sig)
    else:
        # first differences: x_{b+1} - x_b
        for b in range(nbins - 1):
            add_penalty([ibkg(b), ibkg(b+1)], [-1.0, 1.0], lam_bkg)
            add_penalty([isig(b), isig(b+1)], [-1.0, 1.0], lam_sig)

    # -----------------------------
    # optional prior anchor
    # -----------------------------
    if use_prior:
        for b in range(nbins):
            H[ibkg(b)][ibkg(b)] += prior_bkg
            y[ibkg(b)] += prior_bkg * prior[0]

            H[isig(b)][isig(b)] += prior_sig
            y[isig(b)] += prior_sig * prior[1]

    # Solve H x = y
    svd = ROOT.TDecompSVD(H)
    xvec = ROOT.TVectorD(y)
    ok = svd.Solve(xvec)

    if not ok:
        xb = [prior[0]] * nbins
        xs = [prior[1]] * nbins
    else:
        xb = [float(xvec[ibkg(b)]) for b in range(nbins)]
        xs = [float(xvec[isig(b)]) for b in range(nbins)]

    if clamp_nonneg:
        xb = [max(0.0, v) for v in xb]
        xs = [max(0.0, v) for v in xs]

    return xb, xs


def RunUniverseMinimizer_MatrixFit(recipe, data_holders, mc_holders, error_band=None, uni_idx=None, kreg=None):
    """
    Solve for scale factors per bin for all components in recipe, using all regions in recipe.

    data_holders: dict region -> HistHolder (data)
    mc_holders  : dict region -> HistHolder (MC)
    Returns: dict[comp] -> MnvH1D scale histogram (same binning as HIST_TO_FIT ref)
    """
    if kreg is None:
        kreg = BackgroundFitConfig.REGULATION_PARAMETER if hasattr(BackgroundFitConfig, "REGULATION_PARAMETER") else 0.0
    if recipe.kreg is not None:
        # recipe.kreg is "TSVD truncation" in your config; keep that name for later.
        # Here we interpret kreg as Tikhonov strength, not truncation.
        pass

    # Choose a reference hist for binning (use Signal region data)
    # href = _get_region_hist(data_holders[recipe.regions[0]], error_band, uni_idx)
    href = _get_region_hist(data_holders[recipe.regions[0]], use_universe=False)
    
    # define tag early so diagnostics inside loop can use it
    tag = "CV" if error_band is None else f"{error_band}[{uni_idx}]"

    # Output scale hists, one per component
    scales = {comp: href.Clone(f"scale_{comp}") for comp in recipe.components}
    for comp in scales:
        scales[comp].Reset()

    # --- Diagnostics setup ---
    diag = {
        "n_bins_total": 0,
        "n_bins_fallback": 0,
        "n_bins_dropped": 0,   # kept < full rank
        "drop_by_nkept": {},   # {n_kept: count}
        "fallback_bins": [],   # list of (bin, xcenter, total_mc)
        "dropped_bins": [],    # list of (bin, xcenter, n_kept, svals, total_mc)
        "worst_cond_bins": [], # list of (bin, xcenter, smin_kept/smax, n_kept)
    }

    # configurable threshold (events); pick 1.0 as you suggested
    min_total_mc = getattr(recipe, "min_total_mc", None)
    if min_total_mc is None:
        min_total_mc = getattr(BackgroundFitConfig, "MIN_TOTAL_FLOATED_MC", 1.0)

    # For TSVD keep count: interpret recipe.kreg as kmax if set, else keep all
    # kmax = getattr(recipe, "kreg", None)
    kmax_used = None
    tol = 0.02   # or getattr(recipe, "tol", None)

    # Build component sums + fixed sums for each region once
    comp_sums_by_region = {}
    fixed_by_region = {}
    data_by_region = {}
    for reg in recipe.regions:
        comp_sums_by_region[reg] = _build_component_sums(recipe, mc_holders[reg], error_band, uni_idx)
        fixed_by_region[reg] = _build_fixed_sum(recipe, mc_holders[reg], error_band, uni_idx)
        # data_by_region[reg] = _get_region_hist(data_holders[reg], error_band, uni_idx)
        data_by_region[reg] = _get_region_hist(data_holders[reg], use_universe=False)

    # Pseudodata option
    if AnalysisConfig.pseudodata:
        for reg in recipe.regions:
            # overwrite data with MC total CV (for this universe view)
            # safest: use the MC "Total" from holder
            mctot = mc_holders[reg].hists.get("Total", None)
            if mctot is None:
                continue
            mct = mctot
            if error_band is not None:
                if uni_idx is None:
                    mct = mctot.GetVertErrorBand(error_band)
                else:
                    mct = mctot.GetVertErrorBand(error_band).GetHist(uni_idx)

            for b in range(0, data_by_region[reg].GetNbinsX() + 2):
                data_by_region[reg].SetBinContent(b, mct.GetBinContent(b))
                data_by_region[reg].SetBinError(b,   mct.GetBinError(b))

    # Solve bin-by-bin
    nb = href.GetNbinsX()
    # skip_fit_bins = {1}
    for b in range(0, nb + 2):  # include under/overflow
        # if b in skip_fit_bins:
        #     continue
        # b_vec = data - fixed for each region
        b_vec = []
        for reg in recipe.regions:
            bval = data_by_region[reg].GetBinContent(b) - fixed_by_region[reg].GetBinContent(b)
            b_vec.append(bval)

        # A matrix columns (per component): A_cols[k][r] = MC_comp(reg,r,bin b)
        A_cols = []
        for comp in recipe.components:
            col = []
            for reg in recipe.regions:
                col.append(comp_sums_by_region[reg][comp].GetBinContent(b))
            A_cols.append(col)

        # OLD fitting not using regularized SVD
        # x = _solve_scales_per_bin(A_cols, b_vec, kreg=kreg)
        # Example settings:
        # tol=1e-2 keeps modes with s_i >= 1% of smax (often a good start)
        # kmax=recipe.kreg if you want recipe.kreg to mean "keep this many modes"
        # kmax = recipe.kreg if getattr(recipe, "kreg", None) is not None else None
        # x = _solve_scales_per_bin_tsvd(A_cols, b_vec, prior=[1.0]*len(recipe.components), tol=1e-2, kmax=kmax)
        # total floated MC in this bin (sum over regions and floated comps)
        total_mc = 0.0
        for ic in range(len(recipe.components)):
            for ir in range(len(recipe.regions)):
                total_mc += float(A_cols[ic][ir])

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
                A_cols, b_vec,
                prior=[1.0]*len(recipe.components),
                kmax=kmax_used,
                tol=tol,        
                clamp_nonneg=True,
                min_total_mc=None,
                want_diag=True
            )
        diag["n_bins_total"] += 1

        # bin x center (works for 1D)
        xcenter = href.GetXaxis().GetBinCenter(b) if 1 <= b <= href.GetNbinsX() else float("nan")

        if d["used_fallback"]:
            diag["n_bins_fallback"] += 1
            diag["fallback_bins"].append((b, xcenter, d["total_mc"]))

        else:
            # Print singular values if very truncated
            if d["n_kept"] <= 2:
                sv = ",".join(f"{s:.3g}" for s in d["svals"])
                print(f"[TSVD-DIAG {tag}] bin b={b} x={xcenter:.3g} kept={d['n_kept']} s=[{sv}] totalMC={d['total_mc']:.3g}")

            # full rank is min(R,K)
            full_rank = min(len(recipe.regions), len(recipe.components))
            if d["n_kept"] < full_rank:
                diag["n_bins_dropped"] += 1
                diag["drop_by_nkept"][d["n_kept"]] = diag["drop_by_nkept"].get(d["n_kept"], 0) + 1
                diag["dropped_bins"].append((b, xcenter, d["n_kept"], d["svals"], d["total_mc"]))

            # condition-ish metric: smin_kept/smax (smaller => more ill-conditioned)
            in_range = (1 <= b <= href.GetNbinsX())
            if in_range:
                ratio = (d["smin_kept"] / d["smax"]) if (d["smax"] > 0 and d["smin_kept"] > 0) else 0.0
                diag["worst_cond_bins"].append((b, xcenter, ratio, d["n_kept"]))


        for ic, comp in enumerate(recipe.components):
            scales[comp].SetBinContent(b, x[ic])
            scales[comp].SetBinError(b, 0.0)  # you can fill via universe spread later

    # --- Print diagnostics summary (CV or per-universe) ---
    nb = diag["n_bins_total"]

    reg_mode = []
    if kmax_used is not None:
        reg_mode.append(f"kmax={kmax_used}")
    if tol is not None:
        reg_mode.append(f"tol={tol}")
    reg_mode_str = " ".join(reg_mode) if reg_mode else "no_reg"

    print(f"[TSVD-DIAG {tag}] bins={nb} fallback={diag['n_bins_fallback']} dropped={diag['n_bins_dropped']} "
        f"min_total_mc={min_total_mc} {reg_mode_str}")

    if diag["drop_by_nkept"]:
        items = ", ".join([f"kept{nk}:{ct}" for nk, ct in sorted(diag["drop_by_nkept"].items())])
        print(f"[TSVD-DIAG {tag}] dropped breakdown: {items}")

    # Show worst-conditioned bins (lowest smin_kept/smax)
    diag["worst_cond_bins"].sort(key=lambda t: t[2])  # ascending ratio
    for (bb, xc, ratio, nk) in diag["worst_cond_bins"][:5]:
        if not math.isnan(xc):
            print(f"[TSVD-DIAG {tag}] worst bin b={bb} x={xc:.3g} smin/smax={ratio:.3g} n_kept={nk}")
        else:
            print(f"[TSVD-DIAG {tag}] worst bin b={bb} smin/smax={ratio:.3g} n_kept={nk}")

    # Show fallback bins (low MC)
    if diag["n_bins_fallback"] > 0:
        for (bb, xc, tmc) in diag["fallback_bins"][:10]:
            if not math.isnan(xc):
                print(f"[TSVD-DIAG {tag}] fallback bin b={bb} x={xc:.3g} totalMC={tmc:.3g}")
            else:
                print(f"[TSVD-DIAG {tag}] fallback bin b={bb} totalMC={tmc:.3g}")

    return scales





def Evaluate2CompChi2(
    A_cols_by_bin,
    b_vec_by_bin,
    xb,
    xs,
    lam_bkg=1.0,
    lam_sig=1.0,
    curvature=False,
):
    # -------------------------
    # data term
    # -------------------------
    chi2_data = 0.0
    for b in range(len(A_cols_by_bin)):
        a_bkg = A_cols_by_bin[b][0]   # list over regions
        a_sig = A_cols_by_bin[b][1]   # list over regions
        dvec  = b_vec_by_bin[b]       # data-fixed over regions

        for r in range(len(dvec)):
            pred = a_bkg[r] * xb[b] + a_sig[r] * xs[b]
            resid = pred - dvec[r]
            chi2_data += resid * resid

    # -------------------------
    # penalty term
    # -------------------------
    chi2_pen = 0.0

    if curvature:
        for b in range(1, len(xb) - 1):
            db = xb[b+1] - 2.0 * xb[b] + xb[b-1]
            ds = xs[b+1] - 2.0 * xs[b] + xs[b-1]
            chi2_pen += lam_bkg * db * db
            chi2_pen += lam_sig * ds * ds
    else:
        for b in range(len(xb) - 1):
            db = xb[b+1] - xb[b]
            ds = xs[b+1] - xs[b]
            chi2_pen += lam_bkg * db * db
            chi2_pen += lam_sig * ds * ds

    chi2_tot = chi2_data + chi2_pen
    return chi2_data, chi2_pen, chi2_tot

def RunUniverseMinimizer_MatrixFit_2CompSmooth(recipe, data_holders, mc_holders, error_band=None, uni_idx=None):
    """
    Special solver for 2-component recipes:
      - build all bins
      - solve globally across bins with smoothness penalty
    """
    if len(recipe.components) != 2:
        raise RuntimeError("RunUniverseMinimizer_MatrixFit_2CompSmooth requires exactly 2 components")

    # href = _get_region_hist(data_holders[recipe.regions[0]], error_band, uni_idx)
    href = _get_region_hist(data_holders[recipe.regions[0]], use_universe=False)
    tag = "CV" if error_band is None else f"{error_band}[{uni_idx}]"

    scales = {comp: href.Clone(f"scale_{comp}") for comp in recipe.components}
    for comp in scales:
        scales[comp].Reset()

    comp_sums_by_region = {}
    fixed_by_region = {}
    data_by_region = {}

    for reg in recipe.regions:
        comp_sums_by_region[reg] = _build_component_sums(recipe, mc_holders[reg], error_band, uni_idx)
        fixed_by_region[reg] = _build_fixed_sum(recipe, mc_holders[reg], error_band, uni_idx)
        # data_by_region[reg] = _get_region_hist(data_holders[reg], error_band, uni_idx)
        data_by_region[reg] = _get_region_hist(data_holders[reg], use_universe=False)

    if AnalysisConfig.pseudodata:
        for reg in recipe.regions:
            mctot = mc_holders[reg].hists.get("Total", None)
            if mctot is None:
                continue
            mct = mctot if error_band is None else mctot.GetVertErrorBand(error_band).GetHist(uni_idx)
            for b in range(0, data_by_region[reg].GetNbinsX() + 2):
                data_by_region[reg].SetBinContent(b, mct.GetBinContent(b))
                data_by_region[reg].SetBinError(b,   mct.GetBinError(b))

    # collect only real bins; keep under/overflow at 1 later
    A_cols_by_bin = []
    b_vec_by_bin = []
    real_bins = []

    skip_bins = set(getattr(recipe, "skip_fit_bins", [1]))

    min_total_mc = getattr(recipe, "min_total_mc", None)
    if min_total_mc is None:
        min_total_mc = getattr(BackgroundFitConfig, "MIN_TOTAL_FLOATED_MC", 1.0)

    # for b in range(1, href.GetNbinsX() + 1):
    for b in range(1, href.GetNbinsX() + 1):
        if b in skip_bins:
            continue
        
        b_vec = []
        for reg in recipe.regions:
            bval = data_by_region[reg].GetBinContent(b) - fixed_by_region[reg].GetBinContent(b)
            b_vec.append(bval)

        A_cols = []
        total_mc = 0.0
        for comp in recipe.components:
            col = []
            for reg in recipe.regions:
                v = comp_sums_by_region[reg][comp].GetBinContent(b)
                col.append(v)
                total_mc += float(v)
            A_cols.append(col)

        if total_mc < min_total_mc:
            # weak bin: still include, but anchor with near-prior by setting tiny A and d ~ prior
            # simplest is just use actual A,d and let smoothness handle it; if truly empty, use zeros
            pass

        A_cols_by_bin.append(A_cols)
        b_vec_by_bin.append(b_vec)
        real_bins.append(b)


    # curvature = getattr(recipe, "smooth_use_curvature", True)
    curvature = False

    use_original_binwise = getattr(recipe, "use_original_binwise", False)
    use_original_binwise = False

    use_prior = getattr(recipe, "use_prior", False)
    prior_bkg = getattr(recipe, "prior_bkg", 0.2)
    prior_sig = getattr(recipe, "prior_sig", 0.5)

    # ----------------------------------------
    # Lambda scan + plots
    # ----------------------------------------
    lambda_scan = [0.0, 0.1, 0.3, 0.6, 1.0, 3.0, 6.0, 10.0]

    print(f"[LAMBDA-SCAN {tag}] mode={'original' if use_original_binwise else ('curvature' if curvature else 'slope')} use_prior={use_prior}")

    lambda_vals = []
    chi2_data_vals = []
    chi2_pen_vals = []
    chi2_tot_vals = []
    max_bkg_vals = []
    max_sig_vals = []

    for lamb in lambda_scan:
        xb_scan, xs_scan = _solve_scales_all_bins_smooth_2comp(
            A_cols_by_bin,
            b_vec_by_bin,
            prior=(1.0, 1.0),
            lam_bkg=lamb,
            lam_sig=lamb,
            curvature=curvature,
            clamp_nonneg=True,
            use_original_binwise=use_original_binwise,
            use_prior=use_prior,
            prior_bkg=prior_bkg,
            prior_sig=prior_sig,
        )

        chi2_data_scan, chi2_pen_scan, chi2_tot_scan = Evaluate2CompChi2(
            A_cols_by_bin,
            b_vec_by_bin,
            xb_scan,
            xs_scan,
            lam_bkg=lamb,
            lam_sig=lamb,
            curvature=curvature,
        )

        lambda_vals.append(float(lamb))
        chi2_data_vals.append(float(chi2_data_scan))
        chi2_pen_vals.append(float(chi2_pen_scan))
        chi2_tot_vals.append(float(chi2_tot_scan))
        max_bkg_vals.append(float(max(xb_scan)))
        max_sig_vals.append(float(max(xs_scan)))

        print(
            f"[LAMBDA-SCAN {tag}] lam={lamb:.3g} "
            f"chi2_data={chi2_data_scan:.6g} "
            f"chi2_pen={chi2_pen_scan:.6g} "
            f"chi2_tot={chi2_tot_scan:.6g} "
            f"max_{recipe.components[0]}={max(xb_scan):.6g} "
            f"max_{recipe.components[1]}={max(xs_scan):.6g}"
        )

    # --- chi2 graphs ---
    g_data = ROOT.TGraph(len(lambda_vals), array('d', lambda_vals), array('d', chi2_data_vals))
    g_pen  = ROOT.TGraph(len(lambda_vals), array('d', lambda_vals), array('d', chi2_pen_vals))
    g_tot  = ROOT.TGraph(len(lambda_vals), array('d', lambda_vals), array('d', chi2_tot_vals))

    g_pen.SetTitle(f"Lambda scan ({tag});#lambda;#chi^{{2}}")
    g_pen.SetMinimum(1e-4)
    g_data.SetLineWidth(2)
    g_pen.SetLineWidth(2)
    g_tot.SetLineWidth(2)

    g_data.SetLineColor(ROOT.kBlue + 1)
    g_pen.SetLineColor(ROOT.kRed + 1)
    g_tot.SetLineColor(ROOT.kBlack)

    c_lambda = ROOT.TCanvas(f"c_lambda_{tag}", f"c_lambda_{tag}", 800, 600)
    g_pen.Draw("ALP")
    g_data.Draw("LP SAME")
    g_tot.Draw("LP SAME")

    leg = ROOT.TLegend(0.58, 0.28, 0.88, 0.58)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.AddEntry(g_data, "#chi^{2}_{data}", "lp")
    leg.AddEntry(g_pen,  "#chi^{2}_{pen}", "lp")
    leg.AddEntry(g_tot,  "#chi^{2}_{tot}", "lp")
    leg.Draw()

    c_lambda.SetLogy(0)
    c_lambda.Print(f"lambda_scan_chi2_newBinning.png")
    c_lambda.SetLogy(1)
    c_lambda.Print(f"lambda_scan_chi2_newBinning_logy.png")


    # # --- max scale graphs ---
    # g_max_bkg = ROOT.TGraph(len(lambda_vals), array('d', lambda_vals), array('d', max_bkg_vals))
    # g_max_sig = ROOT.TGraph(len(lambda_vals), array('d', lambda_vals), array('d', max_sig_vals))

    # g_max_bkg.SetTitle(f"Lambda scan max scales ({tag});#lambda;max scale factor")
    # g_max_bkg.SetLineWidth(2)
    # g_max_sig.SetLineWidth(2)

    # g_max_bkg.SetLineColor(ROOT.kGreen + 2)
    # g_max_sig.SetLineColor(ROOT.kMagenta + 1)

    # c_max = ROOT.TCanvas(f"c_max_{tag}", f"c_max_{tag}", 800, 600)
    # g_max_bkg.Draw("ALP")
    # g_max_sig.Draw("LP SAME")

    # leg2 = ROOT.TLegend(0.28, 0.75, 0.58, 0.88)
    # leg2.SetFillStyle(0)
    # leg2.SetBorderSize(0)
    # leg2.AddEntry(g_max_bkg, f"max {recipe.components[0]}", "lp")
    # leg2.AddEntry(g_max_sig, f"max {recipe.components[1]}", "lp")
    # leg2.Draw()

    # c_max.Print(f"lambda_scan_maxscale_{tag}.png")




    # -------------------------------------------------
    # Actual fit with the chosen lambda values
    # -------------------------------------------------
    lam_bkg = getattr(recipe, "smooth_lambda_bkg", 1.0)
    lam_sig = getattr(recipe, "smooth_lambda_sig", 1.0)
    xb, xs = _solve_scales_all_bins_smooth_2comp(
        A_cols_by_bin,
        b_vec_by_bin,
        prior=(1.0, 1.0),
        lam_bkg=lam_bkg,
        lam_sig=lam_sig,
        curvature=curvature,
        clamp_nonneg=True,
        use_original_binwise=use_original_binwise,
        use_prior=use_prior,
        prior_bkg=prior_bkg,
        prior_sig=prior_sig,
    )

    chi2_data, chi2_pen, chi2_tot = Evaluate2CompChi2(
        A_cols_by_bin,
        b_vec_by_bin,
        xb,
        xs,
        lam_bkg=lam_bkg,
        lam_sig=lam_sig,
        curvature=curvature,
    )

    comp0, comp1 = recipe.components[0], recipe.components[1]

    # default under/overflow = 1
    for b in range(0, href.GetNbinsX() + 2):
        scales[comp0].SetBinContent(b, 1.0)
        scales[comp1].SetBinContent(b, 1.0)
        scales[comp0].SetBinError(b, 0.0)
        scales[comp1].SetBinError(b, 0.0)

    for i, b in enumerate(real_bins):
        scales[comp0].SetBinContent(b, xb[i])
        scales[comp1].SetBinContent(b, xs[i])

    mode = "original" if use_original_binwise else ("curvature" if curvature else "slope")

    if use_original_binwise:
        print(
            f"[SMOOTH-2COMP {tag}] mode={mode} use_prior={use_prior} "
            f"prior_bkg={prior_bkg} prior_sig={prior_sig}"
        )
    else:
        print(
            f"[SMOOTH-2COMP {tag}] mode={mode} use_prior={use_prior} "
            f"lam_bkg={lam_bkg} lam_sig={lam_sig} "
            f"prior_bkg={prior_bkg} prior_sig={prior_sig}"
        )

    print(f"[SMOOTH-2COMP {tag}] chi2_data={chi2_data:.6g} chi2_pen={chi2_pen:.6g} chi2_tot={chi2_tot:.6g}")
    print(f"[SMOOTH-2COMP {tag}] {comp0}: min={min(xb):.3g} max={max(xb):.3g}")
    print(f"[SMOOTH-2COMP {tag}] {comp1}: min={min(xs):.3g} max={max(xs):.3g}")

    # ---- DEBUG: print bin-by-bin scales ----
    for comp in recipe.components:
        h = scales[comp]
        vals = [(b, h.GetXaxis().GetBinCenter(b), h.GetBinContent(b))
                for b in range(1, h.GetNbinsX()+1)]

        print(f"[SMOOTH-2COMP {tag}] {comp} bins:")
        for bb, xx, vv in vals:
            print(f"  b={bb:2d} x={xx:6.3f} scale={vv:7.3f}")
    # ----------------------------------------

    return scales





# def RunMinimizer_MatrixFit(recipe, data_holders, mc_holders, scale_hists_out):
#     """
#     Fill scale_hists_out (dict[comp] -> MnvH1D with error bands) for:
#       - CV
#       - every vertical error-band universe
#     """
#     # CV
#     cv_scales = RunUniverseMinimizer_MatrixFit(recipe, data_holders, mc_holders, error_band=None, uni_idx=None)
#     WriteScaleToMnvH1D(scale_hists_out, cv_scales, errorband=None, uni_idx=None)
#     print("[MatrixFit] Done with CV scales")

#     any_hist = mc_holders[recipe.regions[0]].GetHist()
#     for error_band in any_hist.GetErrorBandNames():
#         nuni = any_hist.GetVertErrorBand(error_band).GetNHists()
#         for i in range(nuni):
#             uni_scales = RunUniverseMinimizer_MatrixFit(recipe, data_holders, mc_holders, error_band=error_band, uni_idx=i)
#             WriteScaleToMnvH1D(scale_hists_out, uni_scales, errorband=error_band, uni_idx=i)

#         print(f"[MatrixFit] Done with universes for band {error_band} (N={nuni})")

# def RunMinimizer_MatrixFit(recipe, data_holders, mc_holders, scale_hists_out):
#     """
#     Fill scale_hists_out (dict[comp] -> MnvH1D with error bands) for:
#       - CV
#       - every vertical error-band universe
#     """
#     use_smooth_2comp = (len(recipe.components) == 2)

#     if use_smooth_2comp:
#         cv_scales = RunUniverseMinimizer_MatrixFit_2CompSmooth(
#             recipe, data_holders, mc_holders, error_band=None, uni_idx=None
#         )
#     else:
#         cv_scales = RunUniverseMinimizer_MatrixFit(
#             recipe, data_holders, mc_holders, error_band=None, uni_idx=None
#         )

#     WriteScaleToMnvH1D(scale_hists_out, cv_scales, errorband=None, uni_idx=None)
#     print("[MatrixFit] Done with CV scales")

#     any_hist = mc_holders[recipe.regions[0]].GetHist()
#     for error_band in any_hist.GetErrorBandNames():
#         nuni = any_hist.GetVertErrorBand(error_band).GetNHists()
#         for i in range(nuni):
#             if use_smooth_2comp:
#                 uni_scales = RunUniverseMinimizer_MatrixFit_2CompSmooth(
#                     recipe, data_holders, mc_holders, error_band=error_band, uni_idx=i
#                 )
#             else:
#                 uni_scales = RunUniverseMinimizer_MatrixFit(
#                     recipe, data_holders, mc_holders, error_band=error_band, uni_idx=i
#                 )
#             WriteScaleToMnvH1D(scale_hists_out, uni_scales, errorband=error_band, uni_idx=i)

#         print(f"[MatrixFit] Done with universes for band {error_band} (N={nuni})")
def RunMinimizer_MatrixFit(recipe, data_holders, mc_holders, scale_hists_out):
    use_smooth_2comp = (len(recipe.components) == 2)

    # CV
    if use_smooth_2comp:
        cv_scales = RunUniverseMinimizer_MatrixFit_2CompSmooth(
            recipe, data_holders, mc_holders, error_band=None, uni_idx=None
        )
    else:
        cv_scales = RunUniverseMinimizer_MatrixFit(
            recipe, data_holders, mc_holders, error_band=None, uni_idx=None
        )

    WriteScaleToMnvH1D(scale_hists_out, cv_scales)
    print("[MatrixFit] Done with CV scales")

    any_hist = mc_holders[recipe.regions[0]].GetHist()

    for error_band in any_hist.GetErrorBandNames():
        # fill band CV first, unsmoothed
        band_scales = RunUniverseMinimizer_MatrixFit(
            recipe, data_holders, mc_holders, error_band=error_band, uni_idx=None
        )
        WriteScaleToMnvH1D(scale_hists_out, band_scales, errorband=error_band, uni_idx=None)

        # then fill universes, also unsmoothed
        nuni = any_hist.GetVertErrorBand(error_band).GetNHists()
        for i in range(nuni):
            uni_scales = RunUniverseMinimizer_MatrixFit(
                recipe, data_holders, mc_holders, error_band=error_band, uni_idx=i
            )
            WriteScaleToMnvH1D(scale_hists_out, uni_scales, errorband=error_band, uni_idx=i)

        print(f"[MatrixFit] Done with universes for band {error_band} (N={nuni})")


# -------------------------------------------------------------------------
# NEW: apply scales to MC holders using recipe mapping
# -------------------------------------------------------------------------
def ScaleCategories_Recipe(recipe, hist_holder, scale_hists_mapped, prediction=False, signal_comp_name="Signal"):
    """
    Multiply each category histogram by the appropriate scale histogram, based on recipe.cate_to_comp.
    prediction=False: scale all floated comps (including signal comp)
    prediction=True : scale only background comps (leave signal comp untouched)
    """
    for cate, h in hist_holder.hists.items():
        if h is None or cate == "Total":
            continue
        if cate in recipe.fixed_cates:
            continue

        comp = recipe.cate_to_comp(cate)
        if comp is None:
            continue
        if comp not in scale_hists_mapped:
            continue

        if prediction and comp == signal_comp_name:
            continue

        try:
            h.Multiply(h, scale_hists_mapped[comp])
        except Exception:
            print(f"[MatrixFit] Multiply failed for cate={cate} region={hist_holder.sideband} comp={comp}")
            continue


def TuneMC_Recipe(recipe, hist_holder, scale_hists, x_axis=False, y_axis=False, prediction=False, signal_comp_name="Signal"):
    """
    Apply fitted scales to hist_holder.
    If x_axis/y_axis True: first map 1D scale factors onto the hist's axis bins (needed for 2D).
    """
    if x_axis and y_axis:
        return False

    if x_axis or y_axis:
        mapped = MakeComparableMnvHXD(hist_holder.GetHist(), scale_hists, y_axis=y_axis)
        ScaleCategories_Recipe(recipe, hist_holder, mapped, prediction=prediction, signal_comp_name=signal_comp_name)
    else:
        # 1D hist: scales already comparable
        ScaleCategories_Recipe(recipe, hist_holder, scale_hists, prediction=prediction, signal_comp_name=signal_comp_name)

    hist_holder.ResumTotal()
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


def GetScaledDataMC(hist, datafile, mcfile, region, recipe, scale_hists, pot_scale=None, scaled_hist_name=None):
    return GetScaledDataMC_Recipe(
        hist, datafile, mcfile, region,
        recipe=recipe,
        scale_hists=scale_hists,
        pot_scale=pot_scale if pot_scale is not None else 1.0,
        scaled_hist_name=scaled_hist_name,
    )

def BackgroundSubtraction(data_hists, tuned_mc_hists, pred_hists, recipe, signal_comp_name):
    return BackgroundSubtraction_Recipe(
        recipe, data_hists, tuned_mc_hists, pred_hists,
        signal_comp_name=signal_comp_name
    )


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


# -------------------------------------------------------------------------
# LEGACY FUNCTIONS (comment out once matrix-fit is fully used)
# -------------------------------------------------------------------------
# def RunUniverseMinimizer(...): pass
# def RunMinimizer(...): pass
# def TuneMC(...): pass
# def ScaleCategories(...): pass
# def ScaleCategories1D(...): pass
# def BackgroundSubtraction(...): pass
# def GetScaledDataMC(...): pass

def GetSignalCompName(recipe):
    if "Signal" in recipe.components:
        return "Signal"
    if "SignalLike" in recipe.components:
        return "SignalLike"
    raise RuntimeError(f"Cannot determine signal component name from {recipe.components}")


def MakeRatio_Recipe(
    recipe,
    signalHist, sidebandHist,
    normsignalHist, normsidebandHist,
    config,
    signal_comp_name="Signal",
    include_fixed_in_bkg=True,
    flux_band_name="Flux",
):
    # scale
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

    def _sum_background(holder, name):
        out = holder.hists["Total"].Clone(name)
        out.Reset()

        for cate, h in holder.hists.items():
            if h is None or cate == "Total":
                continue

            if cate in recipe.fixed_cates:
                if include_fixed_in_bkg:
                    out.Add(h)
                continue

            comp = recipe.cate_to_comp(cate)
            if comp is None:
                continue

            if comp != signal_comp_name:
                out.Add(h)

        return out

    sig_bkg     = _sum_background(signalHist,     "sig_bkg")
    sid_bkg     = _sum_background(sidebandHist,   "sid_bkg")
    normsig_bkg = _sum_background(normsignalHist, "normsig_bkg")
    normsid_bkg = _sum_background(normsidebandHist,"normsid_bkg")

    # Optional debug: Flux band checks
    def _has_band(h, band):
        try:
            return band in list(h.GetErrorBandNames())
        except Exception:
            return False

    if _has_band(sig_bkg, flux_band_name) and _has_band(sid_bkg, flux_band_name):
        c1 = ROOT.TCanvas()
        sig_bkg.GetVertErrorBand(flux_band_name).DrawAll("hist", True)
        c1.Print(f"{signalHist.plot_name}_post_tuneSigBkg{flux_band_name}.png")

        c1 = ROOT.TCanvas()
        sid_bkg.GetVertErrorBand(flux_band_name).DrawAll("hist", True)
        c1.Print(f"{signalHist.plot_name}_post_tuneSidBkg{flux_band_name}.png")

        # Compare flux error bands pre/post tune (normalized ratio)
        c1 = ROOT.TCanvas()
        fluxerr = sig_bkg.GetVertErrorBand(flux_band_name).GetErrorBand(True, False).Clone()
        normfluxerr = normsig_bkg.GetVertErrorBand(flux_band_name).GetErrorBand(True, False).Clone()
        normfluxerr.Divide(normfluxerr, fluxerr)
        normfluxerr.Draw()
        c1.Print(f"{signalHist.plot_name}_sig_bkgPostTune_{flux_band_name}errband.png")

        c1 = ROOT.TCanvas()
        fluxerr = sid_bkg.GetVertErrorBand(flux_band_name).GetErrorBand(True, False).Clone()
        normfluxerr = normsid_bkg.GetVertErrorBand(flux_band_name).GetErrorBand(True, False).Clone()
        normfluxerr.Divide(normfluxerr, fluxerr)
        normfluxerr.Draw()
        c1.Print(f"{signalHist.plot_name}_sid_bkgPostTune_{flux_band_name}errband.png")

    # Ratio of normalized background shapes: Signal-region-bkg / Sideband-region-bkg
    isig = sig_bkg.Integral()
    isd  = sid_bkg.Integral()
    if isig <= 0 or isd <= 0:
        print(f"[MakeRatio_Recipe] Skip ratio: zero integral (sig={isig}, sid={isd})")
        return

    sig_bkg.Scale(1.0 / isig)
    sid_bkg.Scale(1.0 / isd)

    sig_bkg.Divide(sig_bkg, sid_bkg)
    sig_bkg.GetXaxis().SetTitle("E_{available} + E_{lepton}")
    sig_bkg.GetYaxis().SetTitle("Ratio")
    sig_bkg.SetTitle("Signal/Sideband Background Ratio")

    mnvplotter.DrawMCWithErrorBand(sig_bkg)
    PlotTools.Print(AnalysisConfig.PlotPath("EN4_ratio", "Combined", "bkgratio_matrixfit"))

    c1 = ROOT.TCanvas()
    mnvplotter.DrawErrorSummary(sig_bkg, "TR", True, True, 0)
    c1.Print(f"{signalHist.plot_name}_post_tuneSigBkgErrSummary.png")


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


    # -------------------- WRITE / PLOT SCALE HISTOGRAMS --------------------
    scalefile.cd()

    axis_tag = "EN4" if ("Neutrino" in fit_axis or "EN4" in fit_axis) else "Eel"

    for comp, hist in scale_hists.items():
        # Optional: keep a real physics axis label
        hist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())  # no-op; keeps whatever it has
        # or set explicitly:
        # hist.GetXaxis().SetTitle("E_{estimator} (GeV)" if axis_tag == "EN4" else "E_{e} (GeV)")

        c1 = ROOT.TCanvas()
        mnvplotter.DrawMCWithErrorBand(hist)
        PlotTools.Print(AnalysisConfig.PlotPath(f"{axis_tag}_scales", comp, background_fit_tag), mnvplotter, c1)

        c1 = ROOT.TCanvas()
        mnvplotter.DrawErrorSummary(hist, "TR", True, True, 0)
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

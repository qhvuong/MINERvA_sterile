#!/usr/bin/env python3
"""Fit with a normalized P8 flux prior, then save a fitData-style spectrum/ratio plot.

Install in oscillations/ and run in the usual ROOT/PlotUtils environment.
Uses the production OscillationFitter.DoFit search; --k selects retained modes.
"""
import argparse
import csv
import math
import os
from pathlib import Path
import sys

import numpy as np


class NormalizedFluxPrior:
    """Profile the Gaussian pulls analytically in non-flux-whitened bin space."""

    def __init__(self, nonflux, shifts, scale=1.0):
        if not math.isfinite(scale) or scale <= 0:
            raise ValueError("Flux scale must be positive and finite")
        self.n_universes, self.n_bins = shifts.shape
        self.L = np.linalg.cholesky(nonflux)
        centered = shifts - shifts.mean(axis=0)
        response = centered.T / np.sqrt(self.n_universes)
        white = np.linalg.solve(self.L, response)
        U, s, _ = np.linalg.svd(white, full_matrices=False)
        keep = s > np.finfo(s.dtype).eps * max(white.shape) * s[0]
        self.U, self.s = U[:, keep], scale * s[keep]
        self.K = int(keep.sum())
        self.eigenvalues = self.s**2
        self.flux_cov = response @ response.T

    def fit_whitened(self, y):
        z = self.U.T @ y
        pulls = self.s * z / (1.0 + self.eigenvalues)
        flux_white = self.U @ (self.s * pulls)
        after = y - flux_white
        return float(after @ after), float(pulls @ pulls), flux_white

    def evaluate(self, residual):
        y = np.linalg.solve(self.L, residual)
        residual_chi2, penalty, flux_white = self.fit_whitened(y)
        return residual_chi2 + penalty, residual_chi2, penalty, self.L @ flux_white

    def assert_flux_covariance(self, target):
        diff = np.linalg.solve(self.L, target - self.flux_cov)
        diff = np.linalg.solve(self.L, diff.T).T
        np.testing.assert_allclose(diff, 0.0, atol=1e-8, rtol=0,
                                   err_msg="Normalized universe covariance does not close to ROOT")


class RetainedModeFluxPrior:
    """Normalized modes; optionally solve pulls on a subset and score all bins."""

    def __init__(self, nonflux, shifts, k=None, profile_bins=None, scale=1.0):
        from scipy.linalg import solve_triangular
        self._solve = solve_triangular
        if not math.isfinite(scale) or scale <= 0:
            raise ValueError("Flux scale must be positive and finite")
        self.n_universes, self.n_bins = shifts.shape
        self.L = np.linalg.cholesky(nonflux)
        self.profile_bins = (np.arange(self.n_bins) if profile_bins is None
                             else np.asarray(profile_bins, dtype=int))
        if not len(self.profile_bins):
            raise ValueError("No bins selected for profiling")
        self.L_profile = np.linalg.cholesky(nonflux[np.ix_(self.profile_bins, self.profile_bins)])
        response = (shifts - shifts.mean(axis=0)).T / np.sqrt(self.n_universes)
        white = solve_triangular(self.L_profile, response[self.profile_bins], lower=True)
        U, s, Vh = np.linalg.svd(white, full_matrices=False)
        threshold = np.finfo(s.dtype).eps * max(white.shape) * s[0]
        self.full_K = int(np.count_nonzero(s > threshold))
        self.K = self.full_K if k is None else k
        if not 0 <= self.K <= self.full_K:
            raise ValueError(f"K={self.K} exceeds numerical profiling rank {self.full_K}")
        # Map the same universe-space directions into every scored bin.
        self.modes = scale * response @ Vh[:self.K].T
        ss = scale * s[:self.K]
        self.pull_map = (ss / (1. + ss**2))[:, None] * U[:, :self.K].T
        self.retained_cov = self.modes @ self.modes.T
        self.flux_cov = response @ response.T

    def evaluate(self, residual):
        y = self._solve(self.L_profile, residual[self.profile_bins], lower=True, check_finite=False)
        pulls = self.pull_map @ y
        shift = self.modes @ pulls
        after = self._solve(self.L, residual-shift, lower=True, check_finite=False)
        residual_chi2, penalty = float(after @ after), float(pulls @ pulls)
        return residual_chi2+penalty, residual_chi2, penalty, shift


def vector(hist):
    return np.array([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX()+1)])


def fit_plot_values(*args, **kwargs):
    """Compatibility wrapper; shared implementation lives in tools."""
    from tools.PlotNormalizedFluxPrior import fit_plot_values as evaluate_plot
    return evaluate_plot(*args, **kwargs)


def plot_fit_spectra(*args, **kwargs):
    """Compatibility entry point also used by the mode-convergence study."""
    from tools.PlotNormalizedFluxPrior import PlotNormalizedFluxPrior
    return PlotNormalizedFluxPrior(*args, **kwargs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hist-config-tag", required=True)
    parser.add_argument("--k", type=int, help="Retained normalized modes; default: full available rank")
    parser.add_argument("--exclude-ratios-from-profiling", action="store_true",
                        help="ratios_dropped: solve pulls without ratio bins, but score all bins")
    parser.add_argument("--flux-source", choices=("residual", "mc"), default="residual",
                        help="residual: match full ROOT covariance; mc: earlier study's MC-only prior")
    parser.add_argument("--flux-scale", type=float, default=1.0,
                        help="Fixed multiplier on flux standard deviations (default: 1)")
    parser.add_argument("--spectra-out", type=Path,
                        help="Optional CSV for null/best spectra and fitted residual corrections")
    parser.add_argument("--plots-dir", type=Path,
                        help="Plot output directory (default: oscillations/plots)")
    parser.add_argument("--ntuple-tag", default="MAD", help="Label/prefix for post-fit plots")
    parser.add_argument("--no-plots", action="store_true", help="Skip automatic post-fit plotting")
    args = parser.parse_args()
    if args.k is not None and args.k < 0:
        parser.error("--k must be nonnegative")
    if not math.isfinite(args.flux_scale) or args.flux_scale <= 0:
        parser.error("--flux-scale must be positive and finite")
    osc = Path(__file__).resolve().parent
    tag = args.hist_config_tag
    cfg = osc / f"HIST_CONFIG_{tag}.json"
    rootfile = osc / "rootfiles" / f"NuE_stitched_hists_{tag}.root"
    for path in (cfg, rootfile):
        if not path.is_file():
            parser.error(f"Missing {path}")

    sys.path.insert(0, str(osc.parent))
    os.environ.setdefault("CCNUEROOT", str(osc.parent))
    sys.argv = [sys.argv[0]]  # do not pass these flags to the legacy config parser
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    from scipy.linalg import cho_factor, cho_solve
    from tools.Fitters import OscillationFitter
    from tools.Helper import TMatrix_to_Numpy
    from tools.StitchedHistogram import StitchedHistogram

    h = StitchedHistogram("sample", include_shared_ratio_stats=True)
    h.SetHistConfig(str(cfg))
    h.Load(str(rootfile))
    mh, dh = h.GetMCHistogram(), h.GetDataHistogram()
    d, m = vector(dh), vector(mh)
    C = h.GetCovarianceMatrix(True)
    root_full = h.GetCovarianceMatrix(False)
    A = np.asarray(h.GetAMatrix(), dtype=float)
    if A.shape != (1000, len(d)):
        raise ValueError(f"Expected 1000 P8 shifts across {len(d)} bins, got {A.shape}")

    if args.flux_source == "residual":
        band = dh.GetVertErrorBand("Flux")
        if not band or band.GetNHists() != len(A):
            raise ValueError("Data Flux band does not match the MC Flux band")
        A = np.array([vector(band.GetHist(i)) - d for i in range(len(A))]) - A
        target_flux = root_full - C
    else:
        target_flux = TMatrix_to_Numpy(mh.GetSysErrorMatrix("Flux"))[1:-1, 1:-1]

    scale = args.flux_scale
    # Check the complete 1000-universe prior before any truncation.
    full_prior = NormalizedFluxPrior(C, A, scale=scale)
    full_prior.assert_flux_covariance(target_flux)
    profile_bins = np.arange(len(d))
    if args.exclude_ratios_from_profiling:
        ratio_mask = vector(h.ratio_id)
        if not np.isin(ratio_mask, [0., 1.]).all() or not np.any(ratio_mask == 1.):
            parser.error("This configuration has no valid ratio-bin mask")
        profile_bins = profile_bins[ratio_mask == 0.]
    prior = RetainedModeFluxPrior(C, A, args.k, profile_bins, scale)
    all_profiled = len(profile_bins) == len(d)
    # This covariance identity applies when the same bins solve and score.
    matched = cho_factor(C + prior.retained_cov) if all_profiled else None
    root_factor = cho_factor(root_full)

    def prediction(dm2, ue, umu, utau):
        h.OscillateHistogram(dm2, ue, umu, utau, False, False)
        return vector(h.GetOscillatedHistogram())

    print(f"Flux source: {args.flux_source}; universes={prior.n_universes}; "
          f"profile bins={len(profile_bins)}; scored bins={prior.n_bins}; "
          f"retained K={prior.K}/{prior.full_K}; flux scale={scale:g}; "
          f"lambda relative to unscaled P8={1.0/scale**2:g}")
    print("Unscaled normalized covariance closes to the selected ROOT flux covariance.")

    def report(label, pred):
        r = d-pred
        score, residual, penalty, _ = prior.evaluate(r)
        cov_score = float(r @ cho_solve(matched, r)) if matched is not None else None
        root_score = float(r @ cho_solve(root_factor, r))
        if matched is not None:
            np.testing.assert_allclose(score, cov_score, rtol=1e-9, atol=1e-7)
        if all_profiled and prior.K == prior.full_K and args.flux_source == "residual" and scale == 1:
            np.testing.assert_allclose(score, root_score, rtol=1e-9, atol=1e-7)
        print(f"{label}: profile chi2={score:.8f}, residual={residual:.8f}, "
              f"pull penalty={penalty:.8f}, matched covariance chi2={cov_score}, "
              f"unscaled full ROOT chi2={root_score:.8f}")
        return score

    class NormalizedFitter(OscillationFitter):
        def CalChi2(self, x):
            pred = prediction(100.0*x[0], *x[1:])
            return prior.evaluate(d-pred)[0]

    np.testing.assert_allclose(prediction(0., 0., 0., 0.), m, rtol=1e-7, atol=1e-8)
    null = report("Null", m)
    fitter = NormalizedFitter(h, exclude="", marginalize_flux=False, hist_config=str(cfg))
    # CalChi2 above supplies the normalized profiling objective for every
    # optimizer step; DoFit only provides the production multistart search.
    fit_score, params = fitter.DoFit()
    point = [params[key] for key in ("m", "ue4", "umu4", "utau4")]
    best_pred = prediction(*point)
    best = report("Best fit", best_pred)
    np.testing.assert_allclose(best, fit_score, rtol=1e-9, atol=1e-7)
    print(f"Delta chi2 = {null-best:.8f}")
    print(f"dm2={point[0]:.8g} eV^2, |Ue4|^2={point[1]:.8g}, "
          f"|Umu4|^2={point[2]:.8g}, |Utau4|^2={point[3]:.8g}")
    plots_dir = args.plots_dir if args.plots_dir is not None else osc / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    mode_tag = "ratios_dropped" if args.exclude_ratios_from_profiling else "all_bins"
    output_stem = f"{tag}_{mode_tag}_K{prior.K}_s{scale:g}_{args.flux_source}"
    spectra_out = args.spectra_out or plots_dir / f"{output_stem}_spectra.csv"
    if spectra_out:
        spectra_out.parent.mkdir(parents=True, exist_ok=True)
        null_shift = prior.evaluate(d-m)[3]
        best_shift = prior.evaluate(d-best_pred)[3]
        with spectra_out.open("w", newline="") as output:
            writer = csv.writer(output)
            writer.writerow(("bin", "data", "null_cv", "best_cv", "null_flux_correction",
                             "best_flux_correction", "null_profiled", "best_profiled"))
            for i, values in enumerate(zip(d, m, best_pred, null_shift, best_shift,
                                           m + null_shift, best_pred + best_shift), 1):
                writer.writerow((i, *values))
    import json
    result_path = plots_dir / f"{output_stem}_fit.json"
    result_path.write_text(json.dumps(dict(hist_config_tag=tag, k=prior.K,
        full_rank=prior.full_K, flux_source=args.flux_source, flux_scale=scale,
        n_universes=prior.n_universes, n_profile_bins=len(profile_bins),
        n_score_bins=len(d), exclude_ratios_from_profiling=args.exclude_ratios_from_profiling,
        chi2_null=null, chi2_best=best, delta_chi2=null-best,
        parameters=dict(zip(("dm2", "ue4_sq", "umu4_sq", "utau4_sq"), map(float, point)))), indent=2)+"\n")
    print(f"Saved fit summary: {result_path}\nSaved spectra: {spectra_out}")
    if not args.no_plots:
        plot_fit_spectra(h, prior, d, m, best_pred, params, cfg,
                         plots_dir,
                         f"{tag}_{mode_tag}", args.ntuple_tag, args.flux_source, scale)


if __name__ == "__main__":
    main()

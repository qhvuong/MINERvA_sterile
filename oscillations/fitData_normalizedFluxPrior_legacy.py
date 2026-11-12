#!/usr/bin/env python3
"""Fit all stitched bins with a full-rank, fixed, normalized flux prior."""
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


def vector(hist):
    return np.array([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX()+1)])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hist-config-tag", required=True)
    parser.add_argument("--flux-source", choices=("residual", "mc"), default="residual",
                        help="residual: match full ROOT covariance; mc: earlier study's MC-only prior")
    parser.add_argument("--flux-scale", type=float, default=1.0,
                        help="Fixed multiplier on flux standard deviations (default: 1)")
    parser.add_argument("--spectra-out", type=Path,
                        help="Optional CSV for null/best spectra and fitted residual corrections")
    args = parser.parse_args()
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
    from scipy.linalg import cho_factor, cho_solve, solve_triangular
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
    prior = NormalizedFluxPrior(C, A, scale=scale)
    prior.assert_flux_covariance(target_flux)
    matched = cho_factor(C + scale**2 * target_flux)
    root_factor = cho_factor(root_full)

    def prediction(dm2, ue, umu, utau):
        h.OscillateHistogram(dm2, ue, umu, utau, False, False)
        return vector(h.GetOscillatedHistogram())

    print(f"Flux source: {args.flux_source}; universes={prior.n_universes}; "
          f"fitted bins={prior.n_bins}; full numerical K={prior.K}; flux scale={scale:g}; "
          f"lambda relative to unscaled P8={1.0/scale**2:g}")
    print("Unscaled normalized covariance closes to the selected ROOT flux covariance.")

    def report(label, pred):
        r = d-pred
        score, residual, penalty, _ = prior.evaluate(r)
        cov_score = float(r @ cho_solve(matched, r))
        root_score = float(r @ cho_solve(root_factor, r))
        np.testing.assert_allclose(score, cov_score, rtol=1e-9, atol=1e-7)
        if args.flux_source == "residual" and scale == 1:
            np.testing.assert_allclose(score, root_score, rtol=1e-9, atol=1e-7)
        print(f"{label}: profile chi2={score:.8f}, residual={residual:.8f}, "
              f"pull penalty={penalty:.8f}, matched covariance chi2={cov_score:.8f}, "
              f"unscaled full ROOT chi2={root_score:.8f}")
        return score

    class NormalizedFitter(OscillationFitter):
        def CalChi2(self, x):
            pred = prediction(100.0*x[0], *x[1:])
            y = solve_triangular(prior.L, d-pred, lower=True, check_finite=False)
            residual, penalty, _ = prior.fit_whitened(y)
            return residual + penalty

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
    if args.spectra_out:
        null_shift = prior.evaluate(d-m)[3]
        best_shift = prior.evaluate(d-best_pred)[3]
        with args.spectra_out.open("w", newline="") as output:
            writer = csv.writer(output)
            writer.writerow(("bin", "data", "null_cv", "best_cv", "null_flux_correction",
                             "best_flux_correction", "null_profiled", "best_profiled"))
            for i, values in enumerate(zip(d, m, best_pred, null_shift, best_shift,
                                           m + null_shift, best_pred + best_shift), 1):
                writer.writerow((i, *values))


if __name__ == "__main__":
    main()

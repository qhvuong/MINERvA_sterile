#!/usr/bin/env python3
"""Compare raw and covariance-normalized fixed-A flux profiling for v4.3.1.

Run with --arrays profile_arrays.npz for the NumPy-only null study, or with
--root for covariance closure, fixed-point Asimov tests and production data fits.
"""
import argparse
import csv
import os
from pathlib import Path
import sys

import numpy as np


TAGS = {
    "ratio_all": "predNueel_p8_v4p3p1",
    "direct_all": "predNueel_noRatio_p8_v4p3p1",
}
BASE_TUNE = "v4p3p1"
TUNES = (BASE_TUNE, "v4p3p0", "v4p4p1")


def histogram_tag(case, tune=BASE_TUNE):
    return TAGS[case].rsplit("_", 1)[0] + "_" + tune


NULL = (0., 0., 0., 0.)
PROBES = {
    "low_dm2": (0.2842873055830976, .394, .489, 0.),
    "mu_disappearance": (2., 0., .1, 0.),
    "high_dm2": (18.8770015348253, .05364601633383929, .026969706695506964, .718),
}
LEGACY_NULL = {"ratio_all": 47.17927, "direct_all": 48.71642}


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def profile(C, A, lam, residual):
    """Solve in bin space. Return total, residual, penalty, and MC shift."""
    alpha = np.linalg.solve(C + A.T @ A / lam, residual)
    a = A @ alpha / lam
    shift = A.T @ a
    left = residual - shift
    penalty = lam * (a @ a)
    reschi = left @ np.linalg.solve(C, left)
    total = residual @ alpha
    np.testing.assert_allclose(reschi + penalty, total, rtol=1e-8, atol=1e-8)
    return float(total), float(reschi), float(penalty), shift


def mean_center(A):
    return A - A.mean(axis=0)


def reference_profile(C, A, lam, residual):
    """Independent whitened SVD solve for the centered-basis fit checks."""
    L = np.linalg.cholesky(C)
    B = np.linalg.solve(L, A.T)
    y = np.linalg.solve(L, residual)
    U, singular, _ = np.linalg.svd(B, full_matrices=False)
    pulls = singular * (U.T @ y) / (singular**2 + lam)
    penalty = lam * (pulls @ pulls)
    left = y - U @ (singular * pulls)
    return float(left @ left + penalty), float(penalty)


def settings(count, ns, draws, seed):
    yield "legacy", 50, 1., 0, np.arange(50)
    for n in sorted(set(ns) | {count}):
        yield "normalized", n, float(n), 0, np.arange(n)
        if n < count:
            rng = np.random.default_rng(np.random.SeedSequence([seed, n]))
            for i in range(1, draws + 1):
                yield "normalized", n, float(n), i, np.sort(rng.choice(count, n, replace=False))


def summarize(case, A, C, d, m, args, templates):
    full = mean_center(A)
    ref, _, _, refshift = profile(C, full, len(A), d - m)
    prior = full.T @ full / len(A)
    L = np.linalg.cholesky(C)
    ref_q = {name: profile(C, full, len(A), v - m)[0] for name, v in templates.items()}
    rows, probes, bands = [], [], []
    for mode, n, lam, draw, ix in settings(len(A), args.nprof, args.draws, args.seed):
        B = A[ix] if mode == "legacy" else mean_center(A[ix])
        total, residual, penalty, shift = profile(C, B, lam, d - m)
        change = np.linalg.solve(L, B.T @ B / lam - prior)
        change = np.linalg.solve(L, change.T).T
        meta = dict(case=case, mode=mode, nprof=n, lam=lam,
                    basis="cv_shifts" if mode == "legacy" else "mean_centered", draw=draw)
        rows.append(dict(**meta, chi2_null=total, residual=residual, penalty=penalty,
                         delta_vs_full=total-ref, max_abs_shift_difference=np.max(abs(shift-refshift)),
                         flux_cov_difference_whitened_opnorm=np.linalg.norm(change, 2)))
        if draw == 0:
            bands += [dict(**meta, point="null", bin=i+1, data=d[i], raw_mc=m[i],
                           profiled_mc=m[i]+shift[i]) for i in range(len(d))]
            if mode == "legacy" and args.arrays:
                np.testing.assert_allclose(total, LEGACY_NULL[case], atol=2e-4, rtol=0)
            print(f"{case:10s} {mode:10s} N={n:4d} lambda={lam:7g} "
                  f"null={total:8.3f} residual={residual:8.3f} penalty={penalty:7.3f}")
        for name, v in templates.items():
            q, *_ = profile(C, B, lam, v - m)
            probes.append(dict(**meta, point=name, asimov_null_rejection_chi2=q,
                               full_reference=ref_q[name], delta_vs_full=q-ref_q[name]))
    return rows, probes, bands


def covariance_check(case, lhs, rhs):
    return dict(case=case, relative_frobenius=np.linalg.norm(lhs-rhs)/np.linalg.norm(lhs),
                root_trace=np.trace(lhs), universe_trace=np.trace(rhs),
                max_abs_difference=np.max(abs(lhs-rhs)),
                passed=bool(np.allclose(lhs, rhs, rtol=1e-6, atol=1e-8)))


def compare_conventions(case, target, root_cov, raw):
    checks = []
    for basis, B in (("cv_shifts", raw), ("mean_centered", mean_center(raw))):
        check = covariance_check(case, root_cov, B.T @ B / len(B))
        checks.append(dict(**check, comparison=target, basis=basis))
    return checks


def vec(hist):
    return np.array([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX()+1)])


def root_inputs(args):
    osc = Path(__file__).resolve().parents[1]
    for key in ("USE_RELATIVE_FLUX_A", "USE_HYBRID_EXCLUDED_FLUX_COV"):
        if os.environ.get(key) == "1":
            raise ValueError(f"Unset {key} for the fixed-A, profile-all comparison.")
    sys.path.insert(0, str(osc.parent))
    os.environ.setdefault("CCNUEROOT", str(osc.parent))
    sys.argv = [sys.argv[0]]
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    from tools.StitchedHistogram import StitchedHistogram
    from tools.Helper import TMatrix_to_Numpy
    from tools.Fitters import CheckElasticPredictionInput, OscillationFitter

    for case in args.cases:
        tag = histogram_tag(case, getattr(args, "tune", BASE_TUNE))
        cfg = osc / f"HIST_CONFIG_{tag}.json"
        path = osc / "rootfiles" / f"NuE_stitched_hists_{tag}.root"
        for f in (cfg, path):
            if not f.is_file():
                raise FileNotFoundError(f)
        h = StitchedHistogram("sample")
        h.SetHistConfig(str(cfg))
        h.Load(str(path))
        CheckElasticPredictionInput(h, hist_config=str(cfg), verbose=True)
        mh, dh = h.GetMCHistogram(), h.GetDataHistogram()
        A, C, d, m = h.GetAMatrix(), h.GetCovarianceMatrix(True), vec(dh), vec(mh)
        mc_root = TMatrix_to_Numpy(mh.GetSysErrorMatrix("Flux"))[1:-1, 1:-1]
        checks = compare_conventions(case, "MC Flux band", mc_root, A)
        if dh.HasVertErrorBand("Flux"):
            band = dh.GetVertErrorBand("Flux")
            if band.GetNHists() == len(A):
                residual_shifts = np.array([vec(band.GetHist(i)) - d for i in range(len(A))]) - A
                checks.extend(compare_conventions(case, "stitched residual Flux",
                                                  h.GetCovarianceMatrix(False) - C,
                                                  residual_shifts))
        templates = {}
        for name, point in PROBES.items():
            h.OscillateHistogram(*point, False, False)
            templates[name] = vec(h.GetOscillatedHistogram())
        yield case, A, C, d, m, templates, checks, (h, cfg, OscillationFitter)


def run_fits(case, A, C, d, m, backend, args, fit_rows, spectra):
    h, cfg, fitter_class = backend
    configs = [("legacy", 50, 1.)] + [
        ("normalized", n, float(n)) for n in sorted(set(args.fit_nprof) | {len(A)})
    ]
    from scipy.linalg import cho_factor, cho_solve

    class BinSpaceFitter(fitter_class):
        """Use the production optimizer and prediction with the exact profiled quadratic."""

        def CalChi2(self, x):
            self.hist.OscillateHistogram(100*x[0], *x[1:], False, False)
            r = d - np.asarray(self.hist.GetOscillatedHistogram(), dtype=float)[1:-1]
            return float(r @ cho_solve(self.factor, r))

    for mode, n, lam in configs:
        print(f"\nProduction fit {case} {mode}: N={n}, lambda={lam}", flush=True)
        fitter = BinSpaceFitter(h, exclude="", lam=lam, marginalize_flux=True,
                                profile_n_universes=n, hist_config=str(cfg))
        B = A[:n] if mode == "legacy" else mean_center(A[:n])
        fitter.factor = cho_factor(C + B.T @ B / lam)
        h.OscillateHistogram(*NULL, False, False)

        def check_score(use_osc, residual, total, penalty):
            if mode == "legacy":
                expected = fitter.statistic.Chi2DataMC(marginalize=True, useOsc=use_osc)
            else:
                # Production's raw basis is different. Check the same centered objective
                # by solving for the pulls independently, without changing production code.
                expected = reference_profile(C, B, lam, residual)
            np.testing.assert_allclose([total, penalty], expected, atol=1e-6, rtol=1e-7)

        null, _, null_pen, _ = profile(C, B, lam, d-m)
        check_score(False, d-m, null, null_pen)
        best, res = fitter.DoFit()
        point = (res["m"], res["ue4"], res["umu4"], res["utau4"])
        h.OscillateHistogram(*point, False, False)
        mc = vec(h.GetOscillatedHistogram())
        total, resid, penalty, shift = profile(C, B, lam, d-mc)
        check_score(True, d-mc, total, penalty)
        np.testing.assert_allclose(best, total, atol=1e-6, rtol=1e-7)
        meta = dict(case=case, mode=mode, nprof=n, lam=lam,
                    basis="cv_shifts" if mode == "legacy" else "mean_centered")
        fit_rows.append(dict(**meta, chi2_null=null, null_residual=null-null_pen,
                             null_penalty=null_pen, chi2_bf=total, bf_residual=resid,
                             bf_penalty=penalty, delta_chi2=null-total,
                             dm2=point[0], ue4_sq=point[1], umu4_sq=point[2], utau4_sq=point[3]))
        spectra += [dict(**meta, draw=0, point="bf", bin=i+1, data=d[i], raw_mc=mc[i],
                         profiled_mc=mc[i]+shift[i]) for i in range(len(d))]
        write_csv(args.out / "data_fits.csv", fit_rows)
        write_csv(args.out / "spectra.csv", spectra)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--arrays", type=Path)
    src.add_argument("--root", action="store_true")
    p.add_argument("--out", type=Path, default=Path("normalized_flux_v4p3p1"))
    p.add_argument("--cases", nargs="+", choices=TAGS, default=list(TAGS))
    p.add_argument("--nprof", type=int, nargs="+", default=[50, 100, 500, 1000])
    p.add_argument("--fit-nprof", type=int, nargs="+", default=[50, 100, 500])
    p.add_argument("--draws", type=int, default=100, help="Random subsets per N")
    p.add_argument("--seed", type=int, default=431)
    p.add_argument("--skip-fits", action="store_true")
    args = p.parse_args()
    if min(args.nprof + args.fit_nprof) < 1 or args.draws < 0 or args.seed < 0:
        p.error("Universe counts must be positive; draws and seed must be nonnegative")
    args.out.mkdir(parents=True, exist_ok=True)
    if args.arrays:
        with np.load(args.arrays) as z:
            inputs = [(case, *(z[TAGS[case]+"__"+key] for key in ("A", "C_nonflux", "data", "mc")),
                       {}, [], None) for case in args.cases]
    else:
        inputs = root_inputs(args)
    rows, probes, spectra, covariance, fits, backends = [], [], [], [], [], []
    for case, A, C, d, m, templates, checks, backend in inputs:
        if len(A) < max(50, *args.nprof, *args.fit_nprof):
            p.error(f"Requested more universes than the {len(A)} available for {case}")
        r, q, s = summarize(case, A, C, d, m, args, templates)
        rows.extend(r)
        probes.extend(q)
        spectra.extend(s)
        covariance.extend(checks)
        for check in checks:
            print(case, check["comparison"], check["basis"], "relative difference:", check["relative_frobenius"],
                  "PASS" if check["passed"] else "DIFF")
        if backend:
            backends.append((case, A, C, d, m, backend))
    subset_summary = []
    for case in args.cases:
        for n in sorted(set(args.nprof)):
            random_rows = [r for r in rows if r["case"] == case and r["nprof"] == n and r["draw"]]
            if random_rows:
                lo, median, hi = np.percentile([r["chi2_null"] for r in random_rows], [5, 50, 95])
                spread = np.percentile([abs(r["delta_vs_full"]) for r in random_rows], 95)
                subset_summary.append(dict(case=case, nprof=n, draws=len(random_rows),
                                           chi2_q05=lo, chi2_median=median, chi2_q95=hi,
                                           abs_chi2_difference_vs_full_q95=spread))
    for name, data in (("null_profiles", rows), ("subset_summary", subset_summary),
                       ("asimov_probes", probes),
                       ("spectra", spectra), ("covariance_checks", covariance)):
        write_csv(args.out / (name + ".csv"), data)
    if any(not r["passed"] and r["comparison"] == "MC Flux band"
           and r["basis"] == "mean_centered" for r in covariance):
        p.exit(1, "Mean-centered MC Flux covariance did not close. Check covariance_checks.csv before interpreting normalized fits.\n")
    if not args.skip_fits:
        for spec in backends:
            run_fits(*spec, args, fits, spectra)
    print("Results:", args.out.resolve())


if __name__ == "__main__":
    main()

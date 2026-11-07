#!/usr/bin/env python3
"""Compare the production null flux solve against covariance-whitened SVD.

Run from oscillations/ in the fit environment. No ROOT/JSON inputs are modified,
no oscillation fit runs, and no production convention is changed.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import sys

import numpy as np


def svd_reference(A, C, residual, lam):
    """Solve for a in min ||C^{-1/2}(r-A.T a)||^2 + lam ||a||^2."""
    sym = (C + C.T) / 2
    d, U = np.linalg.eigh(sym)
    scale = np.max(np.abs(d))
    if scale == 0 or np.min(d) < -1e-12 * scale:
        raise ValueError("Solve covariance has no positive scale or is indefinite; eigenvalue range {}".format(
            (float(d[0]), float(d[-1]))))
    keep = d > 1e-15 * scale
    W = U[:, keep].T / np.sqrt(d[keep])[:, None]
    B, y = W @ A.T, W @ residual
    u, s, vt = np.linalg.svd(B, full_matrices=False)
    a = vt.T @ ((s / (s*s + lam)) * (u.T @ y))
    return a, W, s


def vec(h):
    return np.array([h.GetBinContent(i) for i in range(1, h.GetNbinsX()+1)])


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def close(a, b):
    return np.allclose(a, b, rtol=2e-6, atol=1e-7)


def monotonic(values, increasing=True):
    x = np.array(values)
    tol = 2e-6 * np.maximum(1., np.maximum(np.abs(x[1:]), np.abs(x[:-1])))
    return bool(np.all((1 if increasing else -1)*np.diff(x) >= -tol))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tune", default="v4p3p1", choices=["v4p3p0", "v4p3p1", "v4p4p1"])
    p.add_argument("--nprof", nargs="+", type=int, default=[40, 1000])
    p.add_argument("--lambdas", nargs="+", type=float, default=[.01, .1, 1., 10., 100.])
    p.add_argument("--out", type=Path)
    args = p.parse_args()
    ns, ls = sorted(set(args.nprof)), sorted(set(args.lambdas))
    if ns[0] <= 0 or any(not np.isfinite(l) or l <= 0 for l in ls):
        p.error("Nprof and lambda must be positive and finite")
    if os.environ.get("USE_HYBRID_EXCLUDED_FLUX_COV", "0") == "1":
        p.error("Hybrid scoring is enabled. Unset USE_HYBRID_EXCLUDED_FLUX_COV for the nominal same-metric audit.")
    osc = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(osc.parent))
    os.environ.setdefault("CCNUEROOT", str(osc.parent))
    sys.argv = [sys.argv[0]]
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    from tools.StitchedHistogram import StitchedHistogram
    from tools.Fitters import Statistics, GetFluxSolveBinIndices

    out = args.out or Path("null_flux_audit_" + args.tune)
    out.mkdir(parents=True, exist_ok=True)
    rows, binrows, tests, saved = [], [], [], {}
    def record(case, test, ok):
        tests.append(dict(case=case, test=test, passed=bool(ok)))
        print(("PASS" if ok else "FAIL"), case, test)

    for direct in (True, False):
        tag = "predNueel_{}p8_{}".format("noRatio_" if direct else "", args.tune)
        cfg_path = osc / ("HIST_CONFIG_" + tag + ".json")
        root_path = osc / "rootfiles" / ("NuE_stitched_hists_" + tag + ".root")
        print("\n[INPUT]", root_path, "\n[CONFIG]", cfg_path, flush=True)
        cfg = json.loads(cfg_path.read_text())
        h = StitchedHistogram("sample")
        h.SetHistConfig(str(cfg_path))
        h.Load(str(root_path))
        data_h, mc_h = h.GetDataHistogram(), h.GetMCHistogram()
        d, m = vec(data_h), vec(mc_h)
        A0 = np.asarray(h.GetAMatrix(), dtype=float)
        C = np.asarray(h.GetCovarianceMatrix(sansFlux=True), dtype=float)
        if A0.shape[1] != len(m) or C.shape != (len(m), len(m)) or ns[-1] > len(A0):
            raise ValueError("Unexpected A/covariance dimensions or too few universes")
        if not all(np.all(np.isfinite(x)) for x in (d, m, A0, C)):
            raise ValueError("Nonfinite data/MC/flux basis/covariance")
        if not close(C, C.T):
            raise ValueError("Nonsymmetric non-flux covariance")
        names = np.full(len(m), "", dtype=object)
        for name, lim in cfg.items():
            ix = slice(lim["start"], lim["end"]+1)
            if np.any(names[ix] != ""):
                raise ValueError("Overlapping sample ranges in " + str(cfg_path))
            names[ix] = name
        if np.any(names == ""):
            raise ValueError("Incomplete sample ranges in " + str(cfg_path))
        ratio_mask = np.isin(names, ["fhc_ratio", "rhc_ratio"])
        ratio_id = vec(h.ratio_id) > .5
        if not np.array_equal(ratio_mask, ratio_id) or (direct and ratio_mask.any()) or (
                not direct and not {"fhc_ratio", "rhc_ratio"}.issubset(cfg)):
            raise ValueError("ROOT ratioID/config disagree with expected mode for " + tag)

        # Diagnostic only: quantify data-side flux shifts omitted by the MC-only A.
        dband = data_h.GetVertErrorBand("Flux") if data_h.HasVertErrorBand("Flux") else None
        dA = (np.array([vec(dband.GetHist(u)) for u in range(len(A0))]) - d
              if dband and dband.GetNHists() == len(A0) else None)
        for i in range(len(m)):
            binrows.append(dict(tag=tag, bin=i+1, sample=names[i], data=d[i], mc=m[i],
                                nonflux_sigma=np.sqrt(max(C[i,i], 0)),
                                mc_flux_rms=np.sqrt(np.mean(A0[:,i]**2)),
                                data_flux_rms=np.sqrt(np.mean(dA[:,i]**2)) if dA is not None else ""))
        def mask_spec(ratios):
            return {name:list(range(1,lim["end"]-lim["start"]+2)) for name,lim in cfg.items()
                    if (name in ("fhc_ratio", "rhc_ratio")) == ratios}
        cases = ([("direct_all", "", None, None, np.ones(len(m),dtype=bool))] if direct else [
            ("ratio_all", "", None, None, np.ones(len(m),dtype=bool)),
            ("ratios_omitted", "ratio", mask_spec(True), None, ~ratio_mask),
            ("ratios_only", "nonratio", mask_spec(False), "ratio", ratio_mask)])
        for case, exclude, mask, profile_only, expected in cases:
            case_rows = []
            for n in ns:
                for lam in ls:
                    stat = Statistics(h, exclude=exclude, lam=lam, mask_spec=mask,
                                      profile_only=profile_only, profile_n_universes=n,
                                      hist_config=str(cfg_path))
                    solve = np.array(GetFluxSolveBinIndices(
                        h.keys, exclude, profile_only=profile_only, mask_bins=stat.mask_bins,
                        hist_config=str(cfg_path)), dtype=int)
                    score = np.array([i for i in range(len(m)) if i not in stat.mask_bins])
                    want = np.flatnonzero(expected)
                    if not (np.array_equal(solve, want) and np.array_equal(score, want)):
                        raise ValueError("{}: unexpected solve or score bins: {} / {} (expected {})".format(
                            case, solve.tolist(), score.tolist(), want.tolist()))
                    A = A0[:n]
                    As, Cs, r = A[:,solve], C[np.ix_(solve,solve)], (d-m)[solve]
                    ref, W, s = svd_reference(As, Cs, r, lam)
                    chi, penalty = stat.Chi2DataMC(marginalize=True, useOsc=False)
                    ff = stat.GetFluxFitter(useOsc=False)
                    a = np.asarray(ff.GetFluxSolution())
                    _, used = ff.GetProfileAMatrixAndUniverses()
                    if not np.array_equal(used, A):
                        raise ValueError("Production A differs from first N rows of saved basis")
                    prof = m + a @ A
                    if not np.array_equal(solve, score):
                        raise ValueError("Audit requires identical solve/score masks")
                    resid = float(np.linalg.norm(W @ (d-prof)[score])**2)
                    norm2 = float(a @ a)
                    base = float(np.linalg.norm(W @ r)**2)
                    refchi = float(np.linalg.norm(W @ (r-ref @ As))**2 + lam*(ref @ ref))
                    B, y = W @ As.T, W @ r
                    rhs = B.T @ y
                    normal_error = float(np.linalg.norm(B.T @ (B @ a-y)+lam*a) /
                                         max(1., np.linalg.norm(rhs)))
                    ok = (close(a, ref) and close(penalty, lam*norm2)
                          and close(chi, resid+penalty) and close(chi, refchi)
                          and chi <= base + 2e-6*max(1, base))
                    record(case, "N={} lambda={:g}: coefficients, penalty and objective".format(n,lam), ok)
                    row = dict(case=case,nprof=n,lambda_value=lam,solve_bins=len(solve),score_bins=len(score),
                               raw_residual=base,chi2_null=float(chi),resid_null=resid,
                               penalty_null=float(penalty),norm2_null=norm2,norm_a_null=np.sqrt(norm2),
                               chi2_svd=refchi,max_coeff_diff=float(np.max(np.abs(a-ref))),
                               relative_normal_equation_error=normal_error,
                               weighted_svd_rank=int(np.sum(s > np.finfo(float).eps*max(B.shape)*s[0])),
                               effective_flux_cov_trace=float(np.sum(As**2)/lam),
                               empirical_flux_cov_trace=float(np.sum(As**2)/n),
                               negative_scored_prediction_bins=int(np.sum(prof[score]<0)),passed=bool(ok))
                    rows.append(row)
                    case_rows.append(row)
                    prefix = "{}__N{}__lambda{:g}".format(case,n,lam)
                    saved[prefix+"__a"] = a
                    saved[prefix+"__shift"] = prof-m
            for n in ns:
                group = [r for r in case_rows if r["nprof"]==n]
                record(case,"N={}: residual nondecreasing in lambda".format(n),
                       monotonic([r["resid_null"] for r in group]))
                record(case,"N={}: norm squared nonincreasing in lambda".format(n),
                       monotonic([r["norm2_null"] for r in group],False))
                record(case,"N={}: objective nondecreasing in lambda".format(n),
                       monotonic([r["chi2_null"] for r in group]))
            for lam in ls:
                record(case,"lambda={:g}: objective nonincreasing with nested Nprof".format(lam),
                       monotonic([r["chi2_null"] for r in case_rows if r["lambda_value"]==lam],False))
            saved[case+"__solve_bins_zero_based"] = want
        saved[tag+"__A"] = A0
        saved[tag+"__C_nonflux"] = C
        saved[tag+"__data"] = d
        saved[tag+"__mc"] = m
        if dA is not None:
            saved[tag+"__data_flux_shift"] = dA

    write_csv(out/"profile_points.csv", rows)
    write_csv(out/"input_bin_flux_rms.csv", binrows)
    write_csv(out/"checks.csv", tests)
    np.savez_compressed(out/"profile_arrays.npz", **saved)
    failures = sum(not t["passed"] for t in tests)
    print("\n{} failures. Results: {}".format(failures, out.resolve()))
    return int(failures != 0)


if __name__ == "__main__":
    sys.exit(main())

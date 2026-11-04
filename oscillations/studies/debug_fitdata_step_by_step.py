#!/usr/bin/env python3
"""Validate fitData predictions and fixed-A profiling without changing the fitter.

Run --write-config PATH first; enter exact best fits in that JSON when available.
Run --self-test to test the numerical machinery without ROOT.
The full run requires the same ROOT/PlotUtils/config environment as fitData.py.
"""
import argparse
import ast
import contextlib
import csv
import hashlib
import itertools
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback

import numpy as np

PARAMETERS = ("dm2", "ue4", "umu4", "utau4")
DEFAULT_POINTS = {
    "null": [0.0, 0.0, 0.0, 0.0],
    "muon_candidate_2": [2.0, 0.0, 0.10, 0.0],
    "high_dm2_candidate": [39.2, 0.002, 0.12, 0.0],
}


def write_csv(path, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with open(path, "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [jsonable(v) for v in value]
    return value


def dump_json(path, value):
    with open(path, "w") as stream:
        json.dump(jsonable(value), stream, indent=2, allow_nan=False)
        stream.write("\n")


def vector(hist):
    return np.array([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)])


def flux_universes(hist):
    band = hist.GetVertErrorBand("Flux")
    return np.array([vector(band.GetHist(k)) for k in range(band.GetNHists())])


def point(value):
    if isinstance(value, dict):
        # Accept fitData's output dictionary too.
        value = [value.get("dm2", value.get("m")), value["ue4"], value["umu4"], value["utau4"]]
    p = np.array(value, dtype=float)
    if p.shape != (4,) or not np.all(np.isfinite(p)):
        raise ValueError("Point must contain finite [dm2, ue4, umu4, utau4]")
    if p[0] < 0 or np.any(p[1:] < 0) or np.sum(p[1:]) > 1 + 1e-10:
        raise ValueError("Unphysical point: {}".format(p))
    return p


def template_arrays(hist, stitched=True):
    """Return rows=reco bins, columns=L/E bins and actual L/E edges."""
    axis = hist.GetYaxis() if stitched else hist.GetXaxis()
    nle = hist.GetNbinsY() if stitched else hist.GetNbinsX()
    nreco = hist.GetNbinsX() if stitched else hist.GetNbinsY()
    edges = np.array([axis.GetBinLowEdge(i) for i in range(1, nle + 1)] + [axis.GetBinUpEdge(nle)])
    contents = np.array([
        [hist.GetBinContent(r, b) if stitched else hist.GetBinContent(b, r)
         for b in range(1, nle + 1)] for r in range(1, nreco + 1)
    ])
    return contents, edges


def independent_sin_average(contents, edges, dm2):
    """Independent stable midpoint/sinc integral, no project sin_average calls.

    Like the production code, assumes uniform density inside each L/E bin.
    Does not validate the physical units or the input MC event weights.
    """
    k = 1.27 * dm2
    mid = (edges[1:] + edges[:-1]) / 2
    width = np.diff(edges)
    if np.any(width <= 0) or not np.all(np.isfinite(contents)):
        raise ValueError("Invalid L/E template")
    avg = 0.5 * (1 - np.cos(2 * k * mid) * np.sinc(k * width / np.pi))
    norm = contents.sum(axis=1)
    return np.divide(contents @ avg, norm, out=np.zeros_like(norm), where=norm != 0)


def independent_prediction(components, averages, p, elastic_id, ratio_id):
    _, ue, um, ut = p
    se, sm, ss = averages
    pee = 1 - 4 * ue * (1 - ue) * se
    pmm = 1 - 4 * um * (1 - um) * sm
    pme = 4 * ue * um * ss
    pmt = 4 * um * ut * sm
    pet = 4 * ue * ut * se
    e, m, s = components
    nue = e * pee + s * pme
    numu = m * pmm
    tau = elastic_id * (m * pmt + e * pet)
    ratio = np.divide(numu, nue, out=np.zeros_like(numu), where=nue != 0)
    total = np.where(ratio_id > 0.5, ratio, numu + nue) + tau
    return total, {"Pee": pee, "Pmumu": pmm, "Pmue": pme,
                   "Pmutau": pmt, "Petau": pet, "electron": nue,
                   "muon": numu, "tau": tau, "ratio": ratio}


def whitened_lstsq(A, covariance, residual, lam):
    """Solve independently by stacked whitened least squares, not normal equations."""
    if lam <= 0:
        raise ValueError("This audit requires lambda > 0 for a unique ridge solution")
    c = (covariance + covariance.T) / 2
    eig, basis = np.linalg.eigh(c)
    if eig.min() < -1e-10 * max(1., np.max(np.abs(eig))):
        raise ValueError("Non-flux covariance is indefinite")
    # Match numpy pinv's default singular-value cutoff used in the fitter.
    keep = eig > max(np.max(np.abs(eig)), 1e-300) * 1e-15
    W = (basis[:, keep] / np.sqrt(eig[keep])).T
    B = np.vstack((W @ A.T, np.sqrt(lam) * np.eye(A.shape[0])))
    y = np.r_[W @ residual, np.zeros(A.shape[0])]
    a = np.linalg.lstsq(B, y, rcond=None)[0]
    return a, W


class Audit:
    def __init__(self, out):
        self.out = out
        self.checks = []

    def check(self, name, passed, detail="", severity="error"):
        status = "PASS" if passed else ("WARN" if severity == "warning" else "FAIL")
        self.checks.append(dict(check=name, status=status, detail=str(detail)))
        print("[{}] {} {}".format(status, name, detail), flush=True)
        return passed

    def close(self):
        write_csv(self.out / "checks.csv", self.checks)
        return sum(c["status"] == "FAIL" for c in self.checks)


def compare(a, b, rtol=2e-7, atol=1e-8):
    a, b = np.asarray(a), np.asarray(b)
    return bool(np.allclose(a, b, rtol=rtol, atol=atol)), float(np.max(np.abs(a - b), initial=0))


def example_config():
    direct = "predNueel_noRatio_p8_v4p3p0"
    ratios = "predNueel_p8_v4p3p0"
    return {
        "note": "Candidate points are approximate diagnostic benchmarks, NOT measured best fits. Enter exact best_fit dictionaries or use --refit.",
        "lambda": 1.0, "nprof": 40,
        "points": DEFAULT_POINTS,
        "ue_scan": [0, 0.0001, 0.001, 0.003, 0.01, 0.03, 0.06, 0.1],
        "runs": [
            {"name": "direct_all", "tag": direct, "exclude": "", "best_fit": None},
            {"name": "ratio_all", "tag": ratios, "exclude": "", "best_fit": None},
            {"name": "ratio_no_ratio_profile", "tag": ratios, "exclude": "ratio", "best_fit": None},
        ],
    }


class Run:
    def __init__(self, spec, config, repo, audit, modules):
        self.spec, self.name, self.audit, self.modules = spec, spec["name"], audit, modules
        self.out = audit.out / self.name
        self.out.mkdir(parents=True, exist_ok=True)
        self.lam = float(spec.get("lambda", config.get("lambda", 1)))
        self.nprof = int(spec.get("nprof", config.get("nprof", 40)))
        self.exclude = spec.get("exclude", "") or ""
        if self.exclude.lower() == "none":
            self.exclude = ""
        self.profile_only = spec.get("profile_only")
        self.mask_spec = spec.get("mask_spec")
        tag = spec.get("tag")
        def resolve(value, fallback):
            return str((repo / "oscillations" / (value or fallback)).resolve())
        self.hist_config = resolve(spec.get("hist_config"), "HIST_CONFIG_{}.json".format(tag))
        self.root_path = resolve(spec.get("root_file"), "rootfiles/NuE_stitched_hists_{}.root".format(tag))
        with open(self.hist_config) as f:
            self.cfg = json.load(f)
        self.hist = modules["StitchedHistogram"]("sample")
        self.hist.SetHistConfig(self.hist_config)
        self.hist.Load(self.root_path)
        self.mc, self.data = vector(self.hist.GetMCHistogram()), vector(self.hist.GetDataHistogram())
        self.cov = self.hist.GetCovarianceMatrix(sansFlux=True)
        self.full_cov = self.hist.GetCovarianceMatrix(sansFlux=False)
        self.A = np.asarray(self.hist.GetAMatrix(), dtype=float)[:self.nprof].copy()
        if len(self.A) != self.nprof:
            raise ValueError("Requested more profiling universes than available")
        self.ranges = {key: np.arange(v["start"], v["end"] + 1) for key, v in self.cfg.items()}
        self.stat = self.make_stat()
        self.mask = set(self.stat.mask_bins)
        self.keep = np.array([i for i in range(len(self.mc)) if i not in self.mask], dtype=int)
        self.solve = np.array(modules["GetFluxSolveBinIndices"](
            self.hist.keys, self.exclude, profile_only=self.profile_only,
            mask_bins=list(self.mask), hist_config=self.hist_config), dtype=int)
        self.results, self.bin_rows, self.coef_rows = {}, [], []
        self.components = [vector(self.hist.nue_hist), vector(self.hist.numu_hist), vector(self.hist.swap_hist)]
        self.elastic_id, self.ratio_id = vector(self.hist.elastic_id), vector(self.hist.ratio_id)
        self.templates = [template_arrays(t) for t in
                          (self.hist.nue_template, self.hist.numu_template, self.hist.swap_template)]
        self.best_fit = point(spec["best_fit"]) if spec.get("best_fit") is not None else None
        self.meta = dict(name=self.name, root_file=self.root_path, hist_config=self.hist_config,
                         lambda_value=self.lam, nprof=self.nprof, exclude=self.exclude,
                         solve_bins_zero_based=self.solve, score_bins_zero_based=self.keep,
                         covariance_rank=int(np.linalg.matrix_rank(self.cov)),
                         A_sha256=hashlib.sha256(self.A.tobytes()).hexdigest())
        self.bin_info = {}
        for sample, ix in self.ranges.items():
            local = self.hist.mc_hists[sample]
            selected = [i for i in range(1, local.GetNbinsX()+1)
                        if local.GetBinContent(i) > modules['minBinCont']]
            if len(selected) != len(ix):
                continue
            for i, original_bin in zip(ix, selected):
                self.bin_info[int(i)] = dict(
                    source_bin=original_bin,
                    energy_low_GeV=float(local.GetBinLowEdge(original_bin)),
                    energy_high_GeV=float(local.GetBinLowEdge(original_bin)+local.GetBinWidth(original_bin)))
        dump_json(self.out / "setup.json", self.meta)

    def make_stat(self):
        return self.modules["Statistics"](
            self.hist, exclude=self.exclude, lam=self.lam, mask_spec=self.mask_spec,
            profile_only=self.profile_only, profile_n_universes=self.nprof,
            hist_config=self.hist_config)

    def check(self, name, passed, detail="", severity="error"):
        return self.audit.check(self.name + ": " + name, passed, detail, severity)

    def inputs(self):
        n = len(self.mc)
        covered = np.concatenate(list(self.ranges.values()))
        self.check("config covers every bin exactly once", np.array_equal(np.sort(covered), np.arange(n)))
        self.check("finite inputs", all(np.all(np.isfinite(a)) for a in
                                        [self.data, self.mc, self.cov, self.full_cov, self.A]))
        self.check("positive regularization", self.lam > 0)
        self.check("nonnegative nominal MC", np.all(self.mc >= 0))
        print("  solve bins: {} of {}; exclude applies to SOLVE, not scoring".format(len(self.solve), n))
        rows = []
        for name, ix in self.ranges.items():
            local_mc = self.hist.mc_hists.get(name)
            local_data = self.hist.data_hists.get(name)
            if local_mc is None or local_data is None:
                self.check('configured sample exists in ROOT: ' + name, False)
                continue
            selected = [b for b in range(1, local_mc.GetNbinsX()+1)
                        if local_mc.GetBinContent(b) > self.modules['minBinCont']]
            mapping_ok = len(selected) == len(ix)
            self.check('retained bin count matches ROOT: ' + name, mapping_ok,
                       'config={}, ROOT={}'.format(len(ix), len(selected)))
            if mapping_ok:
                for label, local, stitched in [('MC', local_mc, self.mc),
                                               ('data', local_data, self.data)]:
                    expected = np.array([local.GetBinContent(b) for b in selected])
                    ok, delta = compare(stitched[ix], expected)
                    self.check('JSON maps {} to correct {} bins'.format(name, label),
                               ok, 'max abs diff={:.6g}'.format(delta))
            self.check('ratio ID matches JSON sample: ' + name,
                       np.allclose(self.ratio_id[ix], 1. if 'ratio' in name else 0.))
            rows.append(dict(sample=name, nbins=len(ix), mc_sum=self.mc[ix].sum(),
                             data_sum=self.data[ix].sum(),
                             integral_data_mc=float(self.data[ix].sum() / self.mc[ix].sum())
                             if self.mc[ix].sum() != 0 else None))
        write_csv(self.out / "input_samples.csv", rows)
        write_csv(self.out / 'sample_bin_mapping.csv', [dict(global_bin=i+1, **v)
                  for i,v in sorted(self.bin_info.items())])
        for label, cov in [("nonflux", self.cov), ("full", self.full_cov)]:
            eig = np.linalg.eigvalsh((cov + cov.T) / 2)
            self.check(label + " covariance symmetric", np.allclose(cov, cov.T, rtol=1e-8, atol=1e-8))
            self.check(label + " covariance PSD", eig.min() >= -1e-10 * max(1., abs(eig).max()),
                       "min eigenvalue={:.6g}, cond={:.6g}".format(eig.min(), np.linalg.cond(cov)))
        rows = []
        for flavor, (contents, edges), component in zip(("nue", "numu", "swap"), self.templates, self.components):
            norm = contents.sum(axis=1)
            self.check(flavor + " nonnegative template", np.all(contents >= -1e-12))
            self.check(flavor + " nonempty template for nonzero component",
                       np.all(norm[np.abs(component) > 1e-10] > 0))
            self.check(flavor + " nonnegative L/E edges", np.all(edges >= 0),
                       "L/E range {}..{}; verify units are km/GeV (or m/MeV)".format(edges[0], edges[-1]), "warning")
            # Warn on source flow content separately: production drops flow bins.
            source = getattr(self.hist, flavor + "_template")
            flow = np.array([source.GetBinContent(i, 0) + source.GetBinContent(i, source.GetNbinsY()+1)
                             for i in range(1, n+1)])
            self.check(flavor + " no stitched L/E flow content", np.all(np.abs(flow) < 1e-10),
                       "max |underflow+overflow|={:.6g}".format(np.max(np.abs(flow))), "warning")
            for i in range(n):
                rows.append(dict(flavor=flavor, global_bin=i+1, component=component[i],
                                 template_sum=norm[i], flow=flow[i]))
        write_csv(self.out / "template_coverage.csv", rows)
        # Quantify, do not change, data-side flux shifts and universe normalization.
        du = flux_universes(self.hist.GetDataHistogram())
        mu = flux_universes(self.hist.GetMCHistogram())
        if du.shape != mu.shape:
            raise ValueError("Data and MC Flux bands have different shapes")
        dA, mA = du - self.data, mu - self.mc
        residual_A = mA - dA
        rows = []
        for sample, ix in self.ranges.items():
            md, dd = np.linalg.norm(mA[:self.nprof, ix]), np.linalg.norm(dA[:self.nprof, ix])
            rows.append(dict(sample=sample, mc_flux_norm=md, data_flux_norm=dd,
                             data_to_mc_response_norm=dd/md if md else None,
                             residual_response_norm=np.linalg.norm(residual_A[:self.nprof, ix])))
        write_csv(self.out / "data_side_flux_response.csv", rows)
        np.savez_compressed(self.out / "input_matrices.npz", A_used=self.A,
                            A_residual_firstN=residual_A[:self.nprof], covariance_nonflux=self.cov,
                            covariance_full=self.full_cov, flux_covariance=self.full_cov-self.cov,
                            mc=self.mc, data=self.data, solve_bins=self.solve, score_bins=self.keep)

    def averages(self, dm2):
        return [independent_sin_average(c, e, dm2) for c, e in self.templates]

    def oscillation_checks(self, points):
        tests = dict(points)
        tests.update({"zero_mixing_dm2_2": [2, 0, 0, 0],
                      "zero_mixing_dm2_39": [39.2, 0, 0, 0],
                      "zero_dm2_nonzero_mixing": [0, .1, .15, .2],
                      "muon_only": [2, 0, .1, 0],
                      "electron_only": [2, .03, 0, 0],
                      "muon_plus_tau": [2, 0, .1, .2]})
        rows = []
        for label, p in tests.items():
            p = point(p)
            self.hist.OscillateHistogram(*p)
            actual = vector(self.hist.GetOscillatedHistogram())
            avg = self.averages(p[0])
            expected, parts = independent_prediction(self.components, avg, p, self.elastic_id, self.ratio_id)
            ok, delta = compare(actual, expected)
            self.check("independent oscillation reconstruction " + label, ok, "max abs diff={:.5g}".format(delta))
            self.check("finite nonnegative prediction " + label,
                       np.all(np.isfinite(actual)) and np.all(actual >= -1e-9))
            for flavor, vals in zip(("nue", "numu", "swap"), avg):
                self.check("sin2 bounded {} {}".format(label, flavor), np.all((vals >= -1e-10) & (vals <= 1+1e-10)))
            if p[0] == 0 or np.all(p[1:] == 0):
                ok, delta = compare(actual, self.mc)
                self.check("null closure " + label, ok, "max abs diff={:.5g}".format(delta))
            if label in ("muon_only", "electron_only"):
                target = "nue_selection" if label == "muon_only" else "numu_selection"
                for sample, ix in self.ranges.items():
                    if target in sample:
                        ok, delta = compare(actual[ix], self.mc[ix])
                        self.check(label + " leaves " + sample + " unchanged", ok, str(delta))
            for sample, ix in self.ranges.items():
                for local, i in enumerate(ix):
                    row = dict(test=label, sample=sample, bin=local+1, global_bin=i+1,
                               dm2=p[0], ue4=p[1], umu4=p[2], utau4=p[3], nominal=self.mc[i],
                               actual=actual[i], independent=expected[i],
                               survival_ratio=actual[i]/self.mc[i] if self.mc[i] else None,
                               sin2_nue=avg[0][i], sin2_numu=avg[1][i], sin2_swap=avg[2][i])
                    row.update({key: val[i] for key, val in parts.items()})
                    rows.append(row)
            # Local/sub-histogram route vs stitched route, including removed-bin mapping.
            # This also checks template transposition and flavor-component stitching.
            for sample, ix in self.ranges.items():
                local_mc = self.hist.mc_hists[sample]
                physical_bins = [i for i in range(1, local_mc.GetNbinsX()+1)
                                 if local_mc.GetBinContent(i) > self.modules["minBinCont"]]
                if len(physical_bins) != len(ix):
                    self.check("local-to-stitched mapping " + sample, False,
                               "selected local bins do not match config; skip local comparison")
                    continue
                try:
                    with open(os.devnull, "w") as quiet, contextlib.redirect_stdout(quiet):
                        _, sub = self.hist.OscillateSubHistogram(sample, *p)
                    selected = np.array([sub.GetBinContent(i) for i in physical_bins])
                    ok, delta = compare(actual[ix], selected)
                    self.check("subhist vs stitched {} {}".format(label, sample), ok, str(delta))
                except Exception as error:
                    self.check("subhist vs stitched {} {}".format(label, sample), False, repr(error))
        write_csv(self.out / "oscillation_bins.csv", rows)
        self.hist.OscillateHistogram(0., 0., 0., 0.)

    def evaluate(self, label, p, record=True):
        p = point(p)
        self.hist.OscillateHistogram(*p)
        raw = vector(self.hist.GetOscillatedHistogram())
        total, penalty = self.stat.Chi2DataMC(marginalize=True, useOsc=True)
        fitter = self.stat.GetFluxFitter(useOsc=True)
        a = fitter.GetFluxSolution().copy()
        profiled, reported_penalty = fitter.MarginalizeFlux()
        _, used = fitter.GetProfileAMatrixAndUniverses()
        ok, delta = compare(used, self.A)
        self.check("fixed A retained " + label, ok, str(delta)) if record else None
        r = self.data - profiled
        score_cov = self.cov[np.ix_(self.keep, self.keep)]
        W = np.linalg.pinv(score_cov)
        raw_r = (self.data - raw)[self.keep]
        raw_nonflux_chi2 = float(raw_r @ W @ raw_r)
        raw_full_chi2 = float(raw_r @ np.linalg.pinv(self.full_cov[np.ix_(self.keep,self.keep)]) @ raw_r)
        q = np.zeros_like(r)
        q[self.keep] = r[self.keep] * (W @ r[self.keep])
        calc = q.sum() + penalty
        if record:
            ok, delta = compare(total, calc)
            self.check("chi2 reconstruction " + label, ok, str(delta))
            ok, delta = compare(penalty, reported_penalty)
            self.check("penalty reconstruction " + label, ok, str(delta))
            ok, delta = compare(profiled, raw + a @ self.A)
            self.check("profile shift reconstruction " + label, ok, str(delta))
        solve_A = self.A[:, self.solve]
        solve_r = (self.data - raw)[self.solve]
        expected_a, whitening = whitened_lstsq(solve_A, self.cov[np.ix_(self.solve,self.solve)], solve_r, self.lam)
        if record:
            ok, delta = compare(a, expected_a, rtol=2e-6, atol=1e-7)
            self.check("independent whitened solve " + label, ok, str(delta))
        solve_zero = float(np.linalg.norm(whitening @ solve_r)**2)
        solve_objective = float(np.linalg.norm(whitening @ (solve_r - a @ solve_A))**2 + penalty)
        if record:
            self.check("solve improves its own objective " + label, solve_objective <= solve_zero + 1e-7 * max(1.,solve_zero))
        summary = dict(case=self.name, point=label, dm2=p[0], ue4=p[1], umu4=p[2], utau4=p[3],
                       total_chi2=float(total), residual_chi2=float(total-penalty),
                       penalty=float(penalty), coefficient_norm=float(np.linalg.norm(a)),
                       max_abs_coefficient=float(np.abs(a).max()),
                       raw_nonflux_chi2=raw_nonflux_chi2, raw_full_cov_chi2=raw_full_chi2,
                       solve_objective=solve_objective, solve_objective_at_zero=solve_zero)
        for sample, ix in self.ranges.items():
            summary["chi2_"+sample] = float(q[ix].sum())
        result = dict(summary=summary, params=p, coefficients=a, raw=raw,
                      profiled=profiled, shift=a @ self.A, q=q)
        if record:
            self.results[label] = result
            for sample, ix in self.ranges.items():
                for j, i in enumerate(ix):
                    self.bin_rows.append(dict(point=label, sample=sample, bin=j+1, global_bin=i+1,
                                              **self.bin_info.get(int(i), {}),
                                              data=self.data[i], nominal=self.mc[i], raw=raw[i],
                                              profiled=profiled[i], flux_shift=profiled[i]-raw[i],
                                              nonflux_sigma=np.sqrt(max(0.,self.cov[i,i])),
                                              marginal_sigma=np.sqrt(max(0.,self.full_cov[i,i])),
                                              chi2_additive=q[i], scored=i in set(self.keep), solved=i in set(self.solve)))
            for k, ak in enumerate(a):
                self.coef_rows.append(dict(point=label, universe_zero_based=k, coefficient=ak,
                                           penalty_component=self.lam*ak*ak))
            print("  {} {}: chi2={:.6f} residual={:.6f} penalty={:.6f}".format(
                self.name, label, total, total-penalty, penalty), flush=True)
        return result

    def profiling_tests(self):
        # Real loaded covariance/basis, controlled artificial data, never written to input ROOT.
        original = self.hist.GetDataHistogram().Clone()
        rng = np.random.default_rng(5481)
        truth = rng.normal(0, .03, size=self.nprof)
        tests = {"data_equals_nominal": self.mc,
                 "single_universe_injection": self.mc + .1*self.A[0],
                 "multiple_universe_injection": self.mc + truth @ self.A,
                 "constant_2percent_offset": self.mc * 1.02}
        summaries = []
        try:
            for label, data in tests.items():
                toy = original.Clone()
                for i, value in enumerate(data, start=1):
                    toy.SetBinContent(i, float(value))
                self.hist.SetDataHistogram(toy)
                self.data = vector(self.hist.GetDataHistogram())
                result = self.evaluate("toy_" + label, [0,0,0,0], record=False)
                a = result["coefficients"]
                target = (self.data-self.mc)[self.solve]
                aa, _ = whitened_lstsq(self.A[:,self.solve], self.cov[np.ix_(self.solve,self.solve)], target, self.lam)
                ok, delta = compare(a, aa, rtol=2e-6, atol=1e-7)
                self.check("controlled profiling " + label, ok, "max |a-a_LS|={:.6g}".format(delta))
                if label == "data_equals_nominal":
                    self.check("D=MC yields zero chi2 and coefficients",
                               np.max(np.abs(a)) < 1e-10 and abs(result["summary"]["total_chi2"]) < 1e-10)
                B = self.A[:,self.solve]
                V = np.linalg.pinv(self.cov[np.ix_(self.solve,self.solve)])
                grad = 2*((B @ V @ B.T + self.lam*np.eye(self.nprof)) @ a - B @ V @ target)
                scale = max(1., np.linalg.norm(B @ V @ target))
                self.check("normal-equation stationarity " + label, np.linalg.norm(grad)/scale < 1e-6,
                           "relative gradient={:.6g}".format(np.linalg.norm(grad)/scale))
                self.check("toy objective decreases " + label,
                           result["summary"]["solve_objective"] <= result["summary"]["solve_objective_at_zero"] + 1e-7)
                summaries.append(result["summary"])
        finally:
            self.hist.SetDataHistogram(original)
            self.data = vector(self.hist.GetDataHistogram())
            self.stat = self.make_stat()
        write_csv(self.out / "controlled_profile_tests.csv", summaries)
        print("  Injected coefficients need not be recovered exactly: ridge shrinks them; comparison uses the penalized optimum.")

    def results_checks(self, points, ue_scan):
        # Null fitData path (useOsc=False) must agree with useOsc=True at zero mixing.
        nominal_total, nominal_pen = self.stat.Chi2DataMC(marginalize=True, useOsc=False)
        for label, p in points.items():
            self.evaluate(label, p)
        if self.best_fit is not None:
            self.evaluate("best_fit", self.best_fit)
        else:
            self.check("exact best fit supplied", False,
                       "candidate points are evaluated, but are NOT labeled best fits; enter best_fit in JSON or --refit", "warning")
        null = self.results["null"]
        ok, delta = compare([nominal_total, nominal_pen], [null["summary"]["total_chi2"], null["summary"]["penalty"]])
        self.check("fitData null vs oscillated-null chi2", ok, str(delta))
        # Repeat points in reversed order to expose histogram mutation / stale state.
        for label in reversed(list(self.results)):
            original = self.results[label]
            repeat = self.evaluate("repeat_"+label, original["params"], record=False)
            ok, delta = compare([repeat["summary"]["total_chi2"],repeat["summary"]["penalty"]],
                                [original["summary"]["total_chi2"],original["summary"]["penalty"]])
            self.check("order-independent evaluation " + label, ok, str(delta))
            ok, delta = compare(original["coefficients"], repeat["coefficients"])
            self.check("order-independent flux coefficients " + label, ok, str(delta))
        # Explicit Ue scan with dm2, Umu and Utau fixed: a conditional slice, NOT a profile scan.
        scan_rows = []
        anchors = {"muon_candidate_2": point(DEFAULT_POINTS["muon_candidate_2"])}
        if self.best_fit is not None:
            anchors["best_fit"] = self.best_fit
        for anchor, p0 in anchors.items():
            for ue in sorted(set([0., float(p0[1])] + list(map(float,ue_scan)))):
                p = p0.copy(); p[1] = ue
                if p[1:].sum() > 1:
                    continue
                row = self.evaluate("ue_slice", p, record=False)["summary"]
                row["anchor"] = anchor
                scan_rows.append(row)
        write_csv(self.out / "ue_conditional_slice.csv", scan_rows)
        summaries = []
        for label, result in self.results.items():
            row = dict(result["summary"])
            for key, val in null["summary"].items():
                if key in ("total_chi2","residual_chi2","penalty") or key.startswith("chi2_"):
                    row["null_minus_"+key] = val - result["summary"][key]
            summaries.append(row)
        write_csv(self.out / "chi2_breakdown.csv", summaries)
        write_csv(self.out / "spectra.csv", self.bin_rows)
        write_csv(self.out / "flux_coefficients.csv", self.coef_rows)
        np.savez_compressed(self.out / "point_results.npz", **{
            label+"__"+key: r[key] for label,r in self.results.items()
            for key in ("params","coefficients","raw","profiled","shift","q")})

    def refit(self):
        fitter = self.modules["OscillationFitter"](
            self.hist, exclude=self.exclude, lam=self.lam, marginalize_flux=True,
            mask_spec=self.mask_spec, profile_only=self.profile_only,
            profile_n_universes=self.nprof, hist_config=self.hist_config)
        chi, result = fitter.DoFit()
        self.best_fit = point(result)
        dump_json(self.out / "refitted_best_fit.json", dict(best_fit=result, chi2=chi))

    def muon_grid(self):
        rows = []
        for dm in [0., .5, 1., 1.5, 2., 2.5, 3., 5., 10., 20., 39.2, 60., 100.]:
            for um in [0., .02, .05, .1, .15, .2]:
                rows.append(self.evaluate("muon_grid", [dm,0.,um,0.], record=False)["summary"])
        write_csv(self.out / "muon_only_fixed_grid.csv", rows)
        # These are fixed grid values, not a global 2D fit or a confidence contour.


def coefficient_comparison(left, lr, right, rr, label):
    a, b = lr["coefficients"], rr["coefficients"]
    common = sorted(set(left.ranges) & set(right.ranges))
    compatible = a.shape == b.shape and left.lam == right.lam
    for sample in common:
        ia, ib = left.ranges[sample], right.ranges[sample]
        compatible &= (left.A[:,ia].shape == right.A[:,ib].shape and
                       np.allclose(left.A[:,ia],right.A[:,ib],rtol=1e-7,atol=1e-8))
    denom = np.linalg.norm(a)*np.linalg.norm(b)
    row = dict(comparison=label, left=left.name, right=right.name,
               left_point=lr["summary"]["point"], right_point=rr["summary"]["point"],
               shared_basis_agrees=bool(compatible),
               coefficient_l2_difference=float(np.linalg.norm(a-b)) if a.shape == b.shape else None,
               coefficient_cosine=float(a@b/denom) if a.shape==b.shape and denom else None)
    shifts=[]
    for sample in common:
        ia, ib=left.ranges[sample],right.ranges[sample]
        if len(ia)!=len(ib):
            continue
        x,y=lr["shift"][ia],rr["shift"][ib]
        row["shift_l2_difference_"+sample]=float(np.linalg.norm(x-y))
        scale=np.sqrt(np.maximum(np.diag(left.cov)[ia],1e-300))
        row["shift_difference_rms_nonflux_sigma_"+sample]=float(np.sqrt(np.mean(((x-y)/scale)**2)))
        for j,(xx,yy) in enumerate(zip(x,y)):
            shifts.append(dict(comparison=label,left=left.name,right=right.name,sample=sample,
                               bin=j+1,left_shift=xx,right_shift=yy,difference=xx-yy))
    return row, shifts


def compare_runs(runs, out, audit):
    rows, shifts=[],[]
    for run in runs:
        for label in run.results:
            if label == "null":
                continue
            row,s=coefficient_comparison(run,run.results["null"],run,run.results[label],"within_run_null_to_"+label)
            rows.append(row);shifts.extend(s)
    for left,right in itertools.combinations(runs,2):
        labels=sorted(set(left.results)&set(right.results))
        for label in labels:
            row,s=coefficient_comparison(left,left.results[label],right,right.results[label],
                                          "own_best_fits" if label=="best_fit" else "same_point_"+label)
            rows.append(row);shifts.extend(s)
            audit.check("basis comparison {} vs {} {}".format(left.name,right.name,label),
                        row["shared_basis_agrees"],"compare coefficients cautiously if basis differs", "warning")
        if left.root_path==right.root_path and left.nprof==right.nprof:
            audit.check("same-input configs have identical stored A/data: {} vs {}".format(left.name,right.name),
                        np.array_equal(left.A,right.A) and np.array_equal(left.data,right.data))
    write_csv(out/"flux_solution_comparisons.csv",rows)
    write_csv(out/"common_sample_shift_comparisons.csv",shifts)
    return rows


def make_plots(runs, out):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib unavailable; CSV/NPZ diagnostics are complete")
        return
    labels=sorted(set(label for r in runs for label in r.results))
    for label in labels:
        fig,ax=plt.subplots(figsize=(9,4))
        for run in runs:
            if label in run.results:
                ax.plot(run.results[label]["coefficients"], marker=".",label=run.name)
        ax.set(xlabel="Flux universe index (zero based)",ylabel="Fitted coefficient",title=label)
        ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/("coefficients_"+label+".png"),dpi=160);plt.close(fig)
    for run in runs:
        if not run.results:
            continue
        targets=[x for x in ("best_fit","muon_candidate_2","high_dm2_candidate") if x in run.results]
        for label in targets:
            fig,axes=plt.subplots(len(run.ranges),1,figsize=(9,2.5*len(run.ranges)),squeeze=False)
            for ax,(sample,ix) in zip(axes[:,0],run.ranges.items()):
                x=np.arange(1,len(ix)+1)
                null=run.results["null"];r=run.results[label]
                # Counts and ratios have different units; fractional presentation allows shared conventions.
                base=run.mc[ix]
                safe=np.where(base!=0,base,1.)
                ax.errorbar(x,run.data[ix]/safe,yerr=np.sqrt(np.maximum(0,np.diag(run.cov)[ix]))/abs(safe),
                            fmt="k.",label="Data; diagonal non-flux sigma")
                for y,style,text in [(null["raw"][ix],"--","Raw null"),(null["profiled"][ix],"-","Profiled null"),
                                     (r["raw"][ix],"--","Raw "+label),(r["profiled"][ix],"-","Profiled "+label)]:
                    ax.plot(x,y/safe,style,label=text)
                ax.set(title=sample,ylabel="Relative to nominal MC",xlabel="Retained sample bin")
                if all(int(i) in run.bin_info for i in ix):
                    ticklabels=['{:g}-{:g}'.format(run.bin_info[int(i)]['energy_low_GeV'],
                                                run.bin_info[int(i)]['energy_high_GeV']) for i in ix]
                    ax.set_xticks(x)
                    ax.set_xticklabels(ticklabels,rotation=30,ha='right',fontsize=8)
                    ax.set_xlabel('Energy bins (GeV)')
                ax.legend(fontsize=7,ncol=2)
            fig.suptitle(run.name+": "+label);fig.tight_layout();fig.savefig(run.out/("spectra_"+label+".png"),dpi=160);plt.close(fig)


def write_readme(out):
    (out/"READ_ME_FIRST.txt").write_text("""Read diagnostics in this order:
1. checks.csv: FAIL is a failed invariant; WARN needs interpretation. Exit code is nonzero for FAIL.
2. CASE/input_samples.csv and template_coverage.csv: confirm actual fitted samples and components.
3. CASE/oscillation_bins.csv: actual vs independent oscillation prediction, probabilities and L/E averages.
   Both reconstruction paths use the same input templates; agreement does NOT prove input units/weights are physical.
4. CASE/controlled_profile_tests.csv: actual profiling vs independent whitened least squares.
   Injection recovery is compared to the RIDGE optimum, not the injected coefficients.
5. CASE/chi2_breakdown.csv: exact additive sample contributions q_i = r_i*(V^-1 r)_i.
   Sample contributions can be negative. Null-minus columns sum with the penalty to total improvement.
6. CASE/ue_conditional_slice.csv: fix dm2, Umu, Utau and vary Ue, re-solving flux each time.
   This is a conditional slice, NOT a minimization over the other oscillation parameters.
7. flux_solution_comparisons.csv: within-run null/BF, cross-run identical points, and own best fits.
   Compare induced shifts as well as coefficients. Redundant universe directions can make coefficients misleading.
   shared_basis_agrees checks common sample columns, not external flux-universe identity/provenance.
8. common_sample_shift_comparisons.csv and coefficient/spectra PNGs.
9. CASE/data_side_flux_response.csv: quantify the data-side Flux response; no correction is applied.

A is fixed and lambda/Nprof are unchanged. The script forces relative-A and hybrid-covariance
switches OFF and records previous values. --exclude ratio removes ratio bins from the flux
SOLVE ONLY, following fitData; ratios remain in the final residual chi2.
The three default runs use the same common diagnostic points. Approximate candidates are
never claimed as measured best fits. Supply exact best_fit values in config, or use --refit.
If all inputs and predictions pass, a Ue=0 optimum can still be physical: muon-only disappearance
changes CCnumu, ratios and nu+e while direct CCnue stays unchanged; the combined profiled
objective can prefer this even when the inclusive CCnumu integral is nearly correct.
Plot error bars show only covariance diagonals. CSV chi2 uses the full correlated covariance.
No production ROOT files or fitter source are modified. No HIST_CONFIG.json is overwritten.
""")


def self_test(repo):
    """Execute real project classes on fake histograms, without ROOT imports."""
    from scipy import optimize
    ns=dict(np=np,os=os,optimize=optimize,_PRINTED_FLUXFITTER_DIAGNOSTIC=True,
            _PRINTED_FLUX_SVD=set(),_PRINTED_FLUX_WEIGHTED_SVD=set())
    for rel,names in [("tools/Helper.py",{"slicer"}),
                      ("tools/Fitters.py",{"FluxFitter","Statistics"})]:
        tree=ast.parse((repo/rel).read_text())
        nodes=[n for n in tree.body if isinstance(n,(ast.ClassDef,ast.FunctionDef)) and n.name in names]
        exec(compile(ast.Module(body=nodes,type_ignores=[]),str(repo/rel),"exec"),ns)
    class H:
        def __init__(self,values):self.v=np.r_[0.,values,0.]
        def __array__(self,dtype=None,copy=None):return np.asarray(self.v,dtype=dtype)
        def GetNbinsX(self):return len(self.v)-2
    class Toy:
        keys=["all"]
        def __init__(self,A,C,data,mc):self.A,self.C,self.d,self.m=A,C,H(data),H(mc)
        def GetDataHistogram(self):return self.d
        def GetMCHistogram(self):return self.m
        def GetOscillatedHistogram(self):return self.m
        def GetFluxUniverses(self):return self.A+self.m.v[None,1:-1]
        def GetAMatrix(self):return self.A
        def GetCovarianceMatrix(self,sansFlux=False):return self.C if sansFlux else self.C+self.A.T@self.A
    rng=np.random.default_rng(22)
    A=rng.normal(size=(40,6));B=rng.normal(size=(6,6));C=B@B.T+np.eye(6)
    mc=np.arange(6)+100.;data=mc+rng.normal(size=6)
    ns['GetFluxSolveBinIndices']=lambda *a,**kw:list(range(6))
    for subset in [list(range(6)),[0,2,5]]:
        ns['GetFluxSolveBinIndices']=lambda *a,**kw:subset
        for lam in [0.1,1.,40.]:
            toy=Toy(A,C,data,mc)
            stat=ns['Statistics'](toy,exclude='',lam=lam,profile_n_universes=40)
            chi,pen=stat.Chi2DataMC()
            got=stat.GetFluxFitter().GetFluxSolution()
            expected,_=whitened_lstsq(A[:,subset],C[np.ix_(subset,subset)],(data-mc)[subset],lam)
            np.testing.assert_allclose(got,expected,rtol=1e-8,atol=1e-10)
            r=data-mc-got@A
            np.testing.assert_allclose(chi,r@np.linalg.solve(C,r)+lam*got@got,rtol=1e-9)
            if len(subset)==6:
                np.testing.assert_allclose(chi,(data-mc)@np.linalg.solve(C+A.T@A/lam,data-mc),rtol=1e-9)
    # Single nuisance closed form and zero residual.
    ns['GetFluxSolveBinIndices']=lambda *a,**kw:[0]
    toy=Toy(np.array([[2.]]),np.array([[4.]]),[13.],[10.])
    f=ns['FluxFitter'](toy,lam=1,profile_n_universes=1)
    np.testing.assert_allclose(f.GetFluxSolution(),[.75])
    toy.d=H([10.]);f.SetHistogram(toy)
    np.testing.assert_allclose(f.GetFluxSolution(),[0.])
    # Independent sin integrals vs numerical quadrature, including small/large dm2.
    from scipy.integrate import quad
    edges=np.array([0.,.03,.12,.4,.7]);contents=np.array([[2.,3.,1.,4.]])
    for dm in [0.,1e-6,2.,39.2,100.]:
        expected=sum(n*quad(lambda x:np.sin(1.27*dm*x)**2,lo,hi,epsabs=1e-12)[0]/(hi-lo)
                     for n,lo,hi in zip(contents[0],edges[:-1],edges[1:]))/contents.sum()
        np.testing.assert_allclose(independent_sin_average(contents,edges,dm)[0],expected,atol=1e-12)
    # Sanity cases for probabilities, ratio direction, and tau contribution.
    comp=[np.array([100.,100.]),np.array([1000.,1000.]),np.array([20.,20.])]
    av=[np.array([.5,.5])]*3
    total,parts=independent_prediction(comp,av,point([2,0,.1,0]),np.array([0.,1.]),np.array([1.,0.]))
    np.testing.assert_allclose(total,[8.2,920.])
    assert np.allclose(parts['Pee'],1) and np.allclose(parts['Pmumu'],.82)
    print("SELF-TEST PASS: real profile classes, all/subset solves, three lambdas, dense covariance identity,")
    print("single nuisance, zero residual, independent quadrature and pure-muon oscillation/ratio limits.")


def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config',type=Path)
    parser.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[2])
    parser.add_argument('--output',type=Path,default=Path('fitdata_debug'))
    parser.add_argument('--write-config',type=Path)
    parser.add_argument('--self-test',action='store_true')
    parser.add_argument('--stages',default='inputs,oscillation,profile,results,compare',
                        help='Comma separated: inputs,oscillation,profile,results,compare,grid')
    parser.add_argument('--refit',action='store_true',help='Run full production DoFit for every case; potentially slow')
    parser.add_argument('--no-plots',action='store_true')
    args,project_args=parser.parse_known_args()
    repo=args.repo.resolve()
    if args.write_config:
        dump_json(args.write_config,example_config());print('Wrote',args.write_config);return 0
    if args.self_test:
        self_test(repo);return 0
    if not args.config:parser.error('--config is required (or use --write-config / --self-test)')
    stages=set(args.stages.split(','))
    if stages-{'inputs','oscillation','profile','results','compare','grid'}:
        parser.error('Unknown stage')
    if 'compare' in stages:stages.add('results')
    config=json.loads(args.config.read_text())
    specs=config['runs']
    if len({s['name'] for s in specs})!=len(specs):raise ValueError('Run names must be unique')
    for spec in specs:
        if not spec['name'].replace('_','').replace('-','').isalnum():raise ValueError('Use simple alphanumeric run names')
    points={k:point(v) for k,v in config.get('points',DEFAULT_POINTS).items()}
    points['null']=point([0,0,0,0])
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    with open(out/'run.log','w',buffering=1) as logfile:
        class Tee:
            def __init__(self,*streams):self.streams=streams
            def write(self,s):
                for stream in self.streams:stream.write(s)
            def flush(self):
                for stream in self.streams:stream.flush()
        with contextlib.redirect_stdout(Tee(sys.stdout,logfile)),contextlib.redirect_stderr(Tee(sys.stderr,logfile)):
            audit=Audit(out)
            try:
                switches={name:os.environ.get(name) for name in
                          ('USE_RELATIVE_FLUX_A','USE_HYBRID_EXCLUDED_FLUX_COV')}
                for name in switches:os.environ[name]='0'
                os.environ.setdefault('CCNUEROOT',str(repo))
                # Project configs may parse their own remaining CLI arguments.
                sys.argv=[sys.argv[0]]+project_args
                sys.path.insert(0,str(repo))
                import ROOT
                ROOT.gROOT.SetBatch(True)
                ROOT.TH1.AddDirectory(False)
                from tools.StitchedHistogram import StitchedHistogram,minBinCont
                from tools.Fitters import Statistics,OscillationFitter,GetFluxSolveBinIndices
                modules=dict(StitchedHistogram=StitchedHistogram,Statistics=Statistics,
                             OscillationFitter=OscillationFitter,GetFluxSolveBinIndices=GetFluxSolveBinIndices,
                             minBinCont=minBinCont)
                try:commit=subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip()
                except (subprocess.SubprocessError,OSError):commit='unknown'
                dump_json(out/'invocation.json',dict(config=config,repo=str(repo),commit=commit,
                          source_sha256={p:hashlib.sha256((repo/p).read_bytes()).hexdigest() for p in
                                         ('tools/Fitters.py','tools/StitchedHistogram.py','tools/Helper.py')},
                          previous_environment=switches,fixed_A=True,stages=sorted(stages),project_args=project_args))
                print('Fixed-A diagnostic; production source and ROOT inputs remain unchanged. Commit:',commit)
                runs=[]
                for spec in specs:
                    print('\n=== CASE {} ==='.format(spec['name']),flush=True)
                    run=Run(spec,config,repo,audit,modules);runs.append(run)
                    if 'inputs' in stages:run.inputs()
                    if args.refit:run.refit()
                # Evaluate every run's exact BF in every other configuration too.
                # This separates parameter changes from profiling-selection changes.
                for run in runs:
                    if run.best_fit is not None:
                        points['bf_from_'+run.name]=run.best_fit.copy()
                for run in runs:
                    print('\n=== VALIDATE {} ==='.format(run.name),flush=True)
                    osc_points=dict(points)
                    if run.best_fit is not None:osc_points['best_fit']=run.best_fit
                    if 'oscillation' in stages:run.oscillation_checks(osc_points)
                    if 'profile' in stages:run.profiling_tests()
                    if 'results' in stages:run.results_checks(points,config.get('ue_scan',example_config()['ue_scan']))
                    if 'grid' in stages:run.muon_grid()
                if 'compare' in stages:compare_runs(runs,out,audit)
                if not args.no_plots and 'results' in stages:make_plots(runs,out)
                write_readme(out)
            except Exception:
                audit.check('execution completed',False,traceback.format_exc())
                audit.close()
                raise
            failures=audit.close()
            print('\nCompleted: {} failed checks; results in {}'.format(failures,out))
            return 1 if failures else 0


if __name__=='__main__':
    sys.exit(main())

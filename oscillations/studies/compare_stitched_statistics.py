#!/usr/bin/env python3
"""Read-only comparison of old/new stitched inputs and joint-stat covariance.

Run inside the usual ROOT/PlotUtils analysis environment. No ROOT file is
modified and no oscillation minimization is performed.
"""
import argparse
import contextlib
import io
import json
import os
from pathlib import Path
import sys


def summarize_difference(old, new, np, rtol=1e-10, atol=1e-9):
    old, new = np.asarray(old), np.asarray(new)
    if old.shape != new.shape:
        return {"equal": False, "old_shape": list(old.shape), "new_shape": list(new.shape)}
    diff = new - old
    return {
        "equal": bool(np.allclose(old, new, rtol=rtol, atol=atol)),
        "shape": list(old.shape),
        "max_abs_difference": float(np.max(np.abs(diff))) if diff.size else 0.,
        "different_elements": int(np.count_nonzero(~np.isclose(old, new, rtol=rtol, atol=atol))),
        "trace_difference": float(np.trace(diff)) if diff.ndim == 2 and diff.shape[0] == diff.shape[1] else None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-tag", default="predNueel_noRatio_p8_v4p3p1_updatedStat_newEL")
    parser.add_argument("--new-tag", default="LEonly_v4p3p1_jointStat")
    parser.add_argument("--output", type=Path, default=Path("compare_stitched_statistics.json"))
    args = parser.parse_args()
    osc = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(osc.parent))
    os.environ.setdefault("CCNUEROOT", str(osc.parent))
    sys.argv = [sys.argv[0]]  # The legacy config parser must not see these options.
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    import numpy as np
    from tools.StitchedHistogram import StitchedHistogram, minBinCont
    from tools.Helper import TMatrix_to_Numpy

    def vector(hist):
        return np.array([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)])

    def errors2(hist):
        return np.array([hist.GetBinError(i)**2 for i in range(1, hist.GetNbinsX() + 1)])

    def matrix(value):
        return TMatrix_to_Numpy(value)[1:-1, 1:-1]

    def load(tag):
        root_path = osc / "rootfiles" / ("NuE_stitched_hists_" + tag + ".root")
        config_path = osc / ("HIST_CONFIG_" + tag + ".json")
        for path in (root_path, config_path):
            if not path.is_file():
                raise FileNotFoundError(path)
        h = StitchedHistogram("sample", include_shared_ratio_stats=True)
        h.SetHistConfig(str(config_path))
        print("\nLoading:", root_path)
        h.Load(str(root_path))
        with config_path.open() as f:
            cfg = json.load(f)
        # Inventory records what is stored, rather than inferring correction
        # provenance merely from the absence of our new metadata key.
        f = ROOT.TFile.Open(str(root_path), "READ")
        try:
            inventory = [(key.GetName(), key.GetClassName()) for key in f.GetListOfKeys()
                         if any(token in key.GetName().lower() for token in ("stat", "cov"))]
        finally:
            f.Close()
        return h, cfg, inventory

    old, old_cfg, old_inventory = load(args.old_tag)
    new, new_cfg, new_inventory = load(args.new_tag)
    if old_cfg != new_cfg:
        raise RuntimeError("Old/new stitched layouts differ; align sample/bin order before comparing")
    if not hasattr(new, "joint_stat_inputs"):
        raise RuntimeError("Use the updated StitchedHistogram.py with joint-stat support")
    if not new.joint_stat_inputs:
        raise RuntimeError("New file has no joint-stat metadata; cannot perform the ablation")

    report = {"old_tag": args.old_tag, "new_tag": args.new_tag,
              "stored_stat_cov_keys": {"old": old_inventory, "new": new_inventory},
              "comparisons": {}, "null_chi2": {}}

    def compare(label, a, b):
        result = summarize_difference(a, b, np)
        report["comparisons"][label] = result
        print(label + ":", json.dumps(result))
        return result

    old_data, old_mc = vector(old.data_hist), vector(old.mc_hist)
    new_data, new_mc = vector(new.data_hist), vector(new.mc_hist)
    compare("data CV", old_data, new_data)
    compare("MC CV", old_mc, new_mc)
    compare("stored data statistical variances", errors2(old.data_hist), errors2(new.data_hist))
    compare("stored MC statistical variances", errors2(old.mc_hist), errors2(new.mc_hist))
    old_full, new_full = old.GetCovarianceMatrix(False), new.GetCovarianceMatrix(False)
    old_nonflux, new_nonflux = old.GetCovarianceMatrix(True), new.GetCovarianceMatrix(True)
    compare("full covariance", old_full, new_full)
    compare("nonflux covariance", old_nonflux, new_nonflux)
    compare("global flux covariance", old_full - old_nonflux, new_full - new_nonflux)

    residual_histograms = []
    for h in (old, new):
        residual = h.data_hist.Clone()
        residual.Add(h.mc_hist, -1)
        residual_histograms.append(residual)
    ro, rn = residual_histograms
    compare("ROOT independent statistical matrix",
            matrix(ro.GetStatErrorMatrix()), matrix(rn.GetStatErrorMatrix()))
    names = lambda h: set(map(str, h.GetVertErrorBandNames())) | set(map(str, h.GetLatErrorBandNames()))
    old_names, new_names = names(ro), names(rn)
    report["systematic_band_names"] = {"old": sorted(old_names), "new": sorted(new_names)}
    for name in sorted(old_names & new_names):
        compare("systematic covariance: " + name,
                matrix(ro.GetSysErrorMatrix(name)), matrix(rn.GetSysErrorMatrix(name)))
    for side in ("data_hist", "mc_hist"):
        def flux_universes(h):
            band = getattr(h, side).GetVertErrorBand("Flux")
            return np.array([vector(band.GetHist(i)) for i in range(band.GetNHists())])
        compare("Flux universes: " + side, flux_universes(old), flux_universes(new))

    def chi2(label, residual, covariance):
        value = float(residual @ np.linalg.solve(covariance, residual))
        report["null_chi2"][label] = value
        print("{}: {:.12f}".format(label, value))

    old_r, new_r = old_data - old_mc, new_data - new_mc
    print("\n===== NULL CHI2 CROSS CHECKS =====")
    chi2("old residual / old covariance", old_r, old_full)
    chi2("new residual / new covariance", new_r, new_full)
    chi2("new residual / old covariance", new_r, old_full)
    chi2("old residual / new covariance", old_r, new_full)

    # Same NEW ROOT inputs, temporarily use the pre-patch covariance path.
    # No histograms or files are changed, and metadata is restored afterwards.
    saved_joint = new.joint_stat_inputs
    try:
        new.joint_stat_inputs = {}
        with contextlib.redirect_stdout(io.StringIO()):
            new.SetCovarianceMatrices()
        new_without_joint = new.GetCovarianceMatrix(False)
        print("\n===== JOINT-STAT ABLATION (READ ONLY) =====")
        compare("old full vs new with joint step disabled", old_full, new_without_joint)
        compare("observed change vs added joint-stat step",
                new_full - old_full, new_full - new_without_joint)
        report["joint_stat_step_trace"] = float(np.trace(new_full - new_without_joint))
        print("Joint-stat step trace:", report["joint_stat_step_trace"])
        chi2("new residual / new covariance with joint step disabled", new_r, new_without_joint)
    finally:
        new.joint_stat_inputs = saved_joint
        with contextlib.redirect_stdout(io.StringIO()):
            new.SetCovarianceMatrices()

    labels = {}
    for sample, bounds in new_cfg.items():
        hist = new.mc_hists[sample]
        kept = [i for i in range(1, hist.GetNbinsX() + 1) if hist.GetBinContent(i) > minBinCont]
        stitched_indices = list(range(bounds["start"], bounds["end"] + 1))
        if len(kept) != len(stitched_indices):
            raise RuntimeError("Cannot label bins for " + sample)
        labels.update({global_i: (sample, local_i) for global_i, local_i in zip(stitched_indices, kept)})
    delta = new_full - old_full
    print("\n===== LARGEST COVARIANCE CHANGES =====")
    rows = []
    for flat in np.argsort(np.abs(delta).ravel())[::-1][:12]:
        i, j = np.unravel_index(flat, delta.shape)
        if abs(delta[i, j]) < 1e-9:
            continue
        item = {"row": int(i + 1), "column": int(j + 1),
                "row_sample_bin": labels[i], "column_sample_bin": labels[j],
                "old": float(old_full[i, j]), "new": float(new_full[i, j]),
                "difference": float(delta[i, j])}
        rows.append(item)
        print(json.dumps(item))
    report["largest_covariance_changes"] = rows

    # Show whether the older marginal errors already equal the current
    # diagnostics' marginal or residual variance. This is an audit, not an
    # assumption that one format is more complete than the other.
    print("\n===== ELECTRON SOURCE STATISTICS =====")
    electron_rows = []
    for sample, inputs in saved_joint.items():
        for k, cross in enumerate(inputs["cov_data_mc"], start=1):
            vd, vm = inputs["var_data"][k - 1], inputs["var_mc"][k - 1]
            item = {"sample": sample, "bin": k,
                    "old_data_var": float(old.data_hists[sample].GetBinError(k)**2),
                    "old_MC_var": float(old.mc_hists[sample].GetBinError(k)**2),
                    "new_data_var": vd, "new_MC_var": vm,
                    "cov_data_MC": cross,
                    "joint_residual_var": inputs["var_residual"][k - 1],
                    "independent_var": vd + vm}
            electron_rows.append(item)
            print(json.dumps(item))
    report["electron_statistics"] = electron_rows
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print("\nWrote diagnostic report:", args.output)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

import os
import csv
import argparse
import numpy as np

import ROOT
ROOT.TH1.AddDirectory(False)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def get_hist(root_file, hist_name):
    f = ROOT.TFile.Open(root_file, "READ")
    if not f or f.IsZombie():
        raise RuntimeError("Could not open {}".format(root_file))

    h = f.Get(hist_name)
    if not h:
        print("\nAvailable keys in {}:".format(root_file))
        for key in f.GetListOfKeys():
            obj = key.ReadObj()
            print("  {:40s} {}".format(key.GetName(), obj.ClassName()))
        raise RuntimeError("Could not find histogram {} in {}".format(hist_name, root_file))

    h.SetDirectory(0)
    f.Close()
    return h


def hist_to_array(h):
    return np.array([h.GetBinContent(i) for i in range(1, h.GetNbinsX() + 1)], dtype=float)


def get_universe_arrays(h, band_name):
    band = h.GetVertErrorBand(band_name)
    if not band:
        raise RuntimeError("Histogram does not have vertical error band: {}".format(band_name))

    arrs = []
    for i in range(band.GetNHists()):
        arrs.append(hist_to_array(band.GetHist(i)))

    return np.asarray(arrs, dtype=float)


def relative_shift(univ, cv, eps=1e-20):
    return (univ - cv) / np.where(np.abs(cv) > eps, cv, eps)


def safe_corr(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def load_flux_file(path, hist_name, band_name):
    h = get_hist(path, hist_name)
    cv = hist_to_array(h)
    univ = get_universe_arrays(h, band_name)
    rel = np.asarray([relative_shift(univ[i], cv) for i in range(univ.shape[0])])
    return {
        "path": path,
        "hist": h,
        "cv": cv,
        "univ": univ,
        "rel": rel,
    }


def compare_pair(label, A, B, outdir):
    os.makedirs(outdir, exist_ok=True)

    n = min(A["univ"].shape[0], B["univ"].shape[0])

    rows = []
    for i in range(n):
        A_cv = A["cv"]
        B_cv = B["cv"]
        A_u = A["univ"][i]
        B_u = B["univ"][i]
        A_rel = A["rel"][i]
        B_rel = B["rel"][i]

        diff = A_rel - B_rel

        rows.append({
            "universe": i,
            "A_integral_ratio": float(np.sum(A_u) / np.sum(A_cv)),
            "B_integral_ratio": float(np.sum(B_u) / np.sum(B_cv)),
            "A_rel_norm": float(np.linalg.norm(A_rel)),
            "B_rel_norm": float(np.linalg.norm(B_rel)),
            "same_index_corr": safe_corr(A_rel, B_rel),
            "rms_rel_diff": float(np.sqrt(np.mean(diff * diff))),
            "max_abs_rel_diff": float(np.max(np.abs(diff))),
        })

    safe = label.replace(" ", "_").replace("/", "_")

    csv_path = os.path.join(outdir, safe + "_same_index.csv")
    with open(csv_path, "w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print("\n===== {} =====".format(label))
    print("A:", A["path"])
    print("B:", B["path"])
    print("A CV integral:", np.sum(A["cv"]))
    print("B CV integral:", np.sum(B["cv"]))
    print("A/B CV integral ratio:", np.sum(A["cv"]) / np.sum(B["cv"]))
    print("median A universe/CV:", np.median([r["A_integral_ratio"] for r in rows]))
    print("median B universe/CV:", np.median([r["B_integral_ratio"] for r in rows]))
    print("median same-index corr:", np.nanmedian([r["same_index_corr"] for r in rows]))
    print("median RMS rel diff:", np.nanmedian([r["rms_rel_diff"] for r in rows]))
    print("wrote:", csv_path)

    x = np.array([r["universe"] for r in rows], dtype=float)

    plt.figure(figsize=(8, 5.5))
    plt.plot(x, [r["A_integral_ratio"] for r in rows], marker="o", label="A")
    plt.plot(x, [r["B_integral_ratio"] for r in rows], marker="o", label="B")
    plt.xlabel("Universe index")
    plt.ylabel("Universe integral / CV integral")
    plt.title(label + ": integral ratio vs index")
    plt.grid(True, alpha=0.35)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, safe + "_integral_ratio_vs_index.png"), dpi=200)
    plt.close()

    plt.figure(figsize=(8, 5.5))
    plt.plot(x, [r["same_index_corr"] for r in rows], marker="o")
    plt.xlabel("Universe index")
    plt.ylabel("Correlation of relative shifts")
    plt.title(label + ": same-index correlation")
    plt.grid(True, alpha=0.35)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, safe + "_same_index_corr_vs_index.png"), dpi=200)
    plt.close()

    plt.figure(figsize=(8, 5.5))
    plt.plot(x, [r["rms_rel_diff"] for r in rows], marker="o")
    plt.xlabel("Universe index")
    plt.ylabel("RMS relative-shift difference")
    plt.title(label + ": RMS relative-shift difference")
    plt.grid(True, alpha=0.35)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, safe + "_rms_rel_diff_vs_index.png"), dpi=200)
    plt.close()


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--old-produced", required=True)
    parser.add_argument("--new-produced", required=True)
    parser.add_argument("--existing-cvweighted", required=True)
    parser.add_argument("--existing-unweighted", required=True)

    parser.add_argument("--band", default="Flux")

    parser.add_argument(
        "--outdir",
        default="/exp/minerva/data/users/qvuong/flux_studies/plots/compare_old_new_LE1_numu",
    )

    args = parser.parse_args()

    hist_configs = [
        {
            "tag": "cvweighted",
            "produced_hist": "flux_E_cvweighted",
            "existing_file": args.existing_cvweighted,
            "existing_hist": "flux_E_cvweighted",
        },
        {
            "tag": "unweighted",
            "produced_hist": "flux_E_unweighted",
            "existing_file": args.existing_unweighted,
            "existing_hist": "flux_E_unweighted",
        },
    ]

    os.makedirs(args.outdir, exist_ok=True)

    for cfg in hist_configs:
        tag_outdir = os.path.join(args.outdir, cfg["tag"])

        print("\n\n========================================")
        print("Running comparison for:", cfg["tag"])
        print("produced hist:", cfg["produced_hist"])
        print("existing file:", cfg["existing_file"])
        print("existing hist:", cfg["existing_hist"])
        print("========================================")

        oldp = load_flux_file(args.old_produced, cfg["produced_hist"], args.band)
        newp = load_flux_file(args.new_produced, cfg["produced_hist"], args.band)
        exist = load_flux_file(cfg["existing_file"], cfg["existing_hist"], args.band)

        compare_pair("{} old produced vs existing".format(cfg["tag"]), oldp, exist, tag_outdir)
        compare_pair("{} new produced vs existing".format(cfg["tag"]), newp, exist, tag_outdir)
        compare_pair("{} new produced vs old produced".format(cfg["tag"]), newp, oldp, tag_outdir)

    print("\nDone. Output:", args.outdir)


if __name__ == "__main__":
    main()
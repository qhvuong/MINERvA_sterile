#!/usr/bin/env python3

import os
import csv
import argparse
import numpy as np

import ROOT
ROOT.TH1.AddDirectory(False)

try:
    import PlotUtils
except Exception:
    PlotUtils = None

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# DEFAULT_PRODUCED_DIR = "/exp/minerva/data/users/qvuong/flux_studies/producedFluxes_p6_absolute"
DEFAULT_PRODUCED_DIR = "/exp/minerva/data/users/qvuong/flux_studies/producedFluxes_p6"
# DEFAULT_PRODUCED_DIR = "/pnfs/minerva/persistent/users/qvuong/fluxesOLD"
DEFAULT_PRODUCED_P8_DIR = "/exp/minerva/data/users/qvuong/flux_studies/producedFluxes_p8"
# DEFAULT_EXISTING_FLUX_DIR = "/exp/minerva/app/users/qvuong/MAT_AL9/opt/lib/data/flux"
DEFAULT_EXISTING_FLUX_DIR = "/exp/minerva/data/users/qvuong/flux_studies/producedFluxes_p8"
DEFAULT_USED_FLUX_DIR = "/exp/minerva/app/users/qvuong/MAT_AL9/CC-NuE-XSec/custom_plotutils/data/flux"

# DEFAULT_EXISTING_FLUX_DIR = "/exp/minerva/app/users/qvuong/MAT_AL9/opt/lib/data/flux"
# DEFAULT_PRODUCED_DIR = "/exp/minerva/data/users/qvuong/flux_studies"
# DEFAULT_USED_FLUX_DIR = "/exp/minerva/data/users/qvuong/flux_studies/ordered_numu_flux_universes"

PDG_TO_NAME = {
    14: "numu",
    -14: "numubar",
}

def list_keys(root_file):
    f = ROOT.TFile.Open(root_file, "READ")
    if not f or f.IsZombie():
        raise RuntimeError("Could not open {}".format(root_file))

    keys = []
    for key in f.GetListOfKeys():
        obj = key.ReadObj()
        keys.append((key.GetName(), obj.ClassName()))

    f.Close()
    return keys


def get_hist(root_file, hist_name):
    f = ROOT.TFile.Open(root_file, "READ")
    if not f or f.IsZombie():
        raise RuntimeError("Could not open {}".format(root_file))

    h = f.Get(hist_name)
    if not h:
        print("\nAvailable keys in {}:".format(root_file))
        for name, cls in list_keys(root_file):
            print("  {:40s} {}".format(name, cls))
        raise RuntimeError("Could not find histogram {} in {}".format(hist_name, root_file))

    h.SetDirectory(0)
    f.Close()
    return h


def hist_to_array(h):
    return np.array([h.GetBinContent(i) for i in range(1, h.GetNbinsX() + 1)], dtype=float)


def get_band_names(h):
    names = []
    try:
        v = h.GetVertErrorBandNames()
        for i in range(v.size()):
            names.append(str(v[i]))
    except Exception:
        # PyROOT sometimes exposes vector<string> differently.
        try:
            for x in h.GetVertErrorBandNames():
                names.append(str(x))
        except Exception:
            pass
    return names


def get_universe_arrays(h, band_name):
    band = h.GetVertErrorBand(band_name)
    if not band:
        raise RuntimeError("Histogram does not have vertical error band: {}".format(band_name))

    n = band.GetNHists()
    arrs = []

    for i in range(n):
        hu = band.GetHist(i)
        arrs.append(hist_to_array(hu))

    return np.asarray(arrs, dtype=float)


def relative_shift(univ, cv, eps=1e-20):
    return (univ - cv) / np.where(np.abs(cv) > eps, cv, eps)


def norm_shift(univ, cv):
    return np.linalg.norm(univ - cv)


def rel_norm_shift(univ, cv):
    return np.linalg.norm(relative_shift(univ, cv))


def safe_corr(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    if np.std(a) == 0 or np.std(b) == 0:
        return np.nan

    return float(np.corrcoef(a, b)[0, 1])


def compare_universes(prod_arr, exist_arr, cv_prod, cv_exist, max_compare=None):
    """
    Compare universe shifts, using relative shifts to reduce normalization differences.

    If max_compare is set, only compare the first max_compare universes.
    This is useful for comparing new p8 1000-universe flux files against
    older 100-universe files.
    """
    n_prod_total = prod_arr.shape[0]
    n_exist_total = exist_arr.shape[0]

    n_compare = min(n_prod_total, n_exist_total)
    if max_compare is not None:
        n_compare = min(n_compare, max_compare)

    prod_arr = prod_arr[:n_compare]
    exist_arr = exist_arr[:n_compare]

    n_prod = prod_arr.shape[0]
    n_exist = exist_arr.shape[0]

    print("Universe count before truncation:")
    print("  produced total:", n_prod_total)
    print("  existing total:", n_exist_total)
    print("  comparing first:", n_compare)

    prod_rel = np.asarray([relative_shift(prod_arr[i], cv_prod) for i in range(n_prod)])
    exist_rel = np.asarray([relative_shift(exist_arr[j], cv_exist) for j in range(n_exist)])

    # Direct same-index metrics.
    direct_rows = []
    for i in range(min(n_prod, n_exist)):
        dp = prod_rel[i]
        de = exist_rel[i]

        diff = dp - de
        direct_rows.append({
            "universe": i,
            "corr_same_index": safe_corr(dp, de),
            "rms_rel_diff_same_index": float(np.sqrt(np.mean(diff * diff))),
            "max_abs_rel_diff_same_index": float(np.max(np.abs(diff))),
            "prod_rel_norm": float(np.linalg.norm(dp)),
            "exist_rel_norm": float(np.linalg.norm(de)),
            "prod_abs_norm": float(norm_shift(prod_arr[i], cv_prod)),
            "exist_abs_norm": float(norm_shift(exist_arr[i], cv_exist)),
            "prod_integral_ratio": float(np.sum(prod_arr[i]) / np.sum(cv_prod)),
            "exist_integral_ratio": float(np.sum(exist_arr[i]) / np.sum(cv_exist)),
        })

    # Full best-match matrix using RMS relative difference.
    best_rows = []
    rms_matrix = np.zeros((n_prod, n_exist))
    corr_matrix = np.zeros((n_prod, n_exist))

    for i in range(n_prod):
        for j in range(n_exist):
            diff = prod_rel[i] - exist_rel[j]
            rms = np.sqrt(np.mean(diff * diff))
            corr = safe_corr(prod_rel[i], exist_rel[j])

            rms_matrix[i, j] = rms
            corr_matrix[i, j] = corr if np.isfinite(corr) else -999.0

    for i in range(n_prod):
        best_j_rms = int(np.argmin(rms_matrix[i, :]))
        best_j_corr = int(np.argmax(corr_matrix[i, :]))

        best_rows.append({
            "prod_universe": i,
            "best_existing_by_rms": best_j_rms,
            "best_rms_rel_diff": float(rms_matrix[i, best_j_rms]),
            "same_index_rms_rel_diff": float(rms_matrix[i, i]) if i < n_exist else np.nan,
            "best_existing_by_corr": best_j_corr,
            "best_corr": float(corr_matrix[i, best_j_corr]),
            "same_index_corr": float(corr_matrix[i, i]) if i < n_exist else np.nan,
            "is_best_rms_same_index": int(best_j_rms == i),
            "is_best_corr_same_index": int(best_j_corr == i),
        })

    return direct_rows, best_rows, rms_matrix, corr_matrix


def write_csv(rows, path):
    if len(rows) == 0:
        return
    with open(path, "w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def plot_metric(rows, xkey, ykeys, labels, ylabel, title, outpath):
    plt.figure(figsize=(8, 5.5))

    x = np.array([r[xkey] for r in rows], dtype=float)

    for ykey, label in zip(ykeys, labels):
        y = np.array([r[ykey] for r in rows], dtype=float)
        plt.plot(x, y, marker="o", linewidth=1.5, markersize=4, label=label)

    plt.xlabel("Universe index")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, alpha=0.35)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def plot_best_match(best_rows, outpath, mode="rms"):
    plt.figure(figsize=(8, 5.5))

    x = np.array([r["prod_universe"] for r in best_rows], dtype=float)

    if mode == "rms":
        y = np.array([r["best_existing_by_rms"] for r in best_rows], dtype=float)
        title = "Best existing universe by RMS relative-shift difference"
    else:
        y = np.array([r["best_existing_by_corr"] for r in best_rows], dtype=float)
        title = "Best existing universe by relative-shift correlation"

    plt.plot(x, y, marker="o", linestyle="none", markersize=4)
    plt.plot(x, x, linestyle="--", linewidth=1.5, label="same index")
    plt.xlabel("Produced flux universe index")
    plt.ylabel("Best-matching existing flux universe index")
    plt.title(title)
    plt.grid(True, alpha=0.35)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def plot_matrix(mat, outpath, title, vmin=None, vmax=None):
    plt.figure(figsize=(7, 6))
    plt.imshow(mat, origin="lower", aspect="auto", interpolation="nearest", vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.xlabel("Existing flux universe index")
    plt.ylabel("Produced flux universe index")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def summarize_match(label, direct_rows, best_rows):
    n = len(best_rows)
    n_same_rms = sum(r["is_best_rms_same_index"] for r in best_rows)
    n_same_corr = sum(r["is_best_corr_same_index"] for r in best_rows)

    same_corrs = np.array([r["corr_same_index"] for r in direct_rows], dtype=float)
    same_rms = np.array([r["rms_rel_diff_same_index"] for r in direct_rows], dtype=float)

    print("\n===== comparison summary:", label, "=====")
    print("number of produced universes:", n)
    print("same-index is best by RMS:  {} / {}".format(n_same_rms, n))
    print("same-index is best by corr: {} / {}".format(n_same_corr, n))
    print("same-index corr median:", np.nanmedian(same_corrs))
    print("same-index corr min/max:", np.nanmin(same_corrs), np.nanmax(same_corrs))
    print("same-index RMS rel diff median:", np.nanmedian(same_rms))
    print("same-index RMS rel diff min/max:", np.nanmin(same_rms), np.nanmax(same_rms))

    best_rms = np.array([r["best_rms_rel_diff"] for r in best_rows], dtype=float)
    same_rms_from_best = np.array([r["same_index_rms_rel_diff"] for r in best_rows], dtype=float)

    print("best RMS rel diff median:", np.nanmedian(best_rms))
    print("median same-index RMS / best RMS:", np.nanmedian(same_rms_from_best / best_rms))


def process_pair(produced_file, existing_file, hist_name, band_name, label, outdir):
    print("\n\n################################################")
    print("Comparing:", label)
    print("produced:", produced_file)
    print("existing:", existing_file)
    print("hist:", hist_name)
    print("band:", band_name)
    print("################################################")

    h_prod = get_hist(produced_file, hist_name)
    h_exist = get_hist(existing_file, hist_name)

    print("Produced class:", h_prod.ClassName())
    print("Existing class:", h_exist.ClassName())
    print("Produced nbins:", h_prod.GetNbinsX())
    print("Existing nbins:", h_exist.GetNbinsX())
    print("Produced bands:", get_band_names(h_prod))
    print("Existing bands:", get_band_names(h_exist))

    cv_prod = hist_to_array(h_prod)
    cv_exist = hist_to_array(h_exist)

    if len(cv_prod) != len(cv_exist):
        raise RuntimeError(
            "Bin count mismatch: produced {} vs existing {}".format(
                len(cv_prod),
                len(cv_exist),
            )
        )

    prod_arr = get_universe_arrays(h_prod, band_name)
    exist_arr = get_universe_arrays(h_exist, band_name)

    print("Produced universe array:", prod_arr.shape)
    print("Existing universe array:", exist_arr.shape)

    direct_rows, best_rows, rms_matrix, corr_matrix = compare_universes(
        prod_arr,
        exist_arr,
        cv_prod,
        cv_exist,
    )

    os.makedirs(outdir, exist_ok=True)

    safe_label = label.replace(" ", "_").replace("/", "_")

    write_csv(
        direct_rows,
        os.path.join(outdir, "{}_direct_same_index.csv".format(safe_label)),
    )
    write_csv(
        best_rows,
        os.path.join(outdir, "{}_best_match.csv".format(safe_label)),
    )

    summarize_match(label, direct_rows, best_rows)

    # CV comparison.
    rel_cv_diff = relative_shift(cv_prod, cv_exist)
    print("CV integral produced:", np.sum(cv_prod))
    print("CV integral existing:", np.sum(cv_exist))
    print("CV integral ratio produced/existing:", np.sum(cv_prod) / np.sum(cv_exist))
    print("CV max abs rel diff:", np.max(np.abs(rel_cv_diff)))
    print("CV RMS rel diff:", np.sqrt(np.mean(rel_cv_diff * rel_cv_diff)))

    # Plots.
    plot_metric(
        direct_rows,
        "universe",
        ["prod_integral_ratio", "exist_integral_ratio"],
        ["produced", "existing"],
        "Universe integral / CV integral",
        "{}: universe integral ratios".format(label),
        os.path.join(outdir, "{}_integral_ratio_vs_index.png".format(safe_label)),
    )

    plot_metric(
        direct_rows,
        "universe",
        ["prod_rel_norm", "exist_rel_norm"],
        ["produced", "existing"],
        "Norm of relative flux shift",
        "{}: universe relative-shift norms".format(label),
        os.path.join(outdir, "{}_rel_norm_vs_index.png".format(safe_label)),
    )

    plot_metric(
        direct_rows,
        "universe",
        ["corr_same_index"],
        ["same index corr"],
        "Correlation of relative shifts",
        "{}: same-index correlation".format(label),
        os.path.join(outdir, "{}_same_index_corr_vs_index.png".format(safe_label)),
    )

    plot_metric(
        direct_rows,
        "universe",
        ["rms_rel_diff_same_index"],
        ["same index RMS rel diff"],
        "RMS relative-shift difference",
        "{}: same-index RMS difference".format(label),
        os.path.join(outdir, "{}_same_index_rms_diff_vs_index.png".format(safe_label)),
    )

    plot_best_match(
        best_rows,
        os.path.join(outdir, "{}_best_match_by_rms.png".format(safe_label)),
        mode="rms",
    )

    plot_best_match(
        best_rows,
        os.path.join(outdir, "{}_best_match_by_corr.png".format(safe_label)),
        mode="corr",
    )

    plot_matrix(
        rms_matrix,
        os.path.join(outdir, "{}_rms_matrix.png".format(safe_label)),
        "{}: RMS relative-shift difference".format(label),
    )

    plot_matrix(
        corr_matrix,
        os.path.join(outdir, "{}_corr_matrix.png".format(safe_label)),
        "{}: relative-shift correlation".format(label),
        vmin=-1,
        vmax=1,
    )

def process_pair_maybe_different_histnames(
    produced_file,
    existing_file,
    produced_hist_name,
    existing_hist_name,
    band_name,
    label,
    outdir,
    max_compare=None,
):
    print("\n\n################################################")
    print("Comparing:", label)
    print("produced:", produced_file)
    print("existing/used:", existing_file)
    print("produced hist:", produced_hist_name)
    print("existing hist:", existing_hist_name)
    print("band:", band_name)
    print("################################################")

    h_prod = get_hist(produced_file, produced_hist_name)
    h_exist = get_hist(existing_file, existing_hist_name)

    print("Produced class:", h_prod.ClassName())
    print("Existing class:", h_exist.ClassName())
    print("Produced nbins:", h_prod.GetNbinsX())
    print("Existing nbins:", h_exist.GetNbinsX())
    print("Produced bands:", get_band_names(h_prod))
    print("Existing bands:", get_band_names(h_exist))

    cv_prod = hist_to_array(h_prod)
    cv_exist = hist_to_array(h_exist)

    if len(cv_prod) != len(cv_exist):
        raise RuntimeError(
            "Bin count mismatch: produced {} vs existing {}".format(
                len(cv_prod),
                len(cv_exist),
            )
        )

    prod_arr = get_universe_arrays(h_prod, band_name)
    exist_arr = get_universe_arrays(h_exist, band_name)

    print("Produced universe array:", prod_arr.shape)
    print("Existing universe array:", exist_arr.shape)

    direct_rows, best_rows, rms_matrix, corr_matrix = compare_universes(
        prod_arr,
        exist_arr,
        cv_prod,
        cv_exist,
        max_compare=max_compare,
    )

    os.makedirs(outdir, exist_ok=True)

    safe_label = label.replace(" ", "_").replace("/", "_")

    write_csv(
        direct_rows,
        os.path.join(outdir, "{}_direct_same_index.csv".format(safe_label)),
    )
    write_csv(
        best_rows,
        os.path.join(outdir, "{}_best_match.csv".format(safe_label)),
    )

    summarize_match(label, direct_rows, best_rows)

    rel_cv_diff = relative_shift(cv_prod, cv_exist)

    print("CV integral produced:", np.sum(cv_prod))
    print("CV integral existing:", np.sum(cv_exist))
    print("CV integral ratio produced/existing:", np.sum(cv_prod) / np.sum(cv_exist))
    print("CV max abs rel diff:", np.max(np.abs(rel_cv_diff)))
    print("CV RMS rel diff:", np.sqrt(np.mean(rel_cv_diff * rel_cv_diff)))

    plot_metric(
        direct_rows,
        "universe",
        ["prod_integral_ratio", "exist_integral_ratio"],
        ["produced", "existing/used"],
        "Universe integral / CV integral",
        "{}: universe integral ratios".format(label),
        os.path.join(outdir, "{}_integral_ratio_vs_index.png".format(safe_label)),
    )

    plot_metric(
        direct_rows,
        "universe",
        ["prod_rel_norm", "exist_rel_norm"],
        ["produced", "existing/used"],
        "Norm of relative flux shift",
        "{}: universe relative-shift norms".format(label),
        os.path.join(outdir, "{}_rel_norm_vs_index.png".format(safe_label)),
    )

    plot_metric(
        direct_rows,
        "universe",
        ["corr_same_index"],
        ["same index corr"],
        "Correlation of relative shifts",
        "{}: same-index correlation".format(label),
        os.path.join(outdir, "{}_same_index_corr_vs_index.png".format(safe_label)),
    )

    plot_metric(
        direct_rows,
        "universe",
        ["rms_rel_diff_same_index"],
        ["same index RMS rel diff"],
        "RMS relative-shift difference",
        "{}: same-index RMS difference".format(label),
        os.path.join(outdir, "{}_same_index_rms_diff_vs_index.png".format(safe_label)),
    )

    plot_best_match(
        best_rows,
        os.path.join(outdir, "{}_best_match_by_rms.png".format(safe_label)),
        mode="rms",
    )

    plot_best_match(
        best_rows,
        os.path.join(outdir, "{}_best_match_by_corr.png".format(safe_label)),
        mode="corr",
    )

    plot_matrix(
        rms_matrix,
        os.path.join(outdir, "{}_rms_matrix.png".format(safe_label)),
        "{}: RMS relative-shift difference".format(label),
    )

    plot_matrix(
        corr_matrix,
        os.path.join(outdir, "{}_corr_matrix.png".format(safe_label)),
        "{}: relative-shift correlation".format(label),
        vmin=-1,
        vmax=1,
    )

    # Absolute universe integrals.
    prod_abs_integrals = np.sum(prod_arr, axis=1)
    exist_abs_integrals = np.sum(exist_arr, axis=1)

    abs_int_rows = []
    for i in range(min(len(prod_abs_integrals), len(exist_abs_integrals))):
        ratio = prod_abs_integrals[i] / exist_abs_integrals[i] if exist_abs_integrals[i] != 0 else np.nan
        abs_int_rows.append({
            "universe": i,
            "produced_integral": float(prod_abs_integrals[i]),
            "existing_integral": float(exist_abs_integrals[i]),
            "produced_over_existing": float(ratio),
        })

    cv_ratio = np.sum(cv_prod) / np.sum(cv_exist)
    abs_ratios = np.array([r["produced_over_existing"] for r in abs_int_rows], dtype=float)

    print("Absolute universe integral ratio median:", np.nanmedian(abs_ratios))
    print("Absolute universe integral ratio min/max:", np.nanmin(abs_ratios), np.nanmax(abs_ratios))
    print("CV integral ratio produced/existing:", cv_ratio)

    write_csv(
        abs_int_rows,
        os.path.join(outdir, "{}_absolute_universe_integrals.csv".format(safe_label)),
    )

    plot_metric(
        abs_int_rows,
        "universe",
        ["produced_integral", "existing_integral"],
        ["produced", "existing/used"],
        "Absolute universe integral",
        "{}: absolute universe integrals".format(label),
        os.path.join(outdir, "{}_absolute_integral_vs_index.png".format(safe_label)),
    )

    plot_metric(
        abs_int_rows,
        "universe",
        ["produced_over_existing"],
        ["produced/existing"],
        "Produced universe integral / existing universe integral",
        "{}: absolute universe integral ratio".format(label),
        os.path.join(outdir, "{}_absolute_integral_ratio_vs_index.png".format(safe_label)),
    )




def build_paths(args, playlist, pdg):
    species = PDG_TO_NAME[pdg]

    produced = os.path.join(
        args.produced_dir,
        "LE{}_{}.root".format(playlist, species),
    )

    produced_p8 = os.path.join(
        args.produced_p8_dir,
        "LE{}_{}.root".format(playlist, species),
    )

    existing_unweighted = os.path.join(
        args.existing_flux_dir,
        "flux-g4numiv5-pdg{}-minerva{}.root".format(pdg, playlist),
    )

    existing_cvweighted = os.path.join(
        args.existing_flux_dir,
        "flux-gen2thin-pdg{}-minerva{}.root".format(pdg, playlist),
    )

    used_unweighted = os.path.join(
        args.used_flux_dir,
        "flux-g4numiv5-pdg{}-minerva{}.root".format(pdg, playlist),
    )

    used_cvweighted = os.path.join(
        args.used_flux_dir,
        "flux-gen2thin-pdg{}-minerva{}.root".format(pdg, playlist),
    )

    return {
        "species": species,
        "produced": produced,
        "produced_p8": produced_p8,
        "existing_unweighted": existing_unweighted,
        "existing_cvweighted": existing_cvweighted,
        "used_unweighted": used_unweighted,
        "used_cvweighted": used_cvweighted,
    }


def parse_int_list(s):
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def main():
    parser = argparse.ArgumentParser(
        description="Compare produced vs existing/used LE flux universe ordering."
    )

    parser.add_argument(
        "--playlists",
        default="1,5,13",
        help="Comma-separated LE playlist numbers, e.g. 1,5,13.",
    )

    parser.add_argument(
        "--pdgs",
        default="14,-14",
        help="Comma-separated PDG values, e.g. 14,-14.",
    )

    parser.add_argument(
        "--produced-dir",
        default=DEFAULT_PRODUCED_DIR,
        help="Directory containing produced files like LE1_numu.root.",
    )

    parser.add_argument(
        "--existing-flux-dir",
        default=DEFAULT_EXISTING_FLUX_DIR,
        help="Directory containing official opt/lib flux files.",
    )

    parser.add_argument(
        "--used-flux-dir",
        default=DEFAULT_USED_FLUX_DIR,
        help="Directory containing custom_plotutils flux files used by event selection.",
    )

    parser.add_argument(
        "--hist-name",
        default="flux_E_cvweighted",
        help="Default histogram name for produced/cvweighted comparisons.",
    )

    parser.add_argument(
        "--produced-hist-name",
        default=None,
        help="Override histogram name for produced file. Default uses --hist-name.",
    )

    parser.add_argument(
        "--existing-unweighted-hist-name",
        default="flux_E_unweighted",
        help="Histogram name in unweighted files.",
    )

    parser.add_argument(
        "--existing-cvweighted-hist-name",
        default="flux_E_cvweighted",
        help="Histogram name in cvweighted files.",
    )

    parser.add_argument(
        "--band-name",
        default="Flux",
        help="Vertical error band to compare.",
    )

    parser.add_argument(
        "--outdir",
        default="/exp/minerva/data/users/qvuong/flux_studies/plots/check_all_playlists_pdgs",
    )

    parser.add_argument(
        "--produced-p8-dir",
        default=DEFAULT_PRODUCED_P8_DIR,
        help="Directory containing newly produced p8 files like LE1_numu.root.",
    )

    args = parser.parse_args()

    playlists = parse_int_list(args.playlists)
    pdgs = parse_int_list(args.pdgs)

    produced_hist = args.produced_hist_name or args.hist_name

    for playlist in playlists:
        for pdg in pdgs:
            if pdg not in PDG_TO_NAME:
                print("\nWARNING: unsupported PDG:", pdg)
                continue

            paths = build_paths(args, playlist, pdg)

            tag = "LE{}_pdg{}_{}".format(playlist, pdg, paths["species"])
            outdir = os.path.join(args.outdir, tag)

            print("\n\n================================================")
            print("Checking", tag)
            print("produced:", paths["produced"])
            print("produced p8:", paths["produced_p8"])
            print("used cvweighted:", paths["used_cvweighted"])
            print("used unweighted:", paths["used_unweighted"])
            print("existing cvweighted:", paths["existing_cvweighted"])
            print("existing unweighted:", paths["existing_unweighted"])
            print("outdir:", outdir)
            print("================================================")

            # comparisons = [
                # {
                #     "label": "{}_produced_vs_existing_cvweighted".format(tag),
                #     "produced_file": paths["produced"],
                #     "file": paths["existing_cvweighted"],
                #     "produced_hist": produced_hist,
                #     "existing_hist": args.existing_cvweighted_hist_name,
                #     "max_compare": None,
                # },
                # {
                #     "label": "{}_produced_vs_used_customplotutils_cvweighted".format(tag),
                #     "produced_file": paths["produced"],
                #     "file": paths["used_cvweighted"],
                #     "produced_hist": produced_hist,
                #     "existing_hist": args.existing_cvweighted_hist_name,
                #     "max_compare": None,
                # },
                # {
                #     "label": "{}_produced_vs_existing_unweighted".format(tag),
                #     "produced_file": paths["produced"],
                #     "file": paths["existing_unweighted"],
                #     "produced_hist": produced_hist,
                #     "existing_hist": args.existing_unweighted_hist_name,
                #     "max_compare": None,
                # },
                # {
                #     "label": "{}_produced_vs_used_customplotutils_unweighted".format(tag),
                #     "produced_file": paths["produced"],
                #     "file": paths["used_unweighted"],
                #     "produced_hist": produced_hist,
                #     "existing_hist": args.existing_unweighted_hist_name,
                #     "max_compare": None,
                # },
                # {
                #     "label": "{}_produced_p8_first100_vs_existing_cvweighted".format(tag),
                #     "produced_file": paths["produced_p8"],
                #     "file": paths["existing_cvweighted"],
                #     "produced_hist": produced_hist,
                #     "existing_hist": args.existing_cvweighted_hist_name,
                #     "max_compare": 100,
                # },
                # {
                #     "label": "{}_produced_p8_first100_vs_used_customplotutils_cvweighted".format(tag),
                #     "produced_file": paths["produced_p8"],
                #     "file": paths["used_cvweighted"],
                #     "produced_hist": produced_hist,
                #     "existing_hist": args.existing_cvweighted_hist_name,
                #     "max_compare": 100,
                # },
                # {
                #     "label": "{}_produced_p8_first100_vs_existing_unweighted".format(tag),
                #     "produced_file": paths["produced_p8"],
                #     "file": paths["existing_unweighted"],
                #     "produced_hist": produced_hist,
                #     "existing_hist": args.existing_unweighted_hist_name,
                #     "max_compare": 100,
                # },
                # {
                #     "label": "{}_produced_p8_first100_vs_used_customplotutils_unweighted".format(tag),
                #     "produced_file": paths["produced_p8"],
                #     "file": paths["used_unweighted"],
                #     "produced_hist": produced_hist,
                #     "existing_hist": args.existing_unweighted_hist_name,
                #     "max_compare": 100,
                # },
            # ]

            comparisons = [
                {
                    "label": "{}_produced_vs_existing_cvweighted".format(tag),
                    "produced_file": paths["produced"],
                    "file": paths["existing_cvweighted"],
                    "produced_hist": args.existing_cvweighted_hist_name,
                    "existing_hist": args.existing_cvweighted_hist_name,
                    "max_compare": None,
                },
                {
                    "label": "{}_produced_vs_existing_unweighted".format(tag),
                    "produced_file": paths["produced"],
                    "file": paths["existing_unweighted"],
                    "produced_hist": args.existing_unweighted_hist_name,
                    "existing_hist": args.existing_unweighted_hist_name,
                    "max_compare": None,
                },
            ]

            for comp in comparisons:
                try:
                    process_pair_maybe_different_histnames(
                        produced_file=comp["produced_file"],
                        existing_file=comp["file"],
                        produced_hist_name=comp["produced_hist"],
                        existing_hist_name=comp["existing_hist"],
                        band_name=args.band_name,
                        label=comp["label"],
                        outdir=outdir,
                        max_compare=comp["max_compare"],
                    )
                except Exception as e:
                    print("\nWARNING: comparison failed:", comp["label"])
                    print("  produced:", comp["produced_file"])
                    print("  existing/used:", comp["file"])
                    print("  produced hist:", comp["produced_hist"])
                    print("  existing hist:", comp["existing_hist"])
                    print("  max_compare:", comp["max_compare"])
                    print("  error:", e)


    print("\nDone. Wrote outputs under:", args.outdir)


if __name__ == "__main__":
    main()
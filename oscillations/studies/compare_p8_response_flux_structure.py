#!/usr/bin/env python3

import argparse
import os

import ROOT
import numpy as np
import matplotlib.pyplot as plt


ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)


# ============================================================
# STITCHED BINNING
# ============================================================

SAMPLE_BREAKS = [4, 16, 28, 40]

SAMPLE_LABELS_RATIO = [
    r"FHC $\nu+e$",
    r"FHC $\mu/e$",
    r"RHC $\mu/e$",
    r"FHC $\nu_\mu$",
    r"RHC $\bar{\nu}_\mu$",
]

SAMPLE_LABELS_DIRECT = [
    r"FHC $\nu+e$",
    r"FHC $\nu_e$",
    r"RHC $\bar{\nu}_e$",
    r"FHC $\nu_\mu$",
    r"RHC $\bar{\nu}_\mu$",
]


# ============================================================
# ARGUMENTS
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compare stitched Flux covariance structure between "
            "old and new electron-response productions."
        )
    )

    parser.add_argument(
        "--new-ratio-file",
        default="rootfiles/NuE_stitched_hists_prodNueel_p8.root",
    )

    parser.add_argument(
        "--old-ratio-file",
        default="rootfiles/NuE_stitched_hists_prodNueel_p8_oldResponse.root",
    )

    parser.add_argument(
        "--new-direct-file",
        default="rootfiles/NuE_stitched_hists_prodNueel_noRatio_p8.root",
    )

    parser.add_argument(
        "--old-direct-file",
        default="rootfiles/NuE_stitched_hists_prodNueel_noRatio_p8_oldResponse.root",
    )

    parser.add_argument(
        "--hist",
        default="sample_mc",
    )

    parser.add_argument(
        "--error-band",
        default="Flux",
    )

    parser.add_argument(
        "--outdir",
        default="/exp/minerva/data/users/qvuong/surfaces/plots/p8_response_flux_structure",
    )

    return parser.parse_args()


# ============================================================
# ROOT HELPERS
# ============================================================

def open_root(path):
    f = ROOT.TFile.Open(path, "READ")

    if not f or f.IsZombie():
        raise RuntimeError(
            "Could not open {}".format(path)
        )

    return f


def get_hist(root_file, name):
    hist = root_file.Get(name)

    if not hist:
        raise RuntimeError(
            "Could not find histogram '{}' in {}".format(
                name,
                root_file.GetName(),
            )
        )

    return hist


def get_cv(hist):
    return np.array(
        [
            hist.GetBinContent(i)
            for i in range(1, hist.GetNbinsX() + 1)
        ],
        dtype=float,
    )


def get_universes(hist, band_name):
    if not hist.HasVertErrorBand(band_name):
        raise RuntimeError(
            "Histogram '{}' has no vertical error band '{}'. "
            "Available bands: {}".format(
                hist.GetName(),
                band_name,
                list(hist.GetVertErrorBandNames()),
            )
        )

    band = hist.GetVertErrorBand(band_name)

    nuniv = band.GetNHists()
    nbins = hist.GetNbinsX()

    values = np.zeros(
        (nuniv, nbins),
        dtype=float,
    )

    for u in range(nuniv):
        hu = band.GetHist(u)

        for i in range(1, nbins + 1):
            values[u, i - 1] = hu.GetBinContent(i)

    return values


# ============================================================
# ANALYSIS
# ============================================================

def fractional_shifts(universes, cv):
    delta = np.full_like(
        universes,
        np.nan,
        dtype=float,
    )

    good = cv != 0.0

    delta[:, good] = (
        universes[:, good] - cv[good]
    ) / cv[good]

    return delta


def covariance(delta):
    return np.cov(
        delta,
        rowvar=False,
        ddof=1,
    )


def correlation(cov):
    sigma = np.sqrt(np.diag(cov))

    denom = np.outer(
        sigma,
        sigma,
    )

    return np.divide(
        cov,
        denom,
        out=np.zeros_like(cov),
        where=(denom != 0.0),
    )


def eigensystem(cov):
    eigvals, eigvecs = np.linalg.eigh(cov)

    order = np.argsort(eigvals)[::-1]

    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    eigvals[eigvals < 0.0] = 0.0

    total = np.sum(eigvals)

    if total > 0.0:
        eigfrac = eigvals / total
        cumulative = np.cumsum(eigfrac)
    else:
        eigfrac = np.zeros_like(eigvals)
        cumulative = np.zeros_like(eigvals)

    return eigvals, eigvecs, eigfrac, cumulative


def modes_for_fraction(cumulative, target):
    idx = np.where(
        cumulative >= target
    )[0]

    if len(idx) == 0:
        return len(cumulative)

    return int(idx[0] + 1)


def analyze(root_file, hist_name, band_name):
    hist = get_hist(
        root_file,
        hist_name,
    )

    cv = get_cv(hist)

    universes = get_universes(
        hist,
        band_name,
    )

    delta = fractional_shifts(
        universes,
        cv,
    )

    cov = covariance(delta)
    corr = correlation(cov)

    eigvals, eigvecs, eigfrac, cumulative = eigensystem(
        cov
    )

    sigma = np.sqrt(
        np.diag(cov)
    )

    return {
        "cv": cv,
        "universes": universes,
        "delta": delta,
        "cov": cov,
        "corr": corr,
        "sigma": sigma,
        "eigvals": eigvals,
        "eigvecs": eigvecs,
        "eigfrac": eigfrac,
        "cumulative": cumulative,
    }


# ============================================================
# PLOTTING HELPERS
# ============================================================

def draw_sample_boundaries(ax):
    for b in SAMPLE_BREAKS:
        ax.axvline(
            b + 0.5,
            linestyle="--",
            linewidth=1,
            alpha=0.5,
        )


def draw_sample_labels(ax, labels, nbins):
    edges = [0] + SAMPLE_BREAKS + [nbins]

    for i, label in enumerate(labels):
        center = 0.5 * (
            edges[i] + edges[i + 1]
        ) + 0.5

        ax.text(
            center,
            1.01,
            label,
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=9,
        )


def draw_matrix_boundaries(ax):
    for b in SAMPLE_BREAKS:
        pos = b - 0.5

        ax.axvline(
            pos,
            linestyle="--",
            linewidth=1,
            alpha=0.7,
        )

        ax.axhline(
            pos,
            linestyle="--",
            linewidth=1,
            alpha=0.7,
        )


# ============================================================
# NUMERICAL SUMMARY
# ============================================================

def print_summary(label, old, new):
    print()
    print("=" * 100)
    print(label)
    print("=" * 100)

    old_trace = np.trace(
        old["cov"]
    )

    new_trace = np.trace(
        new["cov"]
    )

    print("universes old =", old["universes"].shape[0])
    print("universes new =", new["universes"].shape[0])

    print()
    print("Trace fractional covariance:")
    print("  old = {:.8e}".format(old_trace))
    print("  new = {:.8e}".format(new_trace))

    if old_trace != 0.0:
        print(
            "  new / old = {:.6f}".format(
                new_trace / old_trace
            )
        )

        print(
            "  fractional change = {:+.3f}%".format(
                100.0 * (
                    new_trace / old_trace - 1.0
                )
            )
        )

    print()
    print("Mean fractional sigma:")
    print(
        "  old = {:.6f}".format(
            np.mean(old["sigma"])
        )
    )
    print(
        "  new = {:.6f}".format(
            np.mean(new["sigma"])
        )
    )

    print()
    print("Median fractional sigma:")
    print(
        "  old = {:.6f}".format(
            np.median(old["sigma"])
        )
    )
    print(
        "  new = {:.6f}".format(
            np.median(new["sigma"])
        )
    )

    print()
    print(
        "{:>10s} {:>12s} {:>12s}".format(
            "variance",
            "old modes",
            "new modes",
        )
    )

    for target in [
        0.50,
        0.68,
        0.90,
        0.95,
        0.99,
    ]:
        print(
            "{:9.0f}% {:12d} {:12d}".format(
                100.0 * target,
                modes_for_fraction(
                    old["cumulative"],
                    target,
                ),
                modes_for_fraction(
                    new["cumulative"],
                    target,
                ),
            )
        )


def print_bin_comparison(label, old, new):
    print()
    print("=" * 100)
    print("{} — BIN-BY-BIN FRACTIONAL FLUX UNCERTAINTY".format(label))
    print("=" * 100)

    print(
        "{:>5s} {:>12s} {:>12s} {:>14s} {:>14s}".format(
            "bin",
            "old sigma",
            "new sigma",
            "new-old",
            "new/old",
        )
    )

    for i in range(len(old["sigma"])):
        so = old["sigma"][i]
        sn = new["sigma"][i]

        ratio = np.nan

        if so != 0.0:
            ratio = sn / so

        print(
            "{:5d} {:11.5f}% {:11.5f}% {:+13.5f}% {:14.6f}".format(
                i + 1,
                100.0 * so,
                100.0 * sn,
                100.0 * (sn - so),
                ratio,
            )
        )


# ============================================================
# SIGMA PLOTS
# ============================================================

def plot_sigma_old_new(
    old,
    new,
    labels,
    title,
    outname,
):
    x = np.arange(
        1,
        len(old["sigma"]) + 1,
    )

    fig, ax = plt.subplots(
        figsize=(11, 6)
    )

    ax.plot(
        x,
        100.0 * old["sigma"],
        marker="o",
        linestyle="--",
        label="Old response",
    )

    ax.plot(
        x,
        100.0 * new["sigma"],
        marker="o",
        label="New response",
    )

    draw_sample_boundaries(ax)
    draw_sample_labels(
        ax,
        labels,
        len(old["sigma"]),
    )

    ax.set_xlabel(
        "Stitched bin"
    )

    ax.set_ylabel(
        "Fractional flux uncertainty [%]"
    )

    ax.set_title(title)

    ax.grid(
        True,
        alpha=0.25,
    )

    ax.legend()

    fig.tight_layout()
    fig.savefig(
        outname,
        dpi=200,
    )
    plt.close(fig)


def plot_sigma_change(
    old,
    new,
    labels,
    title,
    outname,
):
    x = np.arange(
        1,
        len(old["sigma"]) + 1,
    )

    diff = 100.0 * (
        new["sigma"] -
        old["sigma"]
    )

    fig, ax = plt.subplots(
        figsize=(11, 6)
    )

    ax.axhline(
        0.0,
        linewidth=1,
    )

    ax.plot(
        x,
        diff,
        marker="o",
    )

    draw_sample_boundaries(ax)
    draw_sample_labels(
        ax,
        labels,
        len(diff),
    )

    ax.set_xlabel(
        "Stitched bin"
    )

    ax.set_ylabel(
        "New - old fractional flux uncertainty [percentage points]"
    )

    ax.set_title(title)

    ax.grid(
        True,
        alpha=0.25,
    )

    fig.tight_layout()
    fig.savefig(
        outname,
        dpi=200,
    )
    plt.close(fig)


# ============================================================
# MATRIX COMPARISONS
# ============================================================

def plot_matrix_difference(
    old_matrix,
    new_matrix,
    title,
    outname,
    cbar_label,
    symmetric=False,
):
    diff = (
        new_matrix -
        old_matrix
    )

    fig, ax = plt.subplots(
        figsize=(9, 8)
    )

    kwargs = {}

    if symmetric:
        vmax = np.nanmax(
            np.abs(diff)
        )

        kwargs["vmin"] = -vmax
        kwargs["vmax"] = vmax

    im = ax.imshow(
        diff,
        origin="lower",
        aspect="auto",
        **kwargs
    )

    cbar = fig.colorbar(
        im,
        ax=ax,
    )

    cbar.set_label(
        cbar_label
    )

    draw_matrix_boundaries(ax)

    ax.set_xlabel(
        "Stitched bin"
    )

    ax.set_ylabel(
        "Stitched bin"
    )

    ax.set_title(title)

    fig.tight_layout()
    fig.savefig(
        outname,
        dpi=200,
    )
    plt.close(fig)


# ============================================================
# EIGENSTRUCTURE
# ============================================================

def plot_eigen_fraction(
    old,
    new,
    title,
    outname,
):
    x = np.arange(
        1,
        len(old["eigfrac"]) + 1,
    )

    fig, ax = plt.subplots(
        figsize=(8, 6)
    )

    ax.plot(
        x,
        old["eigfrac"],
        marker="o",
        linestyle="--",
        label="Old response",
    )

    ax.plot(
        x,
        new["eigfrac"],
        marker="o",
        label="New response",
    )

    ax.set_yscale(
        "log"
    )

    ax.set_xlabel(
        "Flux covariance eigenmode"
    )

    ax.set_ylabel(
        "Fraction of total variance"
    )

    ax.set_title(title)

    ax.grid(
        True,
        alpha=0.25,
    )

    ax.legend()

    fig.tight_layout()
    fig.savefig(
        outname,
        dpi=200,
    )
    plt.close(fig)


def plot_cumulative(
    old,
    new,
    title,
    outname,
):
    x = np.arange(
        1,
        len(old["cumulative"]) + 1,
    )

    fig, ax = plt.subplots(
        figsize=(8, 6)
    )

    ax.plot(
        x,
        old["cumulative"],
        marker="o",
        linestyle="--",
        label="Old response",
    )

    ax.plot(
        x,
        new["cumulative"],
        marker="o",
        label="New response",
    )

    for level in [
        0.90,
        0.95,
        0.99,
    ]:
        ax.axhline(
            level,
            linestyle="--",
            linewidth=1,
        )

    ax.set_xlabel(
        "Number of eigenmodes"
    )

    ax.set_ylabel(
        "Cumulative fraction of flux variance"
    )

    ax.set_ylim(
        0.0,
        1.02,
    )

    ax.set_title(title)

    ax.grid(
        True,
        alpha=0.25,
    )

    ax.legend()

    fig.tight_layout()
    fig.savefig(
        outname,
        dpi=200,
    )
    plt.close(fig)


# ============================================================
# SAVE NUMERICAL RESULTS
# ============================================================

def save_outputs(
    prefix,
    old,
    new,
):
    np.savetxt(
        prefix + "_old_covariance.csv",
        old["cov"],
        delimiter=",",
    )

    np.savetxt(
        prefix + "_new_covariance.csv",
        new["cov"],
        delimiter=",",
    )

    np.savetxt(
        prefix + "_covariance_difference.csv",
        new["cov"] - old["cov"],
        delimiter=",",
    )

    np.savetxt(
        prefix + "_old_correlation.csv",
        old["corr"],
        delimiter=",",
    )

    np.savetxt(
        prefix + "_new_correlation.csv",
        new["corr"],
        delimiter=",",
    )

    np.savetxt(
        prefix + "_correlation_difference.csv",
        new["corr"] - old["corr"],
        delimiter=",",
    )

    table = np.column_stack(
        (
            np.arange(
                1,
                len(old["sigma"]) + 1,
            ),
            old["sigma"],
            new["sigma"],
            new["sigma"] - old["sigma"],
        )
    )

    np.savetxt(
        prefix + "_sigma_comparison.csv",
        table,
        delimiter=",",
        header=(
            "bin,"
            "sigma_old,"
            "sigma_new,"
            "sigma_new_minus_old"
        ),
        comments="",
    )


# ============================================================
# ONE CONFIGURATION
# ============================================================

def compare_configuration(
    name,
    old,
    new,
    labels,
    outdir,
    tag,
):
    print_summary(
        name,
        old,
        new,
    )

    print_bin_comparison(
        name,
        old,
        new,
    )

    prefix = os.path.join(
        outdir,
        tag,
    )

    plot_sigma_old_new(
        old,
        new,
        labels,
        "{}: old vs new response".format(name),
        prefix + "_sigma_old_vs_new.png",
    )

    plot_sigma_change(
        old,
        new,
        labels,
        "{}: response-induced change".format(name),
        prefix + "_sigma_change.png",
    )

    plot_matrix_difference(
        old["cov"],
        new["cov"],
        "{}: new - old Flux covariance".format(name),
        prefix + "_covariance_difference.png",
        "Change in fractional covariance",
        symmetric=True,
    )

    plot_matrix_difference(
        old["corr"],
        new["corr"],
        "{}: new - old Flux correlation".format(name),
        prefix + "_correlation_difference.png",
        "Change in correlation coefficient",
        symmetric=True,
    )

    plot_eigen_fraction(
        old,
        new,
        "{}: Flux eigenvalue spectrum".format(name),
        prefix + "_eigenvalue_fraction.png",
    )

    plot_cumulative(
        old,
        new,
        "{}: cumulative Flux variance".format(name),
        prefix + "_cumulative_variance.png",
    )

    save_outputs(
        prefix,
        old,
        new,
    )


# ============================================================
# MAIN
# ============================================================

def main():
    args = parse_args()

    os.makedirs(
        args.outdir,
        exist_ok=True,
    )

    f_new_ratio = open_root(
        args.new_ratio_file
    )

    f_old_ratio = open_root(
        args.old_ratio_file
    )

    f_new_direct = open_root(
        args.new_direct_file
    )

    f_old_direct = open_root(
        args.old_direct_file
    )

    new_ratio = analyze(
        f_new_ratio,
        args.hist,
        args.error_band,
    )

    old_ratio = analyze(
        f_old_ratio,
        args.hist,
        args.error_band,
    )

    new_direct = analyze(
        f_new_direct,
        args.hist,
        args.error_band,
    )

    old_direct = analyze(
        f_old_direct,
        args.hist,
        args.error_band,
    )

    # Basic consistency checks
    for label, old, new in [
        ("ratio", old_ratio, new_ratio),
        ("direct", old_direct, new_direct),
    ]:
        if len(old["cv"]) != len(new["cv"]):
            raise RuntimeError(
                "{} bin-count mismatch: old={}, new={}".format(
                    label,
                    len(old["cv"]),
                    len(new["cv"]),
                )
            )

        if old["universes"].shape[0] != new["universes"].shape[0]:
            raise RuntimeError(
                "{} universe-count mismatch: old={}, new={}".format(
                    label,
                    old["universes"].shape[0],
                    new["universes"].shape[0],
                )
            )

    compare_configuration(
        "Ratio configuration",
        old_ratio,
        new_ratio,
        SAMPLE_LABELS_RATIO,
        args.outdir,
        "ratio",
    )

    compare_configuration(
        "Direct CCnue configuration",
        old_direct,
        new_direct,
        SAMPLE_LABELS_DIRECT,
        args.outdir,
        "direct",
    )

    f_new_ratio.Close()
    f_old_ratio.Close()
    f_new_direct.Close()
    f_old_direct.Close()

    print()
    print("Done.")
    print(
        "Output directory:",
        args.outdir,
    )


if __name__ == "__main__":
    main()
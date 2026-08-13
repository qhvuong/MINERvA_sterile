#!/usr/bin/env python

import os
import argparse
import ROOT
import PlotUtils
import numpy as np

from matplotlib import pyplot as plt


ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compare flux covariance/correlation structure of "
            "direct CCnue and CC ratio configurations."
        )
    )

    parser.add_argument(
        "--ratio-file",
        required=True,
        help="Stitched ROOT file for CC-ratio configuration",
    )

    parser.add_argument(
        "--direct-file",
        required=True,
        help="Stitched ROOT file for direct-CCnue configuration",
    )

    parser.add_argument(
        "--error-band",
        default="Flux",
    )

    parser.add_argument(
        "--out-prefix",
        default="flux_covariance_comparison",
    )

    return parser.parse_args()


def get_hist(root_file, name):
    hist = root_file.Get(name)

    if not hist:
        raise RuntimeError(
            "Could not find histogram '{}'".format(name)
        )

    return hist


def get_cv_values(hist):
    return np.array([
        hist.GetBinContent(i)
        for i in range(1, hist.GetNbinsX() + 1)
    ], dtype=float)


def get_bin_centers(hist):
    return np.array([
        hist.GetXaxis().GetBinCenter(i)
        for i in range(1, hist.GetNbinsX() + 1)
    ], dtype=float)


def get_universe_values(hist, band_name):
    if not hist.HasVertErrorBand(band_name):
        raise RuntimeError(
            "Histogram '{}' does not contain vertical error band '{}'. "
            "Available bands: {}".format(
                hist.GetName(),
                band_name,
                list(hist.GetVertErrorBandNames()),
            )
        )

    band = hist.GetVertErrorBand(band_name)
    nuniv = band.GetNHists()

    values = np.zeros(
        (nuniv, hist.GetNbinsX()),
        dtype=float,
    )

    for u in range(nuniv):
        h_u = band.GetHist(u)

        for i in range(1, hist.GetNbinsX() + 1):
            values[u, i - 1] = h_u.GetBinContent(i)

    return values


def fractional_shifts(universe_values, cv):
    """
    Return fractional shifts:
        delta[u,i] = (N[u,i] - CV[i]) / CV[i]
    """

    delta = np.full_like(
        universe_values,
        np.nan,
        dtype=float,
    )

    good = cv != 0.0

    delta[:, good] = (
        universe_values[:, good] - cv[good]
    ) / cv[good]

    return delta


def covariance_matrix(delta):
    """
    Covariance of fractional shifts across universes.

    Input:
        delta shape = (Nuniv, Nbin)

    Output:
        covariance shape = (Nbin, Nbin)
    """

    return np.cov(
        delta,
        rowvar=False,
        ddof=1,
    )


def correlation_matrix(cov):
    sigma = np.sqrt(np.diag(cov))

    denom = np.outer(
        sigma,
        sigma,
    )

    corr = np.divide(
        cov,
        denom,
        out=np.zeros_like(cov),
        where=(denom != 0.0),
    )

    return corr


def eigensystem(cov):
    """
    Eigenvalues sorted largest -> smallest.
    """

    eigvals, eigvecs = np.linalg.eigh(cov)

    order = np.argsort(eigvals)[::-1]

    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    eigvals[eigvals < 0.0] = 0.0

    total = np.sum(eigvals)

    if total > 0.0:
        frac = eigvals / total
        cumulative = np.cumsum(frac)
    else:
        frac = np.zeros_like(eigvals)
        cumulative = np.zeros_like(eigvals)

    return eigvals, eigvecs, frac, cumulative


def modes_for_fraction(cumulative, target):
    idx = np.where(cumulative >= target)[0]

    if idx.size == 0:
        return len(cumulative)

    return int(idx[0] + 1)


def analyze_hist(root_file, hist_name, band_name):
    hist = get_hist(
        root_file,
        hist_name,
    )

    cv = get_cv_values(hist)

    centers = get_bin_centers(hist)

    universe_values = get_universe_values(
        hist,
        band_name,
    )

    delta = fractional_shifts(
        universe_values,
        cv,
    )

    cov = covariance_matrix(delta)

    corr = correlation_matrix(cov)

    eigvals, eigvecs, frac, cumulative = eigensystem(
        cov
    )

    sigma = np.sqrt(np.diag(cov))

    return {
        "hist": hist,
        "cv": cv,
        "centers": centers,
        "universes": universe_values,
        "delta": delta,
        "cov": cov,
        "corr": corr,
        "sigma": sigma,
        "eigvals": eigvals,
        "eigvecs": eigvecs,
        "eigfrac": frac,
        "cumulative": cumulative,
    }


def print_summary(label, result):
    print("\n==================================================")
    print(label)
    print("==================================================")

    print("bins =", len(result["sigma"]))

    print("\nFractional flux uncertainty per bin:")

    for i, sigma in enumerate(result["sigma"]):
        print(
            "  bin {:2d}: {:.4f} ({:.2f}%)".format(
                i + 1,
                sigma,
                100.0 * sigma,
            )
        )

    print("\nEigenvalue summary:")

    print(
        "{:>5s} {:>14s} {:>14s} {:>14s}".format(
            "mode",
            "eigenvalue",
            "fraction",
            "cumulative",
        )
    )

    for i in range(len(result["eigvals"])):
        print(
            "{:5d} {:14.6e} {:14.6f} {:14.6f}".format(
                i + 1,
                result["eigvals"][i],
                result["eigfrac"][i],
                result["cumulative"][i],
            )
        )

    for target in [0.50, 0.68, 0.90, 0.95, 0.99]:
        nmodes = modes_for_fraction(
            result["cumulative"],
            target,
        )

        print(
            "Modes needed for {:5.1f}% variance = {}".format(
                100.0 * target,
                nmodes,
            )
        )


def plot_matrix(
    matrix,
    title,
    outname,
    cbar_label,
    vmin=None,
    vmax=None,
):
    fig, ax = plt.subplots(
        figsize=(7, 6)
    )

    im = ax.imshow(
        matrix,
        origin="lower",
        aspect="auto",
        vmin=vmin,
        vmax=vmax,
    )

    cbar = fig.colorbar(
        im,
        ax=ax,
    )

    cbar.set_label(
        cbar_label
    )

    nbins = matrix.shape[0]

    ax.set_xlabel("Bin")
    ax.set_ylabel("Bin")

    ax.set_xticks(
        np.arange(nbins)
    )
    ax.set_yticks(
        np.arange(nbins)
    )

    ax.set_xticklabels(
        np.arange(1, nbins + 1)
    )
    ax.set_yticklabels(
        np.arange(1, nbins + 1)
    )

    ax.set_title(title)

    fig.tight_layout()
    fig.savefig(outname)
    plt.close(fig)


def plot_sigma_comparison(results, outname):
    fig, ax = plt.subplots(
        figsize=(8, 6)
    )

    for label, result in results.items():
        ax.plot(
            np.arange(
                1,
                len(result["sigma"]) + 1,
            ),
            result["sigma"],
            marker="o",
            label=label,
        )

    ax.set_xlabel("Bin")
    ax.set_ylabel(
        "Fractional flux uncertainty"
    )

    ax.grid(
        True,
        alpha=0.25,
    )

    ax.legend()

    fig.tight_layout()
    fig.savefig(outname)
    plt.close(fig)


def plot_eigenvalue_fraction(results, outname):
    fig, ax = plt.subplots(
        figsize=(8, 6)
    )

    for label, result in results.items():
        x = np.arange(
            1,
            len(result["eigfrac"]) + 1,
        )

        ax.plot(
            x,
            result["eigfrac"],
            marker="o",
            label=label,
        )

    ax.set_xlabel(
        "Flux covariance eigenmode"
    )

    ax.set_ylabel(
        "Fraction of total variance"
    )

    ax.set_yscale("log")

    ax.grid(
        True,
        alpha=0.25,
    )

    ax.legend()

    fig.tight_layout()
    fig.savefig(outname)
    plt.close(fig)


def plot_cumulative_variance(results, outname):
    fig, ax = plt.subplots(
        figsize=(8, 6)
    )

    for label, result in results.items():
        x = np.arange(
            1,
            len(result["cumulative"]) + 1,
        )

        ax.plot(
            x,
            result["cumulative"],
            marker="o",
            label=label,
        )

    ax.axhline(
        0.90,
        linestyle="--",
        linewidth=1,
    )

    ax.axhline(
        0.95,
        linestyle="--",
        linewidth=1,
    )

    ax.axhline(
        0.99,
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

    ax.grid(
        True,
        alpha=0.25,
    )

    ax.legend()

    fig.tight_layout()
    fig.savefig(outname)
    plt.close(fig)


def main():
    args = parse_args()

    f_ratio = ROOT.TFile.Open(
        args.ratio_file,
        "READ",
    )

    if not f_ratio or f_ratio.IsZombie():
        raise RuntimeError(
            "Could not open ratio file {}".format(
                args.ratio_file
            )
        )

    f_direct = ROOT.TFile.Open(
        args.direct_file,
        "READ",
    )

    if not f_direct or f_direct.IsZombie():
        raise RuntimeError(
            "Could not open direct file {}".format(
                args.direct_file
            )
        )

    configs = {
        "FHC CC ratio": (
            f_ratio,
            "mc_fhc_ratio",
        ),
        "RHC CC ratio": (
            f_ratio,
            "mc_rhc_ratio",
        ),
        "FHC direct CCnue": (
            f_direct,
            "mc_fhc_nue_selection",
        ),
        "RHC direct CCnue": (
            f_direct,
            "mc_rhc_nue_selection",
        ),
    }

    results = {}

    print("\nRatio configuration:")
    print(" ", args.ratio_file)

    print("\nDirect CCnue configuration:")
    print(" ", args.direct_file)

    print("\nError band:")
    print(" ", args.error_band)

    for label, (root_file, hist_name) in configs.items():

        result = analyze_hist(
            root_file,
            hist_name,
            args.error_band,
        )

        results[label] = result

        print_summary(
            label,
            result,
        )

    # --------------------------------------------------
    # Matrix plots
    # --------------------------------------------------

    for label, result in results.items():
        safe_label = (
            label.lower()
            .replace(" ", "_")
            .replace("/", "_")
        )

        plot_matrix(
            result["cov"],
            "{} flux covariance".format(label),
            "{}_{}_covariance.png".format(
                args.out_prefix,
                safe_label,
            ),
            "Fractional covariance",
        )

        plot_matrix(
            result["corr"],
            "{} flux correlation".format(label),
            "{}_{}_correlation.png".format(
                args.out_prefix,
                safe_label,
            ),
            "Correlation coefficient",
            vmin=-1.0,
            vmax=1.0,
        )

    # --------------------------------------------------
    # Direct comparison plots
    # --------------------------------------------------

    plot_sigma_comparison(
        results,
        args.out_prefix + "_sigma_comparison.png",
    )

    plot_eigenvalue_fraction(
        results,
        args.out_prefix + "_eigenvalue_fraction.png",
    )

    plot_cumulative_variance(
        results,
        args.out_prefix + "_cumulative_variance.png",
    )

    # --------------------------------------------------
    # Save numerical output
    # --------------------------------------------------

    for label, result in results.items():
        safe_label = (
            label.lower()
            .replace(" ", "_")
            .replace("/", "_")
        )

        np.savetxt(
            "{}_{}_covariance.csv".format(
                args.out_prefix,
                safe_label,
            ),
            result["cov"],
            delimiter=",",
        )

        np.savetxt(
            "{}_{}_correlation.csv".format(
                args.out_prefix,
                safe_label,
            ),
            result["corr"],
            delimiter=",",
        )

        eig_table = np.column_stack([
            np.arange(
                1,
                len(result["eigvals"]) + 1,
            ),
            result["eigvals"],
            result["eigfrac"],
            result["cumulative"],
        ])

        np.savetxt(
            "{}_{}_eigenvalues.csv".format(
                args.out_prefix,
                safe_label,
            ),
            eig_table,
            delimiter=",",
            header=(
                "mode,eigenvalue,"
                "fraction,total_cumulative"
            ),
            comments="",
        )

    f_ratio.Close()
    f_direct.Close()

    print("\nSaved comparison plots:")
    print(
        " ",
        args.out_prefix + "_sigma_comparison.png",
    )
    print(
        " ",
        args.out_prefix + "_eigenvalue_fraction.png",
    )
    print(
        " ",
        args.out_prefix + "_cumulative_variance.png",
    )


if __name__ == "__main__":
    main()
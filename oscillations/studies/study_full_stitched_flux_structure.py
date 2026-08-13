#!/usr/bin/env python

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
            "Compare full stitched flux covariance/eigenstructure "
            "between ratio and direct CCnue configurations."
        )
    )

    parser.add_argument(
        "--ratio-file",
        required=True,
        help="Stitched ROOT file for ratio configuration",
    )

    parser.add_argument(
        "--direct-file",
        required=True,
        help="Stitched ROOT file for direct CCnue configuration",
    )

    parser.add_argument(
        "--hist",
        default="sample_mc",
        help="Full stitched MnvH1D [default: sample_mc]",
    )

    parser.add_argument(
        "--error-band",
        default="Flux",
        help="Vertical error band [default: Flux]",
    )

    parser.add_argument(
        "--out-prefix",
        default="full_stitched_flux_structure",
    )

    return parser.parse_args()


def get_hist(root_file, name):
    hist = root_file.Get(name)

    if not hist:
        raise RuntimeError(
            "Could not find histogram '{}'".format(name)
        )

    return hist


def get_cv(hist):
    return np.array([
        hist.GetBinContent(i)
        for i in range(1, hist.GetNbinsX() + 1)
    ], dtype=float)


def get_universes(hist, band_name):
    if not hist.HasVertErrorBand(band_name):
        raise RuntimeError(
            "Histogram '{}' has no vertical error band '{}'. "
            "Available: {}".format(
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

    # numerical cleanup only
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
    where = np.where(
        cumulative >= target
    )[0]

    if len(where) == 0:
        return len(cumulative)

    return int(where[0] + 1)


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

    sigma = np.sqrt(np.diag(cov))

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


def print_summary(label, result):
    print("\n==================================================")
    print(label)
    print("==================================================")

    print("bins =", len(result["cv"]))
    print("universes =", result["universes"].shape[0])

    print(
        "trace fractional covariance =",
        np.trace(result["cov"]),
    )

    print(
        "mean fractional sigma =",
        np.mean(result["sigma"]),
    )

    print(
        "median fractional sigma =",
        np.median(result["sigma"]),
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

    print()

    for target in [
        0.50,
        0.68,
        0.90,
        0.95,
        0.99,
    ]:
        print(
            "Modes needed for {:5.1f}% variance = {}".format(
                100.0 * target,
                modes_for_fraction(
                    result["cumulative"],
                    target,
                ),
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
        figsize=(9, 8)
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

    ax.set_xlabel("Stitched bin")
    ax.set_ylabel("Stitched bin")
    ax.set_title(title)

    fig.tight_layout()
    fig.savefig(outname)
    plt.close(fig)


def plot_sigma_comparison(
    ratio,
    direct,
    outname,
):
    fig, ax = plt.subplots(
        figsize=(10, 6)
    )

    x_ratio = np.arange(
        1,
        len(ratio["sigma"]) + 1,
    )

    x_direct = np.arange(
        1,
        len(direct["sigma"]) + 1,
    )

    ax.plot(
        x_ratio,
        ratio["sigma"],
        marker="o",
        label="Ratio configuration",
    )

    ax.plot(
        x_direct,
        direct["sigma"],
        marker="o",
        label="Direct CCnue configuration",
    )

    ax.set_xlabel("Stitched bin")
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


def plot_eigen_fraction(
    ratio,
    direct,
    outname,
):
    fig, ax = plt.subplots(
        figsize=(8, 6)
    )

    xr = np.arange(
        1,
        len(ratio["eigfrac"]) + 1,
    )

    xd = np.arange(
        1,
        len(direct["eigfrac"]) + 1,
    )

    ax.plot(
        xr,
        ratio["eigfrac"],
        marker="o",
        label="Ratio configuration",
    )

    ax.plot(
        xd,
        direct["eigfrac"],
        marker="o",
        label="Direct CCnue configuration",
    )

    ax.set_yscale("log")

    ax.set_xlabel(
        "Flux covariance eigenmode"
    )

    ax.set_ylabel(
        "Fraction of total variance"
    )

    ax.grid(
        True,
        alpha=0.25,
    )

    ax.legend()

    fig.tight_layout()
    fig.savefig(outname)
    plt.close(fig)


def plot_cumulative(
    ratio,
    direct,
    outname,
):
    fig, ax = plt.subplots(
        figsize=(8, 6)
    )

    xr = np.arange(
        1,
        len(ratio["cumulative"]) + 1,
    )

    xd = np.arange(
        1,
        len(direct["cumulative"]) + 1,
    )

    ax.plot(
        xr,
        ratio["cumulative"],
        marker="o",
        label="Ratio configuration",
    )

    ax.plot(
        xd,
        direct["cumulative"],
        marker="o",
        label="Direct CCnue configuration",
    )

    for y in [
        0.90,
        0.95,
        0.99,
    ]:
        ax.axhline(
            y,
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


def save_eigenvectors(
    result,
    prefix,
    nmodes=10,
):
    nmodes = min(
        nmodes,
        result["eigvecs"].shape[1],
    )

    x = np.arange(
        1,
        result["eigvecs"].shape[0] + 1,
    )

    for k in range(nmodes):

        fig, ax = plt.subplots(
            figsize=(10, 5)
        )

        ax.plot(
            x,
            result["eigvecs"][:, k],
            marker="o",
        )

        ax.axhline(
            0.0,
            linewidth=1,
        )

        ax.set_xlabel(
            "Stitched bin"
        )

        ax.set_ylabel(
            "Eigenvector component"
        )

        ax.set_title(
            "Mode {} — {:.2f}% of flux variance".format(
                k + 1,
                100.0 * result["eigfrac"][k],
            )
        )

        ax.grid(
            True,
            alpha=0.25,
        )

        fig.tight_layout()

        fig.savefig(
            "{}_mode_{:02d}.png".format(
                prefix,
                k + 1,
            )
        )

        plt.close(fig)


def main():
    args = parse_args()

    f_ratio = ROOT.TFile.Open(
        args.ratio_file,
        "READ",
    )

    if not f_ratio or f_ratio.IsZombie():
        raise RuntimeError(
            "Could not open {}".format(
                args.ratio_file
            )
        )

    f_direct = ROOT.TFile.Open(
        args.direct_file,
        "READ",
    )

    if not f_direct or f_direct.IsZombie():
        raise RuntimeError(
            "Could not open {}".format(
                args.direct_file
            )
        )

    ratio = analyze(
        f_ratio,
        args.hist,
        args.error_band,
    )

    direct = analyze(
        f_direct,
        args.hist,
        args.error_band,
    )

    print_summary(
        "FULL STITCHED RATIO CONFIGURATION",
        ratio,
    )

    print_summary(
        "FULL STITCHED DIRECT CCnue CONFIGURATION",
        direct,
    )

    print("\n==================================================")
    print("DIRECT COMPARISON")
    print("==================================================")

    ratio_trace = np.trace(
        ratio["cov"]
    )

    direct_trace = np.trace(
        direct["cov"]
    )

    print(
        "Ratio trace / direct trace =",
        ratio_trace / direct_trace,
    )

    print(
        "Fractional reduction in total variance =",
        1.0 - ratio_trace / direct_trace,
    )

    for target in [
        0.90,
        0.95,
        0.99,
    ]:
        nr = modes_for_fraction(
            ratio["cumulative"],
            target,
        )

        nd = modes_for_fraction(
            direct["cumulative"],
            target,
        )

        print(
            "{}% variance: ratio modes = {}, direct modes = {}".format(
                int(100 * target),
                nr,
                nd,
            )
        )

    # --------------------------------------------------
    # Matrices
    # --------------------------------------------------

    plot_matrix(
        ratio["cov"],
        "Full stitched ratio configuration\nFlux covariance",
        args.out_prefix + "_ratio_covariance.png",
        "Fractional covariance",
    )

    plot_matrix(
        direct["cov"],
        "Full stitched direct CCnue configuration\nFlux covariance",
        args.out_prefix + "_direct_covariance.png",
        "Fractional covariance",
    )

    plot_matrix(
        ratio["corr"],
        "Full stitched ratio configuration\nFlux correlation",
        args.out_prefix + "_ratio_correlation.png",
        "Correlation coefficient",
        vmin=-1.0,
        vmax=1.0,
    )

    plot_matrix(
        direct["corr"],
        "Full stitched direct CCnue configuration\nFlux correlation",
        args.out_prefix + "_direct_correlation.png",
        "Correlation coefficient",
        vmin=-1.0,
        vmax=1.0,
    )

    # --------------------------------------------------
    # Comparisons
    # --------------------------------------------------

    plot_sigma_comparison(
        ratio,
        direct,
        args.out_prefix + "_sigma_comparison.png",
    )

    plot_eigen_fraction(
        ratio,
        direct,
        args.out_prefix + "_eigenvalue_fraction.png",
    )

    plot_cumulative(
        ratio,
        direct,
        args.out_prefix + "_cumulative_variance.png",
    )

    # --------------------------------------------------
    # Leading eigenvectors
    # --------------------------------------------------

    save_eigenvectors(
        ratio,
        args.out_prefix + "_ratio",
        nmodes=10,
    )

    save_eigenvectors(
        direct,
        args.out_prefix + "_direct",
        nmodes=10,
    )

    # --------------------------------------------------
    # Numerical output
    # --------------------------------------------------

    np.savetxt(
        args.out_prefix + "_ratio_covariance.csv",
        ratio["cov"],
        delimiter=",",
    )

    np.savetxt(
        args.out_prefix + "_direct_covariance.csv",
        direct["cov"],
        delimiter=",",
    )

    np.savetxt(
        args.out_prefix + "_ratio_correlation.csv",
        ratio["corr"],
        delimiter=",",
    )

    np.savetxt(
        args.out_prefix + "_direct_correlation.csv",
        direct["corr"],
        delimiter=",",
    )

    f_ratio.Close()
    f_direct.Close()

    print("\nDone.")


if __name__ == "__main__":
    main()
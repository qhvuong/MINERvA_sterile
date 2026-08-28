#!/usr/bin/env python

import argparse
import ROOT
import PlotUtils
import numpy as np

from matplotlib import pyplot as plt


ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

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
SAMPLE_LABELS_COMPARE = [
    r"FHC $\nu+e$",
    r"FHC e-like / $\mu/e$",
    r"RHC e-like / $\mu/e$",
    r"FHC $\nu_\mu$",
    r"RHC $\bar{\nu}_\mu$",
]
def draw_sample_boundaries_1d(ax):
    for b in SAMPLE_BREAKS:
        ax.axvline(
            b + 0.5,
            color="k",
            linestyle="--",
            linewidth=1,
            alpha=0.6,
        )
def draw_sample_boundaries_2d(ax):
    for b in SAMPLE_BREAKS:
        pos = b - 0.5

        ax.axvline(
            pos,
            color="w",
            linestyle="--",
            linewidth=1,
            alpha=0.8,
        )

        ax.axhline(
            pos,
            color="w",
            linestyle="--",
            linewidth=1,
            alpha=0.8,
        )

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compare full stitched flux covariance/eigenstructure "
            "between ratio and direct CCnue configurations."
        )
    )

    # parser.add_argument(
    #     "--ratio-file",
    #     required=True,
    #     help="Stitched ROOT file for ratio configuration",
    # )

    # parser.add_argument(
    #     "--direct-file",
    #     required=True,
    #     help="Stitched ROOT file for direct CCnue configuration",
    # )

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

    parser.add_argument("--total-ratio-file", required=True)
    parser.add_argument("--total-direct-file", required=True)

    parser.add_argument("--ppfx-ratio-file", required=True)
    parser.add_argument("--ppfx-direct-file", required=True)

    parser.add_argument("--beam-ratio-file", required=True)
    parser.add_argument("--beam-direct-file", required=True)

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

    draw_sample_boundaries_2d(ax)

    fig.tight_layout()
    fig.savefig(outname)
    plt.close(fig)


def plot_sigma_comparison(
    ratio,
    direct,
    outname,
    title,
):
    fig, ax = plt.subplots(figsize=(10, 6))

    x_ratio = np.arange(1, len(ratio["sigma"]) + 1)
    x_direct = np.arange(1, len(direct["sigma"]) + 1)

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
    ax.set_ylabel("Fractional flux uncertainty")
    ax.set_title(title)

    ax.grid(True, alpha=0.25)
    draw_sample_boundaries_1d(ax)
    ax.legend()

    fig.tight_layout()
    fig.savefig(outname)
    plt.close(fig)


def plot_ratio_component_comparison(
    total,
    ppfx,
    beam,
    outname,
    sample_labels=None,
):
    fig, ax = plt.subplots(figsize=(10, 6))

    x = np.arange(1, len(total["sigma"]) + 1)

    ax.plot(
        x,
        total["sigma"],
        marker="o",
        label="Total Flux",
    )

    ax.plot(
        x,
        ppfx["sigma"],
        marker="o",
        label="PPFX only",
    )

    ax.plot(
        x,
        beam["sigma"],
        marker="o",
        label="BeamFocus only",
    )

    draw_sample_boundaries_1d(ax)

    if sample_labels is not None:
        sample_edges = [0] + SAMPLE_BREAKS + [len(total["sigma"])]

        for i, label in enumerate(sample_labels):
            center = 0.5 * (
                sample_edges[i] +
                sample_edges[i + 1]
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

    ax.set_xlabel("Stitched bin")
    ax.set_ylabel("Fractional flux uncertainty")

    ax.grid(True, alpha=0.25)
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

        draw_sample_boundaries_1d(ax)

        fig.tight_layout()

        fig.savefig(
            "{}_mode_{:02d}.png".format(
                prefix,
                k + 1,
            )
        )

        plt.close(fig)


def open_root(path):
    f = ROOT.TFile.Open(path, "READ")

    if not f or f.IsZombie():
        raise RuntimeError(
            "Could not open {}".format(path)
        )

    return f


def run_pair_study(
    label,
    ratio,
    direct,
    out_prefix,
):
    print_summary(
        "{} — RATIO CONFIGURATION".format(label),
        ratio,
    )

    print_summary(
        "{} — DIRECT CCnue CONFIGURATION".format(label),
        direct,
    )

    print("\n==================================================")
    print("{} — RATIO VS DIRECT".format(label))
    print("==================================================")

    ratio_trace = np.trace(ratio["cov"])
    direct_trace = np.trace(direct["cov"])

    print("Ratio trace =", ratio_trace)
    print("Direct trace =", direct_trace)

    if direct_trace != 0.0:
        print(
            "Ratio trace / direct trace =",
            ratio_trace / direct_trace,
        )
        print(
            "Fractional reduction in total variance =",
            1.0 - ratio_trace / direct_trace,
        )

    for target in [0.90, 0.95, 0.99]:
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

    # Fractional uncertainty
    plot_sigma_comparison(
        ratio,
        direct,
        out_prefix + "_sigma_comparison.png",
        "{}: ratio vs direct".format(label),
    )

    # Covariance matrices
    plot_matrix(
        ratio["cov"],
        "{} — ratio configuration\nFlux covariance".format(label),
        out_prefix + "_ratio_covariance.png",
        "Fractional covariance",
    )

    plot_matrix(
        direct["cov"],
        "{} — direct CCnue configuration\nFlux covariance".format(label),
        out_prefix + "_direct_covariance.png",
        "Fractional covariance",
    )

    # Correlation matrices
    plot_matrix(
        ratio["corr"],
        "{} — ratio configuration\nFlux correlation".format(label),
        out_prefix + "_ratio_correlation.png",
        "Correlation coefficient",
        vmin=-1.0,
        vmax=1.0,
    )

    plot_matrix(
        direct["corr"],
        "{} — direct CCnue configuration\nFlux correlation".format(label),
        out_prefix + "_direct_correlation.png",
        "Correlation coefficient",
        vmin=-1.0,
        vmax=1.0,
    )

    # Eigenvalue fraction
    plot_eigen_fraction(
        ratio,
        direct,
        out_prefix + "_eigenvalue_fraction.png",
    )

    # Cumulative variance
    plot_cumulative(
        ratio,
        direct,
        out_prefix + "_cumulative_variance.png",
    )

    # Leading eigenvectors
    save_eigenvectors(
        ratio,
        out_prefix + "_ratio",
        nmodes=10,
    )

    save_eigenvectors(
        direct,
        out_prefix + "_direct",
        nmodes=10,
    )

    # Numerical matrices
    np.savetxt(
        out_prefix + "_ratio_covariance.csv",
        ratio["cov"],
        delimiter=",",
    )

    np.savetxt(
        out_prefix + "_direct_covariance.csv",
        direct["cov"],
        delimiter=",",
    )

    np.savetxt(
        out_prefix + "_ratio_correlation.csv",
        ratio["corr"],
        delimiter=",",
    )

    np.savetxt(
        out_prefix + "_direct_correlation.csv",
        direct["corr"],
        delimiter=",",
    )


def print_first_fhc_bin_diagnostic(
    label,
    ratio,
    direct,
    nprint=10,
):
    # --------------------------------------------------
    # Stitched layout:
    #
    # bins  1-4   : FHC nu+e
    # bins  5-16  : FHC nue OR FHC mu/e
    # bins 29-40  : FHC numu
    #
    # numpy indices are therefore:
    #   first FHC nue / mu/e bin = 4
    #   first FHC numu bin       = 28
    # --------------------------------------------------

    idx_e = 4
    idx_ratio = 4
    idx_mu = 28

    delta_e = direct["delta"][:, idx_e]
    delta_mu = direct["delta"][:, idx_mu]
    delta_ratio = ratio["delta"][:, idx_ratio]

    # Linearized expected ratio shift:
    # delta(mu/e) ~ delta_mu - delta_e
    delta_diff = delta_mu - delta_e

    sigma_e = np.std(delta_e, ddof=1)
    sigma_mu = np.std(delta_mu, ddof=1)
    sigma_ratio = np.std(delta_ratio, ddof=1)
    sigma_diff = np.std(delta_diff, ddof=1)

    cov_mu_e = np.cov(
        delta_mu,
        delta_e,
        ddof=1,
    )[0, 1]

    rho_mu_e = np.corrcoef(
        delta_mu,
        delta_e,
    )[0, 1]

    # What sigma_ratio would be from the covariance formula
    sigma_ratio_from_cov = np.sqrt(
        sigma_mu**2
        + sigma_e**2
        - 2.0 * cov_mu_e
    )

    print("\n")
    print("============================================================")
    print("{} : FIRST FHC CC BIN".format(label))
    print("============================================================")

    print("Direct CCnue stitched bin = 5")
    print("Direct CCnumu stitched bin = 29")
    print("mu/e ratio stitched bin    = 5")
    print()

    print("CV values:")
    print(
        "  CCnue  CV = {:.8e}".format(
            direct["cv"][idx_e]
        )
    )
    print(
        "  CCnumu CV = {:.8e}".format(
            direct["cv"][idx_mu]
        )
    )
    print(
        "  mu/e   CV = {:.8e}".format(
            ratio["cv"][idx_ratio]
        )
    )

    print()
    print("Fractional sigmas:")
    print(
        "  sigma_e        = {:.8f}  ({:.4f}%)".format(
            sigma_e,
            100.0 * sigma_e,
        )
    )
    print(
        "  sigma_mu       = {:.8f}  ({:.4f}%)".format(
            sigma_mu,
            100.0 * sigma_mu,
        )
    )
    print(
        "  sigma_mu/e     = {:.8f}  ({:.4f}%)".format(
            sigma_ratio,
            100.0 * sigma_ratio,
        )
    )
    print(
        "  sigma(mu-e)    = {:.8f}  ({:.4f}%)".format(
            sigma_diff,
            100.0 * sigma_diff,
        )
    )
    print(
        "  sigma from cov = {:.8f}  ({:.4f}%)".format(
            sigma_ratio_from_cov,
            100.0 * sigma_ratio_from_cov,
        )
    )

    print()
    print("mu/e relationship:")
    print(
        "  covariance(mu,e) = {:.8e}".format(
            cov_mu_e
        )
    )
    print(
        "  rho(mu,e)        = {:.8f}".format(
            rho_mu_e
        )
    )

    print()
    print("Mean fractional shifts:")
    print(
        "  <delta_e>        = {:+.8e}".format(
            np.mean(delta_e)
        )
    )
    print(
        "  <delta_mu>       = {:+.8e}".format(
            np.mean(delta_mu)
        )
    )
    print(
        "  <delta_mu/e>     = {:+.8e}".format(
            np.mean(delta_ratio)
        )
    )
    print(
        "  <delta_mu-e>     = {:+.8e}".format(
            np.mean(delta_diff)
        )
    )

    print()
    print("First {} universes:".format(nprint))
    print(
        "{:>5s} {:>13s} {:>13s} {:>13s} {:>13s}".format(
            "u",
            "delta_mu",
            "delta_e",
            "mu-e",
            "delta_ratio",
        )
    )

    for u in range(min(nprint, len(delta_e))):
        print(
            "{:5d} {:+13.6e} {:+13.6e} {:+13.6e} {:+13.6e}".format(
                u,
                delta_mu[u],
                delta_e[u],
                delta_diff[u],
                delta_ratio[u],
            )
        )


def print_all_electron_bin_diagnostics(
    total_direct,
    ppfx_direct,
    beam_direct,
):
    """
    Print numerical flux-uncertainty information for every direct
    CCnue and CCnuebar stitched bin.

    Stitched layout:
      bins  5-16  : FHC CCnue
      bins 17-28  : RHC CCnuebar

    numpy indices:
      4-15  : FHC CCnue
      16-27 : RHC CCnuebar
    """

    samples = [
        ("FHC CCnue", 4, 15),
        ("RHC CCnuebar", 16, 27),
    ]

    for sample_name, start, end in samples:

        print("")
        print("================================================================================================================")
        print("{} : ALL BINS".format(sample_name))
        print("================================================================================================================")

        print(
            "{:>5s} {:>6s} {:>14s} "
            "{:>12s} {:>12s} {:>12s} "
            "{:>14s} {:>14s} {:>14s}".format(
                "local",
                "global",
                "CV",
                "Total sigma",
                "PPFX sigma",
                "BF sigma",
                "<delta_BF>",
                "min delta_BF",
                "max delta_BF",
            )
        )

        print("-" * 120)

        for idx in range(start, end + 1):

            local_bin = idx - start + 1
            global_bin = idx + 1

            cv = total_direct["cv"][idx]

            sigma_total = total_direct["sigma"][idx]
            sigma_ppfx = ppfx_direct["sigma"][idx]
            sigma_beam = beam_direct["sigma"][idx]

            delta_beam = beam_direct["delta"][:, idx]

            mean_beam = np.nanmean(delta_beam)
            min_beam = np.nanmin(delta_beam)
            max_beam = np.nanmax(delta_beam)

            print(
                "{:5d} {:6d} {:14.6e} "
                "{:11.5f}% {:11.5f}% {:11.5f}% "
                "{:+14.6e} {:+14.6e} {:+14.6e}".format(
                    local_bin,
                    global_bin,
                    cv,
                    100.0 * sigma_total,
                    100.0 * sigma_ppfx,
                    100.0 * sigma_beam,
                    mean_beam,
                    min_beam,
                    max_beam,
                )
            )

            cv_total = total_direct["cv"][idx]
            cv_ppfx = ppfx_direct["cv"][idx]
            cv_beam = beam_direct["cv"][idx]

            if not (
                np.isclose(cv_total, cv_ppfx, rtol=1e-10, atol=1e-12)
                and
                np.isclose(cv_total, cv_beam, rtol=1e-10, atol=1e-12)
            ):
                print(
                    "WARNING: CV mismatch in global bin {}: "
                    "total={}, ppfx={}, beam={}".format(
                        global_bin,
                        cv_total,
                        cv_ppfx,
                        cv_beam,
                    )
                )

        print("")


def main():
    args = parse_args()

    f_total_ratio = open_root(args.total_ratio_file)
    f_total_direct = open_root(args.total_direct_file)

    f_ppfx_ratio = open_root(args.ppfx_ratio_file)
    f_ppfx_direct = open_root(args.ppfx_direct_file)

    f_beam_ratio = open_root(args.beam_ratio_file)
    f_beam_direct = open_root(args.beam_direct_file)

    total_ratio = analyze(
        f_total_ratio,
        args.hist,
        args.error_band,
    )

    total_direct = analyze(
        f_total_direct,
        args.hist,
        args.error_band,
    )

    ppfx_ratio = analyze(
        f_ppfx_ratio,
        args.hist,
        args.error_band,
    )

    ppfx_direct = analyze(
        f_ppfx_direct,
        args.hist,
        args.error_band,
    )

    beam_ratio = analyze(
        f_beam_ratio,
        args.hist,
        args.error_band,
    )

    beam_direct = analyze(
        f_beam_direct,
        args.hist,
        args.error_band,
    )

    # print_first_fhc_bin_diagnostic(
    #     "TOTAL FLUX",
    #     total_ratio,
    #     total_direct,
    # )

    # print_first_fhc_bin_diagnostic(
    #     "PPFX ONLY",
    #     ppfx_ratio,
    #     ppfx_direct,
    # )

    # print_first_fhc_bin_diagnostic(
    #     "BEAMFOCUS ONLY",
    #     beam_ratio,
    #     beam_direct,
    # )

    print_all_electron_bin_diagnostics(
        total_direct,
        ppfx_direct,
        beam_direct,
    )

    # # --------------------------------------------------
    # # Complete studies for each flux source
    # # --------------------------------------------------

    # run_pair_study(
    #     "Total Flux",
    #     total_ratio,
    #     total_direct,
    #     args.out_prefix + "_total",
    # )

    # run_pair_study(
    #     "PPFX only",
    #     ppfx_ratio,
    #     ppfx_direct,
    #     args.out_prefix + "_ppfx",
    # )

    # run_pair_study(
    #     "BeamFocus only",
    #     beam_ratio,
    #     beam_direct,
    #     args.out_prefix + "_beamfocus",
    # )

    # # --------------------------------------------------
    # # Special three-way ratio comparison
    # # --------------------------------------------------

    # plot_ratio_component_comparison(
    #     total_ratio,
    #     ppfx_ratio,
    #     beam_ratio,
    #     args.out_prefix + "_ratio_total_ppfx_beamfocus.png",
    #     sample_labels=SAMPLE_LABELS_RATIO,
    # )

    # plot_ratio_component_comparison(
    #     total_direct,
    #     ppfx_direct,
    #     beam_direct,
    #     args.out_prefix + "_direct_total_ppfx_beamfocus.png",
    #     sample_labels=SAMPLE_LABELS_DIRECT,
    # )

    f_total_ratio.Close()
    f_total_direct.Close()

    f_ppfx_ratio.Close()
    f_ppfx_direct.Close()

    f_beam_ratio.Close()
    f_beam_direct.Close()

    print("\nDone.")


if __name__ == "__main__":
    main()
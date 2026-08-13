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
        description="Study flux-uncertainty cancellation in CC mu/e ratios."
    )

    parser.add_argument(
        "input",
        help="Stitched ROOT file",
    )

    parser.add_argument(
        "--mu-hist",
        required=True,
        help="CCnumu histogram used as numerator",
    )

    parser.add_argument(
        "--e-hist",
        required=True,
        help="CCnue histogram used as denominator",
    )

    parser.add_argument(
        "--ratio-hist",
        default=None,
        help=(
            "Optional stored mu/e ratio MnvH1D. "
            "If supplied, compare its universe shifts to the ratio "
            "constructed directly from mu/e."
        ),
    )

    parser.add_argument(
        "--error-band",
        default="Flux",
        help="Vertical error band to study [default: Flux]",
    )

    parser.add_argument(
        "--universes",
        nargs="*",
        type=int,
        default=None,
        help=(
            "Specific universes to print/plot in detail. "
            "Example: --universes 0 1 2 10 50"
        ),
    )

    parser.add_argument(
        "--out-prefix",
        default="ratio_flux_cancellation",
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
        uhist = band.GetHist(u)

        for i in range(1, hist.GetNbinsX() + 1):
            values[u, i - 1] = uhist.GetBinContent(i)

    return values


def fractional_shift(universe_values, cv):
    """
    universe_values shape: (Nuniv, Nbin)
    cv shape:              (Nbin,)
    """

    shift = np.full_like(universe_values, np.nan, dtype=float)

    good = cv != 0.0

    shift[:, good] = (
        universe_values[:, good] - cv[good]
    ) / cv[good]

    return shift


def correlation_per_bin(delta_mu, delta_e):
    nbins = delta_mu.shape[1]

    rho = np.full(nbins, np.nan)

    for i in range(nbins):
        x = delta_mu[:, i]
        y = delta_e[:, i]

        good = np.isfinite(x) & np.isfinite(y)

        if np.count_nonzero(good) < 2:
            continue

        if (
            np.std(x[good]) == 0.0
            or np.std(y[good]) == 0.0
        ):
            continue

        rho[i] = np.corrcoef(
            x[good],
            y[good],
        )[0, 1]

    return rho


def main():
    args = parse_args()

    f = ROOT.TFile.Open(args.input, "READ")

    if not f or f.IsZombie():
        raise RuntimeError(
            "Could not open {}".format(args.input)
        )

    h_mu = get_hist(f, args.mu_hist)
    h_e = get_hist(f, args.e_hist)

    if h_mu.GetNbinsX() != h_e.GetNbinsX():
        raise RuntimeError(
            "Muon and electron histograms have different bin counts."
        )

    print("\nInput:")
    print("  file       =", args.input)
    print("  mu hist    =", args.mu_hist)
    print("  e hist     =", args.e_hist)
    print("  error band =", args.error_band)

    # --------------------------------------------------
    # CV
    # --------------------------------------------------

    mu_cv = get_cv_values(h_mu)
    e_cv = get_cv_values(h_e)

    centers = get_bin_centers(h_mu)

    ratio_cv = np.full_like(mu_cv, np.nan)

    good = e_cv != 0.0
    ratio_cv[good] = mu_cv[good] / e_cv[good]

    # --------------------------------------------------
    # Flux universes
    # --------------------------------------------------

    mu_univ = get_universe_values(
        h_mu,
        args.error_band,
    )

    e_univ = get_universe_values(
        h_e,
        args.error_band,
    )

    if mu_univ.shape != e_univ.shape:
        raise RuntimeError(
            "Muon/electron universe arrays have different shapes: "
            "{} vs {}".format(mu_univ.shape, e_univ.shape)
        )

    nuniv, nbins = mu_univ.shape

    print("  universes  =", nuniv)
    print("  bins       =", nbins)

    # --------------------------------------------------
    # Fractional shifts
    # --------------------------------------------------

    delta_mu = fractional_shift(mu_univ, mu_cv)
    delta_e = fractional_shift(e_univ, e_cv)

    # First-order prediction
    delta_ratio_linear = delta_mu - delta_e

    # Exact universe-by-universe mu/e ratio
    ratio_univ = np.full_like(mu_univ, np.nan)

    good = e_univ != 0.0
    ratio_univ[good] = (
        mu_univ[good] / e_univ[good]
    )

    delta_ratio_exact = fractional_shift(
        ratio_univ,
        ratio_cv,
    )

    # --------------------------------------------------
    # Check stored ratio if provided
    # --------------------------------------------------

    delta_ratio_stored = None

    if args.ratio_hist is not None:
        h_ratio = get_hist(
            f,
            args.ratio_hist,
        )

        ratio_stored_cv = get_cv_values(h_ratio)

        ratio_stored_univ = get_universe_values(
            h_ratio,
            args.error_band,
        )

        delta_ratio_stored = fractional_shift(
            ratio_stored_univ,
            ratio_stored_cv,
        )

        print(
            "  stored ratio hist =",
            args.ratio_hist,
        )

    # --------------------------------------------------
    # Summary statistics over universes
    # --------------------------------------------------

    sigma_mu = np.nanstd(
        delta_mu,
        axis=0,
        ddof=1,
    )

    sigma_e = np.nanstd(
        delta_e,
        axis=0,
        ddof=1,
    )

    sigma_ratio = np.nanstd(
        delta_ratio_exact,
        axis=0,
        ddof=1,
    )

    rho = correlation_per_bin(
        delta_mu,
        delta_e,
    )

    uncorr_expectation = np.sqrt(
        sigma_mu**2 + sigma_e**2
    )

    cancellation_factor = np.divide(
        sigma_ratio,
        uncorr_expectation,
        out=np.full_like(sigma_ratio, np.nan),
        where=(uncorr_expectation != 0.0),
    )

    relative_to_larger = np.divide(
        sigma_ratio,
        np.maximum(sigma_mu, sigma_e),
        out=np.full_like(sigma_ratio, np.nan),
        where=(np.maximum(sigma_mu, sigma_e) != 0.0),
    )

    # --------------------------------------------------
    # Print summary
    # --------------------------------------------------

    print("\nPer-bin flux cancellation:")
    print(
        "{:>4s} {:>10s} {:>10s} {:>10s} {:>10s} "
        "{:>10s} {:>10s}".format(
            "bin",
            "sigma_mu",
            "sigma_e",
            "rho",
            "sigma_R",
            "C_uncorr",
            "C_max",
        )
    )

    for i in range(nbins):
        print(
            "{:4d} {:10.4f} {:10.4f} {:10.4f} "
            "{:10.4f} {:10.4f} {:10.4f}".format(
                i + 1,
                sigma_mu[i],
                sigma_e[i],
                rho[i],
                sigma_ratio[i],
                cancellation_factor[i],
                relative_to_larger[i],
            )
        )

    # --------------------------------------------------
    # Detailed selected universes
    # --------------------------------------------------

    if args.universes is not None:
        print("\nSelected universe fractional shifts:")

        for u in args.universes:
            if u < 0 or u >= nuniv:
                print(
                    "Skipping invalid universe",
                    u,
                )
                continue

            print("\nUniverse {}".format(u))
            print(
                "{:>4s} {:>10s} {:>10s} "
                "{:>12s} {:>12s}".format(
                    "bin",
                    "delta_mu",
                    "delta_e",
                    "dR exact",
                    "dmu-de",
                )
            )

            for i in range(nbins):
                print(
                    "{:4d} {:10.4f} {:10.4f} "
                    "{:12.4f} {:12.4f}".format(
                        i + 1,
                        delta_mu[u, i],
                        delta_e[u, i],
                        delta_ratio_exact[u, i],
                        delta_ratio_linear[u, i],
                    )
                )

    # --------------------------------------------------
    # Save full arrays
    # --------------------------------------------------

    np.save(
        args.out_prefix + "_delta_mu.npy",
        delta_mu,
    )

    np.save(
        args.out_prefix + "_delta_e.npy",
        delta_e,
    )

    np.save(
        args.out_prefix + "_delta_ratio_exact.npy",
        delta_ratio_exact,
    )

    # --------------------------------------------------
    # Save summary table
    # --------------------------------------------------

    summary = np.column_stack([
        np.arange(1, nbins + 1),
        centers,
        sigma_mu,
        sigma_e,
        rho,
        sigma_ratio,
        cancellation_factor,
        relative_to_larger,
    ])

    np.savetxt(
        args.out_prefix + "_summary.csv",
        summary,
        delimiter=",",
        header=(
            "bin,center,sigma_mu,sigma_e,rho_mu_e,"
            "sigma_ratio,cancellation_vs_uncorr,"
            "ratio_sigma_vs_max_flavor_sigma"
        ),
        comments="",
    )

    # --------------------------------------------------
    # Plot 1: sigma comparison
    # --------------------------------------------------

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(
        centers,
        sigma_mu,
        marker="o",
        label=r"$\sigma_\mu$",
    )

    ax.plot(
        centers,
        sigma_e,
        marker="s",
        label=r"$\sigma_e$",
    )

    ax.plot(
        centers,
        sigma_ratio,
        marker="^",
        label=r"$\sigma_{\mu/e}$",
    )

    ax.set_xlabel("Bin center")
    ax.set_ylabel("Fractional flux uncertainty")
    ax.legend()
    ax.grid(True, alpha=0.25)

    fig.tight_layout()
    fig.savefig(
        args.out_prefix + "_fractional_uncertainties.png"
    )
    plt.close(fig)

    # --------------------------------------------------
    # Plot 2: mu/e correlation
    # --------------------------------------------------

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(
        centers,
        rho,
        marker="o",
    )

    ax.axhline(0.0, linewidth=1)

    ax.set_ylim(-1.05, 1.05)
    ax.set_xlabel("Bin center")
    ax.set_ylabel(
        r"$\rho(\delta_\mu,\delta_e)$"
    )
    ax.grid(True, alpha=0.25)

    fig.tight_layout()
    fig.savefig(
        args.out_prefix + "_mu_e_correlation.png"
    )
    plt.close(fig)

    # --------------------------------------------------
    # Plot 3: cancellation factor
    # --------------------------------------------------

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(
        centers,
        cancellation_factor,
        marker="o",
    )

    ax.axhline(
        1.0,
        linewidth=1,
        linestyle="--",
    )

    ax.set_xlabel("Bin center")
    ax.set_ylabel(
        r"$\sigma_{\mu/e}/"
        r"\sqrt{\sigma_\mu^2+\sigma_e^2}$"
    )
    ax.grid(True, alpha=0.25)

    fig.tight_layout()
    fig.savefig(
        args.out_prefix + "_cancellation_factor.png"
    )
    plt.close(fig)

    # --------------------------------------------------
    # Stored-ratio closure
    # --------------------------------------------------

    if delta_ratio_stored is not None:
        diff = (
            delta_ratio_stored
            - delta_ratio_exact
        )

        print("\nStored ratio closure:")
        print(
            "  max |stored - constructed| =",
            np.nanmax(np.abs(diff)),
        )

    f.Close()

    print("\nSaved:")
    print(
        " ",
        args.out_prefix + "_summary.csv",
    )
    print(
        " ",
        args.out_prefix + "_fractional_uncertainties.png",
    )
    print(
        " ",
        args.out_prefix + "_mu_e_correlation.png",
    )
    print(
        " ",
        args.out_prefix + "_cancellation_factor.png",
    )


if __name__ == "__main__":
    main()
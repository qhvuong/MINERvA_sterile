#!/usr/bin/env python3

import argparse
import os

import ROOT
import numpy as np
from matplotlib import pyplot as plt

ROOT.TH1.AddDirectory(False)


SPECIES_LATEX = {
    "numu": r"\nu_\mu",
    "numubar": r"\bar{\nu}_\mu",
    "nue": r"\nu_e",
    "nuebar": r"\bar{\nu}_e",
}


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Compare BeamFocus uncertainty in old vs new produced flux files "
            "for multiple LE playlists and neutrino species."
        )
    )

    p.add_argument(
        "--old-dir",
        required=True,
        help="Directory containing old produced flux files",
    )

    p.add_argument(
        "--new-dir",
        required=True,
        help="Directory containing new produced flux files",
    )

    p.add_argument(
        "--playlists",
        nargs="+",
        default=["LE1", "LE5", "LE13"],
        help="Playlists to process",
    )

    p.add_argument(
        "--species",
        nargs="+",
        default=["numu", "numubar", "nue", "nuebar"],
        choices=["numu", "numubar", "nue", "nuebar"],
        help="Species to process",
    )

    p.add_argument(
        "--hist",
        default="flux_E_cvweighted",
        help="Flux histogram name",
    )

    p.add_argument(
        "--band",
        default="Flux_BeamFocus",
        help="Vertical error band name",
    )

    p.add_argument(
        "--old-label",
        default="Old",
        help="Prefix for old-production legend labels",
    )

    p.add_argument(
        "--new-label",
        default="New",
        help="Prefix for new-production legend labels",
    )

    p.add_argument(
        "--output-dir",
        default="beamfocus_old_vs_new",
        help="Directory for output plots",
    )

    p.add_argument(
        "--xmax",
        type=float,
        default=20.0,
        help="Maximum true neutrino energy shown",
    )

    p.add_argument(
        "--ymax",
        type=float,
        default=15.0,
        help="Maximum BeamFocus uncertainty shown",
    )

    return p.parse_args()


def open_hist(path, hname):
    f = ROOT.TFile.Open(path, "READ")

    if not f or f.IsZombie():
        raise RuntimeError(
            "Could not open {}".format(path)
        )

    h = f.Get(hname)

    if not h:
        keys = [
            k.GetName()
            for k in f.GetListOfKeys()
        ]

        raise RuntimeError(
            "Could not find '{}' in {}\n"
            "Available keys include:\n{}".format(
                hname,
                path,
                "\n".join(keys[:50]),
            )
        )

    h.SetDirectory(0)
    f.Close()

    return h


def native_edges(hist):
    ax = hist.GetXaxis()

    return np.array(
        [
            ax.GetBinLowEdge(i)
            for i in range(
                1,
                hist.GetNbinsX() + 1,
            )
        ]
        + [
            ax.GetBinUpEdge(
                hist.GetNbinsX()
            )
        ],
        dtype=float,
    )


def assert_same_binning(h_old, h_new, tol=1e-12):
    e_old = native_edges(h_old)
    e_new = native_edges(h_new)

    if (
        len(e_old) != len(e_new)
        or not np.allclose(
            e_old,
            e_new,
            rtol=0.0,
            atol=tol,
        )
    ):
        raise RuntimeError(
            "Old and new histograms have different binning.\n"
            "Old edges: {}\n"
            "New edges: {}".format(
                e_old,
                e_new,
            )
        )

    return e_old


def get_band(hist, band_name):
    names = [
        str(x)
        for x in hist.GetVertErrorBandNames()
    ]

    if band_name not in names:
        raise RuntimeError(
            "Histogram '{}' has no vertical error band '{}'. "
            "Available: {}".format(
                hist.GetName(),
                band_name,
                ", ".join(names),
            )
        )

    return hist.GetVertErrorBand(
        band_name
    )


def extract_band(hist, band_name):
    nbins = hist.GetNbinsX()

    cv = np.array(
        [
            hist.GetBinContent(i)
            for i in range(
                1,
                nbins + 1,
            )
        ],
        dtype=float,
    )

    band = get_band(
        hist,
        band_name,
    )

    nuniv = band.GetNHists()

    universes = np.zeros(
        (nbins, nuniv),
        dtype=float,
    )

    for u in range(nuniv):
        hu = band.GetHist(u)

        universes[:, u] = [
            hu.GetBinContent(i)
            for i in range(
                1,
                nbins + 1,
            )
        ]

    return cv, universes


def fractional_shifts(cv, universes):
    delta = np.full_like(
        universes,
        np.nan,
        dtype=float,
    )

    good = (
        np.isfinite(cv)
        & (cv != 0.0)
    )

    delta[good, :] = (
        universes[good, :]
        / cv[good, np.newaxis]
        - 1.0
    )

    return delta


def sigma_cv(delta):
    """
    RMS about the CV:
        sqrt(<delta^2>)
    """

    nbins = delta.shape[0]

    sigma = np.full(
        nbins,
        np.nan,
        dtype=float,
    )

    for i in range(nbins):
        x = delta[i]
        x = x[np.isfinite(x)]

        if len(x) == 0:
            continue

        sigma[i] = np.sqrt(
            np.mean(x * x)
        )

    return sigma


def sigma_centered(delta):
    """
    Centered universe spread:
        std(delta, ddof=1)
    """

    nbins = delta.shape[0]

    sigma = np.full(
        nbins,
        np.nan,
        dtype=float,
    )

    for i in range(nbins):
        x = delta[i]
        x = x[np.isfinite(x)]

        if len(x) < 2:
            continue

        sigma[i] = np.std(
            x,
            ddof=1,
        )

    return sigma


def process_one(
    old_path,
    new_path,
    playlist,
    species,
    args,
):
    latex = SPECIES_LATEX[species]

    print()
    print("=" * 90)
    print(
        "{}  {}".format(
            playlist,
            species,
        )
    )
    print("=" * 90)

    print("Old:", old_path)
    print("New:", new_path)

    h_old = open_hist(
        old_path,
        args.hist,
    )

    h_new = open_hist(
        new_path,
        args.hist,
    )

    edges = assert_same_binning(
        h_old,
        h_new,
    )

    centers = 0.5 * (
        edges[:-1]
        + edges[1:]
    )

    cv_old, u_old = extract_band(
        h_old,
        args.band,
    )

    cv_new, u_new = extract_band(
        h_new,
        args.band,
    )

    print(
        "Old universes:",
        u_old.shape[1],
    )

    print(
        "New universes:",
        u_new.shape[1],
    )

    d_old = fractional_shifts(
        cv_old,
        u_old,
    )

    d_new = fractional_shifts(
        cv_new,
        u_new,
    )

    old_cv = (
        100.0
        * sigma_cv(d_old)
    )

    new_cv = (
        100.0
        * sigma_cv(d_new)
    )

    old_centered = (
        100.0
        * sigma_centered(d_old)
    )

    new_centered = (
        100.0
        * sigma_centered(d_new)
    )

    # --------------------------------------------------------
    # Print table
    # --------------------------------------------------------

    print()

    print(
        "{:>15s} "
        "{:>11s} "
        "{:>11s} "
        "{:>11s} "
        "{:>11s}".format(
            "Etrue [GeV]",
            "old CV",
            "new CV",
            "old ctr",
            "new ctr",
        )
    )

    for i in range(
        len(centers)
    ):
        print(
            "{:6.2f}-{:6.2f} "
            "{:10.4f}% "
            "{:10.4f}% "
            "{:10.4f}% "
            "{:10.4f}%".format(
                edges[i],
                edges[i + 1],
                old_cv[i],
                new_cv[i],
                old_centered[i],
                new_centered[i],
            )
        )

    # --------------------------------------------------------
    # Plot centered uncertainty
    # --------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(9, 6)
    )

    ax.plot(
        centers,
        old_centered,
        marker="o",
        label="{} {} ${}$".format(
            args.old_label,
            playlist,
            latex,
        ),
    )

    ax.plot(
        centers,
        new_centered,
        marker="o",
        label="{} {} ${}$".format(
            args.new_label,
            playlist,
            latex,
        ),
    )

    ax.set_xlim(
        0,
        args.xmax,
    )

    ax.set_ylim(
        0,
        args.ymax,
    )

    ax.set_xlabel(
        r"True $E_\nu$ [GeV]"
    )

    ax.set_ylabel(
        "BeamFocus centered uncertainty [%]"
    )

    ax.set_title(
        r"BeamFocus uncertainty: {} ${}$".format(
            playlist,
            latex,
        )
    )

    ax.grid(
        True,
        alpha=0.25,
    )

    ax.legend()

    fig.tight_layout()

    output = os.path.join(
        args.output_dir,
        "beamfocus_old_vs_new_{}_{}.png".format(
            playlist,
            species,
        ),
    )

    fig.savefig(
        output,
        dpi=160,
    )

    plt.close(fig)

    print()
    print("Saved:", output)


def main():
    args = parse_args()

    os.makedirs(
        args.output_dir,
        exist_ok=True,
    )

    n_done = 0
    n_skipped = 0

    for playlist in args.playlists:

        for species in args.species:

            filename = "{}_{}.root".format(
                playlist,
                species,
            )

            old_path = os.path.join(
                args.old_dir,
                filename,
            )

            new_path = os.path.join(
                args.new_dir,
                filename,
            )

            # ------------------------------------------------
            # Check both inputs exist.
            # Skip rather than killing the entire batch.
            # ------------------------------------------------

            if not os.path.isfile(old_path):
                print()
                print(
                    "[SKIP] Missing old file:",
                    old_path,
                )

                n_skipped += 1
                continue

            if not os.path.isfile(new_path):
                print()
                print(
                    "[SKIP] Missing new file:",
                    new_path,
                )

                n_skipped += 1
                continue

            process_one(
                old_path,
                new_path,
                playlist,
                species,
                args,
            )

            n_done += 1

    print()
    print("=" * 90)
    print("Finished")
    print("=" * 90)
    print("Processed :", n_done)
    print("Skipped   :", n_skipped)
    print("Output dir:", args.output_dir)


if __name__ == "__main__":
    main()
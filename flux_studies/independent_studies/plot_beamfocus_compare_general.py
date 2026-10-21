#!/usr/bin/env python3

import argparse
import numpy as np
from matplotlib import pyplot as plt


SPECIES_INDEX = {
    "numu": 0,
    "nue": 1,
    "numubar": 2,
    "nuebar": 3,
}

SPECIES_LATEX = {
    "numu": r"\nu_\mu",
    "nue": r"\nu_e",
    "numubar": r"\bar{\nu}_\mu",
    "nuebar": r"\bar{\nu}_e",
}


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Compare BeamFocus true-E uncertainty shapes between two NPZ files, "
            "typically ME and LE."
        )
    )

    p.add_argument(
        "--me-file",
        required=True,
        help="Input NPZ for ME",
    )

    p.add_argument(
        "--le-file",
        required=True,
        help="Input NPZ for LE",
    )

    p.add_argument(
        "--species",
        required=True,
        choices=["numu", "nue", "numubar", "nuebar"],
        help="Species to compare",
    )

    p.add_argument(
        "--mode",
        required=True,
        choices=["FHC", "RHC", "fhc", "rhc"],
        help="Horn current mode used in plot labels",
    )

    p.add_argument(
        "--definition",
        default="centered",
        choices=["centered", "cv"],
        help="Use centered or cv uncertainty definition",
    )

    p.add_argument(
        "--binning",
        default="LE",
        choices=["LE", "ME", "le", "me"],
        help="Binning key to read from the NPZ files",
    )

    p.add_argument(
        "--scale",
        type=float,
        default=0.50,
        help=(
            "Energy-axis scale applied to ME: "
            "E_target = scale * E_ME"
        ),
    )

    p.add_argument(
        "--me-label",
        default=None,
        help="Optional override for ME legend label",
    )

    p.add_argument(
        "--le-label",
        default=None,
        help="Optional override for LE legend label",
    )

    p.add_argument(
        "--output-prefix",
        default=None,
        help="Prefix for output filenames",
    )

    p.add_argument(
        "--xmax-raw",
        type=float,
        default=15.0,
        help="x max for raw comparison plot",
    )

    p.add_argument(
        "--xmax-scaled",
        type=float,
        default=8.0,
        help="x max for scaled comparison plot",
    )

    p.add_argument(
        "--ymax",
        type=float,
        default=15.0,
        help="y max for plots",
    )

    p.add_argument(
        "--feature-energies",
        type=float,
        nargs="*",
        default=[],
        help="Optional ME feature energies to print after scaling",
    )

    return p.parse_args()


def load_species(path, species, binning, definition):
    d = np.load(path)

    key = binning.lower()
    sigma_key = "{}_sigma_{}".format(
        key,
        definition,
    )
    edges_key = "{}_bin_edges".format(
        key
    )

    if edges_key not in d:
        raise RuntimeError(
            "Missing key '{}' in {}\nAvailable keys:\n{}".format(
                edges_key,
                path,
                "\n".join(sorted(d.files)),
            )
        )

    if sigma_key not in d:
        raise RuntimeError(
            "Missing key '{}' in {}\nAvailable keys:\n{}".format(
                sigma_key,
                path,
                "\n".join(sorted(d.files)),
            )
        )

    edges = d[edges_key]
    centers = 0.5 * (
        edges[:-1] + edges[1:]
    )

    sigma_all = d[sigma_key]

    ispecies = SPECIES_INDEX[species]

    sigma = 100.0 * sigma_all[ispecies]

    return edges, centers, sigma


def main():
    args = parse_args()

    species = args.species
    mode = args.mode.upper()

    latex = SPECIES_LATEX[species]
    species_label = r"${}$".format(latex)

    # --------------------------------------------------------
    # Default legend labels
    # --------------------------------------------------------

    me_label = args.me_label or "{} {} {}".format(
        "ME",
        mode,
        species_label,
    )

    le_label = args.le_label or "{} {} {}".format(
        "LE",
        mode,
        species_label,
    )

    # --------------------------------------------------------
    # Output prefix
    # --------------------------------------------------------

    if args.output_prefix is None:
        output_prefix = "beamfocus_{}_ME_vs_LE_{}".format(
            species,
            mode,
        )
    else:
        output_prefix = args.output_prefix

    scale_tag = "{:.2f}".format(
        args.scale
    ).replace(".", "p")

    # --------------------------------------------------------
    # Load inputs
    # --------------------------------------------------------

    me_edges, me_centers, me_sigma = load_species(
        args.me_file,
        species,
        args.binning,
        args.definition,
    )

    le_edges, le_centers, le_sigma = load_species(
        args.le_file,
        species,
        args.binning,
        args.definition,
    )

    # --------------------------------------------------------
    # Labels
    # --------------------------------------------------------

    if args.definition == "centered":
        ylabel = "BeamFocus centered uncertainty [%]"
        def_label = "centered"
    else:
        ylabel = "BeamFocus CV-RMS uncertainty [%]"
        def_label = "CV-RMS"

    # --------------------------------------------------------
    # Plot 1: raw comparison
    #
    # IMPORTANT:
    # Plot ME first and LE second so matplotlib's automatic
    # colors remain consistent with the scaled plot.
    # --------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(9, 6)
    )

    ax.plot(
        me_centers,
        me_sigma,
        marker="o",
        label=me_label,
    )

    ax.plot(
        le_centers,
        le_sigma,
        marker="o",
        label=le_label,
    )

    ax.set_xlim(
        0,
        args.xmax_raw,
    )

    ax.set_ylim(
        0,
        args.ymax,
    )

    ax.set_xlabel(
        r"True $E_\nu$ [GeV]"
    )

    ax.set_ylabel(
        ylabel
    )

    ax.set_title(
        r"BeamFocus uncertainty comparison for ${}$".format(
            latex
        )
        + "\n{} binning, {}".format(
            args.binning.upper(),
            def_label,
        )
    )

    ax.grid(
        True,
        alpha=0.25,
    )

    ax.legend()

    fig.tight_layout()

    raw_out = (
        output_prefix
        + "_raw.png"
    )

    fig.savefig(
        raw_out,
        dpi=160,
    )

    plt.close(fig)

    # --------------------------------------------------------
    # Plot 2: scaled ME comparison
    #
    # E_LE = scale * E_ME
    #
    # Again plot ME first, then LE, so the colors are identical
    # between the raw and scaled figures.
    # --------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(9, 6)
    )

    scaled_me_label = (
        me_label
        + ", "
        + r"$E_\nu \rightarrow "
        + "{:.2f}".format(
            args.scale
        )
        + r"E_\nu$"
    )

    ax.plot(
        args.scale * me_centers,
        me_sigma,
        marker="o",
        label=scaled_me_label,
    )

    ax.plot(
        le_centers,
        le_sigma,
        marker="o",
        label=le_label,
    )

    ax.set_xlim(
        0,
        args.xmax_scaled,
    )

    ax.set_ylim(
        0,
        args.ymax,
    )

    ax.set_xlabel(
        r"Rescaled true $E_\nu$ [GeV]"
    )

    ax.set_ylabel(
        ylabel
    )

    ax.set_title(
        r"BeamFocus structure mapping for ${}$".format(
            latex
        )
        + "\n{} binning, {}".format(
            args.binning.upper(),
            def_label,
        )
    )

    ax.grid(
        True,
        alpha=0.25,
    )

    ax.legend()

    fig.tight_layout()

    scaled_out = (
        output_prefix
        + "_scaled_{}.png".format(
            scale_tag
        )
    )

    fig.savefig(
        scaled_out,
        dpi=160,
    )

    plt.close(fig)

    # --------------------------------------------------------
    # Print summary
    # --------------------------------------------------------

    print()
    print("Species        :", species)
    print("Mode           :", mode)
    print("Definition     :", args.definition)
    print("Binning        :", args.binning.upper())
    print("Scale          :", args.scale)
    print("ME input       :", args.me_file)
    print("LE input       :", args.le_file)
    print("ME label       :", me_label)
    print("LE label       :", le_label)

    if args.feature_energies:
        print()
        print("Approximate feature mapping:")

        for e_me in args.feature_energies:
            print(
                "  ME {:6.2f} GeV -> target {:6.2f} GeV".format(
                    e_me,
                    args.scale * e_me,
                )
            )

    print()
    print("Saved:")
    print("  {}".format(raw_out))
    print("  {}".format(scaled_out))


if __name__ == "__main__":
    main()
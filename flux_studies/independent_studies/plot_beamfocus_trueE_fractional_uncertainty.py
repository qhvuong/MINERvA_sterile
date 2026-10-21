#!/usr/bin/env python

import argparse
import numpy as np
from matplotlib import pyplot as plt


SPECIES_INDEX = {
    "numu": 0,
    "nue": 1,
    "numubar": 2,
    "nuebar": 3,
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot ME BeamFocus fractional uncertainty versus true E_nu "
            "from the full-stat BeamFocus NPZ output."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="NPZ from study_beamfocus_trueE_all_species.py",
    )

    parser.add_argument(
        "--binning",
        choices=["ME", "LE", "me", "le"],
        default="LE",
        help="True-E binning to plot. Default: LE",
    )

    parser.add_argument(
        "--definition",
        choices=["cv", "centered", "both"],
        default="both",
        help=(
            "Uncertainty definition: "
            "'cv' = RMS relative to nominal CV=1; "
            "'centered' = std about universe mean; "
            "'both' = overlay both on the same plot. "
            "Default: both"
        ),
    )

    parser.add_argument(
        "--output",
        default="beamfocus_trueE_fractional_uncertainty.png",
        help="Output figure",
    )

    parser.add_argument(
        "--xmax",
        type=float,
        default=None,
        help="Optional maximum true neutrino energy to display",
    )

    parser.add_argument(
        "--ymax",
        type=float,
        default=None,
        help="Optional maximum fractional uncertainty in percent",
    )

    return parser.parse_args()


def require_key(data, key):
    if key not in data:
        raise RuntimeError(
            "Required NPZ key '{}' not found.\n"
            "Available keys:\n{}".format(
                key,
                "\n".join(sorted(data.files)),
            )
        )

    return data[key]


def load_definition(data, key, definition):

    # ------------------------------------------------------------
    # Case 1:
    # Final processed NPZ already contains sigma arrays.
    # ------------------------------------------------------------

    sigma_key = "{}_sigma_{}".format(
        key,
        definition,
    )

    if sigma_key in data:
        sigma_species = data[sigma_key]

        sigma_numu = sigma_species[
            SPECIES_INDEX["numu"]
        ]

        sigma_nue = sigma_species[
            SPECIES_INDEX["nue"]
        ]

        sigma_numubar = sigma_species[
            SPECIES_INDEX["numubar"]
        ]

        sigma_nuebar = sigma_species[
            SPECIES_INDEX["nuebar"]
        ]

        sigma_numu_over_nue = require_key(
            data,
            "{}_sigma_{}_numu_over_nue".format(
                key,
                definition,
            ),
        )

        sigma_numubar_over_nuebar = require_key(
            data,
            "{}_sigma_{}_numubar_over_nuebar".format(
                key,
                definition,
            ),
        )

        rho_numu_nue = require_key(
            data,
            "{}_rho_{}_numu_over_nue".format(
                key,
                definition,
            ),
        )

        rho_numubar_nuebar = require_key(
            data,
            "{}_rho_{}_numubar_over_nuebar".format(
                key,
                definition,
            ),
        )

    # ------------------------------------------------------------
    # Case 2:
    # Accumulator checkpoint.
    #
    # Reconstruct:
    #
    #   weights[species, bin, universe]
    #
    # from
    #
    #   universe_sums / cv_sums
    # ------------------------------------------------------------

    else:
        cv_sums = require_key(
            data,
            "{}_cv_sums".format(key),
        )

        universe_sums = require_key(
            data,
            "{}_universe_sums".format(key),
        )

        nspecies, nbins, nuniv = universe_sums.shape

        weights = np.full(
            universe_sums.shape,
            np.nan,
            dtype=np.float64,
        )

        for ispecies in range(nspecies):
            for ibin in range(nbins):

                cv = cv_sums[
                    ispecies,
                    ibin,
                ]

                if not np.isfinite(cv):
                    continue

                if cv <= 0.0:
                    continue

                weights[
                    ispecies,
                    ibin,
                    :,
                ] = (
                    universe_sums[
                        ispecies,
                        ibin,
                        :,
                    ]
                    / cv
                )

        # --------------------------------------------------------
        # Species uncertainties
        # --------------------------------------------------------

        sigma_species = np.full(
            (nspecies, nbins),
            np.nan,
            dtype=np.float64,
        )

        for ispecies in range(nspecies):
            for ibin in range(nbins):

                w = weights[
                    ispecies,
                    ibin,
                    :,
                ]

                valid = np.isfinite(w)

                if not np.any(valid):
                    continue

                w = w[valid]

                if definition == "cv":
                    delta = w - 1.0

                    sigma_species[
                        ispecies,
                        ibin,
                    ] = np.sqrt(
                        np.mean(
                            delta * delta
                        )
                    )

                elif definition == "centered":
                    if len(w) > 1:
                        sigma_species[
                            ispecies,
                            ibin,
                        ] = np.std(
                            w,
                            ddof=1,
                        )

        sigma_numu = sigma_species[
            SPECIES_INDEX["numu"]
        ]

        sigma_nue = sigma_species[
            SPECIES_INDEX["nue"]
        ]

        sigma_numubar = sigma_species[
            SPECIES_INDEX["numubar"]
        ]

        sigma_nuebar = sigma_species[
            SPECIES_INDEX["nuebar"]
        ]

        # --------------------------------------------------------
        # Same-sign ratio diagnostics
        # --------------------------------------------------------

        def calculate_pair(
            numerator_index,
            denominator_index,
        ):
            sigma_ratio = np.full(
                nbins,
                np.nan,
                dtype=np.float64,
            )

            rho = np.full(
                nbins,
                np.nan,
                dtype=np.float64,
            )

            for ibin in range(nbins):

                num = weights[
                    numerator_index,
                    ibin,
                    :,
                ]

                den = weights[
                    denominator_index,
                    ibin,
                    :,
                ]

                valid = (
                    np.isfinite(num)
                    & np.isfinite(den)
                    & (np.abs(den) > 1.0e-12)
                )

                if not np.any(valid):
                    continue

                num = num[valid]
                den = den[valid]

                ratio = num / den

                if definition == "cv":
                    delta_ratio = ratio - 1.0

                    sigma_ratio[ibin] = np.sqrt(
                        np.mean(
                            delta_ratio
                            * delta_ratio
                        )
                    )

                    dnum = num - 1.0
                    dden = den - 1.0

                    snum = np.sqrt(
                        np.mean(
                            dnum * dnum
                        )
                    )

                    sden = np.sqrt(
                        np.mean(
                            dden * dden
                        )
                    )

                    if (
                        snum > 0.0
                        and sden > 0.0
                    ):
                        rho[ibin] = (
                            np.mean(
                                dnum * dden
                            )
                            / (
                                snum
                                * sden
                            )
                        )

                elif definition == "centered":
                    if len(ratio) > 1:
                        sigma_ratio[ibin] = np.std(
                            ratio,
                            ddof=1,
                        )

                        if (
                            np.std(num, ddof=1) > 0.0
                            and
                            np.std(den, ddof=1) > 0.0
                        ):
                            rho[ibin] = np.corrcoef(
                                num,
                                den,
                            )[0, 1]

            return sigma_ratio, rho

        (
            sigma_numu_over_nue,
            rho_numu_nue,
        ) = calculate_pair(
            SPECIES_INDEX["numu"],
            SPECIES_INDEX["nue"],
        )

        (
            sigma_numubar_over_nuebar,
            rho_numubar_nuebar,
        ) = calculate_pair(
            SPECIES_INDEX["numubar"],
            SPECIES_INDEX["nuebar"],
        )

    # ------------------------------------------------------------
    # Convert fractional uncertainties -> percent
    # ------------------------------------------------------------

    return {
        "numu": 100.0 * sigma_numu,
        "nue": 100.0 * sigma_nue,
        "numubar": 100.0 * sigma_numubar,
        "nuebar": 100.0 * sigma_nuebar,

        "numu_over_nue": (
            100.0
            * sigma_numu_over_nue
        ),

        "numubar_over_nuebar": (
            100.0
            * sigma_numubar_over_nuebar
        ),

        "rho_numu_nue": (
            rho_numu_nue
        ),

        "rho_numubar_nuebar": (
            rho_numubar_nuebar
        ),
    }


def print_definition_table(
    definition,
    edges,
    centers,
    values,
):
    print()
    print(
        "Definition : {}".format(
            definition
        )
    )

    print()
    print(
        "{:>15s} "
        "{:>10s} {:>10s} {:>10s} {:>9s} | "
        "{:>10s} {:>10s} {:>10s} {:>9s}".format(
            "E_true [GeV]",
            "numu",
            "nue",
            "numu/nue",
            "rho",
            "numubar",
            "nuebar",
            "mubar/ebar",
            "rho",
        )
    )

    for i in range(len(centers)):
        print(
            "{:6.2f}-{:6.2f} "
            "{:9.3f}% "
            "{:9.3f}% "
            "{:9.3f}% "
            "{:9.4f} | "
            "{:9.3f}% "
            "{:9.3f}% "
            "{:9.3f}% "
            "{:9.4f}".format(
                edges[i],
                edges[i + 1],

                values["numu"][i],
                values["nue"][i],
                values["numu_over_nue"][i],
                values["rho_numu_nue"][i],

                values["numubar"][i],
                values["nuebar"][i],
                values["numubar_over_nuebar"][i],
                values["rho_numubar_nuebar"][i],
            )
        )


def main():
    args = parse_args()

    data = np.load(args.input)

    key = args.binning.lower()

    # ------------------------------------------------------------
    # Read binning
    # ------------------------------------------------------------

    edges = require_key(
        data,
        "{}_bin_edges".format(key),
    )

    centers = 0.5 * (
        edges[:-1]
        + edges[1:]
    )

    # ------------------------------------------------------------
    # Load requested definitions
    # ------------------------------------------------------------

    definitions = []

    if args.definition in ("cv", "both"):
        definitions.append("cv")

    if args.definition in ("centered", "both"):
        definitions.append("centered")

    results = {}

    for definition in definitions:
        results[definition] = load_definition(
            data,
            key,
            definition,
        )

    # ------------------------------------------------------------
    # Print numbers
    # ------------------------------------------------------------

    print()
    print(
        "Input      : {}".format(
            args.input
        )
    )

    print(
        "Binning    : {}".format(
            args.binning.upper()
        )
    )

    for definition in definitions:
        print_definition_table(
            definition,
            edges,
            centers,
            results[definition],
        )

    # ------------------------------------------------------------
    # Plot styling
    # ------------------------------------------------------------

    colors = {
        "mu": "C0",
        "e": "C1",
        "ratio": "C2",
    }

    line_styles = {
        "cv": "-",
        "centered": "--",
    }

    markers = {
        "cv": "o",
        "centered": None,
    }

    label_suffix = {
        "cv": "CV-RMS",
        "centered": "centered",
    }

    if args.definition == "both":
        definition_label = "CV-RMS and mean-centered std"
    elif args.definition == "cv":
        definition_label = "RMS relative to CV"
    else:
        definition_label = "mean-centered standard deviation"


    # ------------------------------------------------------------
    # Build separate output filenames
    # ------------------------------------------------------------

    if args.output.lower().endswith(".png"):
        output_base = args.output[:-4]
    else:
        output_base = args.output

    output_nu = output_base + "_numu_nue.png"
    output_nubar = output_base + "_numubar_nuebar.png"


    # ============================================================
    # Plot 1: numu / nue
    # ============================================================

    fig, ax = plt.subplots(
        figsize=(9, 6),
    )

    for definition in definitions:
        r = results[definition]

        ax.plot(
            centers,
            r["numu"],
            marker="o",
            color=colors["mu"],
            label=r"true $\nu_\mu$",
        )

        ax.plot(
            centers,
            r["nue"],
            marker="o",
            color=colors["e"],
            label=r"true $\nu_e$",
        )

        ax.plot(
            centers,
            r["numu_over_nue"],
            marker="o",
            linewidth=2.5,
            color=colors["ratio"],
            label=r"$\nu_\mu/\nu_e$",
        )

    ax.set_xlabel(
        r"True $E_\nu$ [GeV]"
    )

    ax.set_ylabel(
        "BeamFocus fractional uncertainty [%]"
    )

    ax.set_title(
        "ME BeamFocus uncertainty vs true neutrino energy\n"
        + "{} true-E binning, {}".format(
            args.binning.upper(),
            definition_label,
        )
    )

    ax.grid(
        True,
        alpha=0.25,
    )

    ax.legend(
        fontsize=9,
    )

    if args.xmax is not None:
        ax.set_xlim(
            edges[0],
            args.xmax,
        )

    if args.ymax is not None:
        ax.set_ylim(
            0.0,
            args.ymax,
        )

    fig.tight_layout()

    fig.savefig(
        output_nu,
        dpi=160,
    )

    plt.close(fig)

    print()
    print(
        "Saved:",
        output_nu,
    )


    # ============================================================
    # Plot 2: numubar / nuebar
    # ============================================================

    fig, ax = plt.subplots(
        figsize=(9, 6),
    )

    for definition in definitions:
        r = results[definition]

        ax.plot(
            centers,
            r["numubar"],
            marker="o",
            color=colors["mu"],
            label=r"true $\bar{\nu}_\mu$",
        )

        ax.plot(
            centers,
            r["nuebar"],
            marker="o",
            color=colors["e"],
            label=r"true $\bar{\nu}_e$",
        )

        ax.plot(
            centers,
            r["numubar_over_nuebar"],
            marker="o",
            linewidth=2.5,
            color=colors["ratio"],
            label=r"$\bar{\nu}_\mu/\bar{\nu}_e$",
        )

    ax.set_xlabel(
        r"True $E_\nu$ [GeV]"
    )

    ax.set_ylabel(
        "BeamFocus fractional uncertainty [%]"
    )

    ax.set_title(
        "ME BeamFocus uncertainty vs true neutrino energy\n"
        + "{} true-E binning, {}".format(
            args.binning.upper(),
            definition_label,
        )
    )

    ax.grid(
        True,
        alpha=0.25,
    )

    ax.legend(
        fontsize=9,
    )

    if args.xmax is not None:
        ax.set_xlim(
            edges[0],
            args.xmax,
        )

    if args.ymax is not None:
        ax.set_ylim(
            0.0,
            args.ymax,
        )

    fig.tight_layout()

    fig.savefig(
        output_nubar,
        dpi=160,
    )

    plt.close(fig)

    print(
        "Saved:",
        output_nubar,
    )


if __name__ == "__main__":
    main()
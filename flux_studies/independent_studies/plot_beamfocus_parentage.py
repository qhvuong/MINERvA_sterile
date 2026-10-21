#!/usr/bin/env python3

import argparse
import csv
import os

import numpy as np
from matplotlib import pyplot as plt


SPECIES_LATEX = {
    "numu": r"\nu_\mu",
    "nue": r"\nu_e",
    "numubar": r"\bar{\nu}_\mu",
    "nuebar": r"\bar{\nu}_e",
}

PARENT_LABELS = {
    "all": "All parents",
    "pion": r"$\pi^\pm$ parent",
    "muon": r"$\mu^\pm$ parent",
    "kaon_charged": r"$K^\pm$ parent",
    "kaon_neutral": r"$K^0$ parent",
    "other": "Other parent",
}


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Plot BeamFocus uncertainty versus true neutrino energy "
            "split by immediate neutrino parent category."
        )
    )

    p.add_argument("--input", required=True)
    p.add_argument(
        "--species",
        required=True,
        choices=["numu", "nue", "numubar", "nuebar"],
    )
    p.add_argument(
        "--definition",
        default="centered",
        choices=["centered", "cv"],
    )
    p.add_argument(
        "--parents",
        nargs="+",
        default=[
            "all",
            "pion",
            "muon",
            "kaon_charged",
            "kaon_neutral",
            "other",
        ],
    )
    p.add_argument("--xmin", type=float, default=0.0)
    p.add_argument("--xmax", type=float, default=20.0)
    p.add_argument("--ymax", type=float, default=None)
    p.add_argument("--min-events", type=int, default=0)
    p.add_argument("--title", default=None)
    p.add_argument("--output", default=None)
    p.add_argument("--print-table", action="store_true")

    return p.parse_args()


def load_csv(path):
    rows = []

    with open(path, "r") as f:
        reader = csv.DictReader(f)

        required = {
            "E_low",
            "E_high",
            "species",
            "parent",
            "N",
            "sigma_cv",
            "sigma_centered",
        }

        missing = required - set(reader.fieldnames)

        if missing:
            raise RuntimeError(
                "Missing required columns: {}\nAvailable columns: {}".format(
                    ", ".join(sorted(missing)),
                    ", ".join(reader.fieldnames),
                )
            )

        for row in reader:
            rows.append(
                {
                    "E_low": float(row["E_low"]),
                    "E_high": float(row["E_high"]),
                    "species": row["species"],
                    "parent": row["parent"],
                    "N": int(float(row["N"])),
                    "sigma_cv": float(row["sigma_cv"]),
                    "sigma_centered": float(row["sigma_centered"]),
                }
            )

    return rows


def main():
    args = parse_args()

    rows = load_csv(args.input)

    sigma_key = "sigma_{}".format(args.definition)

    available_species = sorted(
        set(row["species"] for row in rows)
    )

    available_parents = sorted(
        set(row["parent"] for row in rows)
    )

    if args.species not in available_species:
        raise RuntimeError(
            "Species '{}' not found.\nAvailable species: {}".format(
                args.species,
                ", ".join(available_species),
            )
        )

    bad_parents = [
        p for p in args.parents
        if p not in available_parents
    ]

    if bad_parents:
        raise RuntimeError(
            "Parent categories not found: {}\nAvailable parents: {}".format(
                ", ".join(bad_parents),
                ", ".join(available_parents),
            )
        )

    fig, ax = plt.subplots(figsize=(10, 6.5))

    plotted_rows = {}

    for parent in args.parents:

        sub = []

        for row in rows:

            if row["species"] != args.species:
                continue

            if row["parent"] != parent:
                continue

            center = 0.5 * (
                row["E_low"] + row["E_high"]
            )

            if center < args.xmin:
                continue

            if center > args.xmax:
                continue

            if row["N"] < args.min_events:
                continue

            sigma_percent = 100.0 * row[sigma_key]

            if not np.isfinite(sigma_percent):
                continue

            sub.append(
                {
                    "E_low": row["E_low"],
                    "E_high": row["E_high"],
                    "E_center": center,
                    "N": row["N"],
                    "sigma_percent": sigma_percent,
                }
            )

        sub.sort(key=lambda x: x["E_center"])

        plotted_rows[parent] = sub

        if len(sub) == 0:
            print(
                "[WARN] No usable points for parent '{}'".format(
                    parent
                )
            )
            continue

        x = [
            row["E_center"]
            for row in sub
        ]

        y = [
            row["sigma_percent"]
            for row in sub
        ]

        ax.plot(
            x,
            y,
            marker="o",
            label=PARENT_LABELS.get(
                parent,
                parent,
            ),
        )

    latex = SPECIES_LATEX[args.species]

    if args.title is None:
        title = (
            r"BeamFocus uncertainty by parent for ${}$".format(
                latex
            )
            + "\n{} uncertainty".format(
                args.definition
            )
        )
    else:
        title = args.title

    ax.set_title(title)
    ax.set_xlabel(r"True $E_\nu$ [GeV]")

    if args.definition == "centered":
        ax.set_ylabel(
            "BeamFocus centered fractional uncertainty [%]"
        )
    else:
        ax.set_ylabel(
            "BeamFocus RMS uncertainty from CV [%]"
        )

    ax.set_xlim(
        args.xmin,
        args.xmax,
    )

    if args.ymax is not None:
        ax.set_ylim(
            0.0,
            args.ymax,
        )
    else:
        ax.set_ylim(
            bottom=0.0
        )

    ax.grid(
        True,
        alpha=0.25,
    )

    ax.legend()

    fig.tight_layout()

    if args.output is None:
        base = os.path.splitext(
            os.path.basename(args.input)
        )[0]

        args.output = "{}_{}_{}.png".format(
            base,
            args.species,
            args.definition,
        )

    fig.savefig(
        args.output,
        dpi=160,
    )

    plt.close(fig)

    print()
    print("Input      :", args.input)
    print("Species    :", args.species)
    print("Definition :", args.definition)
    print("Parents    :", ", ".join(args.parents))
    print("Saved      :", args.output)

    if args.print_table:

        print()
        print(
            "{:>11s} {:>16s} {:>12s} {:>12s}".format(
                "E [GeV]",
                "parent",
                "N",
                "sigma [%]",
            )
        )

        print("-" * 58)

        for parent in args.parents:

            sub = plotted_rows[parent]

            for row in sub:
                print(
                    "{:5.1f}-{:5.1f} {:>16s} {:12d} {:11.4f}%".format(
                        row["E_low"],
                        row["E_high"],
                        parent,
                        row["N"],
                        row["sigma_percent"],
                    )
                )

            print()


if __name__ == "__main__":
    main()
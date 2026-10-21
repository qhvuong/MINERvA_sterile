#!/usr/bin/env python3

import argparse
import numpy as np
from matplotlib import pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(
        description="Plot LE right-sign muon BeamFocus uncertainty vs true E_nu"
    )
    p.add_argument("--input", required=True, help="Input NPZ file")
    p.add_argument(
        "--mode",
        choices=["FHC", "RHC"],
        required=True,
        help="Beam mode: FHC uses numu, RHC uses numubar",
    )
    p.add_argument(
        "--definition",
        choices=["centered", "cv"],
        default="centered",
        help="Uncertainty definition",
    )
    p.add_argument("--xmax", type=float, default=20.0)
    p.add_argument("--ymax", type=float, default=None)
    p.add_argument(
        "--output",
        required=True,
        help="Output PNG file",
    )
    return p.parse_args()


def main():
    args = parse_args()

    d = np.load(args.input)

    edges = d["le_bin_edges"]
    centers = 0.5 * (edges[:-1] + edges[1:])

    sigma = d[f"le_sigma_{args.definition}"]

    # species order:
    # 0 numu
    # 1 nue
    # 2 numubar
    # 3 nuebar
    if args.mode == "FHC":
        idx = 0
        label = r"LE FHC true $\nu_\mu$"
    else:
        idx = 2
        label = r"LE RHC true $\bar{\nu}_\mu$"

    y = 100.0 * sigma[idx]

    print()
    print("Input      :", args.input)
    print("Mode       :", args.mode)
    print("Definition :", args.definition)
    print()
    print("{:>15s} {:>12s}".format("E_true [GeV]", "uncertainty"))
    for i in range(len(centers)):
        print(
            "{:6.2f}-{:6.2f} {:11.4f}%".format(
                edges[i], edges[i + 1], y[i]
            )
        )

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.plot(
        centers,
        y,
        marker="o",
        label=label,
    )

    ax.set_xlabel(r"True $E_\nu$ [GeV]")
    ax.set_ylabel("BeamFocus fractional uncertainty [%]")

    if args.definition == "centered":
        title_def = "centered uncertainty"
    else:
        title_def = "RMS from CV"

    ax.set_title(f"{label} BeamFocus uncertainty vs true neutrino energy\n{title_def}")
    ax.grid(True, alpha=0.25)
    ax.legend()

    ax.set_xlim(edges[0], args.xmax)
    if args.ymax is not None:
        ax.set_ylim(0.0, args.ymax)

    fig.tight_layout()
    fig.savefig(args.output, dpi=160)
    plt.close(fig)

    print()
    print("Saved:", args.output)


if __name__ == "__main__":
    main()
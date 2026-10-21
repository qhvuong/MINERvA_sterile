#!/usr/bin/env python3

import argparse
import ROOT
import PlotUtils

from tools import Utilities
from tools import SystematicsUniverse


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--playlist",
        required=True,
        help="e.g. le1_p8, le7_p8, le9_p8, le13C_p8",
    )

    parser.add_argument(
        "--ntuple_tag",
        default="MAD",
        help="Ntuple tag, default: MAD",
    )

    args = parser.parse_args()

    print("[MHRW] analysis playlist =", args.playlist)
    print("[MHRW] ntuple tag =", args.ntuple_tag)

    # Build FULL Truth chain
    chain = Utilities.fileChain(
        args.playlist,
        "mc",
        args.ntuple_tag,
        "Truth",
        None,
        None,
    )

    print("[MHRW] tree =", chain.GetTree().GetName())
    print("[MHRW] entries =", chain.GetEntries())

    universes = SystematicsUniverse.GetAllSystematicsUniverses(
        chain,
        False,
        exclude=["all"],
    )

    cv = universes["cv"][0]

    print(
        "[MHRW GENERATE]",
        "tree =", chain.GetTree().GetName(),
        "MAT playlist =", cv.GetPlaylist(),
        "IsTruth =", cv.IsTruth(),
    )

    if chain.GetTree().GetName() != "Truth":
        raise RuntimeError(
            "MHRW generation requires a Truth tree, but got {}".format(
                chain.GetTree().GetName()
            )
        )

    if not cv.IsTruth():
        raise RuntimeError(
            "Universe is not in Truth mode. Refusing to generate MHRW file."
        )

    cv.SetupMHRWeighter()

    print("[MHRW GENERATE] finished")


if __name__ == "__main__":
    main()
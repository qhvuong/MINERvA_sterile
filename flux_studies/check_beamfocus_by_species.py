#!/usr/bin/env python

import argparse
import array
import math
import ROOT


PDG_NAMES = {
    14:  "nu_mu",
    -14: "nubar_mu",
    12:  "nu_e",
    -12: "nubar_e",
}


def load_chain(filelist, tree_name):
    chain = ROOT.TChain(tree_name)

    with open(filelist) as f:
        for line in f:
            filename = line.strip()

            if not filename or filename.startswith("#"):
                continue

            chain.Add(filename)

    return chain


def main():
    parser = argparse.ArgumentParser(
        description="Check P8 LE BeamFocus uncertainties by neutrino species."
    )

    parser.add_argument(
        "--filelist",
        required=True,
        help="Text file containing the input ROOT files."
    )

    parser.add_argument(
        "--tree",
        default="Truth",
        help="Input tree name. Default: Truth"
    )

    parser.add_argument(
        "--n_universes",
        type=int,
        default=1000,
        help="Number of BeamFocus universes. Default: 1000"
    )

    parser.add_argument(
        "--max_events",
        type=int,
        default=10000,
        help="Maximum number of entries to inspect. Use -1 for all entries."
    )

    parser.add_argument(
        "--RHC",
        action="store_true",
        help="Treat the input as RHC. Default is FHC."
    )

    parser.add_argument(
        "--tolerance",
        type=float,
        default=1.0e-6,
        help="Universe range considered nontrivial. Default: 1e-6"
    )

    args = parser.parse_args()

    chain = load_chain(args.filelist, args.tree)

    required_branches = [
        "mc_incoming",
        "mc_hornCurrent_cvweight",
        "mc_wgt_Flux_BeamFocus",
    ]

    for branch_name in required_branches:
        if not chain.GetBranch(branch_name):
            raise RuntimeError(
                "Required branch '{}' was not found.".format(branch_name)
            )

    chain.SetBranchStatus("*", 0)

    for branch_name in required_branches:
        chain.SetBranchStatus(branch_name, 1)

    incoming = array.array("i", [0])
    beam_cv = array.array("d", [0.0])
    beam_weights_int = array.array("i", [0] * args.n_universes)

    chain.SetBranchAddress("mc_incoming", incoming)
    chain.SetBranchAddress("mc_hornCurrent_cvweight", beam_cv)
    chain.SetBranchAddress("mc_wgt_Flux_BeamFocus", beam_weights_int)

    stats = {}

    n_entries = chain.GetEntries()

    if args.max_events >= 0:
        n_entries = min(n_entries, args.max_events)

    for ientry in range(n_entries):
        chain.GetEntry(ientry)

        pdg = incoming[0]
        cv = beam_cv[0]

        if pdg not in stats:
            stats[pdg] = {
                "events": 0,
                "skipped": 0,
                "events_nontrivial": 0,
                "sum_event_rms": 0.0,
                "max_universe_range": 0.0,
                "first_relative_weights": None,
            }

        species = stats[pdg]
        species["events"] += 1

        if math.isnan(cv) or abs(cv) < 1.0e-12:
            species["skipped"] += 1
            continue

        relative_weights = [
            beam_weights_int[u] * 1.0e-7 / cv
            for u in range(args.n_universes)
        ]

        event_mean = sum(relative_weights) / args.n_universes

        event_rms = math.sqrt(
            sum((weight - event_mean) ** 2 for weight in relative_weights)
            / args.n_universes
        )

        event_range = max(relative_weights) - min(relative_weights)

        if event_range > args.tolerance:
            species["events_nontrivial"] += 1

        species["sum_event_rms"] += event_rms
        species["max_universe_range"] = max(
            species["max_universe_range"],
            event_range
        )

        if species["first_relative_weights"] is None:
            species["first_relative_weights"] = relative_weights[:10]

    right_sign_pdg = -14 if args.RHC else 14
    mode_name = "RHC" if args.RHC else "FHC"

    print()
    print("BeamFocus uncertainty check")
    print("===========================")
    print("Mode             :", mode_name)
    print("Right-sign PDG   :", right_sign_pdg)
    print("Entries inspected:", n_entries)
    print("Universes        :", args.n_universes)
    print("Tolerance        :", args.tolerance)
    print()

    ordered_pdgs = [14, -14, 12, -12]

    for pdg in ordered_pdgs:
        if pdg not in stats:
            continue

        species = stats[pdg]
        valid_events = species["events"] - species["skipped"]

        if valid_events > 0:
            mean_event_rms = (
                species["sum_event_rms"] / valid_events
            )

            nontrivial_fraction = (
                float(species["events_nontrivial"])
                / valid_events
            )
        else:
            mean_event_rms = 0.0
            nontrivial_fraction = 0.0

        sign_label = (
            "RIGHT SIGN"
            if pdg == right_sign_pdg
            else "NON-RIGHT SIGN"
        )

        print("{} ({:+d}) [{}]".format(
            PDG_NAMES.get(pdg, "unknown"),
            pdg,
            sign_label
        ))

        print("  Events                    :", species["events"])
        print("  Invalid CV events         :", species["skipped"])
        print("  Events with nonzero spread:", species["events_nontrivial"])
        print("  Nontrivial event fraction : {:.6f}".format(
            nontrivial_fraction
        ))
        print("  Mean RMS across universes : {:.8e}".format(
            mean_event_rms
        ))
        print("  Maximum universe range    : {:.8e}".format(
            species["max_universe_range"]
        ))

        if species["first_relative_weights"] is not None:
            print(
                "  First event, universes 0-9:",
                " ".join(
                    "{:.7f}".format(weight)
                    for weight in species["first_relative_weights"]
                )
            )

        print()

    print("Interpretation")
    print("--------------")
    print(
        "If the colleague's description is correct, the right-sign muon "
        "species should have a visible universe-to-universe RMS and range."
    )
    print(
        "The non-right-sign muon and electron species should have little or "
        "no universe-to-universe spread, even if their common value is not 1."
    )


if __name__ == "__main__":
    main()
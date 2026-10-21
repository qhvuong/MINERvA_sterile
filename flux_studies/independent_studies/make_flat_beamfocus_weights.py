#!/usr/bin/env python

import argparse
import array
import csv
import math
import ROOT
import numpy as np

PDG_NAMES = {
    14: "nu_mu",
    -14: "nubar_mu",
    12: "nu_e",
    -12: "nubar_e",
}

ORDERED_PDGS = [14, -14, 12, -12]


def add_filelist_to_chain(chain, filelist):
    nfiles = 0

    with open(filelist) as f:
        for line in f:
            filename = line.strip()

            if not filename or filename.startswith("#"):
                continue

            chain.Add(filename)
            nfiles += 1

    return nfiles


def passes_compute_flux_fv(vtx):
    x = vtx[0]
    y = vtx[1]
    z = vtx[2]

    if z <= 6117.0 or z >= 8193.0:
        return False

    ax = abs(x)
    ay = abs(y)

    if ax >= 850.0:
        return False

    sqrt3 = math.sqrt(3.0)

    inside_hexagon = (
        ay < 850.0 / sqrt3
        or
        ay < 850.0 * 2.0 / sqrt3 - ax / sqrt3
    )

    return inside_hexagon


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build flat species-dependent BeamFocus universe weights "
            "from P8 ME tuples."
        )
    )

    parser.add_argument(
        "--filelist",
        action="append",
        required=True,
        help=(
            "ME tuple file list. Repeat this option to combine multiple "
            "playlists."
        )
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
        default=-1,
        help="Maximum total entries to process. Default: all entries"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output CSV file"
    )

    parser.add_argument(
        "--mode",
        choices=["FHC", "RHC"],
        required=True,
        help="Horn-current mode"
    )

    parser.add_argument(
        "--accumulator_output",
        help="Save raw accumulator data to an NPZ file."
    )

    args = parser.parse_args()

    chain = ROOT.TChain(args.tree)

    total_files = 0

    for filelist in args.filelist:
        nfiles = add_filelist_to_chain(chain, filelist)
        total_files += nfiles
        print("Added {} files from {}".format(nfiles, filelist))

    print("Total files:", total_files)
    print("Total entries:", chain.GetEntries())

    required_branches = [
        "mc_incoming",
        "mc_current",
        "mc_intType",
        "mc_targetZ",
        "mc_vtx",
        "mc_cvweight_total",
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
    current = array.array("i", [0])
    int_type = array.array("i", [0])
    target_z = array.array("i", [0])
    vtx = array.array("d", [0.0, 0.0, 0.0])

    event_cv = array.array("d", [0.0])
    beam_cv = array.array("d", [0.0])
    beam_weights_int = array.array("i", [0] * args.n_universes)

    chain.SetBranchAddress("mc_incoming", incoming)
    chain.SetBranchAddress("mc_current", current)
    chain.SetBranchAddress("mc_intType", int_type)
    chain.SetBranchAddress("mc_targetZ", target_z)
    chain.SetBranchAddress("mc_vtx", vtx)

    chain.SetBranchAddress("mc_cvweight_total", event_cv)
    chain.SetBranchAddress("mc_hornCurrent_cvweight", beam_cv)
    chain.SetBranchAddress(
        "mc_wgt_Flux_BeamFocus",
        beam_weights_int
    )

    cv_sums = {
        pdg: 0.0
        for pdg in ORDERED_PDGS
    }

    universe_sums = {
        pdg: [0.0] * args.n_universes
        for pdg in ORDERED_PDGS
    }

    event_counts = {
        pdg: 0
        for pdg in ORDERED_PDGS
    }

    skipped_counts = {
        pdg: 0
        for pdg in ORDERED_PDGS
    }

    n_entries = chain.GetEntries()

    if args.max_events >= 0:
        n_entries = min(n_entries, args.max_events)

    for ientry in range(n_entries):
        if ientry % 100000 == 0:
            print("Processing entry {} / {}".format(
                ientry,
                n_entries
            ))

        chain.GetEntry(ientry)

        pdg = incoming[0]

        if pdg not in PDG_NAMES:
            continue

        # Match compute_flux truth selection
        if current[0] != 1:
            continue

        if int_type[0] == 8:
            continue

        if target_z[0] != 6:
            continue

        if not passes_compute_flux_fv(vtx):
            continue

        event_counts[pdg] += 1

        cv = event_cv[0]
        beam_cv_value = beam_cv[0]

        if (
            math.isnan(cv) or
            math.isnan(beam_cv_value) or
            cv < 1.0e-12 or
            abs(beam_cv_value) < 1.0e-12
        ):
            skipped_counts[pdg] += 1
            continue

        cv_sums[pdg] += cv

        for u in range(args.n_universes):
            beam_absolute = beam_weights_int[u] * 1.0e-7
            beam_relative = beam_absolute / beam_cv_value

            universe_sums[pdg][u] += (
                cv * beam_relative
            )

    if args.accumulator_output:
        np.savez(
            args.accumulator_output,
            n_universes=args.n_universes,
            mode=args.mode,
            pdgs=np.array(ORDERED_PDGS, dtype=int),
            cv_sums=np.array(
                [cv_sums[pdg] for pdg in ORDERED_PDGS],
                dtype=float
            ),
            universe_sums=np.array(
                [universe_sums[pdg] for pdg in ORDERED_PDGS],
                dtype=float
            ),
            event_counts=np.array(
                [event_counts[pdg] for pdg in ORDERED_PDGS],
                dtype=np.int64
            ),
            skipped_counts=np.array(
                [skipped_counts[pdg] for pdg in ORDERED_PDGS],
                dtype=np.int64
            ),
        )

        print("Saved accumulators to:", args.accumulator_output)

    flat_weights = {
        pdg: [1.0] * args.n_universes
        for pdg in ORDERED_PDGS
    }

    for pdg in ORDERED_PDGS:
        if cv_sums[pdg] <= 0.0:
            print(
                "WARNING: no valid CV contribution for {}".format(
                    PDG_NAMES[pdg]
                )
            )
            continue

        for u in range(args.n_universes):
            flat_weights[pdg][u] = (
                universe_sums[pdg][u] /
                cv_sums[pdg]
            )

    with open(args.output, "w") as csvfile:
        writer = csv.writer(csvfile)

        writer.writerow([
            "universe",
            "nu_mu",
            "nubar_mu",
            "nu_e",
            "nubar_e",
        ])

        for u in range(args.n_universes):
            writer.writerow([
                u,
                "{:.10f}".format(flat_weights[14][u]),
                "{:.10f}".format(flat_weights[-14][u]),
                "{:.10f}".format(flat_weights[12][u]),
                "{:.10f}".format(flat_weights[-12][u]),
            ])

    print()
    print("Summary")
    print("=======")
    print("Mode:", args.mode)
    print("Entries processed:", n_entries)
    print("Output:", args.output)
    print()

    for pdg in ORDERED_PDGS:
        weights = flat_weights[pdg]

        mean_weight = (
            sum(weights) / args.n_universes
        )

        rms = math.sqrt(
            sum(
                (weight - mean_weight) ** 2
                for weight in weights
            ) / args.n_universes
        )

        print("{} ({:+d})".format(
            PDG_NAMES[pdg],
            pdg
        ))
        print("  Events       :", event_counts[pdg])
        print("  Skipped      :", skipped_counts[pdg])
        print("  CV sum       : {:.8e}".format(cv_sums[pdg]))
        print("  Universe mean: {:.8f}".format(mean_weight))
        print("  Universe RMS : {:.8f}".format(rms))
        print("  Minimum      : {:.8f}".format(min(weights)))
        print("  Maximum      : {:.8f}".format(max(weights)))
        print()

    print("First five CSV rows:")

    for u in range(min(5, args.n_universes)):
        print(
            u,
            "{:.7f}".format(flat_weights[14][u]),
            "{:.7f}".format(flat_weights[-14][u]),
            "{:.7f}".format(flat_weights[12][u]),
            "{:.7f}".format(flat_weights[-12][u])
        )


if __name__ == "__main__":
    main()
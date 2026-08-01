#!/usr/bin/env python3

import argparse
import math
import sys

import ROOT

ROOT.gROOT.SetBatch(True)


def get_leaf(tree, branch_name, required=True):
    leaf = tree.GetLeaf(branch_name)

    if leaf:
        return leaf

    branch = tree.GetBranch(branch_name)

    if branch:
        leaves = branch.GetListOfLeaves()

        if leaves and leaves.GetEntries() == 1:
            return leaves.At(0)

    if required:
        raise RuntimeError(
            f"Could not obtain leaf for branch '{branch_name}'"
        )

    return None


def read_scalar(leaf):
    value = float(leaf.GetValue())

    if not math.isfinite(value):
        raise RuntimeError(f"Non-finite value read from {leaf.GetName()}")

    return value


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare total flux CV weights with PPFX and horn-current "
            "CV weights in a MINERvA AnaTuple."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input ROOT file, including /pnfs/... paths.",
    )

    parser.add_argument(
        "--tree",
        default="MasterAnaDev",
        help="Tree name. Default: MasterAnaDev",
    )

    parser.add_argument(
        "--events",
        type=int,
        default=20,
        help="Number of entries to inspect. Default: 20",
    )

    parser.add_argument(
        "--start-entry",
        type=int,
        default=0,
        help="First entry to inspect. Default: 0",
    )

    args = parser.parse_args()

    root_file = ROOT.TFile.Open(args.input, "READ")

    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open {args.input}")

    tree = root_file.Get(args.tree)

    if not tree:
        raise RuntimeError(
            f"Could not find tree '{args.tree}' in {args.input}"
        )

    total_leaf = get_leaf(tree, "mc_cvweight_totalFlux")
    gen1_total_leaf = get_leaf(
        tree,
        "mc_gen1_cvweight_totalFlux",
        required=False,
    )
    ppfx_leaf = get_leaf(tree, "mc_ppfx1_cvweight")
    horn_leaf = get_leaf(tree, "mc_hornCurrent_cvweight")

    n_entries = int(tree.GetEntries())
    start = max(0, args.start_entry)
    stop = min(n_entries, start + args.events)

    print(f"Input:   {args.input}")
    print(f"Tree:    {args.tree}")
    print(f"Entries: {n_entries}")
    print(f"Range:   [{start}, {stop})")

    print()
    print(
        " entry"
        "          totalFluxCV"
        "               ppfxCV"
        "               hornCV"
        "        ppfx*horn"
        "   total/(ppfx*horn)"
        "       total/ppfx"
        "       total-(ppfx*horn)"
    )

    sum_abs_diff = 0.0
    max_abs_diff = 0.0
    sum_abs_ratio_minus_one = 0.0
    max_abs_ratio_minus_one = 0.0
    valid = 0

    gen1_sum_abs_diff = 0.0
    gen1_max_abs_diff = 0.0
    gen1_valid = 0

    for entry in range(start, stop):
        if tree.GetEntry(entry) <= 0:
            print(f"WARNING: failed to read entry {entry}")
            continue

        total = read_scalar(total_leaf)
        ppfx = read_scalar(ppfx_leaf)
        horn = read_scalar(horn_leaf)

        product = ppfx * horn

        ratio_product = (
            total / product
            if product != 0.0
            else float("nan")
        )

        ratio_ppfx = (
            total / ppfx
            if ppfx != 0.0
            else float("nan")
        )

        difference = total - product

        print(
            f"{entry:6d}"
            f"{total:21.12g}"
            f"{ppfx:21.12g}"
            f"{horn:21.12g}"
            f"{product:21.12g}"
            f"{ratio_product:21.12g}"
            f"{ratio_ppfx:18.12g}"
            f"{difference:24.12g}"
        )

        if math.isfinite(ratio_product):
            abs_diff = abs(difference)
            abs_ratio_minus_one = abs(ratio_product - 1.0)

            sum_abs_diff += abs_diff
            max_abs_diff = max(max_abs_diff, abs_diff)

            sum_abs_ratio_minus_one += abs_ratio_minus_one
            max_abs_ratio_minus_one = max(
                max_abs_ratio_minus_one,
                abs_ratio_minus_one,
            )

            valid += 1

        if gen1_total_leaf:
            gen1_total = read_scalar(gen1_total_leaf)
            gen1_difference = gen1_total - total

            gen1_sum_abs_diff += abs(gen1_difference)
            gen1_max_abs_diff = max(
                gen1_max_abs_diff,
                abs(gen1_difference),
            )

            gen1_valid += 1

            print(
                f"       mc_gen1_cvweight_totalFlux={gen1_total:.12g}"
                f"  gen1-total={gen1_difference:+.12g}"
            )

    print()
    print("=" * 96)
    print("Summary")

    if valid:
        print(
            "Mean |totalFluxCV - ppfxCV*hornCV|: "
            f"{sum_abs_diff / valid:.12g}"
        )
        print(
            "Max  |totalFluxCV - ppfxCV*hornCV|: "
            f"{max_abs_diff:.12g}"
        )
        print(
            "Mean |totalFluxCV/(ppfxCV*hornCV) - 1|: "
            f"{sum_abs_ratio_minus_one / valid:.12g}"
        )
        print(
            "Max  |totalFluxCV/(ppfxCV*hornCV) - 1|: "
            f"{max_abs_ratio_minus_one:.12g}"
        )

    if gen1_total_leaf and gen1_valid:
        print()
        print(
            "Mean |mc_gen1_cvweight_totalFlux "
            "- mc_cvweight_totalFlux|: "
            f"{gen1_sum_abs_diff / gen1_valid:.12g}"
        )
        print(
            "Max  |mc_gen1_cvweight_totalFlux "
            "- mc_cvweight_totalFlux|: "
            f"{gen1_max_abs_diff:.12g}"
        )

    root_file.Close()


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
#!/usr/bin/env python3

import argparse
import glob
import math
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)


REQUIRED_BRANCHES = (
    "mc_wgt_ppfx1_Total",
    "mc_wgt_Flux_BeamFocus",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Inspect PPFX, BeamFocus, and CV weights in a tuple and compare "
            "additive versus multiplicative Flux-weight constructions."
        )
    )

    parser.add_argument(
        "--input",
        nargs="+",
        required=True,
        help="Input ROOT file(s) or glob patterns.",
    )

    parser.add_argument(
        "--tree",
        default=None,
        help=(
            "Tree path/name. If omitted, the script searches for a TTree "
            "containing the required branches."
        ),
    )

    parser.add_argument(
        "--events",
        type=int,
        default=10,
        help="Number of accepted events to print. Default: 10",
    )

    parser.add_argument(
        "--universes",
        type=int,
        default=5,
        help="Number of universes per event to inspect. Default: 5",
    )

    parser.add_argument(
        "--start-entry",
        type=int,
        default=0,
        help="First tree entry to inspect. Default: 0",
    )

    parser.add_argument(
        "--integer-scale",
        type=float,
        default=1.0e-7,
        help=(
            "Scale applied to stored integer universe weights. "
            "Default: 1e-7"
        ),
    )

    parser.add_argument(
        "--total-branch",
        default=None,
        help=(
            "Optional stored total-Flux universe branch. If supplied, the "
            "script tests whether it matches the additive or multiplicative "
            "construction."
        ),
    )

    parser.add_argument(
        "--total-cv-branch",
        default=None,
        help=(
            "Optional CV denominator associated with --total-branch. "
            "If omitted, the total universe values are treated as relative "
            "weights."
        ),
    )

    parser.add_argument(
        "--offset",
        type=int,
        default=0,
        help=(
            "PPFX universe offset. For example, --offset 37 pairs "
            "BeamFocus universe i with PPFX universe i+37. Default: 0"
        ),
    )

    return parser.parse_args()


def expand_input_files(patterns):
    files = []

    for pattern in patterns:
        matches = sorted(glob.glob(pattern))

        if matches:
            files.extend(matches)
        elif os.path.isfile(pattern):
            files.append(pattern)
        else:
            print(
                f"WARNING: no file matched '{pattern}'",
                file=sys.stderr,
            )

    # Preserve order while removing duplicates.
    return list(dict.fromkeys(files))


def iter_root_objects(directory, prefix=""):
    """Recursively yield (path, object) from a ROOT directory."""
    keys = directory.GetListOfKeys()

    if not keys:
        return

    for key in keys:
        name = key.GetName()
        path = f"{prefix}/{name}" if prefix else name
        obj = key.ReadObj()

        if obj.InheritsFrom("TDirectory"):
            yield from iter_root_objects(obj, path)
        else:
            yield path, obj


def find_tree_path(filename):
    root_file = ROOT.TFile.Open(filename, "READ")

    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open {filename}")

    candidates = []

    for path, obj in iter_root_objects(root_file):
        if not obj.InheritsFrom("TTree"):
            continue

        branches = obj.GetListOfBranches()
        branch_names = {
            branches.At(i).GetName()
            for i in range(branches.GetEntries())
        }

        if all(name in branch_names for name in REQUIRED_BRANCHES):
            candidates.append((path, obj.GetEntries()))

    root_file.Close()

    if not candidates:
        raise RuntimeError(
            "Could not find a TTree containing both "
            f"{REQUIRED_BRANCHES[0]} and {REQUIRED_BRANCHES[1]}"
        )

    if len(candidates) > 1:
        print("Found multiple candidate trees:")

        for path, entries in candidates:
            print(f"  {path}: {entries} entries")

        print(f"Using: {candidates[0][0]}")

    return candidates[0][0]


def has_branch(tree, branch_name):
    return bool(branch_name) and bool(tree.GetBranch(branch_name))


def get_leaf(tree, branch_name, required=True):
    if not branch_name:
        return None

    leaf = tree.GetLeaf(branch_name)

    if leaf:
        return leaf

    # Some branches have a leaf name different from the branch name.
    branch = tree.GetBranch(branch_name)

    if branch:
        leaves = branch.GetListOfLeaves()

        if leaves and leaves.GetEntries() == 1:
            return leaves.At(0)

    if required:
        raise RuntimeError(
            f"Could not obtain a leaf for branch '{branch_name}'"
        )

    return None


def scalar_value(leaf, default=1.0):
    if leaf is None:
        return default

    value = float(leaf.GetValue())

    if not math.isfinite(value) or value == 0.0:
        return default

    return value


def universe_value(leaf, index, integer_scale):
    value = float(leaf.GetValue(index))

    # The MINERvA tuple arrays under consideration are stored as integers
    # scaled by 1e7. Adjust --integer-scale if inspecting another format.
    return value * integer_scale


def available_universes(leaf):
    """
    Get the number of valid elements for the current entry.

    GetNdata() is useful for variable-length arrays. GetLen() is retained
    as a fallback.
    """
    ndata = int(leaf.GetNdata())

    if ndata > 0:
        return ndata

    return int(leaf.GetLen())


def print_available_flux_branches(tree):
    print("\nBranches whose names contain 'flux', 'ppfx', or 'horn':")

    branches = tree.GetListOfBranches()

    matched = []

    for i in range(branches.GetEntries()):
        name = branches.At(i).GetName()
        lowered = name.lower()

        if (
            "flux" in lowered
            or "ppfx" in lowered
            or "horn" in lowered
        ):
            matched.append(name)

    for name in sorted(matched):
        print(f"  {name}")


def main():
    args = parse_args()
    files = expand_input_files(args.input)

    if not files:
        raise RuntimeError("No valid input ROOT files were found.")

    tree_path = args.tree or find_tree_path(files[0])

    print(f"Tree: {tree_path}")
    print(f"Input files: {len(files)}")

    chain = ROOT.TChain(tree_path)

    for filename in files:
        added = chain.Add(filename)

        if added == 0:
            print(
                f"WARNING: failed to add {filename}",
                file=sys.stderr,
            )

    n_entries = int(chain.GetEntries())

    if n_entries <= 0:
        raise RuntimeError("The TChain contains no entries.")

    print(f"Entries: {n_entries}")

    print_available_flux_branches(chain)

    ppfx_leaf = get_leaf(chain, "mc_wgt_ppfx1_Total")
    focus_leaf = get_leaf(chain, "mc_wgt_Flux_BeamFocus")

    ppfx_cv_leaf = get_leaf(
        chain,
        "mc_ppfx1_cvweight",
        required=False,
    )

    focus_cv_leaf = get_leaf(
        chain,
        "mc_hornCurrent_cvweight",
        required=False,
    )

    total_leaf = None
    total_cv_leaf = None

    if args.total_branch:
        total_leaf = get_leaf(
            chain,
            args.total_branch,
            required=True,
        )

        if args.total_cv_branch:
            total_cv_leaf = get_leaf(
                chain,
                args.total_cv_branch,
                required=True,
            )

    print("\nCV branches:")
    print(
        "  PPFX CV:  "
        + (
            "mc_ppfx1_cvweight"
            if ppfx_cv_leaf
            else "not found; using 1"
        )
    )
    print(
        "  Focus CV: "
        + (
            "mc_hornCurrent_cvweight"
            if focus_cv_leaf
            else "not found; using 1"
        )
    )

    if total_leaf:
        print(f"  Total universe branch: {args.total_branch}")
        print(
            "  Total CV branch: "
            + (
                args.total_cv_branch
                if total_cv_leaf
                else "not supplied; treating total as relative"
            )
        )
    else:
        print(
            "\nNo stored total-Flux branch was requested. The script will "
            "compare the two candidate formulas, but cannot determine which "
            "one matches an independently stored total universe."
        )

    start = max(0, args.start_entry)
    stop = min(n_entries, start + args.events)

    sum_abs_mult_minus_add = 0.0
    max_abs_mult_minus_add = 0.0
    n_formula_values = 0

    sum_sq_total_minus_mult = 0.0
    sum_sq_total_minus_add = 0.0
    max_total_minus_mult = 0.0
    max_total_minus_add = 0.0
    n_total_values = 0

    printed_events = 0

    for entry in range(start, stop):
        bytes_read = chain.GetEntry(entry)

        if bytes_read <= 0:
            print(f"WARNING: could not read entry {entry}")
            continue

        ppfx_cv = scalar_value(ppfx_cv_leaf, default=1.0)
        focus_cv = scalar_value(focus_cv_leaf, default=1.0)

        n_ppfx = available_universes(ppfx_leaf)
        n_focus = available_universes(focus_leaf)

        if n_ppfx <= 0 or n_focus <= 0:
            print(
                f"Entry {entry}: invalid array lengths "
                f"PPFX={n_ppfx}, Focus={n_focus}"
            )
            continue

        n_common = min(n_ppfx, n_focus)
        n_to_print = min(args.universes, n_common)

        print("\n" + "=" * 88)
        print(
            f"Entry {entry}: "
            f"ppfxCV={ppfx_cv:.9g}, "
            f"focusCV={focus_cv:.9g}, "
            f"Nppfx={n_ppfx}, Nfocus={n_focus}"
        )

        print(
            " final  ppfx  focus      ppfxRel       focusRel"
            "       multiply       additive       mult-add"
        )

        for final_index in range(n_to_print):
            focus_index = final_index
            ppfx_index = (
                final_index + args.offset
            ) % n_ppfx

            ppfx_abs = universe_value(
                ppfx_leaf,
                ppfx_index,
                args.integer_scale,
            )

            focus_abs = universe_value(
                focus_leaf,
                focus_index,
                args.integer_scale,
            )

            ppfx_relative = ppfx_abs / ppfx_cv
            focus_relative = focus_abs / focus_cv

            multiplicative = ppfx_relative * focus_relative
            additive = ppfx_relative + focus_relative - 1.0

            formula_difference = multiplicative - additive

            sum_abs_mult_minus_add += abs(formula_difference)
            max_abs_mult_minus_add = max(
                max_abs_mult_minus_add,
                abs(formula_difference),
            )
            n_formula_values += 1

            line = (
                f"{final_index:6d}"
                f"{ppfx_index:6d}"
                f"{focus_index:7d}"
                f"{ppfx_relative:15.7g}"
                f"{focus_relative:15.7g}"
                f"{multiplicative:15.7g}"
                f"{additive:15.7g}"
                f"{formula_difference:15.7g}"
            )

            if total_leaf:
                n_total = available_universes(total_leaf)

                if final_index < n_total:
                    total_abs = universe_value(
                        total_leaf,
                        final_index,
                        args.integer_scale,
                    )

                    total_cv = scalar_value(
                        total_cv_leaf,
                        default=1.0,
                    )

                    total_relative = total_abs / total_cv

                    residual_mult = (
                        total_relative - multiplicative
                    )
                    residual_add = total_relative - additive

                    sum_sq_total_minus_mult += residual_mult**2
                    sum_sq_total_minus_add += residual_add**2

                    max_total_minus_mult = max(
                        max_total_minus_mult,
                        abs(residual_mult),
                    )

                    max_total_minus_add = max(
                        max_total_minus_add,
                        abs(residual_add),
                    )

                    n_total_values += 1

                    line += (
                        f"  total={total_relative:.7g}"
                        f"  dMult={residual_mult:+.4e}"
                        f"  dAdd={residual_add:+.4e}"
                    )

            print(line)

        printed_events += 1

    print("\n" + "=" * 88)
    print("Summary")

    if n_formula_values:
        print(
            "Mean |multiplicative - additive|: "
            f"{sum_abs_mult_minus_add / n_formula_values:.8g}"
        )
        print(
            "Max  |multiplicative - additive|: "
            f"{max_abs_mult_minus_add:.8g}"
        )

    if n_total_values:
        rms_mult = math.sqrt(
            sum_sq_total_minus_mult / n_total_values
        )
        rms_add = math.sqrt(
            sum_sq_total_minus_add / n_total_values
        )

        print(
            f"\nCompared with stored branch '{args.total_branch}':"
        )
        print(
            f"RMS(total - multiplicative): {rms_mult:.8g}"
        )
        print(
            f"RMS(total - additive):       {rms_add:.8g}"
        )
        print(
            "Max |total - multiplicative|: "
            f"{max_total_minus_mult:.8g}"
        )
        print(
            "Max |total - additive|:       "
            f"{max_total_minus_add:.8g}"
        )

        if rms_mult < rms_add:
            print(
                "\nThe stored total branch is closer to the "
                "multiplicative construction."
            )
        elif rms_add < rms_mult:
            print(
                "\nThe stored total branch is closer to the "
                "additive-shift construction."
            )
        else:
            print(
                "\nThe two constructions agree equally well within "
                "the inspected sample."
            )

    if printed_events == 0:
        raise RuntimeError("No entries were successfully inspected.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
#!/usr/bin/env python3

import argparse
import math
import statistics
import ROOT


SCALE = 1.0e-7


def finite(x):
    try:
        return math.isfinite(float(x))
    except Exception:
        return False


def rms(values):
    values = [float(x) for x in values if finite(x)]

    if not values:
        return float("nan")

    return math.sqrt(sum(x * x for x in values) / len(values))


def median(values):
    values = [float(x) for x in values if finite(x)]

    if not values:
        return float("nan")

    return statistics.median(values)


def correlation(xs, ys):
    pairs = [
        (float(x), float(y))
        for x, y in zip(xs, ys)
        if finite(x) and finite(y)
    ]

    if len(pairs) < 2:
        return float("nan")

    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]

    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)

    covariance = sum(
        (x - mean_x) * (y - mean_y)
        for x, y in zip(xs, ys)
    )

    variance_x = sum((x - mean_x) ** 2 for x in xs)
    variance_y = sum((y - mean_y) ** 2 for y in ys)

    if variance_x <= 0.0 or variance_y <= 0.0:
        return float("nan")

    return covariance / math.sqrt(variance_x * variance_y)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Check whether mc_wgt_Flux_BeamFocus contains absolute "
            "or relative universe weights."
        )
    )

    input_group = parser.add_mutually_exclusive_group(required=True)

    input_group.add_argument(
        "--file",
        help="Single input ROOT tuple",
    )

    input_group.add_argument(
        "--filelist",
        help="Text file containing ROOT tuple paths",
    )

    parser.add_argument(
        "--tree",
        default="Truth",
        help="Tree name, default: Truth",
    )

    parser.add_argument(
        "--max-events",
        type=int,
        default=10000,
        help="Maximum events to inspect, default: 10000",
    )

    parser.add_argument(
        "--n-universes",
        type=int,
        default=1000,
        help="Number of BeamFocus universes, default: 1000",
    )

    parser.add_argument(
        "--pdg",
        type=int,
        default=None,
        help=(
            "Optional incoming-neutrino PDG selection, "
            "for example 14 or -14"
        ),
    )

    parser.add_argument(
        "--print-events",
        type=int,
        default=10,
        help="Number of individual events to print",
    )

    args = parser.parse_args()

    tree = ROOT.TChain(args.tree)

    if args.file:
        added = tree.Add(args.file)

        if added == 0:
            raise RuntimeError(
                "Could not add ROOT file: {}".format(args.file)
            )

        n_files = 1
        input_label = args.file

    else:
        n_files = 0

        with open(args.filelist) as filelist:
            for line in filelist:
                path = line.strip()

                if not path or path.startswith("#"):
                    continue

                added = tree.Add(path)

                if added == 0:
                    print("WARNING: could not add {}".format(path))
                    continue

                n_files += 1

        if n_files == 0:
            raise RuntimeError(
                "No ROOT files were added from {}".format(
                    args.filelist
                )
            )

        input_label = args.filelist

    n_entries = tree.GetEntries()

    if n_entries <= 0:
        raise RuntimeError(
            "TChain '{}' contains no entries from {}".format(
                args.tree,
                input_label,
            )
        )

    required_branches = [
        "mc_wgt_Flux_BeamFocus",
        "mc_hornCurrent_cvweight",
        "mc_incoming",
    ]

    for branch in required_branches:
        if not tree.GetBranch(branch):
            raise RuntimeError(
                "Tree does not contain branch '{}'".format(branch)
            )

    raw_means = []
    cv_values = []
    relative_means = []

    raw_minus_one = []
    raw_minus_cv = []
    relative_minus_one = []

    selected = 0
    skipped = 0
    printed = 0

    n_entries = tree.GetEntries()

    for entry in range(n_entries):
        if selected >= args.max_events:
            break

        tree.GetEntry(entry)

        pdg = int(tree.mc_incoming)

        if args.pdg is not None and pdg != args.pdg:
            continue

        cv = float(tree.mc_hornCurrent_cvweight)

        if not finite(cv) or abs(cv) < 1.0e-12:
            skipped += 1
            continue

        decoded = []

        for universe in range(args.n_universes):
            value = (
                float(tree.mc_wgt_Flux_BeamFocus[universe])
                * SCALE
            )

            if finite(value):
                decoded.append(value)

        if not decoded:
            skipped += 1
            continue

        raw_mean = sum(decoded) / len(decoded)
        relative_mean = raw_mean / cv

        raw_means.append(raw_mean)
        cv_values.append(cv)
        relative_means.append(relative_mean)

        raw_minus_one.append(raw_mean - 1.0)
        raw_minus_cv.append(raw_mean - cv)
        relative_minus_one.append(relative_mean - 1.0)

        if printed < args.print_events:
            universe_rms_raw = math.sqrt(
                sum((x - raw_mean) ** 2 for x in decoded)
                / len(decoded)
            )

            universe_rms_relative = universe_rms_raw / abs(cv)

            print(
                "event={:8d} pdg={:4d} "
                "beamCV={:12.8f} "
                "rawMean={:12.8f} "
                "rawMean/CV={:12.8f} "
                "rawRMS={:12.8f} "
                "relativeRMS={:12.8f}".format(
                    entry,
                    pdg,
                    cv,
                    raw_mean,
                    relative_mean,
                    universe_rms_raw,
                    universe_rms_relative,
                )
            )

            printed += 1

        selected += 1

    if selected == 0:
        raise RuntimeError("No usable events were selected")

    corr_raw_cv = correlation(raw_means, cv_values)

    print("\n[SUMMARY]")
    print("  input                   =", input_label)
    print("  files loaded            =", n_files)
    print("  selected events         =", selected)
    print("  skipped events          =", skipped)
    print("  PDG selection           =", args.pdg)
    print("  universes per event     =", args.n_universes)

    print("\n[CENTRAL VALUES]")
    print(
        "  median BeamFocus CV     = {:.8f}".format(
            median(cv_values)
        )
    )
    print(
        "  median decoded mean     = {:.8f}".format(
            median(raw_means)
        )
    )
    print(
        "  median decoded/CV       = {:.8f}".format(
            median(relative_means)
        )
    )

    print("\n[ABSOLUTE-vs-RELATIVE TESTS]")
    print(
        "  RMS(decoded mean - CV)  = {:.8e}".format(
            rms(raw_minus_cv)
        )
    )
    print(
        "  RMS(decoded mean - 1)   = {:.8e}".format(
            rms(raw_minus_one)
        )
    )
    print(
        "  RMS(decoded/CV - 1)     = {:.8e}".format(
            rms(relative_minus_one)
        )
    )
    print(
        "  corr(decoded mean, CV)  = {:.8f}".format(
            corr_raw_cv
        )
    )

    absolute_score = rms(raw_minus_cv)
    relative_score = rms(raw_minus_one)

    print("\n[INTERPRETATION]")

    if (
        finite(absolute_score)
        and finite(relative_score)
        and absolute_score < 0.5 * relative_score
    ):
        print(
            "  Decoded BeamFocus universes look ABSOLUTE."
        )
        print(
            "  Use: decoded universe / mc_hornCurrent_cvweight"
        )

    elif (
        finite(absolute_score)
        and finite(relative_score)
        and relative_score < 0.5 * absolute_score
    ):
        print(
            "  Decoded BeamFocus universes look RELATIVE."
        )
        print(
            "  Use the decoded universe directly."
        )

    else:
        print(
            "  Result is ambiguous from the simple RMS test."
        )
        print(
            "  Inspect individual events and the correlation."
        )
        print(
            "  Strong decoded-mean/CV agreement and high correlation "
            "support ABSOLUTE weights."
        )



if __name__ == "__main__":
    main()
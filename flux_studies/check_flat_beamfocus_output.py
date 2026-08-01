#!/usr/bin/env python

import argparse
import csv
import math
import sys

import ROOT


def load_csv_weights(csv_path, column):
    rows = {}

    with open(csv_path, "r") as infile:
        reader = csv.DictReader(infile)

        required = {"universe", column}
        missing = required.difference(reader.fieldnames or [])

        if missing:
            raise RuntimeError(
                "CSV is missing required columns: {}".format(
                    ", ".join(sorted(missing))
                )
            )

        for row in reader:
            universe = int(row["universe"])
            rows[universe] = float(row[column])

    if not rows:
        raise RuntimeError("No weights were loaded from the CSV.")

    expected_indices = set(range(len(rows)))

    if set(rows) != expected_indices:
        raise RuntimeError(
            "CSV universe indices are not contiguous from 0 to {}.".format(
                len(rows) - 1
            )
        )

    return [rows[u] for u in range(len(rows))]


def find_mnvh1d(root_file, requested_name=None):
    if requested_name:
        obj = root_file.Get(requested_name)

        if not obj:
            raise RuntimeError(
                "Histogram '{}' was not found.".format(requested_name)
            )

        return requested_name, obj

    candidates = []

    for key in root_file.GetListOfKeys():
        obj = key.ReadObj()

        if obj.InheritsFrom("PlotUtils::MnvH1D"):
            candidates.append((key.GetName(), obj))

    if not candidates:
        raise RuntimeError("No PlotUtils::MnvH1D objects were found.")

    print("Found MnvH1D objects:")

    for index, (name, _) in enumerate(candidates):
        print("  [{}] {}".format(index, name))

    if len(candidates) > 1:
        print(
            "\nNo histogram name was supplied; checking the first object.",
            file=sys.stderr,
        )

    return candidates[0]


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Check that a flat BeamFocus error band matches the selected "
            "CSV column after applying the universe offset."
        )
    )

    parser.add_argument("--root-file", required=True)
    parser.add_argument("--csv", required=True)
    parser.add_argument("--column", default="nu_e")
    parser.add_argument("--offset", type=int, default=509)
    parser.add_argument("--hist-name", default=None)
    parser.add_argument("--band-name", default="Flux_BeamFocus")
    parser.add_argument("--tolerance", type=float, default=1.0e-6)
    parser.add_argument("--print-universes", type=int, default=10)

    args = parser.parse_args()

    weights = load_csv_weights(args.csv, args.column)
    n_csv = len(weights)

    root_file = ROOT.TFile.Open(args.root_file, "READ")

    if not root_file or root_file.IsZombie():
        raise RuntimeError(
            "Could not open ROOT file '{}'.".format(args.root_file)
        )

    hist_name, mnv_hist = find_mnvh1d(
        root_file,
        requested_name=args.hist_name,
    )

    band_names = [
        str(name) for name in mnv_hist.GetVertErrorBandNames()
    ]

    if args.band_name not in band_names:
        raise RuntimeError(
            "Histogram '{}' does not contain band '{}'. Available bands: {}"
            .format(
                hist_name,
                args.band_name,
                ", ".join(band_names),
            )
        )

    band = mnv_hist.GetVertErrorBand(args.band_name)
    n_universes = int(band.GetNHists())

    if n_universes != n_csv:
        raise RuntimeError(
            "ROOT band has {} universes, while CSV has {}.".format(
                n_universes,
                n_csv,
            )
        )

    offset = args.offset % n_universes

    # The MnvH1D central-value histogram.
    cv_hist = mnv_hist.GetCVHistoWithStatError()

    global_max_abs_diff = 0.0
    failed_universes = 0
    checked_bins = 0

    print("\nChecking histogram: {}".format(hist_name))
    print("Band:               {}".format(args.band_name))
    print("CSV column:         {}".format(args.column))
    print("Offset:             {}".format(offset))
    print("Universes:          {}".format(n_universes))
    print("Tolerance:          {:.3e}\n".format(args.tolerance))

    for output_u in range(n_universes):
        source_u = (output_u + offset) % n_universes
        expected = weights[source_u]

        universe_hist = band.GetHist(output_u)

        universe_max_diff = 0.0
        universe_bins = 0

        # Include regular bins only. Underflow and overflow can be added
        # separately if desired.
        for bin_idx in range(1, cv_hist.GetNbinsX() + 1):
            cv_content = cv_hist.GetBinContent(bin_idx)

            if not math.isfinite(cv_content) or abs(cv_content) < 1.0e-15:
                continue

            universe_content = universe_hist.GetBinContent(bin_idx)
            observed = universe_content / cv_content
            difference = abs(observed - expected)

            universe_max_diff = max(universe_max_diff, difference)
            global_max_abs_diff = max(global_max_abs_diff, difference)

            universe_bins += 1
            checked_bins += 1

        passed = universe_max_diff <= args.tolerance

        if not passed:
            failed_universes += 1

        if (
            output_u < args.print_universes
            or not passed
        ):
            print(
                "output u={:4d}  source u={:4d}  "
                "CSV={:.10f}  max|ratio-CSV|={:.3e}  "
                "bins={}  {}".format(
                    output_u,
                    source_u,
                    expected,
                    universe_max_diff,
                    universe_bins,
                    "PASS" if passed else "FAIL",
                )
            )

    print("\n[SUMMARY]")
    print("  checked bins       = {}".format(checked_bins))
    print("  failed universes   = {}".format(failed_universes))
    print(
        "  global max abs diff= {:.6e}".format(
            global_max_abs_diff
        )
    )

    if failed_universes == 0:
        print(
            "  RESULT             = PASS: Flux_BeamFocus matches "
            "the shifted CSV weights."
        )
        return 0

    print(
        "  RESULT             = FAIL: one or more universes do "
        "not match the shifted CSV."
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
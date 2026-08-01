#!/usr/bin/env python

import argparse
import csv
import math
import numpy as np


PDG_NAMES = {
    14: "nu_mu",
    -14: "nubar_mu",
    12: "nu_e",
    -12: "nubar_e",
}

ORDERED_PDGS = [14, -14, 12, -12]


def main():
    parser = argparse.ArgumentParser(
        description="Merge BeamFocus accumulator files."
    )

    parser.add_argument(
        "--input",
        action="append",
        required=True,
        help="Input NPZ accumulator. Repeat for multiple files."
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Final flat-weight CSV."
    )

    args = parser.parse_args()

    total_cv_sums = None
    total_universe_sums = None
    total_event_counts = None
    total_skipped_counts = None
    expected_n_universes = None
    expected_mode = None

    for filename in args.input:
        data = np.load(filename)

        n_universes = int(data["n_universes"])
        mode = str(data["mode"])
        pdgs = list(data["pdgs"].astype(int))

        if pdgs != ORDERED_PDGS:
            raise RuntimeError(
                "Unexpected PDG ordering in {}".format(filename)
            )

        if expected_n_universes is None:
            expected_n_universes = n_universes
            expected_mode = mode

            total_cv_sums = np.zeros_like(data["cv_sums"], dtype=float)
            total_universe_sums = np.zeros_like(
                data["universe_sums"],
                dtype=float
            )
            total_event_counts = np.zeros_like(
                data["event_counts"],
                dtype=np.int64
            )
            total_skipped_counts = np.zeros_like(
                data["skipped_counts"],
                dtype=np.int64
            )
        else:
            if n_universes != expected_n_universes:
                raise RuntimeError(
                    "Universe-count mismatch in {}".format(filename)
                )

            if mode != expected_mode:
                raise RuntimeError(
                    "Horn-mode mismatch in {}".format(filename)
                )

        total_cv_sums += data["cv_sums"]
        total_universe_sums += data["universe_sums"]
        total_event_counts += data["event_counts"]
        total_skipped_counts += data["skipped_counts"]

        print("Merged:", filename)

    flat_weights = np.ones_like(total_universe_sums, dtype=float)

    for ipdg, pdg in enumerate(ORDERED_PDGS):
        if total_cv_sums[ipdg] <= 0.0:
            print(
                "WARNING: no valid CV contribution for {}".format(
                    PDG_NAMES[pdg]
                )
            )
            continue

        flat_weights[ipdg, :] = (
            total_universe_sums[ipdg, :] /
            total_cv_sums[ipdg]
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

        for u in range(expected_n_universes):
            writer.writerow([
                u,
                "{:.10f}".format(flat_weights[0, u]),
                "{:.10f}".format(flat_weights[1, u]),
                "{:.10f}".format(flat_weights[2, u]),
                "{:.10f}".format(flat_weights[3, u]),
            ])

    print()
    print("Mode:", expected_mode)
    print("Output:", args.output)

    for ipdg, pdg in enumerate(ORDERED_PDGS):
        weights = flat_weights[ipdg]
        mean = np.mean(weights)
        rms = math.sqrt(np.mean((weights - mean) ** 2))

        print()
        print("{} ({:+d})".format(PDG_NAMES[pdg], pdg))
        print("  Events       :", total_event_counts[ipdg])
        print("  Skipped      :", total_skipped_counts[ipdg])
        print("  CV sum       : {:.8e}".format(total_cv_sums[ipdg]))
        print("  Universe mean: {:.8f}".format(mean))
        print("  Universe RMS : {:.8f}".format(rms))
        print("  Minimum      : {:.8f}".format(np.min(weights)))
        print("  Maximum      : {:.8f}".format(np.max(weights)))


if __name__ == "__main__":
    main()
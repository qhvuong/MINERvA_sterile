#!/usr/bin/env python3

import argparse
import csv
import glob
import os
import re
import sys


def start_index(path):
    m = re.search(r"_u(\d+)_(\d+)_nullOnly\.csv$", path)
    if m:
        return int(m.group(1))
    return 999999


def combine(pattern, outpath):
    files = sorted(glob.glob(pattern), key=start_index)

    print("Pattern:", pattern)
    print("Found files:", len(files))

    if len(files) == 0:
        print("ERROR: no files found")
        return False

    header = None
    nrows = 0

    os.makedirs(os.path.dirname(outpath), exist_ok=True)

    with open(outpath, "w", newline="") as fout:
        writer = None

        for path in files:
            with open(path, "r") as fin:
                reader = csv.DictReader(fin)

                if header is None:
                    header = reader.fieldnames
                    writer = csv.DictWriter(fout, fieldnames=header)
                    writer.writeheader()
                elif reader.fieldnames != header:
                    raise RuntimeError(f"Header mismatch in {path}")

                local = 0
                for row in reader:
                    writer.writerow(row)
                    nrows += 1
                    local += 1

                print(os.path.basename(path), "rows =", local)

    print("Wrote:", outpath)
    print("Total data rows:", nrows)
    print("Expected if complete: 2200")
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tag", default="prodNueel")
    # parser.add_argument("--mode", required=True, choices=["excludeRatio_maskRatio", "includeAll_maskNone"])
    parser.add_argument("--mode", required=True, choices=["excludeRatio_maskRatio", "includeAll_maskNone", "profileOnlyRatio_maskNonRatio"])
    parser.add_argument("--csv-dir", default="/exp/minerva/data/users/qvuong/surfaces/csvs")
    args = parser.parse_args()

    pattern = os.path.join(
        args.csv_dir,
        f"lambda_scan_fluxUniverseThrows_{args.tag}_{args.mode}_u*_nullOnly.csv",
    )

    outpath = os.path.join(
        args.csv_dir,
        f"lambda_scan_fluxUniverseThrows_{args.tag}_{args.mode}_allChunks_nullOnly.csv",
    )

    ok = combine(pattern, outpath)
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
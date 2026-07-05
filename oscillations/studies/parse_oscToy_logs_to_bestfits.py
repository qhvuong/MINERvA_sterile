#!/usr/bin/env python3

import os
import re
import glob
import argparse
import numpy as np

PATTERN = re.compile(
    r"toy dchi2\s*=\s*([0-9eE+\-.]+).*?"
    r"toy BF\s*=\s*dm2\s*([0-9eE+\-.]+),\s*"
    r"ue4\s*([0-9eE+\-.]+),\s*"
    r"umu4\s*([0-9eE+\-.]+),\s*"
    r"utau4\s*([0-9eE+\-.]+)",
    re.DOTALL,
)

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--log-dir", required=True)
    ap.add_argument("--mode", required=True)
    ap.add_argument("--outdir", default="/exp/minerva/data/users/qvuong/surfaces/npys")
    return ap.parse_args()

def main():
    args = parse_args()

    files = sorted(glob.glob(os.path.join(args.log_dir, "*.log")))
    print("n log files =", len(files))

    rows = []

    for f in files:
        with open(f, "r", errors="replace") as fp:
            txt = fp.read()

        matches = PATTERN.findall(txt)

        if len(matches) == 0:
            print("WARNING: no toy BF matches in", f)
            continue

        for m in matches:
            dchi2, dm2, ue4, umu4, utau4 = map(float, m)
            rows.append([dchi2, dm2, ue4, umu4, utau4])

    rows = np.asarray(rows, dtype=float)

    print("n parsed toys =", len(rows))

    if len(rows) == 0:
        raise RuntimeError("No toys parsed")

    os.makedirs(args.outdir, exist_ok=True)

    csv_out = os.path.join(args.outdir, "toy_bestfits_{}.csv".format(args.mode))
    npy_dchi2 = os.path.join(args.outdir, "asimov_deltachi2s_{}.npy".format(args.mode))
    npy_dm2 = os.path.join(args.outdir, "asimov_best_dm2s_{}.npy".format(args.mode))
    npy_ue4 = os.path.join(args.outdir, "asimov_best_ue4s_{}.npy".format(args.mode))
    npy_umu4 = os.path.join(args.outdir, "asimov_best_umu4s_{}.npy".format(args.mode))
    npy_utau4 = os.path.join(args.outdir, "asimov_best_utau4s_{}.npy".format(args.mode))

    np.savetxt(
        csv_out,
        rows,
        delimiter=",",
        header="dchi2,dm2,ue4,umu4,utau4",
        comments="",
    )

    np.save(npy_dchi2, rows[:, 0])
    np.save(npy_dm2, rows[:, 1])
    np.save(npy_ue4, rows[:, 2])
    np.save(npy_umu4, rows[:, 3])
    np.save(npy_utau4, rows[:, 4])

    print("Saved:", csv_out)
    print("Saved:", npy_dchi2)
    print("Saved:", npy_dm2)
    print("Saved:", npy_ue4)
    print("Saved:", npy_umu4)
    print("Saved:", npy_utau4)

    print("\nSummary:")
    print("dchi2 mean   =", np.mean(rows[:, 0]))
    print("dchi2 median =", np.median(rows[:, 0]))
    print("dchi2 95%    =", np.percentile(rows[:, 0], 95))
    print("dm2 median   =", np.median(rows[:, 1]))
    print("ue4 median   =", np.median(rows[:, 2]))
    print("umu4 median  =", np.median(rows[:, 3]))
    print("utau4 median =", np.median(rows[:, 4]))

if __name__ == "__main__":
    main()
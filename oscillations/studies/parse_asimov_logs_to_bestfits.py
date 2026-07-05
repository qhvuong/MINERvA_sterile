#!/usr/bin/env python3

import os
import re
import glob
import argparse
import numpy as np

NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"

SELECTED_RE = re.compile(
    r"===== selected best fit candidate =====.*?"
    r"best chi2\s*=\s*({num}).*?"
    r"best dm2\s*=\s*({num}).*?"
    r"best Ue4\s*=\s*({num}).*?"
    r"best Umu4\s*=\s*({num}).*?"
    r"best Utau4\s*=\s*({num})".format(num=NUM),
    re.DOTALL,
)

SAVED_RE = re.compile(
    r"Saved:\s*(\S*sample_dchi2s_\S+\.csv)"
)

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--log-dir", required=True)
    ap.add_argument("--text-dir", required=True)
    ap.add_argument("--mode", required=True)
    ap.add_argument("--outdir", default="/exp/minerva/data/users/qvuong/surfaces/npys")
    ap.add_argument("--strict", action="store_true", default=False)
    return ap.parse_args()

def load_dchi2_csv(path):
    arr = np.loadtxt(path, delimiter=",")
    arr = np.asarray(arr, dtype=float).ravel()
    arr = arr[np.isfinite(arr)]
    return arr

def parse_log_file(log_file, text_dir):
    with open(log_file, "r", errors="replace") as fp:
        txt = fp.read()

    # Parse selected BF blocks.
    bf_matches = SELECTED_RE.findall(txt)

    bf_rows = []
    for m in bf_matches:
        best_chi2, dm2, ue4, umu4, utau4 = map(float, m)
        bf_rows.append([best_chi2, dm2, ue4, umu4, utau4])

    # Find corresponding saved dchi2 CSV.
    saved_matches = SAVED_RE.findall(txt)

    if len(saved_matches) == 0:
        return [], "no_saved_csv"

    saved_basename = os.path.basename(saved_matches[-1])
    dchi2_file = os.path.join(text_dir, saved_basename)

    if not os.path.isfile(dchi2_file):
        return [], "missing_text_csv: {}".format(dchi2_file)

    dchi2s = load_dchi2_csv(dchi2_file)

    if len(bf_rows) != len(dchi2s):
        return [], "count_mismatch: nBF={} nDchi2={} file={}".format(
            len(bf_rows), len(dchi2s), os.path.basename(log_file)
        )

    rows = []
    for dchi2, bf in zip(dchi2s, bf_rows):
        best_chi2, dm2, ue4, umu4, utau4 = bf
        rows.append([dchi2, best_chi2, dm2, ue4, umu4, utau4])

    return rows, "ok"

def main():
    args = parse_args()

    logs = sorted(glob.glob(os.path.join(args.log_dir, "*.log")))
    print("n log files =", len(logs))

    all_rows = []
    status_counts = {}

    for log_file in logs:
        rows, status = parse_log_file(log_file, args.text_dir)

        status_counts[status] = status_counts.get(status, 0) + 1

        if status != "ok":
            print("WARNING:", status, "in", log_file)
            if args.strict:
                raise RuntimeError(status)

        all_rows.extend(rows)

    rows = np.asarray(all_rows, dtype=float)

    print("\nStatus counts:")
    for k, v in sorted(status_counts.items()):
        print("  {}: {}".format(k, v))

    print("\nn parsed toys =", len(rows))

    if len(rows) == 0:
        raise RuntimeError("No toys parsed")

    os.makedirs(args.outdir, exist_ok=True)

    csv_out = os.path.join(args.outdir, "toy_bestfits_{}.csv".format(args.mode))

    npy_dchi2 = os.path.join(args.outdir, "asimov_deltachi2s_{}.npy".format(args.mode))
    npy_chi2  = os.path.join(args.outdir, "asimov_best_chi2s_{}.npy".format(args.mode))
    npy_dm2   = os.path.join(args.outdir, "asimov_best_dm2s_{}.npy".format(args.mode))
    npy_ue4   = os.path.join(args.outdir, "asimov_best_ue4s_{}.npy".format(args.mode))
    npy_umu4  = os.path.join(args.outdir, "asimov_best_umu4s_{}.npy".format(args.mode))
    npy_utau4 = os.path.join(args.outdir, "asimov_best_utau4s_{}.npy".format(args.mode))

    np.savetxt(
        csv_out,
        rows,
        delimiter=",",
        header="dchi2,best_chi2,dm2,ue4,umu4,utau4",
        comments="",
    )

    np.save(npy_dchi2, rows[:, 0])
    np.save(npy_chi2,  rows[:, 1])
    np.save(npy_dm2,   rows[:, 2])
    np.save(npy_ue4,   rows[:, 3])
    np.save(npy_umu4,  rows[:, 4])
    np.save(npy_utau4, rows[:, 5])

    print("\nSaved:", csv_out)
    print("Saved:", npy_dchi2)
    print("Saved:", npy_chi2)
    print("Saved:", npy_dm2)
    print("Saved:", npy_ue4)
    print("Saved:", npy_umu4)
    print("Saved:", npy_utau4)

    print("\nSummary:")
    print("dchi2 mean   =", np.mean(rows[:, 0]))
    print("dchi2 median =", np.median(rows[:, 0]))
    print("dchi2 95%    =", np.percentile(rows[:, 0], 95))
    print("dm2 median   =", np.median(rows[:, 2]))
    print("ue4 median   =", np.median(rows[:, 3]))
    print("umu4 median  =", np.median(rows[:, 4]))
    print("utau4 median =", np.median(rows[:, 5]))

if __name__ == "__main__":
    main()
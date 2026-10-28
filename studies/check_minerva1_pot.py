
#!/usr/bin/env python3

import ROOT
from pathlib import Path
from collections import Counter

ROOT.gROOT.SetBatch(True)

BASE = Path("/exp/minerva/app/users/qvuong/MAT_AL9/CC-NuE-XSec")

PLAYLISTS = {
    "Data": BASE / "file_option/playlist_datale1p6.txt",
    "MC": BASE / "file_option/playlist_mcle1p6.txt",
}

P8 = {
    "Data": 9.605330838278293e19,
    "MC": 9.492519777376983e20,
}

MARIA_P6 = {
    "Data": 9.570e19,
    "MC": 9.300e20,
}


def read_playlist(path):
    with open(path) as f:
        return [
            line.strip()
            for line in f
            if line.strip() and not line.lstrip().startswith("#")
        ]


def get_file_pot(path):
    f = ROOT.TFile.Open(path, "READ")

    if not f or f.IsZombie():
        raise RuntimeError("Cannot open ROOT file")

    try:
        meta = f.Get("Meta")
        if not meta:
            raise RuntimeError("Missing Meta tree")

        for branch in ("POT_Used", "POT_Total"):
            if not meta.GetBranch(branch):
                raise RuntimeError(f"Missing {branch} branch")

        used = 0.0
        total = 0.0

        for entry in meta:
            used += float(entry.POT_Used)
            total += float(entry.POT_Total)

        return used, total, meta.GetEntries()
    finally:
        f.Close()


def check_playlist(label, playlist):
    paths = read_playlist(playlist)
    counts = Counter(paths)

    print(f"\n{'=' * 65}")
    print(f"{label}: {playlist}")
    print(f"{'=' * 65}")
    print("Playlist entries:", len(paths))
    print("Unique paths:", len(counts))

    duplicates = {
        path: count
        for path, count in counts.items()
        if count > 1
    }

    if duplicates:
        print("\nDUPLICATED PATHS:")
        for path, count in duplicates.items():
            print(f"  {count}x {path}")

    used_sum = 0.0
    total_sum = 0.0
    failed = []

    for i, path in enumerate(paths, 1):
        try:
            used, total, entries = get_file_pot(path)

            used_sum += used
            total_sum += total

            if i % 100 == 0 or i == len(paths):
                print(
                    f"[{i}/{len(paths)}] "
                    f"Cumulative Used={used_sum:.9e} "
                    f"Cumulative Total={total_sum:.9e}",
                    flush=True,
                )

        except Exception as exc:
            failed.append((path, str(exc)))
            print(f"[FAILED] {path}: {exc}", flush=True)

    if failed:
        print(f"\nWARNING: {len(failed)} files failed!")
        print("The calculated POT is incomplete.")

        for path, reason in failed[:10]:
            print(f"  {path}: {reason}")

    print(f"\n{'-' * 65}")
    print(f"{label} SUMMARY")
    print(f"{'-' * 65}")
    print(f"Files listed:     {len(paths)}")
    print(f"Failed files:     {len(failed)}")
    print(f"Summed POT_Used:  {used_sum:.15e}")
    print(f"Summed POT_Total: {total_sum:.15e}")

    print(f"Our p8 POT:    {P8[label]:.15e}")
    print(f"Maria p6 POT:  {MARIA_P6[label]:.15e}")

    if P8[label]:
        difference = 100 * (used_sum / P8[label] - 1)
        print(f"Difference from p8: {difference:+.6f}%")

    if MARIA_P6[label]:
        difference = 100 * (used_sum / MARIA_P6[label] - 1)
        print(f"Difference from Maria p6: {difference:+.6f}%")

    return used_sum, len(failed)


if __name__ == "__main__":
    data_pot, data_failed = check_playlist("Data", PLAYLISTS["Data"])
    mc_pot, mc_failed = check_playlist("MC", PLAYLISTS["MC"])

    print("\n========== FINAL POT NORMALIZATION ==========")

    if data_failed == 0 and mc_failed == 0 and mc_pot > 0:
        print(f"Our recalculated p6 Data/MC POT: {data_pot / mc_pot:.12f}")
    else:
        print("WARNING: Incomplete POT calculation due to failed files.")

    print(
        "Maria's p6 Data/MC POT:        "
        f"{MARIA_P6['Data'] / MARIA_P6['MC']:.12f}"
    )

    print(
        "Our p8 Data/MC POT:            "
        f"{P8['Data'] / P8['MC']:.12f}"
    )

    print("=============================================")
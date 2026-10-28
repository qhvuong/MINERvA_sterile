
#!/usr/bin/env python3

import argparse
import csv
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

BASE = Path("/exp/minerva/app/users/qvuong/MAT_AL9/CC-NuE-XSec")

PLAYLIST = BASE / "file_option/playlist_mcle1p8.txt"

LE_DIR = BASE / "custom_mhrwKineRenorm"

ME_DIR = Path(
    "/exp/minerva/app/users/qvuong/MAT_AL9/"
    "opt/lib/data/mhrwKineRenorm"
)

LE_FILENAME = (
    "Renorm_Kine_Truth_Minerva1_v10_HD_"
    "withElastics_RwgtNeutronCV_z5991-8408_a850.root"
)

ME_FILENAME = (
    "Renorm_Kine_Truth_MinervaME1A_v10_HD_"
    "withElastics_RwgtNeutronCV_z5991-8408_a850.root"
)

ENERGY_BINS = np.array([
    0, 0.5, 1, 1.5, 2, 2.5, 3, 4,
    5, 6, 8, 10, 12, 15, 20, 30, 50
], dtype=float)


def worker(mode, renorm_dir, output, max_events, max_files):
    import ROOT
    import PlotUtils

    # The analysis config parses sys.argv during import.
    # Temporarily hide this script's diagnostic arguments.
    saved_argv = sys.argv[:]

    try:
        sys.argv = [saved_argv[0]]
        from tools.SystematicsUniverse import CVSystematicUniverse
    finally:
        sys.argv = saved_argv

    ROOT.gROOT.SetBatch(True)

    files = [
        line.strip()
        for line in PLAYLIST.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]

    if max_files > 0:
        files = files[:max_files]

    if not files:
        raise RuntimeError("Empty MC playlist")

    chain = PlotUtils.ChainWrapper("MasterAnaDev")

    for filename in files:
        chain.Add(filename)

    # Configure before evaluating any hadron weights.
    MU = ROOT.PlotUtils.MinervaUniverse

    MU.SetPlaylist("minerva1")
    MU.SetNuEConstraint(False)
    MU.SetAnalysisNuPDG(14)

    MU.SetReadoutVolume("Tracker")
    MU.SetMHRWeightNeutronCVReweight(True)
    MU.SetMHRWeightElastics(True)

    result = MU.SetMHRWeightFilename(
        "Renorm_Kine_Truth",
        str(renorm_dir)
    )
    if result is False:
        raise RuntimeError(
            "SetMHRWeightFilename rejected configuration"
        )

    universe = CVSystematicUniverse(chain, 0)

    total = chain.GetEntries()
    if max_events > 0:
        total = min(total, max_events)

    print(
        f"[{mode}] files={len(files)} entries={total}",
        flush=True
    )

    accepted = 0

    with open(output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["entry", "enu_gev", "hadron_weight"])

        for i in range(total):
            universe.SetEntry(i)

            # True FHC numu CC events.
            # This is not the full reconstructed CCnumu selection.
            if int(universe.mc_incoming) != 14:
                continue

            if int(universe.GetInt("mc_current")) != 1:
                continue

            enu = float(universe.mc_incomingE) / 1000.0

            if not math.isfinite(enu) or enu <= 0:
                continue

            w = float(universe.GetGeantHadronWeight())

            if not math.isfinite(w):
                raise RuntimeError(
                    f"Nonfinite hadron weight: entry={i}, w={w}"
                )

            writer.writerow([i, enu, w])
            accepted += 1

            if accepted <= 10:
                print(
                    f"[{mode}] entry={i} "
                    f"Enu={enu:.4f} GeV w={w:.6f}",
                    flush=True
                )

            if i % 10000 == 0 and i > 0:
                print(
                    f"[{mode}] processed={i}",
                    flush=True
                )

    print(
        f"[{mode}] selected events={accepted}",
        flush=True
    )


def run_worker(mode, directory, output, args):
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker",
        "--mode", mode,
        "--renorm-dir", str(directory),
        "--output", str(output),
        "--max-events", str(args.max_events),
        "--max-files", str(args.max_files),
    ]

    subprocess.run(command, check=True)


def read_weights(filename):
    with open(filename) as f:
        rows = list(csv.DictReader(f))

    entries = np.array(
        [int(r["entry"]) for r in rows], dtype=np.int64
    )
    energies = np.array(
        [float(r["enu_gev"]) for r in rows]
    )
    weights = np.array(
        [float(r["hadron_weight"]) for r in rows]
    )

    return entries, energies, weights


def summarize(energies, le_weights, me_weights, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    centers = 0.5 * (ENERGY_BINS[:-1] + ENERGY_BINS[1:])

    le_means = []
    me_means = []
    counts = []

    summary = outdir / "hadron_cv_energy_summary.csv"

    with open(summary, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "enu_low", "enu_high", "n_events",
            "mean_le", "mean_me",
            "mean_le_minus_unity",
            "mean_me_minus_unity"
        ])

        for lo, hi in zip(ENERGY_BINS[:-1], ENERGY_BINS[1:]):
            mask = (energies >= lo) & (energies < hi)
            n = int(np.count_nonzero(mask))

            if n:
                le = float(np.mean(le_weights[mask]))
                me = float(np.mean(me_weights[mask]))
            else:
                le = np.nan
                me = np.nan

            le_means.append(le)
            me_means.append(me)
            counts.append(n)

            writer.writerow([
                lo, hi, n,
                le, me,
                le - 1, me - 1
            ])

    le_means = np.array(le_means)
    me_means = np.array(me_means)

    print("\n========== HADRON CV WEIGHT SUMMARY ==========")
    print(
        f"{'Enu (GeV)':>14} {'N':>8} "
        f"{'LE':>12} {'ME':>12}"
    )

    for i in range(len(centers)):
        print(
            f"{ENERGY_BINS[i]:5.1f}-"
            f"{ENERGY_BINS[i+1]:5.1f} "
            f"{counts[i]:8d} "
            f"{le_means[i]:12.6f} "
            f"{me_means[i]:12.6f}"
        )

    print("\nOverall means:")
    print("  Unity:", 1.0)
    print("  LE:", np.mean(le_weights))
    print("  ME:", np.mean(me_weights))
    print("  LE nonunity fraction:",
          np.mean(np.abs(le_weights - 1) > 1e-8))
    print("  ME nonunity fraction:",
          np.mean(np.abs(me_weights - 1) > 1e-8))

    plt.figure(figsize=(10, 6))

    plt.axhline(
        1.0, color="black", linestyle="--",
        linewidth=2, label="Unity (no hadron CV)"
    )

    plt.plot(
        centers, le_means,
        "o-", linewidth=2, markersize=5,
        label="LE renormalization (ours)"
    )

    plt.plot(
        centers, me_means,
        "s-", linewidth=2, markersize=5,
        label="ME1A renormalization on LE events"
    )

    plt.xlabel(r"True neutrino energy $E_\nu$ (GeV)")
    plt.ylabel(r"Mean CV hadron weight")
    plt.title("MINERVA1: LE vs ME hadron CV reweight")
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()

    plot = outdir / "hadron_cv_le_vs_me.png"
    plt.savefig(plot, dpi=180)
    plt.close()

    print("\nSaved plot:", plot)
    print("Saved summary:", summary)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--mode", choices=["LE", "ME"])
    parser.add_argument("--renorm-dir")
    parser.add_argument("--output")
    parser.add_argument("--max-events", type=int, default=20000)
    parser.add_argument("--max-files", type=int, default=1)
    parser.add_argument(
        "--outdir", default="studies/hadron_cv_comparison"
    )

    args = parser.parse_args()

    if args.worker:
        worker(
            args.mode, args.renorm_dir,
            args.output, args.max_events, args.max_files
        )
        return

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    le_file = LE_DIR / LE_FILENAME
    me_file = ME_DIR / ME_FILENAME

    for f in (le_file, me_file):
        if not f.is_file():
            raise FileNotFoundError(f)

    le_output = outdir / "events_le.csv"
    me_output = outdir / "events_me.csv"

    # Each subprocess sees one independent renormalization table.
    with tempfile.TemporaryDirectory(prefix="mhrw_compare_") as tmp:
        tmp = Path(tmp)

        le_tmp = tmp / "le"
        me_tmp = tmp / "me"

        le_tmp.mkdir()
        me_tmp.mkdir()

        (le_tmp / LE_FILENAME).symlink_to(le_file)
        (me_tmp / LE_FILENAME).symlink_to(me_file)

        print("\n=== Running LE configuration ===", flush=True)
        run_worker("LE", le_tmp, le_output, args)

        print("\n=== Running ME configuration ===", flush=True)
        run_worker("ME", me_tmp, me_output, args)

    le_entries, le_energy, le_weights = read_weights(le_output)
    me_entries, me_energy, me_weights = read_weights(me_output)

    if not np.array_equal(le_entries, me_entries):
        raise RuntimeError("LE and ME selected different entries")

    if not np.allclose(le_energy, me_energy, rtol=0, atol=1e-8):
        raise RuntimeError("LE and ME energies differ")

    if len(le_entries) == 0:
        raise RuntimeError("No events passed the truth selection")

    summarize(le_energy, le_weights, me_weights, outdir)


if __name__ == "__main__":
    main()

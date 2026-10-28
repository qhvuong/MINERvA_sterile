
#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path

import ROOT
import numpy as np
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

BASE = Path(
    "/exp/minerva/app/users/qvuong/MAT_AL9/CC-NuE-XSec"
)

LE_FILE = (
    BASE / "custom_mhrwKineRenorm"
    / "Renorm_Kine_Truth_Minerva1_v10_HD_"
      "withElastics_RwgtNeutronCV_z5991-8408_a850.root"
)

ME_FILE = Path(
    "/exp/minerva/app/users/qvuong/MAT_AL9/"
    "opt/lib/data/mhrwKineRenorm/"
    "Renorm_Kine_Truth_MinervaME1A_v10_HD_"
    "withElastics_RwgtNeutronCV_z5991-8408_a850.root"
)


def load_histograms(filename):
    f = ROOT.TFile.Open(str(filename), "READ")

    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {filename}")

    histograms = {}

    def visit(directory, prefix=""):
        for key in directory.GetListOfKeys():
            obj = key.ReadObj()
            name = prefix + obj.GetName()

            if obj.InheritsFrom("TDirectory"):
                visit(obj, name + "/")

            elif obj.InheritsFrom("TH2"):
                hist = obj.Clone()
                hist.SetDirectory(0)
                histograms[name] = hist

    visit(f)
    f.Close()

    return histograms


def extract_2d(hist):
    nx = hist.GetNbinsX()
    ny = hist.GetNbinsY()

    xedges = np.array([
        hist.GetXaxis().GetBinLowEdge(i)
        for i in range(1, nx + 2)
    ])

    yedges = np.array([
        hist.GetYaxis().GetBinLowEdge(i)
        for i in range(1, ny + 2)
    ])

    values = np.array([
        [
            hist.GetBinContent(ix, iy)
            for ix in range(1, nx + 1)
        ]
        for iy in range(1, ny + 1)
    ])

    return xedges, yedges, values


def compare_one(name, le, me, outdir, angles):
    le_x, le_y, le_z = extract_2d(le)
    me_x, me_y, me_z = extract_2d(me)


    # The UP histogram stores the average neutron CV weight.
    # PlotUtils divides the raw neutron weight by this factor.
    #
    # Therefore the actual multiplicative correction is 1/R.

    if name.endswith("_2112_up"):

        le_cv_name = name.replace("_up", "_cv")
        me_cv_name = le_cv_name

        le_cv = LE_HISTS[le_cv_name]
        me_cv = ME_HISTS[me_cv_name]

        _, _, le_counts = extract_2d(le_cv)
        _, _, me_counts = extract_2d(me_cv)

        min_entries = 3

        le_z = np.divide(
            1.0,
            le_z,
            out=np.ones_like(le_z),
            where=(le_z != 0) & (le_counts > min_entries)
        )

        me_z = np.divide(
            1.0,
            me_z,
            out=np.ones_like(me_z),
            where=(me_z != 0) & (me_counts > min_entries)
        )

    if (
        le_z.shape != me_z.shape
        or not np.allclose(le_x, me_x)
        or not np.allclose(le_y, me_y)
    ):
        print(f"[SKIP] Different binning: {name}")
        return

    # Histograms are expected to use:
    # X = hadron momentum, Y = hadron angle.
    # Print labels so this assumption can be checked.
    print(f"\nHistogram: {name}")
    print("  X axis:", le.GetXaxis().GetTitle())
    print("  Y axis:", le.GetYaxis().GetTitle())
    print("  Size:", le_z.shape)

    le_integral = np.sum(le_z)
    me_integral = np.sum(me_z)

    print(f"  LE integral: {le_integral:.6g}")
    print(f"  ME integral: {me_integral:.6g}")

    if le_integral <= 0 or me_integral <= 0:
        print("  [SKIP] Non-positive integral")
        return

    print(f"  ME/LE integral: {me_integral/le_integral:.6f}")

    le_norm = le_z / le_integral
    me_norm = me_z / me_integral

    print(
        "  Maximum normalized bin difference:",
        np.max(np.abs(le_norm - me_norm))
    )

    xcenters = 0.5 * (le_x[:-1] + le_x[1:])

    safe_name = re.sub(r"[^A-Za-z0-9_-]", "_", name)

    # --------------------------------------------------
    # 1. LE, ME, UNITY VS HADRON MOMENTUM
    # --------------------------------------------------

    fig, axes = plt.subplots(
        len(angles),
        1,
        figsize=(9, 3.6 * len(angles)),
        squeeze=False
    )

    for row, angle in enumerate(angles):
        ax = axes[row, 0]

        # Choose the ROOT bin containing the requested angle.
        iy_root = le.GetYaxis().FindFixBin(angle)

        iy_root = max(
            1, min(iy_root, le.GetNbinsY())
        )

        iy = iy_root - 1

        actual_lo = le_y[iy]
        actual_hi = le_y[iy + 1]

        ax.axhline(
            1.0,
            color="black",
            linestyle="--",
            linewidth=1.8,
            label="Unity"
        )

        ax.step(
            xcenters,
            le_z[iy],
            where="mid",
            linewidth=2,
            label="LE MINERVA1"
        )

        ax.step(
            xcenters,
            me_z[iy],
            where="mid",
            linewidth=2,
            label="ME1A"
        )

        ax.set_title(
            f"Hadron angle bin: "
            f"{actual_lo:.2f} to {actual_hi:.2f}"
        )

        ax.set_xlabel(
            le.GetXaxis().GetTitle()
            or "Hadron momentum (table units)"
        )

        ax.set_ylabel("Multiplicative neutron CV renormalization")

        ax.legend()
        ax.grid(alpha=0.25)

    fig.suptitle(name, fontsize=11, fontweight="bold")
    fig.tight_layout()

    curve_path = outdir / f"{safe_name}_curves.png"
    fig.savefig(curve_path, dpi=180)
    plt.close(fig)


    fig, axes = plt.subplots(
        len(angles),
        1,
        figsize=(9, 3.6 * len(angles)),
        squeeze=False
    )

    for row, angle in enumerate(angles):
        ax = axes[row, 0]

        iy_root = le.GetYaxis().FindFixBin(angle)
        iy_root = max(1, min(iy_root, le.GetNbinsY()))
        iy = iy_root - 1

        ax.step(
            xcenters,
            le_norm[iy],
            where="mid",
            linewidth=2,
            label="LE normalized"
        )

        ax.step(
            xcenters,
            me_norm[iy],
            where="mid",
            linewidth=2,
            label="ME normalized"
        )

        ax.set_title(
            f"Hadron angle: {le_y[iy]:.1f}–{le_y[iy+1]:.1f} deg"
        )

        ax.set_xlabel("Hadron momentum (MeV)")
        ax.set_ylabel("Fraction of total histogram integral")
        ax.grid(alpha=0.25)
        ax.legend()

    fig.suptitle(
        f"{name}: normalized LE vs ME",
        fontsize=11,
        fontweight="bold"
    )

    fig.tight_layout()

    normalized_path = outdir / f"{safe_name}_normalized.png"
    fig.savefig(normalized_path, dpi=180)
    plt.close(fig)

    print("  Saved:", normalized_path)

    # --------------------------------------------------
    # 2. 2D LE / ME RATIO
    # --------------------------------------------------

    ratio = np.divide(
        le_z,
        me_z,
        out=np.full_like(le_z, np.nan),
        where=(me_z != 0)
    )

    fig, ax = plt.subplots(figsize=(9, 6))

    mesh = ax.pcolormesh(
        le_x,
        le_y,
        np.ma.masked_invalid(ratio),
        shading="auto",
        cmap="coolwarm",
    )

    cb = fig.colorbar(mesh, ax=ax)
    cb.set_label("LE / ME renormalization factor")

    ax.set_xlabel(
        le.GetXaxis().GetTitle()
        or "Hadron momentum (table units)"
    )

    ax.set_ylabel(
        le.GetYaxis().GetTitle()
        or "Hadron angle (table units)"
    )

    ax.set_title(f"LE / ME: {name}")

    fig.tight_layout()

    ratio_path = outdir / f"{safe_name}_ratio2d.png"
    fig.savefig(ratio_path, dpi=180)
    plt.close(fig)

    # --------------------------------------------------
    # 3. CSV OF ALL HISTOGRAM BINS
    # --------------------------------------------------

    csv_path = outdir / f"{safe_name}_values.csv"

    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)

        writer.writerow([
            "momentum_low",
            "momentum_high",
            "angle_low",
            "angle_high",
            "le_factor",
            "me_factor",
            "unity",
            "le_over_me",
            "le_minus_unity",
            "me_minus_unity",
        ])

        for iy in range(le_z.shape[0]):
            for ix in range(le_z.shape[1]):
                writer.writerow([
                    le_x[ix],
                    le_x[ix + 1],
                    le_y[iy],
                    le_y[iy + 1],
                    le_z[iy, ix],
                    me_z[iy, ix],
                    1.0,
                    ratio[iy, ix],
                    le_z[iy, ix] - 1.0,
                    me_z[iy, ix] - 1.0,
                ])

    # Summary: differences are unweighted across table bins.
    valid = (
        np.isfinite(le_z)
        & np.isfinite(me_z)
        & (me_z != 0)
    )

    if np.any(valid):
        difference = le_z[valid] - me_z[valid]

        print(
            "  Mean LE - ME (unweighted table bins):",
            np.mean(difference)
        )

        print(
            "  Max |LE - ME|:",
            np.max(np.abs(difference))
        )

        print(
            "  Max |LE/ME - 1|:",
            np.max(np.abs(ratio[valid] - 1))
        )

    print("  Saved:", curve_path)
    print("  Saved:", ratio_path)
    print("  Saved:", csv_path)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--le-file",
        default=str(LE_FILE)
    )

    parser.add_argument(
        "--me-file",
        default=str(ME_FILE)
    )

    parser.add_argument(
        "--match",
        default=r"^kine_rw_tracker_p_theta_12112_up$",
        help="Regex matching histogram names"
    )

    parser.add_argument(
        "--angles",
        default="15,45,90",
        help="Hadron angles in table Y-axis units"
    )

    parser.add_argument(
        "--max-hists",
        type=int,
        default=4
    )

    parser.add_argument(
        "--list",
        action="store_true",
        help="List common TH2 histograms without plotting"
    )

    parser.add_argument(
        "--outdir",
        default="studies/mhrw_renorm_comparison"
    )

    args = parser.parse_args()

    angles = [
        float(v)
        for v in args.angles.split(",")
    ]

    print("Loading LE:", args.le_file)
    le_hists = load_histograms(args.le_file)

    print("Loading ME:", args.me_file)
    me_hists = load_histograms(args.me_file)


    global LE_HISTS, ME_HISTS

    LE_HISTS = le_hists
    ME_HISTS = me_hists

    common = sorted(
        set(le_hists) & set(me_hists)
    )

    print("\nLE TH2 histograms:", len(le_hists))
    print("ME TH2 histograms:", len(me_hists))
    print("Common TH2 histograms:", len(common))

    if args.list:
        for name in common:
            h = le_hists[name]
            print(
                name,
                "| X:", h.GetXaxis().GetTitle(),
                "| Y:", h.GetYaxis().GetTitle()
            )
        return

    selected = [
        name
        for name in common
        if re.search(args.match, name, re.IGNORECASE)
    ]

    if not selected:
        print("\nNo histograms matched:", args.match)
        print("Run with --list to inspect available names.")
        return

    if args.max_hists > 0:
        selected = selected[:args.max_hists]

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("\nSelected histograms:")
    for name in selected:
        print(" ", name)

    for name in selected:
        compare_one(
            name,
            le_hists[name],
            me_hists[name],
            outdir,
            angles
        )

    print("\nFinished.")


if __name__ == "__main__":
    main()

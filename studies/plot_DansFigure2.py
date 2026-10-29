
#!/usr/bin/env python3
"""Figure-2-style pT panels of Eavail for a chosen p_parallel slice.

Input histogram: Eavail_Pt_Ppar_3D
  X = Eavail
  Y = pT
  Z = p_parallel

Plots triple-differential event densities:
  d^3N / (dEavail dpT dp_parallel)

Numerical densities are multiplied by 1e-5 for display.
Individual MC components are unstacked.

Layout:
  4 + 4 + 1 panels with shared axes and no spacing.
"""

import argparse
from pathlib import Path

import numpy as np
import ROOT
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)


# ============================================================
# CONFIGURATION
# ============================================================

BASE = "/exp/minerva/data/users/qvuong/nu_mu/"

TAG = "ValidateDataMC_CCnumu_DansCuts_v4p4p1_LE_MHRW"

DEFAULT_DATA = BASE + f"kin_dist_dataleFHC_{TAG}_MAD.root"
DEFAULT_MC = BASE + f"kin_dist_mcleFHC_{TAG}_MAD.root"

OUTPUT_DIR = Path("/exp/minerva/data/users/qvuong/studies")

STEM = "Eavail_Pt_Ppar_3D"

SCALE = 1e-5

DAN_LE_PANEL_SCALES = [
    4.0, 1.2, 1.3, 1.4,
    1.5, 1.9, 2.9, 6.7, 18.4
]

SIGNAL = [
    ("QE", "CCNuMuQE", "#27823d"),
    ("Resonant", "CCNuMuDelta", "#3979c4"),
    ("2p2h", "CCNuMu2p2h", "#df9229"),
]

BACKGROUND = [
    "CCNuMuDIS",
    "CCNuMu",
    "CC",
    "Other",
    "NuEElastic",
    "NCDiff",
    "CCNuMuWrongSign",
    "NC",
]

# Plot layout: four + four + one.
PANEL_POSITIONS = [
    (0, 0), (0, 1), (0, 2), (0, 3),
    (1, 0), (1, 1), (1, 2), (1, 3),
    (2, 0),
]

# Lowest occupied row for each column.
BOTTOM_ROW_BY_COLUMN = {
    0: 2,
    1: 1,
    2: 1,
    3: 1,
}


# ============================================================
# HISTOGRAM UTILITIES
# ============================================================

def get_th3(f, name):
    h = f.Get(name)

    if not h or not h.InheritsFrom("TH3"):
        raise RuntimeError(
            f"Missing TH3: {name} in {f.GetName()}"
        )

    return h


def axis_edges(axis):
    return np.array(
        [
            axis.GetBinLowEdge(i)
            for i in range(1, axis.GetNbins() + 2)
        ],
        dtype=float,
    )


def density_slice(h, pt_bin, ppar_bin):
    """Return triple-differential bin density and its statistical error."""

    ex = axis_edges(h.GetXaxis())

    volume = (
        np.diff(ex)
        * h.GetYaxis().GetBinWidth(pt_bin)
        * h.GetZaxis().GetBinWidth(ppar_bin)
    )

    vals = np.array(
        [
            h.GetBinContent(i, pt_bin, ppar_bin)
            for i in range(1, h.GetNbinsX() + 1)
        ],
        dtype=float,
    )

    errs = np.array(
        [
            h.GetBinError(i, pt_bin, ppar_bin)
            for i in range(1, h.GetNbinsX() + 1)
        ],
        dtype=float,
    )

    return vals / volume, errs / volume


# ============================================================
# MAIN
# ============================================================

def main():

    ap = argparse.ArgumentParser(description=__doc__)

    ap.add_argument("--data-file", default=DEFAULT_DATA)
    ap.add_argument("--mc-file", default=DEFAULT_MC)

    ap.add_argument(
        "--mc-scale",
        type=float,
        required=True,
        help="Data POT divided by MC POT",
    )

    ap.add_argument(
        "--ppar-bin",
        type=int,
        default=3,
        help="1-based p_parallel bin; default 3 = 2.5-3.0 GeV/c",
    )

    ap.add_argument(
        "--out",
        default=None,
        help="Output PNG path; default includes analysis tag",
    )

    ap.add_argument(
        "--max-eavail",
        type=float,
        default=0.8,
        help="Maximum Eavail displayed in GeV",
    )

    ap.add_argument(
        "--ymax",
        type=float,
        default=6.0,
        help="Common y-axis maximum (default: 6)",
    )

    ap.add_argument(
        "--panel-scales",
        type=float,
        nargs=9,
        default=None,
        metavar="FACTOR",
        help="Override nine panel display multipliers",
    )

    ap.add_argument(
        "--no-panel-scale",
        action="store_true",
        help="Disable panel scaling",
    )

    args = ap.parse_args()

    if args.no_panel_scale and args.panel_scales is not None:
        ap.error(
            "--no-panel-scale cannot be combined with --panel-scales"
        )

    if args.no_panel_scale:
        args.panel_scales = [1.0] * 9
    elif args.panel_scales is None:
        args.panel_scales = DAN_LE_PANEL_SCALES.copy()

    if (
        args.mc_scale <= 0
        or args.ymax <= 0
        or any(x <= 0 for x in args.panel_scales)
    ):
        ap.error("MC scale, ymax and panel scales must be positive")

    # --------------------------------------------------------
    # Open ROOT files
    # --------------------------------------------------------

    data_file = ROOT.TFile.Open(args.data_file, "READ")
    mc_file = ROOT.TFile.Open(args.mc_file, "READ")

    if (
        not data_file or data_file.IsZombie()
        or not mc_file or mc_file.IsZombie()
    ):
        raise RuntimeError("Cannot open ROOT data/MC files")

    try:

        data = get_th3(data_file, STEM)
        total = get_th3(mc_file, STEM)

        signals = [
            (
                label,
                get_th3(mc_file, STEM + "_" + suffix),
                color,
            )
            for label, suffix, color in SIGNAL
        ]

        backgrounds = [
            get_th3(mc_file, STEM + "_" + suffix)
            for suffix in BACKGROUND
        ]

        e_edges = axis_edges(data.GetXaxis())
        pt_edges = axis_edges(data.GetYaxis())
        par_edges = axis_edges(data.GetZaxis())

        npt = len(pt_edges) - 1

        if npt != 9:
            raise RuntimeError(
                "4-4-1 layout expects nine pT bins"
            )

        if not 1 <= args.ppar_bin < len(par_edges):
            ap.error(
                f"--ppar-bin must be 1 to {len(par_edges) - 1}"
            )

        if (
            args.max_eavail is not None
            and args.max_eavail <= e_edges[0]
        ):
            ap.error("--max-eavail must exceed minimum Eavail")

        # ----------------------------------------------------
        # Validate histogram binning
        # ----------------------------------------------------

        all_mc_hists = (
            [total]
            + [h for _, h, _ in signals]
            + backgrounds
        )

        for h in all_mc_hists:

            comparisons = (
                (e_edges, h.GetXaxis()),
                (pt_edges, h.GetYaxis()),
                (par_edges, h.GetZaxis()),
            )

            for reference, axis in comparisons:

                candidate = axis_edges(axis)

                if (
                    len(reference) != len(candidate)
                    or not np.allclose(
                        reference,
                        candidate,
                        rtol=0,
                        atol=1e-10,
                    )
                ):
                    raise RuntimeError(
                        f"Histogram binning differs: {h.GetName()}"
                    )

        # ----------------------------------------------------
        # Create compact 4-4-1 layout
        # ----------------------------------------------------

        fig = plt.figure(figsize=(15.5, 10.5))

        gs = fig.add_gridspec(
            3,
            4,
            left=0.095,
            right=0.985,
            bottom=0.105,
            top=0.88,
            hspace=0.0,
            wspace=0.0,
        )

        axes = []

        for i, (row, col) in enumerate(PANEL_POSITIONS):

            if i == 0:
                ax = fig.add_subplot(gs[row, col])
            else:
                ax = fig.add_subplot(
                    gs[row, col],
                    sharex=axes[0],
                    sharey=axes[0],
                )

            axes.append(ax)

        e_centers = 0.5 * (
            e_edges[:-1] + e_edges[1:]
        )

        max_mismatch = 0.0
        max_displayed_value = 0.0

        # ----------------------------------------------------
        # Plot each pT bin
        # ----------------------------------------------------

        for j, ax in enumerate(axes, start=1):

            row, col = PANEL_POSITIONS[j - 1]
            factor = args.panel_scales[j - 1]

            # Data
            d, de = density_slice(
                data,
                j,
                args.ppar_bin,
            )

            # Total MC
            t, _ = density_slice(
                total,
                j,
                args.ppar_bin,
            )

            # Scale to display units.
            d = d * SCALE * factor
            de = de * SCALE * factor

            t = t * args.mc_scale * SCALE * factor

            # Individual signal components
            components = []

            for label, h, color in signals:

                v, _ = density_slice(
                    h,
                    j,
                    args.ppar_bin,
                )

                v = (
                    v
                    * args.mc_scale
                    * SCALE
                    * factor
                )

                components.append(
                    (label, v, color)
                )

            # Combined background
            bg = sum(
                (
                    density_slice(
                        h, j, args.ppar_bin
                    )[0]
                    for h in backgrounds
                ),
                np.zeros(len(e_centers)),
            )

            bg = (
                bg
                * args.mc_scale
                * SCALE
                * factor
            )

            # Check component sum
            predicted = bg + sum(
                (v for _, v, _ in components),
                np.zeros_like(bg),
            )

            mismatch = float(
                np.max(np.abs(predicted - t))
            )

            max_mismatch = max(
                max_mismatch,
                mismatch,
            )

            if not np.allclose(
                predicted,
                t,
                rtol=1e-4,
                atol=1e-7,
            ):
                print(
                    f"[WARN] pT panel {j}: "
                    "signal + background differs "
                    "from inclusive MC"
                )

            # ------------------------------------------------
            # Plot curves
            # ------------------------------------------------

            ax.stairs(
                t,
                e_edges,
                color="#bb2a27",
                linewidth=1.8,
                zorder=4,
            )

            for _, v, color in components:
                ax.stairs(
                    v,
                    e_edges,
                    color=color,
                    linewidth=1.5,
                    zorder=3,
                )

            ax.stairs(
                bg,
                e_edges,
                color="#a74491",
                linewidth=1.5,
                zorder=3,
            )

            ax.errorbar(
                e_centers,
                d,
                yerr=de,
                color="black",
                fmt="o",
                ms=3.1,
                capsize=0,
                linewidth=1,
                zorder=5,
            )

            max_displayed_value = max(
                max_displayed_value,
                float(np.max(t)),
                float(np.max(d + de)),
            )

            # ------------------------------------------------
            # pT label INSIDE the panel
            # ------------------------------------------------

            pt_label = (
                rf"${pt_edges[j-1]:g} < p_T < "
                rf"{pt_edges[j]:g}$ GeV/c"
            )

            ax.text(
                0.04,
                0.94,
                pt_label,
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=10,
            )

            # ------------------------------------------------
            # Panel display multiplier
            # ------------------------------------------------

            if not args.no_panel_scale:

                ax.text(
                    0.97,
                    0.80,
                    rf"$\times {factor:g}$",
                    transform=ax.transAxes,
                    ha="right",
                    va="top",
                    fontsize=10,
                )

            # ------------------------------------------------
            # Shared axis configuration
            # ------------------------------------------------

            ax.set_xlim(
                e_edges[0],
                args.max_eavail
                if args.max_eavail is not None
                else e_edges[-1],
            )

            ax.set_ylim(0, args.ymax)

            # Consistent y ticks across panels
            ax.set_yticks(
                np.linspace(0, args.ymax, 6)
            )

            # No repeated axis titles
            ax.set_xlabel("")
            ax.set_ylabel("")

            # Only left column shows y tick labels
            show_y = (col == 0)

            # Only lowest panel in each column
            # shows x tick labels
            show_x = (
                row == BOTTOM_ROW_BY_COLUMN[col]
            )

            ax.tick_params(
                axis="y",
                which="both",
                left=show_y,
                labelleft=show_y,
                labelsize=9,
                direction="in",
            )

            ax.tick_params(
                axis="x",
                which="both",
                bottom=show_x,
                labelbottom=show_x,
                labelsize=9,
                direction="in",
            )

            # Thin borders so adjacent panels
            # visually form one grid.
            for spine in ax.spines.values():
                spine.set_linewidth(0.8)

            ax.grid(
                alpha=0.12,
                axis="y",
            )

            # ------------------------------------------------
            # Diagnostic event counts
            # ------------------------------------------------

            volume = (
                np.diff(e_edges)
                * (
                    pt_edges[j]
                    - pt_edges[j - 1]
                )
                * (
                    par_edges[args.ppar_bin]
                    - par_edges[args.ppar_bin - 1]
                )
            )

            # Undo display scaling before integration.
            nev_data = np.sum(
                (d / (SCALE * factor)) * volume
            )

            nev_mc = np.sum(
                (t / (SCALE * factor)) * volume
            )

            ratio = (
                nev_data / nev_mc
                if nev_mc else float("nan")
            )

            print(
                f"pT {pt_edges[j-1]:g}-{pt_edges[j]:g}: "
                f"data={nev_data:.6g}, "
                f"POT-scaled MC={nev_mc:.6g}, "
                f"D/MC={ratio:.4g}"
            )

        # ----------------------------------------------------
        # Shared figure labels
        # ----------------------------------------------------

        fig.supxlabel(
            r"$E_{\mathrm{avail}}$ (GeV)",
            fontsize=14,
            y=0.055,
        )

        fig.supylabel(
            r"Event rate ($10^{-5}$ events / (GeV$^3$/$c^2$))",
            fontsize=13,
            x=0.025,
        )

        # ----------------------------------------------------
        # Legend in empty bottom-row area
        # ----------------------------------------------------

        legend = [
            Line2D(
                [0], [0],
                color="black",
                marker="o",
                linestyle="None",
                label="Data",
            ),
            Line2D(
                [0], [0],
                color="#bb2a27",
                lw=2,
                label="Total MC",
            ),
        ]

        legend += [
            Line2D(
                [0], [0],
                color=color,
                lw=2,
                label=label,
            )
            for label, _, color in SIGNAL
        ]

        legend.append(
            Line2D(
                [0], [0],
                color="#a74491",
                lw=2,
                label="Background",
            )
        )

        fig.legend(
            handles=legend,
            loc="center",
            bbox_to_anchor=(0.64, 0.20),
            fontsize=12,
            frameon=False,
            ncol=2,
        )

        # ----------------------------------------------------
        # Figure title
        # ----------------------------------------------------

        ppar_low = par_edges[args.ppar_bin - 1]
        ppar_high = par_edges[args.ppar_bin]

        fig.suptitle(
            "Reconstructed available energy in muon "
            "$p_T$ bins\n"
            + rf"${ppar_low:g} < p_{{\parallel}} < "
            + rf"{ppar_high:g}$ GeV/c",
            fontsize=15,
            y=0.96,
        )

        # ----------------------------------------------------
        # Save plots
        # ----------------------------------------------------

        if args.out is None:
            out = (
                OUTPUT_DIR
                / f"dans_figure2_reco_{TAG}.png"
            )
        else:
            out = Path(args.out)

        out.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        fig.savefig(
            out,
            dpi=180,
        )

        fig.savefig(
            out.with_suffix(".pdf"),
        )

        plt.close(fig)

        print(
            "Saved:",
            out,
            "and",
            out.with_suffix(".pdf"),
        )

        print(
            "Maximum displayed category-sum discrepancy:",
            max_mismatch,
        )

        if max_displayed_value > args.ymax:
            print(
                "[WARNING] Some displayed values exceed "
                f"the common y-axis maximum {args.ymax:g}. "
                "These values will be clipped."
            )

    finally:
        data_file.Close()
        mc_file.Close()


if __name__ == "__main__":
    main()

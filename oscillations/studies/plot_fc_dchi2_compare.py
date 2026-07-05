#!/usr/bin/env python3

import os
import argparse
import numpy as np

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

DEFAULT_DATA_DCHI2 = {
    # Nprof=40 values from nominal no-throw fits
    "prodNueel_profiledFlux_includeAll_Nprof40": 10.275,
    "prodNueel_profiledFlux_excluderatio_Nprof40": 42.832,
    "prodNueel_noRatio_profiledFlux_includeAll_Nprof40": 8.650,

    # optional comparison
    "prodNueel_profiledFlux_includeAll_Nprof30": 11.639,
    "prodNueel_profiledFlux_includeAll_Nprof35": 11.580,
    "prodNueel_profiledFlux_includeAll_Nprof50": 9.193,
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Overlay multiple FC toy delta-chi2 distributions."
    )

    parser.add_argument(
        "--input-dir",
        default="/exp/minerva/data/users/qvuong/surfaces/chi2s",
        help="Directory containing asimov_deltachi2s_<mode>.npy files.",
    )

    parser.add_argument(
        "--modes",
        nargs="+",
        required=True,
        help="Mode names, without asimov_deltachi2s_ prefix and .npy suffix.",
    )

    parser.add_argument(
        "--labels",
        nargs="+",
        default=None,
        help="Legend labels. Must match number of modes if provided.",
    )

    parser.add_argument(
        "--data-dchi2",
        nargs="*",
        default=None,
        help=(
            "Optional data dchi2 values. Either one per mode, or omit to use "
            "built-in DEFAULT_DATA_DCHI2 when available. Use 'none' to skip a line."
        ),
    )

    parser.add_argument(
        "--bins",
        type=int,
        default=60,
        help="Number of histogram bins.",
    )

    parser.add_argument(
        "--density",
        action="store_true",
        default=False,
        help="Plot normalized distributions instead of toy counts.",
    )

    parser.add_argument(
        "--xmax",
        type=float,
        default=None,
        help="Optional maximum x value for plotting.",
    )

    parser.add_argument(
        "--outdir",
        default="/exp/minerva/data/users/qvuong/surfaces/plots/fc_compare",
        help="Output plot directory.",
    )

    parser.add_argument(
        "--outname",
        default="FC_dchi2_distribution_compare.png",
        help="Output plot filename.",
    )

    return parser.parse_args()


def load_results(input_dir, mode):
    f = os.path.join(input_dir, "asimov_deltachi2s_{}.npy".format(mode))

    if not os.path.isfile(f):
        raise IOError("Missing FC toy file: {}".format(f))

    arr = np.load(f)
    arr = np.asarray(arr).ravel()
    arr = arr[np.isfinite(arr)]

    if arr.size == 0:
        raise RuntimeError("No finite toy values in {}".format(f))

    return arr, f


def get_data_dchi2_values(args):
    if args.data_dchi2 is None:
        vals = []
        for mode in args.modes:
            vals.append(DEFAULT_DATA_DCHI2.get(mode, None))
        return vals

    if len(args.data_dchi2) != len(args.modes):
        raise ValueError("--data-dchi2 must have same length as --modes")

    vals = []
    for x in args.data_dchi2:
        if str(x).strip().lower() in ["none", "nan", "skip"]:
            vals.append(None)
        else:
            vals.append(float(x))
    return vals


def summarize(label, arr, data_val=None):
    print("\n===== {} =====".format(label))
    print("ntoys  =", arr.size)
    print("min    =", np.min(arr))
    print("max    =", np.max(arr))
    print("mean   =", np.mean(arr))
    print("median =", np.median(arr))

    for cl in [68, 90, 95, 99]:
        print("{}% critical = {}".format(cl, np.percentile(arr, cl)))

    if data_val is not None:
        p_fc = 100.0 * np.mean(arr >= data_val)
        print("data dchi2 =", data_val)
        print("toy fraction >= data =", p_fc, "%")


def main():
    args = parse_args()

    if args.labels is not None and len(args.labels) != len(args.modes):
        raise ValueError("--labels must have same length as --modes")

    labels = args.labels if args.labels is not None else args.modes
    data_vals = get_data_dchi2_values(args)

    loaded = []
    for mode, label, data_val in zip(args.modes, labels, data_vals):
        arr, f = load_results(args.input_dir, mode)
        print("Loaded:", f)
        summarize(label, arr, data_val=data_val)
        loaded.append((mode, label, arr, data_val))

    os.makedirs(args.outdir, exist_ok=True)
    outpath = os.path.join(args.outdir, args.outname)

    # Use common range so all histograms are comparable.
    all_values = np.concatenate([x[2] for x in loaded])
    xmin = 0.0
    xmax = args.xmax if args.xmax is not None else np.nanpercentile(all_values, 99.5)
    hist_range = (xmin, xmax)

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    fig, ax = plt.subplots(figsize=(10, 6))

    case_handles = []

    # Draw histograms and vertical lines.
    for i, (mode, label, arr, data_val) in enumerate(loaded):
        crit95 = np.percentile(arr, 95)
        color = colors[i % len(colors)]

        ax.hist(
            arr,
            bins=args.bins,
            range=hist_range,
            histtype="step",
            linewidth=1.9,
            density=args.density,
            color=color,
        )

        if data_val is not None:
            ax.axvline(
                data_val,
                color=color,
                linestyle="--",
                linewidth=2.0,
                alpha=0.95,
            )

        ax.axvline(
            crit95,
            color=color,
            linestyle=":",
            linewidth=2.3,
            alpha=0.95,
        )

        # Legend only says which color corresponds to which case.
        case_handles.append(
            Line2D(
                [0],
                [0],
                color=color,
                lw=2.0,
                linestyle="-",
                label=label,
            )
        )

    ylabel = "Density" if args.density else "Number of toys"
    ax.set_xlabel(r"Toy $\Delta\chi^2$", fontsize=13)
    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_title(r"FC toy $\Delta\chi^2$ distributions", fontsize=14)
    ax.grid(True, alpha=0.25)

    # Add labels near vertical lines after y-limits are known.
    ymin, ymax = ax.get_ylim()

    for i, (mode, label, arr, data_val) in enumerate(loaded):
        crit95 = np.percentile(arr, 95)
        color = colors[i % len(colors)]

        # Stagger labels vertically to reduce overlap.
        y_data = ymax * (0.92 - 0.10 * (i % 3))
        y_crit = ymax * (0.55 - 0.10 * (i % 3))

        if data_val is not None:
            ax.text(
                data_val,
                y_data,
                "data {:.2f}".format(data_val),
                color=color,
                rotation=90,
                va="top",
                ha="right",
                fontsize=10,
                fontweight="bold",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=2.0),
            )

        ax.text(
            crit95,
            y_crit,
            "95% {:.2f}".format(crit95),
            color=color,
            rotation=90,
            va="top",
            ha="left",
            fontsize=10,
            fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=2.0),
        )

    ax.legend(
        handles=case_handles,
        title="Case",
        fontsize=10,
        title_fontsize=11,
        loc="upper right",
        # bbox_to_anchor=(0.82, 0.98),
        framealpha=0.5,
    )

    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)

    print("\nSaved:", outpath)


if __name__ == "__main__":
    main()
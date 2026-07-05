#!/usr/bin/env python3

import os
import argparse
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot toy dchi2 vs toy best-fit parameters as 2D COLZ-like histograms."
    )

    parser.add_argument(
        "--input-dir",
        default="/exp/minerva/data/users/qvuong/surfaces/npys",
        help="Directory containing asimov_deltachi2s_<mode>.npy and asimov_best_*_<mode>.npy",
    )

    parser.add_argument(
        "--modes",
        nargs="+",
        required=True,
        help="Mode names, e.g. prodNueel_profiledFlux_includeAll_Nprof40",
    )

    parser.add_argument(
        "--labels",
        nargs="+",
        default=None,
        help="Optional plot titles. Must match number of modes if provided.",
    )

    parser.add_argument(
        "--outdir",
        default="/exp/minerva/data/users/qvuong/surfaces/plots/fc_compare",
        help="Output directory",
    )

    parser.add_argument(
        "--outname-prefix",
        default="FC_dchi2_colz",
        help="Prefix for output file names",
    )

    parser.add_argument(
        "--xbins",
        type=int,
        default=50,
        help="Number of x bins in each 2D histogram",
    )

    parser.add_argument(
        "--ybins",
        type=int,
        default=50,
        help="Number of y bins in each 2D histogram",
    )

    parser.add_argument(
        "--ymax",
        type=float,
        default=None,
        help="Optional fixed ymax for dchi2 axis",
    )

    parser.add_argument(
        "--y-percentile-max",
        type=float,
        default=99.5,
        help="If --ymax is not given, use this percentile of dchi2 as ymax",
    )

    parser.add_argument(
        "--logz",
        action="store_true",
        default=False,
        help="Use log color scale",
    )

    parser.add_argument(
        "--log-dm2",
        action="store_true",
        default=False,
        help="Use log x-axis for dm2 panel",
    )

    return parser.parse_args()


def safe_range(x):
    x = np.asarray(x)
    xmin = np.nanmin(x)
    xmax = np.nanmax(x)

    if not np.isfinite(xmin) or not np.isfinite(xmax):
        return (0.0, 1.0)

    if xmin == xmax:
        pad = 1e-6 if xmin == 0 else abs(xmin) * 0.05
        return (xmin - pad, xmax + pad)

    pad = 0.02 * (xmax - xmin)
    return (xmin - pad, xmax + pad)


def load_mode(input_dir, mode):
    f_dchi2 = os.path.join(input_dir, f"asimov_deltachi2s_{mode}.npy")
    f_dm2   = os.path.join(input_dir, f"asimov_best_dm2s_{mode}.npy")
    f_ue4   = os.path.join(input_dir, f"asimov_best_ue4s_{mode}.npy")
    f_umu4  = os.path.join(input_dir, f"asimov_best_umu4s_{mode}.npy")
    f_utau4 = os.path.join(input_dir, f"asimov_best_utau4s_{mode}.npy")

    needed = [f_dchi2, f_dm2, f_ue4, f_umu4, f_utau4]
    for f in needed:
        if not os.path.isfile(f):
            raise IOError(f"Missing file: {f}")

    dchi2 = np.asarray(np.load(f_dchi2)).ravel()
    dm2   = np.asarray(np.load(f_dm2)).ravel()
    ue4   = np.asarray(np.load(f_ue4)).ravel()
    umu4  = np.asarray(np.load(f_umu4)).ravel()
    utau4 = np.asarray(np.load(f_utau4)).ravel()

    n = len(dchi2)
    if not (len(dm2) == len(ue4) == len(umu4) == len(utau4) == n):
        raise RuntimeError(f"Length mismatch in mode {mode}")

    mask = (
        np.isfinite(dchi2)
        & np.isfinite(dm2)
        & np.isfinite(ue4)
        & np.isfinite(umu4)
        & np.isfinite(utau4)
    )

    return {
        "dchi2": dchi2[mask],
        "dm2": dm2[mask],
        "ue4": ue4[mask],
        "umu4": umu4[mask],
        "utau4": utau4[mask],
    }


def summarize(mode, label, arrs):
    dchi2 = arrs["dchi2"]
    dm2   = arrs["dm2"]
    ue4   = arrs["ue4"]
    umu4  = arrs["umu4"]
    utau4 = arrs["utau4"]

    print(f"\n===== {label} =====")
    print("mode          =", mode)
    print("ntoys         =", len(dchi2))
    print("dchi2 mean    =", np.mean(dchi2))
    print("dchi2 median  =", np.median(dchi2))
    print("dchi2 95%     =", np.percentile(dchi2, 95))
    print("dm2 median    =", np.median(dm2))
    print("ue4 median    =", np.median(ue4))
    print("umu4 median   =", np.median(umu4))
    print("utau4 median  =", np.median(utau4))


def make_one_plot(mode, label, arrs, args):
    dchi2 = arrs["dchi2"]
    dm2   = arrs["dm2"]
    ue4   = arrs["ue4"]
    umu4  = arrs["umu4"]
    utau4 = arrs["utau4"]

    if args.ymax is not None:
        ymax = args.ymax
    else:
        ymax = np.percentile(dchi2, args.y_percentile_max)

    ymin = 0.0
    ymask = (dchi2 >= ymin) & (dchi2 <= ymax)

    dchi2_p = dchi2[ymask]
    dm2_p   = dm2[ymask]
    ue4_p   = ue4[ymask]
    umu4_p  = umu4[ymask]
    utau4_p = utau4[ymask]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    axes = axes.flatten()

    panels = [
        (dm2_p,   r"Best-fit $\Delta m^2$ (eV$^2$)", "dm2"),
        (ue4_p,   r"Best-fit $|U_{e4}|^2$", "ue4"),
        (umu4_p,  r"Best-fit $|U_{\mu4}|^2$", "umu4"),
        (utau4_p, r"Best-fit $|U_{\tau4}|^2$", "utau4"),
    ]

    for ax, (x, xlabel, key) in zip(axes, panels):
        xr = safe_range(x)
        yr = (ymin, ymax)

        hist_kwargs = dict(
            bins=[args.xbins, args.ybins],
            range=[xr, yr],
            cmap="viridis",
        )

        if args.logz:
            hist_kwargs["norm"] = LogNorm()

        h = ax.hist2d(x, dchi2_p, **hist_kwargs)

        cb = fig.colorbar(h[3], ax=ax)
        cb.set_label("Toy count")

        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"Toy $\Delta\chi^2$")
        ax.set_title(f"{key} vs dchi2")

        if key == "dm2" and args.log_dm2:
            # only safe if positive
            if np.all(x > 0):
                ax.set_xscale("log")

        ax.grid(alpha=0.2)

    fig.suptitle(label, fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    os.makedirs(args.outdir, exist_ok=True)
    outpath = os.path.join(args.outdir, f"{args.outname_prefix}_{mode}.png")
    fig.savefig(outpath, dpi=200)
    plt.close(fig)

    print("Saved:", outpath)


def main():
    args = parse_args()

    if args.labels is not None and len(args.labels) != len(args.modes):
        raise ValueError("--labels must have same length as --modes")

    labels = args.labels if args.labels is not None else args.modes

    for mode, label in zip(args.modes, labels):
        arrs = load_mode(args.input_dir, mode)
        summarize(mode, label, arrs)
        make_one_plot(mode, label, arrs, args)


if __name__ == "__main__":
    main()
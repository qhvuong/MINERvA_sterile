#!/usr/bin/env python3

import os
import argparse
import glob
import numpy as np

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt


PARAMS = [
    ("dm2",  r"Best-fit $\Delta m^2$ [eV$^2$]"),
    ("ue4",  r"Best-fit $|U_{e4}|^2$"),
    ("umu4", r"Best-fit $|U_{\mu4}|^2$"),
    ("utau4", r"Best-fit $|U_{\tau4}|^2$"),
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot toy delta-chi2 versus fitted oscillation best-fit parameters."
    )

    parser.add_argument(
        "--input-dir",
        default="/exp/minerva/data/users/qvuong/surfaces/npys",
        help="Directory containing merged toy outputs.",
    )

    parser.add_argument(
        "--modes",
        nargs="+",
        required=True,
        help="Mode names, without prefixes/suffixes.",
    )

    parser.add_argument(
        "--labels",
        nargs="+",
        default=None,
        help="Labels for plot titles. Must match number of modes if provided.",
    )

    parser.add_argument(
        "--outdir",
        default="/exp/minerva/data/users/qvuong/surfaces/plots/fc_compare",
        help="Output plot directory.",
    )

    parser.add_argument(
        "--outname-prefix",
        default="FC_dchi2_vs_bestfit",
        help="Prefix for output filenames.",
    )

    parser.add_argument(
        "--alpha",
        type=float,
        default=0.35,
        help="Scatter marker transparency.",
    )

    parser.add_argument(
        "--marker-size",
        type=float,
        default=8.0,
        help="Scatter marker size.",
    )

    parser.add_argument(
        "--xmax-dchi2",
        type=float,
        default=None,
        help="Optional maximum delta-chi2 shown on y-axis.",
    )

    parser.add_argument(
        "--data-dchi2",
        type=float,
        default=None,
        help="Optional horizontal data delta-chi2 line.",
    )

    return parser.parse_args()


def load_csv_table(input_dir, mode):
    """
    Try to load a merged CSV table.

    Expected possible filenames:
      toy_bestfits_<mode>.csv
      asimov_bestfits_<mode>.csv
      sample_dchi2s_<mode>.csv
      <mode>.csv

    Expected columns include some recognizable version of:
      dchi2, dm2, ue4, umu4, utau4

    The script is permissive about column names.
    """
    candidates = [
        os.path.join(input_dir, "toy_bestfits_{}.csv".format(mode)),
        os.path.join(input_dir, "asimov_bestfits_{}.csv".format(mode)),
        os.path.join(input_dir, "sample_dchi2s_{}.csv".format(mode)),
        os.path.join(input_dir, "{}.csv".format(mode)),
    ]

    # Also allow any CSV with the mode name.
    candidates.extend(sorted(glob.glob(os.path.join(input_dir, "*{}*.csv".format(mode)))))

    for f in candidates:
        if not os.path.isfile(f):
            continue

        try:
            arr = np.genfromtxt(f, delimiter=",", names=True, dtype=None, encoding=None)
        except Exception:
            continue

        if arr.size == 0 or arr.dtype.names is None:
            continue

        names = list(arr.dtype.names)
        lower = {name.lower(): name for name in names}

        def find_col(possible):
            for p in possible:
                if p.lower() in lower:
                    return lower[p.lower()]
            return None

        col_dchi2 = find_col(["dchi2", "delta_chi2", "deltachi2", "toy_dchi2"])
        col_dm2   = find_col(["dm2", "m", "best_dm2", "bf_dm2", "delta_m2", "delta_m"])
        col_ue4   = find_col(["ue4", "best_ue4", "bf_ue4", "u_e4", "ue4sq"])
        col_umu4  = find_col(["umu4", "best_umu4", "bf_umu4", "u_mu4", "umu4sq"])
        col_utau4 = find_col(["utau4", "best_utau4", "bf_utau4", "u_tau4", "utau4sq"])

        needed = [col_dchi2, col_dm2, col_ue4, col_umu4, col_utau4]
        if any(x is None for x in needed):
            print("Found CSV but missing needed columns:", f)
            print("columns =", names)
            continue

        # Handle one-row structured arrays.
        arr = np.atleast_1d(arr)

        data = {
            "dchi2": np.asarray(arr[col_dchi2], dtype=float),
            "dm2":   np.asarray(arr[col_dm2], dtype=float),
            "ue4":   np.asarray(arr[col_ue4], dtype=float),
            "umu4":  np.asarray(arr[col_umu4], dtype=float),
            "utau4": np.asarray(arr[col_utau4], dtype=float),
        }

        print("Loaded CSV:", f)
        return clean_data(data)

    return None


def load_npy_arrays(input_dir, mode):
    """
    Fall back to separate .npy arrays.

    Expected:
      asimov_deltachi2s_<mode>.npy
      asimov_best_dm2s_<mode>.npy
      asimov_best_ue4s_<mode>.npy
      asimov_best_umu4s_<mode>.npy
      asimov_best_utau4s_<mode>.npy

    Also accepts:
      best_dm2s_<mode>.npy, best_ue4s_<mode>.npy, etc.
    """

    def find_one(patterns):
        for p in patterns:
            f = os.path.join(input_dir, p.format(mode))
            if os.path.isfile(f):
                return f
        return None

    files = {
        "dchi2": find_one([
            "asimov_deltachi2s_{}.npy",
            "deltachi2s_{}.npy",
            "dchi2s_{}.npy",
        ]),
        "dm2": find_one([
            "asimov_best_dm2s_{}.npy",
            "best_dm2s_{}.npy",
            "bf_dm2s_{}.npy",
        ]),
        "ue4": find_one([
            "asimov_best_ue4s_{}.npy",
            "best_ue4s_{}.npy",
            "bf_ue4s_{}.npy",
        ]),
        "umu4": find_one([
            "asimov_best_umu4s_{}.npy",
            "best_umu4s_{}.npy",
            "bf_umu4s_{}.npy",
        ]),
        "utau4": find_one([
            "asimov_best_utau4s_{}.npy",
            "best_utau4s_{}.npy",
            "bf_utau4s_{}.npy",
        ]),
    }

    missing = [k for k, f in files.items() if f is None]
    if len(missing) > 0:
        raise IOError(
            "Could not find CSV or required .npy arrays for mode {}. Missing: {}".format(
                mode, missing
            )
        )

    data = {}
    for key, f in files.items():
        data[key] = np.asarray(np.load(f), dtype=float).ravel()
        print("Loaded {}: {}".format(key, f))

    return clean_data(data)


def clean_data(data):
    n = min(len(v) for v in data.values())

    for k in data:
        data[k] = np.asarray(data[k], dtype=float).ravel()[:n]

    mask = np.ones(n, dtype=bool)
    for k in data:
        mask &= np.isfinite(data[k])

    for k in data:
        data[k] = data[k][mask]

    if len(data["dchi2"]) == 0:
        raise RuntimeError("No finite rows after cleaning")

    return data


def load_mode(input_dir, mode):
    data = load_csv_table(input_dir, mode)
    if data is not None:
        return data

    return load_npy_arrays(input_dir, mode)


def summarize(label, data):
    print("\n===== {} =====".format(label))
    print("ntoys =", len(data["dchi2"]))

    print("dchi2: mean={:.6g}, median={:.6g}, 95%={:.6g}".format(
        np.mean(data["dchi2"]),
        np.median(data["dchi2"]),
        np.percentile(data["dchi2"], 95),
    ))

    for p, _ in PARAMS:
        print("{}: min={:.6g}, median={:.6g}, max={:.6g}".format(
            p,
            np.min(data[p]),
            np.median(data[p]),
            np.max(data[p]),
        ))


def safe_tag(x):
    return str(x).replace(" ", "_").replace(",", "").replace("/", "-")


def make_one_mode_plot(mode, label, data, args):
    os.makedirs(args.outdir, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharey=True)
    axes = axes.ravel()

    dchi2 = data["dchi2"]

    for ax, (param, xlabel) in zip(axes, PARAMS):
        ax.scatter(
            data[param],
            dchi2,
            s=args.marker_size,
            alpha=args.alpha,
            linewidths=0,
        )

        ax.set_xlabel(xlabel)
        ax.grid(True, alpha=0.25)

        if args.data_dchi2 is not None:
            ax.axhline(
                args.data_dchi2,
                linestyle="--",
                linewidth=1.8,
                label="data Δχ²",
            )

        if args.xmax_dchi2 is not None:
            ax.set_ylim(0, args.xmax_dchi2)

    axes[0].set_ylabel(r"Toy $\Delta\chi^2$")
    axes[2].set_ylabel(r"Toy $\Delta\chi^2$")

    fig.suptitle(
        r"Toy $\Delta\chi^2$ vs fitted BF parameters: {}".format(label),
        fontsize=14,
    )

    if args.data_dchi2 is not None:
        axes[0].legend(framealpha=0.7)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    outpath = os.path.join(
        args.outdir,
        "{}_{}.png".format(args.outname_prefix, safe_tag(mode)),
    )

    fig.savefig(outpath, dpi=200)
    plt.close(fig)

    print("Saved:", outpath)


def make_overlay_param_plots(loaded, args):
    """
    One overlay plot per parameter, comparing modes.
    """
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for param, xlabel in PARAMS:
        fig, ax = plt.subplots(figsize=(9, 6))

        for i, (mode, label, data) in enumerate(loaded):
            color = colors[i % len(colors)]

            ax.scatter(
                data[param],
                data["dchi2"],
                s=args.marker_size,
                alpha=args.alpha,
                linewidths=0,
                label=label,
                color=color,
            )

        ax.set_xlabel(xlabel, fontsize=13)
        ax.set_ylabel(r"Toy $\Delta\chi^2$", fontsize=13)
        ax.set_title(r"Toy $\Delta\chi^2$ vs {}".format(xlabel), fontsize=14)
        ax.grid(True, alpha=0.25)

        if args.xmax_dchi2 is not None:
            ax.set_ylim(0, args.xmax_dchi2)

        ax.legend(framealpha=0.7, fontsize=10)

        fig.tight_layout()

        outpath = os.path.join(
            args.outdir,
            "{}_overlay_{}.png".format(args.outname_prefix, param),
        )

        fig.savefig(outpath, dpi=200)
        plt.close(fig)

        print("Saved:", outpath)


def main():
    args = parse_args()

    if args.labels is not None and len(args.labels) != len(args.modes):
        raise ValueError("--labels must have same length as --modes")

    labels = args.labels if args.labels is not None else args.modes

    loaded = []
    for mode, label in zip(args.modes, labels):
        data = load_mode(args.input_dir, mode)
        summarize(label, data)
        loaded.append((mode, label, data))

    for mode, label, data in loaded:
        make_one_mode_plot(mode, label, data, args)

    if len(loaded) > 1:
        make_overlay_param_plots(loaded, args)


if __name__ == "__main__":
    main()
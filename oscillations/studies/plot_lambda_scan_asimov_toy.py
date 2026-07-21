#!/usr/bin/env python3

import os
import csv
import math
import argparse
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def finite(value):
    try:
        return math.isfinite(float(value))
    except Exception:
        return False


def read_toy_rows(path, toy_index):
    rows = []

    with open(path, "r") as input_file:
        reader = csv.DictReader(input_file)

        for raw_row in reader:
            try:
                row_toy = int(float(raw_row["toy_index"]))
            except Exception:
                continue

            if row_toy != toy_index:
                continue

            row = {}

            for key, value in raw_row.items():
                try:
                    row[key] = float(value)
                except Exception:
                    row[key] = value

            required = [
                "nprof",
                "lambda",
                "resid_null",
                "penalty_null",
                "chi2_null",
                "norm_a_null",
                "max_abs_a_null",
            ]

            if not all(key in row for key in required):
                continue

            if not all(
                finite(row[key])
                for key in required
            ):
                continue

            row["toy_index"] = row_toy
            row["nprof"] = int(row["nprof"])
            row["lambda"] = float(row["lambda"])

            rows.append(row)

    return rows


def group_by_nprof(rows):
    grouped = defaultdict(list)

    for row in rows:
        grouped[int(row["nprof"])].append(row)

    for nprof in grouped:
        grouped[nprof].sort(
            key=lambda row: row["lambda"]
        )

    return grouped


def lambda_is_labeled(lam):
    targets = [
        0.03,
        0.1,
        0.3,
        0.5,
        0.7,
        0.85,
        1.0,
        1.15,
        1.5,
        2.0,
        3.0,
        10.0,
        30.0,
        100.0,
    ]

    return any(
        abs(
            math.log10(lam)
            - math.log10(target)
        ) < 1e-8
        for target in targets
    )


def plot_all_nprof_lcurves(
    grouped,
    toy_index,
    outpath,
):
    """
    One figure containing the L-curve for every Nprof.
    """
    plt.figure(figsize=(10, 7.5))

    for nprof in sorted(grouped):
        rows = grouped[nprof]

        residual = [
            float(row["resid_null"])
            for row in rows
        ]

        penalty = [
            float(row["penalty_null"])
            for row in rows
        ]

        lambdas = [
            float(row["lambda"])
            for row in rows
        ]

        plt.plot(
            residual,
            penalty,
            marker="o",
            linewidth=1.8,
            markersize=4,
            label="Nprof={}".format(nprof),
        )

        # Mark lambda=1 distinctly for every Nprof.
        lambda1_index = min(
            range(len(lambdas)),
            key=lambda index: abs(
                math.log10(lambdas[index])
            ),
        )

        plt.scatter(
            [residual[lambda1_index]],
            [penalty[lambda1_index]],
            marker="x",
            s=55,
        )

    plt.xlabel(r"Residual $\chi^2$")
    plt.ylabel(r"Flux penalty $\lambda a^T a$")
    plt.title(
        "Null pseudo-data L-curves, toy {}\n"
        r"crosses indicate $\lambda=1$".format(
            toy_index
        )
    )
    plt.grid(True, alpha=0.30)
    plt.legend(
        fontsize=8,
        ncol=2,
    )
    plt.tight_layout()
    plt.savefig(
        outpath,
        dpi=200,
    )
    plt.close()

    print("Wrote", outpath)


def plot_one_lcurve_per_nprof(
    grouped,
    toy_index,
    outdir,
):
    """
    Separate, fully labeled L-curve for every Nprof.
    """
    for nprof in sorted(grouped):
        rows = grouped[nprof]

        residual = [
            float(row["resid_null"])
            for row in rows
        ]

        penalty = [
            float(row["penalty_null"])
            for row in rows
        ]

        lambdas = [
            float(row["lambda"])
            for row in rows
        ]

        plt.figure(figsize=(7.5, 6))

        plt.plot(
            residual,
            penalty,
            marker="o",
            linewidth=2,
            markersize=5,
        )

        for x, y, lam in zip(
            residual,
            penalty,
            lambdas,
        ):
            if lambda_is_labeled(lam):
                plt.text(
                    x,
                    y,
                    "{:g}".format(lam),
                    fontsize=7,
                )

        lambda1_index = min(
            range(len(lambdas)),
            key=lambda index: abs(
                math.log10(lambdas[index])
            ),
        )

        plt.scatter(
            [residual[lambda1_index]],
            [penalty[lambda1_index]],
            marker="x",
            s=70,
            label=r"$\lambda=1$",
        )

        plt.xlabel(r"Residual $\chi^2$")
        plt.ylabel(
            r"Flux penalty $\lambda a^T a$"
        )
        plt.title(
            "Null pseudo-data L-curve\n"
            "toy {}, Nprof={}".format(
                toy_index,
                nprof,
            )
        )
        plt.grid(True, alpha=0.30)
        plt.legend(fontsize=9)
        plt.tight_layout()

        outpath = os.path.join(
            outdir,
            "toy{}_lcurve_Nprof{}.png".format(
                toy_index,
                nprof,
            ),
        )

        plt.savefig(
            outpath,
            dpi=200,
        )
        plt.close()

        print("Wrote", outpath)


def plot_metrics_vs_lambda(
    grouped,
    toy_index,
    outdir,
):
    metrics = [
        (
            "chi2_null",
            r"Total profiled null $\chi^2$",
        ),
        (
            "resid_null",
            r"Residual null $\chi^2$",
        ),
        (
            "penalty_null",
            r"Flux penalty $\lambda a^T a$",
        ),
        (
            "norm_a_null",
            r"$||a||$",
        ),
        (
            "max_abs_a_null",
            r"$\max |a_i|$",
        ),
    ]

    for metric, ylabel in metrics:
        plt.figure(figsize=(9, 6))

        for nprof in sorted(grouped):
            rows = grouped[nprof]

            lambdas = [
                float(row["lambda"])
                for row in rows
            ]

            values = [
                float(row[metric])
                for row in rows
            ]

            plt.plot(
                lambdas,
                values,
                marker="o",
                linewidth=1.8,
                markersize=4,
                label="Nprof={}".format(nprof),
            )

        plt.xscale("log")
        plt.axvline(
            1.0,
            linestyle="--",
            linewidth=1.2,
        )
        plt.xlabel(r"$\lambda$")
        plt.ylabel(ylabel)
        plt.title(
            "Null pseudo-data toy {}".format(
                toy_index
            )
        )
        plt.grid(True, alpha=0.30)
        plt.legend(
            fontsize=8,
            ncol=2,
        )
        plt.tight_layout()

        outpath = os.path.join(
            outdir,
            "toy{}_{}_vs_lambda.png".format(
                toy_index,
                metric,
            ),
        )

        plt.savefig(
            outpath,
            dpi=200,
        )
        plt.close()

        print("Wrote", outpath)


def write_lambda1_summary(
    grouped,
    toy_index,
    outpath,
):
    fieldnames = [
        "toy_index",
        "nprof",
        "lambda",
        "chi2_null",
        "resid_null",
        "penalty_null",
        "norm_a_null",
        "max_abs_a_null",
    ]

    output_rows = []

    for nprof in sorted(grouped):
        row = min(
            grouped[nprof],
            key=lambda candidate: abs(
                math.log10(
                    candidate["lambda"]
                )
            ),
        )

        output_rows.append({
            "toy_index": int(toy_index),
            "nprof": int(nprof),
            "lambda": float(row["lambda"]),
            "chi2_null": float(
                row["chi2_null"]
            ),
            "resid_null": float(
                row["resid_null"]
            ),
            "penalty_null": float(
                row["penalty_null"]
            ),
            "norm_a_null": float(
                row["norm_a_null"]
            ),
            "max_abs_a_null": float(
                row["max_abs_a_null"]
            ),
        })

    with open(outpath, "w") as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=fieldnames,
        )

        writer.writeheader()
        writer.writerows(output_rows)

    print("Wrote", outpath)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot the L-curves for one pseudo-data "
            "toy across all Nprof values."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="CSV from lambda_scan_asimov.py.",
    )

    parser.add_argument(
        "--toy",
        type=int,
        default=0,
        help="Toy index to plot. Default: 0.",
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory.",
    )

    args = parser.parse_args()

    if not os.path.exists(args.input):
        raise RuntimeError(
            "Missing input CSV: {}".format(
                args.input
            )
        )

    os.makedirs(
        args.outdir,
        exist_ok=True,
    )

    rows = read_toy_rows(
        args.input,
        args.toy,
    )

    if len(rows) == 0:
        raise RuntimeError(
            "No valid rows found for toy {} in {}".format(
                args.toy,
                args.input,
            )
        )

    grouped = group_by_nprof(rows)

    print("")
    print("===== SINGLE-TOY L-CURVE PLOT SETUP =====")
    print("input        =", args.input)
    print("toy          =", args.toy)
    print("rows         =", len(rows))
    print(
        "Nprof values =",
        sorted(grouped.keys()),
    )
    print(
        "lambda counts=",
        {
            nprof: len(grouped[nprof])
            for nprof in sorted(grouped)
        },
    )
    print("outdir       =", args.outdir)
    print("")

    plot_all_nprof_lcurves(
        grouped,
        args.toy,
        os.path.join(
            args.outdir,
            "toy{}_lcurves_compare_nprof.png".format(
                args.toy
            ),
        ),
    )

    plot_one_lcurve_per_nprof(
        grouped,
        args.toy,
        args.outdir,
    )

    plot_metrics_vs_lambda(
        grouped,
        args.toy,
        args.outdir,
    )

    write_lambda1_summary(
        grouped,
        args.toy,
        os.path.join(
            args.outdir,
            "toy{}_lambda1_nprof_summary.csv".format(
                args.toy
            ),
        ),
    )

    print("")
    print("Done. Output:", args.outdir)


if __name__ == "__main__":
    main()
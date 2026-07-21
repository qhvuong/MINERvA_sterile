#!/usr/bin/env python3

import os
import csv
import math
import argparse
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_LABELS = {
    "direct": r"$\nu+e$ + CC$\nu_e$ + CC$\nu_\mu$, all profiled",
    "excludeRatio": (
        r"$\nu+e$ + CC ratios + CC$\nu_\mu$, "
        r"ratios excluded from profiling"
    ),
    "includeRatio": (
        r"$\nu+e$ + CC ratios + CC$\nu_\mu$, "
        r"all profiled"
    ),
    "onlyRatio": r"Profile only CC ratios",
}


def finite(value):
    try:
        return math.isfinite(float(value))
    except Exception:
        return False


def parse_case_spec(spec):
    """
    Parse:

        case_key=/path/to/file.csv
    """
    if "=" not in spec:
        raise ValueError(
            "Case must be supplied as key=path, got: {}".format(
                spec
            )
        )

    key, path = spec.split("=", 1)

    key = key.strip()
    path = path.strip()

    if key == "":
        raise ValueError(
            "Empty case key in: {}".format(spec)
        )

    if path == "":
        raise ValueError(
            "Empty path in: {}".format(spec)
        )

    return key, path


def read_case_rows(
    path,
    toy_index,
):
    """
    Read one toy from one configuration CSV.

    Returns:
        grouped[nprof] = rows sorted by lambda
    """
    grouped = defaultdict(list)

    with open(path, "r") as input_file:
        reader = csv.DictReader(input_file)

        for raw_row in reader:
            try:
                row_toy = int(
                    float(raw_row["toy_index"])
                )
            except Exception:
                continue

            if row_toy != toy_index:
                continue

            try:
                nprof = int(
                    float(raw_row["nprof"])
                )

                lam = float(
                    raw_row["lambda"]
                )

                residual = float(
                    raw_row["resid_null"]
                )

                penalty = float(
                    raw_row["penalty_null"]
                )

                chi2 = float(
                    raw_row["chi2_null"]
                )

            except Exception:
                continue

            if not all(
                finite(value)
                for value in [
                    lam,
                    residual,
                    penalty,
                    chi2,
                ]
            ):
                continue

            grouped[nprof].append({
                "toy_index": row_toy,
                "nprof": nprof,
                "lambda": lam,
                "resid_null": residual,
                "penalty_null": penalty,
                "chi2_null": chi2,
            })

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


def validate_lambda_grid(
    datasets,
    nprof,
):
    """
    Warn if the cases do not contain identical lambda grids.
    """
    grids = {}

    for key, case in datasets.items():
        if nprof not in case["grouped"]:
            continue

        grids[key] = tuple(
            row["lambda"]
            for row in case["grouped"][nprof]
        )

    if len(grids) <= 1:
        return

    reference_key = next(iter(grids))
    reference_grid = grids[reference_key]

    for key, grid in grids.items():
        if grid != reference_grid:
            print(
                "WARNING: lambda grid differs for "
                "Nprof={} between '{}' and '{}'".format(
                    nprof,
                    reference_key,
                    key,
                )
            )


def plot_cases_for_nprof(
    datasets,
    toy_index,
    nprof,
    outdir,
):
    """
    One L-curve plot for a fixed toy and Nprof.
    Curves are the different sample configurations.
    """
    plt.figure(figsize=(9, 7))

    plotted_any = False

    for key, case in datasets.items():
        grouped = case["grouped"]

        if nprof not in grouped:
            print(
                "WARNING: case '{}' has no Nprof={}".format(
                    key,
                    nprof,
                )
            )
            continue

        rows = grouped[nprof]

        residual = [
            row["resid_null"]
            for row in rows
        ]

        penalty = [
            row["penalty_null"]
            for row in rows
        ]

        lambdas = [
            row["lambda"]
            for row in rows
        ]

        if len(rows) < 2:
            continue

        plotted_any = True

        line, = plt.plot(
            residual,
            penalty,
            marker="o",
            linewidth=2,
            markersize=4.5,
            label=case["label"],
        )

        line_color = line.get_color()

        lambda1_index = min(
            range(len(lambdas)),
            key=lambda index: abs(
                math.log10(
                    lambdas[index]
                )
            ),
        )

        plt.scatter(
            [residual[lambda1_index]],
            [penalty[lambda1_index]],
            marker="x",
            s=65,
            color=line_color,
        )

        # Label lambda points using the same line color.
        for x, y, lam in zip(
            residual,
            penalty,
            lambdas,
        ):
            if lambda_is_labeled(lam):
                plt.annotate(
                    "{:g}".format(lam),
                    xy=(x, y),
                    xytext=(3, 3),
                    textcoords="offset points",
                    fontsize=6.5,
                    color=line_color,
                )

    if not plotted_any:
        plt.close()
        return

    plt.xlabel(r"Residual $\chi^2$")
    plt.ylabel(
        r"Flux penalty $\lambda a^T a$"
    )
    plt.title(
        "Null pseudo-data L-curve comparison\n"
        "toy {}, Nprof={}; crosses indicate "
        r"$\lambda=1$".format(
            toy_index,
            nprof,
        )
    )
    plt.grid(True, alpha=0.30)
    plt.legend(
        fontsize=8,
    )
    plt.tight_layout()

    outpath = os.path.join(
        outdir,
        "toy{}_lcurve_compare_cases_Nprof{}.png".format(
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


def plot_lambda1_points_by_nprof(
    datasets,
    toy_index,
    outdir,
):
    """
    Display only the lambda=1 location for all configurations
    as Nprof changes.
    """
    plt.figure(figsize=(9, 7))

    plotted_any = False

    for key, case in datasets.items():
        nprof_values = []
        residual_values = []
        penalty_values = []

        for nprof in sorted(case["grouped"]):
            rows = case["grouped"][nprof]

            row = min(
                rows,
                key=lambda candidate: abs(
                    math.log10(
                        candidate["lambda"]
                    )
                ),
            )

            nprof_values.append(nprof)
            residual_values.append(
                row["resid_null"]
            )
            penalty_values.append(
                row["penalty_null"]
            )

        if len(nprof_values) == 0:
            continue

        plotted_any = True

        line, = plt.plot(
            residual_values,
            penalty_values,
            marker="o",
            linewidth=2,
            label=case["label"],
        )

        line_color = line.get_color()

        for nprof, x, y in zip(
            nprof_values,
            residual_values,
            penalty_values,
        ):
            plt.annotate(
                str(nprof),
                xy=(x, y),
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=7,
                color=line_color,
            )

    if not plotted_any:
        plt.close()
        return

    plt.xlabel(
        r"Residual $\chi^2$ at $\lambda=1$"
    )
    plt.ylabel(
        r"Flux penalty at $\lambda=1$"
    )
    plt.title(
        "Null pseudo-data toy {}\n".format(toy_index)
        + r"$\lambda=1$ location versus $N_{\rm prof}$"
    )
    plt.grid(True, alpha=0.30)
    plt.legend(fontsize=8)
    plt.tight_layout()

    outpath = os.path.join(
        outdir,
        "toy{}_lambda1_compare_cases_vs_nprof.png".format(
            toy_index
        ),
    )

    plt.savefig(
        outpath,
        dpi=200,
    )
    plt.close()

    print("Wrote", outpath)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare pseudo-data L-curves from multiple "
            "analysis configurations for one toy."
        )
    )

    parser.add_argument(
        "--case",
        action="append",
        required=True,
        help=(
            "Case specification key=/path/to/file.csv. "
            "Repeat once per configuration."
        ),
    )

    parser.add_argument(
        "--toy",
        type=int,
        default=0,
        help="Toy index to compare. Default: 0.",
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory.",
    )

    args = parser.parse_args()

    os.makedirs(
        args.outdir,
        exist_ok=True,
    )

    datasets = {}

    for spec in args.case:
        key, path = parse_case_spec(spec)

        if not os.path.exists(path):
            raise RuntimeError(
                "Missing case file '{}': {}".format(
                    key,
                    path,
                )
            )

        grouped = read_case_rows(
            path=path,
            toy_index=args.toy,
        )

        if len(grouped) == 0:
            raise RuntimeError(
                "No rows found for toy {} in case '{}'".format(
                    args.toy,
                    key,
                )
            )

        label = DEFAULT_LABELS.get(
            key,
            key,
        )

        datasets[key] = {
            "path": path,
            "label": label,
            "grouped": grouped,
        }

        print(
            "Loaded case='{}' toy={} "
            "Nprof values={} path={}".format(
                label,
                args.toy,
                sorted(grouped.keys()),
                path,
            )
        )

    all_nprof = sorted(
        set(
            nprof
            for case in datasets.values()
            for nprof in case["grouped"].keys()
        )
    )

    print("")
    print("===== MULTI-CASE TOY L-CURVE SETUP =====")
    print("toy          =", args.toy)
    print("cases        =", list(datasets.keys()))
    print("Nprof values =", all_nprof)
    print("outdir       =", args.outdir)
    print("")

    for nprof in all_nprof:
        validate_lambda_grid(
            datasets,
            nprof,
        )

        plot_cases_for_nprof(
            datasets=datasets,
            toy_index=args.toy,
            nprof=nprof,
            outdir=args.outdir,
        )

    plot_lambda1_points_by_nprof(
        datasets=datasets,
        toy_index=args.toy,
        outdir=args.outdir,
    )

    print("")
    print("Done. Output:", args.outdir)


if __name__ == "__main__":
    main()
#!/usr/bin/env python3

import os
import csv
import math
import argparse
from collections import defaultdict, Counter

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def finite(value):
    try:
        return math.isfinite(float(value))
    except Exception:
        return False


def read_csv(path):
    rows = []

    with open(path, "r") as input_file:
        reader = csv.DictReader(input_file)

        for raw_row in reader:
            row = {}

            for key, value in raw_row.items():
                try:
                    row[key] = float(value)
                except Exception:
                    row[key] = value

            required = [
                "toy_index",
                "nprof",
                "lambda",
                "resid_null",
                "penalty_null",
                "chi2_null",
            ]

            if not all(key in row for key in required):
                continue

            if not all(
                finite(row[key])
                for key in required
            ):
                continue

            row["toy_index"] = int(row["toy_index"])
            row["nprof"] = int(row["nprof"])
            row["lambda"] = float(row["lambda"])

            rows.append(row)

    return rows


def percentile_summary(values):
    values = np.asarray(
        [
            float(value)
            for value in values
            if finite(value)
        ],
        dtype=float,
    )

    if len(values) == 0:
        return {
            "median": np.nan,
            "p16": np.nan,
            "p84": np.nan,
            "mean": np.nan,
            "std": np.nan,
            "count": 0,
        }

    return {
        "median": float(np.median(values)),
        "p16": float(np.percentile(values, 16.0)),
        "p84": float(np.percentile(values, 84.0)),
        "mean": float(np.mean(values)),
        "std": float(np.std(values)),
        "count": int(len(values)),
    }


def group_rows(rows):
    """
    Return:

        by_nprof_toy[nprof][toy] = rows sorted by lambda
        by_nprof_lambda[nprof][lambda] = rows over toys
    """
    by_nprof_toy = defaultdict(
        lambda: defaultdict(list)
    )

    by_nprof_lambda = defaultdict(
        lambda: defaultdict(list)
    )

    for row in rows:
        nprof = int(row["nprof"])
        toy = int(row["toy_index"])
        lam = float(row["lambda"])

        by_nprof_toy[nprof][toy].append(row)
        by_nprof_lambda[nprof][lam].append(row)

    for nprof in by_nprof_toy:
        for toy in by_nprof_toy[nprof]:
            by_nprof_toy[nprof][toy].sort(
                key=lambda row: row["lambda"]
            )

    return by_nprof_toy, by_nprof_lambda


def safe_name(value):
    return (
        str(value)
        .replace(" ", "_")
        .replace(",", "")
        .replace("/", "-")
        .replace("\\", "")
        .replace("$", "")
        .replace("{", "")
        .replace("}", "")
    )


def get_lambda1_row(rows):
    return min(
        rows,
        key=lambda row: abs(
            math.log10(float(row["lambda"]))
        ),
    )


def median_curve_for_nprof(
    by_lambda,
    xkey,
    ykey,
):
    lambdas = sorted(by_lambda.keys())

    output = []

    for lam in lambdas:
        point_rows = by_lambda[lam]

        x_summary = percentile_summary(
            [
                row[xkey]
                for row in point_rows
            ]
        )

        y_summary = percentile_summary(
            [
                row[ykey]
                for row in point_rows
            ]
        )

        output.append({
            "lambda": float(lam),

            "x_median": x_summary["median"],
            "x_p16": x_summary["p16"],
            "x_p84": x_summary["p84"],

            "y_median": y_summary["median"],
            "y_p16": y_summary["p16"],
            "y_p84": y_summary["p84"],

            "count": min(
                x_summary["count"],
                y_summary["count"],
            ),
        })

    return output


def plot_lcurve_by_nprof(
    by_nprof_toy,
    by_nprof_lambda,
    outdir,
):
    """
    One plot per Nprof:

      - individual toy L-curves in the background;
      - median L-curve in the foreground;
      - selected lambda labels on the median.
    """
    label_lambdas = {
        0.03,
        0.1,
        0.3,
        0.5,
        1.0,
        2.0,
        10.0,
        100.0,
    }

    for nprof in sorted(by_nprof_toy):
        plt.figure(figsize=(8, 6))

        # Individual toy curves.
        for toy, toy_rows in sorted(
            by_nprof_toy[nprof].items()
        ):
            residual = []
            penalty = []

            for row in toy_rows:
                if (
                    finite(row.get("resid_null"))
                    and finite(row.get("penalty_null"))
                ):
                    residual.append(
                        float(row["resid_null"])
                    )
                    penalty.append(
                        float(row["penalty_null"])
                    )

            if len(residual) < 2:
                continue

            plt.plot(
                residual,
                penalty,
                linewidth=0.6,
                alpha=0.10,
            )

        median_curve = median_curve_for_nprof(
            by_nprof_lambda[nprof],
            xkey="resid_null",
            ykey="penalty_null",
        )

        x = [
            point["x_median"]
            for point in median_curve
        ]

        y = [
            point["y_median"]
            for point in median_curve
        ]

        plt.plot(
            x,
            y,
            marker="o",
            linewidth=2.5,
            markersize=5,
            label="Toy median",
        )

        for point in median_curve:
            lam = point["lambda"]

            if any(
                abs(
                    math.log10(lam)
                    - math.log10(target)
                ) < 1e-8
                for target in label_lambdas
            ):
                plt.text(
                    point["x_median"],
                    point["y_median"],
                    "{:g}".format(lam),
                    fontsize=7,
                )

        plt.xlabel(r"Residual $\chi^2$")
        plt.ylabel(r"Flux penalty $\lambda a^T a$")
        plt.title(
            "Null pseudo-data L-curves, "
            r"$N_{\rm prof}={}$".format(nprof)
        )
        plt.grid(True, alpha=0.30)
        plt.legend(fontsize=9)
        plt.tight_layout()

        outpath = os.path.join(
            outdir,
            "toy_lcurve_Nprof{}.png".format(nprof),
        )

        plt.savefig(
            outpath,
            dpi=200,
        )
        plt.close()

        print("Wrote", outpath)


def plot_median_lcurves_all_nprof(
    by_nprof_lambda,
    outdir,
):
    """
    Compare median L-curves for every Nprof.
    """
    plt.figure(figsize=(9, 7))

    for nprof in sorted(by_nprof_lambda):
        curve = median_curve_for_nprof(
            by_nprof_lambda[nprof],
            xkey="resid_null",
            ykey="penalty_null",
        )

        x = [
            point["x_median"]
            for point in curve
        ]

        y = [
            point["y_median"]
            for point in curve
        ]

        plt.plot(
            x,
            y,
            marker="o",
            linewidth=1.8,
            markersize=3.5,
            label="Nprof={}".format(nprof),
        )

        # Mark lambda=1.
        lambda1_point = min(
            curve,
            key=lambda point: abs(
                math.log10(point["lambda"])
            ),
        )

        plt.scatter(
            [lambda1_point["x_median"]],
            [lambda1_point["y_median"]],
            marker="x",
            s=45,
        )

    plt.xlabel(r"Median residual $\chi^2$")
    plt.ylabel(r"Median flux penalty $\lambda a^T a$")
    plt.title(
        "Median null pseudo-data L-curves\n"
        r"cross marks indicate $\lambda=1$"
    )
    plt.grid(True, alpha=0.30)
    plt.legend(fontsize=8, ncol=2)
    plt.tight_layout()

    outpath = os.path.join(
        outdir,
        "median_lcurves_compare_nprof.png",
    )

    plt.savefig(
        outpath,
        dpi=200,
    )
    plt.close()

    print("Wrote", outpath)


def plot_metric_vs_lambda(
    by_nprof_lambda,
    metric,
    ylabel,
    outdir,
):
    """
    For every Nprof, plot the toy median metric versus lambda
    with a 16–84% band.
    """
    plt.figure(figsize=(9, 6))

    for nprof in sorted(by_nprof_lambda):
        lambdas = sorted(
            by_nprof_lambda[nprof].keys()
        )

        medians = []
        lower = []
        upper = []

        for lam in lambdas:
            summary = percentile_summary(
                [
                    row[metric]
                    for row
                    in by_nprof_lambda[nprof][lam]
                ]
            )

            medians.append(summary["median"])
            lower.append(summary["p16"])
            upper.append(summary["p84"])

        lambdas = np.asarray(
            lambdas,
            dtype=float,
        )

        medians = np.asarray(
            medians,
            dtype=float,
        )

        lower = np.asarray(
            lower,
            dtype=float,
        )

        upper = np.asarray(
            upper,
            dtype=float,
        )

        plt.plot(
            lambdas,
            medians,
            marker="o",
            linewidth=1.8,
            markersize=3.5,
            label="Nprof={}".format(nprof),
        )

        plt.fill_between(
            lambdas,
            lower,
            upper,
            alpha=0.08,
        )

    plt.xscale("log")
    plt.xlabel(r"$\lambda$")
    plt.ylabel(ylabel)
    plt.title(
        "Null pseudo-data median and 16–84% toy band"
    )
    plt.grid(True, alpha=0.30)
    plt.legend(fontsize=8, ncol=2)
    plt.tight_layout()

    outpath = os.path.join(
        outdir,
        "{}_vs_lambda_all_nprof.png".format(
            metric
        ),
    )

    plt.savefig(
        outpath,
        dpi=200,
    )
    plt.close()

    print("Wrote", outpath)


def plot_lambda1_metric_vs_nprof(
    by_nprof_toy,
    metric,
    ylabel,
    outdir,
):
    """
    At lambda closest to 1, summarize each metric over toys
    as a function of Nprof.
    """
    nprof_values = []
    medians = []
    lower = []
    upper = []

    for nprof in sorted(by_nprof_toy):
        values = []

        for toy, toy_rows in by_nprof_toy[nprof].items():
            row = get_lambda1_row(toy_rows)

            if finite(row.get(metric)):
                values.append(
                    float(row[metric])
                )

        summary = percentile_summary(values)

        nprof_values.append(nprof)
        medians.append(summary["median"])
        lower.append(summary["p16"])
        upper.append(summary["p84"])

    x = np.asarray(
        nprof_values,
        dtype=float,
    )

    medians = np.asarray(
        medians,
        dtype=float,
    )

    lower = np.asarray(
        lower,
        dtype=float,
    )

    upper = np.asarray(
        upper,
        dtype=float,
    )

    plt.figure(figsize=(8, 5.5))

    plt.plot(
        x,
        medians,
        marker="o",
        linewidth=2,
        label="Toy median",
    )

    plt.fill_between(
        x,
        lower,
        upper,
        alpha=0.20,
        label="16–84%",
    )

    plt.xlabel(
        "Number of profiled flux universes"
    )
    plt.ylabel(ylabel)
    plt.title(
        r"Null pseudo-data at $\lambda=1$"
    )
    plt.grid(True, alpha=0.30)
    plt.legend(fontsize=9)
    plt.tight_layout()

    outpath = os.path.join(
        outdir,
        "lambda1_{}_vs_nprof.png".format(
            metric
        ),
    )

    plt.savefig(
        outpath,
        dpi=200,
    )
    plt.close()

    print("Wrote", outpath)


def estimate_bend_lambda(rows):
    """
    Estimate the L-curve bend using discrete curvature in:

        x = log10(residual chi2)
        y = log10(penalty)

    parameterized by:

        t = log10(lambda)

    Endpoints are not candidates.

    Returns:
        bend_lambda, curvature

    This is a diagnostic estimator, not a formal statistical
    definition of the optimal lambda.
    """
    valid_rows = []

    for row in rows:
        residual = row.get("resid_null")
        penalty = row.get("penalty_null")
        lam = row.get("lambda")

        if not (
            finite(residual)
            and finite(penalty)
            and finite(lam)
        ):
            continue

        residual = float(residual)
        penalty = float(penalty)
        lam = float(lam)

        if (
            residual <= 0.0
            or penalty <= 0.0
            or lam <= 0.0
        ):
            continue

        valid_rows.append({
            "lambda": lam,
            "t": math.log10(lam),
            "x": math.log10(residual),
            "y": math.log10(penalty),
        })

    valid_rows.sort(
        key=lambda point: point["t"]
    )

    if len(valid_rows) < 3:
        return np.nan, np.nan

    t = np.asarray(
        [
            point["t"]
            for point in valid_rows
        ],
        dtype=float,
    )

    x = np.asarray(
        [
            point["x"]
            for point in valid_rows
        ],
        dtype=float,
    )

    y = np.asarray(
        [
            point["y"]
            for point in valid_rows
        ],
        dtype=float,
    )

    dx_dt = np.gradient(x, t)
    dy_dt = np.gradient(y, t)

    d2x_dt2 = np.gradient(dx_dt, t)
    d2y_dt2 = np.gradient(dy_dt, t)

    denominator = (
        dx_dt * dx_dt
        + dy_dt * dy_dt
    ) ** 1.5

    numerator = np.abs(
        dx_dt * d2y_dt2
        - dy_dt * d2x_dt2
    )

    curvature = np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator),
        where=denominator > 1e-14,
    )

    # Avoid selecting a boundary point.
    candidate_indices = np.arange(
        1,
        len(valid_rows) - 1,
    )

    if len(candidate_indices) == 0:
        return np.nan, np.nan

    best_index = candidate_indices[
        np.argmax(
            curvature[candidate_indices]
        )
    ]

    return (
        float(valid_rows[best_index]["lambda"]),
        float(curvature[best_index]),
    )


def calculate_bend_results(by_nprof_toy):
    results = []

    for nprof in sorted(by_nprof_toy):
        for toy in sorted(by_nprof_toy[nprof]):
            bend_lambda, curvature = estimate_bend_lambda(
                by_nprof_toy[nprof][toy]
            )

            results.append({
                "toy_index": int(toy),
                "nprof": int(nprof),
                "bend_lambda": float(bend_lambda),
                "bend_log10_lambda": (
                    float(math.log10(bend_lambda))
                    if finite(bend_lambda)
                    and bend_lambda > 0.0
                    else np.nan
                ),
                "bend_curvature": float(curvature),
            })

    return results


def write_bend_results(
    bend_results,
    outpath,
):
    fieldnames = [
        "toy_index",
        "nprof",
        "bend_lambda",
        "bend_log10_lambda",
        "bend_curvature",
    ]

    with open(outpath, "w") as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=fieldnames,
        )

        writer.writeheader()
        writer.writerows(bend_results)

    print("Wrote", outpath)


def plot_bend_lambda_vs_nprof(
    bend_results,
    outdir,
):
    grouped = defaultdict(list)

    for row in bend_results:
        if (
            finite(row["bend_lambda"])
            and row["bend_lambda"] > 0.0
        ):
            grouped[int(row["nprof"])].append(
                float(row["bend_lambda"])
            )

    nprof_values = []
    medians = []
    lower = []
    upper = []

    for nprof in sorted(grouped):
        # Summarize in log(lambda) space.
        log_values = np.log10(
            np.asarray(
                grouped[nprof],
                dtype=float,
            )
        )

        summary = percentile_summary(
            log_values
        )

        nprof_values.append(nprof)
        medians.append(
            10.0 ** summary["median"]
        )
        lower.append(
            10.0 ** summary["p16"]
        )
        upper.append(
            10.0 ** summary["p84"]
        )

    x = np.asarray(
        nprof_values,
        dtype=float,
    )

    median = np.asarray(
        medians,
        dtype=float,
    )

    lower = np.asarray(
        lower,
        dtype=float,
    )

    upper = np.asarray(
        upper,
        dtype=float,
    )

    plt.figure(figsize=(8, 5.5))

    plt.plot(
        x,
        median,
        marker="o",
        linewidth=2,
        label="Median bend",
    )

    plt.fill_between(
        x,
        lower,
        upper,
        alpha=0.20,
        label="16–84%",
    )

    plt.axhline(
        1.0,
        linestyle="--",
        linewidth=1.5,
        label=r"$\lambda=1$",
    )

    plt.yscale("log")
    plt.xlabel(
        "Number of profiled flux universes"
    )
    plt.ylabel(
        r"Estimated bend $\lambda$"
    )
    plt.title(
        "L-curve bend location across null toys"
    )
    plt.grid(True, alpha=0.30)
    plt.legend(fontsize=9)
    plt.tight_layout()

    outpath = os.path.join(
        outdir,
        "bend_lambda_vs_nprof.png",
    )

    plt.savefig(
        outpath,
        dpi=200,
    )
    plt.close()

    print("Wrote", outpath)


def select_nprof_per_toy(
    bend_results,
):
    """
    For each toy, select the Nprof whose estimated bend is
    closest to lambda=1 in log10(lambda) space.
    """
    by_toy = defaultdict(list)

    for row in bend_results:
        if not finite(
            row["bend_log10_lambda"]
        ):
            continue

        by_toy[int(row["toy_index"])].append(
            row
        )

    selections = []

    for toy, toy_rows in sorted(by_toy.items()):
        selected = min(
            toy_rows,
            key=lambda row: (
                abs(row["bend_log10_lambda"]),
                int(row["nprof"]),
            ),
        )

        selections.append({
            "toy_index": int(toy),
            "selected_nprof": int(
                selected["nprof"]
            ),
            "bend_lambda": float(
                selected["bend_lambda"]
            ),
            "distance_from_lambda1": float(
                abs(
                    selected[
                        "bend_log10_lambda"
                    ]
                )
            ),
        })

    return selections


def write_nprof_selections(
    selections,
    outpath,
):
    fieldnames = [
        "toy_index",
        "selected_nprof",
        "bend_lambda",
        "distance_from_lambda1",
    ]

    with open(outpath, "w") as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=fieldnames,
        )

        writer.writeheader()
        writer.writerows(selections)

    print("Wrote", outpath)


def plot_selected_nprof(
    selections,
    outdir,
):
    counts = Counter(
        int(row["selected_nprof"])
        for row in selections
    )

    nprof_values = sorted(
        counts.keys()
    )

    values = [
        counts[nprof]
        for nprof in nprof_values
    ]

    plt.figure(figsize=(8, 5.5))

    plt.bar(
        [
            str(nprof)
            for nprof in nprof_values
        ],
        values,
    )

    plt.xlabel(
        "Selected number of profiled flux universes"
    )
    plt.ylabel("Number of toys")
    plt.title(
        r"$N_{\rm prof}$ whose L-curve bend is closest to $\lambda=1$"
    )
    plt.grid(
        True,
        axis="y",
        alpha=0.30,
    )
    plt.tight_layout()

    outpath = os.path.join(
        outdir,
        "selected_nprof_distribution.png",
    )

    plt.savefig(
        outpath,
        dpi=200,
    )
    plt.close()

    print("Wrote", outpath)


def write_lambda1_summary(
    by_nprof_toy,
    outpath,
):
    metrics = [
        "chi2_null",
        "resid_null",
        "penalty_null",
        "norm_a_null",
        "max_abs_a_null",
    ]

    fieldnames = [
        "nprof",
        "ntoys",
    ]

    for metric in metrics:
        fieldnames.extend([
            "{}_median".format(metric),
            "{}_p16".format(metric),
            "{}_p84".format(metric),
            "{}_mean".format(metric),
            "{}_std".format(metric),
        ])

    rows = []

    for nprof in sorted(by_nprof_toy):
        row = {
            "nprof": int(nprof),
            "ntoys": int(
                len(by_nprof_toy[nprof])
            ),
        }

        for metric in metrics:
            values = []

            for toy_rows in (
                by_nprof_toy[nprof].values()
            ):
                lambda1_row = get_lambda1_row(
                    toy_rows
                )

                if finite(
                    lambda1_row.get(metric)
                ):
                    values.append(
                        lambda1_row[metric]
                    )

            summary = percentile_summary(
                values
            )

            row[
                "{}_median".format(metric)
            ] = summary["median"]

            row[
                "{}_p16".format(metric)
            ] = summary["p16"]

            row[
                "{}_p84".format(metric)
            ] = summary["p84"]

            row[
                "{}_mean".format(metric)
            ] = summary["mean"]

            row[
                "{}_std".format(metric)
            ] = summary["std"]

        rows.append(row)

    with open(outpath, "w") as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=fieldnames,
        )

        writer.writeheader()
        writer.writerows(rows)

    print("Wrote", outpath)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot null pseudo-data lambda and Nprof "
            "L-curve diagnostics."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="CSV produced by lambda_scan_asimov.py.",
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Directory for plots and summary CSV files.",
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

    rows = read_csv(
        args.input
    )

    if len(rows) == 0:
        raise RuntimeError(
            "No valid rows found in {}".format(
                args.input
            )
        )

    (
        by_nprof_toy,
        by_nprof_lambda,
    ) = group_rows(rows)

    nprof_values = sorted(
        by_nprof_toy.keys()
    )

    toy_values = sorted(
        set(
            row["toy_index"]
            for row in rows
        )
    )

    lambda_values = sorted(
        set(
            row["lambda"]
            for row in rows
        )
    )

    print("")
    print("===== PSEUDO-DATA L-CURVE PLOT SETUP =====")
    print("input         =", args.input)
    print("rows          =", len(rows))
    print("toys          =", len(toy_values))
    print("toy range     =", (
        min(toy_values),
        max(toy_values),
    ))
    print("Nprof values  =", nprof_values)
    print("lambda values =", lambda_values)
    print("outdir        =", args.outdir)
    print("")

    plot_lcurve_by_nprof(
        by_nprof_toy,
        by_nprof_lambda,
        args.outdir,
    )

    plot_median_lcurves_all_nprof(
        by_nprof_lambda,
        args.outdir,
    )

    for metric, ylabel in [
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
    ]:
        plot_metric_vs_lambda(
            by_nprof_lambda,
            metric,
            ylabel,
            args.outdir,
        )

        plot_lambda1_metric_vs_nprof(
            by_nprof_toy,
            metric,
            ylabel,
            args.outdir,
        )

    bend_results = calculate_bend_results(
        by_nprof_toy
    )

    write_bend_results(
        bend_results,
        os.path.join(
            args.outdir,
            "toy_nprof_bend_summary.csv",
        ),
    )

    plot_bend_lambda_vs_nprof(
        bend_results,
        args.outdir,
    )

    selections = select_nprof_per_toy(
        bend_results
    )

    write_nprof_selections(
        selections,
        os.path.join(
            args.outdir,
            "selected_nprof_by_toy.csv",
        ),
    )

    plot_selected_nprof(
        selections,
        args.outdir,
    )

    write_lambda1_summary(
        by_nprof_toy,
        os.path.join(
            args.outdir,
            "lambda1_nprof_toy_summary.csv",
        ),
    )

    print("")
    print("Done. Output:", args.outdir)


if __name__ == "__main__":
    main()
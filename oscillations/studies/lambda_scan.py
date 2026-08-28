import os
import sys
import argparse
import csv
import shutil
import numpy as np

def parse_scan_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hist-config-tag", default=None)
    parser.add_argument("--exclude", default=None)
    parser.add_argument("--mask", default="none")
    parser.add_argument("--out", default=None)
    parser.add_argument(
        "--lambdas",
        # default="0.01,0.03,0.1,0.3,1,3,10,30,100",
        default="0.01,0.03,0.1,0.2,0.3,0.5,0.7,0.85,1,1.15,1.3,1.5,2,3,5,10,30,100",
        help="Comma-separated lambda values"
    )
    parser.add_argument(
        "--skip-bf",
        action="store_true",
        help="Only compute null profiled chi2. Useful for data L-curve checks.",
    )
    parser.add_argument(
        "--profile-only",
        default=None,
        help="Use only selected samples in flux profiling solve. Options: none, ratio.",
    )
    parser.add_argument(
        "--profile-n-universes",
        default=None,
        type=int,
        help="Use only the first N flux universes as profiling nuisance directions.",
    )
    args, remaining = parser.parse_known_args()

    # Hide lambda-scan-only args from AnalysisConfig/project parsers.
    sys.argv = [sys.argv[0]] + remaining

    return args

args = parse_scan_args()

import ROOT
import PlotUtils

from config.AnalysisConfig import AnalysisConfig
from tools.StitchedHistogram import StitchedHistogram
from tools.Fitters import Statistics, OscillationFitter

ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

def write_rows_csv(rows, out_csv):
    if len(rows) == 0:
        return

    fieldnames = list(rows[0].keys())

    with open(out_csv, "w") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

def parse_mask(mask_name, hist_config):
    """
    Predefined masks for stress tests.
    Bin counts are determined from the selected HIST_CONFIG.
    """
    if mask_name in [None, "", "none", "None"]:
        return None

    import json

    with open(hist_config, "r") as f:
        cfg = json.load(f)

    def all_local_bins(sample):
        if sample not in cfg:
            return []

        n_bins = cfg[sample]["end"] - cfg[sample]["start"] + 1
        return list(range(1, n_bins + 1))

    if mask_name == "nonratio":
        return {
            "fhc_elastic": all_local_bins("fhc_elastic"),
            "fhc_numu_selection": all_local_bins("fhc_numu_selection"),
            "rhc_numu_selection": all_local_bins("rhc_numu_selection"),
        }

    if mask_name == "ratio":
        return {
            "fhc_ratio": all_local_bins("fhc_ratio"),
            "rhc_ratio": all_local_bins("rhc_ratio"),
        }

    if mask_name == "bin1":
        return {
            "fhc_ratio": [1],
            "rhc_ratio": [1],
        }

    if mask_name == "bins9_10":
        return {
            "fhc_ratio": [9, 10],
            "rhc_ratio": [9, 10],
        }

    if mask_name == "bins1_9_10":
        return {
            "fhc_ratio": [1, 9, 10],
            "rhc_ratio": [1, 9, 10],
        }

    if mask_name == "fhc_bins1_9_10":
        return {
            "fhc_ratio": [1, 9, 10],
        }

    if mask_name == "rhc_bins1_9_10":
        return {
            "rhc_ratio": [1, 9, 10],
        }

    raise ValueError("Unknown mask_name: {}".format(mask_name))


def get_residual_and_penalty(stat, useOsc=False):
    """
    Calls Chi2DataMC and returns total chi2, residual chi2, penalty, norm(a), max|a|.
    """
    chi2, penalty = stat.Chi2DataMC(
        marginalize=True,
        useOsc=useOsc
    )

    flux_fitter = stat.GetFluxFitter(useOsc=useOsc)
    a = flux_fitter.GetFluxSolution()

    residual = chi2 - penalty
    norm_a = np.linalg.norm(a)
    max_abs_a = np.max(np.abs(a))

    return chi2, residual, penalty, norm_a, max_abs_a


def run_one_lambda(
    sample_histogram,
    lam,
    exclude_samples,
    mask_spec,
    hist_config,
    skip_bf=False,
    profile_only=None,
    profile_n_universes=None,
):
    """
    Run null profiled chi2, and optionally oscillation best fit, for one lambda value.
    """

    # Null at this lambda
    stat_null = Statistics(
        sample_histogram,
        exclude=exclude_samples,
        lam=lam,
        mask_spec=mask_spec,
        profile_only=profile_only,
        profile_n_universes=profile_n_universes,
        hist_config=hist_config,
    )

    chi2_null, resid_null, pen_null, norm_null, max_null = get_residual_and_penalty(
        stat_null,
        useOsc=False
    )

    row = {
        "lambda": lam,
        "log_lambda": np.log10(lam),
        "nprof": profile_n_universes if profile_n_universes is not None else -1,

        "chi2_null": chi2_null,
        "resid_null": resid_null,
        "penalty_null": pen_null,
        "norm_a_null": norm_null,
        "max_abs_a_null": max_null,
    }

    if skip_bf:
        row.update({
            "chi2_bf": np.nan,
            "resid_bf": np.nan,
            "penalty_bf": np.nan,
            "norm_a_bf": np.nan,
            "max_abs_a_bf": np.nan,
            "delta_chi2": np.nan,
            "dm2": np.nan,
            "ue4": np.nan,
            "umu4": np.nan,
            "utau4": np.nan,
        })
        return row

    # Fit BF at this lambda
    fitter = OscillationFitter(
        sample_histogram,
        exclude=exclude_samples,
        lam=lam,
        marginalize_flux=True,
        mask_spec=mask_spec,
        profile_only=profile_only,
        profile_n_universes=profile_n_universes,
        hist_config=hist_config,
    )

    chi2_bf_fit, res = fitter.DoFit()

    # Recompute BF components cleanly at final point
    sample_histogram.OscillateHistogram(
        res["m"],
        res["ue4"],
        res["umu4"],
        res["utau4"]
    )

    stat_bf = Statistics(
        sample_histogram,
        exclude=exclude_samples,
        lam=lam,
        mask_spec=mask_spec,
        profile_only=profile_only,
        profile_n_universes=profile_n_universes,
        hist_config=hist_config,
    )

    chi2_bf, resid_bf, pen_bf, norm_bf, max_bf = get_residual_and_penalty(
        stat_bf,
        useOsc=True
    )

    row.update({
        "chi2_bf": chi2_bf,
        "resid_bf": resid_bf,
        "penalty_bf": pen_bf,
        "norm_a_bf": norm_bf,
        "max_abs_a_bf": max_bf,

        "delta_chi2": chi2_null - chi2_bf,

        "dm2": res["m"],
        "ue4": res["ue4"],
        "umu4": res["umu4"],
        "utau4": res["utau4"],
    })

    return row


def main():
    global args

    ccnueroot = os.environ.get("CCNUEROOT")
    if ccnueroot is None:
        raise RuntimeError("CCNUEROOT is not set")

    plot_tag = args.hist_config_tag
    if plot_tag is None:
        plot_tag = getattr(AnalysisConfig, "hist_config_tag", "default")

    if plot_tag in [None, "", "none"]:
        plot_tag = "default"

    exclude_samples = args.exclude
    if exclude_samples is None:
        exclude_samples = getattr(AnalysisConfig, "exclude", "")

    if exclude_samples is None:
        exclude_samples = ""

    if isinstance(exclude_samples, str):
        if exclude_samples.strip().lower() in ["none", ""]:
            exclude_samples = ""

    filename = "rootfiles/NuE_stitched_hists_{}.root".format(plot_tag)
    file_path = "{}/oscillations/{}".format(ccnueroot, filename)

    hist_config = "HIST_CONFIG_{}.json".format(plot_tag)
    if not os.path.exists(hist_config):
        raise RuntimeError("Missing requested hist config file: {}".format(hist_config))

    mask_spec = parse_mask(args.mask, hist_config)
    
    shutil.copyfile(hist_config, "HIST_CONFIG.json")

    lambda_values = [float(x) for x in args.lambdas.split(",")]

    if args.out is None:
        mask_label = args.mask if args.mask not in [None, ""] else "none"
        exclude_label = exclude_samples if exclude_samples not in [None, ""] else "none"
        out_csv = "lambda_scan_{}_exclude_{}_mask_{}.csv".format(
            plot_tag,
            exclude_label.replace(",", "_"),
            mask_label
        )
    else:
        out_csv = args.out

    out_dir = os.path.dirname(out_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    print("\n===== lambda scan setup =====")
    print("plot_tag       =", plot_tag)
    print("file           =", file_path)
    print("hist_config    =", hist_config)
    print("exclude        =", exclude_samples)
    print("profile_only   =", args.profile_only)
    print("profile_n_universes =", args.profile_n_universes)
    print("mask           =", args.mask)
    print("mask_spec      =", mask_spec)
    print("lambdas        =", lambda_values)
    print("output csv     =", out_csv)
    print("skip_bf        =", args.skip_bf)

    rows = []

    for lam in lambda_values:
        print("\n\n==============================================")
        print("Running lambda =", lam)
        print("==============================================")

        # Reload fresh histogram for each lambda to avoid stale oscillated state.
        sample_histogram = StitchedHistogram("sample")
        sample_histogram.SetHistConfig(hist_config)
        sample_histogram.Load(file_path)

        row = run_one_lambda(
            sample_histogram=sample_histogram,
            lam=lam,
            exclude_samples=exclude_samples,
            mask_spec=mask_spec,
            hist_config=hist_config,
            skip_bf=args.skip_bf,
            profile_only=args.profile_only,
            profile_n_universes=args.profile_n_universes,
        )

        rows.append(row)

        print("\n===== lambda result =====")
        for k in [
            "lambda",
            "log_lambda",
            "nprof",
            "chi2_null",
            "chi2_bf",
            "delta_chi2",
            "resid_null",
            "resid_bf",
            "penalty_null",
            "penalty_bf",
            "norm_a_null",
            "norm_a_bf",
            "max_abs_a_null",
            "max_abs_a_bf",
            "dm2",
            "ue4",
            "umu4",
            "utau4",
        ]:
            print("{:20s} {}".format(k, row[k]))

        write_rows_csv(rows, out_csv)
        print("Wrote partial results to", out_csv)

    write_rows_csv(rows, out_csv)

    print("\n===== final lambda scan table =====")
    for row in rows:
        print(row)
    print("\nSaved:", out_csv)


if __name__ == "__main__":
    main()
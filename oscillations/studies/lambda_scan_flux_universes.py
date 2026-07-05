#!/usr/bin/env python3

import os
import sys
import argparse
import csv
import shutil
import numpy as np
import logging


def parse_scan_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--hist-config-tag", default=None)
    parser.add_argument("--exclude", default=None)
    parser.add_argument("--mask", default="none")
    parser.add_argument("--out", default=None)

    parser.add_argument(
        "--universe-indices",
        default="0,1,2,3,4",
        help='Comma-separated flux universe indices to throw, or "all".',
    )

    parser.add_argument(
        "--lambdas",
        default="0.01,0.03,0.1,0.2,0.3,0.5,0.7,0.85,1,1.15,1.3,1.5,2,3,5,10,30,100",
        # default="0.01,0.1,0.5,0.85,1,1.15,1.3,1.5,2,5,10,30,100",
        help="Comma-separated lambda values",
    )

    parser.add_argument(
        "--poisson",
        action="store_true",
        help="Also Poisson fluctuate the selected flux-universe prediction. Default is pure flux-universe Asimov.",
    )

    parser.add_argument("--seed", type=int, default=12345)

    parser.add_argument(
        "--skip-bf",
        action="store_true",
        help="Only compute null profiled chi2. Useful for quick L-curve checks.",
    )

    parser.add_argument(
        "--profile-only",
        default=None,
        help="Use only selected samples in flux profiling solve. Options: none, ratio.",
    )

    parser.add_argument(
        "--profile-n-universes",
        type=int,
        default=None,
        help="Use only the first N flux universes as profiling basis. Default: use all.",
    )

    parser.add_argument(
        "--nuniverse-throws",
        type=int,
        default=None,
        help="Use thrown flux universe indices 0...(N-1). Overrides --universe-indices.",
    )

    args, remaining = parser.parse_known_args()

    # Hide scan-only args from AnalysisConfig/project parsers.
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


def parse_mask(mask_name):
    if mask_name in [None, "", "none", "None"]:
        return None

    if mask_name == "nonratio":
        return {
            "fhc_elastic": list(range(1, 5)),
            "fhc_numu_selection": list(range(1, 13)),
            "rhc_numu_selection": list(range(1, 13)),
        }

    if mask_name == "ratio":
        return {
            "fhc_ratio": list(range(1, 13)),
            "rhc_ratio": list(range(1, 13)),
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


def parse_universe_indices(spec, n_flux_universes):
    if spec is None or spec == "":
        return list(range(min(5, n_flux_universes)))

    if str(spec).strip().lower() == "all":
        return list(range(n_flux_universes))

    indices = [int(x.strip()) for x in str(spec).split(",") if x.strip() != ""]

    bad = [i for i in indices if i < 0 or i >= n_flux_universes]
    if len(bad) > 0:
        raise RuntimeError(
            "Requested universe indices {} outside valid range 0..{}".format(
                bad,
                n_flux_universes - 1,
            )
        )

    return indices


def overwrite_data_hist_from_stitched_vector(histogram, stitched_vector):
    """
    Replace the stitched data histogram with a pseudo-data vector.

    This follows the same strategy as fitAsimovs_quinn.py:
    start from the existing stitched data histogram, build weights so that
    DivideSingle changes the bin contents to the toy values, then call
    SetDataHistogram().
    """
    stitched_data = histogram.GetDataHistogram()

    toy = np.asarray(stitched_vector, dtype=float)
    if toy.ndim != 1:
        raise RuntimeError("stitched_vector must be 1D")

    if stitched_data.GetNbinsX() != len(toy):
        raise RuntimeError(
            "Toy vector length {} does not match stitched data bins {}".format(
                len(toy),
                stitched_data.GetNbinsX(),
            )
        )

    weights = stitched_data.Clone().GetCVHistoWithStatError()

    for i in range(1, weights.GetNbinsX() + 1):
        if toy[i - 1] != 0:
            weight = stitched_data.GetBinContent(i) / toy[i - 1]
        else:
            weight = stitched_data.GetBinContent(i)

        weights.SetBinContent(i, weight)
        weights.SetBinError(i, 0)

    data_histogram = stitched_data.Clone()
    data_histogram.DivideSingle(data_histogram, weights)

    histogram.SetDataHistogram(data_histogram)


def make_flux_universe_pseudodata(histogram, universe_index, poisson=False, rng=None):
    """
    Build pseudo-data from one stored Flux universe.

    A has shape:
        n_flux_universes x n_bins

    and by construction:
        A[k] = flux_universe_k - CV

    so:
        pseudo = CV + A[k]
    """
    cv = np.array(histogram.GetMCHistogram())[1:-1].astype(float)
    A = histogram.GetAMatrix()

    if universe_index < 0 or universe_index >= A.shape[0]:
        raise RuntimeError(
            "universe_index {} outside A matrix range 0..{}".format(
                universe_index,
                A.shape[0] - 1,
            )
        )

    pseudo_mean = cv + A[universe_index]

    # Stored flux universes should be physical, but protect Poisson from tiny negatives.
    min_before_clip = np.min(pseudo_mean)
    pseudo_mean = np.clip(pseudo_mean, 0.0, None)

    if poisson:
        if rng is None:
            raise RuntimeError("poisson=True requires rng")
        pseudo = rng.poisson(pseudo_mean).astype(float)
    else:
        pseudo = pseudo_mean.astype(float)

    return pseudo, min_before_clip


def get_residual_and_penalty(stat, useOsc=False):
    chi2, penalty = stat.Chi2DataMC(
        marginalize=True,
        useOsc=useOsc,
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
    skip_bf=False,
    profile_only=None,
    profile_n_universes=None,
):
    """
    Run null profiled chi2 and optionally oscillation best fit for one lambda value.
    """

    stat_null = Statistics(
        sample_histogram,
        exclude=exclude_samples,
        lam=lam,
        mask_spec=mask_spec,
        profile_only=profile_only,
        profile_n_universes=profile_n_universes,
    )

    chi2_null, resid_null, pen_null, norm_null, max_null = get_residual_and_penalty(
        stat_null,
        useOsc=False,
    )

    row = {
        "lambda": lam,
        "log_lambda": np.log10(lam),

        "chi2_null": chi2_null,
        "resid_null": resid_null,
        "penalty_null": pen_null,
        "norm_a_null": norm_null,
        "max_abs_a_null": max_null,
    }

    if skip_bf:
        row.update(
            {
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
            }
        )
        return row

    fitter = OscillationFitter(
        sample_histogram,
        exclude=exclude_samples,
        lam=lam,
        marginalize_flux=True,
        mask_spec=mask_spec,
        profile_only=profile_only,
        profile_n_universes=profile_n_universes,
    )

    chi2_bf_fit, res = fitter.DoFit()

    # Recompute BF components cleanly at final point.
    sample_histogram.OscillateHistogram(
        res["m"],
        res["ue4"],
        res["umu4"],
        res["utau4"],
    )

    stat_bf = Statistics(
        sample_histogram,
        exclude=exclude_samples,
        lam=lam,
        mask_spec=mask_spec,
        profile_only=profile_only,
        profile_n_universes=profile_n_universes,
    )

    chi2_bf, resid_bf, pen_bf, norm_bf, max_bf = get_residual_and_penalty(
        stat_bf,
        useOsc=True,
    )

    row.update(
        {
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
        }
    )

    return row


def main():
    global args

    np.random.seed(args.seed)
    rng = np.random.default_rng(args.seed)

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

    mask_spec = parse_mask(args.mask)

    filename = "rootfiles/NuE_stitched_hists_{}.root".format(plot_tag)
    file_path = "{}/oscillations/{}".format(ccnueroot, filename)

    hist_config = "HIST_CONFIG_{}.json".format(plot_tag)
    if not os.path.exists(hist_config):
        raise RuntimeError("Missing requested hist config file: {}".format(hist_config))

    shutil.copyfile(hist_config, "HIST_CONFIG.json")

    lambda_values = [float(x) for x in args.lambdas.split(",")]

    # Load once to inspect flux-universe count.
    source_hist = StitchedHistogram("source")
    source_hist.Load(file_path)

    A = source_hist.GetAMatrix()
    n_flux_universes = A.shape[0]
    n_bins = A.shape[1]

    if args.nuniverse_throws is not None:
        if args.nuniverse_throws <= 0:
            raise RuntimeError("--nuniverse-throws must be positive")

        if args.nuniverse_throws > n_flux_universes:
            raise RuntimeError(
                "--nuniverse-throws={} requested, but only {} flux universes exist".format(
                    args.nuniverse_throws,
                    n_flux_universes,
                )
            )

        universe_indices = list(range(args.nuniverse_throws))
    else:
        universe_indices = parse_universe_indices(
            args.universe_indices,
            n_flux_universes,
        )

    profile_n_universes = args.profile_n_universes

    if profile_n_universes is not None:
        if profile_n_universes <= 0:
            raise RuntimeError("--profile-n-universes must be positive")

        if profile_n_universes > n_flux_universes:
            raise RuntimeError(
                "--profile-n-universes={} requested, but only {} flux universes exist".format(
                    profile_n_universes,
                    n_flux_universes,
                )
            )

    if args.out is None:
        mask_label = args.mask if args.mask not in [None, ""] else "none"
        exclude_label = exclude_samples if exclude_samples not in [None, ""] else "none"
        univ_label = (
            "all"
            if len(universe_indices) == n_flux_universes
            else "u" + "-".join([str(u) for u in universe_indices])
        )
        pois_label = "_poisson" if args.poisson else ""
        bf_label = "_nullOnly" if args.skip_bf else ""
        profile_label = (
            args.profile_only
            if args.profile_only not in [None, "", "none", "None"]
            else "none"
        )
        nprof_label = (
            "allProfUniv"
            if profile_n_universes is None
            else "Nprof{}".format(profile_n_universes)
        )
        out_csv = "lambda_scan_fluxUniverseThrows_{}_exclude_{}_profileOnly_{}_mask_{}_{}_{}{}{}.csv".format(
            plot_tag,
            exclude_label.replace(",", "_"),
            profile_label,
            mask_label,
            nprof_label,
            univ_label,
            pois_label,
            bf_label,
        )
    else:
        out_csv = args.out

    out_dir = os.path.dirname(out_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    print("\n===== flux-universe lambda scan setup =====")
    print("plot_tag             =", plot_tag)
    print("file                 =", file_path)
    print("hist_config          =", hist_config)
    print("exclude              =", exclude_samples)
    print("profile_only         =", args.profile_only)
    print("profile_n_universes =", profile_n_universes)
    print("mask                 =", args.mask)
    print("mask_spec            =", mask_spec)
    print("A shape              =", A.shape)
    print("n_flux_universes     =", n_flux_universes)
    print("n_bins               =", n_bins)
    print("universe_indices     =", universe_indices)
    print("poisson              =", args.poisson)
    print("skip_bf              =", args.skip_bf)
    print("seed                 =", args.seed)
    print("lambdas              =", lambda_values)
    print("output csv           =", out_csv)

    rows = []

    for thrown_universe in universe_indices:
        print("\n\n################################################")
        print("Building pseudo-data from flux universe =", thrown_universe)
        print("################################################")

        # Build pseudo-data from the source histogram.
        pseudo_vector, min_before_clip = make_flux_universe_pseudodata(
            source_hist,
            universe_index=thrown_universe,
            poisson=args.poisson,
            rng=rng,
        )

        cv = np.array(source_hist.GetMCHistogram())[1:-1].astype(float)
        print("CV sum                =", np.sum(cv))
        print("pseudo sum            =", np.sum(pseudo_vector))
        print("pseudo min            =", np.min(pseudo_vector))
        print("pseudo max            =", np.max(pseudo_vector))
        print("min before clip       =", min_before_clip)
        print("sum pseudo-CV         =", np.sum(pseudo_vector - cv))
        print("norm pseudo-CV        =", np.linalg.norm(pseudo_vector - cv))
        print("max abs pseudo-CV     =", np.max(np.abs(pseudo_vector - cv)))

        for lam in lambda_values:
            print("\n\n==============================================")
            print("Thrown universe =", thrown_universe, "lambda =", lam)
            print("==============================================")

            # Reload fresh histogram for each lambda to avoid stale oscillated state.
            sample_histogram = StitchedHistogram("sample")
            sample_histogram.Load(file_path)

            overwrite_data_hist_from_stitched_vector(
                sample_histogram,
                pseudo_vector,
            )

            check_data = np.array(sample_histogram.GetDataHistogram())[1:-1]
            print("check data sum after SetDataHistogram =", np.sum(check_data))
            print("check max |data - pseudo| =", np.max(np.abs(check_data - pseudo_vector)))

            row = run_one_lambda(
                sample_histogram=sample_histogram,
                lam=lam,
                exclude_samples=exclude_samples,
                mask_spec=mask_spec,
                skip_bf=args.skip_bf,
                profile_only=args.profile_only,
                profile_n_universes=profile_n_universes,
            )

            row["thrown_universe"] = thrown_universe
            row["n_flux_universes_total"] = n_flux_universes
            row["profile_n_universes"] = (
                profile_n_universes if profile_n_universes is not None else n_flux_universes
            )
            row["profile_n_universes_label"] = (
                "all" if profile_n_universes is None else str(profile_n_universes)
            )
            row["pseudo_type"] = "stored_flux_universe"
            row["poisson"] = int(args.poisson)
            row["seed"] = args.seed
            row["hist_config_tag"] = plot_tag
            row["exclude"] = exclude_samples if exclude_samples != "" else "none"
            row["mask"] = args.mask
            row["profile_only"] = args.profile_only if args.profile_only not in [None, "", "none", "None"] else "none"

            rows.append(row)

            print("\n===== flux-universe lambda result =====")
            for k in [
                "thrown_universe",
                "lambda",
                "log_lambda",
                "profile_n_universes",
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
                print("{:25s} {}".format(k, row[k]))

            write_rows_csv(rows, out_csv)
            print("Wrote partial results to", out_csv)

    write_rows_csv(rows, out_csv)

    print("\n===== final flux-universe lambda scan table =====")
    for row in rows:
        print(row)

    print("\nSaved:", out_csv)


if __name__ == "__main__":
    main()
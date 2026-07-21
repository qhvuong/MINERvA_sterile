#!/usr/bin/env python3

import os
import sys
import csv
import argparse
import shutil
import numpy as np
import gc
import json

def parse_scan_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run null pseudo-data L-curve scans over lambda and Nprof. "
            "One pseudo-dataset is reused for every lambda and Nprof "
            "within each toy."
        )
    )

    parser.add_argument("--hist-config-tag", default=None)

    parser.add_argument(
        "--toy-source-hist-config-tag",
        default="prodNueel_noRatio_p8",
        help=(
            "Direct CCnue stitched configuration used to generate the "
            "primitive elastic, CCnue, CCnuebar, CCnumu, and CCnumubar "
            "pseudo-data. All analysis configurations should use the "
            "same source tag and seed."
        ),
    )

    parser.add_argument("--exclude", default=None)
    parser.add_argument("--mask", default="none")
    parser.add_argument("--profile-only", default=None)

    parser.add_argument(
        "--lambdas",
        default="0.3,0.5,0.7,0.85,1,1.15,1.3,1.5,2,3,10",
        help="Comma-separated lambda values.",
    )

    parser.add_argument(
        "--nprof-values",
        default="30,40,50,75,100,200,500",
        help="Comma-separated Nprof values.",
    )

    parser.add_argument(
        "--ntoys",
        type=int,
        default=100,
        help="Number of null pseudo-data toys.",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=12345,
        help="Base NumPy random seed.",
    )

    parser.add_argument(
        "--throw-mode",
        choices=[
            "poisson",
            "gaussian_sansflux",
            "asimov",
        ],
        default="poisson",
        help=(
            "Pseudo-data construction. "
            "'poisson' throws each nominal count bin independently; "
            "'gaussian_sansflux' draws a correlated shift using the "
            "sans-flux covariance; "
            "'asimov' uses the nominal MC without fluctuation."
        ),
    )

    parser.add_argument(
        "--out",
        default=None,
        help="Output CSV path.",
    )

    args, remaining = parser.parse_known_args()

    # Hide scan-only arguments from AnalysisConfig/project parsers.
    sys.argv = [sys.argv[0]] + remaining

    return args


args = parse_scan_args()


import ROOT
import PlotUtils

from config.AnalysisConfig import AnalysisConfig
from tools.StitchedHistogram import StitchedHistogram
from tools.Fitters import Statistics

ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)


def parse_csv_floats(spec):
    values = [
        float(value.strip())
        for value in str(spec).split(",")
        if value.strip() != ""
    ]

    if len(values) == 0:
        raise RuntimeError("No floating-point values were provided")

    return values


def parse_csv_ints(spec):
    values = [
        int(value.strip())
        for value in str(spec).split(",")
        if value.strip() != ""
    ]

    if len(values) == 0:
        raise RuntimeError("No integer values were provided")

    return values


def write_rows_csv(rows, out_csv):
    if len(rows) == 0:
        return

    fieldnames = list(rows[0].keys())

    with open(out_csv, "w") as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(rows)


def normalize_exclude(exclude):
    if exclude is None:
        return ""

    if isinstance(exclude, str):
        if exclude.strip().lower() in ["", "none"]:
            return ""

    return exclude


def normalize_profile_only(profile_only):
    if profile_only is None:
        return None

    if isinstance(profile_only, str):
        if profile_only.strip().lower() in ["", "none"]:
            return None

    return profile_only


def parse_mask(mask_name):
    """
    Predefined masks used in the original lambda-scan studies.
    """
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

    raise ValueError(
        "Unknown mask name: {}".format(mask_name)
    )


def covariance_square_root(covariance):
    """
    Return a PSD square root of a symmetric covariance matrix.

    Small negative eigenvalues from numerical precision are clipped to zero.
    """
    covariance = np.asarray(
        covariance,
        dtype=float,
    )

    covariance = 0.5 * (
        covariance + covariance.T
    )

    eigenvalues, eigenvectors = np.linalg.eigh(
        covariance
    )

    minimum_before_clip = float(
        np.min(eigenvalues)
    )

    eigenvalues = np.clip(
        eigenvalues,
        0.0,
        None,
    )

    sqrt_covariance = (
        eigenvectors
        @ np.diag(np.sqrt(eigenvalues))
        @ eigenvectors.T
    )

    return sqrt_covariance, minimum_before_clip


DIRECT_SLICES = {
    "fhc_elastic": slice(0, 4),
    "fhc_nue_selection": slice(4, 16),
    "rhc_nue_selection": slice(16, 28),
    "fhc_numu_selection": slice(28, 40),
    "rhc_numu_selection": slice(40, 52),
}


def make_primitive_null_pseudodata(
    source_histogram,
    throw_mode,
    rng,
    sqrt_sansflux=None,
):
    """
    Generate one toy in the direct primitive sample space:

        fhc_elastic
        fhc_nue_selection
        rhc_nue_selection
        fhc_numu_selection
        rhc_numu_selection

    The same primitive toy can then be transformed into either the
    direct or ratio stitched configuration.
    """
    nominal = np.asarray(
        source_histogram.GetMCHistogram(),
        dtype=float,
    )[1:-1]

    if len(nominal) != 52:
        raise RuntimeError(
            "Primitive direct source must contain 52 bins, got {}".format(
                len(nominal)
            )
        )

    if not np.all(np.isfinite(nominal)):
        raise RuntimeError(
            "Primitive nominal MC contains non-finite values"
        )

    if throw_mode == "asimov":
        pseudo = nominal.copy()

    elif throw_mode == "poisson":
        if np.any(nominal < 0.0):
            bad = np.where(nominal < 0.0)[0]
            raise RuntimeError(
                "Poisson throw requested but primitive nominal bins "
                "are negative: {}".format(bad.tolist())
            )

        pseudo = rng.poisson(
            nominal
        ).astype(float)

    elif throw_mode == "gaussian_sansflux":
        if sqrt_sansflux is None:
            raise RuntimeError(
                "gaussian_sansflux requires the direct-source "
                "sans-flux covariance square root"
            )

        z = rng.normal(
            loc=0.0,
            scale=1.0,
            size=len(nominal),
        )

        pseudo_raw = nominal + sqrt_sansflux @ z

        negative_bins = np.where(pseudo_raw < 0.0)[0]

        if len(negative_bins) > 0:
            print(
                "WARNING: Gaussian primitive toy has {} negative bins "
                "before clipping: global bins {}".format(
                    len(negative_bins),
                    (negative_bins + 1).tolist(),
                )
            )

        pseudo = np.clip(
            pseudo_raw,
            0.0,
            None,
        )

    else:
        raise RuntimeError(
            "Unknown throw mode: {}".format(
                throw_mode
            )
        )

    return nominal, pseudo


def safe_ratio(numerator, denominator, label):
    """
    Construct a finite ratio and reject toys containing zero
    denominator bins.
    """
    numerator = np.asarray(
        numerator,
        dtype=float,
    )

    denominator = np.asarray(
        denominator,
        dtype=float,
    )

    zero_denominator = np.abs(denominator) <= 1e-12

    if np.any(zero_denominator):
        raise RuntimeError(
            "{} has zero denominator in local bins {}".format(
                label,
                (np.where(zero_denominator)[0] + 1).tolist(),
            )
        )

    ratio = np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator, dtype=float),
        where=~zero_denominator,
    )

    if not np.all(np.isfinite(ratio)):
        raise RuntimeError(
            "{} ratio contains non-finite values".format(label)
        )

    return ratio


def convert_primitive_vector_to_target(
    primitive_vector,
    target_has_ratio,
):
    """
    Convert the common primitive direct-space vector to the stitched
    layout used by the target analysis configuration.

    Direct target:
        elastic, nue, nuebar, numu, numubar

    Ratio target:
        elastic, numu/nue, numubar/nuebar, numu, numubar
    """
    primitive_vector = np.asarray(
        primitive_vector,
        dtype=float,
    )

    if len(primitive_vector) != 52:
        raise RuntimeError(
            "Primitive vector must contain 52 bins, got {}".format(
                len(primitive_vector)
            )
        )

    fhc_elastic = primitive_vector[
        DIRECT_SLICES["fhc_elastic"]
    ]

    fhc_nue = primitive_vector[
        DIRECT_SLICES["fhc_nue_selection"]
    ]

    rhc_nue = primitive_vector[
        DIRECT_SLICES["rhc_nue_selection"]
    ]

    fhc_numu = primitive_vector[
        DIRECT_SLICES["fhc_numu_selection"]
    ]

    rhc_numu = primitive_vector[
        DIRECT_SLICES["rhc_numu_selection"]
    ]

    if not target_has_ratio:
        return np.concatenate([
            fhc_elastic,
            fhc_nue,
            rhc_nue,
            fhc_numu,
            rhc_numu,
        ])

    fhc_ratio = safe_ratio(
        fhc_numu,
        fhc_nue,
        "FHC CCnumu/CCnue",
    )

    rhc_ratio = safe_ratio(
        rhc_numu,
        rhc_nue,
        "RHC CCnumubar/CCnuebar",
    )

    return np.concatenate([
        fhc_elastic,
        fhc_ratio,
        rhc_ratio,
        fhc_numu,
        rhc_numu,
    ])


def set_pseudo_histogram_from_vector(
    histogram,
    pseudo_vector,
    toy_index,
):
    """
    Set the StitchedHistogram pseudo histogram from a stitched-bin vector.

    The pseudo histogram is cloned from the nominal MC so it has the
    correct binning and ROOT histogram type.
    """
    pseudo_vector = np.asarray(
        pseudo_vector,
        dtype=float,
    )

    nominal_histogram = histogram.GetMCHistogram()

    if nominal_histogram.GetNbinsX() != len(pseudo_vector):
        raise RuntimeError(
            "Pseudo vector has {} bins, but stitched histogram has {}".format(
                len(pseudo_vector),
                nominal_histogram.GetNbinsX(),
            )
        )

    pseudo_histogram = nominal_histogram.Clone(
        "pseudo_null_toy_{}".format(toy_index)
    )

    pseudo_histogram.SetDirectory(0)

    for bin_number in range(
        1,
        pseudo_histogram.GetNbinsX() + 1,
    ):
        value = float(
            pseudo_vector[bin_number - 1]
        )

        pseudo_histogram.SetBinContent(
            bin_number,
            value,
        )

        # The χ² covariance is handled separately.
        # These ROOT bin errors are not used by Statistics.Chi2DataMC.
        pseudo_histogram.SetBinError(
            bin_number,
            np.sqrt(max(value, 0.0)),
        )

    histogram.SetPseudoHistogram(
        pseudo_histogram
    )


def get_residual_and_penalty(statistic):
    """
    Evaluate the profiled null χ² against the pseudo-data histogram.
    """
    chi2, penalty = statistic.Chi2DataMC(
        marginalize=True,
        usePseudo=True,
        useOsc=False,
    )

    flux_fitter = statistic.GetFluxFitter(
        useOsc=False
    )

    solution = np.array(
        flux_fitter.GetFluxSolution(),
        dtype=float,
        copy=True,
    )

    residual = float(
        chi2 - penalty
    )

    norm_a = float(
        np.linalg.norm(solution)
    )

    max_abs_a = float(
        np.max(np.abs(solution))
    )

    del flux_fitter

    return {
        "chi2_null": float(chi2),
        "resid_null": residual,
        "penalty_null": float(penalty),
        "norm_a_null": norm_a,
        "max_abs_a_null": max_abs_a,
    }


def run_one_point(
    sample_histogram,
    lam,
    nprof,
    exclude_samples,
    mask_spec,
    profile_only,
):
    statistic = Statistics(
        sample_histogram,
        exclude=exclude_samples,
        lam=lam,
        mask_spec=mask_spec,
        profile_only=profile_only,
        profile_n_universes=nprof,
    )

    result = get_residual_and_penalty(
        statistic
    )

    del statistic

    row = {
        "lambda": float(lam),
        "log_lambda": float(np.log10(lam)),
        "nprof": int(nprof),
    }

    row.update(result)

    return row


def main():
    global args

    ccnueroot = os.environ.get(
        "CCNUEROOT"
    )

    if ccnueroot is None:
        raise RuntimeError(
            "CCNUEROOT is not set"
        )

    plot_tag = args.hist_config_tag

    if plot_tag is None:
        plot_tag = getattr(
            AnalysisConfig,
            "hist_config_tag",
            "default",
        )

    if plot_tag in [None, "", "none"]:
        plot_tag = "default"

    exclude_samples = args.exclude

    if exclude_samples is None:
        exclude_samples = getattr(
            AnalysisConfig,
            "exclude",
            "",
        )

    exclude_samples = normalize_exclude(
        exclude_samples
    )

    profile_only = normalize_profile_only(
        args.profile_only
    )

    mask_spec = parse_mask(
        args.mask
    )

    lambda_values = parse_csv_floats(
        args.lambdas
    )

    nprof_values = parse_csv_ints(
        args.nprof_values
    )

    if args.ntoys <= 0:
        raise RuntimeError(
            "--ntoys must be positive"
        )

    if any(value <= 0 for value in nprof_values):
        raise RuntimeError(
            "All Nprof values must be positive"
        )

    # --------------------------------------------------
    # Target configuration: the analysis being scanned.
    # --------------------------------------------------
    target_filename = (
        "rootfiles/"
        "NuE_stitched_hists_{}.root"
    ).format(plot_tag)

    target_file_path = os.path.join(
        ccnueroot,
        "oscillations",
        target_filename,
    )

    target_hist_config = (
        "HIST_CONFIG_{}.json"
    ).format(plot_tag)

    # --------------------------------------------------
    # Primitive source configuration: always direct CCnue.
    # --------------------------------------------------
    toy_source_tag = args.toy_source_hist_config_tag

    source_filename = (
        "rootfiles/"
        "NuE_stitched_hists_{}.root"
    ).format(toy_source_tag)

    source_file_path = os.path.join(
        ccnueroot,
        "oscillations",
        source_filename,
    )

    source_hist_config = (
        "HIST_CONFIG_{}.json"
    ).format(toy_source_tag)

    for label, path in [
        ("target stitched ROOT file", target_file_path),
        ("target HIST_CONFIG", target_hist_config),
        ("toy-source stitched ROOT file", source_file_path),
        ("toy-source HIST_CONFIG", source_hist_config),
    ]:
        if not os.path.exists(path):
            raise RuntimeError(
                "Missing {}: {}".format(label, path)
            )

    # Load the direct primitive source using its own configuration.
    shutil.copyfile(
        source_hist_config,
        "HIST_CONFIG.json",
    )

    primitive_source_histogram = StitchedHistogram(
        "primitive_source"
    )

    primitive_source_histogram.Load(
        source_file_path
    )

    expected_primitive_samples = {
        "fhc_elastic",
        "fhc_nue_selection",
        "rhc_nue_selection",
        "fhc_numu_selection",
        "rhc_numu_selection",
    }

    missing_primitive_samples = (
        expected_primitive_samples
        - set(primitive_source_histogram.keys)
    )

    if missing_primitive_samples:
        raise RuntimeError(
            "Toy source must be the direct configuration. "
            "Missing primitive samples: {}".format(
                sorted(missing_primitive_samples)
            )
        )

    # Now activate and load the target analysis configuration.
    shutil.copyfile(
        target_hist_config,
        "HIST_CONFIG.json",
    )

    target_source_histogram = StitchedHistogram(
        "target_source"
    )

    target_source_histogram.Load(
        target_file_path
    )

    with open(target_hist_config, "r") as config_file:
        target_layout = json.load(config_file)

    has_ratio_samples = (
        "fhc_ratio" in target_layout
        and "rhc_ratio" in target_layout
    )

    has_direct_nue_samples = (
        "fhc_nue_selection" in target_layout
        and "rhc_nue_selection" in target_layout
    )

    if has_ratio_samples == has_direct_nue_samples:
        raise RuntimeError(
            "Could not uniquely identify target stitched layout. "
            "Expected either ratio samples or direct CCnue samples. "
            "HIST_CONFIG keys: {}".format(
                list(target_layout.keys())
            )
        )

    total_flux_universes = (
        target_source_histogram.GetAMatrix().shape[0]
    )

    invalid_nprof = [
        value
        for value in nprof_values
        if value > total_flux_universes
    ]

    if len(invalid_nprof) > 0:
        raise RuntimeError(
            "Requested Nprof values {} exceed the {} available "
            "flux universes".format(
                invalid_nprof,
                total_flux_universes,
            )
        )

    sqrt_sansflux = None
    minimum_sansflux_eigenvalue = np.nan

    if args.throw_mode == "gaussian_sansflux":
        primitive_sansflux_covariance = (
            primitive_source_histogram.GetCovarianceMatrix(
                sansFlux=True
            )
        )

        if primitive_sansflux_covariance.shape != (52, 52):
            raise RuntimeError(
                "Expected primitive direct sans-flux covariance "
                "shape (52, 52), got {}".format(
                    primitive_sansflux_covariance.shape
                )
            )

        (
            sqrt_sansflux,
            minimum_sansflux_eigenvalue,
        ) = covariance_square_root(
            primitive_sansflux_covariance
        )


    if args.out is None:
        exclude_label = (
            exclude_samples
            if exclude_samples != ""
            else "none"
        )

        profile_label = (
            profile_only
            if profile_only is not None
            else "none"
        )

        out_csv = (
            "lambda_scan_asimov_"
            "{}_exclude_{}_profileOnly_{}_mask_{}_"
            "{}_Ntoys{}.csv"
        ).format(
            plot_tag,
            exclude_label.replace(",", "_"),
            profile_label,
            args.mask,
            args.throw_mode,
            args.ntoys,
        )

    else:
        out_csv = args.out

    out_dir = os.path.dirname(
        out_csv
    )

    if out_dir:
        os.makedirs(
            out_dir,
            exist_ok=True,
        )

    print("")
    print("===== NULL PSEUDO-DATA NPROF L-CURVE SCAN =====")
    print("target plot_tag          =", plot_tag)
    print("target file              =", target_file_path)
    print("target hist_config       =", target_hist_config)
    print("toy source tag           =", toy_source_tag)
    print("toy source file          =", source_file_path)
    print("toy source hist_config   =", source_hist_config)
    print(
        "target samples           =",
        target_source_histogram.keys,
    )
    print(
        "primitive source samples =",
        primitive_source_histogram.keys,
    )
    print("target has ratios        =", has_ratio_samples)
    print("exclude                  =", exclude_samples)
    print("profile_only             =", profile_only)
    print("mask                     =", args.mask)
    print("mask_spec                =", mask_spec)
    print("throw_mode               =", args.throw_mode)
    print("ntoys                    =", args.ntoys)
    print("seed                     =", args.seed)
    print("lambdas                  =", lambda_values)
    print("nprof_values             =", nprof_values)
    print("total flux universes     =", total_flux_universes)
    print(
        "minimum primitive sansFlux eigenvalue before clipping =",
        minimum_sansflux_eigenvalue,
    )
    print("output csv               =", out_csv)
    print("")


    # --------------------------------------------------
    # Verify that the direct source exactly reconstructs
    # the target configuration at nominal.
    # --------------------------------------------------
    primitive_nominal = np.asarray(
        primitive_source_histogram.GetMCHistogram(),
        dtype=float,
    )[1:-1]

    rebuilt_target_nominal = convert_primitive_vector_to_target(
        primitive_nominal,
        target_has_ratio=has_ratio_samples,
    )

    stored_target_nominal = np.asarray(
        target_source_histogram.GetMCHistogram(),
        dtype=float,
    )[1:-1]

    if len(rebuilt_target_nominal) != len(stored_target_nominal):
        raise RuntimeError(
            "Rebuilt target nominal has {} bins, but stored target "
            "has {} bins".format(
                len(rebuilt_target_nominal),
                len(stored_target_nominal),
            )
        )

    nominal_difference = (
        rebuilt_target_nominal
        - stored_target_nominal
    )

    max_nominal_abs_difference = float(
        np.max(np.abs(nominal_difference))
    )

    max_nominal_relative_difference = float(
        np.max(
            np.divide(
                np.abs(nominal_difference),
                np.maximum(
                    np.abs(stored_target_nominal),
                    1e-12,
                ),
            )
        )
    )

    print("")
    print("===== MATCHED-TOY NOMINAL CLOSURE =====")
    print("target has ratios        =", has_ratio_samples)
    print(
        "rebuilt target shape    =",
        rebuilt_target_nominal.shape,
    )
    print(
        "stored target shape     =",
        stored_target_nominal.shape,
    )
    print(
        "max absolute difference =",
        max_nominal_abs_difference,
    )
    print(
        "max relative difference =",
        max_nominal_relative_difference,
    )
    print("===== END NOMINAL CLOSURE =====")
    print("")

    if max_nominal_relative_difference > 1e-8:
        raise RuntimeError(
            "Primitive source does not reproduce the target nominal. "
            "Do not run matched toys until this closure is understood."
        )


    rows = []

    for toy_index in range(args.ntoys):
        toy_seed = args.seed + toy_index
        toy_rng = np.random.default_rng(toy_seed)

        (
            primitive_nominal_vector,
            primitive_pseudo_vector,
        ) = make_primitive_null_pseudodata(
            source_histogram=primitive_source_histogram,
            throw_mode=args.throw_mode,
            rng=toy_rng,
            sqrt_sansflux=sqrt_sansflux,
        )

        nominal_vector = convert_primitive_vector_to_target(
            primitive_nominal_vector,
            target_has_ratio=has_ratio_samples,
        )

        pseudo_vector = convert_primitive_vector_to_target(
            primitive_pseudo_vector,
            target_has_ratio=has_ratio_samples,
        )

        print("")
        print("################################################")
        print(
            "Toy {}/{}".format(
                toy_index,
                args.ntoys - 1,
            )
        )
        print("################################################")
        print("toy_index        =", toy_index)
        print("toy_seed         =", toy_seed)
        print("nominal sum      =", np.sum(nominal_vector))
        print("pseudo sum       =", np.sum(pseudo_vector))
        print(
            "sum pseudo-CV    =",
            np.sum(pseudo_vector - nominal_vector),
        )
        print(
            "norm pseudo-CV   =",
            np.linalg.norm(
                pseudo_vector - nominal_vector
            ),
        )
        print(
            "max |pseudo-CV|  =",
            np.max(
                np.abs(
                    pseudo_vector - nominal_vector
                )
            ),
        )

        # This exact pseudo-vector is reused across every Nprof and lambda.
        # Load one fresh histogram for this toy only.
        sample_histogram = StitchedHistogram(
            "sample_toy_{}".format(toy_index)
        )

        shutil.copyfile(
            target_hist_config,
            "HIST_CONFIG.json",
        )

        sample_histogram.Load(
            target_file_path
        )

        set_pseudo_histogram_from_vector(
            sample_histogram,
            pseudo_vector,
            toy_index,
        )

        check_pseudo = np.asarray(
            sample_histogram.GetPseudoHistogram(),
            dtype=float,
        )[1:-1]

        max_pseudo_difference = float(
            np.max(
                np.abs(
                    check_pseudo - pseudo_vector
                )
            )
        )

        if max_pseudo_difference > 1e-10:
            raise RuntimeError(
                "Pseudo histogram assignment failed: "
                "max difference = {}".format(
                    max_pseudo_difference
                )
            )

        for nprof in nprof_values:
            for lam in lambda_values:
                print("")
                print(
                    "toy={} Nprof={} lambda={}".format(
                        toy_index,
                        nprof,
                        lam,
                    )
                )

                row = run_one_point(
                    sample_histogram=sample_histogram,
                    lam=lam,
                    nprof=nprof,
                    exclude_samples=exclude_samples,
                    mask_spec=mask_spec,
                    profile_only=profile_only,
                )

                row.update({
                    "toy_index": int(toy_index),
                    "toy_seed": int(toy_seed),
                    "throw_mode": args.throw_mode,
                    "base_seed": int(args.seed),

                    "hist_config_tag": plot_tag,
                    "toy_source_hist_config_tag": toy_source_tag,
                    "target_has_ratio": int(has_ratio_samples),
                    "matched_primitive_toy": 1,

                    "exclude": (
                        exclude_samples
                        if exclude_samples != ""
                        else "none"
                    ),
                    "profile_only": (
                        profile_only
                        if profile_only is not None
                        else "none"
                    ),
                    "mask": args.mask,

                    "primitive_nominal_sum": float(
                        np.sum(primitive_nominal_vector)
                    ),
                    "primitive_pseudo_sum": float(
                        np.sum(primitive_pseudo_vector)
                    ),
                    "primitive_shift_norm": float(
                        np.linalg.norm(
                            primitive_pseudo_vector
                            - primitive_nominal_vector
                        )
                    ),
                    "nominal_sum": float(
                        np.sum(nominal_vector)
                    ),
                    "pseudo_sum": float(
                        np.sum(pseudo_vector)
                    ),
                    "pseudo_minus_nominal_sum": float(
                        np.sum(
                            pseudo_vector
                            - nominal_vector
                        )
                    ),
                    "pseudo_minus_nominal_norm": float(
                        np.linalg.norm(
                            pseudo_vector
                            - nominal_vector
                        )
                    ),
                    "pseudo_minus_nominal_max_abs": float(
                        np.max(
                            np.abs(
                                pseudo_vector
                                - nominal_vector
                            )
                        )
                    ),
                    "n_flux_universes_total": int(
                        total_flux_universes
                    ),
                })

                rows.append(row)

                print("")
                print("===== toy lambda result =====")
                for key in [
                    "toy_index",
                    "nprof",
                    "lambda",
                    "chi2_null",
                    "resid_null",
                    "penalty_null",
                    "norm_a_null",
                    "max_abs_a_null",
                ]:
                    print(
                        "{:20s} {}".format(
                            key,
                            row[key],
                        )
                    )

        # Save only after the full toy is complete.
        write_rows_csv(
            rows,
            out_csv,
        )

        print(
            "Wrote partial results through toy",
            toy_index,
        )

        del sample_histogram
        gc.collect()

    print("")
    print("===== SCAN COMPLETE =====")
    print("rows written =", len(rows))
    print("saved        =", out_csv)


if __name__ == "__main__":
    main()
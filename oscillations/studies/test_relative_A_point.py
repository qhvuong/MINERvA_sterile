#!/usr/bin/env python3

import os
import sys
import argparse
import shutil
import ROOT


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--hist-config-tag", required=True)
    parser.add_argument("--exclude", default="none")
    parser.add_argument("--profile-only", default=None)
    parser.add_argument("--profile-n-universes", type=int, default=None)
    parser.add_argument("--lam", type=float, default=1.0)

    parser.add_argument("--dm2", type=float, required=True)
    parser.add_argument("--ue4", type=float, required=True)
    parser.add_argument("--umu4", type=float, required=True)
    parser.add_argument("--utau4", type=float, required=True)

    parser.add_argument(
        "--mode",
        choices=["absolute", "relative", "both"],
        default="both",
    )

    args, remaining = parser.parse_known_args()

    # Hide diagnostic-only arguments from AnalysisConfig.
    sys.argv = [sys.argv[0]] + remaining

    return args


args = parse_args()

# Only import project modules after sys.argv has been cleaned.
from tools.StitchedHistogram import StitchedHistogram
from tools.Fitters import Statistics


def normalize_exclude(exclude):
    if exclude is None:
        return ""

    if isinstance(exclude, str):
        if exclude.strip().lower() in ["", "none"]:
            return ""

    return exclude


def run_fixed_point(
    file_path,
    exclude,
    profile_only,
    profile_n_universes,
    lam,
    dm2,
    ue4,
    umu4,
    utau4,
    use_relative_a,
):
    os.environ["USE_RELATIVE_FLUX_A"] = (
        "1" if use_relative_a else "0"
    )

    histogram = StitchedHistogram("sample")
    histogram.Load(file_path)

    histogram.OscillateHistogram(
        dm2,
        ue4,
        umu4,
        utau4,
    )

    statistic = Statistics(
        histogram,
        exclude=exclude,
        lam=lam,
        mask_spec=None,
        profile_only=profile_only,
        profile_n_universes=profile_n_universes,
    )

    chi2, penalty = statistic.Chi2DataMC(
        marginalize=True,
        useOsc=True,
    )

    residual = float(chi2 - penalty)

    flux_fitter = statistic.GetFluxFitter(useOsc=True)
    solution = flux_fitter.GetFluxSolution()

    return {
        "mode": "relative" if use_relative_a else "absolute",
        "chi2": float(chi2),
        "residual": residual,
        "penalty": float(penalty),
        "norm_a": float((solution @ solution) ** 0.5),
        "max_abs_a": float(abs(solution).max()),
    }


def print_result(result):
    print("")
    print("===== {} A RESULT =====".format(result["mode"].upper()))
    print("chi2      =", result["chi2"])
    print("residual  =", result["residual"])
    print("penalty   =", result["penalty"])
    print("|a|       =", result["norm_a"])
    print("max |a_i| =", result["max_abs_a"])


def main():
    global args

    ccnueroot = os.environ.get("CCNUEROOT")
    if ccnueroot is None:
        raise RuntimeError("CCNUEROOT is not set")

    plot_tag = args.hist_config_tag

    file_path = os.path.join(
        ccnueroot,
        "oscillations",
        "rootfiles",
        "NuE_stitched_hists_{}.root".format(plot_tag),
    )

    hist_config = "HIST_CONFIG_{}.json".format(plot_tag)

    if not os.path.exists(hist_config):
        raise RuntimeError(
            "Missing requested hist config file: {}".format(
                hist_config
            )
        )

    if not os.path.exists(file_path):
        raise RuntimeError(
            "Missing stitched ROOT file: {}".format(file_path)
        )

    shutil.copyfile(
        hist_config,
        "HIST_CONFIG.json",
    )

    exclude = normalize_exclude(args.exclude)

    print("\n===== FIXED-POINT RELATIVE-A TEST =====")
    print("file                =", file_path)
    print("hist_config         =", hist_config)
    print("exclude             =", exclude)
    print("profile_only        =", args.profile_only)
    print("profile_n_universes =", args.profile_n_universes)
    print("lambda              =", args.lam)
    print("dm2                 =", args.dm2)
    print("ue4                 =", args.ue4)
    print("umu4                =", args.umu4)
    print("utau4               =", args.utau4)
    print("mode                =", args.mode)

    results = {}

    if args.mode in ["absolute", "both"]:
        results["absolute"] = run_fixed_point(
            file_path=file_path,
            exclude=exclude,
            profile_only=args.profile_only,
            profile_n_universes=args.profile_n_universes,
            lam=args.lam,
            dm2=args.dm2,
            ue4=args.ue4,
            umu4=args.umu4,
            utau4=args.utau4,
            use_relative_a=False,
        )

        print_result(results["absolute"])

    if args.mode in ["relative", "both"]:
        results["relative"] = run_fixed_point(
            file_path=file_path,
            exclude=exclude,
            profile_only=args.profile_only,
            profile_n_universes=args.profile_n_universes,
            lam=args.lam,
            dm2=args.dm2,
            ue4=args.ue4,
            umu4=args.umu4,
            utau4=args.utau4,
            use_relative_a=True,
        )

        print_result(results["relative"])

    if args.mode == "both":
        absolute = results["absolute"]
        relative = results["relative"]

        print("")
        print("===== RELATIVE MINUS ABSOLUTE =====")
        print(
            "delta chi2      =",
            relative["chi2"] - absolute["chi2"],
        )
        print(
            "delta residual  =",
            relative["residual"] - absolute["residual"],
        )
        print(
            "delta penalty   =",
            relative["penalty"] - absolute["penalty"],
        )
        print(
            "delta |a|       =",
            relative["norm_a"] - absolute["norm_a"],
        )
        print(
            "delta max |a_i| =",
            relative["max_abs_a"] - absolute["max_abs_a"],
        )

    print("\n===== END FIXED-POINT TEST =====")


if __name__ == "__main__":
    main()
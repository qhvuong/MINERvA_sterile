#!/usr/bin/env python3

import os
import sys
import argparse
import shutil
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--hist-config-tag", required=True)
    parser.add_argument("--exclude", default="none")
    parser.add_argument("--profile-only", default=None)
    parser.add_argument("--profile-n-universes", type=int, default=50)
    parser.add_argument("--lam", type=float, default=1.0)

    args, remaining = parser.parse_known_args()

    # Hide script-specific arguments before project imports parse sys.argv.
    sys.argv = [sys.argv[0]] + remaining

    return args


args = parse_args()

import ROOT

from tools.StitchedHistogram import StitchedHistogram
from tools.Fitters import Statistics, GetSampleBinRanges

ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)


def normalize_exclude(value):
    if value is None:
        return ""

    if isinstance(value, str) and value.strip().lower() in ["", "none"]:
        return ""

    return value


def print_sample_bin_diagnostics(
    histogram,
    data_vec,
    mc_vec,
    covariance,
    sample_name,
    label,
):
    """
    Print per-bin residual and correlated chi2-like contributions
    for one stitched sample using that sample's covariance sub-block.

    q_i = r_i * (V r)_i
    sum_i q_i = sample sub-block chi2

    Individual q_i values can be negative because of correlations.
    """
    ranges = GetSampleBinRanges(histogram.keys)

    if sample_name not in ranges:
        raise RuntimeError(
            "Sample '{}' not found in stitched histogram ranges".format(
                sample_name
            )
        )

    inds = ranges[sample_name]

    data_sample = np.asarray(data_vec, dtype=float)[inds]
    mc_sample = np.asarray(mc_vec, dtype=float)[inds]

    cov_sample = np.asarray(
        covariance,
        dtype=float,
    )[np.ix_(inds, inds)]

    inv_cov_sample = np.linalg.pinv(cov_sample)

    residual = data_sample - mc_sample
    weighted_residual = inv_cov_sample @ residual
    q = residual * weighted_residual

    chi2_sample = float(residual.T @ inv_cov_sample @ residual)

    print("")
    print(
        "===== PER-BIN DIAGNOSTIC: {} | {} =====".format(
            sample_name,
            label,
        )
    )
    print("sample sub-block chi2 =", chi2_sample)
    print("sum q_i               =", float(np.sum(q)))
    print("")
    print(
        "{:>5s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s}".format(
            "bin",
            "data",
            "prediction",
            "data-pred",
            "(V r)_i",
            "q_i",
        )
    )

    for local_idx in range(len(inds)):
        print(
            "{:5d} {:14.7g} {:14.7g} {:14.7g} "
            "{:14.7g} {:14.7g}".format(
                local_idx + 1,
                data_sample[local_idx],
                mc_sample[local_idx],
                residual[local_idx],
                weighted_residual[local_idx],
                q[local_idx],
            )
        )

    print(
        "===== END PER-BIN DIAGNOSTIC: {} | {} =====".format(
            sample_name,
            label,
        )
    )
    print("")

    return {
        "label": label,
        "sample": sample_name,
        "data": data_sample,
        "prediction": mc_sample,
        "residual": residual,
        "weighted_residual": weighted_residual,
        "q": q,
        "chi2": chi2_sample,
    }


def run_point(
    label,
    file_path,
    dm2,
    ue4,
    umu4,
    utau4,
    exclude,
    profile_only,
    profile_n_universes,
    lam,
    profile_flux,
    use_relative_a,
):
    print("")
    print("############################################################")
    print("POINT:", label)
    print("############################################################")
    print("dm2                 =", dm2)
    print("ue4                 =", ue4)
    print("umu4                =", umu4)
    print("utau4               =", utau4)
    print("profile_flux        =", profile_flux)
    print("use_relative_A      =", use_relative_a)
    print("exclude             =", exclude)
    print("profile_only        =", profile_only)
    print("profile_n_universes =", profile_n_universes)
    print("lambda              =", lam)

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
        marginalize=profile_flux,
        useOsc=True,
    )

    if profile_flux:
        flux_fitter = statistic.GetFluxFitter(useOsc=True)
        solution = flux_fitter.GetFluxSolution()

        norm_a = float(np.linalg.norm(solution))
        max_abs_a = float(np.max(np.abs(solution)))
    else:
        norm_a = 0.0
        max_abs_a = 0.0

    residual = float(chi2 - penalty)

    print("")
    print("===== GLOBAL RESULT: {} =====".format(label))
    print("chi2      =", float(chi2))
    print("residual  =", residual)
    print("penalty   =", float(penalty))
    print("|a|       =", norm_a)
    print("max |a_i| =", max_abs_a)

    # Use the exact prediction and covariance treatment from this point.
    if profile_flux:
        flux_fitter = statistic.GetFluxFitter(useOsc=True)
        mc_vec, _ = flux_fitter.MarginalizeFlux()
        covariance = histogram.GetCovarianceMatrix(sansFlux=True)
    else:
        mc_vec = np.asarray(
            histogram.GetOscillatedHistogram(),
            dtype=float,
        )[1:-1]

        covariance = histogram.GetCovarianceMatrix(sansFlux=False)

    data_vec = np.asarray(
        histogram.GetDataHistogram(),
        dtype=float,
    )[1:-1]

    statistic.DebugChi2BySample(
        data_vec,
        mc_vec,
        covariance,
        label=label,
    )

    sample_bin_results = {}

    for sample_name in [
        "fhc_nue_selection",
        "rhc_nue_selection",
        "fhc_numu_selection",
        "rhc_numu_selection",
    ]:
        if sample_name not in histogram.keys:
            continue

        sample_bin_results[sample_name] = print_sample_bin_diagnostics(
            histogram=histogram,
            data_vec=data_vec,
            mc_vec=mc_vec,
            covariance=covariance,
            sample_name=sample_name,
            label=label,
        )

    return {
        "label": label,
        "chi2": float(chi2),
        "residual": residual,
        "penalty": float(penalty),
        "norm_a": norm_a,
        "max_abs_a": max_abs_a,
        "sample_bins": sample_bin_results,
    }


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

    if not os.path.exists(file_path):
        raise RuntimeError(
            "Missing stitched ROOT file: {}".format(file_path)
        )

    if not os.path.exists(hist_config):
        raise RuntimeError(
            "Missing histogram config: {}".format(hist_config)
        )

    shutil.copyfile(hist_config, "HIST_CONFIG.json")

    exclude = normalize_exclude(args.exclude)

    print("\n===== FIXED-POINT SAMPLE COMPARISON =====")
    print("file                =", file_path)
    print("hist_config         =", hist_config)
    print("exclude             =", exclude)
    print("profile_only        =", args.profile_only)
    print("profile_n_universes =", args.profile_n_universes)
    print("lambda              =", args.lam)

    results = []

    # 1. Profiled null point.
    results.append(
        run_point(
            label="null_profiled",
            file_path=file_path,
            dm2=40.0,
            ue4=0.0,
            umu4=0.0,
            utau4=0.0,
            exclude=exclude,
            profile_only=args.profile_only,
            profile_n_universes=args.profile_n_universes,
            lam=args.lam,
            profile_flux=True,
            use_relative_a=True,
        )
    )

    # --------------------------------------------------
    # 2. Direct-CCnue global best fit
    #
    # Replace these four values with the result from:
    # USE_RELATIVE_FLUX_A=1 python fitData.py ...
    # --------------------------------------------------
    results.append(
        run_point(
            label="direct_globalBF",
            file_path=file_path,

            # REPLACE WITH DIRECT RELATIVE-A GLOBAL BF
            dm2=40.628,
            ue4=0.238,
            umu4=0.01568,
            utau4=0.711,

            exclude=exclude,
            profile_only=args.profile_only,
            profile_n_universes=args.profile_n_universes,
            lam=args.lam,
            profile_flux=True,
            use_relative_a=True,
        )
    )

    # --------------------------------------------------
    # 3. Direct-CCnue lower-mass local minimum
    #
    # Replace these values with an actual local minimum
    # printed by the direct relative-A multistart fit.
    # --------------------------------------------------
    results.append(
        run_point(
            label="direct_lowerBF",
            file_path=file_path,

            # REPLACE WITH DIRECT RELATIVE-A LOWER-MASS BF
            dm2=19.1066,
            ue4=0.238,
            umu4=0.01568,
            utau4=0.711,

            exclude=exclude,
            profile_only=args.profile_only,
            profile_n_universes=args.profile_n_universes,
            lam=args.lam,
            profile_flux=True,
            use_relative_a=True,
        )
    )

    # # 2. Previous no-profile best fit for the ratio configuration.
    # results.append(
    #     run_point(
    #         label="noProfile_BF_dm2_18p854",
    #         file_path=file_path,
    #         dm2=18.85427581855046,
    #         ue4=0.02410474744692575,
    #         umu4=0.03494692651028352,
    #         utau4=0.0,
    #         exclude=exclude,
    #         profile_only=args.profile_only,
    #         profile_n_universes=args.profile_n_universes,
    #         lam=args.lam,
    #         profile_flux=False,
    #         use_relative_a=False,
    #     )
    # )

    # # 3. Relative-A profiled best fit.
    # results.append(
    #     run_point(
    #         label="relativeA_BF_dm2_41p007",
    #         file_path=file_path,
    #         dm2=41.00740619806878,
    #         ue4=0.2137000197823035,
    #         umu4=0.016387800362155338,
    #         utau4=0.6111483546339311,
    #         exclude=exclude,
    #         profile_only=args.profile_only,
    #         profile_n_universes=args.profile_n_universes,
    #         lam=args.lam,
    #         profile_flux=True,
    #         use_relative_a=True,
    #     )
    # )

    # results.append(
    #     run_point(
    #         label="relativeA_at_noProfile_BF_dm2_18p854",
    #         file_path=file_path,
    #         dm2=18.85427581855046,
    #         ue4=0.02410474744692575,
    #         umu4=0.03494692651028352,
    #         utau4=0.0,
    #         exclude=exclude,
    #         profile_only=args.profile_only,
    #         profile_n_universes=args.profile_n_universes,
    #         lam=args.lam,
    #         profile_flux=True,
    #         use_relative_a=True,
    #     )
    # )

    # results.append(
    #     run_point(
    #         label="relativeA_localBF_dm2_10p723",
    #         file_path=file_path,
    #         dm2=10.7232,
    #         ue4=0.24489,
    #         umu4=0.01316,
    #         utau4=0.0,
    #         exclude=exclude,
    #         profile_only=args.profile_only,
    #         profile_n_universes=args.profile_n_universes,
    #         lam=args.lam,
    #         profile_flux=True,
    #         use_relative_a=True,
    #     )
    # )

    # results.append(
    #     run_point(
    #         label="relativeA_localBF_dm2_18p218",
    #         file_path=file_path,
    #         dm2=18.2178,
    #         ue4=0.18250,
    #         umu4=0.01575,
    #         utau4=0.19666,
    #         exclude=exclude,
    #         profile_only=args.profile_only,
    #         profile_n_universes=args.profile_n_universes,
    #         lam=args.lam,
    #         profile_flux=True,
    #         use_relative_a=True,
    #     )
    # )

    def compare_sample_bins(result_a, result_b, sample_name):
        a = result_a["sample_bins"][sample_name]
        b = result_b["sample_bins"][sample_name]

        if len(a["q"]) != len(b["q"]):
            raise RuntimeError(
                "{} bin counts do not match".format(sample_name)
            )

        if not np.allclose(
            a["data"],
            b["data"],
            rtol=0.0,
            atol=1e-12,
        ):
            raise RuntimeError(
                "{} data vectors differ between points".format(sample_name)
            )

        print("")
        print("============================================================")
        print(
            "{} BIN COMPARISON: {} -> {}".format(
                sample_name,
                result_a["label"],
                result_b["label"],
            )
        )
        print("Positive delta_q favors the second point.")
        print("============================================================")

        print(
            "{:>5s} {:>13s} {:>13s} {:>13s} "
            "{:>13s} {:>13s} {:>13s} {:>13s}".format(
                "bin",
                "data",
                "pred_A",
                "pred_B",
                "pred_B-A",
                "q_A",
                "q_B",
                "delta_q",
            )
        )

        delta_q_total = 0.0

        for i in range(len(a["q"])):
            delta_q = a["q"][i] - b["q"][i]
            delta_q_total += delta_q

            print(
                "{:5d} {:13.6g} {:13.6g} {:13.6g} "
                "{:13.6g} {:13.6g} {:13.6g} {:13.6g}".format(
                    i + 1,
                    a["data"][i],
                    a["prediction"][i],
                    b["prediction"][i],
                    b["prediction"][i] - a["prediction"][i],
                    a["q"][i],
                    b["q"][i],
                    delta_q,
                )
            )

        print("")
        print("sum delta_q =", delta_q_total)
        print(
            "sub-block chi2 difference =",
            a["chi2"] - b["chi2"],
        )


    print("")
    print("============================================================")
    print("GLOBAL SUMMARY")
    print("============================================================")

    for result in results:
        print(
            "{:<30s} chi2={:12.6f} residual={:12.6f} "
            "penalty={:10.6f} |a|={:10.6f}".format(
                result["label"],
                result["chi2"],
                result["residual"],
                result["penalty"],
                result["norm_a"],
            )
        )

    result_map = {
        result["label"]: result
        for result in results
    }

    samples_to_compare = [
        "fhc_nue_selection",
        "rhc_nue_selection",
        "fhc_numu_selection",
        "rhc_numu_selection",
    ]

    for sample_name in samples_to_compare:
        if (
            sample_name not in result_map["direct_lowerBF"]["sample_bins"]
            or sample_name not in result_map["direct_globalBF"]["sample_bins"]
        ):
            print(
                "Skipping comparison for missing sample:",
                sample_name,
            )
            continue

        compare_sample_bins(
            result_map["direct_lowerBF"],
            result_map["direct_globalBF"],
            sample_name,
        )


if __name__ == "__main__":
    main()
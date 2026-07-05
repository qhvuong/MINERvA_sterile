#!/usr/bin/env python

import os
import argparse
import numpy as np

from tools.StitchedHistogram import StitchedHistogram
from tools.Helper import sin_average
from tools.Fitters import Statistics

np.set_printoptions(precision=6, linewidth=200)


def safe_integral(h):
    return h.Integral() if h else 0.0


def get_bin_values(h):
    return np.array([h.GetBinContent(i) for i in range(1, h.GetNbinsX() + 1)])


def classify_sample(sample):
    if "elastic" in sample:
        return "elastic"
    if "imd" in sample or "numu" in sample:
        return "numu"
    if "ratio" in sample:
        return "ratio"
    if "nue" in sample:
        return "nue"
    return "other"


def manual_oscillated_bin(
    sample,
    nue,
    numu,
    swap,
    nue_sin,
    numu_sin,
    swap_sin,
    ue4,
    umu4,
    utau4,
):
    """
    Reproduce the sample-level logic in StitchedHistogram.OscillateSubHistogram.

    elastic:
      total = nue*Pee + swap*Pmue + numu*Pmumu + nue*Petau + numu*Pmutau

    numu / imd:
      total = numu*Pmumu

    nue:
      total = nue*Pee + swap*Pmue

    ratio:
      skipped by default in this script
    """
    Pee = 1.0 - 4.0 * ue4 * (1.0 - ue4) * nue_sin
    Pmue = 4.0 * ue4 * umu4 * swap_sin
    Pmumu = 1.0 - 4.0 * umu4 * (1.0 - umu4) * numu_sin
    Petau = 4.0 * ue4 * utau4 * nue_sin
    Pmutau = 4.0 * utau4 * umu4 * numu_sin

    kind = classify_sample(sample)

    if kind == "elastic":
        manual = (
            nue * Pee
            + swap * Pmue
            + numu * Pmumu
            + nue * Petau
            + numu * Pmutau
        )
    elif kind == "numu":
        manual = numu * Pmumu
    elif kind == "nue":
        manual = nue * Pee + swap * Pmue
    else:
        manual = np.nan

    return manual, Pee, Pmue, Pmumu, Petau, Pmutau


def check_one_sample(h, sample, dm2, ue4, umu4, utau4, verbose=False):
    if sample not in h.mc_hists:
        print("\nSkipping missing sample:", sample)
        return None

    if "ratio" in sample:
        print("\nSkipping derived ratio sample:", sample)
        return None

    required_dicts = [
        h.nue_hists,
        h.numu_hists,
        h.swap_hists,
        h.nue_templates,
        h.numu_templates,
        h.swap_templates,
    ]

    if any(sample not in d for d in required_dicts):
        print("\nSkipping sample with missing components/templates:", sample)
        return None

    h_nue = h.nue_hists[sample]
    h_numu = h.numu_hists[sample]
    h_swap = h.swap_hists[sample]

    t_nue = h.nue_templates[sample]
    t_numu = h.numu_templates[sample]
    t_swap = h.swap_templates[sample]

    try:
        _, h_code_total = h.OscillateSubHistogram(sample, dm2, ue4, umu4, utau4)
    except Exception as e:
        print("\nSkipping sample {} because OscillateSubHistogram failed: {}".format(sample, e))
        return None

    manual_integral = 0.0
    code_integral = 0.0
    null_integral = h.mc_hists[sample].Integral()

    max_abs_diff = 0.0
    max_frac_diff = 0.0

    rows = []

    for i in range(1, h_nue.GetNbinsX() + 1):
        nue = h_nue.GetBinContent(i)
        numu = h_numu.GetBinContent(i)
        swap = h_swap.GetBinContent(i)

        # Per-sample templates have L/E on X and reco bin on Y,
        # so use yaxis=False, matching OscillateSubHistogram.
        nue_sin = sin_average(i, dm2, t_nue, False)
        numu_sin = sin_average(i, dm2, t_numu, False)
        swap_sin = sin_average(i, dm2, t_swap, False)

        manual, Pee, Pmue, Pmumu, Petau, Pmutau = manual_oscillated_bin(
            sample,
            nue,
            numu,
            swap,
            nue_sin,
            numu_sin,
            swap_sin,
            ue4,
            umu4,
            utau4,
        )

        code = h_code_total.GetBinContent(i)
        diff = code - manual
        fracdiff = diff / manual if manual not in [0, np.nan] and manual != 0 else 0.0

        manual_integral += manual
        code_integral += code

        max_abs_diff = max(max_abs_diff, abs(diff))
        max_frac_diff = max(max_frac_diff, abs(fracdiff))

        rows.append(
            {
                "bin": i,
                "nue": nue,
                "numu": numu,
                "swap": swap,
                "nue_sin": nue_sin,
                "numu_sin": numu_sin,
                "swap_sin": swap_sin,
                "Pee": Pee,
                "Pmue": Pmue,
                "Pmumu": Pmumu,
                "Petau": Petau,
                "Pmutau": Pmutau,
                "manual": manual,
                "code": code,
                "diff": diff,
                "fracdiff": fracdiff,
            }
        )

    osc_frac_effect = (
        (code_integral - null_integral) / null_integral
        if null_integral != 0
        else 0.0
    )

    print("\n===== {} =====".format(sample))
    print("  type              =", classify_sample(sample))
    print("  null integral      =", null_integral)
    print("  manual osc integral=", manual_integral)
    print("  code osc integral  =", code_integral)
    print("  code - manual      =", code_integral - manual_integral)
    print("  osc - null         =", code_integral - null_integral)
    print("  frac osc effect    =", osc_frac_effect)
    print("  max abs bin diff   =", max_abs_diff)
    print("  max frac bin diff  =", max_frac_diff)
    print("  component integrals:")
    print("    nue  =", safe_integral(h_nue))
    print("    numu =", safe_integral(h_numu))
    print("    swap =", safe_integral(h_swap))

    if verbose:
        print(
            "\n{:>3s} {:>10s} {:>10s} {:>10s} | {:>8s} {:>8s} {:>8s} | {:>8s} {:>8s} {:>8s} | {:>11s} {:>11s} {:>11s} {:>10s}".format(
                "bin",
                "nue",
                "numu",
                "swap",
                "s_nue",
                "s_numu",
                "s_swap",
                "Pee",
                "Pmue",
                "Pmumu",
                "manual",
                "code",
                "diff",
                "fracdiff",
            )
        )

        for r in rows:
            print(
                "{bin:3d} {nue:10.4g} {numu:10.4g} {swap:10.4g} | {nue_sin:8.5f} {numu_sin:8.5f} {swap_sin:8.5f} | {Pee:8.5f} {Pmue:8.5f} {Pmumu:8.5f} | {manual:11.5g} {code:11.5g} {diff:11.4e} {fracdiff:10.3e}".format(
                    **r
                )
            )

    return {
        "sample": sample,
        "kind": classify_sample(sample),
        "null": null_integral,
        "manual": manual_integral,
        "code": code_integral,
        "code_minus_manual": code_integral - manual_integral,
        "osc_minus_null": code_integral - null_integral,
        "frac_osc_effect": osc_frac_effect,
        "max_abs_bin_diff": max_abs_diff,
        "max_frac_bin_diff": max_frac_diff,
        "nue_int": safe_integral(h_nue),
        "numu_int": safe_integral(h_numu),
        "swap_int": safe_integral(h_swap),
    }


def check_global_stitched(h, dm2, ue4, umu4, utau4):
    """
    Check full stitched OscillateHistogram result by sample blocks using
    HIST_CONFIG.json, which stores the actual stitched bin ranges.
    """
    from tools.Helper import GetSliceIndices

    h.OscillateHistogram(dm2, ue4, umu4, utau4)

    h_null = h.GetMCHistogram()
    h_osc = h.GetOscillatedHistogram()

    print("\n\n===== Global stitched block summary =====")

    for sample in h.keys:
        try:
            inds = GetSliceIndices("HIST_CONFIG.json", "", [sample])
        except Exception as e:
            print("  Could not get indices for {}: {}".format(sample, e))
            continue

        if len(inds) == 0:
            continue

        null_sum = 0.0
        osc_sum = 0.0

        for idx0 in inds:
            ibin = idx0 + 1
            null_sum += h_null.GetBinContent(ibin)
            osc_sum += h_osc.GetBinContent(ibin)

        frac = (osc_sum - null_sum) / null_sum if null_sum != 0 else 0.0

        print(
            "  {:20s} bins {:3d}-{:3d}: null={:12.6g} osc={:12.6g} diff={:12.6g} frac={: .5%}".format(
                sample,
                inds[0] + 1,
                inds[-1] + 1,
                null_sum,
                osc_sum,
                osc_sum - null_sum,
                frac,
            )
        )


def check_chi2_profile_path(h, dm2, ue4, umu4, utau4, exclude="ratio,elastic", lam=1):
    print("\n\n===== Chi2 / flux profiling check =====")
    print("exclude =", exclude)
    print("lambda  =", lam)

    # Null/raw
    stat_null = Statistics(h, exclude=exclude, lam=lam)
    chi2_null_raw, pen_null_raw = stat_null.Chi2DataMC(
        marginalize=False,
        useOsc=False,
    )

    # Oscillated/raw
    h.OscillateHistogram(dm2, ue4, umu4, utau4)
    stat_osc = Statistics(h, exclude=exclude, lam=lam)
    chi2_osc_raw, pen_osc_raw = stat_osc.Chi2DataMC(
        marginalize=False,
        useOsc=True,
    )

    # Null/profiled
    stat_null_prof = Statistics(h, exclude=exclude, lam=lam)
    chi2_null_prof, pen_null_prof = stat_null_prof.Chi2DataMC(
        marginalize=True,
        useOsc=False,
    )

    # Oscillated/profiled
    stat_osc_prof = Statistics(h, exclude=exclude, lam=lam)
    chi2_osc_prof, pen_osc_prof = stat_osc_prof.Chi2DataMC(
        marginalize=True,
        useOsc=True,
    )

    print("\nRaw:")
    print("  null chi2       = {:.6f}  penalty = {:.6f}".format(chi2_null_raw, pen_null_raw))
    print("  osc chi2        = {:.6f}  penalty = {:.6f}".format(chi2_osc_raw, pen_osc_raw))
    print("  delta chi2 raw  = {:.6f}".format(chi2_null_raw - chi2_osc_raw))

    print("\nProfiled:")
    print("  null chi2       = {:.6f} = {:.6f} + {:.6f} penalty".format(
        chi2_null_prof, chi2_null_prof - pen_null_prof, pen_null_prof
    ))
    print("  osc chi2        = {:.6f} = {:.6f} + {:.6f} penalty".format(
        chi2_osc_prof, chi2_osc_prof - pen_osc_prof, pen_osc_prof
    ))
    print("  delta chi2 prof = {:.6f}".format(chi2_null_prof - chi2_osc_prof))

    # Optional: show how much profiling moved the MC prediction.
    null_ff = stat_null_prof.GetFluxFitter(useOsc=False)
    osc_ff = stat_osc_prof.GetFluxFitter(useOsc=True)

    if null_ff is not None:
        null_profiled_mc, _ = null_ff.MarginalizeFlux()
        raw_null_mc = np.array(h.GetMCHistogram())[1:-1]
        shift = null_profiled_mc - raw_null_mc
        print("\nNull profiling MC shift:")
        print("  total shift   =", np.sum(shift))
        print("  max abs shift =", np.max(np.abs(shift)))
        print("  first 6 v+e bin shifts =", shift[:6])

    if osc_ff is not None:
        osc_profiled_mc, _ = osc_ff.MarginalizeFlux()
        raw_osc_mc = np.array(h.GetOscillatedHistogram())[1:-1]
        shift = osc_profiled_mc - raw_osc_mc
        print("\nOsc profiling MC shift:")
        print("  total shift   =", np.sum(shift))
        print("  max abs shift =", np.max(np.abs(shift)))
        print("  first 6 v+e bin shifts =", shift[:6])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        default=None,
        help="Stitched ROOT file. Default: $CCNUEROOT/oscillations/rootfiles/NuE_stitched_hists_newNuE.root",
    )
    parser.add_argument("--dm2", type=float, default=5.0)
    parser.add_argument("--ue4", type=float, default=0.05)
    parser.add_argument("--umu4", type=float, default=0.05)
    parser.add_argument("--utau4", type=float, default=0.0)
    parser.add_argument("--profile-exclude", default="ratio,elastic")
    parser.add_argument(
        "--sample",
        default="all",
        help="Sample to check, e.g. fhc_elastic, fhc_numu_selection, fhc_nue_selection, or all.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print bin-by-bin details for each sample.",
    )
    args = parser.parse_args()

    if args.input is None:
        args.input = os.path.join(
            os.environ["CCNUEROOT"],
            "oscillations/rootfiles/NuE_stitched_hists_newNuE.root",
        )

    print("Input file:", args.input)
    print("Oscillation point:")
    print("  dm2   =", args.dm2)
    print("  ue4   =", args.ue4)
    print("  umu4  =", args.umu4)
    print("  utau4 =", args.utau4)

    h = StitchedHistogram("sample")
    h.Load(args.input)

    if args.sample == "all":
        samples = [s for s in h.keys if "ratio" not in s]
    else:
        samples = [args.sample]

    results = []
    for sample in samples:
        res = check_one_sample(
            h,
            sample,
            args.dm2,
            args.ue4,
            args.umu4,
            args.utau4,
            verbose=args.verbose,
        )
        if res is not None:
            results.append(res)

    print("\n\n===== Compact summary =====")
    print(
        "{:<22s} {:>10s} {:>12s} {:>12s} {:>12s} {:>12s}".format(
            "sample",
            "type",
            "null",
            "osc",
            "frac_eff",
            "code-manual",
        )
    )
    for r in results:
        print(
            "{:<22s} {:>10s} {:12.5g} {:12.5g} {:12.5%} {:12.4e}".format(
                r["sample"],
                r["kind"],
                r["null"],
                r["code"],
                r["frac_osc_effect"],
                r["code_minus_manual"],
            )
        )

    check_global_stitched(h, args.dm2, args.ue4, args.umu4, args.utau4)

    check_chi2_profile_path(
        h,
        args.dm2,
        args.ue4,
        args.umu4,
        args.utau4,
        exclude=args.profile_exclude,
        lam=1,
    )

if __name__ == "__main__":
    main()
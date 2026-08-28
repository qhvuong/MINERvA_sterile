#!/usr/bin/env python

import argparse
import json
import ROOT
import PlotUtils
import numpy as np

from tools.PlotLibrary import HistHolder
from tools import Utilities


ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)


SIGNAL_DEFINITION = [
    "CCNuEQE",
    "CCNuEDelta",
    "CCNuEDIS",
    "CCNuE",
    "CCNuE2p2h",
    "CCNuEWrongSign",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Inspect raw CCnue/CCnuebar selection Flux universes "
            "before background fit and stitching."
        )
    )

    parser.add_argument(
        "--config",
        default="stitch_input_files.json",
        help="Stitch input JSON",
    )

    parser.add_argument(
        "--input-set",
        required=True,
        choices=[
            "p8",
            "p8_onlyPPFX",
            "p8_onlyBeamFocus",
        ],
    )

    parser.add_argument(
        "--sample",
        required=True,
        choices=[
            "fhc_ccnue",
            "rhc_ccnuebar",
        ],
    )

    parser.add_argument(
        "--hist",
        default="Biased Neutrino Energy",
        help="HistHolder plot name",
    )

    parser.add_argument(
        "--band",
        default="Flux",
        help="Vertical error band to inspect",
    )

    parser.add_argument(
        "--ntuple-tag",
        default="MAD",
    )

    parser.add_argument(
        "--n-universes-print",
        type=int,
        default=5,
    )

    return parser.parse_args()


def get_config_path(config_file, input_set, sample):
    with open(config_file, "r") as f:
        cfg = json.load(f)

    return cfg[input_set][sample]["mc"]


def load_raw_holder(path, sample, plot_name, ntuple_tag):
    if sample == "fhc_ccnue":
        playlist_name = "CCnue_allSystematics_fullStatsFluxes"
    elif sample == "rhc_ccnuebar":
        playlist_name = "CCnuebar_allSystematics_fullStatsFluxes"
    else:
        raise RuntimeError("Unknown sample {}".format(sample))

    type_path_map = {
        "mc": path,
    }

    _, mc_file, _, _, mc_pot = Utilities.getFilesAndPOTScale(
        playlist_name,
        type_path_map,
        ntuple_tag,
        True,
    )

    if mc_pot is None or mc_pot == 0:
        raise RuntimeError("Could not determine MC POT")

    # No data file here, so just use MC POT as the normalization target.
    stand_pot = mc_pot

    holder = HistHolder(
        plot_name,
        mc_file,
        "Signal",
        True,
        mc_pot,
        stand_pot,
    )

    # Since stand_pot == mc_pot, this should be a scale of 1.
    holder.POTScale(False)

    return mc_file, holder


def get_delta(hist, band_name):
    if not hist.HasVertErrorBand(band_name):
        raise RuntimeError(
            "{} has no band '{}'. Available = {}".format(
                hist.GetName(),
                band_name,
                list(hist.GetVertErrorBandNames()),
            )
        )

    band = hist.GetVertErrorBand(band_name)

    nuniv = band.GetNHists()
    nbins = hist.GetNbinsX()

    cv = np.array(
        [
            hist.GetBinContent(i)
            for i in range(1, nbins + 1)
        ],
        dtype=float,
    )

    universes = np.zeros((nuniv, nbins))

    for u in range(nuniv):
        hu = band.GetHist(u)

        for i in range(1, nbins + 1):
            universes[u, i - 1] = hu.GetBinContent(i)

    delta = np.full_like(universes, np.nan)

    good = cv != 0.0

    delta[:, good] = (
        universes[:, good] - cv[good]
    ) / cv[good]

    return cv, universes, delta


def build_signal_sum(holder, name):
    hsum = None

    for category in SIGNAL_DEFINITION:
        if category not in holder.hists:
            print(
                "WARNING: category {} is missing".format(
                    category
                )
            )
            continue

        h = holder.hists[category]

        if h is None:
            continue

        if hsum is None:
            hsum = h.Clone(name)
            hsum.Reset()
            hsum.SetDirectory(0)

        hsum.Add(h)

    if hsum is None:
        raise RuntimeError("Could not construct signal sum")

    return hsum


def print_cv_composition(holder, signal_sum):
    nbins = signal_sum.GetNbinsX()

    print("")
    print("=" * 150)
    print("CV SIGNAL COMPOSITION")
    print("=" * 150)

    header = "{:>4s} {:>9s} {:>9s} {:>12s}".format(
        "bin",
        "low",
        "high",
        "total",
    )

    for category in SIGNAL_DEFINITION:
        header += " {:>14s}".format(category)

    header += " {:>12s}".format("WS frac")

    print(header)
    print("-" * 150)

    for i in range(1, nbins + 1):

        total = signal_sum.GetBinContent(i)

        values = []

        for category in SIGNAL_DEFINITION:
            if category in holder.hists and holder.hists[category]:
                value = holder.hists[category].GetBinContent(i)
            else:
                value = 0.0

            values.append(value)

        if "CCNuEWrongSign" in holder.hists:
            ws = holder.hists["CCNuEWrongSign"].GetBinContent(i)
        else:
            ws = 0.0

        ws_frac = ws / total if total != 0.0 else np.nan

        line = (
            "{:4d} {:9.3f} {:9.3f} {:12.6g}".format(
                i,
                signal_sum.GetXaxis().GetBinLowEdge(i),
                signal_sum.GetXaxis().GetBinUpEdge(i),
                total,
            )
        )

        for value in values:
            line += " {:14.6g}".format(value)

        line += " {:11.5f}%".format(
            100.0 * ws_frac
        )

        print(line)


def print_flux_summary(hist, band_name, label):
    cv, universes, delta = get_delta(
        hist,
        band_name,
    )

    print("")
    print("=" * 110)
    print("FLUX SUMMARY: {}".format(label))
    print("=" * 110)

    print(
        "{:>4s} {:>9s} {:>9s} {:>12s} "
        "{:>13s} {:>13s} {:>13s} {:>13s}".format(
            "bin",
            "low",
            "high",
            "CV",
            "<delta>",
            "sigma",
            "min",
            "max",
        )
    )

    print("-" * 110)

    for i in range(len(cv)):
        d = delta[:, i]

        print(
            "{:4d} {:9.3f} {:9.3f} {:12.6g} "
            "{:+12.6e} {:12.6e} {:+12.6e} {:+12.6e}".format(
                i + 1,
                hist.GetXaxis().GetBinLowEdge(i + 1),
                hist.GetXaxis().GetBinUpEdge(i + 1),
                cv[i],
                np.nanmean(d),
                np.nanstd(d, ddof=1),
                np.nanmin(d),
                np.nanmax(d),
            )
        )

    return cv, universes, delta


def print_selected_universes(
    hist,
    band_name,
    label,
    nprint=5,
):
    cv, universes, delta = get_delta(
        hist,
        band_name,
    )

    print("")
    print("=" * 110)
    print("SELECTED UNIVERSES: {}".format(label))
    print("=" * 110)

    bins_to_show = [
        0,
        min(5, len(cv) - 1),
        min(8, len(cv) - 1),
        min(10, len(cv) - 1),
        len(cv) - 1,
    ]

    bins_to_show = sorted(set(bins_to_show))

    header = "{:>5s}".format("u")

    for idx in bins_to_show:
        header += " {:>14s}".format(
            "delta_bin{}".format(idx + 1)
        )

    print(header)
    print("-" * len(header))

    for u in range(min(nprint, delta.shape[0])):
        line = "{:5d}".format(u)

        for idx in bins_to_show:
            line += " {:+14.6e}".format(
                delta[u, idx]
            )

        print(line)


def print_flatness_test(hist, band_name, label):
    """
    For each universe, calculate the spread of delta across reco bins.

    A perfectly flat multiplicative weight gives:
        std_bins(delta_u) = 0
    """

    _, _, delta = get_delta(
        hist,
        band_name,
    )

    per_universe_bin_spread = np.nanstd(
        delta,
        axis=1,
        ddof=1,
    )

    print("")
    print("=" * 110)
    print("ENERGY-FLATNESS TEST: {}".format(label))
    print("=" * 110)

    print(
        "mean std across bins = {:.8e}".format(
            np.nanmean(per_universe_bin_spread)
        )
    )

    print(
        "median               = {:.8e}".format(
            np.nanmedian(per_universe_bin_spread)
        )
    )

    print(
        "max                  = {:.8e}".format(
            np.nanmax(per_universe_bin_spread)
        )
    )

    worst = np.nanargmax(per_universe_bin_spread)

    print(
        "worst universe       = {}  std_bins(delta)={:.8e}".format(
            worst,
            per_universe_bin_spread[worst],
        )
    )


def main():
    args = parse_args()

    mc_path = get_config_path(
        args.config,
        args.input_set,
        args.sample,
    )

    print("")
    print("Input set =", args.input_set)
    print("Sample    =", args.sample)
    print("MC file   =", mc_path)
    print("Histogram =", args.hist)
    print("Band      =", args.band)

    mc_file, holder = load_raw_holder(
        mc_path,
        args.sample,
        args.hist,
        args.ntuple_tag,
    )

    print("")
    print("Available categories:")
    for category in holder.hists:
        h = holder.hists[category]

        if h:
            print(
                "  {:25s} integral = {:12.6g}".format(
                    category,
                    h.Integral(),
                )
            )

    signal_sum = build_signal_sum(
        holder,
        "raw_signal_sum",
    )

    # --------------------------------------------------
    # CV composition
    # --------------------------------------------------

    print_cv_composition(
        holder,
        signal_sum,
    )

    # --------------------------------------------------
    # Total signal Flux response
    #
    # This should reproduce what BackgroundFit writes as
    # EN4_predicted_Signal, modulo histogram choice/scaling.
    # --------------------------------------------------

    print_flux_summary(
        signal_sum,
        args.band,
        "TOTAL SIGNAL",
    )

    print_selected_universes(
        signal_sum,
        args.band,
        "TOTAL SIGNAL",
        args.n_universes_print,
    )

    print_flatness_test(
        signal_sum,
        args.band,
        "TOTAL SIGNAL",
    )

    # --------------------------------------------------
    # Individual signal categories
    # --------------------------------------------------

    for category in SIGNAL_DEFINITION:

        if category not in holder.hists:
            continue

        h = holder.hists[category]

        if h is None or h.Integral() == 0:
            continue

        print_flux_summary(
            h,
            args.band,
            category,
        )

        print_flatness_test(
            h,
            args.band,
            category,
        )

    mc_file.Close()

    print("")
    print("Done.")


if __name__ == "__main__":
    main()
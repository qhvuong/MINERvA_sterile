#!/usr/bin/env python3

import argparse
import json
import os
import shutil
import sys

import numpy as np
import ROOT

# from tools.StitchedHistogram import StitchedHistogram


ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)


def parse_dm2_values(text):
    values = []
    for token in text.split(","):
        token = token.strip()
        if not token:
            continue
        value = float(token)
        if value < 0:
            raise argparse.ArgumentTypeError("dm2 values must be non-negative")
        values.append(value)

    if not values:
        raise argparse.ArgumentTypeError("Provide at least one dm2 value")

    return values


def load_hist_config(path):
    with open(path, "r") as handle:
        return json.load(handle)


def get_sample_range(config, sample):
    if sample not in config:
        return None

    start0 = int(config[sample]["start"])
    end0 = int(config[sample]["end"])

    # ROOT stitched bins are one-based.
    return start0 + 1, end0 + 1


def make_ratio_hist(osc_hist, null_hist, name):
    ratio = osc_hist.Clone(name)
    ratio.SetDirectory(0)
    ratio.Divide(ratio, null_hist)
    return ratio


def style_hist(hist, color, line_style):
    hist.SetLineColor(color)
    hist.SetLineWidth(3)
    hist.SetLineStyle(line_style)
    hist.SetMarkerStyle(0)
    hist.SetFillStyle(0)


def draw_sample_panel(
    sample,
    sample_range,
    ratio_hists,
    dm2_values,
    colors,
    line_styles,
    output_path,
):
    start_bin, end_bin = sample_range

    canvas = ROOT.TCanvas(
        "c_{}".format(sample),
        "c_{}".format(sample),
        1000,
        700,
    )
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.05)
    canvas.SetBottomMargin(0.13)
    canvas.SetTopMargin(0.08)

    frame = ratio_hists[0].Clone("frame_{}".format(sample))
    frame.Reset()
    frame.SetDirectory(0)
    frame.SetTitle(
        "{} raw oscillated / raw null;Stitched bin number;Oscillated / null".format(
            sample
        )
    )
    frame.GetXaxis().SetRange(start_bin, end_bin)
    frame.SetMinimum(0.0)
    frame.SetMaximum(1.6)
    frame.Draw("AXIS")

    legend = ROOT.TLegend(0.58, 0.68, 0.93, 0.91)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.031)

    kept = []

    for index, (dm2, hist) in enumerate(zip(dm2_values, ratio_hists)):
        style_hist(
            hist,
            colors[index % len(colors)],
            line_styles[index % len(line_styles)],
        )
        hist.GetXaxis().SetRange(start_bin, end_bin)
        hist.Draw("HIST SAME")
        legend.AddEntry(
            hist,
            "#Delta m^{{2}} = {:g} eV^{{2}}".format(dm2),
            "l",
        )
        kept.append(hist)

    line = ROOT.TLine(
        ratio_hists[0].GetXaxis().GetBinLowEdge(start_bin),
        1.0,
        ratio_hists[0].GetXaxis().GetBinUpEdge(end_bin),
        1.0,
    )
    line.SetLineColor(ROOT.kBlack)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw("SAME")

    legend.Draw()
    canvas.Print(output_path)

    return canvas, frame, legend, line, kept


def print_sample_values(sample, sample_range, ratio_arrays, dm2_values):
    start_bin, end_bin = sample_range
    start0 = start_bin - 1
    end0 = end_bin - 1

    print("")
    print("===== {} raw oscillated / raw null =====".format(sample))

    header = "{:>6s}".format("lbin")
    for dm2 in dm2_values:
        header += " {:>14s}".format("dm2={:g}".format(dm2))
    print(header)

    for global_idx0 in range(start0, end0 + 1):
        local_bin = global_idx0 - start0 + 1
        row = "{:6d}".format(local_bin)

        for array in ratio_arrays:
            row += " {:14.8f}".format(array[global_idx0])

        print(row)


def print_convergence(
    sample,
    sample_range,
    ratio_arrays,
    dm2_values,
):
    if len(dm2_values) < 2:
        return

    start_bin, end_bin = sample_range
    sl = slice(start_bin - 1, end_bin)

    print("")
    print("===== {} dm2-to-dm2 convergence =====".format(sample))

    for previous_index in range(len(dm2_values) - 1):
        dm2_a = dm2_values[previous_index]
        dm2_b = dm2_values[previous_index + 1]

        a = ratio_arrays[previous_index][sl]
        b = ratio_arrays[previous_index + 1][sl]

        frac = np.divide(
            b,
            a,
            out=np.ones_like(b),
            where=np.abs(a) > 1e-12,
        )

        print(
            "{:g} -> {:g}: max |B/A - 1| = {:.8e}, "
            "mean |B/A - 1| = {:.8e}".format(
                dm2_a,
                dm2_b,
                np.max(np.abs(frac - 1.0)),
                np.mean(np.abs(frac - 1.0)),
            )
        )


def print_swap_ratio_check(
    histogram,
    config,
    ue4,
    umu4,
):
    """
    Directly inspect the stored nominal electron and swapped components
    for the FHC/RHC ratio bins.

    At fully averaged oscillations:
        Pee   = 1 - 2*ue4*(1-ue4)
        Pmumu = 1 - 2*umu4*(1-umu4)
        Pmue  = 2*ue4*umu4

    Here ue4 and umu4 mean |Ue4|^2 and |Umu4|^2.
    """

    if histogram.nue_hist is None:
        raise RuntimeError("Loaded stitched file has no nue_hist")

    if histogram.swap_hist is None:
        raise RuntimeError("Loaded stitched file has no swap_hist")

    if histogram.numu_hist is None:
        raise RuntimeError("Loaded stitched file has no numu_hist")

    pee = 1.0 - 2.0 * ue4 * (1.0 - ue4)
    pmumu = 1.0 - 2.0 * umu4 * (1.0 - umu4)
    pmue = 2.0 * ue4 * umu4

    print("")
    print("===== HIGH-DM2 AVERAGED PROBABILITIES =====")
    print("Pee   =", pee)
    print("Pmumu =", pmumu)
    print("Pmue  =", pmue)

    for sample in ["fhc_ratio", "rhc_ratio"]:
        sample_range = get_sample_range(config, sample)

        if sample_range is None:
            continue

        start_bin, end_bin = sample_range

        print("")
        print(
            "===== {} DIRECT Nswap / Ne CHECK =====".format(
                sample
            )
        )

        print(
            "{:>6s} {:>14s} {:>14s} {:>14s} {:>18s}".format(
                "lbin",
                "Ne",
                "Nswap",
                "Nswap/Ne",
                "expected Rosc/R0",
            )
        )

        for global_bin in range(start_bin, end_bin + 1):
            local_bin = global_bin - start_bin + 1

            ne = histogram.nue_hist.GetBinContent(global_bin)
            nswap = histogram.swap_hist.GetBinContent(global_bin)

            if abs(ne) < 1e-12:
                swap_over_nue = float("nan")
                expected = float("nan")
            else:
                swap_over_nue = nswap / ne

                denominator = (
                    pee
                    + pmue * swap_over_nue
                )

                expected = (
                    pmumu / denominator
                    if abs(denominator) > 1e-12
                    else float("nan")
                )

            print(
                "{:6d} {:14.8g} {:14.8g} {:14.8f} {:18.8f}".format(
                    local_bin,
                    ne,
                    nswap,
                    swap_over_nue,
                    expected,
                )
            )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot raw oscillated/raw-null predictions at fixed dm2 values. "
            "No fit and no flux profiling are performed."
        )
    )

    parser.add_argument(
        "--hist-config-tag",
        default="prodNueel_p8",
        help="Tag used by NuE_stitched_hists_<tag>.root and HIST_CONFIG_<tag>.json",
    )
    parser.add_argument(
        "--dm2",
        type=parse_dm2_values,
        default=parse_dm2_values("41.718,100,1000,10000"),
        help="Comma-separated dm2 values in eV^2",
    )
    parser.add_argument("--ue4", type=float, default=0.18939364174315096)
    parser.add_argument("--umu4", type=float, default=0.01763413048189484)
    parser.add_argument("--utau4", type=float, default=0.2574217107021343)
    parser.add_argument(
        "--out-dir",
        default="plots/fixed_oscillation",
        help="Output directory",
    )

    args = parser.parse_args()

    # Importing StitchedHistogram imports AnalysisConfig, whose global parser
    # also tries to parse sys.argv. Hide this diagnostic script's arguments
    # during that import.
    saved_argv = sys.argv[:]

    try:
        sys.argv = [sys.argv[0]]
        from tools.StitchedHistogram import StitchedHistogram
    finally:
        sys.argv = saved_argv

    ccnueroot = os.environ.get("CCNUEROOT")
    if ccnueroot is None:
        raise RuntimeError("CCNUEROOT is not set")

    plot_tag = args.hist_config_tag

    root_path = os.path.join(
        ccnueroot,
        "oscillations",
        "rootfiles",
        "NuE_stitched_hists_{}.root".format(plot_tag),
    )

    hist_config_path = "HIST_CONFIG_{}.json".format(plot_tag)

    if not os.path.exists(root_path):
        raise RuntimeError("Missing stitched ROOT file: {}".format(root_path))

    if not os.path.exists(hist_config_path):
        raise RuntimeError(
            "Missing histogram configuration: {}".format(hist_config_path)
        )

    os.makedirs(args.out_dir, exist_ok=True)

    # StitchedHistogram helper functions expect this exact filename.
    shutil.copyfile(hist_config_path, "HIST_CONFIG.json")
    config = load_hist_config("HIST_CONFIG.json")

    print("")
    print("===== fixed raw oscillation diagnostic =====")
    print("file       =", root_path)
    print("hist config=", hist_config_path)
    print("dm2 values =", args.dm2)
    print("Ue4^2      =", args.ue4)
    print("Umu4^2     =", args.umu4)
    print("Utau4^2    =", args.utau4)
    print("profiling  = OFF")
    print("")

    histogram = StitchedHistogram("fixed_oscillation")
    histogram.Load(root_path)

    print_swap_ratio_check(
        histogram=histogram,
        config=config,
        ue4=args.ue4,
        umu4=args.umu4,
    )

    null_hist = histogram.GetMCHistogram()
    null_array = np.asarray(null_hist, dtype=float)[1:-1]

    osc_hists = []
    ratio_hists = []
    ratio_arrays = []

    for index, dm2 in enumerate(args.dm2):
        histogram.OscillateHistogram(
            dm2,
            args.ue4,
            args.umu4,
            args.utau4,
            False,
            False,
        )

        osc_hist = histogram.GetOscillatedHistogram()
        osc_hist.SetName("osc_dm2_{}".format(index))
        osc_hist.SetDirectory(0)

        ratio_hist = make_ratio_hist(
            osc_hist,
            null_hist,
            "ratio_dm2_{}".format(index),
        )

        osc_hists.append(osc_hist)
        ratio_hists.append(ratio_hist)
        ratio_arrays.append(
            np.asarray(ratio_hist, dtype=float)[1:-1]
        )

    colors = [
        ROOT.kOrange + 7,
        ROOT.kBlue + 1,
        ROOT.kMagenta + 1,
        ROOT.kGreen + 2,
        ROOT.kCyan + 2,
        ROOT.kRed + 1,
    ]
    line_styles = [2, 1, 3, 4, 5, 6]

    samples = [
        "fhc_ratio",
        "rhc_ratio",
    ]

    for sample in samples:
        sample_range = get_sample_range(config, sample)

        if sample_range is None:
            print(
                "WARNING: {} is not present in HIST_CONFIG.json; skipping".format(
                    sample
                )
            )
            continue

        print_sample_values(
            sample,
            sample_range,
            ratio_arrays,
            args.dm2,
        )

        print_convergence(
            sample,
            sample_range,
            ratio_arrays,
            args.dm2,
        )

        output_path = os.path.join(
            args.out_dir,
            "{}_{}_raw_osc_over_null.png".format(
                plot_tag,
                sample,
            ),
        )

        draw_sample_panel(
            sample,
            sample_range,
            ratio_hists,
            args.dm2,
            colors,
            line_styles,
            output_path,
        )

        print("wrote", output_path)

    # Also save the full stitched comparison.
    full_output = os.path.join(
        args.out_dir,
        "{}_all_stitched_raw_osc_over_null.png".format(plot_tag),
    )

    canvas = ROOT.TCanvas("c_full", "c_full", 1200, 700)
    canvas.SetLeftMargin(0.10)
    canvas.SetRightMargin(0.04)
    canvas.SetBottomMargin(0.12)

    frame = ratio_hists[0].Clone("frame_full")
    frame.Reset()
    frame.SetTitle(
        "Raw oscillated / raw null;Stitched bin number;Oscillated / null"
    )
    frame.SetMinimum(0.0)
    frame.SetMaximum(1.6)
    frame.Draw("AXIS")

    legend = ROOT.TLegend(0.68, 0.68, 0.94, 0.91)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)

    for index, (dm2, hist) in enumerate(zip(args.dm2, ratio_hists)):
        style_hist(
            hist,
            colors[index % len(colors)],
            line_styles[index % len(line_styles)],
        )
        hist.Draw("HIST SAME")
        legend.AddEntry(
            hist,
            "#Delta m^{{2}} = {:g} eV^{{2}}".format(dm2),
            "l",
        )

    unity = ROOT.TLine(
        frame.GetXaxis().GetXmin(),
        1.0,
        frame.GetXaxis().GetXmax(),
        1.0,
    )
    unity.SetLineColor(ROOT.kBlack)
    unity.SetLineStyle(2)
    unity.Draw("SAME")

    legend.Draw()
    canvas.Print(full_output)

    print("wrote", full_output)
    print("")
    print("Done.")


if __name__ == "__main__":
    main()
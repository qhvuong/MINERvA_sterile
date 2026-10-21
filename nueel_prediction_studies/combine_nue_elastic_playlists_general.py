#!/usr/bin/env python3

import argparse
import ROOT
from array import array

import PlotUtils.LoadPlotUtilsLib

ROOT.gROOT.SetBatch(True)

HIST_NAMES_DEFAULT = [
    "h_nue_elastic_total",
    "h_nue_elastic_nue",
    "h_nue_elastic_nuebar",
    "h_nue_elastic_numu",
    "h_nue_elastic_numubar",
]


def get_edges(h):
    edges = []

    for i in range(1, h.GetNbinsX() + 1):
        edges.append(h.GetBinLowEdge(i))

    edges.append(
        h.GetBinLowEdge(h.GetNbinsX())
        + h.GetBinWidth(h.GetNbinsX())
    )

    return edges


def check_binning(h1, h2, name):
    if h1.GetNbinsX() != h2.GetNbinsX():
        raise RuntimeError(
            f"{name}: bin-count mismatch: "
            f"{h1.GetNbinsX()} vs {h2.GetNbinsX()}"
        )

    for i in range(1, h1.GetNbinsX() + 2):
        e1 = h1.GetBinLowEdge(i)
        e2 = h2.GetBinLowEdge(i)

        if abs(e1 - e2) > 1.0e-12:
            raise RuntimeError(
                f"{name}: bin-edge mismatch at edge {i}: "
                f"input1={e1}, input2={e2}"
            )


def has_flux_band(h, flux_band_name):
    try:
        return h.HasVertErrorBand(flux_band_name)
    except Exception:
        return False


def make_plain_hist(template, name):
    edges = get_edges(template)

    h = ROOT.TH1D(
        name,
        template.GetTitle(),
        len(edges) - 1,
        array("d", edges),
    )

    h.SetDirectory(0)

    h.GetXaxis().SetTitle(
        template.GetXaxis().GetTitle()
    )

    h.GetYaxis().SetTitle(
        template.GetYaxis().GetTitle()
    )

    return h


def make_mnv_hist(
    template,
    name,
    flux_band_name,
    n_universes,
):
    edges = get_edges(template)

    ROOT.gROOT.cd()

    h = ROOT.PlotUtils.MnvH1D(
        name,
        template.GetTitle(),
        len(edges) - 1,
        array("d", edges),
    )

    ROOT.SetOwnership(h, False)
    h.SetDirectory(0)

    h.GetXaxis().SetTitle(
        template.GetXaxis().GetTitle()
    )

    h.GetYaxis().SetTitle(
        template.GetYaxis().GetTitle()
    )

    h.AddVertErrorBand(
        flux_band_name,
        n_universes,
    )

    return h


def combine_plain_hist(h1, h2, name):
    check_binning(h1, h2, name)

    hout = make_plain_hist(
        h1,
        name,
    )

    for i in range(
        0,
        hout.GetNbinsX() + 2,
    ):
        value = (
            h1.GetBinContent(i)
            + h2.GetBinContent(i)
        )

        hout.SetBinContent(
            i,
            value,
        )

        # Optional: combine stored bin errors in quadrature.
        err1 = h1.GetBinError(i)
        err2 = h2.GetBinError(i)

        hout.SetBinError(
            i,
            (err1 * err1 + err2 * err2) ** 0.5,
        )

    return hout


def combine_mnv_hist(
    h1,
    h2,
    name,
    flux_band_name,
):
    check_binning(
        h1,
        h2,
        name,
    )

    if not has_flux_band(
        h1,
        flux_band_name,
    ):
        raise RuntimeError(
            f"{name}: input1 has no "
            f"'{flux_band_name}' band"
        )

    if not has_flux_band(
        h2,
        flux_band_name,
    ):
        raise RuntimeError(
            f"{name}: input2 has no "
            f"'{flux_band_name}' band"
        )

    b1 = h1.GetVertErrorBand(
        flux_band_name
    )

    b2 = h2.GetVertErrorBand(
        flux_band_name
    )

    n1 = b1.GetNHists()
    n2 = b2.GetNHists()

    if n1 != n2:
        raise RuntimeError(
            f"{name}: universe mismatch: "
            f"input1={n1}, input2={n2}"
        )

    n_universes = n1

    hout = make_mnv_hist(
        h1,
        name,
        flux_band_name,
        n_universes,
    )

    # --------------------------------------------------
    # CV
    # --------------------------------------------------

    for i in range(
        0,
        hout.GetNbinsX() + 2,
    ):
        value = (
            h1.GetBinContent(i)
            + h2.GetBinContent(i)
        )

        hout.SetBinContent(
            i,
            value,
        )

        err1 = h1.GetBinError(i)
        err2 = h2.GetBinError(i)

        hout.SetBinError(
            i,
            (err1 * err1 + err2 * err2) ** 0.5,
        )

    # --------------------------------------------------
    # Flux universes
    # --------------------------------------------------

    bout = hout.GetVertErrorBand(
        flux_band_name
    )

    for u in range(n_universes):

        hu1 = b1.GetHist(u)
        hu2 = b2.GetHist(u)
        huout = bout.GetHist(u)

        for i in range(
            0,
            hout.GetNbinsX() + 2,
        ):
            value = (
                hu1.GetBinContent(i)
                + hu2.GetBinContent(i)
            )

            huout.SetBinContent(
                i,
                value,
            )

    return hout


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Combine two nu-e prediction ROOT files. "
            "Supports both TH1D CV-only files and "
            "MnvH1D files with a Flux vertical error band."
        )
    )

    parser.add_argument(
        "--input1",
        required=True,
        help="First input ROOT file.",
    )

    parser.add_argument(
        "--input2",
        required=True,
        help="Second input ROOT file.",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Combined output ROOT file.",
    )

    parser.add_argument(
        "--flux-band-name",
        default="Flux",
        help="Flux vertical error-band name.",
    )

    parser.add_argument(
        "--hist-names",
        nargs="+",
        default=HIST_NAMES_DEFAULT,
        help="Histogram names to combine.",
    )

    args = parser.parse_args()

    f1 = ROOT.TFile.Open(
        args.input1
    )

    if not f1 or f1.IsZombie():
        raise RuntimeError(
            f"Could not open {args.input1}"
        )

    f2 = ROOT.TFile.Open(
        args.input2
    )

    if not f2 or f2.IsZombie():
        raise RuntimeError(
            f"Could not open {args.input2}"
        )

    print()
    print("Combining nu-e predictions")
    print(f"  input1: {args.input1}")
    print(f"  input2: {args.input2}")
    print(f"  output: {args.output}")
    print()

    output_hists = []

    for name in args.hist_names:

        h1 = f1.Get(name)
        h2 = f2.Get(name)

        if not h1:
            raise RuntimeError(
                f"Could not find {name} "
                f"in {args.input1}"
            )

        if not h2:
            raise RuntimeError(
                f"Could not find {name} "
                f"in {args.input2}"
            )

        h1_has_flux = has_flux_band(
            h1,
            args.flux_band_name,
        )

        h2_has_flux = has_flux_band(
            h2,
            args.flux_band_name,
        )

        if h1_has_flux != h2_has_flux:
            raise RuntimeError(
                f"{name}: inconsistent histogram types: "
                f"input1 has Flux band={h1_has_flux}, "
                f"input2 has Flux band={h2_has_flux}"
            )

        if h1_has_flux:
            hout = combine_mnv_hist(
                h1,
                h2,
                name,
                args.flux_band_name,
            )

            hist_type = "MnvH1D"
        else:
            hout = combine_plain_hist(
                h1,
                h2,
                name,
            )

            hist_type = "TH1D"

        ROOT.SetOwnership(
            hout,
            False,
        )

        output_hists.append(
            hout
        )

        expected = (
            h1.Integral()
            + h2.Integral()
        )

        closure = (
            hout.Integral()
            - expected
        )

        print(
            f"{name:30s} "
            f"[{hist_type}]  "
            f"input1={h1.Integral():10.5f}  "
            f"input2={h2.Integral():10.5f}  "
            f"combined={hout.Integral():10.5f}  "
            f"closure={closure:+.3e}"
        )

    fout = ROOT.TFile(
        args.output,
        "RECREATE",
    )

    for hout in output_hists:
        fout.cd()
        hout.Write()

    fout.Close()
    f1.Close()
    f2.Close()

    print()
    print(
        f"Wrote combined prediction to "
        f"{args.output}"
    )


if __name__ == "__main__":
    main()
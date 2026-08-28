#!/usr/bin/env python3

import ROOT
from array import array

import PlotUtils.LoadPlotUtilsLib


INPUT_LE1 = "/exp/minerva/data/users/qvuong/nueel_prediction_studies/nue_elastic_prediction_TomalakXS_mnv_LE1.root"
INPUT_LE13 = "/exp/minerva/data/users/qvuong/nueel_prediction_studies/nue_elastic_prediction_TomalakXS_mnv_LE13.root"

OUTPUT = "/exp/minerva/data/users/qvuong/nueel_prediction_studies/nue_elastic_prediction_TomalakXS_mnv_FHC.root"

FLUX_BAND = "Flux"

HIST_NAMES = [
    "h_nue_elastic_total",
    "h_nue_elastic_nue",
    "h_nue_elastic_nuebar",
    "h_nue_elastic_numu",
    "h_nue_elastic_numubar",
]


def make_output_hist(template, name, n_universes):
    """
    Make a fresh MnvH1D with the same binning as the input histogram
    and one Flux error band.
    """

    edges = []

    for i in range(1, template.GetNbinsX() + 1):
        edges.append(template.GetBinLowEdge(i))

    edges.append(
        template.GetBinLowEdge(template.GetNbinsX())
        + template.GetBinWidth(template.GetNbinsX())
    )

    # Construct the MnvH1D and its error-band histograms
    # outside the output TFile directory.
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
        FLUX_BAND,
        n_universes,
    )

    return h


def combine_hist(h1, h2, name):
    """
    Combine two MnvH1Ds:

        combined = LE1 + LE13

    No additional POT scaling is applied.

    The same Flux universe index is combined across playlists.
    """

    if h1.GetNbinsX() != h2.GetNbinsX():
        raise RuntimeError(
            f"Binning mismatch for {name}"
        )

    # Check actual bin edges.
    for i in range(1, h1.GetNbinsX() + 2):
        if abs(h1.GetBinLowEdge(i) - h2.GetBinLowEdge(i)) > 1e-12:
            raise RuntimeError(
                f"{name}: bin-edge mismatch at edge {i}: "
                f"LE1={h1.GetBinLowEdge(i)}, "
                f"LE13={h2.GetBinLowEdge(i)}"
            )

    if not h1.HasVertErrorBand(FLUX_BAND):
        raise RuntimeError(
            f"{name}: LE1 histogram has no {FLUX_BAND} band"
        )

    if not h2.HasVertErrorBand(FLUX_BAND):
        raise RuntimeError(
            f"{name}: LE13 histogram has no {FLUX_BAND} band"
        )

    b1 = h1.GetVertErrorBand(FLUX_BAND)
    b2 = h2.GetVertErrorBand(FLUX_BAND)

    n1 = b1.GetNHists()
    n2 = b2.GetNHists()

    if n1 != n2:
        raise RuntimeError(
            f"{name}: universe mismatch: "
            f"LE1={n1}, LE13={n2}"
        )

    n_universes = n1

    hout = make_output_hist(
        h1,
        name,
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

    # --------------------------------------------------
    # Flux universes
    # --------------------------------------------------

    bout = hout.GetVertErrorBand(
        FLUX_BAND
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

    f1 = ROOT.TFile.Open(INPUT_LE1)

    if not f1 or f1.IsZombie():
        raise RuntimeError(
            f"Could not open {INPUT_LE1}"
        )

    f13 = ROOT.TFile.Open(INPUT_LE13)

    if not f13 or f13.IsZombie():
        raise RuntimeError(
            f"Could not open {INPUT_LE13}"
        )

    print()
    print("Combining nu-e predictions")
    print("No additional POT scaling will be applied.")
    print()

    # --------------------------------------------------
    # First build everything in memory.
    # --------------------------------------------------

    ROOT.gROOT.cd()

    output_hists = []

    for name in HIST_NAMES:

        h1 = f1.Get(name)
        h13 = f13.Get(name)

        if not h1:
            raise RuntimeError(
                f"Could not find {name} in {INPUT_LE1}"
            )

        if not h13:
            raise RuntimeError(
                f"Could not find {name} in {INPUT_LE13}"
            )

        hout = combine_hist(
            h1,
            h13,
            name,
        )

        ROOT.SetOwnership(hout, False)
        output_hists.append(hout)

        expected = h1.Integral() + h13.Integral()

        print(
            f"{name:30s} "
            f"LE1={h1.Integral():10.5f}  "
            f"LE13={h13.Integral():10.5f}  "
            f"combined={hout.Integral():10.5f}  "
            f"closure={hout.Integral() - expected:+.3e}"
        )

    # --------------------------------------------------
    # Only now open the output file and write.
    # --------------------------------------------------

    fout = ROOT.TFile(
        OUTPUT,
        "RECREATE",
    )

    for hout in output_hists:
        fout.cd()
        hout.Write()

    fout.Close()

    f1.Close()
    f13.Close()

    print()
    print(
        f"Wrote combined prediction to {OUTPUT}"
    )


if __name__ == "__main__":
    main()
#!/usr/bin/env python3

import argparse
import ROOT
from array import array


# ------------------------------------------------------------------
# Jaewon paper values
# ------------------------------------------------------------------

# Actual published acceptance-corrected counts:
#
# bins   = [0.8-2, 2-3, 3-5, 5-7, 7-9, 9-inf]
# counts = [48.7, 14.4, 20.5, 18.1, 11.9, 21.6]
#
# For VISUALIZATION ONLY, we approximate the plotted 9-20 GeV
# paper bin as 15 events.
#
# The prediction's displayed 9-20 GeV bin is then obtained by
# preserving the prediction/data ratio from the actual 9-inf bin:
#
#   N_pred_plot(9-20)
#       = 15 * N_pred(9-inf) / 21.6
#
# This approximation must NOT be used in the fit.

PAPER_LAST_BIN_ACTUAL = 21.6
PAPER_LAST_BIN_PLOT = 15.0

PLOT_EDGES = [0.8, 2.0, 3.0, 5.0, 7.0, 9.0, 20.0]

PAPER_DATA_COUNTS = {
    "edges": PLOT_EDGES,
    "counts": [48.7, 14.4, 20.5, 18.1, 11.9, PAPER_LAST_BIN_PLOT],
}


def get_hist(root_file, hist_name):
    f = ROOT.TFile.Open(root_file)

    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open file: {root_file}")

    h = f.Get(hist_name)

    if not h:
        raise RuntimeError(
            f"Could not find histogram {hist_name} in {root_file}"
        )

    h.SetDirectory(0)
    f.Close()

    return h


def make_jaewon_style_prediction_hist(
    h,
    name_suffix="_jaewon_plot",
    last_bin_scale=None,
):
    """
    Construct a plotting-only histogram with edges

        [0.8, 2, 3, 5, 7, 9, 20] GeV.

    Bins 1-5 are copied directly from the physical prediction.

    Bin 6 of the physical prediction represents 9-inf
    numerically as 9-100 GeV.

    For visualization only, the final displayed bin is interpreted
    as 9-20 GeV.

    For the total prediction:

        N_plot(9-20)
          = 15 * N_pred(9-inf) / 21.6

    For flavor components, pass the corresponding scale factor so
    the flavor fractions in the last bin are preserved and the stack
    sums exactly to the displayed total.

    This transformation is for plotting only and must never be used
    as a fit prediction.
    """

    h_plot = ROOT.TH1D(
        h.GetName() + name_suffix,
        h.GetTitle(),
        len(PLOT_EDGES) - 1,
        array("d", PLOT_EDGES),
    )

    h_plot.SetDirectory(0)

    # First five bins are unchanged.
    for i in range(1, 6):
        h_plot.SetBinContent(i, h.GetBinContent(i))
        h_plot.SetBinError(i, h.GetBinError(i))

    # Sixth physical bin = 9-inf (numerically 9-100).
    old_last = h.GetBinContent(6)
    old_last_err = h.GetBinError(6)

    if last_bin_scale is None:
        # Total prediction:
        #
        # preserve ratio relative to paper's actual 9-inf count.
        target_last = (
            PAPER_LAST_BIN_PLOT
            * old_last
            / PAPER_LAST_BIN_ACTUAL
        )

        if old_last != 0.0:
            last_bin_scale = target_last / old_last
        else:
            last_bin_scale = 0.0

    h_plot.SetBinContent(
        6,
        old_last * last_bin_scale,
    )

    h_plot.SetBinError(
        6,
        old_last_err * last_bin_scale,
    )

    return h_plot, last_bin_scale


def scale_to_per_2gev(h, name_suffix="_per2gev"):
    """
    Clone histogram and scale each bin to events / 2 GeV:

        N_i_plot = N_i * 2 / bin_width
    """

    h_scaled = h.Clone(h.GetName() + name_suffix)
    h_scaled.SetDirectory(0)

    for i in range(1, h_scaled.GetNbinsX() + 1):
        width = h_scaled.GetBinWidth(i)
        content = h_scaled.GetBinContent(i)
        error = h_scaled.GetBinError(i)

        if width > 0:
            scale = 2.0 / width

            h_scaled.SetBinContent(
                i,
                content * scale,
            )

            h_scaled.SetBinError(
                i,
                error * scale,
            )

    return h_scaled


def make_paper_data_graph():
    edges = PAPER_DATA_COUNTS["edges"]
    counts = PAPER_DATA_COUNTS["counts"]

    n = len(counts)

    g = ROOT.TGraphAsymmErrors(n)

    for i in range(n):
        x_low = edges[i]
        x_high = edges[i + 1]

        width = x_high - x_low
        x = 0.5 * (x_low + x_high)

        y = counts[i] * 2.0 / width

        ex_low = x - x_low
        ex_high = x_high - x

        g.SetPoint(i, x, y)

        # Only horizontal bin widths are shown.
        # Paper uncertainties are not included in this quick overlay.
        g.SetPointError(
            i,
            ex_low,
            ex_high,
            0.0,
            0.0,
        )

    g.SetName("g_paper_data_per2gev")
    g.SetTitle(
        "Paper data;"
        "E_{l}^{true} [GeV];"
        "Events / 2 GeV"
    )

    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.1)
    g.SetLineWidth(2)

    return g


def style_hist(
    h,
    color,
    width=2,
    linestyle=1,
    fill=False,
):
    h.SetLineColor(color)
    h.SetLineWidth(width)
    h.SetLineStyle(linestyle)
    h.SetMarkerStyle(0)

    if fill:
        h.SetFillColor(color)
        h.SetFillStyle(1001)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot nu-e elastic prediction as events / 2 GeV "
            "using a visualization-only Jaewon-style 9-20 GeV "
            "last bin."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help=(
            "Input ROOT file from "
            "make_nue_elastic_from_flux_root*.py."
        ),
    )

    parser.add_argument(
        "--output",
        default="nue_elastic_prediction_per2gev.png",
        help="Output plot file, e.g. pdf/png.",
    )

    parser.add_argument(
        "--total-hist",
        default="h_nue_elastic_total",
        help="Total prediction histogram name.",
    )

    parser.add_argument(
        "--overlay-flavors",
        action="store_true",
        help=(
            "Overlay nue, nuebar, numu, numubar "
            "components if present."
        ),
    )

    parser.add_argument(
        "--overlay-paper-data",
        action="store_true",
        help=(
            "Overlay paper acceptance-corrected event counts, "
            "using the visualization-only ~15-event "
            "9-20 GeV last bin."
        ),
    )

    parser.add_argument(
        "--xmax",
        type=float,
        default=20.0,
        help="Maximum x-axis value to show.",
    )

    parser.add_argument(
        "--logy",
        action="store_true",
        help="Use log y-axis.",
    )

    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    # --------------------------------------------------------------
    # Total prediction
    # --------------------------------------------------------------

    h_total_physical = get_hist(
        args.input,
        args.total_hist,
    )

    physical_last = h_total_physical.GetBinContent(6)

    h_total_plot, last_bin_scale = (
        make_jaewon_style_prediction_hist(
            h_total_physical
        )
    )

    visual_last = h_total_plot.GetBinContent(6)

    print()
    print("Visualization-only last-bin treatment:")
    print(
        f"  physical prediction 9-inf = "
        f"{physical_last:.6f}"
    )
    print(
        f"  paper 9-inf count         = "
        f"{PAPER_LAST_BIN_ACTUAL:.6f}"
    )
    print(
        f"  assumed paper 9-20 count  = "
        f"{PAPER_LAST_BIN_PLOT:.6f}"
    )
    print(
        f"  displayed pred 9-20       = "
        f"{visual_last:.6f}"
    )
    print(
        f"  last-bin scale factor      = "
        f"{last_bin_scale:.6f}"
    )
    print()

    h_total = scale_to_per_2gev(
        h_total_plot
    )

    h_total.GetXaxis().SetTitle(
        "True lepton energy E_{l}^{true} [GeV]"
    )

    h_total.GetYaxis().SetTitle(
        "Events / 2 GeV"
    )

    h_total.SetTitle("")

    style_hist(
        h_total,
        ROOT.kBlack,
        width=3,
    )

    # --------------------------------------------------------------
    # Flavor components
    # --------------------------------------------------------------

    flavor_hists = {}

    if args.overlay_flavors:
        flavor_colors = {
            "numu": ROOT.kBlue + 1,
            "nue": ROOT.kRed + 1,
            "numubar": ROOT.kGreen + 2,
            "nuebar": ROOT.kMagenta + 1,
        }

        for flavor, color in flavor_colors.items():
            hname = f"h_nue_elastic_{flavor}"

            try:
                h_physical = get_hist(
                    args.input,
                    hname,
                )

            except RuntimeError:
                print(
                    f"Warning: could not find "
                    f"{hname}, skipping."
                )
                continue

            # Apply the same last-bin scaling as the total.
            # This preserves the flavor composition.
            h_plot, _ = (
                make_jaewon_style_prediction_hist(
                    h_physical,
                    name_suffix="_jaewon_plot",
                    last_bin_scale=last_bin_scale,
                )
            )

            h = scale_to_per_2gev(
                h_plot
            )

            style_hist(
                h,
                color,
                width=1,
                fill=True,
            )

            flavor_hists[flavor] = h

    # --------------------------------------------------------------
    # Canvas
    # --------------------------------------------------------------

    c = ROOT.TCanvas(
        "c",
        "c",
        900,
        700,
    )

    c.SetLeftMargin(0.13)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.07)

    if args.logy:
        c.SetLogy()

    xmin = 0.0
    xmax = args.xmax

    # --------------------------------------------------------------
    # Determine y range
    # --------------------------------------------------------------

    ymax = 0.0

    for i in range(
        1,
        h_total.GetNbinsX() + 1,
    ):
        x1 = h_total.GetBinLowEdge(i)
        x2 = x1 + h_total.GetBinWidth(i)

        if x2 <= xmin or x1 >= xmax:
            continue

        ymax = max(
            ymax,
            h_total.GetBinContent(i),
        )

    paper_graph = None

    if args.overlay_paper_data:
        paper_graph = make_paper_data_graph()

        for i in range(
            paper_graph.GetN()
        ):
            y = paper_graph.GetPointY(i)

            ymax = max(
                ymax,
                y,
            )

    if args.logy:
        ymin = 0.01
        ymax_plot = ymax * 20.0

    else:
        ymin = 0.0
        ymax_plot = ymax * 1.35

    # --------------------------------------------------------------
    # Draw
    # --------------------------------------------------------------

    frame = c.DrawFrame(
        xmin,
        ymin,
        xmax,
        ymax_plot,
        (
            ";True lepton energy "
            "E_{l}^{true} [GeV];"
            "Events / 2 GeV"
        ),
    )

    stack = ROOT.THStack(
        "stack_components",
        "",
    )

    for flavor in [
        "numu",
        "nue",
        "numubar",
        "nuebar",
    ]:
        if flavor in flavor_hists:
            stack.Add(
                flavor_hists[flavor]
            )

    if len(flavor_hists) > 0:
        stack.Draw("hist same")

    # Total prediction as black outline.
    h_total.Draw("hist same")

    if paper_graph:
        paper_graph.Draw("P same")

    # --------------------------------------------------------------
    # Legend
    # --------------------------------------------------------------

    leg = ROOT.TLegend(
        0.55,
        0.62,
        0.88,
        0.88,
    )

    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)

    leg.AddEntry(
        h_total,
        "Total prediction",
        "l",
    )

    flavor_labels = {
        "numu": "#nu_{#mu} e",
        "nue": "#nu_{e} e",
        "numubar": "#bar{#nu}_{#mu} e",
        "nuebar": "#bar{#nu}_{e} e",
    }

    for flavor in [
        "numu",
        "nue",
        "numubar",
        "nuebar",
    ]:
        if flavor in flavor_hists:
            leg.AddEntry(
                flavor_hists[flavor],
                flavor_labels[flavor],
                "f",
            )

    if paper_graph:
        leg.AddEntry(
            paper_graph,
            "Paper data, acc.-corrected",
            "p",
        )

    leg.Draw()

    # --------------------------------------------------------------
    # Labels
    # --------------------------------------------------------------

    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextSize(0.038)

    label.DrawLatex(
        0.16,
        0.86,
        "#nu e^{-} #rightarrow #nu e^{-}",
    )

    label.DrawLatex(
        0.16,
        0.81,
        "Truth-level prediction",
    )

    label.DrawLatex(
        0.16,
        0.76,
        "Bin contents scaled to events / 2 GeV",
    )

    label.SetTextSize(0.028)

    label.DrawLatex(
        0.16,
        0.71,
        "9-20 GeV bin approximated for visualization only",
    )

    c.SaveAs(
        args.output
    )

    print(
        f"Wrote plot to {args.output}"
    )


if __name__ == "__main__":
    main()
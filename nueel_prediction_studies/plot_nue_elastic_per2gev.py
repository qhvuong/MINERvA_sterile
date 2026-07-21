#!/usr/bin/env python3

import argparse
import ROOT


PAPER_DATA_COUNTS = {
    # Acceptance-corrected MINERvA nu-e event counts from the paper table.
    # Bins: 0.8-2, 2-3, 3-5, 5-7, 7-9, 9-inf
    "edges": [0.8, 2.0, 3.0, 5.0, 7.0, 9.0, 120.0],
    "counts": [48.7, 14.4, 20.5, 18.1, 11.9, 21.6],
}


def get_hist(root_file, hist_name):
    f = ROOT.TFile.Open(root_file)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open file: {root_file}")

    h = f.Get(hist_name)
    if not h:
        raise RuntimeError(f"Could not find histogram {hist_name} in {root_file}")

    h.SetDirectory(0)
    f.Close()
    return h


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
            h_scaled.SetBinContent(i, content * scale)
            h_scaled.SetBinError(i, error * scale)

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

        # Only x-widths here. No paper uncertainty included in this quick overlay.
        ex_low = x - x_low
        ex_high = x_high - x

        g.SetPoint(i, x, y)
        g.SetPointError(i, ex_low, ex_high, 0.0, 0.0)

    g.SetName("g_paper_data_per2gev")
    g.SetTitle("Paper data;E_{l}^{true} [GeV];Events / 2 GeV")
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.1)
    g.SetLineWidth(2)

    return g


def style_hist(h, color, width=2, linestyle=1, fill=False):
    h.SetLineColor(color)
    h.SetLineWidth(width)
    h.SetLineStyle(linestyle)
    h.SetMarkerStyle(0)

    if fill:
        h.SetFillColor(color)
        h.SetFillStyle(1001)


def main():
    parser = argparse.ArgumentParser(
        description="Plot nu-e elastic prediction as events / 2 GeV."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input ROOT file from make_nue_elastic_from_flux_root.py.",
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
        help="Overlay nue, nuebar, numu, numubar components if present.",
    )

    parser.add_argument(
        "--overlay-paper-data",
        action="store_true",
        help="Overlay paper acceptance-corrected event counts, scaled to events / 2 GeV.",
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

    h_total = get_hist(args.input, args.total_hist)
    h_total = scale_to_per_2gev(h_total)

    h_total.GetXaxis().SetTitle("True lepton energy E_{l}^{true} [GeV]")
    h_total.GetYaxis().SetTitle("Predicted events / 2 GeV")
    h_total.SetTitle("")

    style_hist(h_total, ROOT.kBlack, width=3)

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
                h = get_hist(args.input, hname)
            except RuntimeError:
                print(f"Warning: could not find {hname}, skipping.")
                continue

            h = scale_to_per_2gev(h)
            style_hist(h, color, width=1, fill=True)
            flavor_hists[flavor] = h

    c = ROOT.TCanvas("c", "c", 900, 700)
    c.SetLeftMargin(0.13)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.07)

    if args.logy:
        c.SetLogy()

    xmin = 0.0
    xmax = args.xmax

    # Only use visible bins for ymax.
    ymax = 0.0
    for i in range(1, h_total.GetNbinsX() + 1):
        x1 = h_total.GetBinLowEdge(i)
        x2 = x1 + h_total.GetBinWidth(i)

        if x2 <= xmin or x1 >= xmax:
            continue

        ymax = max(ymax, h_total.GetBinContent(i))

    for h in flavor_hists.values():
        ymax = max(ymax, h.GetMaximum())

    paper_graph = None
    if args.overlay_paper_data:
        paper_graph = make_paper_data_graph()
        for i in range(paper_graph.GetN()):
            y = paper_graph.GetPointY(i)
            ymax = max(ymax, y)

    if args.logy:
        h_total.SetMinimum(0.01)
        h_total.SetMaximum(ymax * 20.0)
    else:
        h_total.SetMinimum(0.0)
        h_total.SetMaximum(ymax * 1.35)

    frame = c.DrawFrame(
        xmin,
        h_total.GetMinimum(),
        xmax,
        h_total.GetMaximum(),
        ";True lepton energy E_{l}^{true} [GeV];Predicted events / 2 GeV"
    )

    stack = ROOT.THStack("stack_components", "")

    # Draw largest/most dominant components first or choose your preferred order.
    for flavor in ["numu", "nue", "numubar", "nuebar"]:
        if flavor in flavor_hists:
            stack.Add(flavor_hists[flavor])

    stack.Draw("hist same")

    # Draw total prediction as black outline on top.
    h_total.Draw("hist same")

    if paper_graph:
        paper_graph.Draw("P same")

    leg = ROOT.TLegend(0.55, 0.62, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)

    leg.AddEntry(h_total, "Total prediction", "l")

    flavor_labels = {
        "numu": "#nu_{#mu} e",
        "nue": "#nu_{e} e",
        "numubar": "#bar{#nu}_{#mu} e",
        "nuebar": "#bar{#nu}_{e} e",
    }

    for flavor in ["numu", "nue", "numubar", "nuebar"]:
        if flavor in flavor_hists:
            leg.AddEntry(flavor_hists[flavor], flavor_labels[flavor], "f")

    if paper_graph:
        leg.AddEntry(paper_graph, "Paper data, acc.-corrected", "p")

    leg.Draw()

    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextSize(0.038)
    label.DrawLatex(0.16, 0.86, "#nu e^{-} #rightarrow #nu e^{-}")
    label.DrawLatex(0.16, 0.81, "Truth-level prediction")
    label.DrawLatex(0.16, 0.76, "Bin contents scaled to events / 2 GeV")

    c.SaveAs(args.output)

    print(f"Wrote plot to {args.output}")


if __name__ == "__main__":
    main()
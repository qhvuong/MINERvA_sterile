#!/usr/bin/env python3

import argparse
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


def get_cv_hist(filename, histname, name):
    f = ROOT.TFile.Open(filename)

    if not f or f.IsZombie():
        raise RuntimeError(
            f"Could not open {filename}"
        )

    h = f.Get(histname)

    if not h:
        raise RuntimeError(
            f"Could not find {histname} in {filename}"
        )

    # Works for both TH1D and MnvH1D.
    if hasattr(h, "GetCVHistoWithStatError"):
        out = h.GetCVHistoWithStatError()
        out = out.Clone(name)
    else:
        out = h.Clone(name)

    out.SetDirectory(0)
    f.Close()

    return out


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--lo", required=True)
    parser.add_argument("--ho", required=True)
    parser.add_argument(
        "--hist",
        default="h_nue_elastic_total",
    )
    parser.add_argument(
        "--output",
        default="nue_elastic_LO_vs_Tomalak.png",
    )

    args = parser.parse_args()

    h_lo = get_cv_hist(
        args.lo,
        args.hist,
        "h_lo",
    )

    h_ho = get_cv_hist(
        args.ho,
        args.hist,
        "h_ho",
    )

    if h_lo.GetNbinsX() != h_ho.GetNbinsX():
        raise RuntimeError(
            "LO and HO histograms have different bin counts"
        )

    for i in range(1, h_lo.GetNbinsX() + 1):
        if (
            abs(
                h_lo.GetBinLowEdge(i)
                - h_ho.GetBinLowEdge(i)
            )
            > 1.0e-9
            or abs(
                h_lo.GetBinWidth(i)
                - h_ho.GetBinWidth(i)
            )
            > 1.0e-9
        ):
            raise RuntimeError(
                f"Bin mismatch at bin {i}"
            )

    print()
    print(
        f"LO total      = {h_lo.Integral():.8f}"
    )
    print(
        f"Nominal total = {h_ho.Integral():.8f}"
    )

    if h_lo.Integral() != 0:
        diff = (
            h_ho.Integral()
            / h_lo.Integral()
            - 1.0
        ) * 100.0

        print(
            f"Total change  = {diff:+.4f}%"
        )

    print()
    print(
        "bin   range [GeV]       LO"
        "          nominal       nominal/LO"
    )

    for i in range(1, h_lo.GetNbinsX() + 1):
        low = h_lo.GetBinLowEdge(i)
        high = low + h_lo.GetBinWidth(i)

        lo = h_lo.GetBinContent(i)
        ho = h_ho.GetBinContent(i)

        ratio = ho / lo if lo != 0 else 0.0

        print(
            f"{i:2d}   "
            f"{low:5.1f}-{high:<5.1f}   "
            f"{lo:10.5f}   "
            f"{ho:10.5f}   "
            f"{ratio:10.5f}"
        )

    # ----------------------------------------------------------
    # Plot
    # ----------------------------------------------------------

    canvas = ROOT.TCanvas(
        "c",
        "c",
        800,
        800,
    )

    pad1 = ROOT.TPad(
        "pad1",
        "",
        0.0,
        0.30,
        1.0,
        1.0,
    )

    pad2 = ROOT.TPad(
        "pad2",
        "",
        0.0,
        0.0,
        1.0,
        0.30,
    )

    pad1.SetBottomMargin(0.02)
    pad2.SetTopMargin(0.03)
    pad2.SetBottomMargin(0.34)

    pad1.Draw()
    pad2.Draw()

    # ----------------------------------------------------------
    # Upper panel
    # ----------------------------------------------------------

    pad1.cd()

    h_lo.SetLineColor(ROOT.kBlue + 1)
    h_lo.SetLineWidth(3)
    h_lo.SetLineStyle(2)

    h_ho.SetLineColor(ROOT.kRed + 1)
    h_ho.SetLineWidth(3)

    ymax = max(
        h_lo.GetMaximum(),
        h_ho.GetMaximum(),
    )

    h_lo.SetMaximum(1.30 * ymax)
    h_lo.SetMinimum(0.0)

    h_lo.GetYaxis().SetTitle(
        "Predicted #nu+e events"
    )

    h_lo.GetXaxis().SetLabelSize(0)
    h_lo.GetXaxis().SetTitleSize(0)

    h_lo.Draw("HIST")
    h_ho.Draw("HIST SAME")

    legend = ROOT.TLegend(
        0.50,
        0.70,
        0.88,
        0.87,
    )

    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    legend.AddEntry(
        h_lo,
        "Leading order (tree level)",
        "l",
    )

    legend.AddEntry(
        h_ho,
        "Nominal Tomalak-based",
        "l",
    )

    legend.Draw()

    # ----------------------------------------------------------
    # Lower panel: fractional difference
    # ----------------------------------------------------------

    pad2.cd()

    h_frac = h_ho.Clone(
        "h_frac"
    )

    h_frac.Reset()

    for i in range(
        1,
        h_frac.GetNbinsX() + 1,
    ):
        lo = h_lo.GetBinContent(i)
        ho = h_ho.GetBinContent(i)

        if lo != 0.0:
            frac = (
                ho - lo
            ) / lo
        else:
            frac = 0.0

        h_frac.SetBinContent(
            i,
            frac,
        )

    h_frac.SetLineColor(
        ROOT.kBlack
    )

    h_frac.SetLineWidth(2)

    h_frac.GetYaxis().SetTitle(
        "(Nominal-LO)/LO"
    )

    h_frac.GetXaxis().SetTitle(
        "True electron energy [GeV]"
    )

    h_frac.GetYaxis().SetNdivisions(
        505
    )

    h_frac.GetYaxis().SetTitleSize(
        0.10
    )

    h_frac.GetYaxis().SetLabelSize(
        0.08
    )

    h_frac.GetYaxis().SetTitleOffset(
        0.45
    )

    h_frac.GetXaxis().SetTitleSize(
        0.12
    )

    h_frac.GetXaxis().SetLabelSize(
        0.09
    )

    h_frac.GetXaxis().SetTitleOffset(
        1.05
    )

    max_abs = 0.0

    for i in range(
        1,
        h_frac.GetNbinsX() + 1,
    ):
        max_abs = max(
            max_abs,
            abs(
                h_frac.GetBinContent(i)
            ),
        )

    yrange = max(
        0.02,
        1.25 * max_abs,
    )

    h_frac.SetMinimum(
        -yrange
    )

    h_frac.SetMaximum(
        yrange
    )

    h_frac.Draw("HIST")

    line = ROOT.TLine(
        h_frac.GetXaxis().GetXmin(),
        0.0,
        h_frac.GetXaxis().GetXmax(),
        0.0,
    )

    line.SetLineStyle(2)
    line.Draw()

    canvas.SaveAs(
        args.output
    )

    print()
    print(
        f"Wrote {args.output}"
    )


if __name__ == "__main__":
    main()
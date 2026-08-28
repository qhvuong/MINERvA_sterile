#!/usr/bin/env python3

from array import array
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

OLD_FILE = (
    "/exp/minerva/data/users/qvuong/nueel_prediction_studies/"
    "nue_elastic_prediction_higherOrderXS_mnv_FHC.root"
)

NEW_FILE = (
    "/exp/minerva/data/users/qvuong/nueel_prediction_studies/"
    "nue_elastic_prediction_TomalakXS_mnv_FHC.root"
)

OUTPUT = (
    "/exp/minerva/data/users/qvuong/nueel_prediction_studies/"
    "compare_Valencia_vs_Tomalak_FHC.png"
)

HIST_NAME = "h_nue_elastic_total"


def get_hist(path, name):
    f = ROOT.TFile.Open(path)

    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open {path}")

    h_mnv = f.Get(name)

    if not h_mnv:
        f.Close()
        raise RuntimeError(f"Could not find {name} in {path}")

    # Copy only the CV histogram into a plain TH1D.
    # This avoids MnvH1D / error-band ownership issues when closing the file.
    edges = []

    for i in range(1, h_mnv.GetNbinsX() + 1):
        edges.append(h_mnv.GetBinLowEdge(i))

    edges.append(
        h_mnv.GetBinLowEdge(h_mnv.GetNbinsX())
        + h_mnv.GetBinWidth(h_mnv.GetNbinsX())
    )

    from array import array

    h = ROOT.TH1D(
        name + "_plain",
        h_mnv.GetTitle(),
        h_mnv.GetNbinsX(),
        array("d", edges),
    )

    h.SetDirectory(0)

    for i in range(0, h_mnv.GetNbinsX() + 2):
        h.SetBinContent(i, h_mnv.GetBinContent(i))
        h.SetBinError(i, h_mnv.GetBinError(i))

    f.Close()

    return h


def main():

    h_old = get_hist(
        OLD_FILE,
        HIST_NAME,
    )

    h_new = get_hist(
        NEW_FILE,
        HIST_NAME,
    )

    if h_old.GetNbinsX() != h_new.GetNbinsX():
        raise RuntimeError("Histogram bin-count mismatch")

    for i in range(1, h_old.GetNbinsX() + 2):
        x_old = h_old.GetBinLowEdge(i)
        x_new = h_new.GetBinLowEdge(i)

        if abs(x_old - x_new) > 1e-12:
            raise RuntimeError(
                f"Bin-edge mismatch at edge {i}: "
                f"old={x_old}, new={x_new}"
            )

    print()
    print("Valencia-style vs Tomalak comparison")
    print()

    print(
        f"{'bin':>3s} "
        f"{'range [GeV]':>16s} "
        f"{'old':>12s} "
        f"{'Tomalak':>12s} "
        f"{'diff':>12s} "
        f"{'% diff':>10s}"
    )

    print("-" * 72)

    for i in range(1, h_old.GetNbinsX() + 1):

        xlow = h_old.GetBinLowEdge(i)
        xhigh = xlow + h_old.GetBinWidth(i)

        old = h_old.GetBinContent(i)
        new = h_new.GetBinContent(i)

        diff = new - old

        if old != 0.0:
            frac = 100.0 * diff / old
        else:
            frac = 0.0

        print(
            f"{i:3d} "
            f"{xlow:6.1f}-{xhigh:<6.1f} "
            f"{old:12.6f} "
            f"{new:12.6f} "
            f"{diff:12.6f} "
            f"{frac:9.4f}%"
        )

    old_total = h_old.Integral()
    new_total = h_new.Integral()

    print()
    print(f"Old total     = {old_total:.8f}")
    print(f"Tomalak total = {new_total:.8f}")
    print(f"Difference    = {new_total - old_total:+.8f}")
    print(
        f"Fractional    = "
        f"{100.0 * (new_total - old_total) / old_total:+.6f}%"
    )
    print()

    # ----------------------------------------------------------
    # Ratio histogram
    # ----------------------------------------------------------

    h_ratio = h_new.Clone("h_ratio")
    h_ratio.SetDirectory(0)

    for i in range(1, h_ratio.GetNbinsX() + 1):

        old = h_old.GetBinContent(i)
        new = h_new.GetBinContent(i)

        if old != 0.0:
            h_ratio.SetBinContent(i, new / old)
        else:
            h_ratio.SetBinContent(i, 0.0)

        h_ratio.SetBinError(i, 0.0)

    # ----------------------------------------------------------
    # Canvas
    # ----------------------------------------------------------

    c = ROOT.TCanvas(
        "c",
        "c",
        900,
        850,
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
    pad1.SetLeftMargin(0.13)
    pad1.SetRightMargin(0.05)
    pad1.SetTopMargin(0.07)

    pad2.SetTopMargin(0.03)
    pad2.SetBottomMargin(0.32)
    pad2.SetLeftMargin(0.13)
    pad2.SetRightMargin(0.05)

    pad1.Draw()
    pad2.Draw()

    # ----------------------------------------------------------
    # Top panel
    # ----------------------------------------------------------

    pad1.cd()

    h_old.SetTitle("")
    h_old.GetYaxis().SetTitle("Predicted events")
    h_old.GetYaxis().SetTitleSize(0.055)
    h_old.GetYaxis().SetLabelSize(0.045)
    h_old.GetXaxis().SetLabelSize(0.0)

    h_old.SetLineColor(ROOT.kBlue + 1)
    h_old.SetLineWidth(3)

    h_new.SetLineColor(ROOT.kRed + 1)
    h_new.SetLineWidth(3)
    h_new.SetLineStyle(2)

    ymax = max(
        h_old.GetMaximum(),
        h_new.GetMaximum(),
    )

    h_old.SetMaximum(1.30 * ymax)
    h_old.SetMinimum(0.0)

    h_old.Draw("hist")
    h_new.Draw("hist same")

    leg = ROOT.TLegend(
        0.55,
        0.72,
        0.88,
        0.88,
    )

    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)

    leg.AddEntry(
        h_old,
        "Previous Valencia higher-order",
        "l",
    )

    leg.AddEntry(
        h_new,
        "Tomalak",
        "l",
    )

    leg.Draw()

    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextSize(0.040)

    label.DrawLatex(
        0.16,
        0.86,
        "#nu e^{-} #rightarrow #nu e^{-}",
    )

    label.DrawLatex(
        0.16,
        0.80,
        "Full FHC prediction",
    )

    # ----------------------------------------------------------
    # Bottom ratio panel
    # ----------------------------------------------------------

    pad2.cd()

    h_ratio.SetTitle("")

    h_ratio.GetXaxis().SetTitle(
        "True lepton energy E_{l}^{true} [GeV]"
    )

    h_ratio.GetYaxis().SetTitle(
        "Tomalak / Valencia"
    )

    h_ratio.GetXaxis().SetTitleSize(0.12)
    h_ratio.GetXaxis().SetLabelSize(0.10)

    h_ratio.GetYaxis().SetTitleSize(0.10)
    h_ratio.GetYaxis().SetLabelSize(0.085)
    h_ratio.GetYaxis().SetTitleOffset(0.55)

    h_ratio.SetLineColor(ROOT.kBlack)
    h_ratio.SetLineWidth(2)
    h_ratio.SetMarkerStyle(20)
    h_ratio.SetMarkerSize(0.9)

    h_ratio.SetMinimum(0.95)
    h_ratio.SetMaximum(1.05)

    h_ratio.Draw("hist p")

    line = ROOT.TLine(
        h_ratio.GetXaxis().GetXmin(),
        1.0,
        h_ratio.GetXaxis().GetXmax(),
        1.0,
    )

    line.SetLineStyle(2)
    line.Draw()

    c.SaveAs(OUTPUT)

    print(f"Wrote plot to {OUTPUT}")


if __name__ == "__main__":
    main()
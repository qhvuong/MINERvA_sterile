#!/usr/bin/env python3
import ROOT
import sys
import math

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

try:
    import PlotUtils.LoadPlotUtilsLib
except Exception:
    pass


def to_list(root_vec):
    return [str(root_vec[i]) for i in range(root_vec.size())]


def make_ratio_hist(univ_hist, cv_hist, name, cv_min):
    ratio = univ_hist.Clone(name)
    ratio.SetDirectory(0)
    ratio.Reset()

    for b in range(1, cv_hist.GetNbinsX() + 1):
        cv = cv_hist.GetBinContent(b)
        u = univ_hist.GetBinContent(b)

        if cv > cv_min:
            ratio.SetBinContent(b, u / cv)
        else:
            ratio.SetBinContent(b, 0.0)

        ratio.SetBinError(b, 0.0)

    return ratio


def make_mean_and_rms_hists(ratio_hists, cv_hist, cv_min):
    mean_hist = cv_hist.Clone("mean_ratio")
    mean_hist.SetDirectory(0)
    mean_hist.Reset()

    rms_hist = cv_hist.Clone("rms_ratio")
    rms_hist.SetDirectory(0)
    rms_hist.Reset()

    for b in range(1, cv_hist.GetNbinsX() + 1):
        cv = cv_hist.GetBinContent(b)
        if cv <= cv_min:
            mean_hist.SetBinContent(b, 0.0)
            rms_hist.SetBinContent(b, 0.0)
            continue

        vals = [h.GetBinContent(b) for h in ratio_hists if h.GetBinContent(b) != 0.0]
        if len(vals) == 0:
            mean_hist.SetBinContent(b, 0.0)
            rms_hist.SetBinContent(b, 0.0)
            continue

        mean = sum(vals) / float(len(vals))
        var = sum((x - mean) ** 2 for x in vals) / float(len(vals))
        rms = math.sqrt(var)

        mean_hist.SetBinContent(b, mean)
        mean_hist.SetBinError(b, 0.0)

        rms_hist.SetBinContent(b, rms)
        rms_hist.SetBinError(b, 0.0)

    return mean_hist, rms_hist


def plot_universe_ratios(
    infile="testFlux.root",
    histname="flux_E_cvweighted",
    band="Flux",
    outname_lines="flux_ratio_lines.png",
    outname_mean="flux_ratio_mean.png",
    max_universes=None,
    xmin=0,
    xmax=20,
    ymin=0.7,
    ymax=1.3,
    cv_floor_frac=1e-3,
):
    f = ROOT.TFile.Open(infile)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open file: {infile}")

    h = f.Get(histname)
    if not h:
        f.ls()
        raise RuntimeError(f"Could not find histogram: {histname}")

    print(f"Loaded {histname}")
    print("Class:", h.ClassName())

    band_names = to_list(h.GetVertErrorBandNames())
    print("Vertical error bands:", band_names)

    if band not in band_names:
        raise RuntimeError(f"Band '{band}' not found. Available: {band_names}")

    eb = h.GetVertErrorBand(band)
    nuniv_total = eb.GetNHists()
    nuniv = nuniv_total if max_universes is None else min(max_universes, nuniv_total)
    print(f"Using band {band} with {nuniv}/{nuniv_total} universes")

    cv = h.Clone(f"{histname}_cv_clone")
    cv.SetDirectory(0)

    peak_cv = max(cv.GetBinContent(b) for b in range(1, cv.GetNbinsX() + 1))
    cv_min = cv_floor_frac * peak_cv
    print(f"Using CV threshold = {cv_min:.6g} ({cv_floor_frac} x peak CV)")

    ratio_hists = []
    true_ymin = 999.0
    true_ymax = -999.0

    for i in range(nuniv):
        hu = eb.GetHist(i).Clone(f"{histname}_{band}_univ_{i}_clone")
        hu.SetDirectory(0)

        ratio = make_ratio_hist(hu, cv, f"ratio_{band}_{i}", cv_min)
        ratio_hists.append(ratio)

        for b in range(1, ratio.GetNbinsX() + 1):
            val = ratio.GetBinContent(b)
            if val != 0.0:
                true_ymin = min(true_ymin, val)
                true_ymax = max(true_ymax, val)

    # -------------------------
    # Plot all universe/CV lines
    # -------------------------
    c1 = ROOT.TCanvas("c1", "c1", 1100, 800)
    c1.SetLeftMargin(0.13)
    c1.SetRightMargin(0.05)
    c1.SetBottomMargin(0.12)

    frame = cv.Clone("frame")
    frame.Reset()
    frame.SetDirectory(0)
    frame.SetStats(False)
    frame.SetTitle(f"{histname}: {band} universes / CV;E_{{#nu}} [GeV];Universe / CV")
    frame.GetXaxis().SetRangeUser(xmin, xmax)
    frame.GetYaxis().SetRangeUser(ymin, ymax)
    frame.Draw("AXIS")

    line1 = ROOT.TLine(xmin, 1.0, xmax, 1.0)
    line1.SetLineColor(ROOT.kBlack)
    line1.SetLineStyle(2)
    line1.SetLineWidth(2)
    line1.Draw()

    for r in ratio_hists:
        r.SetLineColorAlpha(ROOT.kBlue, 0.15)
        r.SetLineWidth(1)
        r.Draw("HIST SAME")

    frame.Draw("AXIS SAME")

    leg1 = ROOT.TLegend(0.55, 0.78, 0.90, 0.90)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    if ratio_hists:
        leg1.AddEntry(ratio_hists[0], f"{band} universes / CV", "l")
    leg1.AddEntry(line1, "ratio = 1", "l")
    leg1.Draw()

    c1.Print(outname_lines)
    print(f"Saved {outname_lines}")

    # -------------------------
    # Plot mean ratio and RMS
    # -------------------------
    mean_hist, rms_hist = make_mean_and_rms_hists(ratio_hists, cv, cv_min)

    c2 = ROOT.TCanvas("c2", "c2", 1100, 800)
    c2.SetLeftMargin(0.13)
    c2.SetRightMargin(0.05)
    c2.SetBottomMargin(0.12)

    frame2 = cv.Clone("frame2")
    frame2.Reset()
    frame2.SetDirectory(0)
    frame2.SetStats(False)
    frame2.SetTitle(f"{histname}: mean({band}/CV) and RMS;E_{{#nu}} [GeV];Ratio")
    frame2.GetXaxis().SetRangeUser(xmin, xmax)
    frame2.GetYaxis().SetRangeUser(ymin, ymax)
    frame2.Draw("AXIS")

    line2 = ROOT.TLine(xmin, 1.0, xmax, 1.0)
    line2.SetLineColor(ROOT.kBlack)
    line2.SetLineStyle(2)
    line2.SetLineWidth(2)
    line2.Draw()

    mean_hist.SetLineColor(ROOT.kRed + 1)
    mean_hist.SetLineWidth(3)
    mean_hist.Draw("HIST SAME")

    # Draw mean +/- RMS as thin lines
    plus = mean_hist.Clone("plus")
    minus = mean_hist.Clone("minus")
    plus.SetDirectory(0)
    minus.SetDirectory(0)

    for b in range(1, mean_hist.GetNbinsX() + 1):
        m = mean_hist.GetBinContent(b)
        s = rms_hist.GetBinContent(b)
        if m != 0.0:
            plus.SetBinContent(b, m + s)
            minus.SetBinContent(b, m - s)
        else:
            plus.SetBinContent(b, 0.0)
            minus.SetBinContent(b, 0.0)

    plus.SetLineColor(ROOT.kRed + 1)
    minus.SetLineColor(ROOT.kRed + 1)
    plus.SetLineStyle(2)
    minus.SetLineStyle(2)
    plus.SetLineWidth(2)
    minus.SetLineWidth(2)
    plus.Draw("HIST SAME")
    minus.Draw("HIST SAME")

    frame2.Draw("AXIS SAME")

    leg2 = ROOT.TLegend(0.52, 0.75, 0.90, 0.90)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.AddEntry(mean_hist, "Mean(universe / CV)", "l")
    leg2.AddEntry(plus, "Mean #pm RMS", "l")
    leg2.AddEntry(line2, "ratio = 1", "l")
    leg2.Draw()

    c2.Print(outname_mean)
    print(f"Saved {outname_mean}")

    print(f"Observed ratio range in kept bins: [{true_ymin:.3f}, {true_ymax:.3f}]")

    # do NOT call f.Close() -- ROOT/cppyy ownership can segfault here


if __name__ == "__main__":
    infile = sys.argv[1] if len(sys.argv) > 1 else "testFlux.root"

    # Combined
    plot_universe_ratios(
        infile=infile,
        histname="flux_E_cvweighted",
        band="Flux",
        outname_lines="flux_cvweighted_Flux_ratio_lines.png",
        outname_mean="flux_cvweighted_Flux_ratio_mean.png",
        max_universes=None,
        xmin=0,
        xmax=20,
        ymin=0.2,
        ymax=1.8,
        cv_floor_frac=1e-3,
    )

    # PPFX only
    plot_universe_ratios(
        infile=infile,
        histname="flux_E_cvweighted",
        band="ppfx1_Total",
        outname_lines="flux_cvweighted_ppfx_ratio_lines.png",
        outname_mean="flux_cvweighted_ppfx_ratio_mean.png",
        max_universes=None,
        xmin=0,
        xmax=20,
        ymin=0.2,
        ymax=1.8,
        cv_floor_frac=1e-3,
    )

    # Beam focus only
    plot_universe_ratios(
        infile=infile,
        histname="flux_E_cvweighted",
        band="Flux_BeamFocus",
        outname_lines="flux_cvweighted_beamfocus_ratio_lines.png",
        outname_mean="flux_cvweighted_beamfocus_ratio_mean.png",
        max_universes=None,
        xmin=0,
        xmax=20,
        ymin=0.2,
        ymax=1.8,
        cv_floor_frac=1e-3,
    )
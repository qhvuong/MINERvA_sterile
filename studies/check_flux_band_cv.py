#!/usr/bin/env python

import os
import sys
import ROOT

ROOT.TH1.AddDirectory(False)

DEFAULT_FILE = "/exp/minerva/data/users/qvuong/nu_e/kin_dist_mcleFHC_CCnue_allSystematics_testCVs_MAD.root"

# Start with a few useful histograms. You can add more.
DEFAULT_HISTS = [
    "EN4",
    "EN4_CCNuEQE",
    "EN4_CCNuEDelta",
    "EN4_NCPi0",
    "EN4_dEdX",
    "EN4_dEdX_CCNuEQE",
    "EN4_dEdX_CCNuEDelta",
    "EN4_dEdX_NCPi0",
]

OUTDIR = "band_cv_debug_plots"


def safe_name(s):
    return (
        str(s)
        .replace("/", "_")
        .replace(" ", "_")
        .replace("(", "")
        .replace(")", "")
        .replace("+", "p")
        .replace("-", "m")
    )


def get_band_names(h):
    # Mostly vertical bands here, but this keeps it explicit.
    try:
        return [str(x) for x in h.GetVertErrorBandNames()]
    except Exception:
        return [str(x) for x in h.GetErrorBandNames() if h.HasVertErrorBand(str(x))]


def print_band_cv_check(h, bandname, max_bad_bins=20):
    band = h.GetVertErrorBand(bandname)

    main_int = h.Integral()
    band_int = band.Integral()
    ratio = band_int / main_int if main_int else "NA"

    print("\n  Band:", bandname)
    print("    main integral   =", main_int)
    print("    bandCV integral =", band_int)
    print("    bandCV/main     =", ratio)
    print("    n universes     =", band.GetNHists())

    nbad = 0
    max_abs_diff = 0.0

    for b in range(0, h.GetNbinsX() + 2):
        main = h.GetBinContent(b)
        bcv = band.GetBinContent(b)
        diff = bcv - main
        max_abs_diff = max(max_abs_diff, abs(diff))

        if abs(diff) > 1e-10:
            nbad += 1
            if nbad <= max_bad_bins:
                print(
                    "      DIFF bin {:2d}: main={:.12g} bandCV={:.12g} diff={:.12g} ratio={}".format(
                        b,
                        main,
                        bcv,
                        diff,
                        bcv / main if main else "NA",
                    )
                )

    if nbad == 0:
        print("    bin-by-bin check: OK")
    else:
        print("    bin-by-bin check: {} differing bins, max_abs_diff={:.12g}".format(nbad, max_abs_diff))


def draw_band_cv_and_universes(h, bandname, outdir):
    if not h.HasVertErrorBand(bandname):
        return

    band = h.GetVertErrorBand(bandname)

    h_main = ROOT.TH1D(h)
    h_main.SetDirectory(0)
    h_main.SetStats(0)
    h_main.SetName("{}_{}_mainCV".format(h.GetName(), safe_name(bandname)))
    h_main.SetLineColor(ROOT.kBlack)
    h_main.SetLineWidth(4)
    h_main.SetLineStyle(2)  # dashed

    h_bandcv = ROOT.TH1D(band)
    h_bandcv.SetDirectory(0)
    h_bandcv.SetStats(0)
    h_bandcv.SetName("{}_{}_bandCV".format(h.GetName(), safe_name(bandname)))
    h_bandcv.SetLineColor(ROOT.kBlue)
    h_bandcv.SetLineWidth(2)
    h_bandcv.SetMarkerColor(ROOT.kBlue)
    h_bandcv.SetMarkerStyle(20)
    h_bandcv.SetMarkerSize(1.0)

    univs = []
    ymax = max(h_main.GetMaximum(), h_bandcv.GetMaximum())

    for i in range(band.GetNHists()):
        hu = ROOT.TH1D(band.GetHist(i))
        hu.SetDirectory(0)
        hu.SetStats(0)
        hu.SetName("{}_{}_univ{}".format(h.GetName(), safe_name(bandname), i))
        hu.SetLineColor(ROOT.kGray + 1)
        hu.SetLineWidth(1)
        univs.append(hu)
        ymax = max(ymax, hu.GetMaximum())

    cname = "c_{}_{}".format(h.GetName(), safe_name(bandname))
    c = ROOT.TCanvas(cname, cname, 1200, 900)

    h_main.SetMaximum(1.25 * ymax if ymax > 0 else 1.0)
    h_main.SetMinimum(0.0)
    h_main.Draw("HIST")

    for hu in univs:
        hu.Draw("HIST SAME")

    # band CV as blue line + markers
    h_bandcv.Draw("HIST SAME")
    h_bandcv.Draw("P SAME")

    # redraw main dashed line on top
    h_main.Draw("HIST SAME")

    leg = ROOT.TLegend(0.55, 0.68, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(h_main, "Main CV (dashed)", "l")
    leg.AddEntry(h_bandcv, "Band CV (blue markers)", "lp")
    if univs:
        leg.AddEntry(univs[0], "Universes", "l")
    leg.Draw()

    c.Modified()
    c.Update()

    # Keep PyROOT references alive until SaveAs
    c._h_main = h_main
    c._h_bandcv = h_bandcv
    c._univs = univs
    c._leg = leg

    outname = os.path.join(
        outdir,
        "{}__{}.png".format(safe_name(h.GetName()), safe_name(bandname))
    )
    print("    saving:", outname)
    c.SaveAs(outname)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_FILE

    os.makedirs(OUTDIR, exist_ok=True)

    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError("Could not open file: {}".format(path))

    print("Opened:", path)

    # Use requested hist list, but only keep ones present.
    hist_names = []
    for hname in DEFAULT_HISTS:
        if f.Get(hname):
            hist_names.append(hname)

    print("\nHistograms to check:")
    for hname in hist_names:
        print(" ", hname)

    for hname in hist_names:
        h = f.Get(hname)
        if not h:
            continue

        print("\n" + "=" * 100)
        print("Histogram:", hname)
        print("Title:", h.GetTitle())
        print("Main CV integral:", h.Integral())

        bands = get_band_names(h)
        print("Available vertical bands:")
        for b in bands:
            print(" ", b)

        for bandname in bands:
            if not h.HasVertErrorBand(bandname):
                continue

            print_band_cv_check(h, bandname)
            # draw_band_cv_and_universes(h, bandname, OUTDIR)

    f.Close()
    print("\nDone. Plots are in:", OUTDIR)


if __name__ == "__main__":
    main()
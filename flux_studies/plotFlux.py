#!/usr/bin/env python3
import ROOT
import sys

ROOT.gROOT.SetBatch(True)

# Load PlotUtils if needed
try:
    import PlotUtils.LoadPlotUtilsLib
except Exception:
    pass


def to_list(root_vec):
    return [str(root_vec[i]) for i in range(root_vec.size())]


def plot_flux_universes(
    infile="testFlux.root",
    histname="flux_E_cvweighted",
    band="Flux",
    outname=None,
    max_universes=None,
    xmin=None,
    xmax=None,
    logy=False,
):
    f = ROOT.TFile.Open(infile)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open file: {infile}")

    h = f.Get(histname)
    if not h:
        print("\nAvailable histograms:")
        f.ls()
        raise RuntimeError(f"Could not find histogram: {histname}")

    print(f"Loaded {histname}")
    print("Class:", h.ClassName())

    band_names = to_list(h.GetVertErrorBandNames())
    print("Vertical error bands:", band_names)

    if band not in band_names:
        raise RuntimeError(
            f"Could not find band '{band}'. Available bands are: {band_names}"
        )

    eb = h.GetVertErrorBand(band)
    nuniv_total = eb.GetNHists()
    nuniv = nuniv_total if max_universes is None else min(max_universes, nuniv_total)

    print(f"Using band {band} with {nuniv}/{nuniv_total} universes")

    # cv = h.GetCVHisto().Clone(f"{histname}_cv_clone")
    cv = h.Clone(f"{histname}_cv_clone")
    cv.SetDirectory(0)

    universes = []
    ymax = cv.GetMaximum()
    ymin_pos = None

    for i in range(nuniv):
        hu = eb.GetHist(i).Clone(f"{histname}_{band}_univ_{i}_clone")
        hu.SetDirectory(0)
        universes.append(hu)

        ymax = max(ymax, hu.GetMaximum())

        for b in range(1, hu.GetNbinsX() + 1):
            val = hu.GetBinContent(b)
            if val > 0:
                ymin_pos = val if ymin_pos is None else min(ymin_pos, val)

    if outname is None:
        outname = f"{histname}_{band}_universes.png"

    c = ROOT.TCanvas("c", "c", 1100, 800)
    c.SetLeftMargin(0.13)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(0.12)

    if logy:
        c.SetLogy()

    cv.SetStats(False)
    cv.SetTitle(f"{histname}: CV and {band} universes;E_{{#nu}} [GeV];Flux / m^{{2}} / POT / GeV")
    cv.SetLineColor(ROOT.kBlack)
    cv.SetLineWidth(3)

    if xmin is not None and xmax is not None:
        cv.GetXaxis().SetRangeUser(xmin, xmax)

    if logy:
        ymin = ymin_pos * 0.3 if ymin_pos is not None else 1e-30
        cv.GetYaxis().SetRangeUser(ymin, ymax * 3.0)
    else:
        cv.GetYaxis().SetRangeUser(0.0, ymax * 1.25)

    # Draw CV first to define axes
    cv.Draw("HIST")

    # Draw universes
    for hu in universes:
        hu.SetLineColorAlpha(ROOT.kBlue, 0.18)
        hu.SetLineWidth(1)
        hu.Draw("HIST SAME")

    # Draw CV again on top
    cv.Draw("HIST SAME")

    leg = ROOT.TLegend(0.58, 0.75, 0.90, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(cv, "CV", "l")
    if len(universes) > 0:
        leg.AddEntry(universes[0], f"{band} universes", "l")
    leg.Draw()

    c.Print(outname)
    print(f"Saved {outname}")

    # f.Close()


if __name__ == "__main__":
    infile = sys.argv[1] if len(sys.argv) > 1 else "testFlux.root"

    # Main one you probably want:
    plot_flux_universes(
        infile=infile,
        histname="flux_E_cvweighted",
        band="Flux",
        outname="flux_cvweighted_Flux_universes.png",
        max_universes=None,
        xmin=0,
        xmax=20,
        logy=False,
    )

    # Also useful for checking component bands separately:
    plot_flux_universes(infile, "flux_E_cvweighted", "ppfx1_Total",
                        "flux_cvweighted_ppfx1_Total_universes.png",
                        xmin=0, xmax=20)
    
    plot_flux_universes(infile, "flux_E_cvweighted", "Flux_BeamFocus",
                        "flux_cvweighted_Flux_BeamFocus_universes.png",
                        xmin=0, xmax=20)
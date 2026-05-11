#!/usr/bin/env python3
import ROOT
import sys

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

try:
    import PlotUtils.LoadPlotUtilsLib
except Exception:
    pass


def to_list(v):
    return [str(v[i]) for i in range(v.size())]


def check_flux_bin(infile, histname="flux_E_cvweighted", band="Flux", energy=None, bin_number=None):
    f = ROOT.TFile.Open(infile)
    if not f or f.IsZombie():
        raise RuntimeError("Could not open {}".format(infile))

    h = f.Get(histname)
    if not h:
        print("Contents of file:")
        f.ls()
        raise RuntimeError("Could not find histogram {}".format(histname))

    bands = to_list(h.GetVertErrorBandNames())

    print("\n" + "=" * 100)
    print("File:     {}".format(infile))
    print("Hist:     {}".format(histname))
    print("Class:    {}".format(h.ClassName()))
    print("Bands:    {}".format(bands))

    if band not in bands:
        print("Band {} not found; skipping.".format(band))
        return

    eb = h.GetVertErrorBand(band)
    nuniv = eb.GetNHists()

    if bin_number is not None:
        b = int(bin_number)
    elif energy is not None:
        b = h.GetXaxis().FindBin(float(energy))
    else:
        raise RuntimeError("Need either energy or bin_number")

    xlo = h.GetXaxis().GetBinLowEdge(b)
    xhi = h.GetXaxis().GetBinUpEdge(b)

    cv = h.GetBinContent(b)

    vals = []
    for i in range(nuniv):
        hu = eb.GetHist(i)
        u = hu.GetBinContent(b)
        r = u / cv if cv != 0 else -999
        vals.append((i, u, r))

    min_i, min_u, min_r = min(vals, key=lambda x: x[1])
    max_i, max_u, max_r = max(vals, key=lambda x: x[1])

    n_above = sum(1 for _, _, r in vals if r > 1.0)
    n_below = sum(1 for _, _, r in vals if r < 1.0)
    n_equal = nuniv - n_above - n_below

    print("\nBand: {}".format(band))
    print("Energy requested: {}".format(energy))
    print("Bin {}: [{:.6g}, {:.6g}] GeV".format(b, xlo, xhi))
    print("CV = {:.12g}".format(cv))
    print("min universe = {:.12g}  universe {}  ratio = {:.8g}".format(min_u, min_i, min_r))
    print("max universe = {:.12g}  universe {}  ratio = {:.8g}".format(max_u, max_i, max_r))
    print("universes above CV: {} / {}".format(n_above, nuniv))
    print("universes below CV: {} / {}".format(n_below, nuniv))
    print("universes equal CV: {} / {}".format(n_equal, nuniv))

    if cv < min_u:
        print("RESULT: CV is BELOW all universes in this flux true-energy bin.")
    elif cv > max_u:
        print("RESULT: CV is ABOVE all universes in this flux true-energy bin.")
    else:
        print("RESULT: CV is inside the universe envelope.")

    print("\nSorted universe/CV ratios:")
    print("  universe      ratio        content")
    for i, u, r in sorted(vals, key=lambda x: x[2]):
        print("  {:4d}   {:12.6g}   {:12.6g}".format(i, r, u))

    # Avoid f.Close(); PyROOT/cppyy sometimes segfaults with PlotUtils ownership.


def check_multiple(infile, histname, energies, bands):
    for energy in energies:
        for band in bands:
            check_flux_bin(
                infile=infile,
                histname=histname,
                band=band,
                energy=energy,
            )


if __name__ == "__main__":
    infile = sys.argv[1] if len(sys.argv) > 1 else "custom_plotutils/data/flux/flux-gen2thin-pdg12-minerva1.root"

    histname = "flux_E_cvweighted"

    # Energies chosen to probe the same true-Eν region as the high-EN4 selected events.
    # With the 38-bin LE binning:
    #   14 GeV -> [14, 15]
    #   15 GeV -> [15, 16]
    #   17.1217 GeV -> [17, 18]
    #   24 GeV -> [20, 30]
    energies = [14.0, 15.0, 17.1217, 24.0, 32.0, 43.0]

    bands = ["Flux", "ppfx1_Total", "Flux_BeamFocus"]

    check_multiple(infile, histname, energies, bands)
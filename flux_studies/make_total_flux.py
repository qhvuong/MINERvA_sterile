#!/usr/bin/env python3
import sys
import ROOT

# This is the same mechanism your compute_flux.py uses to get PlotUtils classes
import PlotUtils.LoadPlotUtilsLib

def make_total(infile, outfile):
    fin = ROOT.TFile.Open(infile, "READ")
    if not fin or fin.IsZombie():
        raise RuntimeError(f"Cannot open input file: {infile}")

    h_unw_total = None
    h_cv_total  = None

    keys = fin.GetListOfKeys()
    for k in keys:
        name = k.GetName()

        if name.startswith("flux_E_unweighted_parent"):
            h = fin.Get(name)
            if not h:
                continue
            if h_unw_total is None:
                h_unw_total = h.Clone("flux_E_unweighted")
                h_unw_total.Reset()
            h_unw_total.Add(h)

        if name.startswith("flux_E_cvweighted_parent"):
            h = fin.Get(name)
            if not h:
                continue
            if h_cv_total is None:
                h_cv_total = h.Clone("flux_E_cvweighted")
                h_cv_total.Reset()
            h_cv_total.Add(h)

    if h_unw_total is None or h_cv_total is None:
        raise RuntimeError("Did not find both unweighted and cvweighted parent flux histograms.")

    totalPOT = fin.Get("total_POT")  # TParameter<double>

    fout = ROOT.TFile.Open(outfile, "RECREATE")
    if not fout or fout.IsZombie():
        raise RuntimeError(f"Cannot open output file: {outfile}")

    # Write only what you want
    h_unw_total.Write()
    h_cv_total.Write()
    if totalPOT:
        totalPOT.Write()

    fout.Close()
    fin.Close()
    print(f"Wrote {outfile} with: flux_E_unweighted_total, flux_E_cvweighted_total, total_POT")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: make_total_flux.py input.root output.root")
        sys.exit(1)
    make_total(sys.argv[1], sys.argv[2])

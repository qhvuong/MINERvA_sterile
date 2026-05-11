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


def check_last_n_bins(infile, histname="EN4", band="Flux", n_last=2):
    f = ROOT.TFile.Open(infile)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open {infile}")

    h = f.Get(histname)
    if not h:
        f.ls()
        raise RuntimeError(f"Could not find {histname}")

    bands = to_list(h.GetVertErrorBandNames())
    print("\nHistogram:", histname)
    print("Bands:", bands)

    if band not in bands:
        raise RuntimeError(f"Could not find band {band}")

    eb = h.GetVertErrorBand(band)
    nuniv = eb.GetNHists()

    # MnvH1D itself is the main CV in your binding
    cv = h

    # Find all bins where CV or any universe has nonzero content
    populated_bins = []
    for b in range(1, h.GetNbinsX() + 1):
        cv_val = cv.GetBinContent(b)
        vals = [eb.GetHist(i).GetBinContent(b) for i in range(nuniv)]
        if cv_val != 0 or any(v != 0 for v in vals):
            populated_bins.append(b)

    if not populated_bins:
        print("No populated bins found.")
        return

    bins_to_check = populated_bins[-n_last:]

    print(f"\nChecking last {len(bins_to_check)} populated bins:", bins_to_check)

    for b in bins_to_check:
        xlo = h.GetXaxis().GetBinLowEdge(b)
        xhi = h.GetXaxis().GetBinUpEdge(b)

        cv_val = cv.GetBinContent(b)
        vals = [(i, eb.GetHist(i).GetBinContent(b)) for i in range(nuniv)]

        min_i, min_u = min(vals, key=lambda x: x[1])
        max_i, max_u = max(vals, key=lambda x: x[1])

        ratios = []
        for i, u in vals:
            r = u / cv_val if cv_val != 0 else float("nan")
            ratios.append((i, u, r))

        print("\n" + "-" * 80)
        print(f"Bin {b}: [{xlo}, {xhi}]")
        print(f"  CV = {cv_val:.12g}")
        print(f"  min universe = {min_u:.12g}  universe {min_i}")
        print(f"  max universe = {max_u:.12g}  universe {max_i}")

        if cv_val < min_u:
            print("  RESULT: CV is BELOW the smallest universe in this bin.")
        elif cv_val > max_u:
            print("  RESULT: CV is ABOVE the largest universe in this bin.")
        else:
            print("  RESULT: CV is inside the universe envelope.")

        print("\n  Universe ratios in this bin:")
        print("    universe   content        universe/CV")
        for i, u, r in ratios:
            print(f"    {i:4d}   {u:12.6g}   {r:12.6g}")

        print("\n  Sorted universe/CV ratios:")
        for i, u, r in sorted(ratios, key=lambda x: x[2]):
            print(f"    {i:4d}   ratio={r:12.6g}   content={u:12.6g}")


if __name__ == "__main__":
    infile = sys.argv[1]
    histname = sys.argv[2] if len(sys.argv) > 2 else "EN4"
    band = sys.argv[3] if len(sys.argv) > 3 else "Flux"
    n_last = int(sys.argv[4]) if len(sys.argv) > 4 else 2

    check_last_n_bins(infile, histname, band, n_last)
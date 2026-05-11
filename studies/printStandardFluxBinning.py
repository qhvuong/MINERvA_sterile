#!/usr/bin/env python3
import ROOT
import os

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

try:
    import PlotUtils.LoadPlotUtilsLib
except Exception:
    pass


def get_edges(path, hname="flux_E_cvweighted"):
    path = os.path.expandvars(path)
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open {path}")

    h = f.Get(hname)
    if not h:
        f.ls()
        raise RuntimeError(f"Could not find {hname} in {path}")

    ax = h.GetXaxis()
    edges = [ax.GetBinLowEdge(i) for i in range(1, h.GetNbinsX() + 1)]
    edges.append(ax.GetBinUpEdge(h.GetNbinsX()))
    return edges


numu_le_file = "custom_plotutils/data/flux/flux-gen2thin-pdg14-minerva1.root"
nue_le_file  = "custom_plotutils/data/flux/flux-gen2thin-pdg12-minerva1.root"

numu_edges = get_edges(numu_le_file)
nue_edges = get_edges(nue_le_file)

print("Standard LE numu flux binning:")
print("  file:", numu_le_file)
print("  Number of bins:", len(numu_edges) - 1)
print("  Bin edges:")

for i in range(len(numu_edges) - 1):
    print(f"bin {i+1:3d}: {numu_edges[i]:10.6g}  {numu_edges[i+1]:10.6g}")

print("\nAs Python list:")
print("STANDARD_LE_FLUX_BINNING = [")
for x in numu_edges:
    print(f"    {float(x):.10g},")
print("]")

if len(numu_edges) == len(nue_edges) and all(abs(a-b) < 1e-12 for a, b in zip(numu_edges, nue_edges)):
    print("\nnumu LE and nue LE files have identical binning.")
else:
    print("\nnumu LE and nue LE files do NOT have identical binning.")
    print("  numu bins:", len(numu_edges) - 1)
    print("  nue bins: ", len(nue_edges) - 1)
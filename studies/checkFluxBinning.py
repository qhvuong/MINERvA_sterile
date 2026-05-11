#!/usr/bin/env python3
import ROOT
import sys
import os
import math

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

try:
    import PlotUtils.LoadPlotUtilsLib
except Exception:
    pass


def get_hist(path, hname="flux_E_cvweighted"):
    path = os.path.expandvars(path)
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open {path}")

    h = f.Get(hname)
    if not h:
        print("Contents of", path)
        f.ls()
        raise RuntimeError(f"Could not find {hname} in {path}")

    h.SetDirectory(0)
    return f, h


def axis_edges(h):
    ax = h.GetXaxis()
    return [ax.GetBinLowEdge(i) for i in range(1, h.GetNbinsX() + 1)] + [
        ax.GetBinUpEdge(h.GetNbinsX())
    ]


def print_binning(label, h):
    ax = h.GetXaxis()
    edges = axis_edges(h)

    print("\n" + "=" * 80)
    print(label)
    print("  hist name:", h.GetName())
    print("  class:", h.ClassName())
    print("  nbins:", h.GetNbinsX())
    print("  xmin:", ax.GetXmin())
    print("  xmax:", ax.GetXmax())
    print("  first 10 edges:", ["{:.8g}".format(x) for x in edges[:10]])
    print("  last 10 edges: ", ["{:.8g}".format(x) for x in edges[-10:]])

    uniform = True
    if h.GetNbinsX() > 1:
        bw0 = edges[1] - edges[0]
        for i in range(1, len(edges) - 1):
            if abs((edges[i + 1] - edges[i]) - bw0) > 1e-10:
                uniform = False
                break
        print("  uniform binning:", uniform)
        if uniform:
            print("  bin width:", bw0)


def compare_axes(label_a, h_a, label_b, h_b, tol=1e-10):
    ea = axis_edges(h_a)
    eb = axis_edges(h_b)

    print("\n" + "-" * 80)
    print(f"Compare: {label_a}  vs  {label_b}")

    if h_a.GetNbinsX() != h_b.GetNbinsX():
        print("  RESULT: DIFFERENT number of bins")
        print("   ", label_a, "nbins =", h_a.GetNbinsX())
        print("   ", label_b, "nbins =", h_b.GetNbinsX())
        return False

    bad = []
    for i, (a, b) in enumerate(zip(ea, eb)):
        if abs(a - b) > tol:
            bad.append((i, a, b, a - b))

    if not bad:
        print("  RESULT: SAME bin edges within tolerance")
        return True

    print("  RESULT: DIFFERENT bin edges")
    print("  first few differences:")
    for i, a, b, d in bad[:20]:
        print(f"    edge {i:4d}: {label_a}={a:.12g}, {label_b}={b:.12g}, diff={d:.4g}")
    if len(bad) > 20:
        print(f"    ... and {len(bad)-20} more")
    return False


def main():
    if len(sys.argv) < 3:
        print("""
Usage:
  python checkFluxBinning.py standard1.root produced1.root [standard2.root produced2.root ...]

Example:
  python checkFluxBinning.py \\
    $PLOTUTILSROOT/data/flux/flux-gen2thin-pdg12-minerva1.root custom_plotutils/data/flux/flux-gen2thin-pdg12-minerva1.root \\
    $PLOTUTILSROOT/data/flux/flux-gen2thin-pdg14-minerva1.root custom_plotutils/data/flux/flux-gen2thin-pdg14-minerva1.root
""")
        sys.exit(1)

    hname = "flux_E_cvweighted"

    pairs = []
    args = sys.argv[1:]
    if len(args) % 2 != 0:
        raise RuntimeError("Pass files as pairs: standard produced")

    for i in range(0, len(args), 2):
        pairs.append((args[i], args[i + 1]))

    keep_files = []

    for standard_path, produced_path in pairs:
        fs, hs = get_hist(standard_path, hname)
        fp, hp = get_hist(produced_path, hname)
        keep_files.extend([fs, fp])

        standard_label = "STANDARD: " + standard_path
        produced_label = "PRODUCED: " + produced_path

        print_binning(standard_label, hs)
        print_binning(produced_label, hp)

        compare_axes(standard_label, hs, produced_label, hp)

    print("\nDone.")


if __name__ == "__main__":
    main()
#!/usr/bin/env python
import ROOT
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("--file1", required=True)
ap.add_argument("--file2", required=True)
ap.add_argument("--hist1", required=True)
ap.add_argument("--hist2", required=True)
ap.add_argument("--emin", type=float, default=0.0)
ap.add_argument("--emax", type=float, default=20.0)
args = ap.parse_args()

def get_cv(path, hname):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError("Could not open " + path)

    h = f.Get(hname)
    if not h:
        f.ls()
        raise RuntimeError("Could not find hist {} in {}".format(hname, path))

    if hasattr(h, "GetCVHistoWithStatError"):
        hcv = h.GetCVHistoWithStatError()
    elif hasattr(h, "GetCVHistoWithError"):
        hcv = h.GetCVHistoWithError()
    else:
        hcv = h

    hcv.SetDirectory(0)
    f.Close()
    return hcv

h1 = get_cv(args.file1, args.hist1)
h2 = get_cv(args.file2, args.hist2)

print("[COMPARE]")
print("file1:", args.file1)
print("hist1:", args.hist1)
print("file2:", args.file2)
print("hist2:", args.hist2)
print("{:>5} {:>10} {:>10} {:>16} {:>16} {:>16}".format(
    "bin", "lo", "hi", "file1_CV", "file2_CV", "file2/file1"
))

for b in range(1, min(h1.GetNbinsX(), h2.GetNbinsX()) + 1):
    lo1 = h1.GetBinLowEdge(b)
    hi1 = lo1 + h1.GetBinWidth(b)
    lo2 = h2.GetBinLowEdge(b)
    hi2 = lo2 + h2.GetBinWidth(b)

    if hi1 <= args.emin or lo1 >= args.emax:
        continue

    v1 = h1.GetBinContent(b)
    v2 = h2.GetBinContent(b)
    ratio = v2 / v1 if v1 else -999

    print("{:5d} {:10.4g} {:10.4g} {:16.8e} {:16.8e} {:16.8e}".format(
        b, lo1, hi1, v1, v2, ratio
    ))

    if abs(lo1 - lo2) > 1e-6 or abs(hi1 - hi2) > 1e-6:
        print("      [WARNING] binning mismatch: file1 [{}, {}], file2 [{}, {}]".format(
            lo1, hi1, lo2, hi2
        ))
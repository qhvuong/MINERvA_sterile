#!/usr/bin/env python
import ROOT
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("--p6", required=True)
ap.add_argument("--p8", required=True)
ap.add_argument("--hist", default="flux_E_cvweighted")
ap.add_argument("--emin", type=float, default=0.0)
ap.add_argument("--emax", type=float, default=20.0)
args = ap.parse_args()

def get_cv(path, hname):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError("Could not open " + path)

    h = f.Get(hname)
    if not h:
        raise RuntimeError("Could not find hist {} in {}".format(hname, path))

    # If MnvH1D, get CV hist. If already TH1, use directly.
    if hasattr(h, "GetCVHistoWithStatError"):
        hcv = h.GetCVHistoWithStatError()
    elif hasattr(h, "GetCVHistoWithError"):
        hcv = h.GetCVHistoWithError()
    else:
        hcv = h

    hcv.SetDirectory(0)
    f.Close()
    return hcv

h6 = get_cv(args.p6, args.hist)
h8 = get_cv(args.p8, args.hist)

print("[COMPARE]", args.hist)
print("{:>5} {:>10} {:>10} {:>16} {:>16} {:>16}".format(
    "bin", "lo", "hi", "P6_CV", "P8_CV", "P8/P6"
))

for b in range(1, h6.GetNbinsX() + 1):
    lo = h6.GetBinLowEdge(b)
    hi = lo + h6.GetBinWidth(b)
    if hi <= args.emin or lo >= args.emax:
        continue

    p6 = h6.GetBinContent(b)
    p8 = h8.GetBinContent(b)
    ratio = p8 / p6 if p6 else -999

    print("{:5d} {:10.4g} {:10.4g} {:16.8e} {:16.8e} {:16.8e}".format(
        b, lo, hi, p6, p8, ratio
    ))
#!/usr/bin/env python

import ROOT
import math
import sys

ROOT.TH1.AddDirectory(False)

TAG = "CCnue_allSystematics_testCVs"  # change if needed
BASE = "/exp/minerva/data/users/qvuong/nu_e"

INPUTS = {
    "le1":   f"{BASE}/kin_dist_mcle1_{TAG}_MAD.root",
    "le7":   f"{BASE}/kin_dist_mcle7_{TAG}_MAD.root",
    "le9":   f"{BASE}/kin_dist_mcle9_{TAG}_MAD.root",
    "le13C": f"{BASE}/kin_dist_mcle13C_{TAG}_MAD.root",
}

FHC = f"{BASE}/kin_dist_mcleFHC_{TAG}_MAD.root"

HISTS = ["EN4", "EN4_CCNuEQE", "EN4_CCNuEDelta", "EN4_NCPi0"]
BAND = "Flux"

def clone_hist(h, name):
    out = ROOT.TH1D(h)
    out.SetName(name)
    out.SetDirectory(0)
    return out

def get_hist(path, hname):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError("Cannot open " + path)
    h = f.Get(hname)
    if not h:
        raise RuntimeError("Missing {} in {}".format(hname, path))
    # detach by cloning MnvH1D object itself
    hc = h.Clone("{}_clone".format(hname))
    hc.SetDirectory(0)
    f.Close()
    return hc

def max_abs_diff(h1, h2):
    maxdiff = 0.0
    maxbin = -1
    for b in range(0, h1.GetNbinsX() + 2):
        d = abs(h1.GetBinContent(b) - h2.GetBinContent(b))
        if d > maxdiff:
            maxdiff = d
            maxbin = b
    return maxdiff, maxbin

def check_one(hname, bandname):
    print("\n" + "=" * 100)
    print("Checking", hname, bandname)

    # Sum parent CVs from individual playlist files
    summed_cv = None
    summed_bandcv = None
    summed_univs = []

    for label, path in INPUTS.items():
        h = get_hist(path, hname)

        if summed_cv is None:
            summed_cv = clone_hist(h, "{}_sum_cv".format(hname))
            summed_cv.Reset()

        summed_cv.Add(h)

        if not h.HasVertErrorBand(bandname):
            raise RuntimeError("{} has no band {}".format(path, bandname))

        band = h.GetVertErrorBand(bandname)

        if summed_bandcv is None:
            summed_bandcv = clone_hist(band, "{}_sum_bandcv".format(hname))
            summed_bandcv.Reset()

        summed_bandcv.Add(band)

        if not summed_univs:
            for i in range(band.GetNHists()):
                u = clone_hist(band.GetHist(i), "{}_sum_u{}".format(hname, i))
                u.Reset()
                summed_univs.append(u)

        for i in range(band.GetNHists()):
            summed_univs[i].Add(band.GetHist(i))

        print(
            "{:5s}: CV={:.8g} bandCV={:.8g} bandCV/CV={:.6g}".format(
                label,
                h.Integral(),
                band.Integral(),
                band.Integral() / h.Integral() if h.Integral() else 0.0,
            )
        )

    # Read FHC output
    hfhc = get_hist(FHC, hname)
    bfhc = hfhc.GetVertErrorBand(bandname)

    print("\nFHC:")
    print("  FHC CV integral      =", hfhc.Integral())
    print("  summed CV integral   =", summed_cv.Integral())
    print("  FHC bandCV integral  =", bfhc.Integral())
    print("  summed bandCV int    =", summed_bandcv.Integral())

    d_cv, b_cv = max_abs_diff(hfhc, summed_cv)
    d_bcv, b_bcv = max_abs_diff(bfhc, summed_bandcv)

    print("\nClosure:")
    print("  FHC parent CV vs sum(input CV):       maxdiff={} at bin {}".format(d_cv, b_cv))
    print("  FHC band CV vs sum(input band CV):    maxdiff={} at bin {}".format(d_bcv, b_bcv))

    # Universe closure
    bad = 0
    worst = (0.0, -1, -1)
    for i in range(bfhc.GetNHists()):
        d, b = max_abs_diff(bfhc.GetHist(i), summed_univs[i])
        if d > worst[0]:
            worst = (d, i, b)
        if d > 1e-8:
            bad += 1

    print("  Universe closure bad count:", bad, "/", bfhc.GetNHists())
    print("  Worst universe closure: maxdiff={} universe={} bin={}".format(*worst))

    # Print first few universe integral comparisons
    print("\nFirst few universe integral closure:")
    for i in range(min(10, bfhc.GetNHists())):
        fhc_u = bfhc.GetHist(i).Integral()
        sum_u = summed_univs[i].Integral()
        print(
            "  u{:3d}: FHC={:.8g} sumInputs={:.8g} diff={:.8g} FHC/CV={:.6g}".format(
                i,
                fhc_u,
                sum_u,
                fhc_u - sum_u,
                fhc_u / hfhc.Integral() if hfhc.Integral() else 0.0,
            )
        )

for hname in HISTS:
    check_one(hname, BAND)
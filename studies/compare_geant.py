#!/usr/bin/env python3

import ROOT
ROOT.gROOT.SetBatch(True)

HIST = "EN4"

geant_on = "/exp/minerva/data/users/qvuong/nu_mu/kin_dist_mcle1_Debug_CCnumu_allCuts_geantON_MAD.root"
geant_off = "/exp/minerva/data/users/qvuong/nu_mu/kin_dist_mcle1_Debug_CCnumu_allCuts_geantOFF_MAD.root"


def get_hist(path):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError("Could not open {}".format(path))

    h = f.Get(HIST)
    if not h:
        raise RuntimeError("Could not find {} in {}".format(HIST, path))

    h = h.GetCVHistoWithStatError().Clone()
    h.SetDirectory(0)
    f.Close()

    return h


h_on = get_hist(geant_on)
h_off = get_hist(geant_off)

print()
print("=" * 95)
print(
    "{:<14s} {:>16s} {:>16s} {:>14s} {:>14s}".format(
        "E_nu bin",
        "Geant ON",
        "Geant OFF",
        "ON/OFF",
        "change (%)",
    )
)
print("=" * 95)

for i in range(1, h_on.GetNbinsX() + 1):
    lo = h_on.GetXaxis().GetBinLowEdge(i)
    hi = h_on.GetXaxis().GetBinUpEdge(i)

    n_on = h_on.GetBinContent(i)
    n_off = h_off.GetBinContent(i)

    ratio = n_on / n_off if n_off else 0.0
    change = 100.0 * (ratio - 1.0)

    print(
        "{:5.2f}-{:<5.2f} "
        "{:16.6f} {:16.6f} "
        "{:14.6f} {:+13.3f}%".format(
            lo, hi,
            n_on,
            n_off,
            ratio,
            change,
        )
    )


print("\nEXACT DIFFERENCE CHECK")
print("=" * 95)

max_abs_diff = 0.0

for i in range(0, h_on.GetNbinsX() + 2):
    a = h_on.GetBinContent(i)
    b = h_off.GetBinContent(i)
    diff = a - b

    max_abs_diff = max(max_abs_diff, abs(diff))

    print(
        "bin {:2d}: ON={:20.12f}  OFF={:20.12f}  diff={:+20.12e}".format(
            i, a, b, diff
        )
    )

print("\nMAX ABS DIFF =", max_abs_diff)
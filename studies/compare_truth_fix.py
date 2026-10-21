#!/usr/bin/env python3

import ROOT
ROOT.gROOT.SetBatch(True)

HIST = "EN4"

base = "/exp/minerva/data/users/qvuong/nu_mu/"

paths = {
    "FV": {
        "data": base + "kin_dist_datale1_Debug_CCnumu_allCuts_FV_MAD.root",
        "mc":   base + "kin_dist_mcle1_Debug_CCnumu_allCuts_FV_MAD.root",
    },

    "InverseHasNoBackExitingTracks": {
        "data": base + "kin_dist_datale1_Debug_CCnumu_allCuts_FV_InverseHasNoBackExitingTracks_MAD.root",
        "mc":   base + "kin_dist_mcle1_Debug_CCnumu_allCuts_FV_InverseHasNoBackExitingTracks_MAD.root",
    },

    "Eavail": {
        "data": base + "kin_dist_datale1_Debug_CCnumu_FV_Eavail_MAD.root",
        "mc":   base + "kin_dist_mcle1_Debug_CCnumu_FV_Eavail_MAD.root",
    },

    "Pt": {
        "data": base + "kin_dist_datale1_Debug_CCnumu_FV_Pt_MAD.root",
        "mc":   base + "kin_dist_mcle1_Debug_CCnumu_FV_Pt_MAD.root",
    },

    "Etheta": {
        "data": base + "kin_dist_datale1_Debug_CCnumu_FV_Etheta_MAD.root",
        "mc":   base + "kin_dist_mcle1_Debug_CCnumu_FV_Etheta_MAD.root",
    },

    "LeptonAngle": {
        "data": base + "kin_dist_datale1_Debug_CCnumu_FV_LeptonAngle_MAD.root",
        "mc":   base + "kin_dist_mcle1_Debug_CCnumu_FV_LeptonAngle_MAD.root",
    },
}


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


h = {}

for cut, d in paths.items():
    h[cut] = {
        "data": get_hist(d["data"]),
        "mc": get_hist(d["mc"]),
    }


cuts = [
    "InverseHasNoBackExitingTracks",
    "Eavail",
    "Pt",
    "Etheta",
    "LeptonAngle",
]

h_data_fv = h["FV"]["data"]
h_mc_fv = h["FV"]["mc"]


print()
print("=" * 100)
print("DATA/MC EFFICIENCY RATIO RELATIVE TO FV-ONLY")
print("=" * 100)

header = "{:<13s}".format("E_nu bin")
for cut in cuts:
    header += " {:>15s}".format(cut[:15])

print(header)
print("-" * len(header))


for i in range(1, h_data_fv.GetNbinsX() + 1):

    lo = h_data_fv.GetXaxis().GetBinLowEdge(i)
    hi = h_data_fv.GetXaxis().GetBinUpEdge(i)

    d_fv = h_data_fv.GetBinContent(i)
    m_fv = h_mc_fv.GetBinContent(i)

    line = "{:5.2f}-{:<5.2f}".format(lo, hi)

    for cut in cuts:
        d_cut = h[cut]["data"].GetBinContent(i)
        m_cut = h[cut]["mc"].GetBinContent(i)

        eff_d = d_cut / d_fv if d_fv else 0.0
        eff_m = m_cut / m_fv if m_fv else 0.0

        ratio = eff_d / eff_m if eff_m else 0.0

        line += " {:15.5f}".format(ratio)

    print(line)


print()
print("=" * 100)
print("INDIVIDUAL DATA AND MC EFFICIENCIES")
print("=" * 100)

for cut in cuts:

    print("\n--- {} ---".format(cut))
    print(
        "{:<13s} {:>12s} {:>12s} {:>12s}".format(
            "E_nu bin",
            "Data eff",
            "MC eff",
            "Data/MC",
        )
    )

    for i in range(1, h_data_fv.GetNbinsX() + 1):

        lo = h_data_fv.GetXaxis().GetBinLowEdge(i)
        hi = h_data_fv.GetXaxis().GetBinUpEdge(i)

        d_fv = h_data_fv.GetBinContent(i)
        m_fv = h_mc_fv.GetBinContent(i)

        d_cut = h[cut]["data"].GetBinContent(i)
        m_cut = h[cut]["mc"].GetBinContent(i)

        eff_d = d_cut / d_fv if d_fv else 0.0
        eff_m = m_cut / m_fv if m_fv else 0.0
        ratio = eff_d / eff_m if eff_m else 0.0

        print(
            "{:5.2f}-{:<5.2f} {:12.5f} {:12.5f} {:12.5f}".format(
                lo,
                hi,
                eff_d,
                eff_m,
                ratio,
            )
        )
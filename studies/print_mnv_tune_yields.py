#!/usr/bin/env python3

import ROOT
ROOT.gROOT.SetBatch(True)

SAMPLES = {
    "FHC": ("/exp/minerva/data/users/qvuong/nu_mu", "leFHC", "CCnumu"),
    "RHC": ("/exp/minerva/data/users/qvuong/antinu_mu", "le5", "CCnumubar"),
}

TAGS = [
    "v1_D2_JOINT_LECOH_DiffCOH",
    "v1_D2_LECOH_DiffCOH",
    "v4_D2_MENU1PI_LECOH_DiffCOH",
    "v4_D2_MENU1PI_MECOH_DiffCOH",
    "v4p3p0_D2_MENU1PI_LECOH_DiffCOH",
    "v4p3p0_D2_MENU1PI_MECOH_DiffCOH",
]

REF = "v1_D2_LECOH_DiffCOH"

HISTS = {
    "Flux": "True Neutrino Energy",
    "EN4": "Reconstructed Energy Estimator",
}

def get_hist(f, name):
    h = f.Get(name)
    if not h:
        return None
    if hasattr(h, "GetCVHistoWithStatError"):
        h = h.GetCVHistoWithStatError()
    else:
        h = h.Clone()
    h.SetDirectory(0)
    return h

def get_pot(f):
    tree = f.Get("Meta")
    if not tree:
        return None
    pot = sum(float(e.POT_Used) for e in tree)
    return pot if pot > 0 else None

for mode, (directory, playlist, sample) in SAMPLES.items():

    files = {}
    pots = {}

    for tag in TAGS:
        path = (
            f"{directory}/kin_dist_mc{playlist}_"
            f"Validation_{sample}_{tag}_MAD.root"
        )
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            print("MISSING:", path)
            continue
        files[tag] = f
        pots[tag] = get_pot(f)

    if REF not in files:
        continue

    for name, label in HISTS.items():

        histograms = {}

        for tag, f in files.items():
            h = get_hist(f, name)
            if h is None:
                continue

            if pots[tag] and pots[REF]:
                h.Scale(pots[REF] / pots[tag])

            histograms[tag] = h

        if REF not in histograms:
            continue

        href = histograms[REF]

        print("\n" + "=" * 105)
        print(mode, "|", label)
        print("Reference:", REF)
        print("=" * 105)

        for tag, h in histograms.items():

            print("\nTUNE:", tag)
            print(
                f"{'Energy (GeV)':>17}"
                f"{'Reference':>16}"
                f"{'Tune':>16}"
                f"{'Ratio':>12}"
                f"{'Change (%)':>15}"
                f"{'Delta N':>16}"
            )

            for i in range(1, href.GetNbinsX() + 1):

                low = href.GetXaxis().GetBinLowEdge(i)
                high = href.GetXaxis().GetBinUpEdge(i)

                nref = href.GetBinContent(i)
                ntune = h.GetBinContent(i)

                ratio = ntune / nref if nref > 0 else float("nan")
                change = 100 * (ratio - 1)
                delta = ntune - nref

                print(
                    f"{low:7.2f}-{high:<7.2f}"
                    f"{nref:16.2f}"
                    f"{ntune:16.2f}"
                    f"{ratio:12.4f}"
                    f"{change:15.2f}"
                    f"{delta:16.2f}"
                )

    for f in files.values():
        f.Close()

import os
import ROOT

ROOT.gROOT.SetBatch(True)

from tools import Utilities

base = "/exp/minerva/data/users/qvuong/nu_mu"
prefix = "kin_dist_mcle13C_Debug_CCnumu_v4_oldLowQ2_MAD"

files = {
    "standard": f"{base}/{prefix}_no2p2h.root",
    "special": f"{base}/{prefix}_only2p2h.root",
    "combined": f"{base}/{prefix}.root",
}

pots = {}
hists = {}

for label, path in files.items():
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {path}")

    pots[label] = Utilities.getPOTFromFile(path)

    print(f"\n{label.upper()}: {os.path.basename(path)}")
    print(f"POT = {pots[label]:.8e}")

    for name in ["E_avail_CCNuMu2p2h", "Lepton_Pt_CCNuMu2p2h"]:
        h = f.Get(name)
        if not h:
            print(f"  Missing: {name}")
            continue

        value = sum(h.GetBinContent(i) for i in range(h.GetNcells()))
        hists[label, name] = value

        print(f"  {name}: integral={value:.8f}")

    f.Close()

scale = pots["standard"] / pots["special"]

print(f"\nPOT scale = {scale:.8f}")

for name in ["E_avail_CCNuMu2p2h", "Lepton_Pt_CCNuMu2p2h"]:
    standard = hists.get(("standard", name), 0.0)
    special = hists.get(("special", name), 0.0)
    combined = hists.get(("combined", name), 0.0)

    expected = standard + scale * special

    print(f"\n{name}")
    print(f"  Standard contribution: {standard:.8f}")
    print(f"  Scaled special:        {scale * special:.8f}")
    print(f"  Expected combined:     {expected:.8f}")
    print(f"  Actual combined:       {combined:.8f}")
    print(f"  Difference:            {combined - expected:.8e}")



f = ROOT.TFile.Open(files["combined"])

print("\n=== CCNuMu composition ===")

for key in f.GetListOfKeys():
    name = key.GetName()

    if not name.startswith("E_avail_"):
        continue

    h = key.ReadObj()

    if not h.InheritsFrom("TH1"):
        continue

    integral = sum(
        h.GetBinContent(i)
        for i in range(h.GetNcells())
    )

    print(f"{name:45s} {integral:15.3f}")

f.Close()
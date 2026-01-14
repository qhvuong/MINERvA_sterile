#!/usr/bin/env python3
import os
import glob
import ROOT

ROOT.gROOT.SetBatch(True)

base_dir = "/exp/minerva/app/users/qvuong/custom_plotutils/data/flux"
keep_hists = {"flux_E_cvweighted", "flux_E_unweighted"}

files = glob.glob(os.path.join(base_dir, "*.root"))

for path in files:
    fname = os.path.basename(path)
    print(f"\nProcessing {fname}")

    # Decide which histogram this file SHOULD contain
    if "gen2thin" in fname:
        target = "flux_E_cvweighted"
    elif "g4numiv5" in fname:
        target = "flux_E_unweighted"
    else:
        print("  -> Skipping (not gen2thin or g4numiv5)")
        continue

    f = ROOT.TFile.Open(path, "UPDATE")
    if not f or f.IsZombie():
        print("  -> ERROR: Could not open file")
        continue

    # List all objects in the file
    keys = [k.GetName() for k in f.GetListOfKeys()]

    # Delete everything except the one we want to keep
    for k in keys:
        if k != target:
            print(f"  Deleting: {k}")
            f.Delete(f"{k};*")   # delete all cycles of this object

    # Force rewrite of directory header
    f.Write("", ROOT.TObject.kOverwrite)
    f.Close()

    print(f"  -> Kept only: {target}")

#!/usr/bin/env python3

import glob
import ROOT

ROOT.gROOT.SetBatch(True)

base_dir = "/exp/minerva/data/users/qvuong/nu_e"

pattern = (
    f"{base_dir}/"
    "kin_dist_datale13*_CCnue_newResponse_MAD.root"
)

files = sorted(glob.glob(pattern))

if not files:
    raise RuntimeError(
        f"No files matched: {pattern}"
    )

print(f"Found {len(files)} files\n")

grand_totals = {}

for filename in files:
    print("=" * 100)
    print(filename)

    f = ROOT.TFile.Open(filename)

    if not f or f.IsZombie():
        print("  ERROR: could not open file")
        continue

    meta = f.Get("Meta")

    if not meta:
        print("  ERROR: no Meta tree")
        f.Close()
        continue

    branches = [
        b.GetName()
        for b in meta.GetListOfBranches()
    ]

    pot_branches = [
        b for b in branches
        if "pot" in b.lower()
    ]

    if not pot_branches:
        print("  No POT-like branches found.")
        print("  Available branches:")
        for b in branches:
            print(f"    {b}")
        f.Close()
        continue

    print("  POT-like branches:")
    for b in pot_branches:
        print(f"    {b}")

    print()

    file_sums = {
        b: 0.0
        for b in pot_branches
    }

    nentries = meta.GetEntries()

    print(f"  Meta entries = {nentries}")

    for i in range(nentries):
        meta.GetEntry(i)

        for b in pot_branches:
            try:
                value = getattr(meta, b)

                # Handle scalar-like ROOT objects if needed.
                try:
                    value = float(value)
                except Exception:
                    pass

                print(
                    f"    entry {i:3d}  "
                    f"{b:20s} = {value}"
                )

                if isinstance(value, (int, float)):
                    file_sums[b] += value

            except Exception as e:
                print(
                    f"    Could not read {b}: {e}"
                )

    print()
    print("  File sums:")

    for b, val in file_sums.items():
        print(
            f"    {b:20s} = {val:.12e}"
        )

        grand_totals[b] = (
            grand_totals.get(b, 0.0)
            + val
        )

    f.Close()

    print()

print("=" * 100)
print("SUM OVER ALL LE13 DATA FILES")
print("=" * 100)

for b, val in grand_totals.items():
    print(
        f"{b:20s} = {val:.12e}"
    )
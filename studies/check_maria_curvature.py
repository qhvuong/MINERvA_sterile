
#!/usr/bin/env python3

import argparse
import math
from collections import Counter

import ROOT

ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="Local ROOT file or root:// URL")
parser.add_argument("--max-events", type=int, default=50000)
parser.add_argument("--print-events", type=int, default=20)
args = parser.parse_args()

f = ROOT.TFile.Open(args.input_file)
if not f or f.IsZombie():
    raise RuntimeError("Could not open ROOT file")

# Locate the reco event tree.
trees = []
for key in f.GetListOfKeys():
    obj = key.ReadObj()
    if obj.InheritsFrom("TTree"):
        trees.append(obj)

required = [
    "MasterAnaDev_minos_used_curvature",
    "MasterAnaDev_minos_trk_eqp_qp",
    "MasterAnaDev_minos_trk_qp",
    "MasterAnaDev_nuHelicity",
    "isMinosMatchTrack",
]

tree = None
for candidate in trees:
    if all(candidate.GetBranch(name) for name in required):
        tree = candidate
        break

if tree is None:
    print("Available trees:")
    for candidate in trees:
        print(" ", candidate.GetName())
    print("\nMissing branches by tree:")
    for candidate in trees:
        missing = [name for name in required
                   if not candidate.GetBranch(name)]
        print(candidate.GetName(), missing)
    raise RuntimeError("No tree contains all required branches")

print("Tree:", tree.GetName())
print("Total entries:", tree.GetEntries())

def get_value(name):
    leaf = tree.GetLeaf(name)
    return float(leaf.GetValue(0))

counts = Counter()
curvature_values = []
examples = []

n = min(args.max_events, tree.GetEntries())

for i in range(n):
    tree.GetEntry(i)

    used = get_value("MasterAnaDev_minos_used_curvature")
    eqp = get_value("MasterAnaDev_minos_trk_eqp_qp")
    qp = get_value("MasterAnaDev_minos_trk_qp")
    helicity = get_value("MasterAnaDev_nuHelicity")
    matched = get_value("isMinosMatchTrack")

    counts["all"] += 1
    counts[f"helicity_{helicity:g}"] += 1
    counts[f"used_curvature_{used:g}"] += 1

    if matched != 1:
        continue

    counts["matched"] += 1

    # Work with the FHC helicity population.
    if helicity != 1:
        continue

    counts["matched_FHC"] += 1

    # Our charge selection, including the -1 sentinel.
    our_charge = (qp != -1 and qp < 0)

    if our_charge:
        counts["our_charge_pass"] += 1

    # Exactly reproduce Maria's curvature bypass logic.
    if used != 1:
        maria_pass = True
        significance = None
        counts["range_or_other"] += 1
    else:
        counts["curvature_reco"] += 1

        if eqp == 0 or not math.isfinite(eqp):
            significance = None
            maria_pass = False
            counts["invalid_eqp"] += 1
        else:
            significance = 1.0 / eqp
            curvature_values.append(significance)
            maria_pass = significance <= -5.0

            if significance < 0:
                counts["negative_significance"] += 1
            else:
                counts["positive_significance"] += 1

            if -5 < significance < 0:
                counts["negative_below_5"] += 1

            if significance <= -5:
                counts["negative_above_5"] += 1

    if maria_pass:
        counts["maria_pass"] += 1

    if our_charge and maria_pass:
        counts["both_pass"] += 1
    elif our_charge:
        counts["ours_only"] += 1
    elif maria_pass:
        counts["maria_only"] += 1
    else:
        counts["neither"] += 1

    if len(examples) < args.print_events:
        examples.append((
            i, used, qp, eqp, significance,
            our_charge, maria_pass
        ))

print("\n========== RAW COUNTS ==========")
for key, value in sorted(counts.items()):
    print(f"{key:30s}: {value}")

print("\n========== CUT ACCEPTANCE ==========")
den = counts["matched_FHC"]

if den:
    for key in [
        "curvature_reco",
        "range_or_other",
        "our_charge_pass",
        "maria_pass",
        "both_pass",
        "ours_only",
        "maria_only",
        "neither",
    ]:
        print(
            f"{key:30s}: "
            f"{counts[key]:8d} / {den} "
            f"({100.0 * counts[key] / den:6.2f}%)"
        )

print("\n========== CURVATURE SIGNIFICANCE ==========")
if curvature_values:
    values = sorted(curvature_values)

    for percentile in [0, 1, 5, 25, 50, 75, 95, 99, 100]:
        idx = round((len(values)-1) * percentile / 100)
        print(f"{percentile:3d}th percentile: {values[idx]:12.5g}")

print("\n========== EXAMPLE EVENTS ==========")
print("entry  usedCurv   qp        eqp_qp       1/eqp_qp    ours  maria")

for entry, used, qp, eqp, sig, ours, maria in examples:
    sig_str = f"{sig:12.5g}" if sig is not None else "       N/A  "
    print(
        f"{entry:5d}  {used:8g}  {qp:9.4g} "
        f"{eqp:12.5g} {sig_str}  "
        f"{str(ours):5s} {str(maria):5s}"
    )

f.Close()

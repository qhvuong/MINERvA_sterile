
#!/usr/bin/env python3

import argparse
from collections import Counter

import ROOT

ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser()
parser.add_argument("input_file")
parser.add_argument("--max-events", type=int, default=-1)
parser.add_argument("--show", type=int, default=15)
args = parser.parse_args()

f = ROOT.TFile.Open(args.input_file)
if not f or f.IsZombie():
    raise RuntimeError("Cannot open ROOT file")

tree = f.Get("MasterAnaDev")
if not tree:
    raise RuntimeError("MasterAnaDev tree not found")

required = [
    "has_interaction_vertex",
    "MasterAnaDev_vtx",
]

for name in required:
    if not tree.GetBranch(name):
        raise RuntimeError(f"Missing branch: {name}")

# Match our CCNuMu fiducial configuration
ZMIN = 5980.0
ZMAX = 8422.0
APOTHEM = 850.0

# Exact regular-hexagon apothem, with vertices
# oriented along the Y direction (MINERvA convention).
SQRT3 = 3.0 ** 0.5

def pass_hexagon(x, y):
    x = abs(x)
    y = abs(y)

    if x > APOTHEM:
        return False

    if y > 2.0 * APOTHEM / SQRT3:
        return False

    return y <= (
        2.0 * APOTHEM / SQRT3 - x / SQRT3
    )

counts = Counter()
examples = []

n = tree.GetEntries()
if args.max_events > 0:
    n = min(n, args.max_events)

print(f"Tree: {tree.GetName()}")
print(f"Scanning {n:,} events")

for i in range(n):
    tree.GetEntry(i)

    has_vertex = int(tree.has_interaction_vertex)

    vtx = getattr(tree, "MasterAnaDev_vtx")
    x, y, z = float(vtx[0]), float(vtx[1]), float(vtx[2])

    pass_z = ZMIN <= z <= ZMAX
    pass_xy = pass_hexagon(x, y)
    pass_fid = pass_z and pass_xy
    pass_vertex = has_vertex == 1

    counts["Total"] += 1
    counts[f"has_interaction_vertex={has_vertex}"] += 1

    if pass_vertex:
        counts["Has interaction vertex"] += 1

    if pass_z:
        counts["Pass Vertex_Z"] += 1

    if pass_xy:
        counts["Pass Vertex_Apothem"] += 1

    if pass_fid:
        counts["Pass both fiducial cuts"] += 1

        if pass_vertex:
            counts["Fiducial AND has vertex"] += 1
        else:
            counts["Fiducial BUT NO vertex"] += 1

            if len(examples) < args.show:
                examples.append((i, has_vertex, x, y, z))

    if pass_z and not pass_vertex:
        counts["Pass Z BUT NO vertex"] += 1

    if pass_vertex and not pass_fid:
        counts["Has vertex BUT fails fiducial"] += 1

print("\n========== COUNTS ==========")

keys = [
    "Total",
    "Has interaction vertex",
    "Pass Vertex_Z",
    "Pass Vertex_Apothem",
    "Pass both fiducial cuts",
    "Fiducial AND has vertex",
    "Fiducial BUT NO vertex",
    "Pass Z BUT NO vertex",
    "Has vertex BUT fails fiducial",
]

for key in keys:
    print(f"{key:38s}: {counts[key]:8d}")

print("\n========== VERTEX FLAG VALUES ==========")
for key in sorted(counts):
    if key.startswith("has_interaction_vertex="):
        print(f"{key:38s}: {counts[key]:8d}")

fid = counts["Pass both fiducial cuts"]
missing = counts["Fiducial BUT NO vertex"]

print("\n========== EFFECT ON OUR SELECTION ==========")

if fid:
    print(f"Events passing both fiducial cuts: {fid}")
    print(f"Additional events rejected by Maria: {missing}")
    print(f"Additional rejection: {100.0 * missing / fid:.4f}%")
else:
    print("No events pass the fiducial selection")

print("\n========== EXAMPLE REJECTED EVENTS ==========")
if examples:
    print("Entry  flag      X(mm)       Y(mm)       Z(mm)")
    for i, flag, x, y, z in examples:
        print(
            f"{i:5d}  {flag:4d}  "
            f"{x:10.2f}  {y:10.2f}  {z:10.2f}"
        )
else:
    print("No events pass fiducial cuts while failing the vertex flag.")

f.Close()

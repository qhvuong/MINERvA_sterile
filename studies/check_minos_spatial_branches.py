#!/usr/bin/env python3

import ROOT

f = ROOT.TFile.Open(
    "root://fndcadoor.fnal.gov:1094/"
    "pnfs/fnal.gov/usr/minerva/persistent/DataPreservation/"
    "Low_Energy_Era/p8/StandardMC/Playlist_minerva1/"
    "MasterAnaDev_mc_AnaTuple_run00010200_Playlist.root"
)

t = f.Get("MasterAnaDev")

expr = (
    "MasterAnaDev_t_minos_trk_primFSLepMinosInitProjZ:"
    "MasterAnaDev_t_minos_trk_primFSLepMinosInitProjY:"
    "MasterAnaDev_t_minos_trk_primFSLepMinosInitProjX"
)

cut = "MasterAnaDev_t_minos_trk_primFSLepMinosInitProjZ > -900000"

t.Scan(
    "MasterAnaDev_minos_trk_vtx_x:"
    "MasterAnaDev_minos_trk_vtx_y:"
    "MasterAnaDev_minos_trk_vtx_z",
    "MasterAnaDev_minos_trk_vtx_z > -900000",
    "",
    30
)


from collections import Counter

branch = "MasterAnaDev_minos_trk_vtx_z"

n = t.Draw(
    branch,
    branch + " > 0",
    "goff"
)

values = [t.GetV1()[i] for i in range(n)]

counts = Counter(round(z, 3) for z in values)

print("Selected:", n)
print("Minimum valid z:", min(values))

print("\nFirst 15 z planes:")
for z, count in sorted(counts.items())[:15]:
    print(f"{z:12.3f}   {count}")
#!/usr/bin/env python3

import ROOT

ROOT.gROOT.SetBatch(True)

file_old = (
    "/exp/minerva/data/users/qvuong/elastic_nue/"
    "kin_dist_mcmeFHC_p8_scattering_LElike_6bins_MAD.root"
)

file_new = (
    "/exp/minerva/data/users/qvuong/elastic_nue/"
    "kin_dist_mcme1N_p8_scattering_LElike_6bins_"
    "V1D2_LECOH_DiffCOH_noLowQ2_MAD.root"
)

hist_names = [
    "drawnL_ElepReco_LE_NuEElasticE",
    "drawnL_ElepReco_LE_NuEElasticEBar",
    "drawnL_ElepReco_LE_NuEElasticMu",
    "drawnL_ElepReco_LE_NuEElasticMuBar",
]


def get_edges(axis):
    edges = [
        axis.GetBinLowEdge(i)
        for i in range(1, axis.GetNbins() + 1)
    ]
    edges.append(axis.GetBinUpEdge(axis.GetNbins()))
    return edges


f_old = ROOT.TFile.Open(file_old)
f_new = ROOT.TFile.Open(file_new)

if not f_old or f_old.IsZombie():
    raise RuntimeError("Could not open old file")

if not f_new or f_new.IsZombie():
    raise RuntimeError("Could not open new file")


for hname in hist_names:

    h_old = f_old.Get(hname)
    h_new = f_new.Get(hname)

    if not h_old:
        print("\nMISSING in old file:", hname)
        continue

    if not h_new:
        print("\nMISSING in new file:", hname)
        continue

    print("\n" + "=" * 80)
    print(hname)
    print("=" * 80)

    for axis_name, old_axis, new_axis in [
        ("X", h_old.GetXaxis(), h_new.GetXaxis()),
        ("Y", h_old.GetYaxis(), h_new.GetYaxis()),
    ]:

        old_edges = get_edges(old_axis)
        new_edges = get_edges(new_axis)

        print("\n{} axis".format(axis_name))
        print("  old title :", old_axis.GetTitle())
        print("  new title :", new_axis.GetTitle())

        print("  old nbins :", old_axis.GetNbins())
        print("  new nbins :", new_axis.GetNbins())

        print("  old edges :", old_edges)
        print("  new edges :", new_edges)

        if old_edges == new_edges:
            print("  MATCH")
        else:
            print("  DIFFER")

            n = min(len(old_edges), len(new_edges))

            for i in range(n):
                if old_edges[i] != new_edges[i]:
                    print(
                        "    edge {:2d}: old = {:12.6g}, new = {:12.6g}".format(
                            i,
                            old_edges[i],
                            new_edges[i],
                        )
                    )

            if len(old_edges) != len(new_edges):
                print(
                    "    different number of edges: old={}, new={}".format(
                        len(old_edges),
                        len(new_edges),
                    )
                )


f_old.Close()
f_new.Close()
import ROOT

files = {
    "FHC": "/exp/minerva/app/users/qvuong/MAT_AL9/opt/lib/data/flux/flux-gen2thin-pdg14-minervame1D.root",
    "RHC": "/exp/minerva/app/users/qvuong/MAT_AL9/opt/lib/data/flux/flux-gen2thin-pdg-14-minervame5A.root",
}

hname = "flux_E_cvweighted"

edges = {}

for tag, path in files.items():
    f = ROOT.TFile.Open(path)
    h = f.Get(hname)

    e = [
        h.GetXaxis().GetBinLowEdge(i)
        for i in range(1, h.GetNbinsX() + 1)
    ]
    e.append(h.GetXaxis().GetBinUpEdge(h.GetNbinsX()))

    edges[tag] = e

    print("\n", tag)
    print("nbins =", h.GetNbinsX())
    print("edges =", e)

print("\nSame number of edges:", len(edges["FHC"]) == len(edges["RHC"]))

if len(edges["FHC"]) == len(edges["RHC"]):
    print(
        "Exactly same binning:",
        all(abs(a-b) < 1e-12 for a, b in zip(edges["FHC"], edges["RHC"]))
    )
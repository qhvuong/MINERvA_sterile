import ROOT

normal = ROOT.TFile.Open("/exp/minerva/data/users/qvuong/flux_studies/producedFluxes_p8_ME/ME5A_numubar.root")
offset = ROOT.TFile.Open("/exp/minerva/data/users/qvuong/flux_studies/producedFluxes_p8_ME_offsetPPFX/ME5A_numubar.root")

# Find the first MnvH1D
def find_hist(directory):
    for key in directory.GetListOfKeys():
        obj = key.ReadObj()
        if obj.InheritsFrom("PlotUtils::MnvH1D"):
            return obj
        if obj.InheritsFrom("TDirectory"):
            h = find_hist(obj)
            if h:
                return h
    return None

hn = find_hist(normal)
ho = find_hist(offset)

print("Histogram:", hn.GetName())

for band in ["ppfx1_Total", "Flux_BeamFocus", "Flux"]:
    print("\nBand:", band)

    bn = hn.GetVertErrorBand(band)
    bo = ho.GetVertErrorBand(band)

    for u in [0,1,2,36,37,38]:
        h1 = bn.GetHist(u)
        h2 = bo.GetHist(u)

        diff = 0.
        for b in range(h1.GetNbinsX()+2):
            diff = max(diff, abs(h1.GetBinContent(b)-h2.GetBinContent(b)))

        print(
            f"u={u:3d}"
            f"  normal={h1.Integral(0,h1.GetNbinsX()+1):.8e}"
            f"  offset={h2.Integral(0,h2.GetNbinsX()+1):.8e}"
            f"  maxdiff={diff:.3e}"
        )
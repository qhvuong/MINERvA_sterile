import ROOT
import math

f = ROOT.TFile.Open(
    "/exp/minerva/data/users/qvuong/antinu_mu/kin_dist_mcle5_p8_debug_nonfinite_job6_MAD_6.root"
)

for key in f.GetListOfKeys():
    h = key.ReadObj()

    if not h.InheritsFrom("TH1"):
        continue

    bad_bins = [
        b for b in range(h.GetNcells())
        if not math.isfinite(h.GetBinContent(b))
        or not math.isfinite(h.GetBinError(b))
    ]

    if bad_bins:
        print("[NONFINITE HIST]", h.GetName(), bad_bins)

f.Close()
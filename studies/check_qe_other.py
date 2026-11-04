
import ROOT
import numpy as np

ROOT.gROOT.SetBatch(True)

maria_path = (
    "/exp/minerva/data/users/mmehmood/"
    "2026-10-24-NoNuEConstraint/"
    "MnvTunev1_minerva1_NoNuE_v5_2026-10-27/"
    "Hists_EventSelection_MnvTunev1_minerva1_NoNuE_v5_2026-10-27"
    "_minerva1_nosys_t99_z99_Nu.root"
)

ours_path = (
    "/exp/minerva/data/users/qvuong/nu_mu/"
    "kin_dist_mcle1_ValidateCuts_CCnumu_Inclusive_Curvature_v1_MAD.root"
)

fm = ROOT.TFile.Open(maria_path)
fo = ROOT.TFile.Open(ours_path)

def values(file, name, transpose=False):
    h = file.Get(name)
    if not h:
        raise RuntimeError(f"Missing histogram: {name}")
    a = np.array([
        [h.GetBinContent(ix, iy)
         for iy in range(1, h.GetNbinsY() + 1)]
        for ix in range(1, h.GetNbinsX() + 1)
    ])
    return a.T if transpose else a

m_qe = values(fm, "reco_QE_pZmu_pTmu")
m_other = values(fm, "reco_other_pZmu_pTmu")

# Ours: x=pT, y=ppar, so transpose to Maria's orientation
o_qe = values(fo, "Lepton_Pt_Ppar_CCNuMuQE", True)
o_other = values(fo, "Lepton_Pt_Ppar_CCNuMu", True)

for label, m, o in [
    ("QE", m_qe, o_qe),
    ("Other", m_other, o_other),
    ("QE + Other", m_qe + m_other, o_qe + o_other),
]:
    diff = o - m
    mask = m >= 100
    rel = np.abs(diff[mask] / m[mask]) if np.any(mask) else []

    print(f"\n{label}")
    print(f"Maria integral: {m.sum():.4f}")
    print(f"Ours integral:  {o.sum():.4f}")
    print(f"Difference:     {diff.sum():+.4f}")
    print(f"Sum abs diff:   {np.abs(diff).sum():.4f}")
    if len(rel):
        print(f"Max relative bin diff (Maria >=100): {max(rel)*100:.3f}%")

fm.Close()
fo.Close()

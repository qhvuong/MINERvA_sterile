import ROOT

base = "/exp/minerva/data/users/qvuong"

files = {
    "LE1 DATA": f"{base}/nu_mu/kin_dist_datale1_ValidateCuts_CCnumu_Inclusive_Curvature_v1_MAD.root",
    "LE1 MC": f"{base}/nu_mu/kin_dist_mcle1_ValidateCuts_CCnumu_Inclusive_Curvature_v1_MAD.root",

    "LE5 DATA": f"{base}/antinu_mu/kin_dist_datale5_ValidateCuts_CCnumubar_Inclusive_Curvature_v1_MAD.root",
    "LE5 STANDARD MC": f"{base}/antinu_mu/kin_dist_mcle5_ValidateCuts_CCnumubar_Inclusive_Curvature_v1_MAD_no2p2h.root",
    "LE5 COMBINED MC": f"{base}/antinu_mu/kin_dist_mcle5_ValidateCuts_CCnumubar_Inclusive_Curvature_v1_MAD.root",
}

for label, path in files.items():
    f = ROOT.TFile.Open(path)

    if not f or f.IsZombie():
        print(f"{label}: FILE MISSING OR INVALID")
        continue

    meta = f.Get("Meta")

    if not meta or not meta.GetBranch("POT_Used"):
        print(f"{label}: POT_Used unavailable")
    else:
        pot = sum(float(row.POT_Used) for row in meta)
        print(f"{label:20s} POT = {pot:.12e}  entries = {meta.GetEntries()}")

    f.Close()
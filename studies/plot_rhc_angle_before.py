import ROOT

ROOT.gROOT.SetBatch(True)

base = "/exp/minerva/data/users/qvuong/antinu_mu"
tag = "ValidateCuts_CCnumubar_FV_InverseHasNoBackExitingTracks_DeadTime_LeptonAngle_v1"

fd = ROOT.TFile.Open(f"{base}/kin_dist_datale5_{tag}_MAD.root")
fm = ROOT.TFile.Open(f"{base}/kin_dist_mcle5_{tag}_MAD.root")

data = fd.Get("Lep_Angle").Clone("data_angle")
mc = fm.Get("Lep_Angle").Clone("mc_angle")

data.SetDirectory(0)
mc.SetDirectory(0)

fd.Close()
fm.Close()

# Normalize each histogram to unit area.
if data.Integral() <= 0 or mc.Integral() <= 0:
    raise RuntimeError("Invalid histogram integral")

data.Scale(1.0 / data.Integral())
mc.Scale(1.0 / mc.Integral())

data.SetMarkerStyle(20)
data.SetMarkerColor(ROOT.kBlack)
data.SetLineColor(ROOT.kBlack)
mc.SetLineColor(ROOT.kRed + 1)
mc.SetLineWidth(2)

c = ROOT.TCanvas("c", "Muon angle comparison", 850, 650)

mc.Draw("HIST")
data.Draw("E1 SAME")

mc.GetXaxis().SetTitle("Reconstructed muon angle")
mc.GetYaxis().SetTitle("Normalized events")

leg = ROOT.TLegend(0.63, 0.75, 0.88, 0.88)
leg.AddEntry(data, "RHC Data", "lep")
leg.AddEntry(mc, "RHC MC", "l")
leg.Draw()

c.SaveAs("rhc_lepton_angle_shape.png")
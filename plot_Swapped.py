import ROOT

# File and POT target
file_path = "/exp/minerva/data/users/qvuong/nu_e/kin_dist_mcle1A_p6_SwappedSample_MAD.root"
target_POT = 1.1e21

# Open the ROOT file
f = ROOT.TFile.Open(file_path)

# Get the Meta tree
meta = f.Get("Meta")

# Safely read POT_Used from first entry
meta.GetEntry(0)
POT_Used = getattr(meta, "POT_Used")

print(f"Original POT: {POT_Used:.3e}")

# Scale factor
scale = target_POT / POT_Used
print(f"Scaling histos by: {scale:.3f}")

# Get histograms
hm = f.Get("hm").Clone("hm_scaled")
he = f.Get("he").Clone("he_scaled")
hm.SetStats(0)
he.SetStats(0)
hm.SetMaximum(hm.GetMaximum()*100.0)

hm.Scale(scale)
he.Scale(scale)

# Draw
canvas = ROOT.TCanvas("canvas", "Histos", 800, 600)
canvas.SetLogy()
he.SetLineColor(ROOT.kRed)
hm.SetLineColor(ROOT.kBlue)

he.SetTitle("; E_{#nu} (GeV); Events/yr.POT")
hm.Draw("HIST")
he.Draw("HIST SAME")

legend = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
legend.AddEntry(he, "Truth", "l")
legend.AddEntry(hm, "Swapped Sample", "l")
legend.Draw()

canvas.SaveAs("nu.png")


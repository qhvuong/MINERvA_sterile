import ROOT
import os

# Input ROOT file
filename = "le5_ebar.root"
f = ROOT.TFile.Open(filename)

# Parent PDGs
parent_pdgs = [321, -321, 130, 311, 13, -13, 211, -211]
other_pdg = 999999

PDG_NAMES = {
    211: "pi+",
    -211: "pi-",
    321: "K+",
    -321: "K-",
    130: "KL",
    13: "mu-",
    -13: "mu+",
    311: "K0",
    999999: "Others"
}

colors = {
    211: ROOT.kRed + 1,
    -211: ROOT.kRed - 3,
    321: ROOT.kGreen + 2,
    -321: ROOT.kGreen - 3,  # New color for K-
    130: ROOT.kBlue + 1,
    13:  ROOT.kOrange + 7,
    -13: ROOT.kMagenta + 2,
    311: ROOT.kCyan + 1,
    other_pdg: ROOT.kGray + 2,
}

# Histogram stack
stack = ROOT.THStack("flux_stack", "Stacked Flux by Parent;E_{#nu} (GeV);Flux / m^{2} / POT / GeV")

# Legend
legend = ROOT.TLegend(0.65, 0.6, 0.88, 0.88)
legend.SetBorderSize(0)

# Add histograms
for pdg in parent_pdgs + [other_pdg]:
    hist_name = f"flux_E_cvweighted_parent{pdg}"
    hist = f.Get(hist_name)
    if not hist:
        print(f"Warning: Histogram {hist_name} not found.")
        continue
    if hist.Integral() == 0:
        print(f"Skipping empty histogram: {hist_name}")
        continue
    hist.SetLineColor(colors.get(pdg, ROOT.kBlack))
    hist.SetFillColor(colors.get(pdg, ROOT.kBlack))
    hist.SetMarkerStyle(0)
    stack.Add(hist)

    label = PDG_NAMES.get(pdg, str(pdg))
    legend.AddEntry(hist, label, "f")

# Draw
canvas = ROOT.TCanvas("c", "", 800, 600)
canvas.SetLogy()
stack.SetMinimum(1e-10)  # for example
stack.SetMaximum(1e-6)   # optional
stack.Draw("HIST")
stack.SetTitle("Flux Contribution by Parent (Stacked)")
stack.GetYaxis().SetTitleOffset(1.2)
legend.Draw()

# # Optional: draw total flux
# total = stack.GetStack().Last().Clone("total_flux")
# total.SetLineColor(ROOT.kBlack)
# total.SetLineWidth(2)
# total.SetFillStyle(0)
# total.Draw("SAME HIST")

# Save
canvas.SaveAs("ebarFlux_stackedParent.png")
# canvas.SaveAs("stacked_flux_by_parent.pdf")

import ROOT
import sys

# --- List of input ROOT files ---
nu_files = [
    "$PLOTUTILSROOT/data/flux/flux-g4numiv6-pdg14-minervame1D1M1NWeightedAve.root",     # numerator
    "$PLOTUTILSROOT/data/flux/flux-g4numiv6-pdg12-minervame1D1M1NWeightedAve.root",    # denominator part 1
]

nubar_files = [
    "$PLOTUTILSROOT/data/flux/flux-g4numiv6-pdg-14-minervame1D1M1NWeightedAve.root",     # numerator
    "$PLOTUTILSROOT/data/flux/flux-g4numiv6-pdg-12-minervame1D1M1NWeightedAve.root",    # denominator part 1
]

hist_name = "flux_E_unweighted"
#output_file = "hist_ratio.root"
output_hist_name = "hist_ratio"

names = ["nu", "nubar"]
filenames = [nu_files, nubar_files]
edges = []
ratios = [] 

nu_ratios = []
nubar_ratios = []

nu_edges = []
nubar_edges = []

# --- Load histograms ---
def LoadHistograms(out_name, input_files):
    hists = []
    for fname in input_files:
        f = ROOT.TFile.Open(fname, "READ")
        if not f or f.IsZombie():
            raise RuntimeError(f"Could not open file: {fname}")
        
        h = f.Get(hist_name)
        if not h:
            raise RuntimeError(f"Histogram '{hist_name}' not found in file {fname}")
        
        h.SetDirectory(0)  # Detach from file so it doesn't get deleted
        f.Close()
        hists.append(h)

    # --- Compute ratio ---
    hist_ratio = hists[0].Clone(output_hist_name)
    hist_ratio.Divide(hists[0], hists[1])
    
    # --- Get bin edges and values ---
    xaxis = hist_ratio.GetXaxis()
    nbins = xaxis.GetNbins()
    
    bin_edges = [xaxis.GetBinLowEdge(1 + i) for i in range(nbins)]
    bin_edges.append(xaxis.GetBinUpEdge(nbins))  # add final upper edge

    ratio = [hist_ratio.GetBinContent(i + 1) for i in range(nbins)]

    edges.append(bin_edges)
    ratios.append(ratio)

    # --- Write to txt ---
    with open(f"{out_name}_ratio.txt", "w") as out:
        out.write(str(bin_edges) + "\n")
        out.write(str(ratio) + "\n")

    # --- Save to new file ---
    out = ROOT.TFile.Open(f"{out_name}_ratio.root", "RECREATE")
    hist_ratio.Write()
    out.Close()



#print(f"Saved ratio histogram to {output_file}")

for i in range(2):
    LoadHistograms(names[i], filenames[i])

print(edges[0] == edges[1])

import ROOT, math, os
import PlotUtils

fname = "rootfiles/NuE_stitched_hists_prodNueel_noRatio.root"
# fname = "rootfiles/NuE_stitched_hists_prodNueel.root"

f = ROOT.TFile.Open(fname)
if not f or f.IsZombie():
    raise RuntimeError("Could not open " + fname)

data = f.Get("sample_data")
mc   = f.Get("sample_mc")

if not data or not mc:
    raise RuntimeError("Could not get sample_data/sample_mc")

h_test = data.Clone("h_test_data_minus_mc")
h_test.Add(mc, -1)

print("file:", fname)
print("nbins:", data.GetNbinsX())

print("\n{:>4s} {:>12s} {:>12s} {:>12s} | {:>12s} {:>12s} {:>12s} | {:>12s} {:>12s} {:>12s}".format(
    "bin",
    "data", "data_err", "sqrt(data)",
    "mc", "mc_err", "sqrt(mc)",
    "data-mc", "h_err", "quad_err"
))

for i in range(1, data.GetNbinsX()+1):
    d  = data.GetBinContent(i)
    de = data.GetBinError(i)

    m  = mc.GetBinContent(i)
    me = mc.GetBinError(i)

    r  = h_test.GetBinContent(i)
    re = h_test.GetBinError(i)

    sd = math.sqrt(abs(d)) if d >= 0 else float("nan")
    sm = math.sqrt(abs(m)) if m >= 0 else float("nan")
    qe = math.sqrt(de*de + me*me)

    print("{:4d} {:12.5g} {:12.5g} {:12.5g} | {:12.5g} {:12.5g} {:12.5g} | {:12.5g} {:12.5g} {:12.5g}".format(
        i, d, de, sd, m, me, sm, r, re, qe
    ))

# f.Close()
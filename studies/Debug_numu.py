#!/usr/bin/env python3

import ROOT
import math

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# ============================================================
# INPUT
# ============================================================

data_fv_path = (
    "/exp/minerva/data/users/qvuong/nu_mu/"
    "kin_dist_dataleFHC_Debug_CCnumu_onlyFV_MAD.root"
)

mc_fv_path = (
    "/exp/minerva/data/users/qvuong/nu_mu/"
    "kin_dist_mcleFHC_Debug_CCnumu_onlyFV_MAD.root"
)

data_sel_path = (
    "/exp/minerva/data/users/qvuong/nu_mu/"
    "kin_dist_dataleFHC_Debug_CCnumu_FV_InverseHasNoBackExitingTracks_MAD.root"
)

mc_sel_path = (
    "/exp/minerva/data/users/qvuong/nu_mu/"
    "kin_dist_mcleFHC_Debug_CCnumu_FV_InverseHasNoBackExitingTracks_MAD.root"
)

HIST_NAME = "EN4"


# ============================================================
# HELPERS
# ============================================================

def open_file(path):
    f = ROOT.TFile.Open(path)

    if not f or f.IsZombie():
        raise RuntimeError("Could not open:\n{}".format(path))

    return f


def get_cv_hist(root_file, hist_name, new_name):
    """
    Get MnvH1D and extract ordinary TH1D CV histogram.
    """
    h_mnv = root_file.Get(hist_name)

    if not h_mnv:
        raise RuntimeError(
            "Could not find histogram '{}' in {}".format(
                hist_name,
                root_file.GetName()
            )
        )

    # Extract CV histogram as an ordinary ROOT TH1D
    h = h_mnv.GetCVHistoWithStatError().Clone(new_name)
    h.SetDirectory(0)

    return h


# ============================================================
# OPEN FILES
# ============================================================

f_data_fv = open_file(data_fv_path)
f_mc_fv   = open_file(mc_fv_path)

f_data_sel = open_file(data_sel_path)
f_mc_sel   = open_file(mc_sel_path)


# ============================================================
# GET CV HISTOGRAMS
# ============================================================

h_data_fv = get_cv_hist(
    f_data_fv,
    HIST_NAME,
    "h_data_fv"
)

h_mc_fv = get_cv_hist(
    f_mc_fv,
    HIST_NAME,
    "h_mc_fv"
)

h_data_sel = get_cv_hist(
    f_data_sel,
    HIST_NAME,
    "h_data_sel"
)

h_mc_sel = get_cv_hist(
    f_mc_sel,
    HIST_NAME,
    "h_mc_sel"
)


# ============================================================
# SANITY CHECK: ARE NUMERATOR AND DENOMINATOR ACTUALLY DIFFERENT?
# ============================================================

print("\n============================================================")
print("RAW HISTOGRAM CHECK")
print("============================================================")

print("\nIntegrated EN4 yields:")
print("Data FV        =", h_data_fv.Integral())
print("Data selected  =", h_data_sel.Integral())
print("MC FV          =", h_mc_fv.Integral())
print("MC selected    =", h_mc_sel.Integral())

print("\nBin contents BEFORE efficiency calculation:")
print(
    "{:>4s} {:>10s} {:>14s} {:>14s} {:>14s} {:>14s}".format(
        "bin",
        "Enu",
        "Data FV",
        "Data sel",
        "MC FV",
        "MC sel",
    )
)

print("-" * 80)

identical_data = True
identical_mc = True

for i in range(1, h_data_fv.GetNbinsX() + 1):

    enu = h_data_fv.GetBinCenter(i)

    data_fv = h_data_fv.GetBinContent(i)
    data_sel = h_data_sel.GetBinContent(i)

    mc_fv = h_mc_fv.GetBinContent(i)
    mc_sel = h_mc_sel.GetBinContent(i)

    if abs(data_fv - data_sel) > 1e-9:
        identical_data = False

    if abs(mc_fv - mc_sel) > 1e-9:
        identical_mc = False

    print(
        "{:4d} {:10.3f} {:14.3f} {:14.3f} {:14.3f} {:14.3f}".format(
            i,
            enu,
            data_fv,
            data_sel,
            mc_fv,
            mc_sel,
        )
    )


print("\nHistogram comparison:")

if identical_data:
    print("WARNING: Data FV and Data selected EN4 histograms are IDENTICAL.")
else:
    print("Data FV and selected histograms are different.")

if identical_mc:
    print("WARNING: MC FV and MC selected EN4 histograms are IDENTICAL.")
else:
    print("MC FV and selected histograms are different.")


# ============================================================
# CALCULATE EFFICIENCIES MANUALLY
# ============================================================

h_eff_data = h_data_fv.Clone("h_eff_data")
h_eff_mc = h_mc_fv.Clone("h_eff_mc")
h_eff_ratio = h_data_fv.Clone("h_eff_ratio")

h_eff_data.Reset()
h_eff_mc.Reset()
h_eff_ratio.Reset()

print("\n============================================================")
print("BIN-BY-BIN SELECTION EFFICIENCY")
print("============================================================")

print(
    "{:>4s} {:>10s} {:>12s} {:>12s} {:>12s}".format(
        "bin",
        "Enu",
        "Data eff",
        "MC eff",
        "Data/MC",
    )
)

print("-" * 60)

for i in range(1, h_data_fv.GetNbinsX() + 1):

    enu = h_data_fv.GetBinCenter(i)

    n_data_fv = h_data_fv.GetBinContent(i)
    n_data_sel = h_data_sel.GetBinContent(i)

    n_mc_fv = h_mc_fv.GetBinContent(i)
    n_mc_sel = h_mc_sel.GetBinContent(i)

    # ----------------------------
    # Data efficiency
    # ----------------------------

    if n_data_fv > 0.0:
        eff_data = n_data_sel / n_data_fv

        # approximate binomial uncertainty
        variance = eff_data * (1.0 - eff_data) / n_data_fv
        err_data = math.sqrt(max(0.0, variance))
    else:
        eff_data = 0.0
        err_data = 0.0

    # ----------------------------
    # MC efficiency
    # ----------------------------

    if n_mc_fv > 0.0:
        eff_mc = n_mc_sel / n_mc_fv

        variance = eff_mc * (1.0 - eff_mc) / n_mc_fv
        err_mc = math.sqrt(max(0.0, variance))
    else:
        eff_mc = 0.0
        err_mc = 0.0

    # ----------------------------
    # Data / MC efficiency
    # ----------------------------

    if eff_mc > 0.0:
        ratio = eff_data / eff_mc
    else:
        ratio = 0.0

    h_eff_data.SetBinContent(i, eff_data)
    h_eff_data.SetBinError(i, err_data)

    h_eff_mc.SetBinContent(i, eff_mc)
    h_eff_mc.SetBinError(i, err_mc)

    h_eff_ratio.SetBinContent(i, ratio)

    print(
        "{:4d} {:10.3f} {:12.4f} {:12.4f} {:12.4f}".format(
            i,
            enu,
            eff_data,
            eff_mc,
            ratio,
        )
    )


# ============================================================
# INTEGRATED EFFICIENCY
# ============================================================

data_fv_total = h_data_fv.Integral()
data_sel_total = h_data_sel.Integral()

mc_fv_total = h_mc_fv.Integral()
mc_sel_total = h_mc_sel.Integral()

data_eff_total = (
    data_sel_total / data_fv_total
    if data_fv_total > 0 else 0.0
)

mc_eff_total = (
    mc_sel_total / mc_fv_total
    if mc_fv_total > 0 else 0.0
)

ratio_total = (
    data_eff_total / mc_eff_total
    if mc_eff_total > 0 else 0.0
)

print("\n============================================================")
print("INTEGRATED SELECTION EFFICIENCY")
print("============================================================")

print("Data efficiency =", data_eff_total)
print("MC efficiency   =", mc_eff_total)
print("Data / MC       =", ratio_total)


# ============================================================
# EFFICIENCY PLOT
# ============================================================

c = ROOT.TCanvas("c", "c", 800, 700)

h_eff_data.SetLineColor(ROOT.kBlack)
h_eff_data.SetMarkerColor(ROOT.kBlack)
h_eff_data.SetMarkerStyle(20)

h_eff_mc.SetLineColor(ROOT.kRed + 1)
h_eff_mc.SetMarkerColor(ROOT.kRed + 1)
h_eff_mc.SetMarkerStyle(24)

h_eff_data.SetTitle(
    "FHC CC#nu_{#mu} back-exiting-track selection efficiency;"
    "Reconstructed E_{#nu} [GeV];"
    "Selection efficiency"
)

h_eff_data.SetMinimum(0.0)
h_eff_data.SetMaximum(1.05)

h_eff_data.Draw("E1")
h_eff_mc.Draw("E1 SAME")

leg = ROOT.TLegend(0.62, 0.72, 0.88, 0.88)
leg.SetBorderSize(0)

leg.AddEntry(
    h_eff_data,
    "Data",
    "lep"
)

leg.AddEntry(
    h_eff_mc,
    "MC",
    "lep"
)

leg.Draw()

c.SaveAs(
    "CCnumu_BackExitingTrack_efficiency.png"
)


# ============================================================
# DATA / MC EFFICIENCY RATIO
# ============================================================

c2 = ROOT.TCanvas("c2", "c2", 800, 700)

h_eff_ratio.SetMarkerStyle(20)
h_eff_ratio.SetMarkerColor(ROOT.kBlack)
h_eff_ratio.SetLineColor(ROOT.kBlack)

h_eff_ratio.SetTitle(
    "FHC CC#nu_{#mu} back-exiting-track efficiency ratio;"
    "Reconstructed E_{#nu} [GeV];"
    "#epsilon_{Data}/#epsilon_{MC}"
)

h_eff_ratio.SetMinimum(0.5)
h_eff_ratio.SetMaximum(1.5)

h_eff_ratio.Draw("E1")

line = ROOT.TLine(
    h_eff_ratio.GetXaxis().GetXmin(),
    1.0,
    h_eff_ratio.GetXaxis().GetXmax(),
    1.0,
)

line.SetLineStyle(2)
line.Draw()

c2.SaveAs(
    "CCnumu_BackExitingTrack_efficiency_DataMC.png"
)
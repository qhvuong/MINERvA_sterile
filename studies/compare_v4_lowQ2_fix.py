#!/usr/bin/env python3

import os
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# ============================================================
# CONFIG
# ============================================================

HIST_NAME = "EN4"

BASE_DIR = "/exp/minerva/data/users/qvuong/nu_mu"
OUT_DIR = "/exp/minerva/data/users/qvuong/studies"

OLD_TAG = "DebugCuts_CCnumu_oldFV_MnvTunesV4"
FIXED_TAG = "Debug_CCnumu_v4_oldLowQ2"

CONFIGS = [
    {
        "tag": OLD_TAG,
        "label": "V4: old LowQ2 implementation",
        "color": ROOT.kRed + 1,
        "style": 2,
        "width": 4,
    },
    {
        "tag": FIXED_TAG,
        "label": "V4: corrected LowQ2 implementation",
        "color": ROOT.kBlue + 1,
        "style": 1,
        "width": 4,
    },
]


# ============================================================
# PATHS
# ============================================================

def data_path(tag):
    return os.path.join(
        BASE_DIR,
        "kin_dist_dataleFHC_{}_MAD.root".format(tag)
    )


def mc_path(tag):
    return os.path.join(
        BASE_DIR,
        "kin_dist_mcleFHC_{}_MAD.root".format(tag)
    )


# ============================================================
# HELPERS
# ============================================================

def get_pot(root_file):
    meta = root_file.Get("Meta")

    if not meta:
        raise RuntimeError(
            "Could not find Meta tree in {}".format(root_file.GetName())
        )

    pot = 0.0

    for entry in meta:
        pot += entry.POT_Used

    return pot


def get_hist(root_file, name):
    h_mnv = root_file.Get(name)

    if not h_mnv:
        raise RuntimeError(
            "Could not find histogram {} in {}".format(
                name,
                root_file.GetName()
            )
        )

    # CV only
    h = h_mnv.GetCVHistoWithStatError()

    h.SetName(
        "{}_cv_{}".format(
            name,
            os.path.basename(root_file.GetName()).replace(".root", "")
        )
    )

    h.SetDirectory(0)

    return h


def integral(h):
    # Include underflow and overflow
    return h.Integral(0, h.GetNbinsX() + 1)


def make_ratio(num, den, name):
    ratio = num.Clone(name)
    ratio.SetDirectory(0)
    ratio.Divide(num, den, 1.0, 1.0, "")
    return ratio


# ============================================================
# LOAD FILES
# ============================================================

results = {}

for cfg in CONFIGS:

    tag = cfg["tag"]

    print("")
    print("Loading:", tag)
    print("  Data:", data_path(tag))
    print("  MC  :", mc_path(tag))

    fd = ROOT.TFile.Open(data_path(tag))
    fm = ROOT.TFile.Open(mc_path(tag))

    if not fd or fd.IsZombie():
        raise RuntimeError(
            "Could not open {}".format(data_path(tag))
        )

    if not fm or fm.IsZombie():
        raise RuntimeError(
            "Could not open {}".format(mc_path(tag))
        )

    data_pot = get_pot(fd)
    mc_pot = get_pot(fm)

    h_data = get_hist(fd, HIST_NAME)
    h_mc = get_hist(fm, HIST_NAME)

    pot_scale = data_pot / mc_pot
    h_mc.Scale(pot_scale)

    results[tag] = {
        "cfg": cfg,
        "data_pot": data_pot,
        "mc_pot": mc_pot,
        "pot_scale": pot_scale,
        "data": h_data,
        "mc": h_mc,
        "data_yield": integral(h_data),
        "mc_yield": integral(h_mc),
    }

    fd.Close()
    fm.Close()


# ============================================================
# NUMERICAL SUMMARY
# ============================================================

old = results[OLD_TAG]
fixed = results[FIXED_TAG]

old_mc = old["mc_yield"]
fixed_mc = fixed["mc_yield"]

old_data = old["data_yield"]
fixed_data = fixed["data_yield"]

print("")
print("======================================================================")
print("CCnumu V4 LowQ2 implementation comparison: {}".format(HIST_NAME))
print("======================================================================")
print("")

print(
    "{:<38s} {:>14s} {:>14s} {:>14s}".format(
        "Configuration",
        "Data",
        "MC",
        "Data/MC",
    )
)

print("-" * 85)

for cfg in CONFIGS:

    tag = cfg["tag"]
    r = results[tag]

    data_mc = (
        r["data_yield"] / r["mc_yield"]
        if r["mc_yield"] else 0.0
    )

    print(
        "{:<38s} {:>14.3f} {:>14.3f} {:>14.6f}".format(
            cfg["label"],
            r["data_yield"],
            r["mc_yield"],
            data_mc,
        )
    )


# ============================================================
# DIRECT EFFECT OF LOWQ2 FIX
# ============================================================

print("")
print("===== EFFECT OF CORRECTING LOWQ2 IMPLEMENTATION =====")

mc_shift = (
    (fixed_mc / old_mc - 1.0) * 100.0
    if old_mc else 0.0
)

print(
    "Corrected V4 / old V4 MC = {:.6f}".format(
        fixed_mc / old_mc if old_mc else 0.0
    )
)

print(
    "MC normalization shift   = {:+.3f}%".format(
        mc_shift
    )
)

print(
    "Old V4 Data/MC           = {:.6f}".format(
        old_data / old_mc if old_mc else 0.0
    )
)

print(
    "Corrected V4 Data/MC     = {:.6f}".format(
        fixed_data / fixed_mc if fixed_mc else 0.0
    )
)


# ============================================================
# POT CHECK
# ============================================================

print("")
print("===== POT SUMMARY =====")

for cfg in CONFIGS:

    r = results[cfg["tag"]]

    print(
        "{:<38s} data POT = {:.8e}   MC POT = {:.8e}   scale = {:.8g}".format(
            cfg["label"],
            r["data_pot"],
            r["mc_pot"],
            r["pot_scale"],
        )
    )


# ============================================================
# DATA CONSISTENCY CHECK
# ============================================================

print("")
print("===== DATA CONSISTENCY CHECK =====")

old_data_hist = old["data"]
fixed_data_hist = fixed["data"]

max_diff = 0.0

for i in range(0, old_data_hist.GetNbinsX() + 2):

    diff = abs(
        old_data_hist.GetBinContent(i)
        - fixed_data_hist.GetBinContent(i)
    )

    max_diff = max(max_diff, diff)

print("Maximum data bin difference =", max_diff)


# ============================================================
# BIN-BY-BIN MC CHANGE
# ============================================================

print("")
print("===== BIN-BY-BIN CORRECTED / OLD MC =====")
print("")

print(
    "{:>4s} {:>12s} {:>12s} {:>14s} {:>14s} {:>12s}".format(
        "bin",
        "low edge",
        "high edge",
        "old MC",
        "fixed MC",
        "fixed/old",
    )
)

print("-" * 78)

h_old = old["mc"]
h_fixed = fixed["mc"]

for i in range(1, h_old.GetNbinsX() + 1):

    lo = h_old.GetXaxis().GetBinLowEdge(i)
    hi = h_old.GetXaxis().GetBinUpEdge(i)

    y_old = h_old.GetBinContent(i)
    y_fixed = h_fixed.GetBinContent(i)

    ratio = y_fixed / y_old if y_old else 0.0

    print(
        "{:4d} {:12.4f} {:12.4f} {:14.4f} {:14.4f} {:12.6f}".format(
            i,
            lo,
            hi,
            y_old,
            y_fixed,
            ratio,
        )
    )


# ============================================================
# PLOT 1:
# DATA + OLD/FIXED V4
# BOTTOM = DATA / MC
# ============================================================

canvas = ROOT.TCanvas("c", "c", 1100, 950)

pad1 = ROOT.TPad(
    "pad1", "pad1",
    0.0, 0.30, 1.0, 1.0
)

pad2 = ROOT.TPad(
    "pad2", "pad2",
    0.0, 0.00, 1.0, 0.30
)

pad1.SetBottomMargin(0.02)
pad1.SetLeftMargin(0.13)
pad1.SetRightMargin(0.05)

pad2.SetTopMargin(0.03)
pad2.SetBottomMargin(0.34)
pad2.SetLeftMargin(0.13)
pad2.SetRightMargin(0.05)

pad1.Draw()
pad2.Draw()


# ---------------- TOP ----------------

pad1.cd()

h_data_plot = old_data_hist.Clone("data_plot")
h_data_plot.SetDirectory(0)

h_data_plot.SetMarkerStyle(20)
h_data_plot.SetMarkerSize(1.0)
h_data_plot.SetLineColor(ROOT.kBlack)
h_data_plot.SetLineWidth(2)

h_data_plot.GetYaxis().SetTitle("Events")
h_data_plot.GetYaxis().SetTitleSize(0.055)
h_data_plot.GetYaxis().SetLabelSize(0.045)
h_data_plot.GetYaxis().SetTitleOffset(1.10)

h_data_plot.GetXaxis().SetLabelSize(0.0)
h_data_plot.GetXaxis().SetTitleSize(0.0)

max_y = max(
    h_data_plot.GetMaximum(),
    h_old.GetMaximum(),
    h_fixed.GetMaximum(),
)

h_data_plot.SetMaximum(1.30 * max_y)
h_data_plot.SetMinimum(0.0)

h_data_plot.Draw("E1")

for cfg in CONFIGS:

    h = results[cfg["tag"]]["mc"]

    h.SetLineColor(cfg["color"])
    h.SetLineStyle(cfg["style"])
    h.SetLineWidth(cfg["width"])
    h.SetMarkerSize(0)

    h.Draw("HIST SAME")

h_data_plot.Draw("E1 SAME")


legend = ROOT.TLegend(
    0.46, 0.66,
    0.92, 0.88
)

legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextSize(0.034)

legend.AddEntry(
    h_data_plot,
    "Data",
    "lep"
)

for cfg in CONFIGS:

    legend.AddEntry(
        results[cfg["tag"]]["mc"],
        cfg["label"],
        "l"
    )

legend.Draw()


latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextSize(0.045)

latex.DrawLatex(
    0.15,
    0.92,
    "LE FHC CC#nu_{#mu}"
)


# ---------------- BOTTOM: DATA/MC ----------------

pad2.cd()

ratio_hists = []

for i, cfg in enumerate(CONFIGS):

    tag = cfg["tag"]

    ratio = make_ratio(
        old_data_hist,
        results[tag]["mc"],
        "data_mc_{}".format(tag)
    )

    ratio.SetLineColor(cfg["color"])
    ratio.SetLineStyle(cfg["style"])
    ratio.SetLineWidth(cfg["width"])
    ratio.SetMarkerSize(0)

    ratio.GetYaxis().SetTitle("Data / MC")
    ratio.GetXaxis().SetTitle("Energy (GeV)")

    ratio.GetYaxis().SetTitleSize(0.10)
    ratio.GetYaxis().SetLabelSize(0.085)
    ratio.GetYaxis().SetTitleOffset(0.52)

    ratio.GetXaxis().SetTitleSize(0.12)
    ratio.GetXaxis().SetLabelSize(0.095)
    ratio.GetXaxis().SetTitleOffset(1.08)

    ratio.GetYaxis().SetNdivisions(505)

    ratio.SetMinimum(0.70)
    ratio.SetMaximum(1.20)

    if i == 0:
        ratio.Draw("HIST")
    else:
        ratio.Draw("HIST SAME")

    ratio_hists.append(ratio)


line = ROOT.TLine(
    old_data_hist.GetXaxis().GetXmin(),
    1.0,
    old_data_hist.GetXaxis().GetXmax(),
    1.0
)

line.SetLineStyle(2)
line.SetLineWidth(2)
line.Draw("SAME")


outname = os.path.join(
    OUT_DIR,
    "CCnumu_EN4_V4_old_vs_corrected_LowQ2.png"
)

canvas.Print(outname)

print("")
print("Saved:", outname)


# ============================================================
# PLOT 2:
# DIRECT FIXED / OLD MC RATIO
# ============================================================

cratio = ROOT.TCanvas(
    "cratio",
    "cratio",
    1100,
    800
)

mc_ratio = make_ratio(
    h_fixed,
    h_old,
    "fixed_over_old"
)

mc_ratio.SetLineColor(ROOT.kBlue + 1)
mc_ratio.SetLineWidth(4)
mc_ratio.SetMarkerSize(0)

mc_ratio.GetYaxis().SetTitle(
    "Corrected V4 / old V4"
)

mc_ratio.GetXaxis().SetTitle(
    "Energy (GeV)"
)

mc_ratio.GetYaxis().SetTitleSize(0.050)
mc_ratio.GetYaxis().SetLabelSize(0.045)
mc_ratio.GetYaxis().SetTitleOffset(1.15)

mc_ratio.GetXaxis().SetTitleSize(0.050)
mc_ratio.GetXaxis().SetLabelSize(0.045)

mc_ratio.SetMinimum(0.75)
mc_ratio.SetMaximum(1.10)

mc_ratio.Draw("HIST")


ratio_line = ROOT.TLine(
    h_old.GetXaxis().GetXmin(),
    1.0,
    h_old.GetXaxis().GetXmax(),
    1.0
)

ratio_line.SetLineStyle(2)
ratio_line.SetLineWidth(2)
ratio_line.Draw("SAME")


latex_ratio = ROOT.TLatex()
latex_ratio.SetNDC()
latex_ratio.SetTextSize(0.040)

latex_ratio.DrawLatex(
    0.15,
    0.92,
    "LE FHC CC#nu_{#mu}: LowQ2 implementation effect"
)


ratio_outname = os.path.join(
    OUT_DIR,
    "CCnumu_EN4_V4_corrected_over_old_LowQ2.png"
)

cratio.Print(ratio_outname)

print("Saved:", ratio_outname)


# ============================================================
# PLOT 3:
# SHAPE-ONLY CORRECTED / OLD
# ============================================================

cshape = ROOT.TCanvas(
    "cshape",
    "cshape",
    1100,
    800
)

old_shape = h_old.Clone("old_shape")
fixed_shape = h_fixed.Clone("fixed_shape")

old_shape.SetDirectory(0)
fixed_shape.SetDirectory(0)

old_int = integral(old_shape)
fixed_int = integral(fixed_shape)

if old_int:
    old_shape.Scale(1.0 / old_int)

if fixed_int:
    fixed_shape.Scale(1.0 / fixed_int)


shape_ratio = make_ratio(
    fixed_shape,
    old_shape,
    "shape_fixed_over_old"
)

shape_ratio.SetLineColor(ROOT.kBlue + 1)
shape_ratio.SetLineWidth(4)
shape_ratio.SetMarkerSize(0)

shape_ratio.GetYaxis().SetTitle(
    "Normalized corrected / old V4"
)

shape_ratio.GetXaxis().SetTitle(
    "Energy (GeV)"
)

shape_ratio.GetYaxis().SetTitleSize(0.050)
shape_ratio.GetYaxis().SetLabelSize(0.045)
shape_ratio.GetYaxis().SetTitleOffset(1.15)

shape_ratio.GetXaxis().SetTitleSize(0.050)
shape_ratio.GetXaxis().SetLabelSize(0.045)

shape_ratio.SetMinimum(0.85)
shape_ratio.SetMaximum(1.15)

shape_ratio.Draw("HIST")


shape_line = ROOT.TLine(
    old_shape.GetXaxis().GetXmin(),
    1.0,
    old_shape.GetXaxis().GetXmax(),
    1.0
)

shape_line.SetLineStyle(2)
shape_line.SetLineWidth(2)
shape_line.Draw("SAME")


latex_shape = ROOT.TLatex()
latex_shape.SetNDC()
latex_shape.SetTextSize(0.040)

latex_shape.DrawLatex(
    0.15,
    0.92,
    "LE FHC CC#nu_{#mu}: LowQ2 shape effect"
)


shape_outname = os.path.join(
    OUT_DIR,
    "CCnumu_EN4_V4_old_vs_corrected_LowQ2_shape.png"
)

cshape.Print(shape_outname)

print("Saved:", shape_outname)
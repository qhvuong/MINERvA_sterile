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

CONFIGS = [
    {
        "tag": "DebugCuts_CCnumu_oldFV_MnvTunesV1LECOH",
        "label": "V1 + LE COH",
        "color": ROOT.kBlack,
        "style": 1,
        "width": 4,
    },
    {
        "tag": "CCnumu_V1_LECOH_JOINT",
        "label": "V1 + LE COH + JOINT",
        "color": ROOT.kBlue + 1,
        "style": 2,
        "width": 4,
    },
    {
        "tag": "DebugCuts_CCnumu_oldFV_MnvTunesV1LowQ2",
        "label": "V1 + MENU1PI",
        "color": ROOT.kRed + 1,
        "style": 3,
        "width": 4,
    },
    {
        "tag": "CCnumu_V1D2_LECOH_DiffCOH_noLowQ2",
        "label": "V1 + D2 + LE COH + DiffCOH, no LowQ2",
        "color": ROOT.kGreen + 2,
        "style": 4,
        "width": 4,
    },
    {
        "tag": "DebugCuts_CCnumu_oldFV_MnvTunesV4",
        "label": "Current / V4",
        "color": ROOT.kMagenta + 1,
        "style": 5,
        "width": 4,
    },
]

REF_TAG = "DebugCuts_CCnumu_oldFV_MnvTunesV1LECOH"


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
                name, root_file.GetName()
            )
        )

    # CV-only comparison: convert MnvH1D -> ordinary TH1D.
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
    # Include underflow and overflow, same convention as previous script.
    return h.Integral(0, h.GetNbinsX() + 1)


def make_ratio(num, den, name):
    ratio = num.Clone(name)
    ratio.SetDirectory(0)

    ratio.Divide(
        num,
        den,
        1.0,
        1.0,
        ""
    )

    return ratio


# ============================================================
# LOAD EVERYTHING
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
        raise RuntimeError("Could not open {}".format(data_path(tag)))

    if not fm or fm.IsZombie():
        raise RuntimeError("Could not open {}".format(mc_path(tag)))

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

ref = results[REF_TAG]

ref_data_yield = ref["data_yield"]
ref_mc_yield = ref["mc_yield"]

print("")
print("====================================================================")
print("CCnumu tune comparison: {}".format(HIST_NAME))
print("Reference = V1 + LE COH")
print("====================================================================")
print("")

print(
    "{:<46s} {:>14s} {:>14s} {:>12s} {:>18s}".format(
        "Configuration",
        "Data",
        "MC",
        "Data/MC",
        "vs V1+LECOH [%]",
    )
)

print("-" * 110)

for cfg in CONFIGS:

    tag = cfg["tag"]
    r = results[tag]

    data_yield = r["data_yield"]
    mc_yield = r["mc_yield"]

    data_mc = data_yield / mc_yield if mc_yield else 0.0

    shift_vs_ref = (
        (mc_yield / ref_mc_yield - 1.0) * 100.0
        if ref_mc_yield else 0.0
    )

    print(
        "{:<46s} {:>14.3f} {:>14.3f} {:>12.6f} {:>18.3f}".format(
            cfg["label"],
            data_yield,
            mc_yield,
            data_mc,
            shift_vs_ref,
        )
    )


# ============================================================
# DIRECT JOINT EFFECT
# ============================================================

joint_tag = "CCnumu_V1_LECOH_JOINT"
old_tag = "CCnumu_V1D2_LECOH_DiffCOH_noLowQ2"

joint_mc = results[joint_tag]["mc_yield"]
old_mc = results[old_tag]["mc_yield"]

print("")
print("===== TUNE EFFECTS RELATIVE TO V1 + LE COH =====")

print(
    "JOINT effect = {:+.3f}%".format(
        (joint_mc / ref_mc_yield - 1.0) * 100.0
    )
)

print(
    "D2 + DiffCOH, no LowQ2 effect = {:+.3f}%".format(
        (old_mc / ref_mc_yield - 1.0) * 100.0
    )
)

print(
    "JOINT vs old no-LowQ2 config = {:+.3f}%".format(
        (joint_mc / old_mc - 1.0) * 100.0
    )
)


# ============================================================
# CHECK THAT DATA IS IDENTICAL
# ============================================================

print("")
print("===== DATA CONSISTENCY CHECK =====")

ref_data = ref["data"]

for cfg in CONFIGS:

    tag = cfg["tag"]

    if tag == REF_TAG:
        continue

    h = results[tag]["data"]

    max_diff = 0.0

    for i in range(0, h.GetNbinsX() + 2):

        diff = abs(
            h.GetBinContent(i) -
            ref_data.GetBinContent(i)
        )

        max_diff = max(max_diff, diff)

    print(
        "{:<44s} max bin difference = {}".format(
            cfg["label"],
            max_diff
        )
    )


# ============================================================
# PLOT
# ============================================================

canvas = ROOT.TCanvas("c", "c", 1100, 950)

pad1 = ROOT.TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0)
pad2 = ROOT.TPad("pad2", "pad2", 0.0, 0.00, 1.0, 0.30)

pad1.SetBottomMargin(0.02)
pad1.SetLeftMargin(0.13)
pad1.SetRightMargin(0.05)

pad2.SetTopMargin(0.03)
pad2.SetBottomMargin(0.34)
pad2.SetLeftMargin(0.13)
pad2.SetRightMargin(0.05)

pad1.Draw()
pad2.Draw()


# ============================================================
# TOP: SPECTRA
# ============================================================

pad1.cd()

h_data_plot = ref_data.Clone("data_plot")

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
    [h_data_plot.GetMaximum()] +
    [results[cfg["tag"]]["mc"].GetMaximum() for cfg in CONFIGS]
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


legend = ROOT.TLegend(0.48, 0.57, 0.92, 0.88)

legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextSize(0.034)

legend.AddEntry(h_data_plot, "Data", "lep")

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


# ============================================================
# BOTTOM: DATA / MC
# ============================================================

pad2.cd()

ratio_hists = []

for i, cfg in enumerate(CONFIGS):

    tag = cfg["tag"]

    ratio = make_ratio(
        ref_data,
        results[tag]["mc"],
        "ratio_{}".format(tag)
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
    ref_data.GetXaxis().GetXmin(),
    1.0,
    ref_data.GetXaxis().GetXmax(),
    1.0
)

line.SetLineStyle(2)
line.SetLineWidth(2)
line.Draw("SAME")


# ============================================================
# SHAPE-ONLY COMPARISON
# Normalize every MC histogram to unit area, then divide by
# normalized V1 + LECOH.
# ============================================================

cshape = ROOT.TCanvas("cshape", "cshape", 1100, 800)

ref_shape = results[REF_TAG]["mc"].Clone("ref_shape")
ref_shape.SetDirectory(0)

ref_int = integral(ref_shape)

if ref_int:
    ref_shape.Scale(1.0 / ref_int)


shape_ratios = []

for i, cfg in enumerate(CONFIGS):

    tag = cfg["tag"]

    h = results[tag]["mc"].Clone(
        "shape_{}".format(tag)
    )
    h.SetDirectory(0)

    h_int = integral(h)

    if h_int:
        h.Scale(1.0 / h_int)

    ratio = h.Clone(
        "shape_ratio_{}".format(tag)
    )

    ratio.Divide(ref_shape)

    ratio.SetLineColor(cfg["color"])
    ratio.SetLineStyle(cfg["style"])
    ratio.SetLineWidth(cfg["width"])
    ratio.SetMarkerSize(0)

    ratio.GetYaxis().SetTitle(
        "Normalized MC / V1+LECOH"
    )
    ratio.GetXaxis().SetTitle(
        "Energy (GeV)"
    )

    ratio.GetYaxis().SetTitleSize(0.050)
    ratio.GetYaxis().SetLabelSize(0.045)
    ratio.GetYaxis().SetTitleOffset(1.15)

    ratio.GetXaxis().SetTitleSize(0.050)
    ratio.GetXaxis().SetLabelSize(0.045)

    ratio.SetMinimum(0.80)
    ratio.SetMaximum(1.20)

    if i == 0:
        ratio.Draw("HIST")
    else:
        ratio.Draw("HIST SAME")

    shape_ratios.append(ratio)


line_shape = ROOT.TLine(
    ref_shape.GetXaxis().GetXmin(),
    1.0,
    ref_shape.GetXaxis().GetXmax(),
    1.0
)

line_shape.SetLineStyle(2)
line_shape.SetLineWidth(2)
line_shape.Draw("SAME")


leg_shape = ROOT.TLegend(0.49, 0.63, 0.91, 0.88)

leg_shape.SetBorderSize(0)
leg_shape.SetFillStyle(0)
leg_shape.SetTextSize(0.032)

for i, cfg in enumerate(CONFIGS):
    leg_shape.AddEntry(
        shape_ratios[i],
        cfg["label"],
        "l"
    )

leg_shape.Draw()


latex_shape = ROOT.TLatex()
latex_shape.SetNDC()
latex_shape.SetTextSize(0.040)

latex_shape.DrawLatex(
    0.15,
    0.92,
    "LE FHC CC#nu_{#mu}: shape only"
)


shape_outname = (
    "/exp/minerva/data/users/qvuong/studies/CCnumu_EN4_MnvTune_shape_only_comparison.png"
)

cshape.Print(shape_outname)

print("Saved:", shape_outname)



# ============================================================
# SAVE
# ============================================================

outname = "/exp/minerva/data/users/qvuong/studies/CCnumu_EN4_V1LECOH_JOINT_comparison.png"

canvas.Print(outname)

print("")
print("Saved:", outname)
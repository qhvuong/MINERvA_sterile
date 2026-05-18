import os
import copy
import ROOT
import PlotUtils
import numpy as np

from tools.StitchedHistogram import *
from tools.Helper import *

from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS
from config.AnalysisConfig import AnalysisConfig
from tools import Utilities
from tools.PlotLibrary import HistHolder

ccnueroot = os.environ.get("CCNUEROOT")

MNVPLOTTER = PlotUtils.MnvPlotter()
MNVPLOTTER.error_summary_group_map.clear()
for k, v in CONSOLIDATED_ERROR_GROUPS.items():
    vec = ROOT.vector("std::string")()
    for vs in v:
        vec.push_back(vs)
    MNVPLOTTER.error_summary_group_map[k] = vec

errsToRemove = ["LowQ2Pi","elETracker"]
ROOT.TH1.AddDirectory(False)


def clone_total(holder, name):
    h = holder.hists["Total"].Clone(name)
    h.SetDirectory(0)
    return h

def get_hist_checked(root_file, names, label):
    """
    Try several possible histogram names and return the first one found.
    """
    for name in names:
        h = root_file.Get(name)
        if h:
            h = h.Clone(label)
            h.SetDirectory(0)
            print("Loaded {} from {}".format(label, name))
            return h

    print("\nAvailable keys containing Eel/EReco/Elastic:")
    for key in root_file.GetListOfKeys():
        kname = key.GetName()
        if ("Eel" in kname) or ("EReco" in kname) or ("Elastic" in kname):
            print("  ", kname)

    raise RuntimeError("Could not find {}. Tried: {}".format(label, names))

def sum_hists_from_file(root_file, hist_names, out_name):
    hsum = None

    for hist_name in hist_names:
        h = root_file.Get(hist_name)
        if not h:
            print("WARNING: missing", hist_name)
            continue

        h = h.Clone(out_name + "_" + hist_name)
        h.SetDirectory(0)

        if hsum is None:
            hsum = h.Clone(out_name)
            hsum.SetDirectory(0)
            hsum.Reset()

        hsum.Add(h)

    if hsum is None:
        raise RuntimeError("Could not build {}".format(out_name))

    return hsum

def check_bad_bins(h, label):
    print("\nChecking", label)
    print("Integral:", h.Integral())
    for i in range(0, h.GetNbinsX() + 2):
        v = h.GetBinContent(i)
        e = h.GetBinError(i)
        if not np.isfinite(v) or not np.isfinite(e):
            print("  BAD bin", i, "content =", v, "error =", e)

def print_hist_bins(h, label, max_bins=20):
    print("\n===== {} =====".format(label))
    print("Integral =", h.Integral())
    print("Nbins    =", h.GetNbinsX())
    for i in range(1, min(h.GetNbinsX(), max_bins) + 1):
        print(
            "  bin {:2d}: low={:8.3f} high={:8.3f} content={:12.6g} error={:12.6g}".format(
                i,
                h.GetBinLowEdge(i),
                h.GetBinLowEdge(i) + h.GetBinWidth(i),
                h.GetBinContent(i),
                h.GetBinError(i),
            )
        )

def print_pot_scale_check(label, data_pot, mc_pot, standPOT):
    print("\n===== POT SCALE CHECK: {} =====".format(label))
    print("  data_pot =", data_pot)
    print("  mc_pot   =", mc_pot)
    print("  standPOT =", standPOT)

    if mc_pot is None or mc_pot == 0:
        print("  MC scale = BAD")
    else:
        print("  MC scale = standPOT / mc_pot =", standPOT / mc_pot)

    if data_pot is None or data_pot == 0:
        print("  Data scale = no explicit scaling expected")
    else:
        print("  Data scale = standPOT / data_pot =", standPOT / data_pot)

def print_integral_change(label, h_raw, h_scaled):
    print("\n===== INTEGRAL CHANGE: {} =====".format(label))
    print("  raw integral    =", h_raw.Integral())
    print("  scaled integral =", h_scaled.Integral())

    if h_raw.Integral() != 0:
        print("  scaled/raw      =", h_scaled.Integral() / h_raw.Integral())
    else:
        print("  scaled/raw      = undefined")

def compare_1d_2d_template(label, h_1d, h_2d, project_axis="x"):
    """
    Compare a 1D spectrum integral to the integral of a 2D template.
    For your templates, total TH2 integral should usually match the 1D
    component integral after POT scaling, unless content is outside the
    template axis range or under/overflow is treated differently.
    """
    print("\n===== 1D vs 2D TEMPLATE CHECK: {} =====".format(label))
    print("  1D integral =", h_1d.Integral())
    print("  2D integral =", h_2d.Integral())

    diff = h_1d.Integral() - h_2d.Integral()
    print("  1D - 2D     =", diff)

    if h_1d.Integral() != 0:
        print("  frac diff   =", diff / h_1d.Integral())

    # Check under/overflow-inclusive integral if available.
    try:
        nx = h_2d.GetNbinsX()
        ny = h_2d.GetNbinsY()
        total_with_flow = h_2d.Integral(0, nx + 1, 0, ny + 1)
        print("  2D integral including under/overflow =", total_with_flow)
        print("  1D - 2D incl flow =", h_1d.Integral() - total_with_flow)
    except Exception as e:
        print("  Could not compute 2D flow-inclusive integral:", e)




def load_nue_elastic_fhc():
    type_path_map = {
        "data": "/exp/minerva/data/users/qvuong/elastic_nue/kin_dist_dataleFHC_NuE_allSystematics_fixedMnvTunes_MAD.root",
        "mc":   "/exp/minerva/data/users/qvuong/elastic_nue/kin_dist_mcleFHC_NuE_allSystematics_fixedMnvTunes_MAD.root",
    }

    data_file, mc_file, pot_scale, data_pot, mc_pot = Utilities.getFilesAndPOTScale(
        "NuE_allSystematics_fixedMnvTunes", type_path_map, "MAD", True
    )

    standPOT = data_pot if data_pot is not None else mc_pot
    binwidthScale = getattr(AnalysisConfig, "binwidth", False)

    print_pot_scale_check("FHC NuE elastic", data_pot, mc_pot, standPOT)

    tuned_file = ROOT.TFile.Open(
        "/exp/minerva/data/users/qvuong/elastic_nue/bkgfit_leFHC_nueElastic_matrix_NuE_allSystematics_fixedMnvTunes_MAD.root"
    )

    if not tuned_file or tuned_file.IsZombie():
        raise RuntimeError("Could not open FHC NuE elastic bkgfit file.")

    h_mc = tuned_file.Get("Eel_predicted_Signal")
    h_data = tuned_file.Get("Eel_data_bkgSubbed")

    if not h_mc or not h_data:
        raise RuntimeError("Could not load NuE elastic tuned histograms from bkgfit file.")

    h_mc = h_mc.Clone("fhc_elastic_mc")
    h_data = h_data.Clone("fhc_elastic_data")
    h_mc.SetDirectory(0)
    h_data.SetDirectory(0)

    print_hist_bins(h_mc,   "FHC elastic tuned MC from bkgfit")
    print_hist_bins(h_data, "FHC elastic bkg-subtracted data from bkgfit")

    # -------------------------------------------------
    # L/E template from original elastic MC production
    # -------------------------------------------------
    template_holder = HistHolder(
        "Reco Lepton Energy vs L/E",
        mc_file,
        "Signal",
        True,
        mc_pot,
        standPOT
    )

    h_template_raw = clone_total(template_holder, "fhc_elastic_template_total_raw")
    print_hist_bins(
        h_template_raw.ProjectionX("fhc_elastic_template_raw_projx"),
        "FHC elastic template RAW projection before POTScale"
    )

    template_holder.POTScale(binwidthScale)

    h_template_total = clone_total(template_holder, "fhc_elastic_template_total")
    h_template_total.SetDirectory(0)

    print_hist_bins(
        h_template_total.ProjectionX("fhc_elastic_template_scaled_projx"),
        "FHC elastic template after POTScale projection"
    )

    print_integral_change(
        "FHC elastic total template",
        h_template_raw,
        h_template_total
    )

    # -------------------------------------------------
    # 1D flavor components for oscillation bookkeeping
    # -------------------------------------------------
    h_electron = get_hist_checked(
        mc_file,
        ["Eel_NuEElasticE"],
        "electron_fhc_elastic"
    )

    h_muon = get_hist_checked(
        mc_file,
        ["Eel_NuEElasticMu"],
        "muon_fhc_elastic"
    )

    if mc_pot is None or mc_pot == 0:
        raise RuntimeError("Bad mc_pot for FHC NuE elastic: {}".format(mc_pot))

    scale = standPOT / mc_pot

    print("\nFHC elastic manual POT info:")
    print("  data_pot =", data_pot)
    print("  mc_pot   =", mc_pot)
    print("  standPOT =", standPOT)
    print("  scale    =", scale)

    h_electron_raw = h_electron.Clone("electron_fhc_elastic_raw")
    h_muon_raw = h_muon.Clone("muon_fhc_elastic_raw")
    h_electron_raw.SetDirectory(0)
    h_muon_raw.SetDirectory(0)

    print_hist_bins(h_electron_raw, "FHC elastic electron RAW before manual POT scale")
    print_hist_bins(h_muon_raw,     "FHC elastic muon RAW before manual POT scale")

    h_electron.Scale(scale)
    h_muon.Scale(scale)

    print_hist_bins(h_electron, "FHC elastic electron after manual POT scale")
    print_hist_bins(h_muon,     "FHC elastic muon after manual POT scale")

    print_integral_change("FHC elastic electron flavor", h_electron_raw, h_electron)
    print_integral_change("FHC elastic muon flavor",     h_muon_raw,     h_muon)

    # Match electron+muon components to tuned prediction.
    h_sum = h_electron.Clone("fhc_elastic_flavor_sum")
    h_sum.SetDirectory(0)
    h_sum.Add(h_muon)

    if h_sum.Integral() > 0:
        flavor_scale = h_mc.Integral() / h_sum.Integral()
        print("FHC elastic flavor scale to tuned prediction =", flavor_scale)

        h_electron.Scale(flavor_scale)
        h_muon.Scale(flavor_scale)

        print_hist_bins(h_electron, "FHC elastic electron after flavor_scale")
        print_hist_bins(h_muon,     "FHC elastic muon after flavor_scale")
    else:
        raise RuntimeError("FHC elastic electron+muon flavor sum is zero.")

    # -------------------------------------------------
    # 2D flavor templates
    # -------------------------------------------------
    h_template_electron = get_hist_checked(
        mc_file,
        ["ElepReco_LE_NuEElasticE"],
        "fhc_elastic_template_electron"
    )
    h_template_electron.Scale(scale)

    h_template_muon = get_hist_checked(
        mc_file,
        ["ElepReco_LE_NuEElasticMu"],
        "fhc_elastic_template_muon"
    )
    h_template_muon.Scale(scale)

    compare_1d_2d_template(
        "FHC elastic electron",
        h_electron,
        h_template_electron
    )

    compare_1d_2d_template(
        "FHC elastic muon",
        h_muon,
        h_template_muon
    )

    print("\nFHC NuE elastic check:")
    print("  tuned mc integral       =", h_mc.Integral())
    print("  tuned data integral     =", h_data.Integral())
    print("  electron component int  =", h_electron.Integral())
    print("  muon component int      =", h_muon.Integral())
    print("  electron+muon int       =", h_electron.Integral() + h_muon.Integral())
    print("  template electron int   =", h_template_electron.Integral())
    print("  template muon int       =", h_template_muon.Integral())

    return {
        "mc": h_mc,
        "data": h_data,
        "electron": h_electron,
        "muon": h_muon,
        "template_electron": h_template_electron,
        "template_muon": h_template_muon,
    }

def load_ccnue_fhc():
    type_path_map = {
        "data": "/exp/minerva/data/users/qvuong/nu_e/kin_dist_dataleFHC_CCnue_allSystematics_fixedMnvTunes_MAD.root",
        "mc":   "/exp/minerva/data/users/qvuong/nu_e/kin_dist_mcleFHC_CCnue_allSystematics_fixedMnvTunes_MAD.root",
    }

    data_file, mc_file, pot_scale, data_pot, mc_pot = Utilities.getFilesAndPOTScale(
        "CCnue_allSystematics_fixedMnvTunes", type_path_map, "MAD", True
    )

    standPOT = data_pot if data_pot is not None else mc_pot
    binwidthScale = getattr(AnalysisConfig, "binwidth", False)

    print_pot_scale_check("FHC CCnue", data_pot, mc_pot, standPOT)

    tuned_file = ROOT.TFile.Open(
        "/exp/minerva/data/users/qvuong/nu_e/bkgfit_leFHC_N4_tune_CCnue_allSystematics_fixedMnvTunes_MAD.root"
    )

    if not tuned_file or tuned_file.IsZombie():
        raise RuntimeError("Could not open FHC CCnue bkgfit file.")

    h_mc = tuned_file.Get("EN4_predicted_Signal")
    h_data = tuned_file.Get("EN4_data_bkgSubbed")

    if not h_mc or not h_data:
        raise RuntimeError("Could not load CCnue tuned histograms from bkgfit file.")

    h_mc = h_mc.Clone("fhc_ccnue_mc")
    h_data = h_data.Clone("fhc_ccnue_data")
    h_mc.SetDirectory(0)
    h_data.SetDirectory(0)

    print_hist_bins(h_mc,   "FHC CCnue tuned MC from bkgfit")
    print_hist_bins(h_data, "FHC CCnue bkg-subtracted data from bkgfit")

    # Nominal nue L/E template.
    template_file = ROOT.TFile.Open(
        "/exp/minerva/data/users/qvuong/nu_e/kin_dist_mcleFHC_CCnue_allSystematics_fixedMnvTunes_MAD.root"
    )

    template_holder = HistHolder(
        "Reco Energy vs L/E",
        template_file,
        "Signal",
        True,
        mc_pot,
        standPOT
    )

    h_template_raw = clone_total(template_holder, "fhc_ccnue_template_raw")
    template_holder.POTScale(binwidthScale)
    h_template = clone_total(template_holder, "fhc_ccnue_template")
    h_template.SetDirectory(0)

    print_integral_change("FHC CCnue template", h_template_raw, h_template)

    # Swapped sample.
    swap_type_path_map = {
        "mc": "/exp/minerva/data/users/qvuong/nu_e_swapped/kin_dist_mcleFHC_CCnueswap_allSystematics_fixedMnvTunes_MAD.root"
    }

    _, swap_mc_file, _, _, swap_mc_pot = Utilities.getFilesAndPOTScale(
        "CCnueswap_allSystematics_fixedMnvTunes", swap_type_path_map, "MAD", True
    )

    print_pot_scale_check("FHC CCnue swap", None, swap_mc_pot, standPOT)

    swap_file = ROOT.TFile.Open(
        "/exp/minerva/data/users/qvuong/nu_e_swapped/kin_dist_mcleFHC_CCnueswap_allSystematics_fixedMnvTunes_MAD.root"
    )

    swap_template_holder = HistHolder(
        "Reco Energy vs L/E",
        swap_file,
        "Signal",
        True,
        swap_mc_pot,
        standPOT
    )

    swap_hist_holder = HistHolder(
        "Biased Neutrino Energy",
        swap_file,
        "Signal",
        True,
        swap_mc_pot,
        standPOT
    )

    h_swap_raw = clone_total(swap_hist_holder, "fhc_ccnue_swap_raw")
    h_swap_template_raw = clone_total(swap_template_holder, "fhc_ccnue_swap_template_raw")

    swap_hist_holder.POTScale(binwidthScale)
    swap_template_holder.POTScale(binwidthScale)

    h_swap = clone_total(swap_hist_holder, "fhc_ccnue_swap")
    h_swap_template = clone_total(swap_template_holder, "fhc_ccnue_swap_template")

    h_swap.SetDirectory(0)
    h_swap_template.SetDirectory(0)

    print_integral_change("FHC CCnue swap 1D", h_swap_raw, h_swap)
    print_integral_change("FHC CCnue swap template", h_swap_template_raw, h_swap_template)

    compare_1d_2d_template(
        "FHC CCnue nominal",
        h_mc,
        h_template
    )

    compare_1d_2d_template(
        "FHC CCnue swap",
        h_swap,
        h_swap_template
    )

    return {
        "mc": h_mc,
        "data": h_data,
        "swap": h_swap,
        "template_nue": h_template,
        "template_swap": h_swap_template,
    }


def load_ccnuebar_rhc():
    type_path_map = {
        "data": "/exp/minerva/data/users/qvuong/antinu_e/kin_dist_datale5_CCnuebar_allSystematics_fixedMnvTunes_MAD.root",
        "mc":   "/exp/minerva/data/users/qvuong/antinu_e/kin_dist_mcle5_CCnuebar_allSystematics_fixedMnvTunes_MAD.root",
    }

    data_file, mc_file, pot_scale, data_pot, mc_pot = Utilities.getFilesAndPOTScale(
        "CCnuebar_allSystematics_fixedMnvTunes", type_path_map, "MAD", True
    )

    standPOT = data_pot if data_pot is not None else mc_pot
    binwidthScale = getattr(AnalysisConfig, "binwidth", False)

    print_pot_scale_check("RHC CCnuebar", data_pot, mc_pot, standPOT)

    tuned_file = ROOT.TFile.Open(
        "/exp/minerva/data/users/qvuong/antinu_e/bkgfit_le5_N4_tune_CCnuebar_allSystematics_fixedMnvTunes_MAD.root"
    )

    if not tuned_file or tuned_file.IsZombie():
        raise RuntimeError("Could not open RHC CCnuebar bkgfit file.")

    h_mc = tuned_file.Get("EN4_predicted_Signal")
    h_data = tuned_file.Get("EN4_data_bkgSubbed")

    if not h_mc or not h_data:
        raise RuntimeError("Could not load CCnuebar tuned histograms from bkgfit file.")

    h_mc = h_mc.Clone("rhc_ccnuebar_mc")
    h_data = h_data.Clone("rhc_ccnuebar_data")
    h_mc.SetDirectory(0)
    h_data.SetDirectory(0)

    print_hist_bins(h_mc,   "RHC CCnuebar tuned MC from bkgfit")
    print_hist_bins(h_data, "RHC CCnuebar bkg-subtracted data from bkgfit")

    # Nominal nuebar L/E template.
    template_file = ROOT.TFile.Open(
        "/exp/minerva/data/users/qvuong/antinu_e/kin_dist_mcle5_CCnuebar_allSystematics_fixedMnvTunes_MAD.root"
    )

    template_holder = HistHolder(
        "Reco Energy vs L/E",
        template_file,
        "Signal",
        True,
        mc_pot,
        standPOT
    )

    h_template_raw = clone_total(template_holder, "rhc_ccnuebar_template_raw")
    template_holder.POTScale(binwidthScale)
    h_template = clone_total(template_holder, "rhc_ccnuebar_template")
    h_template.SetDirectory(0)

    print_integral_change("RHC CCnuebar template", h_template_raw, h_template)

    # Swapped sample.
    swap_type_path_map = {
        "mc": "/exp/minerva/data/users/qvuong/antinu_e_swapped/kin_dist_mcle5_CCnuebarswap_allSystematics_fixedMnvTunes_MAD.root"
    }

    _, swap_mc_file, _, _, swap_mc_pot = Utilities.getFilesAndPOTScale(
        "CCnuebarswap_allSystematics_fixedMnvTunes", swap_type_path_map, "MAD", True
    )

    print_pot_scale_check("RHC CCnuebar swap", None, swap_mc_pot, standPOT)

    swap_file = ROOT.TFile.Open(
        "/exp/minerva/data/users/qvuong/antinu_e_swapped/kin_dist_mcle5_CCnuebarswap_allSystematics_fixedMnvTunes_MAD.root"
    )

    swap_template_holder = HistHolder(
        "Reco Energy vs L/E",
        swap_file,
        "Signal",
        True,
        swap_mc_pot,
        standPOT
    )

    swap_hist_holder = HistHolder(
        "Biased Neutrino Energy",
        swap_file,
        "Signal",
        True,
        swap_mc_pot,
        standPOT
    )

    h_swap_raw = clone_total(swap_hist_holder, "rhc_ccnuebar_swap_raw")
    h_swap_template_raw = clone_total(swap_template_holder, "rhc_ccnuebar_swap_template_raw")

    swap_hist_holder.POTScale(binwidthScale)
    swap_template_holder.POTScale(binwidthScale)

    h_swap = clone_total(swap_hist_holder, "rhc_ccnuebar_swap")
    h_swap_template = clone_total(swap_template_holder, "rhc_ccnuebar_swap_template")

    h_swap.SetDirectory(0)
    h_swap_template.SetDirectory(0)

    print_integral_change("RHC CCnuebar swap 1D", h_swap_raw, h_swap)
    print_integral_change("RHC CCnuebar swap template", h_swap_template_raw, h_swap_template)

    compare_1d_2d_template(
        "RHC CCnuebar nominal",
        h_mc,
        h_template
    )

    compare_1d_2d_template(
        "RHC CCnuebar swap",
        h_swap,
        h_swap_template
    )

    return {
        "mc": h_mc,
        "data": h_data,
        "swap": h_swap,
        "template_nue": h_template,
        "template_swap": h_swap_template,
    }


def load_ccnumu_fhc():
    type_path_map = {
        "data": "/exp/minerva/data/users/qvuong/nu_mu/kin_dist_dataleFHC_CCnumu_allSystematics_fixedMnvTunes_MAD.root",
        "mc":   "/exp/minerva/data/users/qvuong/nu_mu/kin_dist_mcleFHC_CCnumu_allSystematics_fixedMnvTunes_MAD.root",
    }

    data_file, mc_file, pot_scale, data_pot, mc_pot = Utilities.getFilesAndPOTScale(
        "CCnumu_allSystematics_fixedMnvTunes", type_path_map, "MAD", True
    )

    standPOT = data_pot if data_pot is not None else mc_pot
    binwidthScale = getattr(AnalysisConfig, "binwidth", False)

    print_pot_scale_check("FHC CCnumu", data_pot, mc_pot, standPOT)

    mc_holder = HistHolder(
        "Biased Neutrino Energy",
        mc_file,
        "Signal",
        True,
        mc_pot,
        standPOT
    )

    data_holder = HistHolder(
        "Biased Neutrino Energy",
        data_file,
        "Signal",
        False,
        data_pot,
        standPOT
    )

    template_holder = HistHolder(
        "Reco Energy vs L/E",
        mc_file,
        "Signal",
        True,
        mc_pot,
        standPOT
    )

    h_mc_raw = clone_total(mc_holder, "fhc_ccnumu_mc_raw")
    h_data_raw = data_holder.GetHist().Clone("fhc_ccnumu_data_raw")
    h_template_raw = clone_total(template_holder, "fhc_ccnumu_template_raw")

    h_data_raw.SetDirectory(0)

    mc_holder.POTScale(binwidthScale)
    data_holder.POTScale(binwidthScale)
    template_holder.POTScale(binwidthScale)

    h_mc = clone_total(mc_holder, "fhc_ccnumu_mc")
    h_data = data_holder.GetHist().Clone("fhc_ccnumu_data")
    h_template = clone_total(template_holder, "fhc_ccnumu_template")

    h_data.SetDirectory(0)
    h_template.SetDirectory(0)

    print_integral_change("FHC CCnumu MC", h_mc_raw, h_mc)
    print_integral_change("FHC CCnumu data", h_data_raw, h_data)
    print_integral_change("FHC CCnumu template", h_template_raw, h_template)

    compare_1d_2d_template(
        "FHC CCnumu",
        h_mc,
        h_template
    )

    return {
        "mc": h_mc,
        "data": h_data,
        "template_numu": h_template,
    }


def load_ccnumubar_rhc():
    type_path_map = {
        "data": "/exp/minerva/data/users/qvuong/antinu_mu/kin_dist_datale5_CCnumubar_allSystematics_fixedMnvTunes_MAD.root",
        "mc":   "/exp/minerva/data/users/qvuong/antinu_mu/kin_dist_mcle5_CCnumubar_allSystematics_fixedMnvTunes_MAD.root",
    }

    data_file, mc_file, pot_scale, data_pot, mc_pot = Utilities.getFilesAndPOTScale(
        "CCnumubar_allSystematics_fixedMnvTunes", type_path_map, "MAD", True
    )

    standPOT = data_pot if data_pot is not None else mc_pot
    binwidthScale = getattr(AnalysisConfig, "binwidth", False)

    print_pot_scale_check("RHC CCnumubar", data_pot, mc_pot, standPOT)

    mc_holder = HistHolder(
        "Biased Neutrino Energy",
        mc_file,
        "Signal",
        True,
        mc_pot,
        standPOT
    )

    data_holder = HistHolder(
        "Biased Neutrino Energy",
        data_file,
        "Signal",
        False,
        data_pot,
        standPOT
    )

    template_holder = HistHolder(
        "Reco Energy vs L/E",
        mc_file,
        "Signal",
        True,
        mc_pot,
        standPOT
    )

    h_mc_raw = clone_total(mc_holder, "rhc_ccnumubar_mc_raw")
    h_data_raw = data_holder.GetHist().Clone("rhc_ccnumubar_data_raw")
    h_template_raw = clone_total(template_holder, "rhc_ccnumubar_template_raw")

    h_data_raw.SetDirectory(0)

    mc_holder.POTScale(binwidthScale)
    data_holder.POTScale(binwidthScale)
    template_holder.POTScale(binwidthScale)

    h_mc = clone_total(mc_holder, "rhc_ccnumubar_mc")
    h_data = data_holder.GetHist().Clone("rhc_ccnumubar_data")
    h_template = clone_total(template_holder, "rhc_ccnumubar_template")

    h_mc.SetDirectory(0)
    h_data.SetDirectory(0)
    h_template.SetDirectory(0)

    print_integral_change("RHC CCnumubar MC", h_mc_raw, h_mc)
    print_integral_change("RHC CCnumubar data", h_data_raw, h_data)
    print_integral_change("RHC CCnumubar template", h_template_raw, h_template)

    compare_1d_2d_template(
        "RHC CCnumubar",
        h_mc,
        h_template
    )

    return {
        "mc": h_mc,
        "data": h_data,
        "template_numu": h_template,
    }

if __name__ == "__main__":
    binwidthScale = getattr(AnalysisConfig, "binwidth", False)
    ratio_mode = getattr(AnalysisConfig, "ratio", False)
    exclude_samples = getattr(AnalysisConfig, "exclude", [])

    print("Loading FHC CCnue...")
    fhc_ccnue = load_ccnue_fhc()

    print("Loading RHC CCnuebar...")
    rhc_ccnuebar = load_ccnuebar_rhc()

    print("Loading FHC CCnumu...")
    fhc_ccnumu = load_ccnumu_fhc()

    print("Loading RHC CCnumubar...")
    rhc_ccnumubar = load_ccnumubar_rhc()

    print("Loading FHC nue elastic...")
    fhc_elastic = load_nue_elastic_fhc()

    print("\n===== INPUT CHECK =====")
    print("fhc_ccnue mc integral   :", fhc_ccnue["mc"].Integral())
    print("fhc_ccnue data integral :", fhc_ccnue["data"].Integral())
    print("fhc_ccnue swap integral :", fhc_ccnue["swap"].Integral())
    print("rhc_ccnuebar mc integral   :", rhc_ccnuebar["mc"].Integral())
    print("rhc_ccnuebar data integral :", rhc_ccnuebar["data"].Integral())
    print("rhc_ccnuebar swap integral :", rhc_ccnuebar["swap"].Integral())
    print("fhc_ccnumu mc integral  :", fhc_ccnumu["mc"].Integral())
    print("fhc_ccnumu data integral:", fhc_ccnumu["data"].Integral())
    print("rhc_ccnumubar mc integral  :", rhc_ccnumubar["mc"].Integral())
    print("rhc_ccnumubar data integral:", rhc_ccnumubar["data"].Integral())
    print("fhc_elastic mc integral       :", fhc_elastic["mc"].Integral())
    print("fhc_elastic data integral     :", fhc_elastic["data"].Integral())
    print("fhc_elastic electron integral :", fhc_elastic["electron"].Integral())
    print("fhc_elastic muon integral     :", fhc_elastic["muon"].Integral())

    sample_histogram = StitchedHistogram("sample")

    # Add regular samples first, especially numu samples
    sample_histogram.AddHistograms("fhc_nue_selection", fhc_ccnue["mc"], fhc_ccnue["data"])
    sample_histogram.AddHistograms("rhc_nue_selection", rhc_ccnuebar["mc"], rhc_ccnuebar["data"])

    sample_histogram.AddHistograms("fhc_numu_selection", fhc_ccnumu["mc"], fhc_ccnumu["data"])
    sample_histogram.AddHistograms("rhc_numu_selection", rhc_ccnumubar["mc"], rhc_ccnumubar["data"])

    # Now swapped samples can safely look up fhc/rhc numu integrals
    sample_histogram.AddSwappedSample("fhc_nue_selection", fhc_ccnue["swap"])
    sample_histogram.AddSwappedSample("rhc_nue_selection", rhc_ccnuebar["swap"])

    # Now nue elastic scattering
    sample_histogram.AddScatteringFlavors("electron_fhc_elastic", fhc_elastic["electron"])
    sample_histogram.AddScatteringFlavors("muon_fhc_elastic", fhc_elastic["muon"])

    sample_histogram.AddHistograms("fhc_elastic", fhc_elastic["mc"], fhc_elastic["data"])



    # Add all templates before making ratios / copying
    sample_histogram.AddTemplates(
        "fhc_elastic",
        nue=fhc_elastic["template_electron"],
        numu=fhc_elastic["template_muon"],
        swap=fhc_elastic["template_muon"],
    )
    sample_histogram.AddTemplates(
        "fhc_numu_selection",
        numu=fhc_ccnumu["template_numu"],
    )
    sample_histogram.AddTemplates(
        "rhc_numu_selection",
        numu=rhc_ccnumubar["template_numu"],
    )
    sample_histogram.AddTemplates(
        "fhc_nue_selection",
        nue=fhc_ccnue["template_nue"],
        swap=fhc_ccnue["template_swap"],
    )
    sample_histogram.AddTemplates(
        "rhc_nue_selection",
        nue=rhc_ccnuebar["template_nue"],
        swap=rhc_ccnuebar["template_swap"],
    )

    sample_histogram.CleanErrorBands(errsToRemove)
    stitched = copy.deepcopy(sample_histogram)

    if ratio_mode:
        print("\nRatio mode is ON")
        print("  Making FHC flavor ratio")
        stitched.MakeRatio("fhc")
        print("  Making RHC flavor ratio")
        stitched.MakeRatio("rhc")
    else:
        print("\nRatio mode is OFF: using direct CCnue/CCnuebar/CCnumu/CCnumubar samples")

    # Apply broad keyword exclusions, e.g. --exclude nue
    stitched.ApplyExclusion(exclude_samples)

    # Force final stitched order.
    if ratio_mode:
        desired_order = [
            "fhc_elastic",
            "fhc_ratio",
            "rhc_ratio",
            "fhc_numu_selection",
            "rhc_numu_selection",
        ]
    else:
        desired_order = [
            "fhc_elastic",
            "fhc_nue_selection",
            "rhc_nue_selection",
            "fhc_numu_selection",
            "rhc_numu_selection",
        ]

    stitched.stitchKeys = [
        k for k in desired_order
        if k in stitched.stitchKeys
    ]

    print("\nSamples to be stitched:")
    for k in stitched.stitchKeys:
        print(" ", k)

    def check_bad_bins(h, label):
        print("\nChecking", label)
        print("Integral:", h.Integral())
        for i in range(0, h.GetNbinsX() + 2):
            v = h.GetBinContent(i)
            e = h.GetBinError(i)
            if not np.isfinite(v) or not np.isfinite(e):
                print("  BAD bin", i, "content =", v, "error =", e)

    check_bad_bins(rhc_ccnumubar["mc"], "rhc_ccnumubar mc")
    check_bad_bins(rhc_ccnumubar["data"], "rhc_ccnumubar data")

    print("About to stitch...")
    stitched.Stitch()
    print("Finished stitch.")

    print("\n===== STITCHED CHECK =====")
    print("stitched data integral:", stitched.data_hist.Integral())
    print("stitched mc integral  :", stitched.mc_hist.Integral())
    print("stitched swap integral:", stitched.swap_hist.Integral())
    print("stitched nbins        :", stitched.mc_hist.GetNbinsX())

    np.savetxt("mc_cv.csv", np.array(stitched.mc_hist)[1:-1], delimiter=",")
    np.savetxt("data_cv.csv", np.array(stitched.data_hist)[1:-1], delimiter=",")

    outroot = "{}/oscillations/rootfiles/NuE_stitched_hists.root".format(ccnueroot)
    stitched.Write(outroot)
    print("Wrote stitched file to", outroot)

    c = ROOT.TCanvas("c", "c", 900, 700)
    c.SetLogy()
    stitched.mc_hist.SetLineWidth(2)
    stitched.mc_hist.Draw("HIST")
    stitched.data_hist.Draw("E1 SAME")
    c.Print("stitched_test.png")

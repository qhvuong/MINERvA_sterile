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

errsToRemove = ["LowQ2Pi"]
ROOT.TH1.AddDirectory(False)


def clone_total(holder, name):
    h = holder.hists["Total"].Clone(name)
    h.SetDirectory(0)
    return h

def load_nue_elastic_fhc():
    type_path_map = {
        "data": "/exp/minerva/data/users/qvuong/elastic_nue/kin_dist_dataleFHC_NuE_allSystematics_newFlux_MAD.root",
        "mc":   "/exp/minerva/data/users/qvuong/elastic_nue/kin_dist_mcleFHC_NuE_allSystematics_newFlux_MAD.root",
    }

    data_file, mc_file, pot_scale, data_pot, mc_pot = Utilities.getFilesAndPOTScale(
        "NuEElastic", type_path_map, "MAD", True
    )

    standPOT = data_pot if data_pot is not None else mc_pot

    # Main 1D observable for the stitched elastic sample.
    # Change "Electron Energy" if your elastic sample uses a different PlotLibrary name.
    mc_holder = HistHolder("Electron Energy", mc_file, "Signal", True, mc_pot, standPOT)
    data_holder = HistHolder("Electron Energy", data_file, "Signal", False, data_pot, standPOT)

    # L/E templates.
    template_holder = HistHolder("Reco Energy vs L/E", mc_file, "Signal", True, mc_pot, standPOT)

    binwidthScale = getattr(AnalysisConfig, "binwidth", False)

    mc_holder.POTScale(binwidthScale)
    data_holder.POTScale(binwidthScale)
    template_holder.POTScale(binwidthScale)

    h_mc = clone_total(mc_holder, "fhc_elastic_mc")
    h_data = data_holder.GetHist().Clone("fhc_elastic_data")
    h_template = clone_total(template_holder, "fhc_elastic_template")

    # Flavor pieces for oscillation bookkeeping.
    # These category names may need to match your actual elastic file.
    h_electron = mc_file.Get("Electron_Scattering")
    h_muon     = mc_file.Get("Muon_Scattering")

    if not h_electron:
        h_electron = mc_file.Get("ENueElastic")
    if not h_muon:
        h_muon = mc_file.Get("ENumuElastic")

    if not h_electron or not h_muon:
        raise RuntimeError("Could not find electron/muon flavor components for FHC elastic sample")

    h_electron = h_electron.Clone("electron_fhc_elastic")
    h_muon = h_muon.Clone("muon_fhc_elastic")

    h_data.SetDirectory(0)
    h_template.SetDirectory(0)
    h_electron.SetDirectory(0)
    h_muon.SetDirectory(0)

    return {
        "mc": h_mc,
        "data": h_data,
        "electron": h_electron,
        "muon": h_muon,
        "template_electron": h_template,
        "template_muon": h_template.Clone("fhc_elastic_template_muon"),
    }

def load_ccnue_fhc():
    type_path_map = {
        "data": "/exp/minerva/data/users/qvuong/nu_e/kin_dist_dataleFHC_CCnue_allSystematics_newFlux_MAD.root",
        "mc":   "/exp/minerva/data/users/qvuong/nu_e/kin_dist_mcleFHC_CCnue_allSystematics_newFlux_MAD.root",
    }

    data_file, mc_file, pot_scale, data_pot, mc_pot = Utilities.getFilesAndPOTScale(
        "CCnue_allSystematics_newFlux", type_path_map, "MAD", True
    )
    standPOT = data_pot if data_pot is not None else mc_pot

    tuned_file = ROOT.TFile.Open(
        "/exp/minerva/data/users/qvuong/nu_e/bkgfit_leFHC_N4_tune_CCnue_allSystematics_newFlux_MAD.root"
    )
    h_mc = tuned_file.Get("EN4_predicted_Signal")
    h_data = tuned_file.Get("EN4_data_bkgSubbed")

    if not h_mc or not h_data:
        raise RuntimeError("Could not load CCnue tuned histograms from bkgfit file.")

    h_mc = h_mc.Clone("fhc_ccnue_mc")
    h_data = h_data.Clone("fhc_ccnue_data")
    h_mc.SetDirectory(0)
    h_data.SetDirectory(0)

    template_file = ROOT.TFile.Open(
        "/exp/minerva/data/users/qvuong/nu_e/kin_dist_mcleFHC_CCnue_allSystematics_newFlux_MAD.root"
    )
    template_holder = HistHolder("Reco Energy vs L/E", template_file, "Signal", True, mc_pot, standPOT)

    swap_type_path_map = {
        "mc": "/exp/minerva/data/users/qvuong/nu_e_swapped/kin_dist_mcleFHC_CCnueswap_allSystematics_newFlux_MAD.root"
    }
    _, swap_mc_file, _, _, swap_mc_pot = Utilities.getFilesAndPOTScale(
        "CCnueswap_allSystematics", swap_type_path_map, "MAD", True
    )

    swap_file = ROOT.TFile.Open(
        "/exp/minerva/data/users/qvuong/nu_e_swapped/kin_dist_mcleFHC_CCnueswap_allSystematics_newFlux_MAD.root"
    )
    swap_template_holder = HistHolder("Reco Energy vs L/E", swap_file, "Signal", True, swap_mc_pot, standPOT)
    swap_hist_holder = HistHolder("Biased Neutrino Energy", swap_file, "Signal", True, swap_mc_pot, standPOT)

    swap_hist_holder.POTScale(getattr(AnalysisConfig, "binwidth", False))
    h_swap = clone_total(swap_hist_holder, "fhc_ccnue_swap")
    h_swap.SetDirectory(0)

    h_template = clone_total(template_holder, "fhc_ccnue_template")
    h_swap_template = clone_total(swap_template_holder, "fhc_ccnue_swap_template")

    return {
        "mc": h_mc,
        "data": h_data,
        "swap": h_swap,
        "template_nue": h_template,
        "template_swap": h_swap_template,
    }


def load_ccnuebar_rhc():
    type_path_map = {
        "data": "/exp/minerva/data/users/qvuong/antinu_e/kin_dist_datale5_CCnuebar_allSystematics_newFlux_MAD.root",
        "mc":   "/exp/minerva/data/users/qvuong/antinu_e/kin_dist_mcle5_CCnuebar_allSystematics_newFlux_MAD.root",
    }

    data_file, mc_file, pot_scale, data_pot, mc_pot = Utilities.getFilesAndPOTScale(
        "CCnuebar_allSystematics", type_path_map, "MAD", True
    )
    standPOT = data_pot if data_pot is not None else mc_pot

    tuned_file = ROOT.TFile.Open(
        "/exp/minerva/data/users/qvuong/antinu_e/bkgfit_le5_N4_tune_CCnuebar_allSystematics_newFlux_MAD.root"
    )
    h_mc = tuned_file.Get("EN4_predicted_Signal")
    h_data = tuned_file.Get("EN4_data_bkgSubbed")

    if not h_mc or not h_data:
        raise RuntimeError("Could not load CCnuebar tuned histograms from bkgfit file.")

    h_mc = h_mc.Clone("rhc_ccnuebar_mc")
    h_data = h_data.Clone("rhc_ccnuebar_data")
    h_mc.SetDirectory(0)
    h_data.SetDirectory(0)

    template_file = ROOT.TFile.Open(
        "/exp/minerva/data/users/qvuong/antinu_e/kin_dist_mcle5_CCnuebar_allSystematics_newFlux_MAD.root"
    )
    template_holder = HistHolder("Reco Energy vs L/E", template_file, "Signal", True, mc_pot, standPOT)

    swap_type_path_map = {
        "mc": "/exp/minerva/data/users/qvuong/antinu_e_swapped/kin_dist_mcle5_CCnuebarswap_allSystematics_newFlux_MAD.root"
    }
    _, swap_mc_file, _, _, swap_mc_pot = Utilities.getFilesAndPOTScale(
        "CCnuebarswap_allSystematics", swap_type_path_map, "MAD", True
    )

    swap_file = ROOT.TFile.Open(
        "/exp/minerva/data/users/qvuong/antinu_e_swapped/kin_dist_mcle5_CCnuebarswap_allSystematics_newFlux_MAD.root"
    )
    swap_template_holder = HistHolder("Reco Energy vs L/E", swap_file, "Signal", True, swap_mc_pot, standPOT)
    swap_hist_holder = HistHolder("Biased Neutrino Energy", swap_file, "Signal", True, swap_mc_pot, standPOT)

    swap_hist_holder.POTScale(getattr(AnalysisConfig, "binwidth", False))
    h_swap = clone_total(swap_hist_holder, "rhc_ccnuebar_swap")
    h_swap.SetDirectory(0)

    h_template = clone_total(template_holder, "rhc_ccnuebar_template")
    h_swap_template = clone_total(swap_template_holder, "rhc_ccnuebar_swap_template")

    return {
        "mc": h_mc,
        "data": h_data,
        "swap": h_swap,
        "template_nue": h_template,
        "template_swap": h_swap_template,
    }


def load_ccnumu_fhc():
    type_path_map = {
        "data": "/exp/minerva/data/users/qvuong/nu_mu/kin_dist_dataleFHC_CCnumu_allSystematics_newFlux_MAD.root",
        "mc":   "/exp/minerva/data/users/qvuong/nu_mu/kin_dist_mcleFHC_CCnumu_allSystematics_newFlux_MAD.root",
    }

    data_file, mc_file, pot_scale, data_pot, mc_pot = Utilities.getFilesAndPOTScale(
        "CCnumu_allSystematics", type_path_map, "MAD", True
    )
    standPOT = data_pot if data_pot is not None else mc_pot

    mc_holder = HistHolder("Biased Neutrino Energy", mc_file, "Signal", True, mc_pot, standPOT)
    data_holder = HistHolder("Biased Neutrino Energy", data_file, "Signal", False, data_pot, standPOT)
    template_holder = HistHolder("Reco Energy vs L/E", mc_file, "Signal", True, mc_pot, standPOT)

    binwidthScale = getattr(AnalysisConfig, "binwidth", False)
    mc_holder.POTScale(binwidthScale)
    data_holder.POTScale(binwidthScale)
    template_holder.POTScale(binwidthScale)

    # For first-stage stitching validation, use Total directly.
    h_mc = clone_total(mc_holder, "fhc_ccnumu_mc")
    h_data = data_holder.GetHist().Clone("fhc_ccnumu_data")
    h_template = clone_total(template_holder, "fhc_ccnumu_template")

    h_data.SetDirectory(0)
    h_template.SetDirectory(0)

    return {
        "mc": h_mc,
        "data": h_data,
        "template_numu": h_template,
    }

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

def load_ccnumubar_rhc():
    type_path_map = {
        "data": "/exp/minerva/data/users/qvuong/antinu_mu/kin_dist_datale5_CCnumubar_allSystematics_newFlux_MAD.root",
        "mc":   "/exp/minerva/data/users/qvuong/antinu_mu/kin_dist_mcle5_CCnumubar_allSystematics_newFlux_MAD.root",
    }

    data_file, mc_file, pot_scale, data_pot, mc_pot = Utilities.getFilesAndPOTScale(
        "CCnumubar_allSystematics", type_path_map, "MAD", True
    )

    standPOT = data_pot if data_pot is not None else mc_pot

    if mc_pot is None or mc_pot == 0:
        raise RuntimeError("Bad mc_pot for RHC CCnumubar: {}".format(mc_pot))

    scale = standPOT / mc_pot

    print("\nPOT DEBUG: rhc_ccnumubar")
    print("  data_pot =", data_pot)
    print("  mc_pot   =", mc_pot)
    print("  standPOT =", standPOT)
    print("  scale    =", scale)

    # Use CCnumu-style exclusive categories.
    # Do NOT include EN4_CC or EN4_NC here unless you verify they are exclusive.
    ccnumu_1d_components = [
        "EN4_CCNuMuWrongSign",
        "EN4_CCNuMuQE",
        "EN4_CCNuMuDelta",
        "EN4_CCNuMuDIS",
        "EN4_CCNuMu2p2h",
        "EN4_CCNuMu",
        "EN4_NCDiff",
        "EN4_NuEElastic",
        "EN4_Other",
    ]

    ccnumu_2d_components = [
        "EReco_LE_CCNuMuWrongSign",
        "EReco_LE_CCNuMuQE",
        "EReco_LE_CCNuMuDelta",
        "EReco_LE_CCNuMuDIS",
        "EReco_LE_CCNuMu2p2h",
        "EReco_LE_CCNuMu",
        "EReco_LE_NCDiff",
        "EReco_LE_NuEElastic",
        "EReco_LE_Other",
    ]

    # Build MC total manually from CCnumu categories.
    h_mc = sum_hists_from_file(
        mc_file,
        ccnumu_1d_components,
        "rhc_ccnumubar_mc"
    )
    h_mc.Scale(scale)

    # Data likely only has inclusive EN4, so use it directly.
    h_data = data_file.Get("EN4")
    if not h_data:
        raise RuntimeError("Could not find EN4 in RHC CCnumubar data file.")

    h_data = h_data.Clone("rhc_ccnumubar_data")
    h_data.SetDirectory(0)

    # Build L/E template manually too, avoiding CCnue HistHolder categories.
    h_template = sum_hists_from_file(
        mc_file,
        ccnumu_2d_components,
        "rhc_ccnumubar_template"
    )
    h_template.Scale(scale)

    check_bad_bins(h_mc, "rhc_ccnumubar mc FINAL")
    check_bad_bins(h_data, "rhc_ccnumubar data FINAL")

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

    sample_histogram = StitchedHistogram("sample")

    # Add regular samples first, especially numu samples
    sample_histogram.AddHistograms("fhc_nue_selection", fhc_ccnue["mc"], fhc_ccnue["data"])
    sample_histogram.AddHistograms("rhc_nue_selection", rhc_ccnuebar["mc"], rhc_ccnuebar["data"])

    sample_histogram.AddHistograms("fhc_numu_selection", fhc_ccnumu["mc"], fhc_ccnumu["data"])
    sample_histogram.AddHistograms("rhc_numu_selection", rhc_ccnumubar["mc"], rhc_ccnumubar["data"])

    # Now swapped samples can safely look up fhc/rhc numu integrals
    sample_histogram.AddSwappedSample("fhc_nue_selection", fhc_ccnue["swap"])
    sample_histogram.AddSwappedSample("rhc_nue_selection", rhc_ccnuebar["swap"])

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
    sample_histogram.AddTemplates(
        "fhc_numu_selection",
        numu=fhc_ccnumu["template_numu"],
    )
    sample_histogram.AddTemplates(
        "rhc_numu_selection",
        numu=rhc_ccnumubar["template_numu"],
    )

    sample_histogram.CleanErrorBands(errsToRemove)
    stitched = copy.deepcopy(sample_histogram)

    if ratio_mode:
        if "fhc" not in exclude_samples:
            stitched.MakeRatio("fhc")

    stitched.ApplyExclusion(exclude_samples)

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

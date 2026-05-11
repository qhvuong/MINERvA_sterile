import os
import copy
import ROOT
import PlotUtils
import numpy as np

from Tools.Histogram import *
from Tools.Helper import *

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


def load_ccnue_fhc():
    type_path_map = {
        "data": "/exp/minerva/data/users/qvuong/nu_e/kin_dist_dataleFHC_CCnue_allSystematics_newFlux_MAD.root",
        "mc":   "/exp/minerva/data/users/qvuong/nu_e/kin_dist_mcleFHC_CCnue_allSystematics_newFlux_MAD.root",
    }

    data_file, mc_file, pot_scale, data_pot, mc_pot = Utilities.getFilesAndPOTScale(
        "CCnue_allSystematics_testCVs", type_path_map, "MAD", True
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

def load_ccnumubar_rhc():
    type_path_map = {
        "data": "/exp/minerva/data/users/qvuong/antinu_mu/kin_dist_datale5_CCnumubar_allSystematics_newFlux_MAD.root",
        "mc":   "/exp/minerva/data/users/qvuong/antinu_mu/kin_dist_mcle5_CCnumubar_allSystematics_newFlux_MAD.root",
    }

    data_file, mc_file, pot_scale, data_pot, mc_pot = Utilities.getFilesAndPOTScale(
        "CCnumubar_allSystematics", type_path_map, "MAD", True
    )
    standPOT = data_pot if data_pot is not None else mc_pot

    mc_holder = HistHolder("Biased Neutrino Energy", mc_file, "Signal", True, mc_pot, standPOT)
    data_holder = HistHolder("Biased Neutrino Energy", data_file, "Signal", False, data_pot, standPOT)
    template_holder = HistHolder("Reco Energy vs L/E", mc_file, "Signal", True, mc_pot, standPOT)

    binwidthScale = getattr(AnalysisConfig, "binwidth", False)
    mc_holder.POTScale(binwidthScale)
    data_holder.POTScale(binwidthScale)
    template_holder.POTScale(binwidthScale)

    h_mc = clone_total(mc_holder, "rhc_ccnumubar_mc")
    h_data = data_holder.GetHist().Clone("rhc_ccnumubar_data")
    h_template = clone_total(template_holder, "rhc_ccnumubar_template")

    h_data.SetDirectory(0)
    h_template.SetDirectory(0)

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

    sample_histogram.AddSwappedSample("fhc_nue_selection", fhc_ccnue["swap"])
    sample_histogram.AddSwappedSample("rhc_nue_selection", rhc_ccnuebar["swap"])

    sample_histogram.AddHistograms("fhc_nue_selection", fhc_ccnue["mc"], fhc_ccnue["data"])
    sample_histogram.AddHistograms("rhc_nue_selection", rhc_ccnuebar["mc"], rhc_ccnuebar["data"])

    sample_histogram.AddHistograms("fhc_numu_selection", fhc_ccnumu["mc"], fhc_ccnumu["data"])
    sample_histogram.AddHistograms("rhc_numu_selection", rhc_ccnumubar["mc"], rhc_ccnumubar["data"])

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

    # c2 = ROOT.TCanvas("c2", "c2", 900, 700)
    # stitched.swap_hist.SetLineWidth(2)
    # stitched.swap_hist.Draw("HIST")
    # c2.Print("fhc_ccnue_swap_stitched_test.png")

    # print(f"Wrote {outroot}, fhc_ccnue_ccnumu_stitched_test.png, and fhc_ccnue_swap_stitched_test.png")
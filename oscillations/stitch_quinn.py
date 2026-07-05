import os
import copy
import ROOT
import PlotUtils
import numpy as np
import shutil
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

PAPER_ELASTIC_FILE = (
    "/exp/minerva/data/users/qvuong/elastic_nue/"
    "paper_nue_elastic_decomposed_with2D.root"
)

LELIKE_4BIN_ELASTIC_TEMPLATE_FILE = (
    "/exp/minerva/data/users/qvuong/elastic_nue/"
    "kin_dist_mcmeFHC_p6_scattering_LElike_MAD.root"
)

# POT normalization
# Published MINERvA nu+e paper exposure
PAPER_NUE_ELASTIC_POT = 3.43e20

# FHC exposure used by your current FHC CCnue/CCnumu samples
FHC_ANALYSIS_POT = 3.323050142731963e20

def get_pot_from_meta(fin):
    meta = fin.Get("Meta")
    if not meta:
        raise RuntimeError("Input file has no Meta tree.")

    branch_names = [b.GetName() for b in meta.GetListOfBranches()]

    candidates = [
        "POT_Used",
        "POTUsed",
        "used_POT",
        "UsedPOT",
        "POT",
        "pot",
        "POT_Total",
        "total_pot",
    ]

    for name in candidates:
        if name in branch_names:
            total = 0.0
            for i in range(meta.GetEntries()):
                meta.GetEntry(i)
                total += float(getattr(meta, name))

            print("Read POT from Meta branch {} = {}".format(name, total))
            return total

    raise RuntimeError(
        "Could not find POT branch in Meta. Available branches = {}".format(branch_names)
    )

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

def scale_tmatrix(mat, scale):
    """
    Scale a covariance matrix by scale^2.
    Counts scale as scale, covariance scales as scale^2.
    """
    if mat is None:
        return None

    for i in range(mat.GetNrows()):
        for j in range(mat.GetNcols()):
            mat[i][j] = mat[i][j] * scale * scale

    return mat


def scale_hist_list(hists, scale):
    """
    Scale ROOT histograms by a common POT factor.
    This scales contents and bin errors.
    """
    for h in hists:
        if h is not None:
            h.Scale(scale)


def exclude_has(exclude, token):
    """
    Check comma-separated AnalysisConfig.exclude safely.
    Example: exclude='nue,elastic'
    """
    if exclude is None:
        return False
    if isinstance(exclude, str):
        parts = [x.strip() for x in exclude.split(",") if x.strip()]
    else:
        parts = list(exclude)
    return token in parts

# This is my LE FHC nue elastic sample
# This is my LE FHC nue elastic sample, using production 1D normalization
# and FHC-scaled LElike 2D templates from PAPER_ELASTIC_FILE.
def load_nue_elastic_fhc():
    type_path_map = {
        "data": "/exp/minerva/data/users/qvuong/elastic_nue/kin_dist_dataleFHC_NuE_allSystematics_fullStatsFluxes_MAD.root",
        "mc":   "/exp/minerva/data/users/qvuong/elastic_nue/kin_dist_mcleFHC_NuE_allSystematics_fullStatsFluxes_MAD.root",
    }

    data_file, mc_file, pot_scale, data_pot, mc_pot = Utilities.getFilesAndPOTScale(
        "NuE_allSystematics_fullStatsFluxes", type_path_map, "MAD", True
    )

    standPOT = data_pot if data_pot is not None else mc_pot
    print_pot_scale_check("FHC NuE elastic", data_pot, mc_pot, standPOT)

    tuned_file = ROOT.TFile.Open(
        "/exp/minerva/data/users/qvuong/elastic_nue/bkgfit_leFHC_nueElastic_matrix_NuE_allSystematics_fullStatsFluxes_MAD.root"
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
    # 1D flavor components from production MC
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

    print("\nFHC elastic production flavor POT info:")
    print("  data_pot =", data_pot)
    print("  mc_pot   =", mc_pot)
    print("  standPOT =", standPOT)
    print("  scale    =", scale)

    h_electron_raw = h_electron.Clone("electron_fhc_elastic_raw")
    h_muon_raw = h_muon.Clone("muon_fhc_elastic_raw")
    h_electron_raw.SetDirectory(0)
    h_muon_raw.SetDirectory(0)

    h_electron.Scale(scale)
    h_muon.Scale(scale)

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
    else:
        raise RuntimeError("FHC elastic electron+muon flavor sum is zero.")

    # -------------------------------------------------
    # 2D flavor templates from already-FHC-scaled LElike templates
    # stored in PAPER_ELASTIC_FILE.
    # -------------------------------------------------
    template_file_2d = ROOT.TFile.Open(LELIKE_4BIN_ELASTIC_TEMPLATE_FILE)
    if not template_file_2d or template_file_2d.IsZombie():
        raise RuntimeError(
            "Could not open LElike 4-bin elastic template file: {}".format(
                LELIKE_4BIN_ELASTIC_TEMPLATE_FILE
            )
        )

    print("\nLoading FHC elastic 2D templates from LElike 4-bin file")
    print("  file =", LELIKE_4BIN_ELASTIC_TEMPLATE_FILE)

    h2_nue = get_hist_checked(
        template_file_2d,
        ["drawnL_ElepReco_LE_NuEElasticE"],
        "fhc_elastic_template_nue",
    )
    h2_nuebar = get_hist_checked(
        template_file_2d,
        ["drawnL_ElepReco_LE_NuEElasticEBar"],
        "fhc_elastic_template_nuebar",
    )
    h2_numu = get_hist_checked(
        template_file_2d,
        ["drawnL_ElepReco_LE_NuEElasticMu"],
        "fhc_elastic_template_numu",
    )
    h2_numubar = get_hist_checked(
        template_file_2d,
        ["drawnL_ElepReco_LE_NuEElasticMuBar"],
        "fhc_elastic_template_numubar",
    )

    template_pot = get_pot_from_meta(template_file_2d)

    if template_pot is None or template_pot == 0:
        raise RuntimeError(
            "Bad POT for LElike 4-bin template file: {}".format(template_pot)
        )

    template_pot_scale = standPOT / template_pot

    print("\nLElike 4-bin template POT scaling:")
    print("  template POT =", template_pot)
    print("  target POT   =", standPOT)
    print("  scale        =", template_pot_scale)

    h2_nue_raw = h2_nue.Clone("fhc_elastic_template_nue_raw")
    h2_nuebar_raw = h2_nuebar.Clone("fhc_elastic_template_nuebar_raw")
    h2_numu_raw = h2_numu.Clone("fhc_elastic_template_numu_raw")
    h2_numubar_raw = h2_numubar.Clone("fhc_elastic_template_numubar_raw")

    for h in [h2_nue_raw, h2_nuebar_raw, h2_numu_raw, h2_numubar_raw]:
        h.SetDirectory(0)

    for h in [h2_nue, h2_nuebar, h2_numu, h2_numubar]:
        h.Scale(template_pot_scale)

    print_integral_change("FHC elastic 2D nue template", h2_nue_raw, h2_nue)
    print_integral_change("FHC elastic 2D nuebar template", h2_nuebar_raw, h2_nuebar)
    print_integral_change("FHC elastic 2D numu template", h2_numu_raw, h2_numu)
    print_integral_change("FHC elastic 2D numubar template", h2_numubar_raw, h2_numubar)

    h_template_electron = h2_nue.Clone("fhc_elastic_template_electron")
    h_template_electron.SetDirectory(0)
    h_template_electron.Add(h2_nuebar)

    h_template_muon = h2_numu.Clone("fhc_elastic_template_muon")
    h_template_muon.SetDirectory(0)
    h_template_muon.Add(h2_numubar)

    template_file_2d.Close()
    tuned_file.Close()

    compare_1d_2d_template(
        "FHC elastic electron production",
        h_electron,
        h_template_electron
    )

    compare_1d_2d_template(
        "FHC elastic muon production",
        h_muon,
        h_template_muon
    )

    print("\nFHC NuE elastic production check:")
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
        "template_nue": h2_nue,
        "template_nuebar": h2_nuebar,
        "template_numu": h2_numu,
        "template_numubar": h2_numubar,
        "cov": None,
    }

# This is Jaewon's result
def load_nue_elastic_fhc_jaewon():
    """
    Load FHC nu+e elastic sample from the paper-normalized decomposition file.

    The returned objects are:
      mc/data:
        total paper-normalized 1D spectra

      electron/muon:
        nue+nuebar and numu+numubar components

      template_electron/template_muon:
        corresponding reco Ee vs true L/E templates

      cov:
        published 6x6 paper covariance matrix
    """

    paper_file = ROOT.TFile.Open(PAPER_ELASTIC_FILE)
    if not paper_file or paper_file.IsZombie():
        raise RuntimeError("Could not open paper decomposed elastic file: {}".format(PAPER_ELASTIC_FILE))

    print("\n===== Loading FHC nu+e elastic from paper decomposition =====")
    print("  file =", PAPER_ELASTIC_FILE)

    h_data = get_hist_checked(
        paper_file,
        ["paper_total_nue_elastic"],
        "fhc_elastic_data_paper",
    )

    h_mc = get_hist_checked(
        paper_file,
        ["paper_decomp_sum"],
        "fhc_elastic_mc_paper_decomp_sum",
    )

    # # Use the paper total diagonal errors on the nominal/reference total.
    # # The full bin-to-bin covariance is stored separately in cov.
    # for i in range(0, h_mc.GetNbinsX() + 2):
    #     h_mc.SetBinError(i, h_data.GetBinError(i))

    h_nue = get_hist_checked(
        paper_file,
        ["paper_decomp_nue"],
        "fhc_elastic_component_nue",
    )
    h_nuebar = get_hist_checked(
        paper_file,
        ["paper_decomp_nuebar"],
        "fhc_elastic_component_nuebar",
    )
    h_numu = get_hist_checked(
        paper_file,
        ["paper_decomp_numu"],
        "fhc_elastic_component_numu",
    )
    h_numubar = get_hist_checked(
        paper_file,
        ["paper_decomp_numubar"],
        "fhc_elastic_component_numubar",
    )

    h_electron = h_nue.Clone("electron_fhc_elastic")
    h_electron.SetDirectory(0)
    h_electron.Add(h_nuebar)

    h_muon = h_numu.Clone("muon_fhc_elastic")
    h_muon.SetDirectory(0)
    h_muon.Add(h_numubar)

    h2_nue = get_hist_checked(
        paper_file,
        ["source_fhc_2d_nue"],
        "fhc_elastic_template_nue",
    )
    h2_nuebar = get_hist_checked(
        paper_file,
        ["source_fhc_2d_nuebar"],
        "fhc_elastic_template_nuebar",
    )
    h2_numu = get_hist_checked(
        paper_file,
        ["source_fhc_2d_numu"],
        "fhc_elastic_template_numu",
    )
    h2_numubar = get_hist_checked(
        paper_file,
        ["source_fhc_2d_numubar"],
        "fhc_elastic_template_numubar",
    )

    h_template_electron = h2_nue.Clone("fhc_elastic_template_electron")
    h_template_electron.SetDirectory(0)
    h_template_electron.Add(h2_nuebar)

    h_template_muon = h2_numu.Clone("fhc_elastic_template_muon")
    h_template_muon.SetDirectory(0)
    h_template_muon.Add(h2_numubar)

    cov = paper_file.Get("paper_covariance_matrix")
    if cov:
        cov = cov.Clone("fhc_elastic_paper_covariance_matrix")
        print("\nLoaded paper covariance matrix.")
        print("  Nrows =", cov.GetNrows())
        print("  Ncols =", cov.GetNcols())
    else:
        print("\nWARNING: paper_covariance_matrix not found in paper file.")
        cov = None

    # print_hist_bins(h_data, "FHC elastic paper total data")
    # print_hist_bins(h_mc, "FHC elastic paper decomposed MC sum")
    # print_hist_bins(h_electron, "FHC elastic electron = nue+nuebar")
    # print_hist_bins(h_muon, "FHC elastic muon = numu+numubar")

    # # -------------------------------------------------
    # # POT-scale paper nu+e result to the FHC exposure
    # # used by the other FHC samples in the stitched fit.
    # # -------------------------------------------------
    # pot_scale = PAPER_TO_FHC_POT_SCALE

    # print("\nScaling paper nu+e elastic sample to FHC analysis POT")
    # print("  paper POT        =", PAPER_NUE_ELASTIC_POT)
    # print("  FHC analysis POT =", FHC_ANALYSIS_POT)
    # print("  scale            =", pot_scale)

    # scale_hist_list(
    #     [
    #         h_data,
    #         h_mc,
    #         h_nue,
    #         h_nuebar,
    #         h_numu,
    #         h_numubar,
    #         h_electron,
    #         h_muon,
    #         h_template_electron,
    #         h_template_muon,
    #         h2_nue,
    #         h2_nuebar,
    #         h2_numu,
    #         h2_numubar,
    #     ],
    #     pot_scale,
    # )

    # # The covariance is in event-count units^2, so scale by POT^2.
    # cov = scale_tmatrix(cov, pot_scale)

    print("\nPaper decomposition file is already scaled to FHC analysis POT.")
    # print("  target FHC POT =", FHC_ANALYSIS_POT)
    # print("  No additional POT scaling applied in stitch_quinn.py.")

    # Make sure the nominal/reference total uses the same diagonal errors as data.
    # This is just for display/default histogram errors; full covariance is separate.
    for i in range(0, h_mc.GetNbinsX() + 2):
        h_mc.SetBinError(i, h_data.GetBinError(i))

    # Closure check must be built AFTER scaling.
    h_flavor_sum = h_electron.Clone("fhc_elastic_flavor_sum_check")
    h_flavor_sum.SetDirectory(0)
    h_flavor_sum.Add(h_muon)

    print_hist_bins(h_data, "FHC elastic paper total data, already FHC POT scaled")
    print_hist_bins(h_mc, "FHC elastic paper decomposed MC sum, already FHC POT scaled")
    print_hist_bins(h_electron, "FHC elastic electron = nue+nuebar, already FHC POT scaled")
    print_hist_bins(h_muon, "FHC elastic muon = numu+numubar, already FHC POT scaled")

    print("\nFHC elastic closure check:")
    print("  paper data integral    =", h_data.Integral())
    print("  paper MC integral      =", h_mc.Integral())
    print("  electron integral      =", h_electron.Integral())
    print("  muon integral          =", h_muon.Integral())
    print("  electron+muon integral =", h_flavor_sum.Integral())
    print("  e+mu - paper data      =", h_flavor_sum.Integral() - h_data.Integral())

    compare_1d_2d_template(
        "FHC elastic electron paper-decomp",
        h_electron,
        h_template_electron,
    )
    compare_1d_2d_template(
        "FHC elastic muon paper-decomp",
        h_muon,
        h_template_muon,
    )

    paper_file.Close()

    return {
        "mc": h_mc,
        "data": h_data,
        "electron": h_electron,
        "muon": h_muon,
        "template_electron": h_template_electron,
        "template_muon": h_template_muon,
        "nue": h_nue,
        "nuebar": h_nuebar,
        "numu": h_numu,
        "numubar": h_numubar,
        "template_nue": h2_nue,
        "template_nuebar": h2_nuebar,
        "template_numu": h2_numu,
        "template_numubar": h2_numubar,
        "cov": cov,
    }


def load_ccnue_fhc():
    type_path_map = {
        "data": "/exp/minerva/data/users/qvuong/nu_e/kin_dist_dataleFHC_CCnue_allSystematics_fullStatsFluxes_MAD.root",
        "mc":   "/exp/minerva/data/users/qvuong/nu_e/kin_dist_mcleFHC_CCnue_allSystematics_fullStatsFluxes_MAD.root",
    }

    data_file, mc_file, pot_scale, data_pot, mc_pot = Utilities.getFilesAndPOTScale(
        "CCnue_allSystematics_fullStatsFluxes", type_path_map, "MAD", True
    )

    standPOT = data_pot if data_pot is not None else mc_pot
    binwidthScale = getattr(AnalysisConfig, "binwidth", False)

    print_pot_scale_check("FHC CCnue", data_pot, mc_pot, standPOT)

    tuned_file = ROOT.TFile.Open(
        "/exp/minerva/data/users/qvuong/nu_e/bkgfit_leFHC_N4_tune_CCnue_allSystematics_fullStatsFluxes_MAD.root"
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
        "/exp/minerva/data/users/qvuong/nu_e/kin_dist_mcleFHC_CCnue_allSystematics_fullStatsFluxes_MAD.root"
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
        "mc": "/exp/minerva/data/users/qvuong/nu_e_swapped/kin_dist_mcleFHC_CCnueswap_allSystematics_fullStatsFluxes_MAD.root"
    }

    _, swap_mc_file, _, _, swap_mc_pot = Utilities.getFilesAndPOTScale(
        "CCnueswap_allSystematics_fullStatsFluxes", swap_type_path_map, "MAD", True
    )

    print_pot_scale_check("FHC CCnue swap", None, swap_mc_pot, standPOT)

    swap_file = ROOT.TFile.Open(
        "/exp/minerva/data/users/qvuong/nu_e_swapped/kin_dist_mcleFHC_CCnueswap_allSystematics_fullStatsFluxes_MAD.root"
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
        "data": "/exp/minerva/data/users/qvuong/antinu_e/kin_dist_datale5_CCnuebar_allSystematics_fullStatsFluxes_MAD.root",
        "mc":   "/exp/minerva/data/users/qvuong/antinu_e/kin_dist_mcle5_CCnuebar_allSystematics_fullStatsFluxes_MAD.root",
    }

    data_file, mc_file, pot_scale, data_pot, mc_pot = Utilities.getFilesAndPOTScale(
        "CCnuebar_allSystematics_fullStatsFluxes", type_path_map, "MAD", True
    )

    standPOT = data_pot if data_pot is not None else mc_pot
    binwidthScale = getattr(AnalysisConfig, "binwidth", False)

    print_pot_scale_check("RHC CCnuebar", data_pot, mc_pot, standPOT)

    tuned_file = ROOT.TFile.Open(
        "/exp/minerva/data/users/qvuong/antinu_e/bkgfit_le5_N4_tune_CCnuebar_allSystematics_fullStatsFluxes_MAD.root"
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
        "/exp/minerva/data/users/qvuong/antinu_e/kin_dist_mcle5_CCnuebar_allSystematics_fullStatsFluxes_MAD.root"
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
        "mc": "/exp/minerva/data/users/qvuong/antinu_e_swapped/kin_dist_mcle5_CCnuebarswap_allSystematics_fullStatsFluxes_MAD.root"
    }

    _, swap_mc_file, _, _, swap_mc_pot = Utilities.getFilesAndPOTScale(
        "CCnuebarswap_allSystematics_fullStatsFluxes", swap_type_path_map, "MAD", True
    )

    print_pot_scale_check("RHC CCnuebar swap", None, swap_mc_pot, standPOT)

    swap_file = ROOT.TFile.Open(
        "/exp/minerva/data/users/qvuong/antinu_e_swapped/kin_dist_mcle5_CCnuebarswap_allSystematics_fullStatsFluxes_MAD.root"
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
        "data": "/exp/minerva/data/users/qvuong/nu_mu/kin_dist_dataleFHC_CCnumu_allSystematics_fullStatsFluxes_MAD.root",
        "mc":   "/exp/minerva/data/users/qvuong/nu_mu/kin_dist_mcleFHC_CCnumu_allSystematics_fullStatsFluxes_MAD.root",
    }

    data_file, mc_file, pot_scale, data_pot, mc_pot = Utilities.getFilesAndPOTScale(
        "CCnumu_allSystematics_fullStatsFluxes", type_path_map, "MAD", True
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
        "data": "/exp/minerva/data/users/qvuong/antinu_mu/kin_dist_datale5_CCnumubar_allSystematics_fullStatsFluxes_MAD.root",
        "mc":   "/exp/minerva/data/users/qvuong/antinu_mu/kin_dist_mcle5_CCnumubar_allSystematics_fullStatsFluxes_MAD.root",
    }

    data_file, mc_file, pot_scale, data_pot, mc_pot = Utilities.getFilesAndPOTScale(
        "CCnumubar_allSystematics_fullStatsFluxes", type_path_map, "MAD", True
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

    hist_config_tag = getattr(AnalysisConfig, "hist_config_tag", "default")

    if hist_config_tag in [None, "", "none"]:
        hist_config_tag = "default"

    hist_config_out = "HIST_CONFIG_{}.json".format(hist_config_tag)
    OUTROOT = "{}/oscillations/rootfiles/NuE_stitched_hists_{}.root".format(
        ccnueroot,
        hist_config_tag,
    )

    print("hist_config_tag =", hist_config_tag)
    print("hist_config_out =", hist_config_out)
    print("OUTROOT         =", OUTROOT)

    print("Loading FHC CCnue...")
    fhc_ccnue = load_ccnue_fhc()

    print("Loading RHC CCnuebar...")
    rhc_ccnuebar = load_ccnuebar_rhc()

    print("Loading FHC CCnumu...")
    fhc_ccnumu = load_ccnumu_fhc()

    print("Loading RHC CCnumubar...")
    rhc_ccnumubar = load_ccnumubar_rhc()

    elastic_source = getattr(AnalysisConfig, "elastic_source", "production")

    print("Loading FHC nue elastic...")
    print("elastic_source =", elastic_source)

    if elastic_source in ["paper", "jaewon"]:
        fhc_elastic = load_nue_elastic_fhc_jaewon()
    elif elastic_source in ["production", "prod", "mine"]:
        fhc_elastic = load_nue_elastic_fhc()
    else:
        raise RuntimeError("Unknown elastic_source: {}".format(elastic_source))

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

    elastic_excluded = exclude_has(exclude_samples, "elastic")

    if (not elastic_excluded) and fhc_elastic.get("cov") is not None:
        stitched.external_covariances["fhc_elastic"] = fhc_elastic["cov"].Clone(
            "fhc_elastic_paper_covariance_matrix"
        )
        print("Registered external covariance for fhc_elastic before stitching.")
    else:
        print("Not registering fhc_elastic external covariance.")
        print("  elastic excluded =", elastic_excluded)
        print("  cov is None      =", fhc_elastic.get("cov") is None)

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

    if os.path.exists("HIST_CONFIG.json"):
        shutil.copyfile("HIST_CONFIG.json", hist_config_out)
        print("Copied HIST_CONFIG.json ->", hist_config_out)
    else:
        print("WARNING: HIST_CONFIG.json was not produced by Stitch().")

    print("\n===== STITCHED CHECK =====")
    print("stitched data integral:", stitched.data_hist.Integral())
    print("stitched mc integral  :", stitched.mc_hist.Integral())
    print("stitched swap integral:", stitched.swap_hist.Integral())
    print("stitched nbins        :", stitched.mc_hist.GetNbinsX())

    # np.savetxt("mc_cv.csv", np.array(stitched.mc_hist)[1:-1], delimiter=",")
    # np.savetxt("data_cv.csv", np.array(stitched.data_hist)[1:-1], delimiter=",")

    csv_dir = "{}/oscillations/csvs/{}".format(ccnueroot, hist_config_tag)

    if not os.path.exists(csv_dir):
        os.makedirs(csv_dir)

    mc_cv = np.array(stitched.mc_hist)[1:-1]
    data_cv = np.array(stitched.data_hist)[1:-1]

    cov_full = stitched.GetCovarianceMatrix(False)
    cov_sans = stitched.GetCovarianceMatrix(True)
    cov_flux = cov_full - cov_sans

    np.savetxt(
        "{}/mc_cv.csv".format(csv_dir),
        mc_cv,
        delimiter=",",
        header="mc_cv",
        comments="",
    )

    np.savetxt(
        "{}/data_cv.csv".format(csv_dir),
        data_cv,
        delimiter=",",
        header="data_cv",
        comments="",
    )

    np.savetxt(
        "{}/cov_full.csv".format(csv_dir),
        cov_full,
        delimiter=",",
    )

    np.savetxt(
        "{}/cov_sansFlux.csv".format(csv_dir),
        cov_sans,
        delimiter=",",
    )

    np.savetxt(
        "{}/cov_flux.csv".format(csv_dir),
        cov_flux,
        delimiter=",",
    )

    np.savetxt(
        "{}/inv_cov_full.csv".format(csv_dir),
        stitched.GetInverseCovarianceMatrix(False),
        delimiter=",",
    )

    np.savetxt(
        "{}/inv_cov_sansFlux.csv".format(csv_dir),
        stitched.GetInverseCovarianceMatrix(True),
        delimiter=",",
    )

    np.savetxt(
        "{}/A_flux_universes_minus_cv.csv".format(csv_dir),
        stitched.GetAMatrix(),
        delimiter=",",
    )

    print("Wrote stitched CSV outputs to", csv_dir)

    stitched.Write(OUTROOT)
    print("Wrote stitched file to", OUTROOT)

    # c = ROOT.TCanvas("c", "c", 900, 700)
    # c.SetLogy()
    # stitched.mc_hist.SetLineWidth(2)
    # stitched.mc_hist.Draw("HIST")
    # stitched.data_hist.Draw("E1 SAME")
    # c.Print("stitched_test.png")

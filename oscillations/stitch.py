import os
import sys
import argparse
import copy
import shutil
import json
import numpy as np
from array import array

SIGNAL_DEFINITION = [
    "CCNuEQE",
    "CCNuEDelta",
    "CCNuEDIS",
    "CCNuE",
    "CCNuE2p2h",
    "CCNuEWrongSign",
]

# POT normalization
# Published MINERvA nu+e paper exposure
PAPER_NUE_ELASTIC_POT = 3.43e20

# Full FHC data exposure:
# Minerva1 + 7 + 9 + 13A + 13B + 13C + 13D + 13E
FHC_ANALYSIS_POT = 3.331982991676675e20

# Parse stitch-only args before ANY analysis/config imports see sys.argv.
_stitch_parser = argparse.ArgumentParser(add_help=False)
_stitch_parser.add_argument(
    "--stitch-input-config",
    "--stitch_input_config",
    dest="stitch_input_config",
    default="stitch_input_files.json",
)
_stitch_parser.add_argument(
    "--stitch-input-set",
    "--stitch_input_set",
    dest="stitch_input_set",
    default="p6",
    choices=["p6", "p8", "p8_onlyPPFX", "p8_onlyBeamFocus"],
)

_stitch_args, _remaining_argv = _stitch_parser.parse_known_args()
sys.argv = [sys.argv[0]] + _remaining_argv

import ROOT
import PlotUtils

from tools.StitchedHistogram import *
from tools.Helper import *

from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS
from config.AnalysisConfig import AnalysisConfig
from tools import Utilities
from tools.PlotLibrary import HistHolder

AnalysisConfig.stitch_input_config = _stitch_args.stitch_input_config
AnalysisConfig.stitch_input_set = _stitch_args.stitch_input_set

ccnueroot = os.environ.get("CCNUEROOT")

def load_stitch_input_config():
    input_config_path = getattr(
        AnalysisConfig,
        "stitch_input_config",
        "stitch_input_files.json",
    )

    input_set = getattr(
        AnalysisConfig,
        "stitch_input_set",
        "p6",
    )

    with open(input_config_path, "r") as f:
        all_inputs = json.load(f)

    if input_set not in all_inputs:
        raise RuntimeError(
            "stitch_input_set={} not found in {}. Available sets = {}".format(
                input_set,
                input_config_path,
                sorted(all_inputs.keys()),
            )
        )

    print("Using stitch input config:", input_config_path)
    print("Using stitch input set   :", input_set)

    return all_inputs[input_set]


STITCH_INPUTS = load_stitch_input_config()


def stitch_path(sample, key):
    try:
        path = STITCH_INPUTS[sample][key]
    except KeyError:
        raise RuntimeError(
            "Missing stitch input path for sample='{}', key='{}'".format(
                sample, key
            )
        )

    print("[INPUT] {}.{} = {}".format(sample, key, path))
    return path


MNVPLOTTER = PlotUtils.MnvPlotter()
MNVPLOTTER.error_summary_group_map.clear()
for k, v in CONSOLIDATED_ERROR_GROUPS.items():
    vec = ROOT.vector("std::string")()
    for vs in v:
        vec.push_back(vs)
    MNVPLOTTER.error_summary_group_map[k] = vec

errsToRemove = ["LowQ2Pi","elETracker"]
ROOT.TH1.AddDirectory(False)

JAEWON_INPUT_FILE = "jaewon_nue_elastic_data.json"

LELIKE_6BIN_ELASTIC_TEMPLATE_FILE = (
    "/exp/minerva/data/users/qvuong/elastic_nue/kin_dist_mcmeFHC_p8_scattering_LElike_6bins_MAD.root"
)

PREDICTED_ELASTIC_FILE = (
    "/exp/minerva/data/users/qvuong/nueel_prediction_studies/nue_elastic_prediction_higherOrderXS_mnv_FHC.root"
)

JAEWON_ELEP_EDGES = [0.8, 2.0, 3.0, 5.0, 7.0, 9.0, 100.0]


def check_1d_edges(h, expected_edges, label, tol=1e-8):
    if h.GetNbinsX() != len(expected_edges) - 1:
        raise RuntimeError(
            "{} has {} bins, expected {}".format(
                label,
                h.GetNbinsX(),
                len(expected_edges) - 1,
            )
        )

    actual = [
        h.GetXaxis().GetBinLowEdge(i)
        for i in range(1, h.GetNbinsX() + 1)
    ]
    actual.append(h.GetXaxis().GetBinUpEdge(h.GetNbinsX()))

    for i, (a, e) in enumerate(zip(actual, expected_edges)):
        if abs(a - e) > tol:
            raise RuntimeError(
                "{} edge {} mismatch: actual={} expected={}".format(
                    label,
                    i,
                    a,
                    e,
                )
            )

def check_1d_component_closure(
    h_total,
    components,
    label,
    tol=1e-8,
    raise_on_failure=False,
):
    """
    Check whether a list of 1D component histograms reproduces h_total
    bin by bin.

    The check includes regular bins only. It reports both absolute and
    fractional differences and optionally raises on failed closure.
    """
    if not components:
        raise RuntimeError(
            "{}: no component histograms were provided".format(label)
        )

    n_bins = h_total.GetNbinsX()

    for component in components:
        if component.GetNbinsX() != n_bins:
            raise RuntimeError(
                "{}: bin-count mismatch: total has {} bins, "
                "component {} has {} bins".format(
                    label,
                    n_bins,
                    component.GetName(),
                    component.GetNbinsX(),
                )
            )

    print("\n===== 1D COMPONENT CLOSURE: {} =====".format(label))

    max_abs_diff = 0.0
    max_frac_diff = 0.0
    failed_bins = []

    for i in range(1, n_bins + 1):
        total = h_total.GetBinContent(i)
        component_values = [
            component.GetBinContent(i)
            for component in components
        ]
        component_sum = sum(component_values)

        diff = component_sum - total
        frac = (
            diff / total
            if abs(total) > 1e-12
            else float("nan")
        )

        max_abs_diff = max(max_abs_diff, abs(diff))

        if np.isfinite(frac):
            max_frac_diff = max(max_frac_diff, abs(frac))

        print(
            "bin {:2d}: total={:12.6g} "
            "components={} sum={:12.6g} "
            "diff={:12.6g} frac={:12.6g}".format(
                i,
                total,
                ["{:.6g}".format(value) for value in component_values],
                component_sum,
                diff,
                frac,
            )
        )

        scale = max(1.0, abs(total))

        if abs(diff) > tol * scale:
            failed_bins.append(i)

    print("max absolute difference =", max_abs_diff)
    print("max fractional difference =", max_frac_diff)
    print("failed bins =", failed_bins)

    if failed_bins and raise_on_failure:
        raise RuntimeError(
            "{} failed component closure in bins {}".format(
                label,
                failed_bins,
            )
        )

    return failed_bins

def check_flux_universe_component_closure(
    h_total,
    components,
    label,
    band_name="Flux",
    tol=1e-8,
    max_universes_to_print=5,
    raise_on_failure=False,
):
    """
    Check whether component histograms reproduce h_total bin by bin
    for every universe in a vertical error band.

    Expected:
        sum_k component_k^(u)(i) == total^(u)(i)

    for every universe u and regular bin i.
    """

    print(
        "\n===== FLUX UNIVERSE COMPONENT CLOSURE: {} =====".format(
            label
        )
    )

    # Make sure the requested band exists everywhere.
    hists = [h_total] + list(components)

    for h in hists:
        band_names = [str(x) for x in h.GetVertErrorBandNames()]

        if band_name not in band_names:
            raise RuntimeError(
                "{} does not contain vertical error band '{}'. "
                "Available bands = {}".format(
                    h.GetName(),
                    band_name,
                    band_names,
                )
            )

    total_band = h_total.GetVertErrorBand(band_name)
    component_bands = [
        h.GetVertErrorBand(band_name)
        for h in components
    ]

    n_universes = total_band.GetNHists()

    print("band       =", band_name)
    print("universes  =", n_universes)

    for component, band in zip(components, component_bands):
        if band.GetNHists() != n_universes:
            raise RuntimeError(
                "{} has {} {} universes; total has {}".format(
                    component.GetName(),
                    band.GetNHists(),
                    band_name,
                    n_universes,
                )
            )

    failed = []
    max_abs_diff = 0.0
    max_frac_diff = 0.0

    for u in range(n_universes):
        h_total_u = total_band.GetHist(u)
        component_u = [
            band.GetHist(u)
            for band in component_bands
        ]

        universe_failed = False
        universe_max_abs = 0.0
        universe_max_frac = 0.0

        for i in range(1, h_total.GetNbinsX() + 1):
            total = h_total_u.GetBinContent(i)

            component_values = [
                h.GetBinContent(i)
                for h in component_u
            ]

            component_sum = sum(component_values)

            diff = component_sum - total

            frac = (
                diff / total
                if abs(total) > 1e-12
                else float("nan")
            )

            universe_max_abs = max(
                universe_max_abs,
                abs(diff),
            )

            max_abs_diff = max(
                max_abs_diff,
                abs(diff),
            )

            if np.isfinite(frac):
                universe_max_frac = max(
                    universe_max_frac,
                    abs(frac),
                )

                max_frac_diff = max(
                    max_frac_diff,
                    abs(frac),
                )

            scale = max(1.0, abs(total))

            if abs(diff) > tol * scale:
                universe_failed = True

                failed.append(
                    {
                        "universe": u,
                        "bin": i,
                        "total": total,
                        "components": component_values,
                        "sum": component_sum,
                        "diff": diff,
                        "frac": frac,
                    }
                )

        if (
            u < max_universes_to_print
            or universe_failed
        ):
            print(
                "univ {:4d}: max_abs_diff={:12.6g} "
                "max_frac_diff={:12.6g} {}".format(
                    u,
                    universe_max_abs,
                    universe_max_frac,
                    "FAIL" if universe_failed else "OK",
                )
            )

    print("")
    print("overall max absolute difference =", max_abs_diff)
    print("overall max fractional difference =", max_frac_diff)
    print("failed universe/bin pairs =", len(failed))

    if failed:
        print("\nFirst failed entries:")

        for entry in failed[:10]:
            print(
                "  univ {:4d} bin {:2d}: "
                "total={:12.6g} components={} "
                "sum={:12.6g} diff={:12.6g} "
                "frac={:12.6g}".format(
                    entry["universe"],
                    entry["bin"],
                    entry["total"],
                    [
                        "{:.6g}".format(x)
                        for x in entry["components"]
                    ],
                    entry["sum"],
                    entry["diff"],
                    entry["frac"],
                )
            )

        if raise_on_failure:
            raise RuntimeError(
                "{} failed {} universe closure in {} "
                "universe/bin pairs".format(
                    label,
                    band_name,
                    len(failed),
                )
            )

    else:
        print(
            "{} universes close exactly within tolerance.".format(
                band_name
            )
        )

    return failed

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

def clone_category_sum(hist_holder, categories, name):
    output = hist_holder.hists["Total"].Clone(name)
    output.Reset()
    output.SetDirectory(0)

    used_categories = []

    for category in categories:
        if category not in hist_holder.hists:
            raise RuntimeError(
                "Missing required signal category '{}' while building {}. "
                "Available categories: {}".format(
                    category,
                    name,
                    sorted(hist_holder.hists.keys()),
                )
            )

        output.Add(hist_holder.hists[category])
        used_categories.append(category)

    print("")
    print("===== SIGNAL CATEGORY SUM: {} =====".format(name))
    print("used categories =", used_categories)
    print("integral        =", output.Integral())

    return output

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
    # type_path_map = {
    #     "data": "/exp/minerva/data/users/qvuong/elastic_nue/kin_dist_dataleFHC_NuE_allSystematics_fullStatsFluxes_MAD.root",
    #     "mc":   "/exp/minerva/data/users/qvuong/elastic_nue/kin_dist_mcleFHC_NuE_allSystematics_fullStatsFluxes_MAD.root",
    # }
    type_path_map = {
        "data": stitch_path("fhc_elastic", "data"),
        "mc":   stitch_path("fhc_elastic", "mc"),
    }

    data_file, mc_file, pot_scale, data_pot, mc_pot = Utilities.getFilesAndPOTScale(
        "NuE_allSystematics_fullStatsFluxes", type_path_map, "MAD", True
    )

    standPOT = data_pot if data_pot is not None else mc_pot
    print_pot_scale_check("FHC NuE elastic", data_pot, mc_pot, standPOT)

    # tuned_file = ROOT.TFile.Open(
    #     "/exp/minerva/data/users/qvuong/elastic_nue/bkgfit_leFHC_nueElastic_matrix_NuE_allSystematics_fullStatsFluxes_MAD.root"
    # )
    tuned_file = ROOT.TFile.Open(
        stitch_path("fhc_elastic", "bkgfit")
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

    # -------------------------------------------------
    # Check whether the production flavor decomposition
    # reproduces the tuned elastic signal bin by bin.
    # Do not modify the components yet.
    # -------------------------------------------------
    h_sum = h_electron.Clone("fhc_elastic_flavor_sum")
    h_sum.SetDirectory(0)
    h_sum.Add(h_muon)

    print("\n===== FHC ELASTIC FLAVOR INTEGRAL CLOSURE =====")
    print("tuned MC integral       =", h_mc.Integral())
    print("electron integral       =", h_electron.Integral())
    print("muon integral           =", h_muon.Integral())
    print("electron + muon         =", h_sum.Integral())
    print(
        "flavor sum - tuned MC  =",
        h_sum.Integral() - h_mc.Integral(),
    )

    elastic_closure_failed_bins = check_1d_component_closure(
        h_mc,
        [h_electron, h_muon],
        "FHC elastic tuned MC vs electron+muon",
        tol=1e-8,
        raise_on_failure=False,
    )

    if elastic_closure_failed_bins:
        print(
            "\nWARNING: FHC elastic flavor components do not close "
            "to the tuned MC bin by bin."
        )
        print(
            "No correction has been applied. Inspect the differences "
            "before deciding whether to normalize the components."
        )
    else:
        print(
            "\nFHC elastic flavor components close to the tuned MC "
            "bin by bin. No additional normalization is needed."
        )

    # -------------------------------------------------
    # Check the same electron+muon decomposition in every
    # Flux universe.
    # -------------------------------------------------
    elastic_flux_closure_failures = (
        check_flux_universe_component_closure(
            h_mc,
            [h_electron, h_muon],
            "FHC elastic tuned MC vs electron+muon",
            band_name="Flux",
            tol=1e-8,
            max_universes_to_print=5,
            raise_on_failure=False,
        )
    )

    if elastic_flux_closure_failures:
        print(
            "\nWARNING: FHC elastic electron+muon components do not "
            "reproduce the tuned MC Flux universes exactly."
        )
        print(
            "Do not modify them yet; inspect the size/pattern of the "
            "differences first."
        )
    else:
        print(
            "\nFHC elastic electron+muon components also close to the "
            "tuned MC for every Flux universe."
        )

    # -------------------------------------------------
    # 2D flavor templates used only to obtain the
    # conditional L/E distribution in each reco-energy bin.
    # The path is selected through stitch_input_files.json.
    # Overall template normalization cancels in the
    # relative oscillation weighting.
    # -------------------------------------------------
    template_2d_path = stitch_path(
        "fhc_elastic",
        "template_2d",
    )

    template_file_2d = ROOT.TFile.Open(
        template_2d_path
    )

    if not template_file_2d or template_file_2d.IsZombie():
        raise RuntimeError(
            "Could not open LElike 4-bin elastic template file: {}".format(
                template_2d_path
            )
        )

    print("\nLoading FHC elastic 2D templates from LElike 4-bin file")
    print("  file =", template_2d_path)

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



def normalize_template_reco_slices(
    h2,
    h1_target,
    out_name,
    reco_axis="y",
):
    """
    Normalize each reconstructed-energy slice of h2 so that its
    projection equals h1_target bin by bin.

    reco_axis="y":
        x = true L/E_nu
        y = reconstructed E_lep

    reco_axis="x":
        x = reconstructed E_lep
        y = true L/E_nu
    """
    h2_out = h2.Clone(out_name)
    h2_out.SetDirectory(0)

    if reco_axis == "y":
        n_reco = h2_out.GetNbinsY()
        n_osc = h2_out.GetNbinsX()
    elif reco_axis == "x":
        n_reco = h2_out.GetNbinsX()
        n_osc = h2_out.GetNbinsY()
    else:
        raise RuntimeError(
            "Unknown reco_axis: {}".format(reco_axis)
        )

    if n_reco != h1_target.GetNbinsX():
        raise RuntimeError(
            "{} reco bins do not match {} target bins".format(
                n_reco,
                h1_target.GetNbinsX(),
            )
        )

    for i_reco in range(1, n_reco + 1):
        target = h1_target.GetBinContent(i_reco)

        current = 0.0

        for j_osc in range(1, n_osc + 1):
            if reco_axis == "y":
                current += h2_out.GetBinContent(
                    j_osc,
                    i_reco,
                )
            else:
                current += h2_out.GetBinContent(
                    i_reco,
                    j_osc,
                )

        if current == 0.0:
            if abs(target) > 1e-12:
                raise RuntimeError(
                    "{}: reco bin {} has zero template content "
                    "but target prediction is {}".format(
                        out_name,
                        i_reco,
                        target,
                    )
                )
            continue

        scale = target / current

        for j_osc in range(1, n_osc + 1):
            if reco_axis == "y":
                xbin = j_osc
                ybin = i_reco
            else:
                xbin = i_reco
                ybin = j_osc

            content = h2_out.GetBinContent(xbin, ybin)
            error = h2_out.GetBinError(xbin, ybin)

            h2_out.SetBinContent(
                xbin,
                ybin,
                content * scale,
            )
            h2_out.SetBinError(
                xbin,
                ybin,
                error * abs(scale),
            )

    return h2_out


def check_template_projection(
    h2,
    h1_target,
    label,
    reco_axis="y",
    tol=1e-8,
):
    if reco_axis == "y":
        # x = oscillation variable (true L/E)
        # y = reconstructed energy
        #
        # Use regular x bins only, matching
        # normalize_template_reco_slices().
        projection = h2.ProjectionY(
            label + "_reco_projection",
            1,
            h2.GetNbinsX(),
        )

    elif reco_axis == "x":
        # y = oscillation variable
        # x = reconstructed energy
        #
        # Use regular y bins only.
        projection = h2.ProjectionX(
            label + "_reco_projection",
            1,
            h2.GetNbinsY(),
        )

    else:
        raise RuntimeError(
            "Unknown reco_axis: {}".format(reco_axis)
        )

    projection.SetDirectory(0)

    if projection.GetNbinsX() != h1_target.GetNbinsX():
        raise RuntimeError(
            "{} projection has {} bins; target has {}".format(
                label,
                projection.GetNbinsX(),
                h1_target.GetNbinsX(),
            )
        )

    print(
        "\n===== TEMPLATE PROJECTION CHECK: {} =====".format(
            label
        )
    )

    for i in range(1, h1_target.GetNbinsX() + 1):
        target = h1_target.GetBinContent(i)
        actual = projection.GetBinContent(i)
        diff = actual - target

        print(
            "bin {:2d}: target={:12.6g} "
            "projection={:12.6g} diff={:12.6g}".format(
                i,
                target,
                actual,
                diff,
            )
        )

        scale = max(1.0, abs(target))

        if abs(diff) > tol * scale:
            raise RuntimeError(
                "{} projection closure failed in bin {}: "
                "target={}, actual={}".format(
                    label,
                    i,
                    target,
                    actual,
                )
            )

def check_2d_reco_edges(
    h2,
    expected_edges,
    label,
    reco_axis="y",
    tol=1e-8,
):
    axis = h2.GetYaxis() if reco_axis == "y" else h2.GetXaxis()

    if axis.GetNbins() != len(expected_edges) - 1:
        raise RuntimeError(
            "{} reco axis has {} bins, expected {}".format(
                label,
                axis.GetNbins(),
                len(expected_edges) - 1,
            )
        )

    actual_edges = [
        axis.GetBinLowEdge(i)
        for i in range(1, axis.GetNbins() + 1)
    ]
    actual_edges.append(axis.GetBinUpEdge(axis.GetNbins()))

    for i, (actual, expected) in enumerate(
        zip(actual_edges, expected_edges)
    ):
        if abs(actual - expected) > tol:
            raise RuntimeError(
                "{} edge {} mismatch: actual={}, expected={}".format(
                    label,
                    i,
                    actual,
                    expected,
                )
            )



def load_jaewon_data_json(path):
    """
    Load Jaewon's published nu+e elastic data and non-flux covariance
    directly from JSON.

    Returns:
        h_data : TH1D
        cov    : TMatrixD
        config : decoded JSON dictionary
    """

    with open(path, "r") as f:
        config = json.load(f)

    edges = config["binning"]["analysis_edges"]
    data = config["data"]
    errors = config["data_errors"]
    covariance = config["covariance"]

    n_bins = len(data)

    if len(edges) != n_bins + 1:
        raise RuntimeError(
            "Jaewon JSON has {} data bins but {} edges".format(
                n_bins,
                len(edges),
            )
        )

    if len(errors) != n_bins:
        raise RuntimeError(
            "Jaewon JSON data_errors length {} does not match "
            "{} data bins".format(
                len(errors),
                n_bins,
            )
        )

    if len(covariance) != n_bins:
        raise RuntimeError(
            "Jaewon JSON covariance has {} rows, expected {}".format(
                len(covariance),
                n_bins,
            )
        )

    for i, row in enumerate(covariance):
        if len(row) != n_bins:
            raise RuntimeError(
                "Jaewon JSON covariance row {} has {} entries, "
                "expected {}".format(
                    i,
                    len(row),
                    n_bins,
                )
            )

    # -------------------------------------------------
    # Data histogram
    # -------------------------------------------------
    h_data = PlotUtils.MnvH1D(
        "fhc_elastic_data",
        "Published #nu-e elastic data",
        n_bins,
        array("d", edges),
    )
    h_data.SetDirectory(0)

    for i in range(n_bins):
        h_data.SetBinContent(
            i + 1,
            float(data[i]),
        )
        h_data.SetBinError(
            i + 1,
            float(errors[i]),
        )

    # -------------------------------------------------
    # Published non-flux covariance
    # -------------------------------------------------
    cov = ROOT.TMatrixD(
        n_bins,
        n_bins,
    )

    for i in range(n_bins):
        for j in range(n_bins):
            cov[i][j] = float(
                covariance[i][j]
            )

    print("\n===== JAEWON NUE ELASTIC INPUT =====")
    print("file =", path)
    print("POT  =", config["pot"])
    print("edges =", edges)
    print("data  =", data)
    print("errors =", errors)

    print("\nCovariance:")
    for i in range(n_bins):
        print(
            "  " +
            " ".join(
                "{:12.6g}".format(cov[i][j])
                for j in range(n_bins)
            )
        )

    return h_data, cov, config


# This is my flux-folded v+e prediction and Jaewon's result
def load_nue_elastic_fhc_prediction():
    """
    Construct the FHC nu+e elastic sample using:

    data/covariance:
        Jaewon's published 6-bin result from JSON,
        rescaled from the paper POT to the full FHC analysis POT.
        The covariance is the published non-flux covariance.

    nominal MC and flavor decomposition:
        Our flux-folded higher-order P8 nu+e prediction,
        already normalized to the full FHC analysis data POT.

    oscillation templates:
        High-statistics nu+e 2D templates providing the true-L/E
        distribution within each reconstructed electron-energy bin.

    The 2D templates are normalized reco-bin by reco-bin to the
    corresponding flux-folded flavor prediction.
    """

    # =================================================
    # Jaewon published data + non-flux covariance
    # =================================================

    h_data, cov, jaewon_config = load_jaewon_data_json(
        JAEWON_INPUT_FILE
    )

    jaewon_pot = float(
        jaewon_config["pot"]
    )

    print("\n===== JAEWON POT CHECK =====")
    print("JSON POT     =", jaewon_pot)
    print("Expected POT =", PAPER_NUE_ELASTIC_POT)

    if abs(
        jaewon_pot - PAPER_NUE_ELASTIC_POT
    ) > 1e-6 * PAPER_NUE_ELASTIC_POT:
        raise RuntimeError(
            "Jaewon JSON POT mismatch: "
            "JSON={} expected={}".format(
                jaewon_pot,
                PAPER_NUE_ELASTIC_POT,
            )
        )

    # =================================================
    # Scale Jaewon's published data and covariance
    # from the paper POT to the FHC analysis POT.
    #
    # Counts/errors:  x r
    # Covariance:     x r^2
    # =================================================

    jaewon_pot_scale = FHC_ANALYSIS_POT / jaewon_pot

    print("\n===== JAEWON POT RESCALING =====")
    print("paper POT        =", jaewon_pot)
    print("analysis POT     =", FHC_ANALYSIS_POT)
    print("data scale       =", jaewon_pot_scale)
    print("covariance scale =", jaewon_pot_scale ** 2)

    print("data integral before =", h_data.Integral())

    # ROOT histogram scaling scales both contents and bin errors.
    h_data.Scale(jaewon_pot_scale)

    # Published covariance is absolute covariance in events^2.
    scale_tmatrix(
        cov,
        jaewon_pot_scale,
    )

    print("\nScaled Jaewon covariance:")
    for i in range(cov.GetNrows()):
        print(
            "  " +
            " ".join(
                "{:12.6g}".format(cov[i][j])
                for j in range(cov.GetNcols())
            )
        )

    print("data integral after  =", h_data.Integral())

    check_1d_edges(
        h_data,
        JAEWON_ELEP_EDGES,
        "Jaewon data",
    )

    print_hist_bins(
        h_data,
        "Jaewon published nu+e data",
    )

    # =================================================
    # New flux-folded prediction
    # =================================================

    prediction_file = ROOT.TFile.Open(
        PREDICTED_ELASTIC_FILE
    )

    if (
        not prediction_file
        or prediction_file.IsZombie()
    ):
        raise RuntimeError(
            "Could not open predicted elastic file: {}".format(
                PREDICTED_ELASTIC_FILE
            )
        )

    h_mc = get_hist_checked(
        prediction_file,
        ["h_nue_elastic_total"],
        "fhc_elastic_mc",
    )

    h_nue = get_hist_checked(
        prediction_file,
        ["h_nue_elastic_nue"],
        "fhc_elastic_component_nue",
    )

    h_nuebar = get_hist_checked(
        prediction_file,
        ["h_nue_elastic_nuebar"],
        "fhc_elastic_component_nuebar",
    )

    h_numu = get_hist_checked(
        prediction_file,
        ["h_nue_elastic_numu"],
        "fhc_elastic_component_numu",
    )

    h_numubar = get_hist_checked(
        prediction_file,
        ["h_nue_elastic_numubar"],
        "fhc_elastic_component_numubar",
    )

    # =================================================
    # Check prediction binning
    # =================================================

    for label, hist in [
        ("total prediction", h_mc),
        ("nue prediction", h_nue),
        ("nuebar prediction", h_nuebar),
        ("numu prediction", h_numu),
        ("numubar prediction", h_numubar),
    ]:
        check_1d_edges(
            hist,
            JAEWON_ELEP_EDGES,
            label,
        )

    # =================================================
    # Check P8 Flux error bands
    # =================================================

    print("\n===== NUE ELASTIC FLUX BAND CHECK =====")

    prediction_hists = [
        ("total", h_mc),
        ("nue", h_nue),
        ("nuebar", h_nuebar),
        ("numu", h_numu),
        ("numubar", h_numubar),
    ]

    for label, hist in prediction_hists:

        band_names = [
            str(x)
            for x in hist.GetVertErrorBandNames()
        ]

        if "Flux" not in band_names:
            raise RuntimeError(
                "{} prediction has no Flux band. "
                "Available bands = {}".format(
                    label,
                    band_names,
                )
            )

        n_flux = (
            hist
            .GetVertErrorBand("Flux")
            .GetNHists()
        )

        print(
            "{:8s}: {} Flux universes".format(
                label,
                n_flux,
            )
        )

        if n_flux != 1000:
            raise RuntimeError(
                "{} prediction has {} Flux universes; "
                "expected 1000".format(
                    label,
                    n_flux,
                )
            )

    # =================================================
    # Build electron- and muon-flavor groups
    #
    # electron:
    #     nue + nuebar
    #
    # muon:
    #     numu + numubar
    # =================================================

    h_electron = h_nue.Clone(
        "electron_fhc_elastic"
    )
    h_electron.SetDirectory(0)
    h_electron.Add(h_nuebar)

    h_muon = h_numu.Clone(
        "muon_fhc_elastic"
    )
    h_muon.SetDirectory(0)
    h_muon.Add(h_numubar)

    # =================================================
    # CV component closure
    # =================================================

    h_prediction_sum = h_electron.Clone(
        "fhc_elastic_prediction_sum_check"
    )
    h_prediction_sum.SetDirectory(0)
    h_prediction_sum.Add(h_muon)

    print(
        "\n===== NEW NUE ELASTIC PREDICTION CLOSURE ====="
    )

    print("total prediction =", h_mc.Integral())
    print("nue              =", h_nue.Integral())
    print("nuebar           =", h_nuebar.Integral())
    print("numu             =", h_numu.Integral())
    print("numubar          =", h_numubar.Integral())
    print("electron         =", h_electron.Integral())
    print("muon             =", h_muon.Integral())
    print("flavor sum       =", h_prediction_sum.Integral())
    print(
        "flavor sum-total =",
        h_prediction_sum.Integral()
        - h_mc.Integral(),
    )

    print(
        "\n===== PREDICTION BIN-BY-BIN CLOSURE ====="
    )

    for i in range(
        1,
        h_mc.GetNbinsX() + 1,
    ):

        total = h_mc.GetBinContent(i)

        flavor_sum = (
            h_prediction_sum.GetBinContent(i)
        )

        diff = flavor_sum - total

        print(
            "bin {:2d}: total={:12.6g} "
            "flavor_sum={:12.6g} "
            "diff={:12.6g}".format(
                i,
                total,
                flavor_sum,
                diff,
            )
        )

        if (
            abs(diff)
            > 1e-8 * max(1.0, abs(total))
        ):
            raise RuntimeError(
                "Prediction flavor closure failed "
                "in bin {}".format(i)
            )

    # =================================================
    # Flux-universe component closure
    #
    # These should close by construction because the
    # total prediction was made from the same four
    # flavor predictions universe by universe.
    # =================================================

    check_flux_universe_component_closure(
        h_mc,
        [
            h_nue,
            h_nuebar,
            h_numu,
            h_numubar,
        ],
        "FHC elastic total vs four flavors",
        band_name="Flux",
        tol=1e-8,
        max_universes_to_print=5,
        raise_on_failure=True,
    )

    check_flux_universe_component_closure(
        h_mc,
        [
            h_electron,
            h_muon,
        ],
        "FHC elastic total vs electron+muon",
        band_name="Flux",
        tol=1e-8,
        max_universes_to_print=5,
        raise_on_failure=True,
    )

    # =================================================
    # Load high-statistics 2D oscillation templates
    # =================================================

    template_file = ROOT.TFile.Open(
        LELIKE_6BIN_ELASTIC_TEMPLATE_FILE
    )

    if (
        not template_file
        or template_file.IsZombie()
    ):
        raise RuntimeError(
            "Could not open 6-bin elastic template file: {}".format(
                LELIKE_6BIN_ELASTIC_TEMPLATE_FILE
            )
        )

    h2_nue_raw = get_hist_checked(
        template_file,
        ["drawnL_ElepReco_LE_NuEElasticE"],
        "fhc_elastic_template_nue_raw",
    )

    h2_nuebar_raw = get_hist_checked(
        template_file,
        ["drawnL_ElepReco_LE_NuEElasticEBar"],
        "fhc_elastic_template_nuebar_raw",
    )

    h2_numu_raw = get_hist_checked(
        template_file,
        ["drawnL_ElepReco_LE_NuEElasticMu"],
        "fhc_elastic_template_numu_raw",
    )

    h2_numubar_raw = get_hist_checked(
        template_file,
        ["drawnL_ElepReco_LE_NuEElasticMuBar"],
        "fhc_elastic_template_numubar_raw",
    )

    # =================================================
    # Identify reconstructed-energy axis
    # =================================================

    print("\nElastic template axes:")
    print(
        "x =",
        h2_nue_raw.GetXaxis().GetTitle(),
    )
    print(
        "y =",
        h2_nue_raw.GetYaxis().GetTitle(),
    )

    # Current template definition:
    #
    # x = true L/E
    # y = reconstructed electron energy
    #
    reco_axis = "y"

    # =================================================
    # Match numerical final-bin edge
    #
    # The published final bin is physically 9-inf.
    #
    # Older templates used 120 GeV as the numerical
    # endpoint, while the new P8 flux prediction uses
    # 100 GeV.  We only change the numerical axis edge;
    # no bin content is changed.
    # =================================================

    for label, hist in [
        ("nue template", h2_nue_raw),
        ("nuebar template", h2_nuebar_raw),
        ("numu template", h2_numu_raw),
        ("numubar template", h2_numubar_raw),
    ]:
        check_2d_reco_edges(
            hist,
            JAEWON_ELEP_EDGES,
            label,
            reco_axis,
        )

    # =================================================
    # Normalize each reco-energy slice of each flavor's
    # L/E template to the new flux-folded prediction.
    # =================================================

    h2_nue = normalize_template_reco_slices(
        h2_nue_raw,
        h_nue,
        "fhc_elastic_template_nue",
        reco_axis,
    )

    h2_nuebar = normalize_template_reco_slices(
        h2_nuebar_raw,
        h_nuebar,
        "fhc_elastic_template_nuebar",
        reco_axis,
    )

    h2_numu = normalize_template_reco_slices(
        h2_numu_raw,
        h_numu,
        "fhc_elastic_template_numu",
        reco_axis,
    )

    h2_numubar = normalize_template_reco_slices(
        h2_numubar_raw,
        h_numubar,
        "fhc_elastic_template_numubar",
        reco_axis,
    )

    # =================================================
    # Check that each normalized template projects back
    # exactly to its 1D prediction.
    # =================================================

    check_template_projection(
        h2_nue,
        h_nue,
        "nue",
        reco_axis,
    )

    check_template_projection(
        h2_nuebar,
        h_nuebar,
        "nuebar",
        reco_axis,
    )

    check_template_projection(
        h2_numu,
        h_numu,
        "numu",
        reco_axis,
    )

    check_template_projection(
        h2_numubar,
        h_numubar,
        "numubar",
        reco_axis,
    )

    # =================================================
    # Group 2D templates for the oscillation machinery
    # =================================================

    h_template_electron = h2_nue.Clone(
        "fhc_elastic_template_electron"
    )
    h_template_electron.SetDirectory(0)
    h_template_electron.Add(h2_nuebar)

    h_template_muon = h2_numu.Clone(
        "fhc_elastic_template_muon"
    )
    h_template_muon.SetDirectory(0)
    h_template_muon.Add(h2_numubar)

    # Final grouped-template projection checks.
    check_template_projection(
        h_template_electron,
        h_electron,
        "electron flavor",
        reco_axis,
    )

    check_template_projection(
        h_template_muon,
        h_muon,
        "muon flavor",
        reco_axis,
    )

    print("\n===== NUE ELASTIC EXPOSURE CHECK =====")
    print("Jaewon original POT =", jaewon_pot)
    print("analysis POT        =", FHC_ANALYSIS_POT)
    print("Jaewon scale        =", jaewon_pot_scale)
    print("prediction POT      =", FHC_ANALYSIS_POT)

    # =================================================
    # Final diagnostic summary
    # =================================================

    print("\n===== FINAL FHC NUE ELASTIC INPUT =====")

    print(
        "Jaewon data integral      =",
        h_data.Integral(),
    )

    print(
        "P8 prediction integral    =",
        h_mc.Integral(),
    )

    print(
        "electron component        =",
        h_electron.Integral(),
    )

    print(
        "muon component            =",
        h_muon.Integral(),
    )

    print(
        "electron + muon            =",
        h_electron.Integral()
        + h_muon.Integral(),
    )

    print(
        "electron template integral =",
        h_template_electron.Integral(),
    )

    print(
        "muon template integral     =",
        h_template_muon.Integral(),
    )

    print(
        "external covariance size   = "
        "{} x {}".format(
            cov.GetNrows(),
            cov.GetNcols(),
        )
    )

    # Histograms were cloned/detached above, so the files
    # can now safely be closed.
    prediction_file.Close()
    template_file.Close()

    return {
        "mc": h_mc,
        "data": h_data,

        "electron": h_electron,
        "muon": h_muon,

        "nue": h_nue,
        "nuebar": h_nuebar,
        "numu": h_numu,
        "numubar": h_numubar,

        "template_electron": h_template_electron,
        "template_muon": h_template_muon,

        "template_nue": h2_nue,
        "template_nuebar": h2_nuebar,
        "template_numu": h2_numu,
        "template_numubar": h2_numubar,

        "cov": cov,
    }


def load_ccnue_fhc():
    type_path_map = {
        "data": stitch_path("fhc_ccnue", "data"),
        "mc":   stitch_path("fhc_ccnue", "mc"),
    }

    data_file, mc_file, pot_scale, data_pot, mc_pot = Utilities.getFilesAndPOTScale(
        "CCnue_allSystematics_fullStatsFluxes", type_path_map, "MAD", True
    )

    standPOT = data_pot if data_pot is not None else mc_pot
    binwidthScale = getattr(AnalysisConfig, "binwidth", False)

    print_pot_scale_check("FHC CCnue", data_pot, mc_pot, standPOT)

    tuned_file = ROOT.TFile.Open(
        stitch_path("fhc_ccnue", "bkgfit")
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
        stitch_path("fhc_ccnue", "mc")
    )

    template_holder = HistHolder(
        "Reco Energy vs L/E",
        template_file,
        "Signal",
        True,
        mc_pot,
        standPOT
    )

    # Full nominal template before scaling, for diagnostics.
    h_template_total_raw = clone_total(
        template_holder,
        "fhc_ccnue_template_total_raw",
    )

    # Apply POT/bin-width scaling once.
    template_holder.POTScale(binwidthScale)

    # Full nominal template after scaling, for diagnostics.
    h_template_total = clone_total(
        template_holder,
        "fhc_ccnue_template_total",
    )

    # Physics template: signal categories only.
    h_template = clone_category_sum(
        template_holder,
        SIGNAL_DEFINITION,
        "fhc_ccnue_template_signal",
    )

    h_template.SetDirectory(0)

    print_integral_change(
        "FHC CCnue template total",
        h_template_total_raw,
        h_template_total,
    )

    print("")
    print("===== FHC NOMINAL TEMPLATE TOTAL VS SIGNAL =====")
    print("total selected integral =", h_template_total.Integral())
    print("signal-only integral    =", h_template.Integral())

    swap_type_path_map = {
        "mc": stitch_path("fhc_ccnue", "swap_mc")
    }

    _, swap_mc_file, _, _, swap_mc_pot = Utilities.getFilesAndPOTScale(
        "CCnueswap_allSystematics_fullStatsFluxes", swap_type_path_map, "MAD", True
    )

    print_pot_scale_check("FHC CCnue swap", None, swap_mc_pot, standPOT)

    swap_file = ROOT.TFile.Open(
        stitch_path("fhc_ccnue", "swap_mc")
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

    # Diagnostic full selected swap before scaling.
    h_swap_total_raw = clone_total(
        swap_hist_holder,
        "fhc_ccnue_swap_total_raw",
    )

    h_swap_template_total_raw = clone_total(
        swap_template_holder,
        "fhc_ccnue_swap_template_total_raw",
    )

    swap_hist_holder.POTScale(binwidthScale)
    swap_template_holder.POTScale(binwidthScale)

    # Diagnostic full selected swap after scaling.
    h_swap_total = clone_total(
        swap_hist_holder,
        "fhc_ccnue_swap_total",
    )

    h_swap_template_total = clone_total(
        swap_template_holder,
        "fhc_ccnue_swap_template_total",
    )

    # Physics objects: signal categories only.
    h_swap = clone_category_sum(
        swap_hist_holder,
        SIGNAL_DEFINITION,
        "fhc_ccnue_swap_signal",
    )

    h_swap_template = clone_category_sum(
        swap_template_holder,
        SIGNAL_DEFINITION,
        "fhc_ccnue_swap_template_signal",
    )

    print("")
    print("===== FHC SWAP TOTAL VS SIGNAL =====")
    print("total selected integral =", h_swap_total.Integral())
    print("signal-only integral    =", h_swap.Integral())

    for i in range(1, h_swap.GetNbinsX() + 1):
        total = h_swap_total.GetBinContent(i)
        signal = h_swap.GetBinContent(i)

        fraction = signal / total if abs(total) > 1e-12 else float("nan")

        print(
            "bin {:2d}: total={:12.6g} signal={:12.6g} "
            "signal/total={:10.6f}".format(
                i,
                total,
                signal,
                fraction,
            )
        )

    h_swap.SetDirectory(0)
    h_swap_template.SetDirectory(0)

    print_integral_change(
        "FHC CCnue swap total 1D",
        h_swap_total_raw,
        h_swap_total,
    )

    print_integral_change(
        "FHC CCnue swap total template",
        h_swap_template_total_raw,
        h_swap_template_total,
    )

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
        "data": stitch_path("rhc_ccnuebar", "data"),
        "mc":   stitch_path("rhc_ccnuebar", "mc"),
    }

    data_file, mc_file, pot_scale, data_pot, mc_pot = Utilities.getFilesAndPOTScale(
        "CCnuebar_allSystematics_fullStatsFluxes", type_path_map, "MAD", True
    )

    standPOT = data_pot if data_pot is not None else mc_pot
    binwidthScale = getattr(AnalysisConfig, "binwidth", False)

    print_pot_scale_check("RHC CCnuebar", data_pot, mc_pot, standPOT)

    tuned_file = ROOT.TFile.Open(
        stitch_path("rhc_ccnuebar", "bkgfit")
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

    template_file = ROOT.TFile.Open(
        stitch_path("rhc_ccnuebar", "mc")
    )

    template_holder = HistHolder(
        "Reco Energy vs L/E",
        template_file,
        "Signal",
        True,
        mc_pot,
        standPOT
    )

    # Full nominal template before scaling, for diagnostics.
    h_template_total_raw = clone_total(
        template_holder,
        "rhc_ccnuebar_template_total_raw",
    )

    # Apply POT/bin-width scaling once.
    template_holder.POTScale(binwidthScale)

    # Full nominal template after scaling, for diagnostics.
    h_template_total = clone_total(
        template_holder,
        "rhc_ccnuebar_template_total",
    )

    # Physics template: signal categories only.
    h_template = clone_category_sum(
        template_holder,
        SIGNAL_DEFINITION,
        "rhc_ccnuebar_template_signal",
    )

    h_template.SetDirectory(0)

    print_integral_change(
        "RHC CCnuebar template total",
        h_template_total_raw,
        h_template_total,
    )

    print("")
    print("===== RHC NOMINAL TEMPLATE TOTAL VS SIGNAL =====")
    print("total selected integral =", h_template_total.Integral())
    print("signal-only integral    =", h_template.Integral())

    swap_type_path_map = {
        "mc": stitch_path("rhc_ccnuebar", "swap_mc")
    }

    _, swap_mc_file, _, _, swap_mc_pot = Utilities.getFilesAndPOTScale(
        "CCnuebarswap_allSystematics_fullStatsFluxes", swap_type_path_map, "MAD", True
    )

    print_pot_scale_check("RHC CCnuebar swap", None, swap_mc_pot, standPOT)

    swap_file = ROOT.TFile.Open(
        stitch_path("rhc_ccnuebar", "swap_mc")
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

    # Full selected swap before POT scaling.
    h_swap_total_raw = clone_total(
        swap_hist_holder,
        "rhc_ccnuebar_swap_total_raw",
    )

    h_swap_template_total_raw = clone_total(
        swap_template_holder,
        "rhc_ccnuebar_swap_template_total_raw",
    )

    # Apply POT/bin-width scaling exactly once.
    swap_hist_holder.POTScale(binwidthScale)
    swap_template_holder.POTScale(binwidthScale)

    # Full selected swap after scaling, for diagnostics.
    h_swap_total = clone_total(
        swap_hist_holder,
        "rhc_ccnuebar_swap_total",
    )

    h_swap_template_total = clone_total(
        swap_template_holder,
        "rhc_ccnuebar_swap_template_total",
    )

    # Signal-only physics objects.
    h_swap = clone_category_sum(
        swap_hist_holder,
        SIGNAL_DEFINITION,
        "rhc_ccnuebar_swap_signal",
    )

    h_swap_template = clone_category_sum(
        swap_template_holder,
        SIGNAL_DEFINITION,
        "rhc_ccnuebar_swap_template_signal",
    )

    h_swap.SetDirectory(0)
    h_swap_template.SetDirectory(0)

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
        "data": stitch_path("fhc_ccnumu", "data"),
        "mc":   stitch_path("fhc_ccnumu", "mc"),
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
        "data": stitch_path("rhc_ccnumubar", "data"),
        "mc":   stitch_path("rhc_ccnumubar", "mc"),
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

    print("\n===== RHC CCnumubar 2p2h keys =====")
    for key in mc_file.GetListOfKeys():
        name = key.GetName()
        if "2p2h" in name.lower():
            obj = mc_file.Get(name)
            if obj and hasattr(obj, "Integral"):
                print(
                    "{:<60s} integral = {}".format(
                        name,
                        obj.Integral(),
                    )
                )
            else:
                print(name)


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

    elastic_source = getattr(
        AnalysisConfig,
        "elastic_source",
        "production",
    )

    print("Loading FHC nue elastic...")
    print("elastic_source =", elastic_source)

    if elastic_source in ["paper", "jaewon", "prediction"]:
        fhc_elastic = load_nue_elastic_fhc_prediction()

    elif elastic_source in ["production", "prod", "mine"]:
        fhc_elastic = load_nue_elastic_fhc()

    else:
        raise RuntimeError(
            "Unknown elastic_source: {}".format(
                elastic_source
            )
        )

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

    input_set = getattr(AnalysisConfig, "stitch_input_set", "p6")

    n_flux_universes = {
        "p6": 100,
        "p8": 1000,
        "p8_onlyPPFX": 1000,
        "p8_onlyBeamFocus": 1000,
    }[input_set]

    sample_histogram.SetNFluxUniverses(n_flux_universes)

    print(
        "Stitch input set = {}, using {} flux universes".format(
            input_set,
            n_flux_universes,
        )
    )

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

    # csv_dir = "{}/oscillations/csvs/{}".format(ccnueroot, hist_config_tag)

    # if not os.path.exists(csv_dir):
    #     os.makedirs(csv_dir)

    # mc_cv = np.array(stitched.mc_hist)[1:-1]
    # data_cv = np.array(stitched.data_hist)[1:-1]

    # cov_full = stitched.GetCovarianceMatrix(False)
    # cov_sans = stitched.GetCovarianceMatrix(True)
    # cov_flux = cov_full - cov_sans

    # np.savetxt(
    #     "{}/mc_cv.csv".format(csv_dir),
    #     mc_cv,
    #     delimiter=",",
    #     header="mc_cv",
    #     comments="",
    # )

    # np.savetxt(
    #     "{}/data_cv.csv".format(csv_dir),
    #     data_cv,
    #     delimiter=",",
    #     header="data_cv",
    #     comments="",
    # )

    # np.savetxt(
    #     "{}/cov_full.csv".format(csv_dir),
    #     cov_full,
    #     delimiter=",",
    # )

    # np.savetxt(
    #     "{}/cov_sansFlux.csv".format(csv_dir),
    #     cov_sans,
    #     delimiter=",",
    # )

    # np.savetxt(
    #     "{}/cov_flux.csv".format(csv_dir),
    #     cov_flux,
    #     delimiter=",",
    # )

    # np.savetxt(
    #     "{}/inv_cov_full.csv".format(csv_dir),
    #     stitched.GetInverseCovarianceMatrix(False),
    #     delimiter=",",
    # )

    # np.savetxt(
    #     "{}/inv_cov_sansFlux.csv".format(csv_dir),
    #     stitched.GetInverseCovarianceMatrix(True),
    #     delimiter=",",
    # )

    # np.savetxt(
    #     "{}/A_flux_universes_minus_cv.csv".format(csv_dir),
    #     stitched.GetAMatrix(),
    #     delimiter=",",
    # )

    # print("Wrote stitched CSV outputs to", csv_dir)

    stitched.Write(OUTROOT)
    print("Wrote stitched file to", OUTROOT)

    # c = ROOT.TCanvas("c", "c", 900, 700)
    # c.SetLogy()
    # stitched.mc_hist.SetLineWidth(2)
    # stitched.mc_hist.Draw("HIST")
    # stitched.data_hist.Draw("E1 SAME")
    # c.Print("stitched_{}.png".format(hist_config_tag))

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
_stitch_parser.add_argument("--elastic-era", choices=["LE", "ME"], default=None)

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
    # choices=["p6", "p8", "p8_onlyPPFX", "p8_onlyBeamFocus", "p8_oldResponse"],
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

PREDICTED_ELASTIC_FILE = (
    "/exp/minerva/data/users/qvuong/nueel_prediction_studies/nue_elastic_prediction_TomalakXS_fixedTotalNorm_newElectronResponse_mnv_FHC.root"
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
        if not np.isfinite(target) or target < 0:
            raise RuntimeError("Invalid target prediction in reco bin {}".format(i_reco))

        current = 0.0

        for j_osc in range(1, n_osc + 1):
            value = (h2_out.GetBinContent(j_osc, i_reco)
                     if reco_axis == "y"
                     else h2_out.GetBinContent(i_reco, j_osc))
            if not np.isfinite(value) or value < 0:
                raise RuntimeError(
                    "{} has invalid template content in reco bin {}".format(out_name, i_reco)
                )
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


def load_nue_elastic_prediction(era="LE", beam="FHC"):
    """Load a flux-folded prediction, measured spectrum, and L/E templates.

    LE: Jaewon's JSON measurement and drawnL_ElepReco_LE_* templates.
    ME: paper_data_{FHC,RHC} and ElepReco_LE_* templates.
    The 2D inputs supply CV conditional L/E shapes, not energy migration.
    Prediction files retain their paired Flux universes.
    """
    era, beam = str(era).upper(), str(beam).upper()
    if era not in ("LE", "ME") or beam not in ("FHC", "RHC"):
        raise ValueError("Invalid elastic era/beam: {}/{}".format(era, beam))
    if era == "LE" and beam != "FHC":
        raise ValueError("The current LE measurement supports FHC only")

    sample = beam.lower() + "_elastic"
    cfg = STITCH_INPUTS[sample]
    edges = JAEWON_ELEP_EDGES
    categories = {
        "nue": "NuEElasticE", "nuebar": "NuEElasticEBar",
        "numu": "NuEElasticMu", "numubar": "NuEElasticMuBar",
    }

    def open_input(path):
        fin = ROOT.TFile.Open(path, "READ")
        if not fin or fin.IsZombie():
            raise RuntimeError("Cannot open elastic input: " + path)
        if fin.TestBit(ROOT.TFile.kRecovered):
            fin.Close()
            raise RuntimeError("Recovered elastic input: " + path)
        print("[INPUT]", path)
        return fin

    def positive_scale(key, default):
        value = float(cfg.get(key, default))
        if not np.isfinite(value) or value <= 0:
            raise ValueError("{} must be finite and positive".format(key))
        return value

    prediction_path = (
        cfg.get("prediction", PREDICTED_ELASTIC_FILE)
        if era == "LE" else stitch_path(sample, "prediction")
    )
    prefix = cfg.get(
        "prediction_hist_prefix",
        "h_nue_elastic_" if era == "LE" else "h_nue_elastic_reco_",
    )
    if prefix not in ("h_nue_elastic_", "h_nue_elastic_reco_"):
        raise ValueError("Unsupported prediction_hist_prefix: " + prefix)
    if era == "ME" and prefix == "h_nue_elastic_":
        if cfg.get("allow_truth_prediction") is not True:
            raise RuntimeError(
                "ME h_nue_elastic_* are truth-energy predictions. "
                "Supply response-folded h_nue_elastic_reco_*; "
                "allow_truth_prediction=true permits a diagnostic only."
            )
        print("WARNING: ME truth-energy diagnostic; no energy response applied.")

    prediction_scale = positive_scale("prediction_pot_scale", 1.0)
    prediction_file = open_input(prediction_path)
    joint_cov = None
    try:
        prediction = {
            flavor: get_hist_checked(
                prediction_file, [prefix + flavor],
                sample + "_" + ("mc" if flavor == "total" else flavor),
            )
            for flavor in ("total",) + tuple(categories)
        }
        for hist in prediction.values():
            hist.Scale(prediction_scale)

        if era == "LE":
            h_data, cov, measurement = load_jaewon_data_json(
                cfg.get("measurement_json", JAEWON_INPUT_FILE)
            )
            paper_pot = float(measurement["pot"])
            if not np.isclose(paper_pot, PAPER_NUE_ELASTIC_POT, rtol=1e-6, atol=0):
                raise RuntimeError("Unexpected LE measurement POT: {}".format(paper_pot))
            data_scale = positive_scale(
                "data_pot_scale", FHC_ANALYSIS_POT / paper_pot
            )
            h_data.Scale(data_scale)
            scale_tmatrix(cov, data_scale)
        else:
            data_path = cfg.get("measurement_file", prediction_path)
            data_file = open_input(data_path)
            try:
                raw_data = get_hist_checked(
                    data_file, ["paper_data_" + beam], sample + "_data_raw"
                )
                # The published ME objects are TH1D, while stitched data use MnvH1D.
                h_data = PlotUtils.MnvH1D(raw_data)
                h_data.SetName(sample + "_data")
                h_data.SetDirectory(0)
                matrix = data_file.Get("cov_nonflux_ME_FHC_RHC")
                if not matrix or matrix.GetNrows() != 12 or matrix.GetNcols() != 12:
                    raise RuntimeError("Missing/invalid joint 12x12 ME measurement covariance")
                joint_cov = matrix.Clone(sample + "_joint_cov_raw")
            finally:
                data_file.Close()
            data_scale = positive_scale("data_pot_scale", 1.0)
            h_data.Scale(data_scale)
            cov = None
    finally:
        prediction_file.Close()

    h_mc = prediction.pop("total")
    for label, hist in [("data", h_data), ("total", h_mc)] + list(prediction.items()):
        check_1d_edges(hist, edges, sample + "_" + label)
        values = [hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)]
        if not np.all(np.isfinite(values)) or np.any(np.array(values) < 0):
            raise RuntimeError("Nonfinite/negative elastic spectrum: " + label)

    for label, hist in [("total", h_mc)] + list(prediction.items()):
        if "Flux" not in [str(x) for x in hist.GetVertErrorBandNames()]:
            raise RuntimeError(sample + " " + label + " has no Flux band")
        if hist.GetVertErrorBand("Flux").GetNHists() != 1000:
            raise RuntimeError(sample + " " + label + " must have 1000 Flux universes")
    check_1d_component_closure(
        h_mc, list(prediction.values()), sample + " CV flavor closure",
        raise_on_failure=True,
    )
    check_flux_universe_component_closure(
        h_mc, list(prediction.values()), sample + " Flux flavor closure",
        raise_on_failure=True,
    )

    template_file = open_input(stitch_path(sample, "prediction_template_2d"))
    templates = {}
    try:
        template_prefix = "drawnL_ElepReco_LE_" if era == "LE" else "ElepReco_LE_"
        for flavor, category in categories.items():
            raw = get_hist_checked(
                template_file, [template_prefix + category],
                sample + "_template_" + flavor + "_raw",
            )
            if raw.GetDimension() != 2:
                raise RuntimeError(raw.GetName() + " is not two-dimensional")
            # Oscillation templates use CV shapes only. Do not retain unscaled
            # MnvH2D universe bands while normalizing its CV slices.
            if hasattr(raw, "GetCVHistoWithStatError"):
                raw = raw.GetCVHistoWithStatError()
                raw.SetDirectory(0)
            check_2d_reco_edges(raw, edges, sample + "_" + flavor, reco_axis="y")
            templates[flavor] = normalize_template_reco_slices(
                raw, prediction[flavor], sample + "_template_" + flavor,
                reco_axis="y",
            )
            check_template_projection(
                templates[flavor], prediction[flavor],
                sample + "_" + flavor, reco_axis="y",
            )
    finally:
        template_file.Close()

    def combine(first, second, name):
        hist = first.Clone(name)
        hist.SetDirectory(0)
        hist.Add(second)
        return hist

    electron = combine(prediction["nue"], prediction["nuebar"], "electron_" + sample)
    muon = combine(prediction["numu"], prediction["numubar"], "muon_" + sample)
    template_electron = combine(
        templates["nue"], templates["nuebar"], sample + "_template_electron"
    )
    template_muon = combine(
        templates["numu"], templates["numubar"], sample + "_template_muon"
    )
    check_template_projection(template_electron, electron, sample + "_electron")
    check_template_projection(template_muon, muon, sample + "_muon")
    print("Elastic era/beam:", era, beam)
    print("Measurement scale:", data_scale, "prediction scale:", prediction_scale)
    print_hist_bins(h_data, sample + " measured data")
    print_hist_bins(h_mc, sample + " predicted MC")

    return dict(
        mc=h_mc, data=h_data, electron=electron, muon=muon,
        template_electron=template_electron, template_muon=template_muon,
        cov=cov, joint_cov=joint_cov, data_scale=data_scale,
        **prediction,
        **{"template_" + flavor: hist for flavor, hist in templates.items()}
    )


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

    elastic_era = str(
        _stitch_args.elastic_era or STITCH_INPUTS.get("elastic_era", "LE")
    ).upper()
    if elastic_era not in ("LE", "ME"):
        raise ValueError("elastic_era must be LE or ME")
    # The legacy production-normalized loader lives in stitch_legacy.py.
    elastic_beams = ("FHC",) if elastic_era == "LE" else ("FHC", "RHC")
    if elastic_era == "ME" and not getattr(
        StitchedHistogram, "ME_ELASTIC_JOINT_COVARIANCE_SUPPORTED", False
    ):
        raise RuntimeError(
            "ME stitching requires the joint ME elastic covariance update "
            "in tools/StitchedHistogram.py. Refusing to omit FHC/RHC correlations."
        )
    elastic_samples = {
        beam.lower() + "_elastic": load_nue_elastic_prediction(elastic_era, beam)
        for beam in elastic_beams
    }

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
    for sample, payload in elastic_samples.items():
        for component in ("mc", "data", "electron", "muon"):
            print(sample, component, "integral:", payload[component].Integral())

    sample_histogram = StitchedHistogram("sample")

    input_set = getattr(AnalysisConfig, "stitch_input_set", "p6")

    if input_set == "p6":
        n_flux_universes = 100
    elif input_set.startswith(("p8", "me_p8")):
        n_flux_universes = 1000
    else:
        raise RuntimeError(
            "Unknown stitch_input_set for flux universes: {}".format(input_set)
        )

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

    for sample, payload in elastic_samples.items():
        sample_histogram.AddScatteringFlavors("electron_" + sample, payload["electron"])
        sample_histogram.AddScatteringFlavors("muon_" + sample, payload["muon"])
        sample_histogram.AddHistograms(sample, payload["mc"], payload["data"])
        sample_histogram.AddTemplates(
            sample,
            nue=payload["template_electron"],
            numu=payload["template_muon"],
            swap=payload["template_muon"],
        )

    # Add CC templates before making ratios / copying.
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

    if not elastic_excluded:
        if elastic_era == "ME":
            fhc = elastic_samples["fhc_elastic"]
            rhc = elastic_samples["rhc_elastic"]
            cov = fhc["joint_cov"].Clone("me_elastic_paper_covariance_matrix")
            for i in range(12):
                for j in range(12):
                    if not np.isclose(
                        cov[i][j], rhc["joint_cov"][i][j], rtol=1e-10, atol=1e-12
                    ):
                        raise RuntimeError("FHC/RHC joint ME measurement covariance differs")
            # Published ordering: FHC six bins followed by RHC six bins.
            scales = [fhc["data_scale"]] * 6 + [rhc["data_scale"]] * 6
            for i in range(12):
                for j in range(12):
                    cov[i][j] *= scales[i] * scales[j]
            stitched.external_covariances["me_elastic"] = cov
        else:
            cov = elastic_samples["fhc_elastic"]["cov"]
            stitched.external_covariances["fhc_elastic"] = cov.Clone(
                "fhc_elastic_paper_covariance_matrix"
            )

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
    elastic_order = [name for name in ("fhc_elastic", "rhc_elastic")
                     if name in elastic_samples]
    if ratio_mode:
        desired_order = elastic_order + [
            "fhc_ratio",
            "rhc_ratio",
            "fhc_numu_selection",
            "rhc_numu_selection",
        ]
    else:
        desired_order = elastic_order + [
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


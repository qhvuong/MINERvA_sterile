#!/usr/bin/env python3

import argparse
import math
import os

import ROOT


ROOT.TH1.AddDirectory(False)


FILES = {
    "fhc": {
        "nominal": (
            "/exp/minerva/data/users/qvuong/nu_e/"
            "kin_dist_mcleFHC_p8Tuples_CCnue_updatedFluxes_MAD.root"
        ),
        "bkgfit": (
            "/exp/minerva/data/users/qvuong/nu_e/"
            "bkgfit_leFHC_N4_tune_p8Tuples_CCnue_updatedFluxes_MAD.root"
        ),
        "swap": (
            "/exp/minerva/data/users/qvuong/nu_e_swapped/"
            "kin_dist_mcleFHC_p8Tuples_CCnueswap_updatedFluxes_MAD.root"
        ),
        "data_pot": 3.331982991676676e20,
        "mc_pot": 2.2166978413209162e21,
        "expected_bkgfit_integral": 2565.0640736705986,
        "expected_bkgfit_bin1": 101.15006,
        "expected_swap_bin1": 417.12764,
    },
    "rhc": {
        "nominal": (
            "/exp/minerva/data/users/qvuong/antinu_e/"
            "kin_dist_mcle5_p8Tuples_CCnuebar_updatedFluxes_MAD.root"
        ),
        "bkgfit": (
            "/exp/minerva/data/users/qvuong/antinu_e/"
            "bkgfit_le5_N4_tune_p8Tuples_CCnuebar_updatedFluxes_MAD.root"
        ),
        "swap": (
            "/exp/minerva/data/users/qvuong/antinu_e_swapped/"
            "kin_dist_mcle5_p8Tuples_CCnuebarswap_updatedFluxes_MAD.root"
        ),
        # Replace with the RHC data POT printed by stitch.py if needed.
        "data_pot": None,
        "mc_pot": 9.513217104576512e20,
        "expected_bkgfit_integral": None,
        "expected_bkgfit_bin1": 16.202755,
        "expected_swap_bin1": 37.351749,
    },
}


def open_file(path):
    root_file = ROOT.TFile.Open(path)

    if not root_file or root_file.IsZombie():
        raise RuntimeError("Could not open {}".format(path))

    return root_file


def clone_hist(root_file, name, clone_name):
    obj = root_file.Get(name)

    if not obj:
        return None

    if not obj.InheritsFrom("TH1"):
        return None

    hist = obj.Clone(clone_name)
    hist.SetDirectory(0)

    return hist


def integral(hist):
    return hist.Integral(1, hist.GetNbinsX())


def print_hist_summary(label, hist):
    print("")
    print(label)
    print("  name     =", hist.GetName())
    print("  class    =", hist.ClassName())
    print("  nbins    =", hist.GetNbinsX())
    print("  integral =", integral(hist))

    for index in range(1, hist.GetNbinsX() + 1):
        print(
            "  bin {:2d}: [{:7.3f}, {:7.3f}] "
            "content={:14.8g} error={:14.8g}".format(
                index,
                hist.GetXaxis().GetBinLowEdge(index),
                hist.GetXaxis().GetBinUpEdge(index),
                hist.GetBinContent(index),
                hist.GetBinError(index),
            )
        )


def list_histograms(path):
    root_file = open_file(path)
    results = []

    for key in root_file.GetListOfKeys():
        obj = key.ReadObj()

        if not obj or not obj.InheritsFrom("TH1"):
            continue

        hist = obj.Clone(
            "diag_{}_{}".format(
                os.path.basename(path).replace(".root", ""),
                key.GetName(),
            )
        )
        hist.SetDirectory(0)

        results.append(
            {
                "name": key.GetName(),
                "class": key.GetClassName(),
                "hist": hist,
                "nbins": hist.GetNbinsX(),
                "integral": integral(hist),
                "bin1": hist.GetBinContent(1),
            }
        )

    root_file.Close()

    return results


def candidate_score(
    item,
    expected_integral=None,
    expected_bin1=None,
    expected_nbins=12,
):
    score = 0.0

    if item["nbins"] != expected_nbins:
        return float("inf")

    if expected_integral is not None:
        denominator = max(abs(expected_integral), 1.0)
        score += abs(item["integral"] - expected_integral) / denominator

    if expected_bin1 is not None:
        denominator = max(abs(expected_bin1), 1.0)
        score += abs(item["bin1"] - expected_bin1) / denominator

    return score


def print_bkgfit_candidates(cfg):
    print("")
    print("=" * 100)
    print("BACKGROUND-FIT HISTOGRAM SEARCH")
    print("file =", cfg["bkgfit"])

    items = list_histograms(cfg["bkgfit"])

    ranked = sorted(
        items,
        key=lambda item: candidate_score(
            item,
            expected_integral=cfg["expected_bkgfit_integral"],
            expected_bin1=cfg["expected_bkgfit_bin1"],
        ),
    )

    print("")
    print(
        "{:4s} {:45s} {:8s} {:16s} {:16s} {:14s}".format(
            "rank",
            "name",
            "nbins",
            "integral",
            "bin1",
            "score",
        )
    )

    for rank, item in enumerate(ranked[:20], start=1):
        score = candidate_score(
            item,
            expected_integral=cfg["expected_bkgfit_integral"],
            expected_bin1=cfg["expected_bkgfit_bin1"],
        )

        print(
            "{:4d} {:45s} {:8d} {:16.8g} {:16.8g} {:14.8g}".format(
                rank,
                item["name"],
                item["nbins"],
                item["integral"],
                item["bin1"],
                score,
            )
        )

    best = ranked[0]

    print("")
    print("Best matching background-fit histogram:")
    print("  key      =", best["name"])
    print("  integral =", best["integral"])
    print("  bin 1    =", best["bin1"])

    return best["name"], best["hist"]


def inspect_nominal_and_swap(cfg):
    nominal_file = open_file(cfg["nominal"])
    swap_file = open_file(cfg["swap"])

    nominal_en4 = clone_hist(
        nominal_file,
        "EN4",
        "nominal_EN4",
    )

    nominal_signal = clone_hist(
        nominal_file,
        "EN4_CCNuE",
        "nominal_EN4_CCNuE",
    )

    swap_en4 = clone_hist(
        swap_file,
        "EN4",
        "swap_EN4",
    )

    swap_signal = clone_hist(
        swap_file,
        "EN4_CCNuE",
        "swap_EN4_CCNuE",
    )

    nominal_file.Close()
    swap_file.Close()

    for label, hist in [
        ("Raw nominal EN4", nominal_en4),
        ("Raw nominal EN4_CCNuE", nominal_signal),
        ("Raw swapped EN4", swap_en4),
        ("Raw swapped EN4_CCNuE", swap_signal),
    ]:
        if hist:
            print_hist_summary(label, hist)

    return nominal_en4, nominal_signal, swap_en4, swap_signal


def compare_pipeline(
    cfg,
    bkgfit_hist,
    swap_en4,
):
    print("")
    print("=" * 100)
    print("PIPELINE CHECK")

    if cfg["data_pot"] is None:
        print("No data POT configured; skipping explicit POT scaling.")
        return

    pot_scale = cfg["data_pot"] / cfg["mc_pot"]

    print("data POT =", cfg["data_pot"])
    print("MC POT   =", cfg["mc_pot"])
    print("POT scale=", pot_scale)

    swap_scaled = swap_en4.Clone("swap_EN4_pot_scaled")
    swap_scaled.SetDirectory(0)
    swap_scaled.Scale(pot_scale)

    print("")
    print(
        "{:>6s} {:>16s} {:>16s} {:>16s}".format(
            "bin",
            "tuned Nue",
            "scaled Nswap",
            "Nswap/Nue",
        )
    )

    for index in range(1, bkgfit_hist.GetNbinsX() + 1):
        nue = bkgfit_hist.GetBinContent(index)
        swap = swap_scaled.GetBinContent(index)

        ratio = (
            swap / nue
            if abs(nue) > 1e-12
            else float("nan")
        )

        print(
            "{:6d} {:16.8g} {:16.8g} {:16.8f}".format(
                index,
                nue,
                swap,
                ratio,
            )
        )

    print("")
    print("First-bin closure checks:")
    print(
        "  tuned Nue bin 1       =",
        bkgfit_hist.GetBinContent(1),
    )
    print(
        "  scaled Nswap bin 1    =",
        swap_scaled.GetBinContent(1),
    )
    print(
        "  scaled Nswap/Nue bin1 =",
        swap_scaled.GetBinContent(1)
        / bkgfit_hist.GetBinContent(1),
    )

    if cfg["expected_bkgfit_bin1"] is not None:
        print(
            "  expected tuned Nue    =",
            cfg["expected_bkgfit_bin1"],
        )

    if cfg["expected_swap_bin1"] is not None:
        print(
            "  expected scaled Nswap =",
            cfg["expected_swap_bin1"],
        )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--mode",
        choices=["fhc", "rhc"],
        default="fhc",
    )

    parser.add_argument(
        "--bkgfit-hist",
        default=None,
        help=(
            "Explicit background-fit histogram key. "
            "When omitted, search for the object matching stitch.py output."
        ),
    )

    args = parser.parse_args()
    cfg = FILES[args.mode]

    print("")
    print("=" * 100)
    print("{} STITCHING DIAGNOSTIC".format(args.mode.upper()))

    print("")
    print("nominal =", cfg["nominal"])
    print("bkgfit  =", cfg["bkgfit"])
    print("swap    =", cfg["swap"])

    nominal_en4, nominal_signal, swap_en4, swap_signal = (
        inspect_nominal_and_swap(cfg)
    )

    if args.bkgfit_hist:
        bkgfit_file = open_file(cfg["bkgfit"])

        bkgfit_hist = clone_hist(
            bkgfit_file,
            args.bkgfit_hist,
            "selected_bkgfit_hist",
        )

        bkgfit_file.Close()

        if not bkgfit_hist:
            raise RuntimeError(
                "Could not load {} from {}".format(
                    args.bkgfit_hist,
                    cfg["bkgfit"],
                )
            )

        selected_key = args.bkgfit_hist
    else:
        selected_key, bkgfit_hist = print_bkgfit_candidates(cfg)

    print("")
    print("Selected bkgfit key =", selected_key)
    print_hist_summary(
        "Selected tuned Nue histogram",
        bkgfit_hist,
    )

    if swap_en4 is None:
        raise RuntimeError("Missing EN4 in swapped file")

    compare_pipeline(
        cfg,
        bkgfit_hist,
        swap_en4,
    )


if __name__ == "__main__":
    main()
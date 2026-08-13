#!/usr/bin/env python3

import argparse
import os

import ROOT


ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetOptStat(0)


BASE_DIR = "/exp/minerva/data/users/qvuong"


# Replace these defaults with the exact updatedFluxes files used by stitch.py.
SAMPLES = {
    "fhc": {
        "nue": (
            "/exp/minerva/data/users/qvuong/nu_e/"
            "bkgfit_leFHC_N4_tune_p8Tuples_CCnue_updatedFluxes_MAD.root"
        ),
        "swap": (
            "/exp/minerva/data/users/qvuong/nu_e_swapped/"
            "kin_dist_mcleFHC_p8Tuples_CCnueswap_updatedFluxes_MAD.root"
        ),
        "label": "FHC",
    },
    "rhc": {
        "nue": (
            "/exp/minerva/data/users/qvuong/antinu_e/"
            "bkgfit_le5_N4_tune_p8Tuples_CCnuebar_updatedFluxes_MAD.root"
        ),
        "swap": (
            "/exp/minerva/data/users/qvuong/antinu_e_swapped/"
            "kin_dist_mcle5_p8Tuples_CCnuebarswap_updatedFluxes_MAD.root"
        ),
        "label": "RHC",
    },
}


def get_hist(path, hname, clone_name):
    root_file = ROOT.TFile.Open(path)

    if not root_file or root_file.IsZombie():
        raise RuntimeError("Could not open {}".format(path))

    hist = root_file.Get(hname)

    if not hist:
        print("")
        print("Could not find exact histogram:", hname)
        print("File:", path)
        print("Available matching keys:")

        matches = []

        for key in root_file.GetListOfKeys():
            name = key.GetName()

            if hname.lower() in name.lower():
                matches.append(name)

        if matches:
            for name in matches:
                print("  ", name)
        else:
            print("  No matching keys found")

        root_file.Close()

        raise RuntimeError(
            "Missing histogram {} in {}".format(hname, path)
        )

    hist = hist.Clone(clone_name)
    hist.SetDirectory(0)

    root_file.Close()

    return hist


def make_ratio(numerator, denominator, name):
    ratio = numerator.Clone(name)
    ratio.SetDirectory(0)

    ratio.Divide(
        numerator,
        denominator,
        1.0,
        1.0,
        "",
    )

    return ratio


def make_bin_number_hist(source, name):
    nbins = source.GetNbinsX()

    result = ROOT.TH1D(
        name,
        name,
        nbins,
        0.5,
        nbins + 0.5,
    )
    result.SetDirectory(0)

    for index in range(1, nbins + 1):
        result.SetBinContent(
            index,
            source.GetBinContent(index),
        )
        result.SetBinError(
            index,
            source.GetBinError(index),
        )

    return result


def print_bin_values(label, nue_hist, swap_hist, ratio_hist):
    print("")
    print("===== {} Nswap / Nue =====".format(label))

    print(
        "{:>6s} {:>16s} {:>16s} {:>16s} "
        "{:>14s} {:>14s}".format(
            "bin",
            "E low",
            "E high",
            "Nue",
            "Nswap",
            "Nswap/Nue",
        )
    )

    for index in range(1, ratio_hist.GetNbinsX() + 1):
        low_edge = ratio_hist.GetXaxis().GetBinLowEdge(index)
        high_edge = ratio_hist.GetXaxis().GetBinUpEdge(index)

        nue = nue_hist.GetBinContent(index)
        swap = swap_hist.GetBinContent(index)
        ratio = ratio_hist.GetBinContent(index)

        print(
            "{:6d} {:16.8g} {:16.8g} {:16.8g} "
            "{:14.8g} {:14.8f}".format(
                index,
                low_edge,
                high_edge,
                nue,
                swap,
                ratio,
            )
        )


def draw_ratio(
    ratio,
    title,
    xtitle,
    outname,
    ymin,
    ymax,
):
    canvas = ROOT.TCanvas(
        "c_{}".format(
            os.path.basename(outname).replace(".png", "")
        ),
        "c",
        900,
        700,
    )

    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.05)
    canvas.SetBottomMargin(0.13)
    canvas.SetTopMargin(0.08)

    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetMarkerColor(ROOT.kBlack)
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerSize(1.1)
    ratio.SetLineWidth(2)

    ratio.SetTitle(
        "{};{};N_{{swap}} / N_{{e}}".format(
            title,
            xtitle,
        )
    )

    ratio.GetYaxis().SetRangeUser(ymin, ymax)
    ratio.Draw("E1")

    unity = ROOT.TLine(
        ratio.GetXaxis().GetXmin(),
        1.0,
        ratio.GetXaxis().GetXmax(),
        1.0,
    )
    unity.SetLineColor(ROOT.kRed + 1)
    unity.SetLineStyle(2)
    unity.SetLineWidth(2)
    unity.Draw("SAME")

    canvas.Print(outname)
    canvas.Print(outname.replace(".png", ".pdf"))

    print("Wrote", outname)
    print("Wrote", outname.replace(".png", ".pdf"))


def resolve_sample_paths(args, sample_key):
    cfg = dict(SAMPLES[sample_key])

    if sample_key == "fhc":
        if args.fhc_nue:
            cfg["nue"] = args.fhc_nue

        if args.fhc_swap:
            cfg["swap"] = args.fhc_swap

    elif sample_key == "rhc":
        if args.rhc_nue:
            cfg["nue"] = args.rhc_nue

        if args.rhc_swap:
            cfg["swap"] = args.rhc_swap

    return cfg


def get_pot_scale(sample_key, args):
    """
    Both samples should be scaled to one common target POT.

    Since the common target cancels in Nswap/Nue,

        (Nswap_raw / POTswap) / (Nue_raw / POTnue)
        = (Nswap_raw / Nue_raw) * POTnue / POTswap.

    We therefore scale:
        Nue  -> Nue
        swap -> swap * POTnue/POTswap
    """

    if sample_key == "fhc":
        nue_pot = args.fhc_nue_pot
        swap_pot = args.fhc_swap_pot
    else:
        nue_pot = args.rhc_nue_pot
        swap_pot = args.rhc_swap_pot

    if nue_pot is None and swap_pot is None:
        return 1.0

    if nue_pot is None or swap_pot is None:
        raise RuntimeError(
            "For {}, provide both nominal and swapped MC POT values".format(
                sample_key
            )
        )

    if nue_pot <= 0.0 or swap_pot <= 0.0:
        raise RuntimeError("MC POT values must be positive")

    return nue_pot / swap_pot


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot the swapped-electron / nominal-electron MC ratio "
            "using the updated-flux selection files."
        )
    )

    parser.add_argument(
        "--hist",
        default="EN4",
        help="Histogram name in both nominal and swapped files",
    )

    parser.add_argument(
        "--outdir",
        default="swap_nue_ratio_plots",
    )

    parser.add_argument(
        "--samples",
        default="all",
        help="fhc, rhc, or all",
    )

    parser.add_argument(
        "--ymin",
        type=float,
        default=0.0,
    )

    parser.add_argument(
        "--ymax",
        type=float,
        default=5.0,
    )

    parser.add_argument("--fhc-nue")
    parser.add_argument("--fhc-swap")
    parser.add_argument("--rhc-nue")
    parser.add_argument("--rhc-swap")

    parser.add_argument(
        "--fhc-nue-pot",
        type=float,
        help="Generated POT represented by the FHC nominal CCnue MC",
    )
    parser.add_argument(
        "--fhc-swap-pot",
        type=float,
        help="Generated POT represented by the FHC swapped MC",
    )
    parser.add_argument(
        "--rhc-nue-pot",
        type=float,
        help="Generated POT represented by the RHC nominal CCnuebar MC",
    )
    parser.add_argument(
        "--rhc-swap-pot",
        type=float,
        help="Generated POT represented by the RHC swapped MC",
    )

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    if args.samples == "all":
        sample_keys = ["fhc", "rhc"]
    else:
        sample_keys = [
            item.strip()
            for item in args.samples.split(",")
            if item.strip()
        ]

    for sample_key in sample_keys:
        if sample_key not in SAMPLES:
            raise RuntimeError(
                "Unknown sample '{}'. Choose fhc, rhc, or all.".format(
                    sample_key
                )
            )

    output_root_path = os.path.join(
        args.outdir,
        "swap_nue_ratios.root",
    )

    output_root = ROOT.TFile(
        output_root_path,
        "RECREATE",
    )

    for sample_key in sample_keys:
        cfg = resolve_sample_paths(args, sample_key)

        print("")
        print("===== {} =====".format(sample_key))
        print("nominal Nue:", cfg["nue"])
        print("swapped:    ", cfg["swap"])

        nue_hist = get_hist(
            cfg["nue"],
            args.hist,
            "{}_nue_{}".format(sample_key, args.hist),
        )

        swap_hist = get_hist(
            cfg["swap"],
            args.hist,
            "{}_swap_{}".format(sample_key, args.hist),
        )

        if nue_hist.GetNbinsX() != swap_hist.GetNbinsX():
            raise RuntimeError(
                "{} nominal and swapped histograms have "
                "different bin counts: {} versus {}".format(
                    sample_key,
                    nue_hist.GetNbinsX(),
                    swap_hist.GetNbinsX(),
                )
            )

        pot_scale = get_pot_scale(sample_key, args)

        print("swap POT correction =", pot_scale)

        swap_hist.Scale(pot_scale)

        ratio_energy = make_ratio(
            swap_hist,
            nue_hist,
            "{}_swap_over_nue_energy".format(sample_key),
        )

        ratio_bin = make_bin_number_hist(
            ratio_energy,
            "{}_swap_over_nue_bin_number".format(sample_key),
        )

        print_bin_values(
            cfg["label"],
            nue_hist,
            swap_hist,
            ratio_energy,
        )

        title = "{} updated-flux swapped / nominal CC#nu_{{e}}".format(
            cfg["label"]
        )

        draw_ratio(
            ratio_energy,
            title,
            "E_{#nu}^{reco} [GeV]",
            os.path.join(
                args.outdir,
                "{}_swap_over_nue_energy.png".format(sample_key),
            ),
            args.ymin,
            args.ymax,
        )

        draw_ratio(
            ratio_bin,
            title,
            "Analysis bin number",
            os.path.join(
                args.outdir,
                "{}_swap_over_nue_bin_number.png".format(sample_key),
            ),
            args.ymin,
            args.ymax,
        )

        output_root.cd()

        nue_hist.Write(
            "{}_nue_{}".format(sample_key, args.hist)
        )
        swap_hist.Write(
            "{}_swap_{}_pot_scaled".format(
                sample_key,
                args.hist,
            )
        )
        ratio_energy.Write()
        ratio_bin.Write()

    output_root.Close()

    print("")
    print("Wrote", output_root_path)


if __name__ == "__main__":
    main()
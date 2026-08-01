#!/usr/bin/env python3

import argparse
import math
import os
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPaintTextFormat(".2f")


def find_mnvh1d(root_file, requested_name=None):
    """Return a requested MnvH1D, or list candidates when no name is given."""
    if requested_name:
        hist = root_file.Get(requested_name)
        if not hist:
            raise RuntimeError(
                f"Could not find histogram '{requested_name}' "
                f"in {root_file.GetName()}"
            )
        return hist

    candidates = []

    for key in root_file.GetListOfKeys():
        obj = key.ReadObj()

        if obj.InheritsFrom("PlotUtils::MnvH1D"):
            candidates.append(obj.GetName())

    if not candidates:
        raise RuntimeError(
            f"No PlotUtils::MnvH1D objects found in {root_file.GetName()}"
        )

    print(f"\nMnvH1D objects in {root_file.GetName()}:")
    for name in candidates:
        print(f"  {name}")

    raise RuntimeError(
        "Supply the histogram name using --hist-a and --hist-b."
    )


def get_universe_histograms(hist, band_name):
    """Get universe TH1 histograms from one vertical error band."""
    if not hist.HasVertErrorBand(band_name):
        available = [
            str(name)
            for name in hist.GetVertErrorBandNames()
        ]

        raise RuntimeError(
            f"Histogram '{hist.GetName()}' does not contain vertical band "
            f"'{band_name}'. Available bands: {available}"
        )

    band = hist.GetVertErrorBand(band_name)
    nuniv = band.GetNHists()

    return [band.GetHist(u) for u in range(nuniv)]


def make_cross_correlation(hist_a, hist_b, band_name):
    """
    Construct the cross-correlation matrix using

        C_ij = (1/N) sum_u
               (phi^u_A,i - phi^CV_A,i)
               (phi^u_B,j - phi^CV_B,j)

        rho_ij = C_ij / sqrt(C_A,ii * C_B,jj)

    Universe u in histogram A is paired with universe u in histogram B.
    """

    cv_a = hist_a.GetCVHistoWithStatError()
    cv_b = hist_b.GetCVHistoWithStatError()

    universes_a = get_universe_histograms(hist_a, band_name)
    universes_b = get_universe_histograms(hist_b, band_name)

    if len(universes_a) != len(universes_b):
        raise RuntimeError(
            f"Universe-count mismatch: "
            f"A has {len(universes_a)}, B has {len(universes_b)}"
        )

    nuniv = len(universes_a)
    nbins_a = cv_a.GetNbinsX()
    nbins_b = cv_b.GetNbinsX()

    print(f"\nUsing error band: {band_name}")
    print(f"Matched universes: {nuniv}")
    print(f"A energy bins:     {nbins_a}")
    print(f"B energy bins:     {nbins_b}")

    # Absolute deviations from the CV, matching the formula on the slide.
    deviations_a = [
        [
            universes_a[u].GetBinContent(i + 1)
            - cv_a.GetBinContent(i + 1)
            for i in range(nbins_a)
        ]
        for u in range(nuniv)
    ]

    deviations_b = [
        [
            universes_b[u].GetBinContent(j + 1)
            - cv_b.GetBinContent(j + 1)
            for j in range(nbins_b)
        ]
        for u in range(nuniv)
    ]

    print("\nFirst 5 universe fractional shifts in selected bins:")

    test_bins_a = [1, min(5, nbins_a), min(10, nbins_a)]
    test_bins_b = [1, min(5, nbins_b), min(10, nbins_b)]

    for u in range(min(5, nuniv)):
        print(f"\nUniverse {u}")

        for i in test_bins_a:
            cv = cv_a.GetBinContent(i)
            varied = universes_a[u].GetBinContent(i)
            frac = varied / cv - 1.0 if cv != 0.0 else 0.0

            print(
                f"  A bin {i:3d}: "
                f"CV={cv:.6e}, universe={varied:.6e}, "
                f"fractional shift={frac:+.6e}"
            )

        for j in test_bins_b:
            cv = cv_b.GetBinContent(j)
            varied = universes_b[u].GetBinContent(j)
            frac = varied / cv - 1.0 if cv != 0.0 else 0.0

            print(
                f"  B bin {j:3d}: "
                f"CV={cv:.6e}, universe={varied:.6e}, "
                f"fractional shift={frac:+.6e}"
            )

    # Diagonal covariance terms for each file separately.
    variance_a = [
        sum(
            deviations_a[u][i] ** 2
            for u in range(nuniv)
        ) / nuniv
        for i in range(nbins_a)
    ]

    variance_b = [
        sum(
            deviations_b[u][j] ** 2
            for u in range(nuniv)
        ) / nuniv
        for j in range(nbins_b)
    ]

    from array import array

    x_edges = [
        cv_a.GetXaxis().GetBinLowEdge(i + 1)
        for i in range(nbins_a)
    ]
    x_edges.append(cv_a.GetXaxis().GetBinUpEdge(nbins_a))

    y_edges = [
        cv_b.GetXaxis().GetBinLowEdge(j + 1)
        for j in range(nbins_b)
    ]
    y_edges.append(cv_b.GetXaxis().GetBinUpEdge(nbins_b))

    corr = ROOT.TH2D(
        "energy_correlation",
        "",
        nbins_a,
        array("d", x_edges),
        nbins_b,
        array("d", y_edges),
    )

    for i in range(nbins_a):
        for j in range(nbins_b):
            covariance = sum(
                deviations_a[u][i] * deviations_b[u][j]
                for u in range(nuniv)
            ) / nuniv

            denominator = math.sqrt(
                variance_a[i] * variance_b[j]
            )

            rho = covariance / denominator if denominator > 0.0 else 0.0

            corr.SetBinContent(i + 1, j + 1, rho)

    corr.SetMinimum(-0.4)
    corr.SetMaximum(0.4)

    zero_var_a = [
        i + 1 for i, value in enumerate(variance_a)
        if value <= 0.0
    ]

    zero_var_b = [
        j + 1 for j, value in enumerate(variance_b)
        if value <= 0.0
    ]

    if zero_var_a:
        print(f"WARNING: zero-variance bins in A: {zero_var_a}")

    if zero_var_b:
        print(f"WARNING: zero-variance bins in B: {zero_var_b}")

    values = [
        corr.GetBinContent(i + 1, j + 1)
        for i in range(nbins_a)
        for j in range(nbins_b)
    ]

    print(
        f"Correlation range: "
        f"min={min(values):+.6f}, max={max(values):+.6f}"
    )

    return corr


def main():
    parser = argparse.ArgumentParser(
        description="Plot a cross-file flux energy correlation matrix."
    )

    parser.add_argument("--file-a", required=True)
    parser.add_argument("--file-b", required=True)

    parser.add_argument(
        "--hist-a",
        default=None,
        help="MnvH1D name in file A",
    )
    parser.add_argument(
        "--hist-b",
        default=None,
        help="MnvH1D name in file B",
    )

    parser.add_argument(
        "--band",
        default="Flux",
        help="Vertical error-band name; default: Flux",
    )

    parser.add_argument(
        "--x-title",
        default="FHC neutrino energy [GeV]",
    )
    parser.add_argument(
        "--y-title",
        default="RHC neutrino energy [GeV]",
    )

    parser.add_argument(
        "--title",
        default="FHC–RHC flux energy correlation",
    )

    parser.add_argument(
        "--output",
        default="fhc_rhc_energy_correlation.png",
    )

    args = parser.parse_args()

    file_a = ROOT.TFile.Open(args.file_a, "READ")
    file_b = ROOT.TFile.Open(args.file_b, "READ")

    if not file_a or file_a.IsZombie():
        raise RuntimeError(f"Could not open {args.file_a}")

    if not file_b or file_b.IsZombie():
        raise RuntimeError(f"Could not open {args.file_b}")

    hist_a = find_mnvh1d(file_a, args.hist_a)
    hist_b = find_mnvh1d(file_b, args.hist_b)

    corr = make_cross_correlation(hist_a, hist_b, args.band)

    canvas = ROOT.TCanvas("canvas", "canvas", 1000, 850)
    canvas.SetRightMargin(0.15)
    canvas.SetLeftMargin(0.13)
    canvas.SetBottomMargin(0.12)

    corr.SetTitle(args.title)
    corr.GetXaxis().SetTitle(args.x_title)
    corr.GetYaxis().SetTitle(args.y_title)
    corr.GetZaxis().SetTitle("Correlation coefficient")

    corr.GetXaxis().CenterTitle()
    corr.GetYaxis().CenterTitle()
    corr.GetZaxis().CenterTitle()

    corr.Draw("COLZ")

    # For a small number of energy bins, use:
    # corr.Draw("COLZ TEXT")

    canvas.SaveAs(args.output)

    root_output = os.path.splitext(args.output)[0] + ".root"
    output_file = ROOT.TFile.Open(root_output, "RECREATE")
    corr.Write()
    output_file.Close()

    print(f"\nSaved plot:   {args.output}")
    print(f"Saved matrix: {root_output}")

    file_a.Close()
    file_b.Close()


if __name__ == "__main__":
    main()
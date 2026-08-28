#!/usr/bin/env python3

import argparse
import math
import ROOT
from array import array

# -----------------------------
# Constants
# -----------------------------

GF = 1.1663787e-5          # GeV^-2
ME = 0.00051099895        # electron mass [GeV]
GEV2_TO_CM2 = 0.389379e-27
SIN2THETAW_DEFAULT = 0.23129


def get_hist(root_file, hist_name):
    f = ROOT.TFile.Open(root_file)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open file: {root_file}")

    h = f.Get(hist_name)
    if not h:
        raise RuntimeError(f"Could not find histogram {hist_name} in {root_file}")

    # Detach from file so it survives after file close
    h.SetDirectory(0)
    f.Close()
    return h


def couplings(flavor, sin2thetaW):
    """
    Tree-level couplings for neutrino-electron elastic scattering.
    """

    if flavor == "nue" or flavor == "nuebar":
        gL = 0.5 + sin2thetaW
        gR = sin2thetaW
    elif flavor == "numu" or flavor == "numubar":
        gL = -0.5 + sin2thetaW
        gR = sin2thetaW
    else:
        raise ValueError(f"Unknown flavor: {flavor}")

    return gL, gR


def te_max(Enu):
    """
    Maximum electron kinetic energy for a neutrino with energy Enu.

    Enu, Te in GeV.
    """

    return 2.0 * Enu * Enu / (ME + 2.0 * Enu)


def dsigma_dTe(Enu, Te, flavor, sin2thetaW):
    """
    d sigma / d Te for nu-e elastic scattering.

    Input:
      Enu [GeV]
      Te  [GeV]

    Output:
      cm^2 / GeV
    """

    if Enu <= 0.0:
        return 0.0

    if Te < 0.0 or Te > te_max(Enu):
        return 0.0

    gL, gR = couplings(flavor, sin2thetaW)

    # For antineutrinos, swap gL and gR
    if flavor.endswith("bar"):
        gL, gR = gR, gL

    y = Te / Enu

    bracket = (
        gL * gL
        + gR * gR * (1.0 - y) * (1.0 - y)
        - gL * gR * ME * Te / (Enu * Enu)
    )

    prefactor = 2.0 * GF * GF * ME / math.pi

    return prefactor * bracket * GEV2_TO_CM2


def make_output_hist_from_edges(name, title, edges):
    h = ROOT.TH1D(name, title, len(edges) - 1, array("d", edges))
    h.GetXaxis().SetTitle("True lepton energy E_{l}^{true} [GeV]")
    h.GetYaxis().SetTitle("Predicted #nu-e events")
    h.SetDirectory(0)
    return h


def fold_one_flux(
    h_flux,
    h_out,
    flavor,
    target_electrons,
    pot,
    sin2thetaW,
    flux_bin_integrated,
    flux_area_unit,
):
    """
    Fold one flux histogram into electron kinetic-energy distribution.

    Expected flux units:
      Usually flux is something like nu / m^2 / POT / GeV,
      or nu / cm^2 / POT / GeV.

    Cross section here is in cm^2.
    Therefore:
      if flux_area_unit == "m2", convert flux from /m^2 to /cm^2 by dividing by 1e4.
      if flux_area_unit == "cm2", leave as-is.
    """

    if flux_area_unit == "m2":
        area_factor = 1.0 / 1.0e4
    elif flux_area_unit == "cm2":
        area_factor = 1.0
    else:
        raise ValueError("flux_area_unit must be either 'm2' or 'cm2'")

    n_Enu_bins = h_flux.GetNbinsX()
    n_Te_bins = h_out.GetNbinsX()

    total = 0.0

    for j in range(1, n_Enu_bins + 1):
        Enu_low = h_flux.GetBinLowEdge(j)
        dEnu = h_flux.GetBinWidth(j)
        Enu = Enu_low + 0.5 * dEnu

        flux = h_flux.GetBinContent(j)

        if Enu <= 0.0 or dEnu <= 0.0 or flux == 0.0:
            continue

        # If flux is differential per GeV, multiply by dEnu.
        # If histogram content is already bin-integrated, do not.
        if flux_bin_integrated:
            flux_bin = flux
        else:
            flux_bin = flux * dEnu

        # Apply POT and area conversion.
        flux_bin *= pot
        flux_bin *= area_factor

        for i in range(1, n_Te_bins + 1):
            El_low = h_out.GetBinLowEdge(i)
            El_high = El_low + h_out.GetBinWidth(i)

            # The paper bins are wide, especially the last bin,
            # so internally subdivide the lepton-energy bin.
            n_sub = max(20, int((El_high - El_low) / 0.05))
            dEl = (El_high - El_low) / n_sub

            dN_bin = 0.0

            for k in range(n_sub):
                El = El_low + (k + 0.5) * dEl

                # Cross section is written in terms of electron kinetic energy.
                Te = El - ME

                if Te <= 0.0:
                    continue

                if Te > te_max(Enu):
                    continue

                dsig = dsigma_dTe(
                    Enu=Enu,
                    Te=Te,
                    flavor=flavor,
                    sin2thetaW=sin2thetaW,
                )

                # dTe = dEl since Te = El - me.
                dN_bin += flux_bin * dsig * dEl * target_electrons

            old = h_out.GetBinContent(i)
            h_out.SetBinContent(i, old + dN_bin)

            total += dN_bin

    return total


def main():
    parser = argparse.ArgumentParser(
        description="Fold MINERvA/DUNE-style flux ROOT histograms into a nu-e elastic prediction."
    )

    parser.add_argument("--nue", default=None, help="ROOT file containing nue flux.")
    parser.add_argument("--nuebar", default=None, help="ROOT file containing nuebar flux.")
    parser.add_argument("--numu", default=None, help="ROOT file containing numu flux.")
    parser.add_argument("--numubar", default=None, help="ROOT file containing numubar flux.")

    parser.add_argument(
        "--hist-name",
        default="flux_E_cvweighted",
        help="Flux histogram name inside each ROOT file.",
    )

    parser.add_argument(
        "--output",
        default="nue_elastic_prediction.root",
        help="Output ROOT file.",
    )

    parser.add_argument("--Te-min", type=float, default=0.0)
    parser.add_argument("--Te-max", type=float, default=20.0)
    parser.add_argument("--n-Te-bins", type=int, default=200)

    parser.add_argument(
        "--target-electrons",
        type=float,
        default=1.0,
        help="Number of target electrons. Use 1 if your flux/event-rate input already includes this normalization.",
    )

    parser.add_argument(
        "--pot",
        type=float,
        default=1.0,
        help="POT normalization. Use 1 if the flux is already POT-scaled.",
    )

    parser.add_argument(
        "--sin2thetaW",
        type=float,
        default=SIN2THETAW_DEFAULT,
        help="sin^2(theta_W).",
    )

    parser.add_argument(
        "--flux-bin-integrated",
        action="store_true",
        help="Use if flux histogram contents are already integrated over Enu bin width.",
    )

    parser.add_argument(
        "--flux-area-unit",
        choices=["cm2", "m2"],
        default="m2",
        help="Area unit in the flux denominator. Cross section is cm^2. MINERvA fluxes are often per m^2, so default is m2.",
    )

    args = parser.parse_args()

    PAPER_LEPTON_ENERGY_EDGES = [0.8, 2.0, 3.0, 5.0, 7.0, 9.0, 120.0]

    input_files = {
        "nue": args.nue,
        "nuebar": args.nuebar,
        "numu": args.numu,
        "numubar": args.numubar,
    }

    h_total = make_output_hist_from_edges(
        "h_nue_elastic_total",
        "Predicted #nu-e elastic distribution;E_{l}^{true} [GeV];Events",
        PAPER_LEPTON_ENERGY_EDGES,
    )

    h_by_flavor = {}
    totals = {}

    for flavor, filename in input_files.items():
        if filename is None:
            continue

        print(f"Reading {flavor} flux from {filename}")

        h_flux = get_hist(filename, args.hist_name)

        h_pred = make_output_hist_from_edges(
            f"h_nue_elastic_{flavor}",
            f"Predicted #nu-e elastic from {flavor};E_{{l}}^{{true}} [GeV];Events",
            PAPER_LEPTON_ENERGY_EDGES,
        )

        total_flavor = fold_one_flux(
            h_flux=h_flux,
            h_out=h_pred,
            flavor=flavor,
            target_electrons=args.target_electrons,
            pot=args.pot,
            sin2thetaW=args.sin2thetaW,
            flux_bin_integrated=args.flux_bin_integrated,
            flux_area_unit=args.flux_area_unit,
        )

        h_total.Add(h_pred)
        h_by_flavor[flavor] = h_pred
        totals[flavor] = total_flavor

    fout = ROOT.TFile(args.output, "RECREATE")
    h_total.Write()

    for flavor, h in h_by_flavor.items():
        h.Write()

    fout.Close()

    print()
    print(f"Wrote output to {args.output}")
    print()
    print("Predicted totals:")
    for flavor in ["nue", "nuebar", "numu", "numubar"]:
        if flavor in totals:
            print(f"  {flavor:7s}: {totals[flavor]:.8e}")

    print(f"  {'total':7s}: {h_total.Integral():.8e}")


if __name__ == "__main__":
    main()
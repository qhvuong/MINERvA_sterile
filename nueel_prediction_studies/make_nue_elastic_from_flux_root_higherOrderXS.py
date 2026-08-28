#!/usr/bin/env python3

import argparse
import math
import ROOT
from array import array

# -----------------------------
# Constants
# -----------------------------

GF = 1.1663787e-5          # Fermi constant [GeV^-2]
ME = 0.00051099895        # electron mass [GeV]
ALPHA_EM = 1.0 / 137.035999084
GEV2_TO_CM2 = 0.389379e-27

# One-loop chiral couplings used in the MINERvA
# neutrino-electron scattering treatment (Park et al.).
#
# C_LL depends on the incoming neutrino flavor.
# C_LR is common to nu_e and nu_mu.
CLL_NUE = 0.7276
CLL_NUMU = -0.2730
CLR = 0.2334


def get_hist(root_file, hist_name):
    f = ROOT.TFile.Open(root_file)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open file: {root_file}")

    h = f.Get(hist_name)
    if not h:
        raise RuntimeError(f"Could not find histogram {hist_name} in {root_file}")

    h.SetDirectory(0)
    f.Close()
    return h


def dilog(x):
    """
    Real dilogarithm Li_2(x) for 0 <= x <= 1.

    Li_2(x) = sum_{k=1}^\infty x^k/k^2.

    For x > 1/2, use
      Li_2(x) = pi^2/6 - ln(x) ln(1-x) - Li_2(1-x)
    to improve convergence.
    """
    if x < 0.0 or x > 1.0:
        raise ValueError(f"dilog implementation requires 0 <= x <= 1; got {x}")

    if x == 0.0:
        return 0.0

    if x == 1.0:
        return math.pi * math.pi / 6.0

    if x > 0.5:
        return (
            math.pi * math.pi / 6.0
            - math.log(x) * math.log1p(-x)
            - dilog(1.0 - x)
        )

    total = 0.0
    term = x
    k = 1

    while True:
        add = term / (k * k)
        total += add

        if abs(add) < 1.0e-15:
            break

        k += 1
        term *= x

        if k > 100000:
            raise RuntimeError("dilog series failed to converge")

    return total


def get_one_loop_couplings(flavor):
    """
    One-loop dimensionless chiral couplings C_LL and C_LR.

    Values:
      nu_e, nubar_e   : C_LL =  0.7276
      nu_mu, nubar_mu : C_LL = -0.2730
      all             : C_LR =  0.2334

    Neutrino versus antineutrino is handled in dsigma_dTe_one_loop()
    by interchanging the leading left/right pieces, not by changing
    these coupling values.
    """
    if flavor in ("nue", "nuebar"):
        cLL = CLL_NUE
    elif flavor in ("numu", "numubar"):
        cLL = CLL_NUMU
    else:
        raise ValueError(f"Unknown flavor: {flavor}")

    return cLL, CLR


def te_max(Enu):
    """
    Maximum electron kinetic energy for elastic nu-e scattering.

    Enu, Te in GeV.
    """
    return 2.0 * Enu * Enu / (ME + 2.0 * Enu)


def x1_radiative(Enu, y):
    """
    O(alpha) radiative correction X1 for the electron-energy spectrum.
    """
    log_y = math.log(y)
    log_1my = math.log1p(-y)
    log_scale = math.log(2.0 * Enu / ME)
    log_ratio = math.log(1.0 / y - 1.0)
    li2_y = dilog(y)

    return (
        (1.0 / 12.0)
        * (6.0 * y + 12.0 * log_1my - 6.0 * log_y - 5.0)
        * log_scale
        - 0.5 * li2_y
        + y * y / 24.0
        - 11.0 * y / 12.0
        - 0.5 * log_ratio * log_ratio
        + y * log_y
        - (1.0 / 12.0) * (6.0 * y + 23.0) * log_1my
        + math.pi * math.pi / 12.0
        - 47.0 / 36.0
    )


def x2_radiative(Enu, y):
    """
    O(alpha) radiative correction X2 for the electron-energy spectrum.
    """
    one_minus_y = 1.0 - y
    omy2 = one_minus_y * one_minus_y

    log_y = math.log(y)
    log_1my = math.log1p(-y)
    log_scale = math.log(2.0 * Enu / ME)
    li2_y = dilog(y)

    term1 = (
        -4.0 * y * y
        + (-6.0 * y * y + 6.0 * y - 3.0) * log_y
        + 11.0 * y
        + 6.0 * omy2 * log_1my
        - 7.0
    ) * log_scale / (6.0 * omy2)

    term2 = (
        (-y * y + y - 0.5)
        * (li2_y + log_y * log_y - math.pi * math.pi / 6.0)
        / omy2
    )

    term3 = (
        (4.0 * y * y + 2.0 * y - 3.0)
        * log_y
        / (4.0 * omy2)
    )

    term4 = -(31.0 - 49.0 * y) / (72.0 * one_minus_y)

    term5 = (
        (10.0 * y - 7.0)
        * log_1my
        / (6.0 * one_minus_y)
    )

    term6 = log_1my * (log_y - 0.5 * log_1my)

    return term1 + term2 + term3 + term4 + term5 + term6


def x3_radiative(Enu, y):
    """
    O(alpha) radiative correction X3 for the electron-energy spectrum.

    Written directly in terms of y = Te / Enu.
    """
    Te = y * Enu
    root = math.sqrt(Te * (2.0 * ME + Te))
    Ee = ME + Te

    # Physical Te > 0 guarantees root > 0.
    log_arg_inner = (root + Ee) / ME
    inner = Ee * math.log(log_arg_inner) / root - 1.0

    log_arg_outer = -ME / (root + Ee) + 1.0 - y

    if log_arg_outer <= 0.0:
        raise ValueError(
            "Non-positive logarithm argument in X3: "
            f"Enu={Enu}, y={y}, arg={log_arg_outer}"
        )

    return math.log(log_arg_outer) * inner


def dsigma_dTe_one_loop(Enu, Te, flavor):
    """
    One-loop electroweak + O(alpha) radiatively corrected
    neutrino-electron elastic differential cross section.

    Inputs:
      Enu [GeV]
      Te  [GeV]

    Output:
      d sigma / d Te [cm^2 / GeV]

    The calculation is performed using y = Te / Enu and
    d sigma / dTe = (1 / Enu) d sigma / dy.
    """
    if Enu <= 0.0:
        return 0.0

    if Te <= 0.0 or Te > te_max(Enu):
        return 0.0

    y = Te / Enu

    # The analytic O(alpha) expressions contain endpoint logarithms.
    # With midpoint integration we should never evaluate exactly at
    # y = 0 or y = 1, but enforce that here explicitly.
    if not (0.0 < y < 1.0):
        return 0.0

    cLL, cLR = get_one_loop_couplings(flavor)

    X1 = x1_radiative(Enu, y)
    X2 = x2_radiative(Enu, y)
    X3 = x3_radiative(Enu, y)

    alpha_over_pi = ALPHA_EM / math.pi

    corr1 = 1.0 + alpha_over_pi * X1
    corr2 = 1.0 + alpha_over_pi * X2
    corr3 = 1.0 + alpha_over_pi * X3

    # Exact Mandelstam s for an electron initially at rest.
    s = ME * ME + 2.0 * ME * Enu

    if flavor in ("nue", "numu"):
        bracket = (
            cLL * cLL * corr1
            + cLR * cLR * (1.0 - y) * (1.0 - y) * corr2
            - cLL * cLR * ME * y / Enu * corr3
        )
    elif flavor in ("nuebar", "numubar"):
        bracket = (
            cLR * cLR * corr1
            + cLL * cLL * (1.0 - y) * (1.0 - y) * corr2
            - cLL * cLR * ME * y / Enu * corr3
        )
    else:
        raise ValueError(f"Unknown flavor: {flavor}")

    dsigma_dy = (GF * GF * s / math.pi) * bracket

    # y = Te / Enu  ->  dy/dTe = 1/Enu
    dsigma_dTe = dsigma_dy / Enu

    return dsigma_dTe * GEV2_TO_CM2


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
    flux_bin_integrated,
    flux_area_unit,
):
    """
    Fold one MINERvA flux histogram into an electron-energy distribution.

    Standard MINERvA flux input:
      nu / m^2 / POT / GeV

    The cross section here is in cm^2, so for flux per m^2
    we multiply by 1e-4 to convert m^-2 -> cm^-2.
    """
    if flux_area_unit == "m2":
        area_factor = 1.0 / 1.0e4
    elif flux_area_unit == "cm2":
        area_factor = 1.0
    else:
        raise ValueError("flux_area_unit must be either 'm2' or 'cm2'")

    n_Enu_bins = h_flux.GetNbinsX()
    n_El_bins = h_out.GetNbinsX()

    total = 0.0

    for j in range(1, n_Enu_bins + 1):
        Enu_low = h_flux.GetBinLowEdge(j)
        dEnu = h_flux.GetBinWidth(j)
        Enu = Enu_low + 0.5 * dEnu

        flux = h_flux.GetBinContent(j)

        if Enu <= 0.0 or dEnu <= 0.0 or flux == 0.0:
            continue

        if flux_bin_integrated:
            flux_bin = flux
        else:
            flux_bin = flux * dEnu

        flux_bin *= pot
        flux_bin *= area_factor

        for i in range(1, n_El_bins + 1):
            El_low = h_out.GetBinLowEdge(i)
            El_high = El_low + h_out.GetBinWidth(i)

            # Wide output bins are numerically integrated with
            # sub-bins in total electron energy.
            n_sub = max(20, int((El_high - El_low) / 0.05))
            dEl = (El_high - El_low) / n_sub

            dN_bin = 0.0

            for k in range(n_sub):
                El = El_low + (k + 0.5) * dEl

                # Te = Ee - me; dTe = dEe.
                Te = El - ME

                if Te <= 0.0 or Te > te_max(Enu):
                    continue

                dsig = dsigma_dTe_one_loop(
                    Enu=Enu,
                    Te=Te,
                    flavor=flavor,
                )

                dN_bin += flux_bin * dsig * dEl * target_electrons

            old = h_out.GetBinContent(i)
            h_out.SetBinContent(i, old + dN_bin)
            total += dN_bin

    return total


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Fold MINERvA flux ROOT histograms into a nu-e elastic "
            "prediction using a one-loop electroweak + O(alpha) "
            "radiative-correction approximation."
        )
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
        default="nue_elastic_prediction_higherOrderXS.root",
        help="Output ROOT file.",
    )

    parser.add_argument(
        "--target-electrons",
        type=float,
        default=1.0,
        help=(
            "Number of target electrons. "
            "Use the target-electron normalization corresponding "
            "to the desired MINERvA fiducial volume."
        ),
    )

    parser.add_argument(
        "--pot",
        type=float,
        default=1.0,
        help="POT normalization. Standard MINERvA flux is per POT.",
    )

    parser.add_argument(
        "--flux-bin-integrated",
        action="store_true",
        help=(
            "Use only if flux histogram contents are already integrated "
            "over Enu bin width. Do NOT use for standard flux_E_cvweighted."
        ),
    )

    parser.add_argument(
        "--flux-area-unit",
        choices=["cm2", "m2"],
        default="m2",
        help=(
            "Area unit in the flux denominator. "
            "Standard MINERvA flux is per m^2."
        ),
    )

    args = parser.parse_args()

    PAPER_LEPTON_ENERGY_EDGES = [0.8, 2.0, 3.0, 5.0, 7.0, 9.0, 100.0]

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
            (
                f"Predicted #nu-e elastic from {flavor};"
                "E_{l}^{true} [GeV];Events"
            ),
            PAPER_LEPTON_ENERGY_EDGES,
        )

        total_flavor = fold_one_flux(
            h_flux=h_flux,
            h_out=h_pred,
            flavor=flavor,
            target_electrons=args.target_electrons,
            pot=args.pot,
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
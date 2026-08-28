#!/usr/bin/env python3

import argparse
import math
import ROOT
from array import array

import PlotUtils.LoadPlotUtilsLib

ME = 0.00051099895
GEV2_TO_CM2 = 0.389379e-27

# Tomalak effective couplings used in Maria's DSigmaDY_Tomalak().
# These are dimensionful (GeV^-2), so the normalization below is
# m_e E_nu / (4 pi), exactly as in Maria's implementation.
CLL_TOMALAK_NUE = 2.39818e-5
CLL_TOMALAK_NUMU = -0.90084e-5
CLR_TOMALAK = 0.76911e-5

PAPER_LEPTON_ENERGY_EDGES = [0.8, 2.0, 3.0, 5.0, 7.0, 9.0, 100.0]


def get_mnv_hist(root_file, hist_name):
    """
    Open a ROOT file and return both the file handle and the MnvH1D.

    Keep the ROOT file open for as long as the MnvH1D and its
    error-band universe histograms are being accessed.
    """

    f = ROOT.TFile.Open(root_file)

    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open file: {root_file}")

    h = f.Get(hist_name)

    if not h:
        f.Close()
        raise RuntimeError(
            f"Could not find histogram {hist_name} in {root_file}"
        )

    # Do NOT Clone(), SetDirectory(0), or close the file here.
    ROOT.SetOwnership(h, False)

    return f, h


def make_plain_hist(name, title, edges):
    h = ROOT.TH1D(name, title, len(edges) - 1, array("d", edges))
    h.SetDirectory(0)
    return h


def make_output_mnvh1d(name, title, edges, flux_band_name, n_universes):
    h = ROOT.PlotUtils.MnvH1D(
        name,
        title,
        len(edges) - 1,
        array("d", edges),
    )
    h.SetDirectory(0)

    if n_universes > 0:
        h.AddVertErrorBand(flux_band_name, n_universes)

    h.GetXaxis().SetTitle("True lepton energy E_{l}^{true} [GeV]")
    h.GetYaxis().SetTitle("Predicted #nu-e events")
    return h


def copy_plain_to_hist(h_plain, h_dest):
    if h_plain.GetNbinsX() != h_dest.GetNbinsX():
        raise RuntimeError("Histogram bin-count mismatch")

    for i in range(0, h_dest.GetNbinsX() + 2):
        h_dest.SetBinContent(i, h_plain.GetBinContent(i))
        h_dest.SetBinError(i, h_plain.GetBinError(i))


def _read_two_column_csv(path):
    xs = []
    ys = []

    with open(path, "r", encoding="utf-8-sig") as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = [x.strip() for x in line.split(",")]

            if len(fields) < 2:
                continue

            if fields[0] == "" or fields[1] == "":
                continue

            try:
                x = float(fields[0])
                y = float(fields[1])
            except ValueError:
                continue

            xs.append(x)
            ys.append(y)

    if len(xs) < 7:
        raise RuntimeError(
            f"Need at least 7 points for a pol6 fit; got {len(xs)} from {path}"
        )

    return xs, ys


def _fit_pol6_csv(path, name):
    """
    Reproduce Maria's TGraph(...)->Fit("pol6") parameterization.
    The Q option only suppresses ROOT fit output; it does not change the fit.
    """
    xs, ys = _read_two_column_csv(path)
    xarr = array("d", xs)
    yarr = array("d", ys)

    graph = ROOT.TGraph(len(xs), xarr, yarr)
    func = ROOT.TF1(name, "pol6", min(xs), max(xs))
    status = graph.Fit(func, "Q")

    if int(status) != 0:
        raise RuntimeError(
            f"pol6 fit failed for {path} with ROOT fit status {int(status)}"
        )

    return [func.GetParameter(i) for i in range(7)]


class RadCorrTomalak:
    """Python port of Maria's RadCorrectionTomalak.cxx."""

    def __init__(self, f1_csv, f10_csv, label):
        self.label = label

        self.a = _fit_pol6_csv(
            f1_csv,
            f"tomalak_pol6_{label}_f1",
        )
        self.b = _fit_pol6_csv(
            f10_csv,
            f"tomalak_pol6_{label}_f10",
        )

        print(f"Tomalak pol6 coefficients for {label}:")
        print(
            "  E_nu = 1 GeV :",
            " ".join(f"{x:.12e}" for x in self.a),
        )
        print(
            "  E_nu = 10 GeV:",
            " ".join(f"{x:.12e}" for x in self.b),
        )

    def get(self, Enu, y):
        if Enu <= 0.0:
            return 0.0

        log_interp = math.log(Enu) / math.log(10.0)

        rad_percent = 0.0
        y_power = 1.0
        for ak, bk in zip(self.a, self.b):
            coeff = ak + (bk - ak) * log_interp
            rad_percent += coeff * y_power
            y_power *= y

        return rad_percent / 100.0


def get_tomalak_couplings(flavor):
    if flavor in ("nue", "nuebar"):
        cLL = CLL_TOMALAK_NUE
    elif flavor in ("numu", "numubar"):
        cLL = CLL_TOMALAK_NUMU
    else:
        raise ValueError(f"Unknown flavor: {flavor}")
    return cLL, CLR_TOMALAK


def te_max(Enu):
    return 2.0 * Enu * Enu / (ME + 2.0 * Enu)


def dsigma_dTe_tomalak(Enu, Te, flavor, rad_corr):
    """Maria's DSigmaDY_Tomalak() converted to d sigma / d T_e."""
    if Enu <= 0.0:
        return 0.0
    if Te <= 0.0 or Te > te_max(Enu):
        return 0.0

    y = Te / Enu
    if not (0.0 < y < 1.0):
        return 0.0

    cLL, cLR = get_tomalak_couplings(flavor)
    factor = ME * Enu / (4.0 * math.pi)

    if flavor in ("nue", "numu"):
        bracket = (
            cLL * cLL
            + cLR * cLR * (1.0 - y) * (1.0 - y)
            - cLL * cLR * ME * y / Enu
        )
    elif flavor in ("nuebar", "numubar"):
        bracket = (
            cLR * cLR
            + cLL * cLL * (1.0 - y) * (1.0 - y)
            - cLR * cLL * ME * y / Enu
        )
    else:
        raise ValueError(f"Unknown flavor: {flavor}")

    dsigma_dy = factor * bracket
    delta = rad_corr.get(Enu, y)
    dsigma_dy *= (1.0 + delta)

    return (dsigma_dy / Enu) * GEV2_TO_CM2


def fold_one_flux(
    h_flux,
    h_out,
    flavor,
    target_electrons,
    pot,
    flux_bin_integrated,
    flux_area_unit,
    rad_corr,
):
    if flux_area_unit == "m2":
        area_factor = 1.0 / 1.0e4
    elif flux_area_unit == "cm2":
        area_factor = 1.0
    else:
        raise ValueError("flux_area_unit must be 'm2' or 'cm2'")

    total = 0.0

    for j in range(1, h_flux.GetNbinsX() + 1):
        Enu_low = h_flux.GetBinLowEdge(j)
        dEnu = h_flux.GetBinWidth(j)
        Enu = Enu_low + 0.5 * dEnu
        flux = h_flux.GetBinContent(j)

        if Enu <= 0.0 or dEnu <= 0.0 or flux == 0.0:
            continue

        flux_bin = flux if flux_bin_integrated else flux * dEnu
        flux_bin *= pot
        flux_bin *= area_factor

        for i in range(1, h_out.GetNbinsX() + 1):
            El_low = h_out.GetBinLowEdge(i)
            El_high = El_low + h_out.GetBinWidth(i)

            n_sub = max(20, int((El_high - El_low) / 0.05))
            dEl = (El_high - El_low) / n_sub
            dN_bin = 0.0

            for k in range(n_sub):
                El = El_low + (k + 0.5) * dEl
                Te = El - ME

                if Te <= 0.0 or Te > te_max(Enu):
                    continue

                dsig = dsigma_dTe_tomalak(
                    Enu=Enu,
                    Te=Te,
                    flavor=flavor,
                    rad_corr=rad_corr,
                )

                dN_bin += flux_bin * dsig * dEl * target_electrons

            h_out.SetBinContent(i, h_out.GetBinContent(i) + dN_bin)
            total += dN_bin

    return total


def fold_cv_and_flux_universes(
    h_flux_mnv,
    flavor,
    target_electrons,
    pot,
    flux_bin_integrated,
    flux_area_unit,
    flux_band_name,
    rad_corr,
):
    if not h_flux_mnv.HasVertErrorBand(flux_band_name):
        available = [str(x) for x in h_flux_mnv.GetVertErrorBandNames()]
        raise RuntimeError(
            f"Input flux histogram has no vertical error band "
            f"'{flux_band_name}'. Available bands: {available}"
        )

    band = h_flux_mnv.GetVertErrorBand(flux_band_name)
    n_universes = band.GetNHists()

    h_cv = make_plain_hist(
        f"h_tmp_{flavor}_cv",
        "",
        PAPER_LEPTON_ENERGY_EDGES,
    )

    fold_one_flux(
        h_flux=h_flux_mnv,
        h_out=h_cv,
        flavor=flavor,
        target_electrons=target_electrons,
        pot=pot,
        flux_bin_integrated=flux_bin_integrated,
        flux_area_unit=flux_area_unit,
        rad_corr=rad_corr,
    )

    universe_hists = []

    for u in range(n_universes):
        h_flux_u = band.GetHist(u)

        h_u = make_plain_hist(
            f"h_tmp_{flavor}_flux_u{u}",
            "",
            PAPER_LEPTON_ENERGY_EDGES,
        )

        fold_one_flux(
            h_flux=h_flux_u,
            h_out=h_u,
            flavor=flavor,
            target_electrons=target_electrons,
            pot=pot,
            flux_bin_integrated=flux_bin_integrated,
            flux_area_unit=flux_area_unit,
            rad_corr=rad_corr,
        )

        universe_hists.append(h_u)

    return h_cv, universe_hists


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build a Tomalak-radiatively-corrected nu-e elastic MnvH1D "
            "prediction from MINERvA fluxes while preserving Flux universes."
        )
    )

    parser.add_argument("--nue", default=None)
    parser.add_argument("--nuebar", default=None)
    parser.add_argument("--numu", default=None)
    parser.add_argument("--numubar", default=None)

    parser.add_argument(
        "--hist-name",
        default="flux_E_cvweighted",
        help="Input MnvH1D flux histogram name.",
    )

    parser.add_argument(
        "--flux-band-name",
        default="Flux",
        help="Vertical flux error band to propagate.",
    )

    parser.add_argument(
        "--output",
        default="nue_elastic_prediction_TomalakXS_mnv.root",
    )


    parser.add_argument(
        "--tomalak-dir",
        default=(
            "/exp/minerva/data/users/qvuong/"
            "nueel_prediction_studies/tomalak_inputs"
        ),
        help="Directory containing the frozen Tomalak correction CSV files.",
    )

    parser.add_argument(
        "--target-electrons",
        type=float,
        default=1.98e30,
    )

    parser.add_argument(
        "--pot",
        type=float,
        default=3.331982991676675e20,
    )

    parser.add_argument(
        "--flux-bin-integrated",
        action="store_true",
        help=(
            "Do NOT use for standard MINERvA flux_E_cvweighted, "
            "which is differential per GeV."
        ),
    )

    parser.add_argument(
        "--flux-area-unit",
        choices=["cm2", "m2"],
        default="m2",
    )

    args = parser.parse_args()

    tomalak_dir = args.tomalak_dir.rstrip("/")

    rad_corr = {
        "numu": RadCorrTomalak(
            f"{tomalak_dir}/f1_muonNeutrinos.csv",
            f"{tomalak_dir}/f10_muonNeutrinos.csv",
            "numu",
        ),

        "numubar": RadCorrTomalak(
            f"{tomalak_dir}/f1_muonAntineutrinos.csv",
            f"{tomalak_dir}/f10_muonAntineutrinos.csv",
            "numubar",
        ),

        "nue": RadCorrTomalak(
            f"{tomalak_dir}/f1_electronNeutrinos.csv",
            f"{tomalak_dir}/f10_electronNeutrinos.csv",
            "nue",
        ),

        "nuebar": RadCorrTomalak(
            f"{tomalak_dir}/f1_electronAntineutrinos.csv",
            f"{tomalak_dir}/f10_electronAntineutrinos.csv",
            "nuebar",
        ),
    }

    input_files = {
        "nue": args.nue,
        "nuebar": args.nuebar,
        "numu": args.numu,
        "numubar": args.numubar,
    }

    active_flavors = [
        flavor for flavor, filename in input_files.items()
        if filename is not None
    ]

    if not active_flavors:
        raise RuntimeError("No input flux files were supplied")

    flavor_cv = {}
    flavor_universes = {}
    n_universes = None

    for flavor in active_flavors:
        filename = input_files[flavor]
        print(f"Reading {flavor} flux from {filename}")

        f_flux, h_flux = get_mnv_hist(
            filename,
            args.hist_name,
        )

        if not h_flux.HasVertErrorBand(args.flux_band_name):
            available = [str(x) for x in h_flux.GetVertErrorBandNames()]
            raise RuntimeError(
                f"{filename}: histogram {args.hist_name} has no "
                f"'{args.flux_band_name}' vertical error band. "
                f"Available bands: {available}"
            )

        this_n_universes = (
            h_flux.GetVertErrorBand(args.flux_band_name).GetNHists()
        )

        if n_universes is None:
            n_universes = this_n_universes
            print(
                f"Using {n_universes} universes from "
                f"'{args.flux_band_name}'"
            )
        elif this_n_universes != n_universes:
            raise RuntimeError(
                "Flux universe count mismatch: "
                f"{flavor} has {this_n_universes}, "
                f"expected {n_universes}"
            )

        h_cv, h_univs = fold_cv_and_flux_universes(
            h_flux_mnv=h_flux,
            flavor=flavor,
            target_electrons=args.target_electrons,
            pot=args.pot,
            flux_bin_integrated=args.flux_bin_integrated,
            flux_area_unit=args.flux_area_unit,
            flux_band_name=args.flux_band_name,
            rad_corr=rad_corr[flavor],
        )

        flavor_cv[flavor] = h_cv
        flavor_universes[flavor] = h_univs

        print(
            f"  {flavor:7s} CV total = "
            f"{h_cv.Integral():.8e}"
        )

        # Now all information from this input MnvH1D has been copied
        # into our own temporary TH1Ds, so it is safe to close.
        f_flux.Close()

    h_by_flavor = {}

    for flavor in active_flavors:
        h_out = make_output_mnvh1d(
            f"h_nue_elastic_{flavor}",
            f"Predicted #nu-e elastic from {flavor};E_{{l}}^{{true}} [GeV];Events",
            PAPER_LEPTON_ENERGY_EDGES,
            args.flux_band_name,
            n_universes,
        )

        copy_plain_to_hist(flavor_cv[flavor], h_out)

        out_band = h_out.GetVertErrorBand(args.flux_band_name)
        for u in range(n_universes):
            copy_plain_to_hist(
                flavor_universes[flavor][u],
                out_band.GetHist(u),
            )

        h_by_flavor[flavor] = h_out

    h_total = make_output_mnvh1d(
        "h_nue_elastic_total",
        "Predicted #nu-e elastic distribution;E_{l}^{true} [GeV];Events",
        PAPER_LEPTON_ENERGY_EDGES,
        args.flux_band_name,
        n_universes,
    )

    for flavor in active_flavors:
        for i in range(0, h_total.GetNbinsX() + 2):
            h_total.SetBinContent(
                i,
                h_total.GetBinContent(i)
                + flavor_cv[flavor].GetBinContent(i),
            )

    total_band = h_total.GetVertErrorBand(args.flux_band_name)

    for u in range(n_universes):
        h_total_u = total_band.GetHist(u)

        for flavor in active_flavors:
            h_flavor_u = flavor_universes[flavor][u]

            for i in range(0, h_total.GetNbinsX() + 2):
                h_total_u.SetBinContent(
                    i,
                    h_total_u.GetBinContent(i)
                    + h_flavor_u.GetBinContent(i),
                )

    fout = ROOT.TFile(args.output, "RECREATE")
    h_total.Write()

    for flavor in active_flavors:
        h_by_flavor[flavor].Write()

    fout.Write()
    fout.Close()

    print()
    print(f"Wrote output to {args.output}")
    print()
    print("Predicted CV totals:")

    for flavor in ["nue", "nuebar", "numu", "numubar"]:
        if flavor in flavor_cv:
            print(f"  {flavor:7s}: {flavor_cv[flavor].Integral():.8e}")

    print(f"  {'total':7s}: {h_total.Integral():.8e}")
    print()
    print(
        f"Preserved vertical error band '{args.flux_band_name}' "
        f"with {n_universes} universes."
    )
    print(
        "The total uses the same universe index u across all flavors."
    )


if __name__ == "__main__":
    main()
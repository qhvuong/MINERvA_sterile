#!/usr/bin/env python

import argparse
import os

import numpy as np
import ROOT


def make_hist(name, title, values, nbins=80, xmin=None, xmax=None):
    values = np.asarray(values)
    values = values[np.isfinite(values)]

    if len(values) == 0:
        return None

    if xmin is None or xmax is None:
        vmin = np.min(values)
        vmax = np.max(values)

        # Give the distribution a little padding.
        width = vmax - vmin

        if width <= 0.0:
            width = 0.01

        xmin = vmin - 0.10 * width
        xmax = vmax + 0.10 * width

    h = ROOT.TH1D(
        name,
        title,
        nbins,
        xmin,
        xmax
    )

    h.SetDirectory(0)

    for value in values:
        h.Fill(value)

    return h


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Plot true-E BeamFocus universe distributions for "
            "numu, nue, and numu/nue using the NPZ output from "
            "study_beamfocus_trueE_cancellation.py."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input NPZ file"
    )

    parser.add_argument(
        "--output_dir",
        default="beamfocus_trueE_plots",
        help="Directory for plots"
    )

    parser.add_argument(
        "--nbins",
        type=int,
        default=80,
        help="Number of histogram bins for universe distributions"
    )

    parser.add_argument(
        "--emin",
        type=float,
        default=None,
        help="Only plot true-E bins starting at or above this value"
    )

    parser.add_argument(
        "--emax",
        type=float,
        default=None,
        help="Only plot true-E bins ending at or below this value"
    )

    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    os.makedirs(args.output_dir, exist_ok=True)

    # ------------------------------------------------------------
    # Load NPZ
    # ------------------------------------------------------------

    data = np.load(args.input)

    bin_edges = data["bin_edges"]
    weights_numu = data["weights_numu"]
    weights_nue = data["weights_nue"]

    n_bins = len(bin_edges) - 1

    print("Input:", args.input)
    print("Energy bins:", n_bins)
    print("Numu shape:", weights_numu.shape)
    print("Nue shape :", weights_nue.shape)

    # ------------------------------------------------------------
    # Output ROOT file containing the distributions too
    # ------------------------------------------------------------

    root_output = os.path.join(
        args.output_dir,
        "beamfocus_trueE_distributions.root"
    )

    fout = ROOT.TFile.Open(
        root_output,
        "RECREATE"
    )

    # ------------------------------------------------------------
    # Loop over true-E bins
    # ------------------------------------------------------------

    for ibin in range(n_bins):

        elo = bin_edges[ibin]
        ehi = bin_edges[ibin + 1]

        if args.emin is not None and elo < args.emin:
            continue

        if args.emax is not None and ehi > args.emax:
            continue

        mu = weights_numu[ibin, :]
        e = weights_nue[ibin, :]

        valid = (
            np.isfinite(mu)
            & np.isfinite(e)
            & (np.abs(e) > 1.0e-12)
        )

        mu = mu[valid]
        e = e[valid]

        if len(mu) == 0:
            print(
                "Skipping {:.1f}-{:.1f} GeV: no valid universes".format(
                    elo,
                    ehi
                )
            )
            continue

        # --------------------------------------------------------
        # Quantities we want to inspect
        # --------------------------------------------------------

        # Individual BeamFocus shifts
        dmu = mu - 1.0
        de = e - 1.0

        # Actual analysis-ratio direction
        ratio = mu / e
        dratio = ratio - 1.0

        sigma_mu = np.sqrt(
            np.mean(dmu ** 2)
        )

        sigma_e = np.sqrt(
            np.mean(de ** 2)
        )

        sigma_ratio = np.sqrt(
            np.mean(dratio ** 2)
        )

        correlation = np.corrcoef(
            dmu,
            de
        )[0, 1]

        print(
            "{:.1f}-{:.1f} GeV : "
            "sig_mu={:.3f}%  "
            "sig_e={:.3f}%  "
            "rho={:.4f}  "
            "sig_mu/e={:.3f}%".format(
                elo,
                ehi,
                100.0 * sigma_mu,
                100.0 * sigma_e,
                correlation,
                100.0 * sigma_ratio
            )
        )

        # --------------------------------------------------------
        # Use percent shifts in the plots
        # --------------------------------------------------------

        dmu_percent = 100.0 * dmu
        de_percent = 100.0 * de
        dratio_percent = 100.0 * dratio

        # Common range for numu/nue comparison.
        all_species = np.concatenate(
            [
                dmu_percent,
                de_percent
            ]
        )

        xmin_species = np.min(all_species)
        xmax_species = np.max(all_species)

        width_species = xmax_species - xmin_species

        if width_species <= 0.0:
            width_species = 1.0

        xmin_species -= 0.10 * width_species
        xmax_species += 0.10 * width_species

        # --------------------------------------------------------
        # Numu distribution
        # --------------------------------------------------------

        h_mu = make_hist(
            "beamfocus_numu_Ebin{:02d}".format(ibin),
            (
                "BeamFocus true #nu_{{#mu}}, "
                "{:.1f} < E_{{#nu}}^{{true}} < {:.1f} GeV;"
                "BeamFocus shift [%];Universes"
            ).format(elo, ehi),
            dmu_percent,
            args.nbins,
            xmin_species,
            xmax_species
        )

        # --------------------------------------------------------
        # Nue distribution
        # --------------------------------------------------------

        h_e = make_hist(
            "beamfocus_nue_Ebin{:02d}".format(ibin),
            (
                "BeamFocus true #nu_{{e}}, "
                "{:.1f} < E_{{#nu}}^{{true}} < {:.1f} GeV;"
                "BeamFocus shift [%];Universes"
            ).format(elo, ehi),
            de_percent,
            args.nbins,
            xmin_species,
            xmax_species
        )

        # --------------------------------------------------------
        # Numu / nue distribution
        # --------------------------------------------------------

        h_ratio = make_hist(
            "beamfocus_numu_over_nue_Ebin{:02d}".format(ibin),
            (
                "BeamFocus #nu_{{#mu}}/#nu_{{e}}, "
                "{:.1f} < E_{{#nu}}^{{true}} < {:.1f} GeV;"
                "W_{{#nu_{{#mu}}}}/W_{{#nu_{{e}}}} - 1 [%];"
                "Universes"
            ).format(elo, ehi),
            dratio_percent,
            args.nbins
        )

        # --------------------------------------------------------
        # Plot 1:
        # numu and nue universe distributions overlaid
        # --------------------------------------------------------

        c1 = ROOT.TCanvas(
            "c_species_{}".format(ibin),
            "",
            900,
            700
        )

        h_mu.SetLineColor(ROOT.kBlue + 1)
        h_mu.SetLineWidth(3)

        h_e.SetLineColor(ROOT.kRed + 1)
        h_e.SetLineWidth(3)

        ymax = max(
            h_mu.GetMaximum(),
            h_e.GetMaximum()
        )

        h_mu.SetMaximum(1.20 * ymax)

        h_mu.SetTitle(
            (
                "BeamFocus universes, "
                "{:.1f} < E_{{#nu}}^{{true}} < {:.1f} GeV;"
                "BeamFocus shift [%];Universes"
            ).format(elo, ehi)
        )

        h_mu.Draw("hist")
        h_e.Draw("hist same")

        line = ROOT.TLine(
            0.0,
            0.0,
            0.0,
            1.20 * ymax
        )

        line.SetLineStyle(2)
        line.Draw()

        legend = ROOT.TLegend(
            0.60,
            0.72,
            0.88,
            0.88
        )

        legend.SetBorderSize(0)
        legend.SetFillStyle(0)

        legend.AddEntry(
            h_mu,
            "#nu_{#mu}",
            "l"
        )

        legend.AddEntry(
            h_e,
            "#nu_{e}",
            "l"
        )

        legend.Draw()

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.035)

        latex.DrawLatex(
            0.17,
            0.84,
            "#sigma_{#mu} = %.2f%%" % (
                100.0 * sigma_mu
            )
        )

        latex.DrawLatex(
            0.17,
            0.79,
            "#sigma_{e} = %.2f%%" % (
                100.0 * sigma_e
            )
        )

        latex.DrawLatex(
            0.17,
            0.74,
            "#rho(#nu_{#mu},#nu_{e}) = %.3f" % (
                correlation
            )
        )

        species_png = os.path.join(
            args.output_dir,
            "beamfocus_species_E_{:.1f}_{:.1f}.png".format(
                elo,
                ehi
            )
        )

        c1.SaveAs(species_png)

        # --------------------------------------------------------
        # Plot 2:
        # cancellation in numu / nue
        # --------------------------------------------------------

        c2 = ROOT.TCanvas(
            "c_ratio_{}".format(ibin),
            "",
            900,
            700
        )

        h_ratio.SetLineColor(ROOT.kBlack)
        h_ratio.SetLineWidth(3)

        h_ratio.Draw("hist")

        ratio_ymax = h_ratio.GetMaximum()

        line_ratio = ROOT.TLine(
            0.0,
            0.0,
            0.0,
            1.10 * ratio_ymax
        )

        line_ratio.SetLineStyle(2)
        line_ratio.Draw()

        latex_ratio = ROOT.TLatex()
        latex_ratio.SetNDC()
        latex_ratio.SetTextSize(0.04)

        latex_ratio.DrawLatex(
            0.18,
            0.84,
            "#sigma(#nu_{#mu}/#nu_{e}) = %.2f%%" % (
                100.0 * sigma_ratio
            )
        )

        ratio_png = os.path.join(
            args.output_dir,
            "beamfocus_numu_over_nue_E_{:.1f}_{:.1f}.png".format(
                elo,
                ehi
            )
        )

        c2.SaveAs(ratio_png)

        # --------------------------------------------------------
        # Save histograms
        # --------------------------------------------------------

        fout.cd()

        h_mu.Write()
        h_e.Write()
        h_ratio.Write()

    fout.Close()

    print()
    print("Saved plots to:", args.output_dir)
    print("Saved distributions to:", root_output)


if __name__ == "__main__":
    main()
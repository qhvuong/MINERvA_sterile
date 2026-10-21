#!/usr/bin/env python3

import os
import ROOT
from PlotUtils import MnvPlotter

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# ============================================================
# INPUT
# ============================================================

flux_dir = "../custom_plotutils/data/flux"

# MINERvA LE playlist mapping
playlist_map = {
    "FHC": ["minerva1", "minerva13"],
    "RHC": ["minerva5"],
}

species = {
    "numu": {
        "pdg": 14,
        "label": "#nu_{#mu}",
        "color": ROOT.kBlack,
        "style": 1,
    },
    "numubar": {
        "pdg": -14,
        "label": "#bar{#nu}_{#mu}",
        "color": ROOT.kBlue + 1,
        "style": 1,
    },
    "nue": {
        "pdg": 12,
        "label": "#nu_{e}",
        "color": ROOT.kRed + 1,
        "style": 1,
    },
    "nuebar": {
        "pdg": -12,
        "label": "#bar{#nu}_{e}",
        "color": ROOT.kPink + 1,
        "style": 1,
    },
}

hist_name = "flux_E_cvweighted"


pot_map = {
    "minerva1":  9.605330838278293e19,
    "minerva13": 1.2512606809221975e20,
}

# ============================================================
# LOAD HISTOGRAM
# ============================================================

open_files = []


def get_flux_hist(pdg, playlist):

    filename = os.path.join(
        flux_dir,
        f"flux-gen2thin-pdg{pdg}-{playlist}.root"
    )

    print(f"Opening: {filename}")

    f = ROOT.TFile.Open(filename, "READ")

    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open {filename}")

    # Prevent PyROOT from trying to destroy the TFile at Python shutdown.
    ROOT.SetOwnership(f, False)

    open_files.append(f)

    h = f.Get(hist_name)

    if not h:
        raise RuntimeError(
            f"Could not find {hist_name} in {filename}"
        )

    # Histogram belongs to the TFile, not Python.
    ROOT.SetOwnership(h, False)

    return h


def get_combined_flux_hist(pdg, playlists):

    total_pot = sum(
        pot_map[p]
        for p in playlists
    )

    h_combined = None

    for playlist in playlists:

        h = get_flux_hist(
            pdg,
            playlist,
        )

        weight = pot_map[playlist] / total_pot

        print(
            f"  {playlist}: "
            f"POT = {pot_map[playlist]:.6e}, "
            f"weight = {weight:.6f}"
        )

        if h_combined is None:

            h_combined = h.Clone(
                f"h_flux_{pdg}_{'_'.join(playlists)}"
            )

            h_combined.SetDirectory(0)
            ROOT.SetOwnership(h_combined, False)

            h_combined.Scale(weight)

        else:

            h_combined.Add(
                h,
                weight,
            )

    return h_combined


def get_mode_flux_hist(pdg, mode):

    playlists = playlist_map[mode]

    if len(playlists) == 1:
        return get_flux_hist(
            pdg,
            playlists[0],
        )

    return get_combined_flux_hist(
        pdg,
        playlists,
    )

# ============================================================
# PLOT
# ============================================================

def make_flux_plot(mode):

    # --------------------------------------------------------
    # MINERvA plotter
    # --------------------------------------------------------

    plotter = MnvPlotter()

    # Apply MINERvA-style global settings
    plotter.ApplyStyle(7)

    # --------------------------------------------------------
    # Load fluxes
    # --------------------------------------------------------

    hists = {}

    for name, info in species.items():

        h = get_mode_flux_hist(
            info["pdg"],
            mode,
        )

        h.SetLineColor(info["color"])
        h.SetLineStyle(info["style"])
        h.SetLineWidth(3)

        h.SetMarkerStyle(0)

        hists[name] = h

        print(
            f"{mode:3s} {name:7s}: "
            f"Integral = {h.Integral():.6e}"
        )

    # --------------------------------------------------------
    # Canvas
    # --------------------------------------------------------

    c = ROOT.TCanvas(
        f"c_flux_{mode}",
        "",
        900,
        700,
    )

    c.SetLeftMargin(0.14)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(0.13)
    c.SetTopMargin(0.08)

    c.SetLogy()

    # --------------------------------------------------------
    # Axis range
    # --------------------------------------------------------

    xmin = 0.0
    xmax = 20.0

    ymax = 0.0
    ymin_positive = float("inf")

    for h in hists.values():

        bin_min = h.GetXaxis().FindBin(xmin)
        bin_max = h.GetXaxis().FindBin(xmax - 1e-6)

        for i in range(bin_min, bin_max + 1):

            content = h.GetBinContent(i)

            if content > 0.0:
                ymax = max(ymax, content)
                ymin_positive = min(ymin_positive, content)

    # Log-y plotting range
    ymin = 0.3 * ymin_positive
    ymax_plot = 3.0 * ymax

    print(
        f"{mode} plotting range: "
        f"{ymin:.3e} -- {ymax_plot:.3e}"
    )

    h0 = hists["numu"]

    h0.SetTitle("")

    h0.GetXaxis().SetTitle("Energy (GeV)")
    h0.GetYaxis().SetTitle(
        "Neutrinos / m^{2} / P.O.T. / GeV"
    )

    h0.GetXaxis().CenterTitle()
    h0.GetYaxis().CenterTitle()

    h0.GetXaxis().SetTitleSize(0.050)
    h0.GetYaxis().SetTitleSize(0.050)

    h0.GetXaxis().SetLabelSize(0.040)
    h0.GetYaxis().SetLabelSize(0.040)

    h0.GetXaxis().SetTitleOffset(1.10)
    h0.GetYaxis().SetTitleOffset(1.25)

    h0.GetXaxis().SetRangeUser(xmin, xmax)

    # Must be positive for log-y
    h0.SetMinimum(ymin)
    h0.SetMaximum(ymax_plot)

    # --------------------------------------------------------
    # Draw
    # --------------------------------------------------------

    draw_order = [
        "numu",
        "numubar",
        "nue",
        "nuebar",
    ]

    hists[draw_order[0]].Draw("HIST")

    for name in draw_order[1:]:
        hists[name].Draw("HIST SAME")

    # --------------------------------------------------------
    # Legend
    # --------------------------------------------------------

    legend = ROOT.TLegend(
        0.67,
        0.64,
        0.89,
        0.88,
    )

    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.040)

    for name in draw_order:

        legend.AddEntry(
            hists[name],
            species[name]["label"],
            "l",
        )

    legend.Draw()

    # --------------------------------------------------------
    # MINERvA label
    # --------------------------------------------------------

    plotter.AddHistoTitle(
        f"MINERvA {mode} Flux",
        0.05,
    )

    # --------------------------------------------------------
    # Redraw axes and save
    # --------------------------------------------------------

    c.RedrawAxis()

    output_base = f"minerva_flux_{mode}"

    c.Print(output_base + ".png")
    c.Print(output_base + ".pdf")

    print()
    print(f"Saved {output_base}.png")
    print(f"Saved {output_base}.pdf")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    make_flux_plot("FHC")
    make_flux_plot("RHC")
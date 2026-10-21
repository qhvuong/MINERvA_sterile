#!/usr/bin/env python

import argparse
import array
import csv
import math

import ROOT
import numpy as np


PDG_NUMU = 14
PDG_NUE = 12

# This is LE flux binning
LE_FLUX_BINNING = [
    0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,
    5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5,
    10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
    20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0,
]

# This is ME flux binning
ME_FLUX_BINNING = [
    0.0, 0.5, 1.0, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 
    2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 
    3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 
    4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 
    5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 
    6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 
    7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 
    8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 
    9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.5, 
    11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 
    16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 21.0, 
    22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 35.0, 
    40.0, 45.0, 50.0, 55.0, 60.0, 70.0, 80.0, 90.0, 100.0 
]



def add_filelist_to_chain(chain, filelist):
    nfiles = 0

    with open(filelist) as f:
        for line in f:
            filename = line.strip()

            if not filename or filename.startswith("#"):
                continue

            chain.Add(filename)
            nfiles += 1

    return nfiles


def passes_compute_flux_fv(vtx):
    x = vtx[0]
    y = vtx[1]
    z = vtx[2]

    if z <= 6117.0 or z >= 8193.0:
        return False

    ax = abs(x)
    ay = abs(y)

    if ax >= 850.0:
        return False

    sqrt3 = math.sqrt(3.0)

    inside_hexagon = (
        ay < 850.0 / sqrt3
        or
        ay < 850.0 * 2.0 / sqrt3 - ax / sqrt3
    )

    return inside_hexagon


def make_hist(name, title, bin_edges):
    h = ROOT.TH1D(
        name,
        title,
        len(bin_edges) - 1,
        array.array("d", bin_edges)
    )
    h.SetDirectory(0)
    return h


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Study ME BeamFocus uncertainty versus true neutrino energy "
            "for true numu and nue, including cross-flavor cancellation."
        )
    )

    parser.add_argument(
        "--filelist",
        action="append",
        required=True,
        help="ME tuple file list. Repeat to combine multiple playlists."
    )

    parser.add_argument(
        "--tree",
        default="Truth",
        help="Input tree name. Default: Truth"
    )

    parser.add_argument(
        "--n_universes",
        type=int,
        default=1000,
        help="Number of BeamFocus universes. Default: 1000"
    )

    parser.add_argument(
        "--max_events",
        type=int,
        default=-1,
        help="Maximum total entries to process. Default: all"
    )

    parser.add_argument(
        "--emin",
        type=float,
        default=0.0,
        help="Minimum true neutrino energy in GeV. Default: 0"
    )

    parser.add_argument(
        "--emax",
        type=float,
        default=20.0,
        help="Maximum true neutrino energy in GeV. Default: 20"
    )

    parser.add_argument(
        "--bin_width",
        type=float,
        default=0.5,
        help="True-E bin width in GeV. Default: 0.5"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output ROOT file"
    )

    parser.add_argument(
        "--csv",
        help="Optional summary CSV output"
    )

    parser.add_argument(
        "--npz",
        help="Optional NPZ output containing all binned universe weights"
    )

    args = parser.parse_args()

    # ------------------------------------------------------------
    # True-E binning
    # ------------------------------------------------------------

    n_bins = int(round((args.emax - args.emin) / args.bin_width))

    bin_edges = [
        args.emin + i * args.bin_width
        for i in range(n_bins + 1)
    ]

    print("True-E binning:")
    print("  Emin      =", args.emin)
    print("  Emax      =", args.emax)
    print("  bin width =", args.bin_width)
    print("  N bins    =", n_bins)

    # ------------------------------------------------------------
    # Build chain
    # ------------------------------------------------------------

    chain = ROOT.TChain(args.tree)

    total_files = 0

    for filelist in args.filelist:
        nfiles = add_filelist_to_chain(chain, filelist)
        total_files += nfiles
        print("Added {} files from {}".format(nfiles, filelist))

    print("Total files:", total_files)
    print("Total entries:", chain.GetEntries())

    # ------------------------------------------------------------
    # Required branches
    # ------------------------------------------------------------

    required_branches = [
        "mc_incoming",
        "mc_incomingE",
        "mc_current",
        "mc_intType",
        "mc_targetZ",
        "mc_vtx",
        "mc_cvweight_total",
        "mc_hornCurrent_cvweight",
        "mc_wgt_Flux_BeamFocus",
    ]

    for branch_name in required_branches:
        if not chain.GetBranch(branch_name):
            raise RuntimeError(
                "Required branch '{}' was not found.".format(branch_name)
            )

    chain.SetBranchStatus("*", 0)

    for branch_name in required_branches:
        chain.SetBranchStatus(branch_name, 1)

    # ------------------------------------------------------------
    # Branch buffers
    # ------------------------------------------------------------

    incoming = array.array("i", [0])
    incoming_e = array.array("d", [0.0])

    current = array.array("i", [0])
    int_type = array.array("i", [0])
    target_z = array.array("i", [0])

    vtx = array.array("d", [0.0, 0.0, 0.0])

    event_cv = array.array("d", [0.0])
    beam_cv = array.array("d", [0.0])

    beam_weights_int = array.array(
        "i",
        [0] * args.n_universes
    )

    chain.SetBranchAddress("mc_incoming", incoming)
    chain.SetBranchAddress("mc_incomingE", incoming_e)

    chain.SetBranchAddress("mc_current", current)
    chain.SetBranchAddress("mc_intType", int_type)
    chain.SetBranchAddress("mc_targetZ", target_z)

    chain.SetBranchAddress("mc_vtx", vtx)

    chain.SetBranchAddress("mc_cvweight_total", event_cv)
    chain.SetBranchAddress("mc_hornCurrent_cvweight", beam_cv)

    chain.SetBranchAddress(
        "mc_wgt_Flux_BeamFocus",
        beam_weights_int
    )

    # ------------------------------------------------------------
    # Accumulators
    #
    # species index:
    #   0 = numu
    #   1 = nue
    #
    # universe_sums[species, Ebin, universe]
    # ------------------------------------------------------------

    cv_sums = np.zeros(
        (2, n_bins),
        dtype=np.float64
    )

    universe_sums = np.zeros(
        (2, n_bins, args.n_universes),
        dtype=np.float64
    )

    event_counts = np.zeros(
        (2, n_bins),
        dtype=np.int64
    )

    skipped_counts = np.zeros(
        (2, n_bins),
        dtype=np.int64
    )

    # ------------------------------------------------------------
    # Event loop
    # ------------------------------------------------------------

    n_entries = chain.GetEntries()

    if args.max_events >= 0:
        n_entries = min(n_entries, args.max_events)

    for ientry in range(n_entries):

        if ientry % 100000 == 0:
            print(
                "Processing entry {} / {}".format(
                    ientry,
                    n_entries
                )
            )

        chain.GetEntry(ientry)

        pdg = incoming[0]

        if pdg == PDG_NUMU:
            ispecies = 0
        elif pdg == PDG_NUE:
            ispecies = 1
        else:
            continue

        # Match compute_flux truth selection
        if current[0] != 1:
            continue

        if int_type[0] == 8:
            continue

        if target_z[0] != 6:
            continue

        if not passes_compute_flux_fv(vtx):
            continue

        enu_gev = incoming_e[0] * 1.0e-3

        if enu_gev < args.emin or enu_gev >= args.emax:
            continue

        ibin = int(
            (enu_gev - args.emin) / args.bin_width
        )

        if ibin < 0 or ibin >= n_bins:
            continue

        event_counts[ispecies, ibin] += 1

        cv = event_cv[0]
        beam_cv_value = beam_cv[0]

        if (
            not math.isfinite(cv)
            or not math.isfinite(beam_cv_value)
            or cv < 1.0e-12
            or abs(beam_cv_value) < 1.0e-12
        ):
            skipped_counts[ispecies, ibin] += 1
            continue

        cv_sums[ispecies, ibin] += cv

        for u in range(args.n_universes):

            beam_absolute = (
                beam_weights_int[u] * 1.0e-7
            )

            beam_relative = (
                beam_absolute / beam_cv_value
            )

            universe_sums[ispecies, ibin, u] += (
                cv * beam_relative
            )

    # ------------------------------------------------------------
    # Convert accumulators into W_s,u(E)
    # ------------------------------------------------------------

    weights = np.full(
        (2, n_bins, args.n_universes),
        np.nan,
        dtype=np.float64
    )

    for ispecies in range(2):
        for ibin in range(n_bins):

            if cv_sums[ispecies, ibin] <= 0.0:
                continue

            weights[ispecies, ibin, :] = (
                universe_sums[ispecies, ibin, :]
                / cv_sums[ispecies, ibin]
            )

    w_numu = weights[0]
    w_nue = weights[1]

    # ------------------------------------------------------------
    # Per-bin diagnostics
    # ------------------------------------------------------------

    sigma_numu = np.full(n_bins, np.nan)
    sigma_nue = np.full(n_bins, np.nan)

    rho_numu_nue = np.full(n_bins, np.nan)

    sigma_ratio = np.full(n_bins, np.nan)
    mean_ratio = np.full(n_bins, np.nan)

    for ibin in range(n_bins):

        mu = w_numu[ibin]
        e = w_nue[ibin]

        valid = (
            np.isfinite(mu)
            & np.isfinite(e)
            & (np.abs(mu) > 1.0e-12)
        )

        if not np.any(valid):
            continue

        mu = mu[valid]
        e = e[valid]

        dmu = mu - 1.0
        de = e - 1.0

        # Systematic RMS relative to CV = 1
        sigma_numu[ibin] = np.sqrt(
            np.mean(dmu * dmu)
        )

        sigma_nue[ibin] = np.sqrt(
            np.mean(de * de)
        )

        # Correlation of CV-referenced systematic shifts
        denom = (
            sigma_numu[ibin]
            * sigma_nue[ibin]
        )

        if denom > 0.0:
            rho_numu_nue[ibin] = (
                np.mean(dmu * de)
                / denom
            )

        # Direct cancellation in nue / numu
        ratio = e / mu

        mean_ratio[ibin] = np.mean(ratio)

        sigma_ratio[ibin] = np.sqrt(
            np.mean(
                (ratio - 1.0) ** 2
            )
        )

    # ------------------------------------------------------------
    # ROOT output
    # ------------------------------------------------------------

    h_sigma_numu = make_hist(
        "beamfocus_sigma_numu",
        ";True E_{#nu} [GeV];BeamFocus fractional uncertainty",
        bin_edges
    )

    h_sigma_nue = make_hist(
        "beamfocus_sigma_nue",
        ";True E_{#nu} [GeV];BeamFocus fractional uncertainty",
        bin_edges
    )

    h_rho = make_hist(
        "beamfocus_rho_numu_nue",
        ";True E_{#nu} [GeV];#rho(#nu_{#mu}, #nu_{e})",
        bin_edges
    )

    h_sigma_ratio = make_hist(
        "beamfocus_sigma_nue_over_numu",
        ";True E_{#nu} [GeV];Residual #nu_{e}/#nu_{#mu} BeamFocus uncertainty",
        bin_edges
    )

    h_mean_ratio = make_hist(
        "beamfocus_mean_nue_over_numu",
        ";True E_{#nu} [GeV];Mean W_{#nu_{e}}/W_{#nu_{#mu}}",
        bin_edges
    )

    h_events_numu = make_hist(
        "event_count_numu",
        ";True E_{#nu} [GeV];Selected true #nu_{#mu} events",
        bin_edges
    )

    h_events_nue = make_hist(
        "event_count_nue",
        ";True E_{#nu} [GeV];Selected true #nu_{e} events",
        bin_edges
    )

    for ibin in range(n_bins):

        root_bin = ibin + 1

        if np.isfinite(sigma_numu[ibin]):
            h_sigma_numu.SetBinContent(
                root_bin,
                sigma_numu[ibin]
            )

        if np.isfinite(sigma_nue[ibin]):
            h_sigma_nue.SetBinContent(
                root_bin,
                sigma_nue[ibin]
            )

        if np.isfinite(rho_numu_nue[ibin]):
            h_rho.SetBinContent(
                root_bin,
                rho_numu_nue[ibin]
            )

        if np.isfinite(sigma_ratio[ibin]):
            h_sigma_ratio.SetBinContent(
                root_bin,
                sigma_ratio[ibin]
            )

        if np.isfinite(mean_ratio[ibin]):
            h_mean_ratio.SetBinContent(
                root_bin,
                mean_ratio[ibin]
            )

        h_events_numu.SetBinContent(
            root_bin,
            event_counts[0, ibin]
        )

        h_events_nue.SetBinContent(
            root_bin,
            event_counts[1, ibin]
        )

    fout = ROOT.TFile.Open(
        args.output,
        "RECREATE"
    )

    if not fout or fout.IsZombie():
        raise RuntimeError(
            "Could not create output ROOT file: {}".format(
                args.output
            )
        )

    fout.cd()

    h_sigma_numu.Write()
    h_sigma_nue.Write()
    h_rho.Write()
    h_sigma_ratio.Write()
    h_mean_ratio.Write()
    h_events_numu.Write()
    h_events_nue.Write()

    fout.Close()

    print("Wrote ROOT output:", args.output)

    # ------------------------------------------------------------
    # Optional NPZ
    # ------------------------------------------------------------

    if args.npz:
        np.savez(
            args.npz,
            bin_edges=np.array(bin_edges),
            weights_numu=w_numu,
            weights_nue=w_nue,
            sigma_numu=sigma_numu,
            sigma_nue=sigma_nue,
            rho_numu_nue=rho_numu_nue,
            mean_ratio=mean_ratio,
            sigma_ratio=sigma_ratio,
            event_counts_numu=event_counts[0],
            event_counts_nue=event_counts[1],
            skipped_counts_numu=skipped_counts[0],
            skipped_counts_nue=skipped_counts[1],
        )

        print("Wrote NPZ output:", args.npz)

    # ------------------------------------------------------------
    # Optional CSV
    # ------------------------------------------------------------

    if args.csv:

        with open(args.csv, "w") as csvfile:

            writer = csv.writer(csvfile)

            writer.writerow([
                "E_low",
                "E_high",
                "N_numu",
                "N_nue",
                "sigma_numu",
                "sigma_nue",
                "rho_numu_nue",
                "mean_nue_over_numu",
                "sigma_nue_over_numu",
            ])

            for ibin in range(n_bins):

                writer.writerow([
                    "{:.6f}".format(bin_edges[ibin]),
                    "{:.6f}".format(bin_edges[ibin + 1]),
                    int(event_counts[0, ibin]),
                    int(event_counts[1, ibin]),
                    "{:.10e}".format(sigma_numu[ibin]),
                    "{:.10e}".format(sigma_nue[ibin]),
                    "{:.10e}".format(rho_numu_nue[ibin]),
                    "{:.10e}".format(mean_ratio[ibin]),
                    "{:.10e}".format(sigma_ratio[ibin]),
                ])

        print("Wrote CSV output:", args.csv)

    # ------------------------------------------------------------
    # Print summary
    # ------------------------------------------------------------

    print()
    print("True-E BeamFocus summary")
    print("========================")

    for ibin in range(n_bins):

        if (
            event_counts[0, ibin] == 0
            and event_counts[1, ibin] == 0
        ):
            continue

        print(
            "{:.1f}-{:.1f} GeV : "
            "Nmu={:8d} Ne={:8d} "
            "sig_mu={:7.3f}% "
            "sig_e={:7.3f}% "
            "rho={:7.4f} "
            "sig_e/mu={:7.3f}%".format(
                bin_edges[ibin],
                bin_edges[ibin + 1],
                int(event_counts[0, ibin]),
                int(event_counts[1, ibin]),
                100.0 * sigma_numu[ibin]
                    if np.isfinite(sigma_numu[ibin])
                    else float("nan"),
                100.0 * sigma_nue[ibin]
                    if np.isfinite(sigma_nue[ibin])
                    else float("nan"),
                rho_numu_nue[ibin],
                100.0 * sigma_ratio[ibin]
                    if np.isfinite(sigma_ratio[ibin])
                    else float("nan"),
            )
        )


if __name__ == "__main__":
    main()
#!/usr/bin/env python

import os
import argparse
import array
import csv
import math

import ROOT
import numpy as np


SPECIES = [
    ("numu", 14, "#nu_{#mu}"),
    ("nue", 12, "#nu_{e}"),
    ("numubar", -14, "#bar{#nu}_{#mu}"),
    ("nuebar", -12, "#bar{#nu}_{e}"),
]
SPECIES_INDEX = {pdg: i for i, (_, pdg, _) in enumerate(SPECIES)}

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


def save_accumulator_checkpoint(
    accum,
    binnings,
    output_npz,
    completed_filelists,
    checkpoint_index=None,
):
    payload = {
        "species_names": np.array([x[0] for x in SPECIES]),
        "species_pdgs": np.array([x[1] for x in SPECIES], dtype=int),
        "completed_filelists": np.array(completed_filelists),
    }

    for tag, edges in binnings.items():
        key = tag.lower()
        a = accum[tag]

        payload["{}_bin_edges".format(key)] = np.array(edges)
        payload["{}_cv_sums".format(key)] = a["cv_sums"]
        payload["{}_universe_sums".format(key)] = a["universe_sums"]
        payload["{}_event_counts".format(key)] = a["event_counts"]
        payload["{}_skipped_counts".format(key)] = a["skipped_counts"]

    base = output_npz
    if base.lower().endswith(".npz"):
        base = base[:-4]

    if checkpoint_index is None:
        checkpoint_path = base + "_accumulators.npz"
    else:
        last_filelist = os.path.basename(completed_filelists[-1])
        if last_filelist.endswith(".txt"):
            last_filelist = last_filelist[:-4]

        checkpoint_path = (
            "{}_accumulators_{:02d}_{}.npz".format(
                base,
                checkpoint_index,
                last_filelist,
            )
        )

    # Write atomically: don't leave a corrupted checkpoint if the write dies.
    tmp_path = checkpoint_path + ".tmp.npz"

    np.savez(
        tmp_path,
        **payload
    )

    os.replace(
        tmp_path,
        checkpoint_path,
    )

    print("Wrote cumulative accumulator checkpoint:", checkpoint_path)

    return checkpoint_path


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Study ME BeamFocus uncertainty versus true neutrino energy "
            "for true numu, nue, numubar, and nuebar, including "
            "same-sign cross-flavor cancellation."
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
        help=(
            "Maximum total entries to process across all playlists. "
            "Default: all"
        )
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
    # True-E binnings
    # ------------------------------------------------------------

    binnings = {
        "ME": ME_FLUX_BINNING,
        "LE": LE_FLUX_BINNING,
    }

    for tag, edges in binnings.items():
        print(
            "{} true-E binning: {} bins, {}--{} GeV".format(
                tag,
                len(edges) - 1,
                edges[0],
                edges[-1],
            )
        )

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

    # ------------------------------------------------------------
    # Cumulative accumulators
    #
    # species order:
    #   0 numu
    #   1 nue
    #   2 numubar
    #   3 nuebar
    #
    # These arrays persist across ALL playlists.
    # ------------------------------------------------------------

    accum = {}

    for tag, edges in binnings.items():
        n_bins = len(edges) - 1

        accum[tag] = {
            "cv_sums": np.zeros(
                (4, n_bins),
                dtype=np.float64,
            ),
            "universe_sums": np.zeros(
                (4, n_bins, args.n_universes),
                dtype=np.float64,
            ),
            "event_counts": np.zeros(
                (4, n_bins),
                dtype=np.int64,
            ),
            "skipped_counts": np.zeros(
                (4, n_bins),
                dtype=np.int64,
            ),
        }

    # ------------------------------------------------------------
    # Global bookkeeping
    # ------------------------------------------------------------

    total_files = 0
    total_entries_available = 0
    total_entries_processed = 0

    total_fv_seen = 0
    total_fv_pass = 0

    completed_filelists = []

    # ------------------------------------------------------------
    # Process one playlist/filelist at a time
    # ------------------------------------------------------------

    for ifilelist, filelist in enumerate(
        args.filelist,
        start=1,
    ):

        # Respect --max_events as a GLOBAL limit across all playlists.
        if (
            args.max_events >= 0
            and total_entries_processed >= args.max_events
        ):
            print()
            print(
                "Reached global --max_events limit of {} entries.".format(
                    args.max_events
                )
            )
            break

        # --------------------------------------------------------
        # Build chain for THIS playlist only
        # --------------------------------------------------------

        chain = ROOT.TChain(args.tree)

        nfiles = add_filelist_to_chain(
            chain,
            filelist,
        )

        total_files += nfiles

        playlist_entries_available = chain.GetEntries()
        total_entries_available += playlist_entries_available

        print()
        print("=" * 100)
        print(
            "Starting playlist {}/{}: {}".format(
                ifilelist,
                len(args.filelist),
                filelist,
            )
        )
        print("Files:", nfiles)
        print("Entries available:", playlist_entries_available)
        print("=" * 100)

        # --------------------------------------------------------
        # Check required branches for THIS chain
        # --------------------------------------------------------

        for branch_name in required_branches:
            if not chain.GetBranch(branch_name):
                raise RuntimeError(
                    "Required branch '{}' was not found in playlist {}.".format(
                        branch_name,
                        filelist,
                    )
                )

        chain.SetBranchStatus("*", 0)

        for branch_name in required_branches:
            chain.SetBranchStatus(
                branch_name,
                1,
            )

        # --------------------------------------------------------
        # Branch buffers for THIS chain
        # --------------------------------------------------------

        incoming = array.array(
            "i",
            [0],
        )

        incoming_e = array.array(
            "d",
            [0.0],
        )

        current = array.array(
            "i",
            [0],
        )

        int_type = array.array(
            "i",
            [0],
        )

        target_z = array.array(
            "i",
            [0],
        )

        vtx = array.array(
            "d",
            [0.0, 0.0, 0.0],
        )

        event_cv = array.array(
            "d",
            [0.0],
        )

        beam_cv = array.array(
            "d",
            [0.0],
        )

        beam_weights_int = array.array(
            "i",
            [0] * args.n_universes,
        )

        chain.SetBranchAddress(
            "mc_incoming",
            incoming,
        )

        chain.SetBranchAddress(
            "mc_incomingE",
            incoming_e,
        )

        chain.SetBranchAddress(
            "mc_current",
            current,
        )

        chain.SetBranchAddress(
            "mc_intType",
            int_type,
        )

        chain.SetBranchAddress(
            "mc_targetZ",
            target_z,
        )

        chain.SetBranchAddress(
            "mc_vtx",
            vtx,
        )

        chain.SetBranchAddress(
            "mc_cvweight_total",
            event_cv,
        )

        chain.SetBranchAddress(
            "mc_hornCurrent_cvweight",
            beam_cv,
        )

        chain.SetBranchAddress(
            "mc_wgt_Flux_BeamFocus",
            beam_weights_int,
        )

        # --------------------------------------------------------
        # Determine how much of THIS playlist to process
        # --------------------------------------------------------

        playlist_n_entries = playlist_entries_available

        if args.max_events >= 0:
            remaining = (
                args.max_events
                - total_entries_processed
            )

            playlist_n_entries = min(
                playlist_n_entries,
                remaining,
            )

        print(
            "Entries to process from this playlist:",
            playlist_n_entries,
        )

        playlist_fv_seen = 0
        playlist_fv_pass = 0

        # --------------------------------------------------------
        # Event loop for THIS playlist
        # --------------------------------------------------------

        for ientry in range(playlist_n_entries):

            if ientry % 100000 == 0:
                print(
                    "Playlist {}/{}: processing entry {} / {} "
                    "(global processed so far: {})".format(
                        ifilelist,
                        len(args.filelist),
                        ientry,
                        playlist_n_entries,
                        total_entries_processed + ientry,
                    )
                )

            chain.GetEntry(ientry)

            pdg = incoming[0]

            if pdg not in SPECIES_INDEX:
                continue

            ispecies = SPECIES_INDEX[pdg]

            # Match compute_flux truth selection.
            if current[0] != 1:
                continue

            if int_type[0] == 8:
                continue

            if target_z[0] != 6:
                continue

            playlist_fv_seen += 1
            total_fv_seen += 1

            if not passes_compute_flux_fv(vtx):
                continue

            playlist_fv_pass += 1
            total_fv_pass += 1

            enu_gev = incoming_e[0] * 1.0e-3

            cv = event_cv[0]
            beam_cv_value = beam_cv[0]

            valid_weight = (
                math.isfinite(cv)
                and math.isfinite(beam_cv_value)
                and cv >= 1.0e-12
                and abs(beam_cv_value) >= 1.0e-12
            )

            # Build all relative BeamFocus universe weights once/event.
            beam_relative = None

            if valid_weight:
                beam_relative = np.empty(
                    args.n_universes,
                    dtype=np.float64,
                )

                for u in range(args.n_universes):
                    beam_absolute = (
                        beam_weights_int[u]
                        * 1.0e-7
                    )

                    beam_relative[u] = (
                        beam_absolute
                        / beam_cv_value
                    )

            # Fill the exact same selected event into both
            # ME and LE true-E binnings.
            for tag, edges in binnings.items():

                if (
                    enu_gev < edges[0]
                    or enu_gev >= edges[-1]
                ):
                    continue

                ibin = int(
                    np.searchsorted(
                        edges,
                        enu_gev,
                        side="right",
                    )
                    - 1
                )

                if (
                    ibin < 0
                    or ibin >= len(edges) - 1
                ):
                    continue

                a = accum[tag]

                a["event_counts"][
                    ispecies,
                    ibin,
                ] += 1

                if not valid_weight:
                    a["skipped_counts"][
                        ispecies,
                        ibin,
                    ] += 1
                    continue

                a["cv_sums"][
                    ispecies,
                    ibin,
                ] += cv

                a["universe_sums"][
                    ispecies,
                    ibin,
                    :,
                ] += (
                    cv
                    * beam_relative
                )

        # --------------------------------------------------------
        # Playlist completed
        # --------------------------------------------------------

        total_entries_processed += playlist_n_entries

        print()
        print(
            "Finished playlist {}/{}: {}".format(
                ifilelist,
                len(args.filelist),
                filelist,
            )
        )

        print(
            "Playlist FV candidates after CC/intType/target cuts:",
            playlist_fv_seen,
        )

        print(
            "Playlist FV passing:",
            playlist_fv_pass,
        )

        print(
            "Cumulative entries processed:",
            total_entries_processed,
        )

        print(
            "Cumulative FV candidates:",
            total_fv_seen,
        )

        print(
            "Cumulative FV passing:",
            total_fv_pass,
        )

        # This is a true completed playlist only if we processed
        # every available entry in it.
        playlist_completed_fully = (
            playlist_n_entries
            == playlist_entries_available
        )

        if playlist_completed_fully:
            completed_filelists.append(
                filelist
            )

        # --------------------------------------------------------
        # Cumulative checkpoint after this playlist
        # --------------------------------------------------------

        if args.npz:

            # If --max_events stopped us in the middle of a playlist,
            # still save the accumulated work, but identify that this
            # checkpoint includes a partial final playlist.
            checkpoint_filelists = list(
                completed_filelists
            )

            if not playlist_completed_fully:
                checkpoint_filelists.append(
                    filelist + " [PARTIAL]"
                )

            save_accumulator_checkpoint(
                accum=accum,
                binnings=binnings,
                output_npz=args.npz,
                completed_filelists=checkpoint_filelists,
                checkpoint_index=ifilelist,
            )

        # If this playlist was truncated by --max_events,
        # the global limit has been reached.
        if not playlist_completed_fully:
            print(
                "Stopped partway through this playlist because "
                "--max_events was reached."
            )
            break

    # ------------------------------------------------------------
    # Final cumulative bookkeeping
    # ------------------------------------------------------------

    print()
    print("=" * 100)
    print("Finished event processing")
    print("=" * 100)

    print("Total files added:", total_files)
    print("Total entries available in visited playlists:", total_entries_available)
    print("Total entries processed:", total_entries_processed)
    print("Total FV candidates:", total_fv_seen)
    print("Total FV passing:", total_fv_pass)

    print(
        "Fully completed playlists: {}/{}".format(
            len(completed_filelists),
            len(args.filelist),
        )
    )

    # ------------------------------------------------------------
    # Final canonical accumulator checkpoint
    # ------------------------------------------------------------

    if args.npz:

        # For your normal full-stat run all playlists will be complete,
        # so completed_filelists is the complete ordered playlist list.
        #
        # Save one canonical cumulative file in addition to the
        # numbered per-playlist checkpoints.
        final_filelists = list(
            completed_filelists
        )

        # This should only matter for a --max_events test.
        if (
            len(final_filelists) == 0
            and len(args.filelist) > 0
        ):
            final_filelists = [
                args.filelist[0] + " [PARTIAL]"
            ]

        save_accumulator_checkpoint(
            accum=accum,
            binnings=binnings,
            output_npz=args.npz,
            completed_filelists=final_filelists,
            checkpoint_index=None,
        )

    # ------------------------------------------------------------
    # Convert accumulators and write outputs
    # ------------------------------------------------------------

    results = {}

    for tag, edges in binnings.items():
        a = accum[tag]
        n_bins = len(edges) - 1

        weights = np.full(
            (4, n_bins, args.n_universes), np.nan, dtype=np.float64
        )

        # Keep both uncertainty definitions:
        #   sigma_cv       : RMS displacement from nominal CV (=1)
        #   sigma_centered : standard deviation about universe mean
        # Also save the mean displacement from CV.
        sigma_cv = np.full((4, n_bins), np.nan, dtype=np.float64)
        sigma_centered = np.full((4, n_bins), np.nan, dtype=np.float64)
        mean_shift = np.full((4, n_bins), np.nan, dtype=np.float64)

        for ispecies in range(4):
            for ibin in range(n_bins):
                if a["cv_sums"][ispecies, ibin] <= 0.0:
                    continue

                weights[ispecies, ibin, :] = (
                    a["universe_sums"][ispecies, ibin, :]
                    / a["cv_sums"][ispecies, ibin]
                )

                w = weights[ispecies, ibin, :]
                valid = np.isfinite(w)
                if not np.any(valid):
                    continue

                w = w[valid]
                d = w - 1.0

                mean_shift[ispecies, ibin] = np.mean(d)
                sigma_cv[ispecies, ibin] = np.sqrt(np.mean(d * d))

                if len(w) > 1:
                    sigma_centered[ispecies, ibin] = np.std(
                        w,
                        ddof=1,
                    )

        # Same-sign flavor cancellation diagnostics:
        #   nue/numu and nuebar/numubar.
        pair_defs = {
            "numu_over_nue": (SPECIES_INDEX[14], SPECIES_INDEX[12]),
            "numubar_over_nuebar": (SPECIES_INDEX[-14], SPECIES_INDEX[-12]),
        }

        pairs = {}

        for pair_name, (inum, iden) in pair_defs.items():
            # Save both correlation conventions too:
            #   rho_cv       : CV-centered, paired with RMS-to-CV sigmas
            #   rho_centered : ordinary Pearson correlation
            rho_cv = np.full(n_bins, np.nan)
            rho_centered = np.full(n_bins, np.nan)

            mean_ratio = np.full(n_bins, np.nan)
            mean_ratio_shift = np.full(n_bins, np.nan)
            sigma_ratio_cv = np.full(n_bins, np.nan)
            sigma_ratio_centered = np.full(n_bins, np.nan)

            for ibin in range(n_bins):
                num = weights[inum, ibin, :]
                den = weights[iden, ibin, :]

                valid = (
                    np.isfinite(num)
                    & np.isfinite(den)
                    & (np.abs(den) > 1.0e-12)
                )
                if not np.any(valid):
                    continue

                num = num[valid]
                den = den[valid]

                dnum = num - 1.0
                dden = den - 1.0

                snum_cv = np.sqrt(np.mean(dnum * dnum))
                sden_cv = np.sqrt(np.mean(dden * dden))

                if snum_cv > 0.0 and sden_cv > 0.0:
                    rho_cv[ibin] = (
                        np.mean(dnum * dden)
                        / (snum_cv * sden_cv)
                    )

                if len(num) > 1:
                    rho_centered[ibin] = np.corrcoef(
                        num,
                        den,
                    )[0, 1]

                ratio = num / den
                delta_ratio = ratio - 1.0

                mean_ratio[ibin] = np.mean(ratio)
                mean_ratio_shift[ibin] = np.mean(delta_ratio)
                sigma_ratio_cv[ibin] = np.sqrt(
                    np.mean(delta_ratio * delta_ratio)
                )

                if len(ratio) > 1:
                    sigma_ratio_centered[ibin] = np.std(
                        ratio,
                        ddof=1,
                    )

            pairs[pair_name] = {
                "rho_cv": rho_cv,
                "rho_centered": rho_centered,
                "mean_ratio": mean_ratio,
                "mean_ratio_shift": mean_ratio_shift,
                "sigma_ratio_cv": sigma_ratio_cv,
                "sigma_ratio_centered": sigma_ratio_centered,
            }

        results[tag] = {
            "weights": weights,
            "sigma_cv": sigma_cv,
            "sigma_centered": sigma_centered,
            "mean_shift": mean_shift,
            "pairs": pairs,
            "event_counts": a["event_counts"],
            "skipped_counts": a["skipped_counts"],
        }

    # ------------------------------------------------------------
    # ROOT output
    # ------------------------------------------------------------

    fout = ROOT.TFile.Open(args.output, "RECREATE")

    if not fout or fout.IsZombie():
        raise RuntimeError(
            "Could not create output ROOT file: {}".format(args.output)
        )

    for tag, edges in binnings.items():
        r = results[tag]
        fout.mkdir(tag)
        fout.cd(tag)

        for ispecies, (name, pdg, latex) in enumerate(SPECIES):
            h_sigma_cv = make_hist(
                "beamfocus_sigma_cv_{}".format(name),
                ";True E_{#nu} [GeV];BeamFocus RMS uncertainty from CV",
                edges,
            )
            h_sigma_centered = make_hist(
                "beamfocus_sigma_centered_{}".format(name),
                ";True E_{#nu} [GeV];BeamFocus centered fractional uncertainty",
                edges,
            )
            h_mean_shift = make_hist(
                "beamfocus_mean_shift_{}".format(name),
                ";True E_{#nu} [GeV];Mean BeamFocus shift from CV",
                edges,
            )
            h_events = make_hist(
                "event_count_{}".format(name),
                ";True E_{{#nu}} [GeV];Selected true {} events".format(latex),
                edges,
            )

            for ibin in range(len(edges) - 1):
                rb = ibin + 1

                if np.isfinite(r["sigma_cv"][ispecies, ibin]):
                    h_sigma_cv.SetBinContent(
                        rb, r["sigma_cv"][ispecies, ibin]
                    )

                if np.isfinite(r["sigma_centered"][ispecies, ibin]):
                    h_sigma_centered.SetBinContent(
                        rb, r["sigma_centered"][ispecies, ibin]
                    )

                if np.isfinite(r["mean_shift"][ispecies, ibin]):
                    h_mean_shift.SetBinContent(
                        rb, r["mean_shift"][ispecies, ibin]
                    )

                h_events.SetBinContent(
                    rb, r["event_counts"][ispecies, ibin]
                )

            h_sigma_cv.Write()
            h_sigma_centered.Write()
            h_mean_shift.Write()
            h_events.Write()

        for pair_name, pr in r["pairs"].items():
            h_rho_cv = make_hist(
                "beamfocus_rho_cv_{}".format(pair_name),
                ";True E_{#nu} [GeV];CV-centered BeamFocus correlation",
                edges,
            )
            h_rho_centered = make_hist(
                "beamfocus_rho_centered_{}".format(pair_name),
                ";True E_{#nu} [GeV];Pearson BeamFocus correlation",
                edges,
            )
            h_sigma_ratio_cv = make_hist(
                "beamfocus_sigma_cv_{}".format(pair_name),
                ";True E_{#nu} [GeV];Residual BeamFocus RMS uncertainty from CV",
                edges,
            )
            h_sigma_ratio_centered = make_hist(
                "beamfocus_sigma_centered_{}".format(pair_name),
                ";True E_{#nu} [GeV];Residual centered BeamFocus uncertainty",
                edges,
            )
            h_mean_ratio = make_hist(
                "beamfocus_mean_{}".format(pair_name),
                ";True E_{#nu} [GeV];Mean universe ratio",
                edges,
            )
            h_mean_ratio_shift = make_hist(
                "beamfocus_mean_shift_{}".format(pair_name),
                ";True E_{#nu} [GeV];Mean universe ratio shift from CV",
                edges,
            )

            for ibin in range(len(edges) - 1):
                rb = ibin + 1

                if np.isfinite(pr["rho_cv"][ibin]):
                    h_rho_cv.SetBinContent(rb, pr["rho_cv"][ibin])
                if np.isfinite(pr["rho_centered"][ibin]):
                    h_rho_centered.SetBinContent(
                        rb, pr["rho_centered"][ibin]
                    )
                if np.isfinite(pr["sigma_ratio_cv"][ibin]):
                    h_sigma_ratio_cv.SetBinContent(
                        rb, pr["sigma_ratio_cv"][ibin]
                    )
                if np.isfinite(pr["sigma_ratio_centered"][ibin]):
                    h_sigma_ratio_centered.SetBinContent(
                        rb, pr["sigma_ratio_centered"][ibin]
                    )
                if np.isfinite(pr["mean_ratio"][ibin]):
                    h_mean_ratio.SetBinContent(rb, pr["mean_ratio"][ibin])
                if np.isfinite(pr["mean_ratio_shift"][ibin]):
                    h_mean_ratio_shift.SetBinContent(
                        rb, pr["mean_ratio_shift"][ibin]
                    )

            h_rho_cv.Write()
            h_rho_centered.Write()
            h_sigma_ratio_cv.Write()
            h_sigma_ratio_centered.Write()
            h_mean_ratio.Write()
            h_mean_ratio_shift.Write()

        fout.cd()

    fout.Close()
    print("Wrote ROOT output:", args.output)

    # ------------------------------------------------------------
    # Optional NPZ
    # ------------------------------------------------------------

    if args.npz:
        payload = {
            "species_names": np.array([x[0] for x in SPECIES]),
            "species_pdgs": np.array([x[1] for x in SPECIES], dtype=int),
        }

        for tag, edges in binnings.items():
            key = tag.lower()
            r = results[tag]

            payload["{}_bin_edges".format(key)] = np.array(edges)
            payload["{}_weights".format(key)] = r["weights"]
            payload["{}_sigma_cv".format(key)] = r["sigma_cv"]
            payload["{}_sigma_centered".format(key)] = r["sigma_centered"]
            payload["{}_mean_shift".format(key)] = r["mean_shift"]
            payload["{}_event_counts".format(key)] = r["event_counts"]
            payload["{}_skipped_counts".format(key)] = r["skipped_counts"]

            for pair_name, pr in r["pairs"].items():
                payload["{}_rho_cv_{}".format(key, pair_name)] = pr["rho_cv"]
                payload["{}_rho_centered_{}".format(key, pair_name)] = (
                    pr["rho_centered"]
                )
                payload["{}_mean_{}".format(key, pair_name)] = pr["mean_ratio"]
                payload["{}_mean_shift_{}".format(key, pair_name)] = (
                    pr["mean_ratio_shift"]
                )
                payload["{}_sigma_cv_{}".format(key, pair_name)] = (
                    pr["sigma_ratio_cv"]
                )
                payload["{}_sigma_centered_{}".format(key, pair_name)] = (
                    pr["sigma_ratio_centered"]
                )

        np.savez(args.npz, **payload)
        print("Wrote NPZ output:", args.npz)

    # ------------------------------------------------------------
    # Optional CSV: one file per binning
    # ------------------------------------------------------------

    if args.csv:
        if args.csv.lower().endswith(".csv"):
            csv_base = args.csv[:-4]
        else:
            csv_base = args.csv

        for tag, edges in binnings.items():
            csv_path = "{}_{}.csv".format(csv_base, tag)
            r = results[tag]

            with open(csv_path, "w") as csvfile:
                writer = csv.writer(csvfile)

                writer.writerow([
                    "E_low", "E_high",
                    "N_numu", "N_nue", "N_numubar", "N_nuebar",

                    "mean_shift_numu", "mean_shift_nue",
                    "mean_shift_numubar", "mean_shift_nuebar",

                    "sigma_cv_numu", "sigma_cv_nue",
                    "sigma_cv_numubar", "sigma_cv_nuebar",

                    "sigma_centered_numu", "sigma_centered_nue",
                    "sigma_centered_numubar", "sigma_centered_nuebar",

                    "rho_cv_numu_nue",
                    "rho_centered_numu_nue",
                    "mean_shift_numu_over_nue",
                    "sigma_cv_numu_over_nue",
                    "sigma_centered_numu_over_nue",

                    "rho_cv_numubar_nuebar",
                    "rho_centered_numubar_nuebar",
                    "mean_shift_numubar_over_nuebar",
                    "sigma_cv_numubar_over_nuebar",
                    "sigma_centered_numubar_over_nuebar",
                ])

                pnu = r["pairs"]["numu_over_nue"]
                pbar = r["pairs"]["numubar_over_nuebar"]

                for ibin in range(len(edges) - 1):
                    writer.writerow([
                        "{:.6f}".format(edges[ibin]),
                        "{:.6f}".format(edges[ibin + 1]),

                        *[
                            int(r["event_counts"][sidx, ibin])
                            for sidx in range(4)
                        ],

                        *[
                            "{:.10e}".format(r["mean_shift"][sidx, ibin])
                            for sidx in range(4)
                        ],

                        *[
                            "{:.10e}".format(r["sigma_cv"][sidx, ibin])
                            for sidx in range(4)
                        ],

                        *[
                            "{:.10e}".format(
                                r["sigma_centered"][sidx, ibin]
                            )
                            for sidx in range(4)
                        ],

                        "{:.10e}".format(pnu["rho_cv"][ibin]),
                        "{:.10e}".format(pnu["rho_centered"][ibin]),
                        "{:.10e}".format(pnu["mean_ratio_shift"][ibin]),
                        "{:.10e}".format(pnu["sigma_ratio_cv"][ibin]),
                        "{:.10e}".format(
                            pnu["sigma_ratio_centered"][ibin]
                        ),

                        "{:.10e}".format(pbar["rho_cv"][ibin]),
                        "{:.10e}".format(pbar["rho_centered"][ibin]),
                        "{:.10e}".format(pbar["mean_ratio_shift"][ibin]),
                        "{:.10e}".format(pbar["sigma_ratio_cv"][ibin]),
                        "{:.10e}".format(
                            pbar["sigma_ratio_centered"][ibin]
                        ),
                    ])

            print("Wrote CSV output:", csv_path)

    # ------------------------------------------------------------
    # Print summaries
    # ------------------------------------------------------------

    for tag, edges in binnings.items():
        r = results[tag]

        print()
        print("{}-binned BeamFocus summary".format(tag))
        print("=" * 100)

        for ibin in range(len(edges) - 1):
            if np.sum(r["event_counts"][:, ibin]) == 0:
                continue

            sig_cv = [
                100.0 * r["sigma_cv"][sidx, ibin]
                if np.isfinite(r["sigma_cv"][sidx, ibin])
                else float("nan")
                for sidx in range(4)
            ]

            sig_centered = [
                100.0 * r["sigma_centered"][sidx, ibin]
                if np.isfinite(r["sigma_centered"][sidx, ibin])
                else float("nan")
                for sidx in range(4)
            ]

            pnu = r["pairs"]["numu_over_nue"]
            pbar = r["pairs"]["numubar_over_nuebar"]

            print(
                "{:.1f}-{:.1f} GeV : "
                "Nmu={:8d} Ne={:8d} Nmubar={:8d} Nebar={:8d} | "
                "CV-RMS: mu={:7.3f}% e={:7.3f}% "
                "mubar={:7.3f}% ebar={:7.3f}% | "
                "centered: mu={:7.3f}% e={:7.3f}% "
                "mubar={:7.3f}% ebar={:7.3f}% | "
                "mu/e: rhoCV={:7.4f} rhoP={:7.4f} "
                "sigCV={:7.3f}% sigCtr={:7.3f}% | "
                "mubar/ebar: rhoCV={:7.4f} rhoP={:7.4f} "
                "sigCV={:7.3f}% sigCtr={:7.3f}%".format(
                    edges[ibin],
                    edges[ibin + 1],
                    *[
                        int(r["event_counts"][sidx, ibin])
                        for sidx in range(4)
                    ],
                    *sig_cv,
                    *sig_centered,
                    pnu["rho_cv"][ibin],
                    pnu["rho_centered"][ibin],
                    100.0 * pnu["sigma_ratio_cv"][ibin]
                    if np.isfinite(pnu["sigma_ratio_cv"][ibin])
                    else float("nan"),
                    100.0 * pnu["sigma_ratio_centered"][ibin]
                    if np.isfinite(pnu["sigma_ratio_centered"][ibin])
                    else float("nan"),
                    pbar["rho_cv"][ibin],
                    pbar["rho_centered"][ibin],
                    100.0 * pbar["sigma_ratio_cv"][ibin]
                    if np.isfinite(pbar["sigma_ratio_cv"][ibin])
                    else float("nan"),
                    100.0 * pbar["sigma_ratio_centered"][ibin]
                    if np.isfinite(pbar["sigma_ratio_centered"][ibin])
                    else float("nan"),
                )
            )


if __name__ == "__main__":
    main()

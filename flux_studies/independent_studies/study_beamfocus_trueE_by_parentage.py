#!/usr/bin/env python

import os
import argparse
import array
import csv
import math
from collections import Counter

import ROOT
import numpy as np


SPECIES = [
    ("numu", 14, "#nu_{#mu}"),
    ("nue", 12, "#nu_{e}"),
    ("numubar", -14, "#bar{#nu}_{#mu}"),
    ("nuebar", -12, "#bar{#nu}_{e}"),
]
SPECIES_INDEX = {pdg: i for i, (_, pdg, _) in enumerate(SPECIES)}

# "all" is an inclusive cross-check against study_beamfocus_trueE_all_species.py.
# The remaining categories are mutually exclusive.
PARENT_CATEGORIES = [
    "all",
    "pion",
    "muon",
    "kaon_charged",
    "kaon_neutral",
    "other",
]
PARENT_INDEX = {name: i for i, name in enumerate(PARENT_CATEGORIES)}

# Common possible names. If none is found, pass --parent-branch explicitly.
PARENT_BRANCH_CANDIDATES = [
    "mc_fr_nuParentID",
    "mc_parent_pdg",
    "mc_parentPDG",
    "mc_parent",
    "parent_pdg",
    "parentPDG",
]

# LE flux binning
LE_FLUX_BINNING = [
    0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,
    5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5,
    10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
    20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0,
]

# ME flux binning
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
    40.0, 45.0, 50.0, 55.0, 60.0, 70.0, 80.0, 90.0, 100.0,
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

    return (
        ay < 850.0 / sqrt3
        or
        ay < 850.0 * 2.0 / sqrt3 - ax / sqrt3
    )


def classify_parent(parent_pdg):
    """
    Classify the IMMEDIATE neutrino parent by PDG.

    pion:
        pi+/- (|PDG| = 211)

    muon:
        mu+/- (|PDG| = 13)

    kaon_charged:
        K+/- (|PDG| = 321)

    kaon_neutral:
        K0 / K0L / K0S (311, 130, 310)

    other:
        everything else

    We intentionally keep charged and neutral kaons separate because that
    can be useful for intrinsic nue/nuebar parentage. They can be combined
    later if desired.
    """
    apdg = abs(int(parent_pdg))

    if apdg == 211:
        return "pion"

    if apdg == 13:
        return "muon"

    if apdg == 321:
        return "kaon_charged"

    if apdg in (130, 310, 311):
        return "kaon_neutral"

    return "other"


def find_parent_branch(chain, requested=None):
    if requested:
        if not chain.GetBranch(requested):
            raise RuntimeError(
                "Requested parent branch '{}' was not found.".format(requested)
            )
        return requested

    found = [
        name for name in PARENT_BRANCH_CANDIDATES
        if chain.GetBranch(name)
    ]

    if len(found) == 1:
        return found[0]

    if len(found) > 1:
        raise RuntimeError(
            "Multiple possible parent branches found: {}. "
            "Please choose one with --parent-branch.".format(found)
        )

    parent_like = []
    branches = chain.GetListOfBranches()
    for i in range(branches.GetEntries()):
        name = branches.At(i).GetName()
        if "parent" in name.lower():
            parent_like.append(name)

    raise RuntimeError(
        "Could not auto-detect the parent-PDG branch. "
        "Branches containing 'parent' are: {}. "
        "Rerun with --parent-branch BRANCHNAME.".format(parent_like)
    )


def make_hist(name, title, bin_edges):
    h = ROOT.TH1D(
        name,
        title,
        len(bin_edges) - 1,
        array.array("d", bin_edges),
    )
    h.SetDirectory(0)
    return h


def save_accumulator_checkpoint(
    accum,
    binnings,
    output_npz,
    completed_filelists,
    parent_branch,
    checkpoint_index=None,
):
    payload = {
        "species_names": np.array([x[0] for x in SPECIES]),
        "species_pdgs": np.array([x[1] for x in SPECIES], dtype=int),
        "parent_categories": np.array(PARENT_CATEGORIES),
        "parent_branch": np.array(parent_branch),
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

    base = output_npz[:-4] if output_npz.lower().endswith(".npz") else output_npz

    if checkpoint_index is None:
        checkpoint_path = base + "_accumulators.npz"
    else:
        last_filelist = os.path.basename(completed_filelists[-1])
        if last_filelist.endswith(".txt"):
            last_filelist = last_filelist[:-4]

        checkpoint_path = "{}_accumulators_{:02d}_{}.npz".format(
            base,
            checkpoint_index,
            last_filelist,
        )

    tmp_path = checkpoint_path + ".tmp.npz"
    np.savez(tmp_path, **payload)
    os.replace(tmp_path, checkpoint_path)

    print("Wrote cumulative accumulator checkpoint:", checkpoint_path)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Study BeamFocus uncertainty versus true neutrino energy, "
            "separated by neutrino flavor and immediate parent PDG."
        )
    )

    parser.add_argument(
        "--filelist",
        action="append",
        required=True,
        help="Tuple file list. Repeat to combine multiple playlists.",
    )

    parser.add_argument(
        "--tree",
        default="Truth",
        help="Input tree name. Default: Truth",
    )

    parser.add_argument(
        "--parent-branch",
        default=None,
        help=(
            "Branch containing the immediate neutrino-parent PDG. "
            "If omitted, try to auto-detect it."
        ),
    )

    parser.add_argument(
        "--n_universes",
        type=int,
        default=1000,
        help="Number of BeamFocus universes. Default: 1000",
    )

    parser.add_argument(
        "--max_events",
        type=int,
        default=-1,
        help="Maximum total entries across all playlists. Default: all",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output ROOT file",
    )

    parser.add_argument(
        "--csv",
        help="Optional summary CSV basename",
    )

    parser.add_argument(
        "--npz",
        help="Optional NPZ output",
    )

    args = parser.parse_args()

    if args.n_universes != 1000:
        raise RuntimeError(
            "This script expects mc_wgt_Flux_BeamFocus to be int[1000]. "
            "Use --n_universes 1000 unless you have verified the tuple branch "
            "has a different physical length."
        )

    binnings = {
        "ME": ME_FLUX_BINNING,
        "LE": LE_FLUX_BINNING,
    }

    # ------------------------------------------------------------
    # Accumulator dimensions:
    #
    # [neutrino species, parent category, true-E bin]
    # [neutrino species, parent category, true-E bin, universe]
    #
    # Event contributes to:
    #   parent="all"
    # and exactly one exclusive parent category.
    # ------------------------------------------------------------

    accum = {}

    for tag, edges in binnings.items():
        n_bins = len(edges) - 1

        accum[tag] = {
            "cv_sums": np.zeros(
                (len(SPECIES), len(PARENT_CATEGORIES), n_bins),
                dtype=np.float64,
            ),
            "universe_sums": np.zeros(
                (
                    len(SPECIES),
                    len(PARENT_CATEGORIES),
                    n_bins,
                    args.n_universes,
                ),
                dtype=np.float64,
            ),
            "event_counts": np.zeros(
                (len(SPECIES), len(PARENT_CATEGORIES), n_bins),
                dtype=np.int64,
            ),
            "skipped_counts": np.zeros(
                (len(SPECIES), len(PARENT_CATEGORIES), n_bins),
                dtype=np.int64,
            ),
        }

    total_files = 0
    total_entries_available = 0
    total_entries_processed = 0
    total_fv_seen = 0
    total_fv_pass = 0

    completed_filelists = []
    resolved_parent_branch = None

    # Raw immediate-parent PDG bookkeeping, after the full truth selection.
    parent_pdg_counts = {
        name: Counter()
        for name, _, _ in SPECIES
    }

    for ifilelist, filelist in enumerate(args.filelist, start=1):

        if (
            args.max_events >= 0
            and total_entries_processed >= args.max_events
        ):
            break

        chain = ROOT.TChain(args.tree)
        nfiles = add_filelist_to_chain(chain, filelist)

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

        parent_branch = find_parent_branch(
            chain,
            requested=args.parent_branch,
        )

        if resolved_parent_branch is None:
            resolved_parent_branch = parent_branch
            print("Using parent-PDG branch:", resolved_parent_branch)
        elif parent_branch != resolved_parent_branch:
            raise RuntimeError(
                "Parent branch changed between playlists: '{}' vs '{}'".format(
                    resolved_parent_branch,
                    parent_branch,
                )
            )

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
            resolved_parent_branch,
        ]

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
            chain.SetBranchStatus(branch_name, 1)

        incoming = array.array("i", [0])
        incoming_e = array.array("d", [0.0])
        current = array.array("i", [0])
        int_type = array.array("i", [0])
        target_z = array.array("i", [0])
        vtx = array.array("d", [0.0, 0.0, 0.0])
        event_cv = array.array("d", [0.0])
        beam_cv = array.array("d", [0.0])
        beam_weights_int = array.array("i", [0] * args.n_universes)
        parent_pdg = array.array("i", [0])

        chain.SetBranchAddress("mc_incoming", incoming)
        chain.SetBranchAddress("mc_incomingE", incoming_e)
        chain.SetBranchAddress("mc_current", current)
        chain.SetBranchAddress("mc_intType", int_type)
        chain.SetBranchAddress("mc_targetZ", target_z)
        chain.SetBranchAddress("mc_vtx", vtx)
        chain.SetBranchAddress("mc_cvweight_total", event_cv)
        chain.SetBranchAddress("mc_hornCurrent_cvweight", beam_cv)
        chain.SetBranchAddress("mc_wgt_Flux_BeamFocus", beam_weights_int)
        chain.SetBranchAddress(resolved_parent_branch, parent_pdg)

        playlist_n_entries = playlist_entries_available

        if args.max_events >= 0:
            remaining = args.max_events - total_entries_processed
            playlist_n_entries = min(playlist_n_entries, remaining)

        playlist_fv_seen = 0
        playlist_fv_pass = 0

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

            species_name = SPECIES[ispecies][0]
            ppdg = int(parent_pdg[0])

            parent_pdg_counts[species_name][ppdg] += 1

            parent_category = classify_parent(ppdg)
            iparent = PARENT_INDEX[parent_category]
            iall = PARENT_INDEX["all"]

            enu_gev = incoming_e[0] * 1.0e-3

            cv = event_cv[0]
            beam_cv_value = beam_cv[0]

            valid_weight = (
                math.isfinite(cv)
                and math.isfinite(beam_cv_value)
                and cv >= 1.0e-12
                and abs(beam_cv_value) >= 1.0e-12
            )

            beam_relative = None

            if valid_weight:
                beam_relative = np.empty(
                    args.n_universes,
                    dtype=np.float64,
                )

                for u in range(args.n_universes):
                    beam_absolute = beam_weights_int[u] * 1.0e-7
                    beam_relative[u] = beam_absolute / beam_cv_value

            for tag, edges in binnings.items():

                if enu_gev < edges[0] or enu_gev >= edges[-1]:
                    continue

                ibin = int(
                    np.searchsorted(
                        edges,
                        enu_gev,
                        side="right",
                    ) - 1
                )

                if ibin < 0 or ibin >= len(edges) - 1:
                    continue

                a = accum[tag]

                # Inclusive "all" plus one exclusive parent category.
                for ip in (iall, iparent):

                    a["event_counts"][ispecies, ip, ibin] += 1

                    if not valid_weight:
                        a["skipped_counts"][ispecies, ip, ibin] += 1
                        continue

                    a["cv_sums"][ispecies, ip, ibin] += cv
                    a["universe_sums"][ispecies, ip, ibin, :] += (
                        cv * beam_relative
                    )

        total_entries_processed += playlist_n_entries

        print()
        print("Finished playlist:", filelist)
        print("Playlist FV candidates:", playlist_fv_seen)
        print("Playlist FV passing:", playlist_fv_pass)
        print("Cumulative entries processed:", total_entries_processed)

        playlist_completed_fully = (
            playlist_n_entries == playlist_entries_available
        )

        if playlist_completed_fully:
            completed_filelists.append(filelist)

        if args.npz:
            checkpoint_filelists = list(completed_filelists)

            if not playlist_completed_fully:
                checkpoint_filelists.append(filelist + " [PARTIAL]")

            save_accumulator_checkpoint(
                accum=accum,
                binnings=binnings,
                output_npz=args.npz,
                completed_filelists=checkpoint_filelists,
                parent_branch=resolved_parent_branch,
                checkpoint_index=ifilelist,
            )

        if not playlist_completed_fully:
            break

    print()
    print("=" * 100)
    print("Finished event processing")
    print("=" * 100)
    print("Total files added:", total_files)
    print("Total entries available in visited playlists:", total_entries_available)
    print("Total entries processed:", total_entries_processed)
    print("Total FV candidates:", total_fv_seen)
    print("Total FV passing:", total_fv_pass)
    print("Parent-PDG branch:", resolved_parent_branch)

    print()
    print("Most common immediate parent PDGs after selection:")
    for species_name, _, _ in SPECIES:
        print()
        print(species_name)
        for pdg, count in parent_pdg_counts[species_name].most_common(15):
            print("  {:8d} : {:12d}".format(pdg, count))

    if args.npz:
        final_filelists = list(completed_filelists)

        if len(final_filelists) == 0 and len(args.filelist) > 0:
            final_filelists = [args.filelist[0] + " [PARTIAL]"]

        save_accumulator_checkpoint(
            accum=accum,
            binnings=binnings,
            output_npz=args.npz,
            completed_filelists=final_filelists,
            parent_branch=resolved_parent_branch,
            checkpoint_index=None,
        )

    # ------------------------------------------------------------
    # Convert accumulators to per-bin universe-average weights.
    #
    # Dimensions:
    #   weights/sigma/... = [species, parent, E-bin, (universe)]
    # ------------------------------------------------------------

    results = {}

    for tag, edges in binnings.items():
        a = accum[tag]
        n_bins = len(edges) - 1

        weights = np.full(
            (
                len(SPECIES),
                len(PARENT_CATEGORIES),
                n_bins,
                args.n_universes,
            ),
            np.nan,
            dtype=np.float64,
        )

        sigma_cv = np.full(
            (len(SPECIES), len(PARENT_CATEGORIES), n_bins),
            np.nan,
            dtype=np.float64,
        )

        sigma_centered = np.full_like(
            sigma_cv,
            np.nan,
        )

        mean_shift = np.full_like(
            sigma_cv,
            np.nan,
        )

        for ispecies in range(len(SPECIES)):
            for iparent in range(len(PARENT_CATEGORIES)):
                for ibin in range(n_bins):

                    cvsum = a["cv_sums"][ispecies, iparent, ibin]
                    if cvsum <= 0.0:
                        continue

                    weights[ispecies, iparent, ibin, :] = (
                        a["universe_sums"][ispecies, iparent, ibin, :]
                        / cvsum
                    )

                    w = weights[ispecies, iparent, ibin, :]
                    valid = np.isfinite(w)

                    if not np.any(valid):
                        continue

                    w = w[valid]
                    d = w - 1.0

                    mean_shift[ispecies, iparent, ibin] = np.mean(d)
                    sigma_cv[ispecies, iparent, ibin] = np.sqrt(
                        np.mean(d * d)
                    )

                    if len(w) > 1:
                        sigma_centered[ispecies, iparent, ibin] = np.std(
                            w,
                            ddof=1,
                        )

        results[tag] = {
            "weights": weights,
            "sigma_cv": sigma_cv,
            "sigma_centered": sigma_centered,
            "mean_shift": mean_shift,
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

        for ispecies, (species_name, _, _) in enumerate(SPECIES):

            for iparent, parent_name in enumerate(PARENT_CATEGORIES):

                stem = "{}_{}".format(species_name, parent_name)

                h_sigma_cv = make_hist(
                    "beamfocus_sigma_cv_{}".format(stem),
                    ";True E_{#nu} [GeV];BeamFocus RMS uncertainty from CV",
                    edges,
                )

                h_sigma_centered = make_hist(
                    "beamfocus_sigma_centered_{}".format(stem),
                    ";True E_{#nu} [GeV];BeamFocus centered fractional uncertainty",
                    edges,
                )

                h_mean_shift = make_hist(
                    "beamfocus_mean_shift_{}".format(stem),
                    ";True E_{#nu} [GeV];Mean BeamFocus shift from CV",
                    edges,
                )

                h_events = make_hist(
                    "event_count_{}".format(stem),
                    ";True E_{#nu} [GeV];Selected events",
                    edges,
                )

                for ibin in range(len(edges) - 1):
                    rb = ibin + 1

                    if np.isfinite(
                        r["sigma_cv"][ispecies, iparent, ibin]
                    ):
                        h_sigma_cv.SetBinContent(
                            rb,
                            r["sigma_cv"][ispecies, iparent, ibin],
                        )

                    if np.isfinite(
                        r["sigma_centered"][ispecies, iparent, ibin]
                    ):
                        h_sigma_centered.SetBinContent(
                            rb,
                            r["sigma_centered"][ispecies, iparent, ibin],
                        )

                    if np.isfinite(
                        r["mean_shift"][ispecies, iparent, ibin]
                    ):
                        h_mean_shift.SetBinContent(
                            rb,
                            r["mean_shift"][ispecies, iparent, ibin],
                        )

                    h_events.SetBinContent(
                        rb,
                        r["event_counts"][ispecies, iparent, ibin],
                    )

                h_sigma_cv.Write()
                h_sigma_centered.Write()
                h_mean_shift.Write()
                h_events.Write()

        fout.cd()

    fout.Close()
    print("Wrote ROOT output:", args.output)

    # ------------------------------------------------------------
    # NPZ output
    # ------------------------------------------------------------

    if args.npz:
        payload = {
            "species_names": np.array([x[0] for x in SPECIES]),
            "species_pdgs": np.array([x[1] for x in SPECIES], dtype=int),
            "parent_categories": np.array(PARENT_CATEGORIES),
            "parent_branch": np.array(resolved_parent_branch),
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

        np.savez(args.npz, **payload)
        print("Wrote NPZ output:", args.npz)

    # ------------------------------------------------------------
    # CSV output: one file per binning.
    #
    # Long format is much easier to inspect:
    #   E_low,E_high,species,parent,N,...,sigma...
    # ------------------------------------------------------------

    if args.csv:
        csv_base = args.csv[:-4] if args.csv.lower().endswith(".csv") else args.csv

        for tag, edges in binnings.items():
            csv_path = "{}_{}.csv".format(csv_base, tag)
            r = results[tag]

            with open(csv_path, "w") as csvfile:
                writer = csv.writer(csvfile)

                writer.writerow([
                    "E_low",
                    "E_high",
                    "species",
                    "parent",
                    "N",
                    "N_skipped",
                    "mean_shift",
                    "sigma_cv",
                    "sigma_centered",
                ])

                for ispecies, (species_name, _, _) in enumerate(SPECIES):
                    for iparent, parent_name in enumerate(PARENT_CATEGORIES):
                        for ibin in range(len(edges) - 1):

                            writer.writerow([
                                "{:.6f}".format(edges[ibin]),
                                "{:.6f}".format(edges[ibin + 1]),
                                species_name,
                                parent_name,
                                int(
                                    r["event_counts"][
                                        ispecies,
                                        iparent,
                                        ibin,
                                    ]
                                ),
                                int(
                                    r["skipped_counts"][
                                        ispecies,
                                        iparent,
                                        ibin,
                                    ]
                                ),
                                "{:.10e}".format(
                                    r["mean_shift"][
                                        ispecies,
                                        iparent,
                                        ibin,
                                    ]
                                ),
                                "{:.10e}".format(
                                    r["sigma_cv"][
                                        ispecies,
                                        iparent,
                                        ibin,
                                    ]
                                ),
                                "{:.10e}".format(
                                    r["sigma_centered"][
                                        ispecies,
                                        iparent,
                                        ibin,
                                    ]
                                ),
                            ])

            print("Wrote CSV output:", csv_path)

    # ------------------------------------------------------------
    # Compact text summary, useful around the focusing region.
    # ------------------------------------------------------------

    for tag, edges in binnings.items():
        r = results[tag]

        print()
        print("{}-binned parentage summary".format(tag))
        print("=" * 100)

        for ispecies, (species_name, _, _) in enumerate(SPECIES):

            print()
            print(species_name)

            for iparent, parent_name in enumerate(PARENT_CATEGORIES):

                total_n = int(
                    np.sum(
                        r["event_counts"][ispecies, iparent, :]
                    )
                )

                if total_n == 0:
                    continue

                print(
                    "  {:14s} total selected = {:12d}".format(
                        parent_name,
                        total_n,
                    )
                )


if __name__ == "__main__":
    main()

#!/usr/bin/env python

"""
  eventSelection.py:
  The executalbe to perform event selection

"""

import os
import ROOT
import string
import cProfile
import signal
import sys
from itertools import chain
import math

#start loading my modules
from tools import Utilities
from config.PlotConfig import HISTS_TO_MAKE
from config.AnalysisConfig import AnalysisConfig
from tools.SystematicsUniverse import GetAllSystematicsUniverses
from tools.SystematicsUniverse import get_flux_ratio_me_to_le
from tools.EventClassification import EventClassifier
from tools.KinematicsCalculator import KinematicsCalculator
from tools.MyHistograms import MakePlotProcessors
from tools import TruthTools
from config import SystematicsConfig

#def timeout_handler(signum,frame):
#    print "some steps takes forever to finish, I am not going to wait."
#    raise Exception("time out")

ROOT.TH1.AddDirectory(False)




def plotRecoKin(mc, chainwrapper, outfile):
    """The main code of event selection."""

    from collections import Counter

    kin_cal = KinematicsCalculator(
        correct_beam_angle=True,
        correct_MC_energy_scale=False,
        calc_true=mc,
        is_pc=AnalysisConfig.is_pc
    )

    eventClassifier = EventClassifier(
        classifiers=["Reco", "Truth"] if mc else ["Reco"],
        use_kin_cuts=True,
        use_sideband=AnalysisConfig.sidebands
    )

    universes = GetAllSystematicsUniverses(
        chainwrapper,
        not mc,
        AnalysisConfig.is_pc,
        AnalysisConfig.exclude_universes
    )

    for univ in chain.from_iterable(iter(universes.values())):
        univ.LoadTools(kin_cal, eventClassifier)

    Plots = preparePlots(universes, mc)

    nEvents = chainwrapper.GetEntries()

    print(f"Total number of events RECO: {nEvents}")

    if AnalysisConfig.testing and nEvents > 2000:
        nEvents = 2000

    print("plotRecoKin, mc", mc)

    # --------------------------------------------------
    # Initialize CV 2p2h diagnostic counters
    # --------------------------------------------------
    interaction_counts = Counter()
    target_counts = Counter()

    if mc:
        from tools.SystematicsUniverse import CVPythonUniverse

        if hasattr(CVPythonUniverse, "_debug_genie_counts"):
            delattr(CVPythonUniverse, "_debug_genie_counts")

        # Reset the per-sample 2p2h diagnostic counters.
        if hasattr(CVPythonUniverse, "_debug_2p2h_counts"):
            delattr(CVPythonUniverse, "_debug_2p2h_counts")

        print(
            f"[2P2H AUDIT] Starting CV diagnostics: "
            f"{nEvents} events",
            flush=True
        )

    setAlarm = AnalysisConfig.grid

    # --------------------------------------------------
    # Main event loop
    # --------------------------------------------------
    for counter in range(nEvents):

        if counter % 10000 == 0:
            print(counter, flush=True)

            if setAlarm:
                signal.alarm(900)

        flux_debug_this_event = 0
        selected_event_for_debug = False

        for universe in chain.from_iterable(iter(universes.values())):

            universe.SetEntry(counter)

            # ------------------------------------------
            # CV-only 2p2h diagnostics
            # ------------------------------------------
            if mc and universe.ShortName() == "cv":

                int_type = int(universe.mc_intType)
                target_z = int(universe.mc_targetZ)

                interaction_counts[int_type] += 1
                target_counts[(int_type, target_z)] += 1

                # Print a few raw MEC examples.
                if (
                    int_type == 8
                    and interaction_counts[8] <= 10
                ):
                    print(
                        "[2P2H MEC EVENT]",
                        f"entry={counter}",
                        f"mc_intType={int_type}",
                        f"mc_targetZ={target_z}",
                        f"mc_Q2={universe.mc_Q2}",
                        f"mc_w={universe.mc_w}",
                        flush=True
                    )

                # Compare MAT with the existing Q2/W guard.
                universe.Debug2p2hCV()
                universe.DebugGenieCV()

            universe.__dict__.pop(
                "_cached_new_true_l_over_e",
                None
            )

            universe.ResetWeight()

            if (
                mc
                and AnalysisConfig.skip_2p2h
                and universe.mc_intType == 8
            ):
                continue

            # Recalculate kinematics and classification
            # for non-vertical universes.
            if not universe.IsVerticalOnly():
                kin_cal.CalculateKinematics(universe)
                eventClassifier.Classify(universe)

            if (
                eventClassifier.side_band is not None
                or eventClassifier.is_true_signal
            ):
                for entry in Plots:
                    entry.Process(universe)

    # --------------------------------------------------
    # Final diagnostic summary
    # --------------------------------------------------
    if mc:
        counts = getattr(
            CVPythonUniverse,
            "_debug_2p2h_counts",
            {}
        )

        print("\n" + "=" * 55, flush=True)
        print("CV LOW-RECOIL 2P2H AUDIT", flush=True)
        print("=" * 55, flush=True)

        print(f"Total MC entries processed: {nEvents}")

        print("\nInteraction type counts:")
        for int_type, count in sorted(interaction_counts.items()):
            print(
                f"  mc_intType={int_type:3d}: {count}",
                flush=True
            )

        print("\nInteraction type / target-Z counts:")
        for (int_type, target_z), count in sorted(
            target_counts.items()
        ):
            print(
                f"  mc_intType={int_type:3d}, "
                f"mc_targetZ={target_z:3d}: {count}",
                flush=True
            )

        total_mec = interaction_counts.get(8, 0)

        eligible_mec = sum(
            count
            for (int_type, target_z), count
            in target_counts.items()
            if int_type == 8 and target_z >= 2
        )

        print("\nMEC summary:")
        print(f"  All MEC events:             {total_mec}")
        print(f"  MEC on nuclear targets:     {eligible_mec}")
        print(
            f"  MEC on hydrogen/invalid Z: {total_mec - eligible_mec}"
        )

        print("\nMAT 2p2h weight diagnostic:")

        if not counts:
            print("  No eligible MEC events reached the diagnostic.")
        else:
            for key, value in sorted(counts.items()):
                print(f"  {key:24s}: {value}")

            if counts["total_mec"] > 0:
                fraction = (
                    counts["changed_weight"]
                    / counts["total_mec"]
                )

                print(
                    f"\n  Fraction with changed weights: "
                    f"{fraction:.6%}"
                )

            print(
                "  Difference in sums of 2p2h factors:",
                counts["sum_mat"] - counts["sum_current"]
            )

        print("=" * 55 + "\n", flush=True)

        print("\n" + "=" * 55)
        print("CV GENIE REWEIGHT AUDIT")
        print("=" * 55)

        genie_counts = getattr(
            CVPythonUniverse, "_debug_genie_counts", {}
        )

        for key, value in sorted(genie_counts.items()):
            print(f"{key:26s}: {value}")

        if genie_counts and genie_counts["total"] > 0:
            n = genie_counts["total"]

            print("\nMean GENIE factors:")
            print(
                "V1 NonResPi only:",
                genie_counts["sum_v1_nonrespi"] / n
            )
            print(
                "V1 + Deuterium:",
                genie_counts["sum_v1_deuterium"] / n
            )
            print(
                "Actual MAT GENIE:",
                genie_counts["sum_mat"] / n
            )

        print("=" * 55, flush=True)


    signal.alarm(0)

    outfile.cd()

    for entry in Plots:
        entry.Finalize()



def plotTruthKin(chainwrapper,outfile):
    kin_cal = KinematicsCalculator(correct_beam_angle=True, correct_MC_energy_scale=False, calc_true = True, calc_reco = False)
    eventClassifier = EventClassifier(classifiers=["Truth"],use_kin_cuts=True, use_sideband=[])
    universes = GetAllSystematicsUniverses(chainwrapper, False)
    for univ in chain.from_iterable(iter(universes.values())):
        univ.LoadTools(kin_cal,eventClassifier)
    nEvents = chainwrapper.GetEntries()
    print(f"Total number of events TRUTH: ", {nEvents})
    Plots = prepareTruthPlots(universes)
    print("plotTruthKin")
    if AnalysisConfig.testing and nEvents > 1000:
        nEvents = 1000
    for counter in range(nEvents):
        #half an hour for a event, should be much more than needed unless stuck in I/O
        if counter %100000 == 0:
            print(counter)

        for universe in chain.from_iterable(iter(universes.values())):
            universe.SetEntry(counter)
            universe.ResetWeight()

            #only update kin_cal & eventClassifier when universe in not vertical only.
            if not universe.IsVerticalOnly():
                kin_cal.CalculateKinematics(universe)
                
                reco_before = eventClassifier.counter[0]
                true_before = eventClassifier.counter[1]
                eventClassifier.Classify(universe)

            if eventClassifier.is_true_signal:
                for entry in Plots:
                    entry.Process(universe)

    outfile.cd()
    for entry in Plots:
        entry.Finalize()

def preparePlots(universes,mc):
    # make a bunch of Plot Processor, grouped by signal/sideband
    plots=set([])

    for entry in HISTS_TO_MAKE:
        if (isinstance(entry,str) and entry.startswith("True Signal")):
            continue
        settings = {"key":entry,"region":AnalysisConfig.sidebands,"mc":mc}
        plots.update(MakePlotProcessors(**settings))

    # add the errorband map to plot processor.
    for plot in plots:
        plot.AddErrorBands(universes)

    return plots


def prepareTruthPlots(universes):
    plots = []
    for entry in HISTS_TO_MAKE:
        if not (isinstance(entry, str) and entry.startswith("True Signal")):
            continue
        # Pass plot name string only: MakePlotProcessors calls TranslateSettings internally.
        # True Signal plots must use truth_signal_tags (includes ignore_selection) in PlotLibrary.
        settings = {"key": entry, "region": ["Signal"], "mc": True}
        plots.extend(MakePlotProcessors(**settings))

    for entry in plots:
        entry.AddErrorBands(universes)

    return plots

def CopyMetaTreeToOutPutFile(outfile):
    metatree = Utilities.fileChain(AnalysisConfig.playlist,st,AnalysisConfig.ntuple_tag,"Meta",AnalysisConfig.count[0],AnalysisConfig.count[1])
    if AnalysisConfig.is_pc: #PC samples are on top of existing sample. Dont count more POT
        return None
    raw_metatree=metatree.GetChain()
    raw_metatree.SetBranchStatus("*",0)
    for _ in ["POT_Used"]:
        raw_metatree.SetBranchStatus(_,1)
    outfile.cd()
    copiedtree = raw_metatree.CopyTree("")
    del metatree
    copiedtree.Write()


if __name__ == "__main__":

    Reco = AnalysisConfig.run_reco
    Truth = AnalysisConfig.truth
    POT_cal = AnalysisConfig.POT_cal
    print("playlist %s running ---------" % AnalysisConfig.playlist)
    for st in AnalysisConfig.data_types:
        print(st)
        outputSelectionHistogram = AnalysisConfig.SelectionHistoPath(AnalysisConfig.playlist,"data" in st,True)
        output_file = ROOT.TFile.Open(outputSelectionHistogram,"RECREATE")
        CopyMetaTreeToOutPutFile(output_file)

        if st=="mc" and Truth:
            print("selecting truth")
            plotTruthKin(Utilities.fileChain(AnalysisConfig.playlist,"mc",AnalysisConfig.ntuple_tag,"Truth",AnalysisConfig.count[0],AnalysisConfig.count[1]),output_file)
        if Reco :
            print("selecting reco")
            plotRecoKin(st=="mc", Utilities.fileChain(AnalysisConfig.playlist,st,AnalysisConfig.ntuple_tag,None,AnalysisConfig.count[0],AnalysisConfig.count[1]), output_file)
        output_file.Close()
        print("selection is done for ", st, AnalysisConfig.playlist)
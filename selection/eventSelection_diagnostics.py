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
    """Main event selection, with final EM-scale diagnostic."""

    kin_cal = KinematicsCalculator(
        correct_beam_angle=True,
        correct_MC_energy_scale=False,
        calc_true=mc,
        is_pc=AnalysisConfig.is_pc,
    )

    eventClassifier = EventClassifier(
        classifiers=["Reco", "Truth"] if mc else ["Reco"],
        use_kin_cuts=True,
        use_sideband=AnalysisConfig.sidebands,
    )

    universes = GetAllSystematicsUniverses(
        chainwrapper,
        not mc,
        AnalysisConfig.is_pc,
        AnalysisConfig.exclude_universes,
    )

    # Keep only CV and regional electron-energy systematics
    # for this diagnostic.
    if mc and abs(SystematicsConfig.AnaNuPDG) == 12:
        required = ("cv", "elE_ECAL", "elE_HCAL")

        assert all(name in universes for name in required), (
            f"Missing required universes. Available: {list(universes)}"
        )

        universes = {
            name: universes[name]
            for name in required
        }

        print(
            "[EM AUDIT] Active universes:",
            {name: len(univs) for name, univs in universes.items()}
        )

    if not mc:
        universes = {"cv": universes["cv"]}

        print(
            "[DATA EM AUDIT] Active universes:",
            {name: len(univs) for name, univs in universes.items()}
        )

    for univ in chain.from_iterable(universes.values()):
        univ.LoadTools(kin_cal, eventClassifier)

    Plots = preparePlots(universes, mc)

    nEvents = chainwrapper.GetEntries()
    print(f"Total number of events RECO: {nEvents}")

    if AnalysisConfig.testing and nEvents > 5000:
        nEvents = 5000

    print("plotRecoKin, mc", mc)

    setAlarm = AnalysisConfig.grid

    # ====================================================
    # Final regional electron-energy diagnostic
    # ====================================================
    audit_em = (
        mc and abs(SystematicsConfig.AnaNuPDG) == 12
    )

    regions = ("ECAL", "HCAL")
    regional_names = {
        "elE_ECAL": "ECAL",
        "elE_HCAL": "HCAL",
    }

    stats = {
        "total": 0,
        "invalid_ecal_all": 0,
        "invalid_hcal_all": 0,
        "invalid_either_all": 0,
        "cv_selected": 0,
        "cv_selected_invalid_ecal": 0,
        "cv_selected_invalid_hcal": 0,
        "cv_selected_invalid_either": 0,
    }

    reg_stats = {
        (region, sigma): {
            "selected": 0,
            "selected_invalid": 0,
            "selected_invalid_either": 0,
            "selected_not_cv": 0,
            "selected_not_cv_invalid": 0,
            "cv_selected_not_region": 0,
            "cv_selected_not_region_invalid": 0,
        }
        for region in regions
        for sigma in (-1, 1)
    }

    sample_messages = []
    max_samples = 12

    def invalid_energy(value):
        return not math.isfinite(value) or value < 0

    # ====================================================
    # Event loop
    # ====================================================
    for counter in range(nEvents):

        if counter % 10000 == 0:
            print(counter)
            if setAlarm:
                signal.alarm(900)

        # Read CV regional energies once per event.
        if audit_em:
            cv = universes["cv"][0]
            cv.SetEntry(counter)

            e_ecal = cv.GetVecElem(
                "prong_ECALCalibE", 0
            )
            e_hcal = cv.GetVecElem(
                "prong_HCALCalibE", 0
            )

            bad_ecal = invalid_energy(e_ecal)
            bad_hcal = invalid_energy(e_hcal)
            bad_either = bad_ecal or bad_hcal

            stats["total"] += 1
            stats["invalid_ecal_all"] += int(bad_ecal)
            stats["invalid_hcal_all"] += int(bad_hcal)
            stats["invalid_either_all"] += int(bad_either)

        # Per-event selection decisions.
        cv_selected = False
        regional_results = {}

        # ====================================================
        # Systematic universe loop
        # ====================================================
        for universe in chain.from_iterable(
            universes.values()
        ):
            universe.SetEntry(counter)

            # Original event-selection workflow.
            universe.__dict__.pop(
                "_cached_new_true_l_over_e", None
            )
            universe.ResetWeight()

            if (
                mc
                and AnalysisConfig.skip_2p2h
                and universe.mc_intType == 8
            ):
                continue

            if not universe.IsVerticalOnly():
                kin_cal.CalculateKinematics(universe)
                eventClassifier.Classify(universe)

                if not mc and universe is universes["cv"][0]:
                    selected = (
                        eventClassifier.side_band is not None
                        or eventClassifier.is_true_signal
                    )

                    if selected:
                        e_ecal = universe.GetVecElem(
                            "prong_ECALCalibE", 0
                        )

                        if not math.isfinite(e_ecal) or e_ecal < 0:
                            print(
                                "[DATA SELECTED INVALID ECAL]",
                                f"event={counter}",
                                f"ECALCalibE={e_ecal}",
                                f"EMShift={universe.GetEMEnergyShift()}",
                            )

                selected = (
                    eventClassifier.side_band is not None
                    or eventClassifier.is_true_signal
                )

                if audit_em:

                    # CV selection decision.
                    if universe is universes["cv"][0]:
                        cv_selected = selected

                        if selected:
                            stats["cv_selected"] += 1

                            if bad_ecal:
                                stats[
                                    "cv_selected_invalid_ecal"
                                ] += 1

                            if bad_hcal:
                                stats[
                                    "cv_selected_invalid_hcal"
                                ] += 1

                            if bad_either:
                                stats[
                                    "cv_selected_invalid_either"
                                ] += 1

                    # Regional systematic decisions.
                    name = universe.ShortName()

                    if name in regional_names:
                        region = regional_names[name]
                        sigma = int(universe.GetSigma())

                        key = (region, sigma)

                        calib_e = (
                            e_ecal
                            if region == "ECAL"
                            else e_hcal
                        )

                        bad_region = invalid_energy(
                            calib_e
                        )

                        regional_results[key] = (
                            selected,
                            bad_region,
                        )

                        if selected:
                            reg_stats[key]["selected"] += 1

                            if bad_region:
                                reg_stats[key][
                                    "selected_invalid"
                                ] += 1

                            if bad_either:
                                reg_stats[key][
                                    "selected_invalid_either"
                                ] += 1

            # Original histogram processing.
            if (
                eventClassifier.side_band is not None
                or eventClassifier.is_true_signal
            ):
                for entry in Plots:
                    entry.Process(universe)

        # ====================================================
        # Compare regional selection against CV
        # ====================================================
        if audit_em:
            for key, result in regional_results.items():
                selected, bad_region = result
                region, sigma = key

                if selected and not cv_selected:
                    reg_stats[key][
                        "selected_not_cv"
                    ] += 1

                    if bad_region:
                        reg_stats[key][
                            "selected_not_cv_invalid"
                        ] += 1

                elif cv_selected and not selected:
                    reg_stats[key][
                        "cv_selected_not_region"
                    ] += 1

                    if bad_region:
                        reg_stats[key][
                            "cv_selected_not_region_invalid"
                        ] += 1

                # Keep a few concrete examples.
                if (
                    bad_region
                    and (selected or cv_selected)
                    and len(sample_messages) < max_samples
                ):
                    sample_messages.append(
                        (
                            counter,
                            region,
                            sigma,
                            e_ecal,
                            e_hcal,
                            cv_selected,
                            selected,
                        )
                    )

    signal.alarm(0)

    # ====================================================
    # Final diagnostic report
    # ====================================================
    if audit_em:

        print("\n" + "=" * 65)
        print("FINAL ELECTRON REGIONAL ENERGY-SCALE AUDIT")
        print("=" * 65)

        print(f"Events processed: {stats['total']}")

        print("\n--- All reconstructed entries ---")
        print(
            "Invalid ECAL:",
            stats["invalid_ecal_all"]
        )
        print(
            "Invalid HCAL:",
            stats["invalid_hcal_all"]
        )
        print(
            "Invalid either:",
            stats["invalid_either_all"]
        )

        print("\n--- CV selection ---")
        print(
            "CV selected:",
            stats["cv_selected"]
        )
        print(
            "Selected with invalid ECAL:",
            stats["cv_selected_invalid_ecal"]
        )
        print(
            "Selected with invalid HCAL:",
            stats["cv_selected_invalid_hcal"]
        )
        print(
            "Selected with invalid either:",
            stats["cv_selected_invalid_either"]
        )

        print("\n--- Regional systematic universes ---")

        for key, s in sorted(reg_stats.items()):
            region, sigma = key

            print(
                f"\n{region} sigma={sigma:+d}"
            )
            print(
                "  Selected:",
                s["selected"]
            )
            print(
                "  Selected with invalid region energy:",
                s["selected_invalid"]
            )
            print(
                "  Selected with invalid either energy:",
                s["selected_invalid_either"]
            )
            print(
                "  Entered selection relative to CV:",
                s["selected_not_cv"]
            )
            print(
                "    With invalid region energy:",
                s["selected_not_cv_invalid"]
            )
            print(
                "  Left selection relative to CV:",
                s["cv_selected_not_region"]
            )
            print(
                "    With invalid region energy:",
                s["cv_selected_not_region_invalid"]
            )

        print("\n--- Example invalid-energy events ---")

        if not sample_messages:
            print(
                "No invalid regional-energy events "
                "selected in CV or regional universes."
            )
        else:
            for (
                event,
                region,
                sigma,
                e_ecal,
                e_hcal,
                cv_sel,
                region_sel,
            ) in sample_messages:
                print(
                    f"event={event}",
                    f"region={region}",
                    f"sigma={sigma:+d}",
                    f"ECAL={e_ecal:.3f}",
                    f"HCAL={e_hcal:.3f}",
                    f"CV_selected={cv_sel}",
                    f"region_selected={region_sel}",
                )

        print("=" * 65 + "\n")

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
        if Reco :
            print("selecting reco")
            plotRecoKin(st=="mc", Utilities.fileChain(AnalysisConfig.playlist,st,AnalysisConfig.ntuple_tag,None,AnalysisConfig.count[0],AnalysisConfig.count[1]), output_file)
            print("done selection reco")
        if st=="mc" and Truth:
            print("selecting truth")
            plotTruthKin(Utilities.fileChain(AnalysisConfig.playlist,"mc",AnalysisConfig.ntuple_tag,"Truth",AnalysisConfig.count[0],AnalysisConfig.count[1]),output_file)
        output_file.Close()
        print("selection is done for ", st, AnalysisConfig.playlist)
        # if st=="mc" and Truth:
        #     print("selecting truth")
        #     plotTruthKin(Utilities.fileChain(AnalysisConfig.playlist,"mc",AnalysisConfig.ntuple_tag,"Truth",AnalysisConfig.count[0],AnalysisConfig.count[1]),output_file)
        # if Reco :
        #     print("selecting reco")
        #     plotRecoKin(st=="mc", Utilities.fileChain(AnalysisConfig.playlist,st,AnalysisConfig.ntuple_tag,None,AnalysisConfig.count[0],AnalysisConfig.count[1]), output_file)
        # output_file.Close()
        # print("selection is done for ", st, AnalysisConfig.playlist)
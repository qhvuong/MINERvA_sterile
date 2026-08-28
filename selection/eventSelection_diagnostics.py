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

#def timeout_handler(signum,frame):
#    print "some steps takes forever to finish, I am not going to wait."
#    raise Exception("time out")

ROOT.TH1.AddDirectory(False)

def plotRecoKin(mc, chainwrapper, outfile):
    """ The main code of event selection """

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

    # ============================================================
    # DIAGNOSTIC SETUP ONLY
    # ============================================================

    cv_universe = universes.get("cv", [None])[0]
    flux_universe0 = universes.get("Flux", [None])[0]

    debug_same_entry_count = 0
    DEBUG_SAME_ENTRY_MAX = 5

    # ============================================================
    # ORIGINAL CODE
    # ============================================================

    for univ in chain.from_iterable(iter(universes.values())):
        univ.LoadTools(kin_cal, eventClassifier)

    Plots = preparePlots(universes, mc)
    nEvents = chainwrapper.GetEntries()
    debug_prints = 0

    print(f"Total number of events RECO: ", {nEvents})

    if AnalysisConfig.testing and nEvents > 1500:
        nEvents = 1500

    print("plotRecoKin, mc ", mc)

    setAlarm = AnalysisConfig.grid

    for counter in range(nEvents):

        # 1/4 hour for 10k event, should be more than needed unless stuck in I/O
        if counter % 10000 == 0:
            print(counter)

            if setAlarm:
                signal.alarm(900)

        flux_debug_this_event = 0
        selected_event_for_debug = False

        # ========================================================
        # DIAGNOSTIC ONLY:
        # Compare CV and Flux universe 0 on the exact same
        # true CCNuEQE event.
        # ========================================================

        if (
            mc
            and cv_universe is not None
            and flux_universe0 is not None
            and debug_same_entry_count < DEBUG_SAME_ENTRY_MAX
        ):

            # ----------------------------------------------------
            # CV
            # ----------------------------------------------------

            cv_universe.SetEntry(counter)

            cv_universe.__dict__.pop(
                "_cached_new_true_l_over_e",
                None
            )

            cv_universe.ResetWeight()

            cv_pdg = int(cv_universe.mc_incoming)
            cv_enu = float(cv_universe.mc_incomingE) * 1.0e-3
            cv_inttype = int(cv_universe.mc_intType)

            cv_run = (
                int(cv_universe.mc_run)
                if hasattr(cv_universe, "mc_run")
                else -999
            )

            cv_subrun = (
                int(cv_universe.mc_subrun)
                if hasattr(cv_universe, "mc_subrun")
                else -999
            )

            cv_evt = (
                int(cv_universe.mc_nthEvtInFile)
                if hasattr(cv_universe, "mc_nthEvtInFile")
                else -999
            )

            # Recalculate only for diagnostic classification
            kin_cal.CalculateKinematics(cv_universe)
            eventClassifier.Classify(cv_universe)

            cv_truth_class = eventClassifier.truth_class

            if cv_truth_class == "CCNuEQE":

                cv_en4 = float(
                    kin_cal.reco_E_lep
                    + kin_cal.reco_visE
                )

                cv_weight = float(
                    cv_universe.GetWeight(False)
                )

                # Print individual pieces making up CV GetWeight()
                cv_universe.DebugStandardWeightPieces(
                    "CV entry={}".format(counter)
                )

                # ------------------------------------------------
                # Flux universe 0
                # ------------------------------------------------

                flux_universe0.SetEntry(counter)

                flux_universe0.__dict__.pop(
                    "_cached_new_true_l_over_e",
                    None
                )

                flux_universe0.ResetWeight()

                flux_pdg = int(flux_universe0.mc_incoming)
                flux_enu = (
                    float(flux_universe0.mc_incomingE)
                    * 1.0e-3
                )

                flux_inttype = int(
                    flux_universe0.mc_intType
                )

                flux_run = (
                    int(flux_universe0.mc_run)
                    if hasattr(flux_universe0, "mc_run")
                    else -999
                )

                flux_subrun = (
                    int(flux_universe0.mc_subrun)
                    if hasattr(flux_universe0, "mc_subrun")
                    else -999
                )

                flux_evt = (
                    int(flux_universe0.mc_nthEvtInFile)
                    if hasattr(flux_universe0, "mc_nthEvtInFile")
                    else -999
                )

                # Explicit recalc for diagnostic comparison
                kin_cal.CalculateKinematics(flux_universe0)
                eventClassifier.Classify(flux_universe0)

                flux_truth_class = eventClassifier.truth_class

                flux_en4 = float(
                    kin_cal.reco_E_lep
                    + kin_cal.reco_visE
                )

                flux_weight = float(
                    flux_universe0.GetWeight(False)
                )

                # Print individual pieces making up Flux0 GetWeight()
                flux_universe0.DebugStandardWeightPieces(
                    "FLUX0 entry={}".format(counter)
                )

                weight_ratio = (
                    flux_weight / cv_weight
                    if cv_weight != 0.0
                    else float("nan")
                )

                print(
                    "[SAME_ENTRY_CHECK]"
                    " counter={}"
                    " cvRun={}"
                    " fluxRun={}"
                    " cvSubrun={}"
                    " fluxSubrun={}"
                    " cvEvt={}"
                    " fluxEvt={}"
                    " cvPDG={}"
                    " fluxPDG={}"
                    " cvIntType={}"
                    " fluxIntType={}"
                    " cvEnu={:.6f}"
                    " fluxEnu={:.6f}"
                    " cvEN4={:.6f}"
                    " fluxEN4={:.6f}"
                    " cvClass={}"
                    " fluxClass={}"
                    " cvW={:.12g}"
                    " fluxW={:.12g}"
                    " fluxW/cvW={:.12g}".format(
                        counter,
                        cv_run,
                        flux_run,
                        cv_subrun,
                        flux_subrun,
                        cv_evt,
                        flux_evt,
                        cv_pdg,
                        flux_pdg,
                        cv_inttype,
                        flux_inttype,
                        cv_enu,
                        flux_enu,
                        cv_en4,
                        flux_en4,
                        cv_truth_class,
                        flux_truth_class,
                        cv_weight,
                        flux_weight,
                        weight_ratio
                    )
                )

                debug_same_entry_count += 1

        # ========================================================
        # ORIGINAL UNIVERSE LOOP -- UNCHANGED
        # ========================================================

        for universe in chain.from_iterable(iter(universes.values())):

            universe.SetEntry(counter)

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

            # only update kin_cal & eventClassifier when universe
            # is not vertical only.
            if not universe.IsVerticalOnly():
                kin_cal.CalculateKinematics(universe)
                eventClassifier.Classify(universe)

            if (
                eventClassifier.side_band is not None
                or eventClassifier.is_true_signal
            ):
                for entry in Plots:
                    entry.Process(universe)

    signal.alarm(0)

    outfile.cd()

    for entry in Plots:
        # print entry.histwrapper.name
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
# def prepareTruthPlots(universes):
#     plots=[]
#     for entry in HISTS_TO_MAKE:
#         if not (isinstance(entry,str) and entry.startswith("True Signal")):
#             continue
#         settings = {"key":entry,"region":"Signal","mc":True}
#         plots.extend(MakePlotProcessors(**settings))

#     for entry in plots:
#         entry.AddErrorBands(universes)

#     return plots

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
        # print("DEBUG: about to open output file")
        output_file = ROOT.TFile.Open(outputSelectionHistogram,"RECREATE")

        # print("DEBUG outputSelectionHistogram =", outputSelectionHistogram)
        # print("DEBUG AnalysisConfig.count =", AnalysisConfig.count)
        
        # print("DEBUG: about to copy meta tree")
        CopyMetaTreeToOutPutFile(output_file)
        # CopyMetaTreeToOutPutFile(output_file, st)

        if st=="mc" and Truth:
            print("selecting truth")
            plotTruthKin(Utilities.fileChain(AnalysisConfig.playlist,"mc",AnalysisConfig.ntuple_tag,"Truth",AnalysisConfig.count[0],AnalysisConfig.count[1]),output_file)
            # print("DEBUG: finished truth loop")
        if Reco :
            print("selecting reco")
            plotRecoKin(st=="mc", Utilities.fileChain(AnalysisConfig.playlist,st,AnalysisConfig.ntuple_tag,None,AnalysisConfig.count[0],AnalysisConfig.count[1]), output_file)
            # print("DEBUG: finished reco loop")
        output_file.Close()
        print("selection is done for ", st, AnalysisConfig.playlist)
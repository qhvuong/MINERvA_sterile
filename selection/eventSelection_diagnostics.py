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



# Replace plotRecoKin in selection/eventSelection.py with this function.
# Uses the imports already present in eventSelection.py; extra imports are local.
def plotRecoKin(mc, chainwrapper, outfile):
    """Original selection flow, with bounded CV swapped-flux diagnostics."""
    from config import SystematicsConfig
    from tools import SystematicsUniverse as syst

    debug_limit = int(os.environ.get("SWAP_FLUX_DEBUG_N", "20"))
    if debug_limit < 0:
        raise ValueError("SWAP_FLUX_DEBUG_N must be nonnegative")
    swapped = bool(getattr(SystematicsConfig, "USE_SWAPPED", False))
    debug_enabled = mc and swapped and not AnalysisConfig.is_pc and debug_limit > 0
    debug_prints = 0
    bad_calls = 0
    bad_weights = 0

    print("[SWAP DEBUG] config={} USE_SWAPPED={} MC={} is_pc={} max_CV_events={}".format(
        SystematicsConfig.__file__, swapped, mc, AnalysisConfig.is_pc, debug_limit), flush=True)
    if debug_enabled and not callable(getattr(syst, "get_swapped_flux_ratio", None)):
        raise RuntimeError(
            "tools.SystematicsUniverse.get_swapped_flux_ratio is missing. "
            "Import the helper there, or define it there, before running this diagnostic.")

    kin_cal = KinematicsCalculator(
        correct_beam_angle=True, correct_MC_energy_scale=False,
        calc_true=mc, is_pc=AnalysisConfig.is_pc)
    eventClassifier = EventClassifier(
        classifiers=["Reco", "Truth"] if mc else ["Reco"],
        use_kin_cuts=True, use_sideband=AnalysisConfig.sidebands)
    universes = GetAllSystematicsUniverses(
        chainwrapper, not mc, AnalysisConfig.is_pc, AnalysisConfig.exclude_universes)
    for univ in chain.from_iterable(iter(universes.values())):
        univ.LoadTools(kin_cal, eventClassifier)

    Plots = preparePlots(universes, mc)
    nEvents = chainwrapper.GetEntries()
    print("Total number of events RECO:", nEvents)
    if AnalysisConfig.testing and nEvents > 5000:
        nEvents = 5000
    print("plotRecoKin, mc", mc)
    setAlarm = AnalysisConfig.grid
    for counter in range(nEvents):
        if counter % 10000 == 0:
            print(counter)
            if setAlarm:
                signal.alarm(900)
        for universe in chain.from_iterable(iter(universes.values())):
            universe.SetEntry(counter)
            universe.__dict__.pop("_cached_new_true_l_over_e", None)
            universe.ResetWeight()
            if mc and AnalysisConfig.skip_2p2h and universe.mc_intType == 8:
                continue
            if not universe.IsVerticalOnly():
                kin_cal.CalculateKinematics(universe)
                eventClassifier.Classify(universe)

            if eventClassifier.side_band is not None or eventClassifier.is_true_signal:
                if (debug_enabled and debug_prints < debug_limit
                        and str(universe.ShortName()).lower() == "cv"):
                    source = int(universe.mc_incoming)
                    target = {12: 14, -12: -14, 14: 12, -14: -12}.get(source)
                    if target is None:
                        raise ValueError("Unexpected incoming neutrino PDG: {}".format(source))
                    energy = float(universe.mc_incomingE) * 1e-3
                    ratio_helper = syst.get_swapped_flux_ratio
                    expected_ratio = float(ratio_helper(universe))
                    cv_flux = float(universe.GetFluxAndCVWeight(energy, target))
                    calls = []

                    # Observe calls from the ACTUAL weight implementation.
                    # This wrapper returns the helper's original value unchanged.
                    def traced_ratio(*args, **kwargs):
                        value = ratio_helper(*args, **kwargs)
                        calls.append(float(value))
                        return value

                    syst.get_swapped_flux_ratio = traced_ratio
                    try:
                        total_weight = float(universe.GetWeight())
                    finally:
                        syst.get_swapped_flux_ratio = ratio_helper

                    expected_calls = 1
                    calls_ok = (len(calls) == expected_calls and
                                math.isclose(calls[0], expected_ratio,
                                             rel_tol=1e-12, abs_tol=0.0))
                    finite = math.isfinite(total_weight) and math.isfinite(cv_flux)
                    bad_calls += int(not calls_ok)
                    bad_weights += int(not finite)
                    debug_prints += 1
                    print("\n[SWAP DEBUG {}/{}] entry={} playlist={} trueE={:.8g} GeV".format(
                        debug_prints, debug_limit, counter,
                        universe.GetPlaylist(), energy), flush=True)
                    print("  source PDG={} -> flux PDG={} sideband={} truth_class={}".format(
                        source, target, eventClassifier.side_band,
                        getattr(eventClassifier, "truth_class", None)))
                    print("  extra generated target/source ratio = {:.12g}".format(expected_ratio))
                    print("  mapped-PDG CV flux weight           = {:.12g}".format(cv_flux))
                    print("  expected combined flux factor       = {:.12g}".format(
                        cv_flux * expected_ratio))
                    print("  helper calls during actual GetWeight = {} values={} [{}]".format(
                        len(calls), calls, "PASS" if calls_ok else "CHECK IMPLEMENTATION"))
                    print("  actual total event weight           = {:.12g} [{}]".format(
                        total_weight, "FINITE" if finite else "INVALID"), flush=True)
                    # GetWeight caches the actual weight just computed. Plot
                    # processors use that same result; do not multiply again here.

                for entry in Plots:
                    entry.Process(universe)

    signal.alarm(0)
    if debug_enabled:
        print("\n[SWAP DEBUG SUMMARY] inspected_CV_events={} helper_call_failures={} "
              "nonfinite_weight_events={}".format(debug_prints, bad_calls, bad_weights), flush=True)
        if debug_prints == 0:
            print("[SWAP DEBUG] No eligible CV events reached the diagnostic; "
                  "this run has not validated the implementation.", flush=True)
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
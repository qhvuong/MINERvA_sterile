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
from collections import defaultdict
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

def _new_wstat():
    return defaultdict(lambda: {
        "n": 0,
        "sumw": 0.0,
        "sumw2": 0.0,
        "maxw": 0.0,
        "max_entry": -1,
    })


def _fill_wstat(stats, energy_gev, weight, entry, bin_width=0.5):
    if energy_gev is None:
        return

    if not math.isfinite(energy_gev):
        return

    if not math.isfinite(weight):
        print("[WSTAT_BAD_WEIGHT] entry={} E={} weight={}".format(entry, energy_gev, weight))
        return

    ibin = int(energy_gev / bin_width)

    s = stats[ibin]
    s["n"] += 1
    s["sumw"] += weight
    s["sumw2"] += weight * weight

    if abs(weight) > abs(s["maxw"]):
        s["maxw"] = weight
        s["max_entry"] = entry


def _print_wstat(stats, label, bin_width=0.5):
    print("[WSTAT] ===== {} =====".format(label))

    for ibin, s in sorted(stats.items()):
        lo = ibin * bin_width
        hi = lo + bin_width

        sumw = s["sumw"]
        sumw2 = s["sumw2"]
        err = math.sqrt(sumw2) if sumw2 >= 0.0 else float("nan")
        relerr = err / sumw if sumw != 0.0 else 0.0
        neff = sumw * sumw / sumw2 if sumw2 > 0.0 else 0.0

        print(
            "[WSTAT] {} E=[{:.1f},{:.1f}) n={} sumw={:.8g} err={:.8g} "
            "relerr={:.4g} neff={:.4g} maxw={:.8g} max_entry={}".format(
                label,
                lo,
                hi,
                s["n"],
                sumw,
                err,
                relerr,
                neff,
                s["maxw"],
                s["max_entry"],
            )
        )


def _is_cv_universe(universe):
    try:
        return universe.ShortName() == "cv"
    except Exception:
        return False


def _safe_reco_lepton_energy_gev(universe):
    try:
        return universe.LeptonEnergy() * 1e-3
    except Exception:
        return None


def _safe_true_energy_gev(universe):
    try:
        return universe.mc_incomingE * 1e-3
    except Exception:
        return None


def _diagnose_weight_pieces(universe, label, entry, energy_gev, weight, threshold=5.0):
    if 5.5 <= energy_gev < 6.5:
        if abs(weight) > threshold:
            print(
                "[HIGH_WEIGHT_6GEV] label={} entry={} E={:.6g} weight={:.8g}".format(
                    label, entry, energy_gev, weight
                )
            )
            try:
                universe.DebugStandardWeightPieces(label)
            except Exception as e:
                print("[HIGH_WEIGHT_6GEV] DebugStandardWeightPieces failed:", e)

def plotRecoKin(mc, chainwrapper, outfile):
    """ The main code of event selection """
    kin_cal = KinematicsCalculator(correct_beam_angle=True, correct_MC_energy_scale=False, calc_true = mc, is_pc = AnalysisConfig.is_pc)
    eventClassifier = EventClassifier(classifiers=["Reco","Truth"] if mc else ["Reco"], use_kin_cuts=True, use_sideband = AnalysisConfig.sidebands)
    universes = GetAllSystematicsUniverses(chainwrapper, not mc, AnalysisConfig.is_pc, AnalysisConfig.exclude_universes)
    # universes = GetAllSystematicsUniverses(
    #     chainwrapper,
    #     not mc,
    #     AnalysisConfig.is_pc,
    #     AnalysisConfig.exclude_universes,
    #     AnalysisConfig.playlist,
    # )

    # cv_universe = universes.get("cv", [None])[0]
    # flux_universes = universes.get("Flux", [])

    for univ in chain.from_iterable(iter(universes.values())):
        univ.LoadTools(kin_cal,eventClassifier)

    Plots = preparePlots(universes,mc)
    nEvents = chainwrapper.GetEntries()
    debug_prints = 0
    reco_trueE_wstats = _new_wstat()
    reco_recoLepE_wstats = _new_wstat()
    print(f"Total number of events RECO: ", {nEvents})
    if AnalysisConfig.testing and nEvents > 10000:
        nEvents = 10000
    print("plotRecoKin, mc ",mc)
    setAlarm = AnalysisConfig.grid
    for counter in range(nEvents):
        #1/4 hour for 10k event, should be more than needed unless stuck in I/O
        if counter %10000 == 0:
            print(counter)
            if setAlarm:
                signal.alarm(900)
        flux_debug_this_event = 0
        selected_event_for_debug = False

        for universe in chain.from_iterable(iter(universes.values())):
            universe.SetEntry(counter)
            universe._debug_entry = counter
            universe.__dict__.pop("_cached_new_true_l_over_e", None)
            universe.ResetWeight()
            if mc and AnalysisConfig.skip_2p2h and universe.mc_intType==8:
                continue

            #only update kin_cal & eventClassifier when universe in not vertical only.
            if not universe.IsVerticalOnly():
                kin_cal.CalculateKinematics(universe)
                eventClassifier.Classify(universe)

            # if eventClassifier.side_band is not None or eventClassifier.is_true_signal:
            #     for entry in Plots:
            #         entry.Process(universe)
            if eventClassifier.side_band is not None or eventClassifier.is_true_signal:

                if mc and _is_cv_universe(universe):
                    w = universe.GetWeight()
                    trueE = _safe_true_energy_gev(universe)
                    recoLepE = _safe_reco_lepton_energy_gev(universe)

                    _fill_wstat(reco_trueE_wstats, trueE, w, counter)
                    _fill_wstat(reco_recoLepE_wstats, recoLepE, w, counter)

                    if trueE is not None:
                        _diagnose_weight_pieces(
                            universe,
                            "RECOCUT_TRUEE_6GEV",
                            counter,
                            trueE,
                            w,
                            threshold=5.0,
                        )

                    if recoLepE is not None:
                        _diagnose_weight_pieces(
                            universe,
                            "RECOCUT_RECOLEPE_6GEV",
                            counter,
                            recoLepE,
                            w,
                            threshold=5.0,
                        )

                for entry in Plots:
                    entry.Process(universe)

    if mc:
        _print_wstat(reco_trueE_wstats, "RECO_SELECTED_CV_TRUE_E")
        _print_wstat(reco_recoLepE_wstats, "RECO_SELECTED_CV_RECO_LEPTON_E")
    signal.alarm(0)
    outfile.cd()
    for entry in Plots:
        #print entry.histwrapper.name
        entry.Finalize()

def plotTruthKin(chainwrapper,outfile):
    kin_cal = KinematicsCalculator(correct_beam_angle=True, correct_MC_energy_scale=False, calc_true = True, calc_reco = False)
    eventClassifier = EventClassifier(classifiers=["Truth"],use_kin_cuts=True, use_sideband=[])
    universes = GetAllSystematicsUniverses(chainwrapper, False)
    # universes = GetAllSystematicsUniverses(
    #     chainwrapper,
    #     False,
    #     AnalysisConfig.is_pc,
    #     AnalysisConfig.exclude_universes,
    #     AnalysisConfig.playlist,
    # )

    for univ in chain.from_iterable(iter(universes.values())):
        univ.LoadTools(kin_cal,eventClassifier)
    nEvents = chainwrapper.GetEntries()
    print(f"Total number of events TRUTH: ", {nEvents})
    Plots = prepareTruthPlots(universes)
    print("plotTruthKin")
    truth_trueE_wstats = _new_wstat()
    if AnalysisConfig.testing and nEvents > 10000:
        nEvents = 10000
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

            # if eventClassifier.is_true_signal:
            #     for entry in Plots:
            #         entry.Process(universe)
            if eventClassifier.is_true_signal:

                if _is_cv_universe(universe):
                    w = universe.GetWeight()
                    trueE = _safe_true_energy_gev(universe)

                    _fill_wstat(truth_trueE_wstats, trueE, w, counter)

                    if trueE is not None:
                        _diagnose_weight_pieces(
                            universe,
                            "TRUTH_SIGNAL_TRUEE_6GEV",
                            counter,
                            trueE,
                            w,
                            threshold=5.0,
                        )

                for entry in Plots:
                    entry.Process(universe)

    _print_wstat(truth_trueE_wstats, "TRUTH_SIGNAL_CV_TRUE_E")
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
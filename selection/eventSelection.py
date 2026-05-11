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
    kin_cal = KinematicsCalculator(correct_beam_angle=True, correct_MC_energy_scale=False, calc_true = mc, is_pc = AnalysisConfig.is_pc)
    eventClassifier = EventClassifier(classifiers=["Reco","Truth"] if mc else ["Reco"], use_kin_cuts=True, use_sideband = AnalysisConfig.sidebands)
    universes = GetAllSystematicsUniverses(chainwrapper, not mc, AnalysisConfig.is_pc, AnalysisConfig.exclude_universes)

    cv_universe = universes.get("cv", [None])[0]
    flux_universes = universes.get("Flux", [])

    for univ in chain.from_iterable(iter(universes.values())):
        univ.LoadTools(kin_cal,eventClassifier)

    Plots = preparePlots(universes,mc)
    nEvents = chainwrapper.GetEntries()
    debug_prints = 0
    
    def get_reco_en4_gev(u):
        """
        Reconstruct the EN4 variable used in the histogram:
        EN4 = E_e + E_avail
        assuming ElectronEnergy() and AvailableEnergy() return MeV.
        """
        return (u.ElectronEnergy() + u.AvailableEnergy()) * 1e-3


    def get_en4_debug_bin(en4):
        bins = [0.0, 2.0, 2.5, 3.0, 3.8, 4.6, 5.5, 6.5, 8.0, 10.0, 12.5, 16.0, 20.0]
        for i in range(len(bins) - 1):
            if bins[i] <= en4 < bins[i + 1]:
                return i + 1
        return None


    def debug_flux_cv_vs_universe(counter, cv_universe, flux_universe, label=""):
        pdg0 = cv_universe.mc_incoming
        pdg = pdg0

        true_enu = cv_universe.mc_incomingE * 1e-3
        reco_en4 = get_reco_en4_gev(cv_universe)
        en4_bin = get_en4_debug_bin(reco_en4)

        cv_flux = cv_universe.GetFluxAndCVWeight(true_enu, pdg)
        univ_flux = flux_universe.GetFluxAndCVWeight(true_enu, pdg)

        flux_ratio = univ_flux / cv_flux if cv_flux != 0 else -999

        cv_total = cv_universe.GetWeight(False)
        univ_total = flux_universe.GetWeight(False)
        total_ratio = univ_total / cv_total if cv_total != 0 else -999

        print(
            "[FLUX_LAST2BIN_CHECK] {} counter={} "
            "EN4={:.6g} EN4bin={} trueEv={:.6g} "
            "pdg0={} pdgUsed={} intType={} current={} "
            "cvFlux={:.8g} univFlux={:.8g} fluxRatio={:.8g} "
            "cvTotalWeight={:.8g} univTotalWeight={:.8g} totalRatio={:.8g} "
            "cvShort={} univShort={} univSigma={}".format(
                label,
                counter,
                reco_en4,
                en4_bin,
                true_enu,
                pdg0,
                pdg,
                cv_universe.mc_intType,
                cv_universe.mc_current,
                cv_flux,
                univ_flux,
                flux_ratio,
                cv_total,
                univ_total,
                total_ratio,
                cv_universe.ShortName(),
                flux_universe.ShortName(),
                flux_universe.GetSigma(),
            )
        )

    def debug_mnvtune_universe(counter, cv_universe, universe, label=""):
        true_enu = cv_universe.mc_incomingE * 1e-3
        reco_en4 = get_reco_en4_gev(cv_universe)
        en4_bin = get_en4_debug_bin(reco_en4)

        cv_total = cv_universe.GetWeight(False)
        univ_total = universe.GetWeight(False)
        total_ratio = univ_total / cv_total if cv_total != 0 else -999

        # IMPORTANT:
        # Do not use hasattr/getattr here because CVPythonUniverse.__getattr__
        # will try to read missing names as TTree branches.
        if "iweight" in universe.__dict__:
            extra_id = "iweight={}".format(universe.__dict__["iweight"])
        elif "universe_number" in universe.__dict__:
            extra_id = "universe_number={}".format(universe.__dict__["universe_number"])
        else:
            extra_id = "extraID=NA"

        print(
            "[MNVTUNE_EVT] {} counter={} "
            "histEN4={:.6g} EN4bin={} trueEv={:.6g} "
            "truth_class={} sideband={} intType={} current={} Q2={:.6g} W={:.6g} "
            "cvTotal={:.8g} univTotal={:.8g} totalRatio={:.8g} "
            "short={} class={} sigma={} {}".format(
                label,
                counter,
                reco_en4,
                en4_bin,
                true_enu,
                cv_universe.classifier.truth_class,
                cv_universe.classifier.side_band,
                cv_universe.mc_intType,
                cv_universe.mc_current,
                cv_universe.mc_Q2 / 1e6 if "mc_Q2" in cv_universe.LeafGetters else cv_universe.mc_Q2 / 1e6,
                cv_universe.mc_w / 1e3 if "mc_w" in cv_universe.LeafGetters else cv_universe.mc_w / 1e3,
                cv_total,
                univ_total,
                total_ratio,
                universe.ShortName(),
                type(universe).__name__,
                universe.GetSigma(),
                extra_id,
            )
        )

    debug_bands = {
        "LowQ2Pi",
        "fsi_weight",
        "Low_Recoil_2p2h_Tune",
        "SuSA_Valencia_Weight",
        "MK_model",
    }

    debug_categories = {
        "CCNuEDelta",
        "CCNuEDIS",
        "CCNuEQE",
        "NCPi0",
        "CCPi0",
    }

    mnvtune_debug_prints = 0
    max_mnvtune_debug_prints = 300

    print(f"Total number of events RECO: ", {nEvents})
    if AnalysisConfig.testing and nEvents > 1000:
        nEvents = 1000
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
            universe.__dict__.pop("_cached_new_true_l_over_e", None)
            universe.ResetWeight()
            if mc and AnalysisConfig.skip_2p2h and universe.mc_intType==8:
                continue

            #only update kin_cal & eventClassifier when universe in not vertical only.
            if not universe.IsVerticalOnly():
                kin_cal.CalculateKinematics(universe)
                eventClassifier.Classify(universe)

            if eventClassifier.side_band is not None or eventClassifier.is_true_signal:
                # if (
                #     mc
                #     and mnvtune_debug_prints < max_mnvtune_debug_prints
                #     and universe.ShortName() in debug_bands
                #     and eventClassifier.truth_class in debug_categories
                # ):
                #     debug_mnvtune_universe(counter, cv_universe, universe, "before_fill")
                #     if "DebugExtraWeight" in universe.__class__.__dict__:
                #         universe.DebugExtraWeight("before_fill")
                #     mnvtune_debug_prints += 1
                for entry in Plots:
                    entry.Process(universe)

    signal.alarm(0)
    outfile.cd()
    for entry in Plots:
        #print entry.histwrapper.name
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
    plots=[]
    for entry in HISTS_TO_MAKE:
        if not (isinstance(entry,str) and entry.startswith("True Signal")):
            continue
        settings = {"key":entry,"region":"Signal","mc":True}
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
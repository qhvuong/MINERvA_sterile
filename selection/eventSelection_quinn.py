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

#start loading my modules
from tools import Utilities
print("Loaded Utilities from:", Utilities.__file__)
from config.PlotConfig import HISTS_TO_MAKE
from config.AnalysisConfig import AnalysisConfig
from tools.SystematicsUniverse import GetAllSystematicsUniverses
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
    for univ in chain.from_iterable(iter(universes.values())):
        univ.LoadTools(kin_cal,eventClassifier)

    Plots = preparePlots(universes,mc)
    nEvents = chainwrapper.GetEntries()
    print(f"Total number of events RECO: ", {nEvents})
    if AnalysisConfig.testing and nEvents > 1000:
        nEvents = 1000
    # nEvents = 10

    # def is_flux(u):
    #     try:
    #         return u.ShortName() == "Flux"
    #     except Exception:
    #         return False

    setAlarm = AnalysisConfig.grid
    for counter in range(nEvents):
        #1/4 hour for 10k event, should be more than needed unless stuck in I/O
        if counter %10000 == 0:
            print(counter)
            if setAlarm:
                signal.alarm(900)

        for universe in chain.from_iterable(iter(universes.values())):
            universe.SetEntry(counter)
            # print("counter: ", counter)
            # print("universe name: ", type(universe).__name__)
            # universe.ResetWeight()
            if mc and AnalysisConfig.skip_2p2h and universe.mc_intType==8:
                continue

            #only update kin_cal & eventClassifier when universe in not vertical only.
            if not universe.IsVerticalOnly():
                kin_cal.CalculateKinematics(universe)
                eventClassifier.Classify(universe)

            if eventClassifier.side_band is not None or eventClassifier.is_true_signal:
                for entry in Plots:
                    # if counter == 0 and type(universe).__name__ == "Flux":
                    # if type(universe).__name__ == "FluxUniverse":
                    #     print(" calling Process for FluxUniverse: ", entry, counter)
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
    if AnalysisConfig.testing and nEvents > 1000:
        nEvents = 1000
    for counter in range(nEvents):
        #half an hour for a event, should be much more than needed unless stuck in I/O
        if counter %100000 == 0:
            print(counter)

        for universe in chain.from_iterable(iter(universes.values())):
            universe.SetEntry(counter)
            # universe.ResetWeight()

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
        print("HISTS_TO_MAKE", entry)
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
       
        output_file = ROOT.TFile.Open(outputSelectionHistogram,"RECREATE")
        
        CopyMetaTreeToOutPutFile(output_file)
        if st=="mc" and Truth:
            print("selecting truth")
            plotTruthKin(Utilities.fileChain(AnalysisConfig.playlist,"mc",AnalysisConfig.ntuple_tag,"Truth",AnalysisConfig.count[0],AnalysisConfig.count[1]),output_file)
        if Reco :
            print("selecting reco")
            plotRecoKin(st=="mc", Utilities.fileChain(AnalysisConfig.playlist,st,AnalysisConfig.ntuple_tag,None,AnalysisConfig.count[0],AnalysisConfig.count[1]), output_file)
            print("done selection reco")
        output_file.Close()
        print("selection is done for ", st, AnalysisConfig.playlist)

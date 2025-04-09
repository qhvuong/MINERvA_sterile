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
import numpy as np
import ast
import array

#start loading my modules
from tools import Utilities
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
    #print("AnalysisConfig.is_pc = ", AnalysisConfig.is_pc)
    kin_cal = KinematicsCalculator(correct_beam_angle=True, correct_MC_energy_scale=False, calc_true = mc, is_pc = AnalysisConfig.is_pc)


    eventClassifier = EventClassifier(classifiers=["Reco","Truth"] if mc else ["Reco"], use_kin_cuts=True, use_sideband = AnalysisConfig.sidebands)
    
    universes = GetAllSystematicsUniverses(chainwrapper, not mc, AnalysisConfig.is_pc, AnalysisConfig.exclude_universes)
    
    for univ in chain.from_iterable(iter(universes.values())):
        #print("univ: ", univ)
        univ.LoadTools(kin_cal,eventClassifier)
    
    Plots = preparePlots(universes,mc)
    nEvents = chainwrapper.GetEntries()
    print(f"Total number of events RECO: ", {nEvents})
    if AnalysisConfig.testing and nEvents > 1000:
        nEvents = 1000

    #cv_universe = list(universes["cv"])[0]

    setAlarm = AnalysisConfig.grid
    for counter in range(nEvents):
        #1/4 hour for 10k event, should be more than needed unless stuck in I/O
        if counter %10000 == 0:
            print(counter)
            if setAlarm:
                signal.alarm(900)
        
        for universe in chain.from_iterable(iter(universes.values())):
            """
            if AnalysisConfig.flavor_swap_type in ["nue_swapped", "nuebar_swapped"]:
                if universe.ShortName() != "cv":
                    continue
            """
            #print("universe: ", universe.values())
            universe.SetEntry(counter)
            if mc and AnalysisConfig.skip_2p2h and universe.mc_intType==8:
                continue
        
            #only update kin_cal & eventClassifier when universe in not vertical only.
            if not universe.IsVerticalOnly():
                kin_cal.CalculateKinematics(universe)
                eventClassifier.Classify(universe)
                #TruthTools.Print(universe)
                #file1 = open("background_debug.txt","a")
                #if eventClassifier.side_band == "Signal" and eventClassifier.truth_class in ["CCNuEAntiNu","CCNuEQE","CCNuEDelta","CCNuEDIS","CCNuE2p2h","CCNuE"]:
                #    file1.write((TruthTools.Print(universe))+"\n")
            """
            if universe.ShortName() == "cv" and eventClassifier.is_true_signal:
                print("🔹 MC incoming E [GeV]:", universe.mc_incomingE / 1000.)
                print("🔸 event.GetWeight():", universe.GetWeight())
            """
        
            if eventClassifier.side_band is not None or eventClassifier.is_true_signal:
                for entry in Plots:
                    #print("Plots entry: ", entry)
                    entry.Process(universe)

            #if universe.ShortName() =="cv" and eventClassifier.is_true_signal:
            #    print universe.kin_cal.true_visE, universe.kin_cal.true_q3, universe.GetWeight()
            #    raw_input("...")

        
        
    signal.alarm(0)
    outfile.cd()
    for entry in Plots:
        #print(entry.histwrapper.name)
        entry.Finalize()
    

    #eventClassifier.GetStatTree().Write("",ROOT.TObject.kOverwrite)
    #print("counter for reco,truth: ", eventClassifier.counter)

def plotTruthKin(chainwrapper,outfile):
    kin_cal = KinematicsCalculator(correct_beam_angle=True, correct_MC_energy_scale=False, calc_true = True, calc_reco = False)
    eventClassifier = EventClassifier(classifiers=["Truth"],use_kin_cuts=True, use_sideband=[])
    universes = GetAllSystematicsUniverses(chainwrapper, False)
    for univ in chain.from_iterable(iter(universes.values())):
        univ.LoadTools(kin_cal,eventClassifier)

    nEvents = chainwrapper.GetEntries()
    print(f"Total number of events TRUTH: ", {nEvents})
    #output_file = ROOT.TFile.Open(outname,"RECREATE")
    Plots = prepareTruthPlots(universes)
    if AnalysisConfig.testing and nEvents > 1000:
        nEvents = 1000
    for counter in range(nEvents):
        #half an hour for a event, should be much more than needed unless stuck in I/O
        if counter %100000 == 0:
            print(counter)

        for universe in chain.from_iterable(iter(universes.values())):
            universe.SetEntry(counter)

            #only update kin_cal & eventClassifier when universe in not vertical only.
            if not universe.IsVerticalOnly():
                kin_cal.CalculateKinematics(universe)
                
                reco_before = eventClassifier.counter[0]
                true_before = eventClassifier.counter[1]
                eventClassifier.Classify(universe)
                #if eventClassifier.counter[0] > reco_before:
                #    print("event %i in reco" % counter)
                #    TruthTools.Print(universe)
                #if eventClassifier.counter[1] > true_before:
                #    print("event %i in truth" % counter)
                #    TruthTools.Print(universe)

            if eventClassifier.is_true_signal:
                for entry in Plots:
                    entry.Process(universe)

    outfile.cd()
    for entry in Plots:
        entry.Finalize()

    #print("counter for reco,truth: ", eventClassifier.counter)

def preparePlots(universes,mc):
    # make a bunch of Plot Processor, grouped by signal/sideband
    plots=set([])
    #for region in ["Signal"]+AnalysisConfig.sidebands:

    for entry in HISTS_TO_MAKE:
        print("HISTS_TO_MAKE entry ", entry)
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
        print("HISTS_TO_MAKE entry", entry)
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
        print("AnalysisConfig.data_types: ", st)
        outputSelectionHistogram = AnalysisConfig.SelectionHistoPath(AnalysisConfig.playlist,"data" in st,True)
       
        output_file = ROOT.TFile.Open(outputSelectionHistogram,"RECREATE")
        
        CopyMetaTreeToOutPutFile(output_file)
        
        if Reco :
            print("selecting reco playlist: %s" % AnalysisConfig.playlist)
            #cProfile.run('plotRecoKin(st=="mc", Utilities.fileChain(AnalysisConfig.playlist,st,AnalysisConfig.ntuple_tag,None,AnalysisConfig.count[0],AnalysisConfig.count[1]), output_file)')
            plotRecoKin(st=="mc", Utilities.fileChain(AnalysisConfig.playlist,st,AnalysisConfig.ntuple_tag,None,AnalysisConfig.count[0],AnalysisConfig.count[1]), output_file)
            print("done selection reco")
        
        if st=="mc" and Truth:
            print("selecting truth")
            plotTruthKin(Utilities.fileChain(AnalysisConfig.playlist,"mc",AnalysisConfig.ntuple_tag,"Truth",AnalysisConfig.count[0],AnalysisConfig.count[1]),output_file)
        output_file.Close()
        print("selection is done for ", st, AnalysisConfig.playlist)
        
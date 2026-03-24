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
    print("plotRecoKin, mc ",mc)
    setAlarm = AnalysisConfig.grid
    for counter in range(nEvents):
        #1/4 hour for 10k event, should be more than needed unless stuck in I/O
        if counter %10000 == 0:
            print(counter)
            if setAlarm:
                signal.alarm(900)


        for universe in chain.from_iterable(iter(universes.values())):
            universe.SetEntry(counter)
            universe.ResetWeight()
            if mc and AnalysisConfig.skip_2p2h and universe.mc_intType==8:
                continue

            #only update kin_cal & eventClassifier when universe in not vertical only.
            if not universe.IsVerticalOnly():
                kin_cal.CalculateKinematics(universe)
                eventClassifier.Classify(universe)

            # # if mc and universe.GetSigma() == 0 and debug_prints < max_debug_prints:
            # if mc and universe.GetSigma() == 0 and eventClassifier.is_reco_signal == True:
            #     p = universe.LeptonP3D()
            #     px, py, pz = p.X(), p.Y(), p.Z()
            #     pt = math.sqrt(px*px + py*py)

            #     print(f"PX={px:.1f}  PY={py:.1f}  PZ={pz:.1f}  PT={pt:.1f}")
                # # debug_prints += 1
                # cone_outside_e = universe.ConeOutsideE()
                # neighborhood_e = universe.NeighborhoodE(_debug=True)
                # print(
                #     "entry={} sideband={} reco_signal={} true_signal={} ConeOutsideE value={:.6g} NeighborhoodE value={:.6g}".format(
                #         counter,
                #         eventClassifier.side_band,
                #         eventClassifier.is_reco_signal,
                #         eventClassifier.is_true_signal,
                #         cone_outside_e,
                #         neighborhood_e,
                #     )
                # )

            # # --- put this near the top of plotRecoKin (before the event loop) ---
            # debug_prints = 0
            # max_debug_prints = 50   # change to whatever you want

            # # --- inside your event loop, after eventClassifier.Classify(universe) ---
            # # Only print/debug for events that PASS your selection logic
            # if (eventClassifier.side_band is not None or eventClassifier.is_true_signal):

            #     # Print only for CV universe (and only up to N events)
            #     if universe.GetSigma() == 0 and debug_prints < max_debug_prints:
            #         debug_prints += 1

            #         print("\n" + "="*80)
            #         print(f"[post-cut] entry={counter} univ={type(universe).__name__}  sideband={eventClassifier.side_band}  is_true_signal={eventClassifier.is_true_signal}")

            #         # Optional detailed reco electron-theta debug (uses prong_part_E)
            #         if hasattr(universe, "DebugElectronTheta"):
            #             universe.DebugElectronTheta(
            #                 prong=0,
            #                 entry_tag=f"entry={counter}, univ={type(universe).__name__}"
            #             )

            #         # --- reco theta (what your selection is effectively using) ---
            #         try:
            #             reco_theta = universe.ElectronTheta()
            #             # print(f"[reco] theta_rot = {reco_theta:.6f} rad")
            #         except Exception as e:
            #             reco_theta = None
            #             print(f"[reco] ElectronTheta() failed: {e}")

            #         # --- truth vs reco theta (MC only) using kin_cal truth already computed ---
            #         if mc:
            #             n = universe.GetInt("mc_nFSPart")

            #             # collect PDGs
            #             pdgs = [int(universe.GetVecElem("mc_FSPartPDG", i)) for i in range(n)]

            #             # pick a final-state charged lepton if present
            #             lep_pdgs = [p for p in pdgs if abs(p) in (11, 13, 15)]
            #             prim_lep_pdg = lep_pdgs[0] if lep_pdgs else None

            #             print(f"[truth] mc_nFSPart={n} prim_lep_pdg={prim_lep_pdg}  (lep candidates={lep_pdgs})")
            #             # optional: also print first few FS PDGs so you see the ordering
            #             print(f"[truth] first FS PDGs: {pdgs}")

            #             px0 = universe.GetVecElem("mc_primFSLepton", 0)
            #             py0 = universe.GetVecElem("mc_primFSLepton", 1)
            #             pz0 = universe.GetVecElem("mc_primFSLepton", 2)
            #             E0  = universe.GetVecElem("mc_primFSLepton", 3)
            #             print(f"[truth] (px,py,pz,E)=({px0:.3g}, {py0:.3g}, {pz0:.3g}, {E0:.3g})")

            #             true_theta = getattr(kin_cal, "true_theta_lep_rad", None)
            #             true_Etheta2 = E0*true_theta**2
            #             print(f"[truth] true Etheta2={true_Etheta2:.6g} MeV")

            #             if true_theta is None:
            #                 print("[truth] kin_cal.true_theta_lep_rad missing (truth theta not available)")
            #             else:
            #                 if reco_theta is not None:
            #                     print(f"[theta check] reco={reco_theta:.6f} true={true_theta:.6f}  (reco-true)={reco_theta-true_theta:+.6f}")
            #                 else:
            #                     print(f"[truth] true_theta = {true_theta:.6f} rad (reco missing)")

            #         # --- prong-level info (pid + prong0 4-vector) ---
            #         try:
            #             n_prongs = universe.GetInt("n_prongs")
            #             pid0 = universe.GetVecElem("prong_part_pid", 0) if n_prongs > 0 else None
            #             pid1 = universe.GetVecElem("prong_part_pid", 1) if n_prongs > 1 else None
            #             print(f"[prongs] n_prongs={n_prongs} pid0={pid0} pid1={pid1}")

            #             # if n_prongs > 0:
            #             #     px0 = universe.GetVecElem("prong_part_E", 0, 0)
            #             #     py0 = universe.GetVecElem("prong_part_E", 0, 1)
            #             #     pz0 = universe.GetVecElem("prong_part_E", 0, 2)
            #             #     E0  = universe.GetVecElem("prong_part_E", 0, 3)
            #             #     print(f"[prong0] (px,py,pz,E)=({px0:.3g}, {py0:.3g}, {pz0:.3g}, {E0:.3g})")
            #         except Exception as e:
            #             print(f"[prongs] failed to read prong branches: {e}")

            #         # EMLikeTrackScore lives in prong_part_score[n_prongs]/D
            #         try:
            #             n_prongs = universe.GetInt("n_prongs")
            #             score0 = universe.GetVecElem("prong_part_score", 0) if n_prongs > 0 else None
            #             print(f"[score] n_prongs={n_prongs} score0={score0}")
            #         except Exception as e:
            #             print(f"[score] failed to read prong_part_score: {e}")

            #         # --- Vertex-dependent nu direction, done in BEAM COORDS (matches ElectronP3D) ---
            #         vx = universe.GetVecElem("vtx", 0)  # mm
            #         vy = universe.GetVecElem("vtx", 1)  # mm
            #         vz = universe.GetVecElem("vtx", 2)  # mm

            #         Lmm = 900e3  # 900 m in mm (tuneable)

            #         # Rotate vertex into the same "beam frame" as ElectronP3D() uses
            #         vtx_det = ROOT.TVector3(vx, vy, vz)
            #         vtx_beam = ROOT.TVector3(vx, vy, vz)
            #         vtx_beam.RotateX(SystematicsConfig.BEAM_ANGLE)

            #         # Neutrino direction in beam frame:
            #         # source ~ (0,0,z_beam - Lmm)  -> direction to vertex is (x_beam, y_beam, Lmm)
            #         nu_dir_beam = ROOT.TVector3(vtx_beam.X(), vtx_beam.Y(), Lmm).Unit()

            #         # Lepton direction in beam frame (ElectronP3D already rotated)
            #         p = universe.ElectronP3D()  # ROOT.Math.XYZVector in beam frame
            #         lep_dir_beam = ROOT.TVector3(p.X(), p.Y(), p.Z()).Unit()

            #         theta_vtxcorr = lep_dir_beam.Angle(nu_dir_beam)  # radians

            #         E_GeV = universe.ElectronEnergy()/1e3
            #         Etheta2_vtxcorr_MeV = (E_GeV * theta_vtxcorr**2) * 1e3

            #         print(f"[vtxcorr] vtx_det(mm)=({vx:.1f},{vy:.1f},{vz:.1f})  vtx_beam(mm)=({vtx_beam.X():.1f},{vtx_beam.Y():.1f},{vtx_beam.Z():.1f})  L(m)={Lmm/1e3:.0f}")
            #         print(f"[vtxcorr] theta={theta_vtxcorr:.6f} rad  Etheta2={Etheta2_vtxcorr_MeV:.3f} MeV")

            if eventClassifier.side_band is not None or eventClassifier.is_true_signal:
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
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
    """
    Main reconstructed-event selection.

    Optional muon-angle validation:
      Compare MAT GetThetamu()/GetPhimu() against the angles
      calculated from MuonP3D_beam() using the same MC event.

    Enable with:
        DEBUG_MUON_ANGLES=1

    This diagnostic does not change event selection, weights,
    systematic universes, or histogram filling.
    """
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

    for univ in chain.from_iterable(iter(universes.values())):
        univ.LoadTools(kin_cal, eventClassifier)

    Plots = preparePlots(universes, mc)

    nEvents = chainwrapper.GetEntries()
    print("Total number of events RECO:", nEvents)

    if AnalysisConfig.testing and nEvents > 2000:
        nEvents = 2000

    print("plotRecoKin, mc", mc)

    setAlarm = AnalysisConfig.grid

    # ============================================================
    # Optional CV muon-angle comparison
    # ============================================================

    from config import SystematicsConfig

    debug_muon_angles = (
        os.environ.get("DEBUG_MUON_ANGLES", "0") == "1"
        and mc
        and abs(SystematicsConfig.AnaNuPDG) == 14
    )

    max_debug_events = 5000
    max_print_outliers = 10

    angle_stats = {
        "n": 0,
        "skipped": 0,
        "sum_dtheta": 0.0,
        "sum_dphi": 0.0,
        "sum_dtheta2": 0.0,
        "sum_dphi2": 0.0,
        "max_dtheta": 0.0,
        "max_dphi": 0.0,
        "outliers": [],
    }

    if debug_muon_angles:
        print("\n========== MUON ANGLE VALIDATION ==========")
        print("MAT: GetThetamu(), GetPhimu()")
        print("Vector: MuonP3D_beam().Theta(), .Phi()")
        print("All reported angles are in degrees.")
        print("Comparisons use the CV universe only.")
        print("Maximum comparisons:", max_debug_events)
        print("===========================================\n")

    cv_universe = universes["cv"][0] if debug_muon_angles else None

    def compare_muon_angles(univ, counter):
        """
        Compare MAT angles to momentum-derived angles.

        Ignore invalid momenta and nonfinite angular values.
        This is a diagnostic, not an event-selection cut.
        """
        stats = angle_stats

        if stats["n"] >= max_debug_events:
            return

        p_det = univ.MuonP3D_det()
        p_beam = univ.MuonP3D_beam()

        # Only compare forward-going muons.
        # MAT's phi3D() uses angular slopes, which can differ by
        # 180 degrees from the momentum-vector phi for backward tracks.
        if p_beam.Z() <= 0:
            return

        components = (
            float(p_det.X()),
            float(p_det.Y()),
            float(p_det.Z()),
            float(p_beam.X()),
            float(p_beam.Y()),
            float(p_beam.Z()),
        )

        if (
            not all(math.isfinite(x) for x in components)
            or not math.isfinite(float(p_beam.R()))
            or p_beam.R() <= 0
        ):
            stats["skipped"] += 1
            return

        theta_mat = float(univ.GetThetamu())
        phi_mat = float(univ.GetPhimu())

        theta_vec = float(p_beam.Theta())
        phi_vec = float(p_beam.Phi())

        angles = (theta_mat, phi_mat, theta_vec, phi_vec)

        if not all(math.isfinite(x) for x in angles):
            stats["skipped"] += 1
            return

        if not (
            0.0 <= theta_mat <= math.pi
            and 0.0 <= theta_vec <= math.pi
        ):
            stats["skipped"] += 1
            return

        dtheta = math.degrees(theta_vec - theta_mat)

        # Wrap phi difference into [-pi, pi].
        # This avoids false ~360-degree disagreements.
        dphi_rad = math.atan2(
            math.sin(phi_vec - phi_mat),
            math.cos(phi_vec - phi_mat),
        )
        dphi = math.degrees(dphi_rad)

        stats["n"] += 1
        stats["sum_dtheta"] += dtheta
        stats["sum_dphi"] += dphi
        stats["sum_dtheta2"] += dtheta * dtheta
        stats["sum_dphi2"] += dphi * dphi

        stats["max_dtheta"] = max(
            stats["max_dtheta"], abs(dtheta)
        )
        stats["max_dphi"] = max(
            stats["max_dphi"], abs(dphi)
        )

        result = {
            "entry": counter,
            "theta_mat": math.degrees(theta_mat),
            "theta_vec": math.degrees(theta_vec),
            "phi_mat": math.degrees(phi_mat),
            "phi_vec": math.degrees(phi_vec),
            "dtheta": dtheta,
            "dphi": dphi,
            "p_det": components[:3],
            "p_beam": components[3:],
            "pmu_mat": float(univ.GetPmu()),
        }

        # Keep events with the largest combined angular difference.
        outliers = stats["outliers"]
        outliers.append(result)

        outliers.sort(
            key=lambda r: max(abs(r["dtheta"]), abs(r["dphi"])),
            reverse=True,
        )

        if len(outliers) > max_print_outliers:
            outliers.pop()

    # ============================================================
    # Standard reconstructed-event loop
    # ============================================================

    for counter in range(nEvents):

        if counter % 10000 == 0:
            print(counter)

            if setAlarm:
                signal.alarm(900)

        flux_debug_this_event = 0
        selected_event_for_debug = False

        for universe in chain.from_iterable(iter(universes.values())):

            universe.SetEntry(counter)
            universe.__dict__.pop(
                "_cached_new_true_l_over_e", None
            )
            universe.ResetWeight()

            # Diagnostic before classification so that the angle
            # comparison is independent of selection decisions.
            if (
                debug_muon_angles
                and universe is cv_universe
                and angle_stats["n"] < max_debug_events
            ):
                compare_muon_angles(universe, counter)

            if (
                mc
                and AnalysisConfig.skip_2p2h
                and universe.mc_intType == 8
            ):
                continue

            # Recalculate kinematics and classification only
            # for universes with lateral shifts.
            if not universe.IsVerticalOnly():
                kin_cal.CalculateKinematics(universe)
                eventClassifier.Classify(universe)

            if (
                eventClassifier.side_band is not None
                or eventClassifier.is_true_signal
            ):
                for entry in Plots:
                    entry.Process(universe)

    # ============================================================
    # Muon-angle diagnostic summary
    # ============================================================

    if debug_muon_angles:
        stats = angle_stats
        n = stats["n"]

        print("\n========== MUON ANGLE COMPARISON ==========")
        print("Valid comparisons:", n)
        print("Invalid entries skipped:", stats["skipped"])

        if n > 0:
            print("\n--- Angular differences (degrees) ---")

            print(
                "Delta theta:",
                "mean =", stats["sum_dtheta"] / n,
                "RMS =", math.sqrt(stats["sum_dtheta2"] / n),
                "max abs =", stats["max_dtheta"],
            )

            print(
                "Delta phi:",
                "mean =", stats["sum_dphi"] / n,
                "RMS =", math.sqrt(stats["sum_dphi2"] / n),
                "max abs =", stats["max_dphi"],
            )

            print("\n--- Largest angular discrepancies ---")

            for r in stats["outliers"]:
                print(
                    "\nEntry:", r["entry"],
                    "\n  theta MAT:", r["theta_mat"],
                    "\n  theta vec:", r["theta_vec"],
                    "\n  dtheta:", r["dtheta"],
                    "\n  phi MAT:", r["phi_mat"],
                    "\n  phi vec:", r["phi_vec"],
                    "\n  dphi:", r["dphi"],
                    "\n  p det:", r["p_det"],
                    "\n  p beam:", r["p_beam"],
                    "\n  |p| MAT:", r["pmu_mat"],
                )

        else:
            print("WARNING: No valid muon comparisons!")

        print("\n===========================================\n")

    # ============================================================
    # Standard histogram finalization
    # ============================================================

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
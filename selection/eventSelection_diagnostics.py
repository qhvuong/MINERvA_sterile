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

    Optional diagnostics:
      DEBUG_MUON_ANGLES=1
          Compare MAT GetThetamu()/GetPhimu() with the new
          MuonTheta()/MuonPhi() methods, including backward muons.

      DEBUG_BEAM_ANGLE=1
          Verify that beam-angle systematic universes shift
          the reconstructed muon angle and preserve |p|.

    Neither diagnostic modifies selection, weights, or histograms.
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
    # Optional muon-angle diagnostics
    # ============================================================

    from config import SystematicsConfig

    is_muon_mc = (
        mc and abs(SystematicsConfig.AnaNuPDG) == 14
    )

    debug_muon_angles = (
        is_muon_mc
        and os.environ.get("DEBUG_MUON_ANGLES", "0") == "1"
    )

    debug_beam_angle = (
        is_muon_mc
        and os.environ.get("DEBUG_BEAM_ANGLE", "0") == "1"
    )

    max_debug_events = 5000
    max_print_outliers = 10

    def wrapped_delta_phi(phi1, phi2):
        """Wrapped difference in radians, within [-pi, pi]."""
        return math.atan2(
            math.sin(phi1 - phi2),
            math.cos(phi1 - phi2),
        )

    # CV vs MAT statistics, separated by forward/backward tracks.
    def new_angle_stats():
        return {
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

    angle_stats = {
        "forward": new_angle_stats(),
        "backward": new_angle_stats(),
    }

    cv_universe = universes["cv"][0] if is_muon_mc else None

    # Identify beam-angle universes through the systematic factory,
    # rather than PyROOT's class-name representation.
    beam_universes = universes.get("beam_angle", [])

    if debug_beam_angle:
        print("Available systematic groups:", list(universes.keys()))
        print("Number of beam-angle universes:", len(beam_universes))

        for i, u in enumerate(beam_universes):
            print(
                "Beam universe", i,
                "sigma =", u.GetSigma(),
                "type =", type(u),
            )

    def compare_muon_angles(univ, counter):
        """
        Compare new momentum-derived angles against MAT.

        Include both forward and backward reconstructed muons.
        """
        p_det = univ.MuonP3D_det()
        p_beam = univ.MuonP3D_beam()

        components = (
            float(p_det.X()),
            float(p_det.Y()),
            float(p_det.Z()),
            float(p_beam.X()),
            float(p_beam.Y()),
            float(p_beam.Z()),
        )

        direction = (
            "forward" if p_beam.Z() >= 0 else "backward"
        )
        stats = angle_stats[direction]

        if stats["n"] >= max_debug_events:
            return

        if (
            not all(math.isfinite(x) for x in components)
            or not math.isfinite(float(p_beam.R()))
            or p_beam.R() <= 0
        ):
            stats["skipped"] += 1
            return

        theta_mat = float(univ.GetThetamu())
        phi_mat = float(univ.GetPhimu())

        # New functions added to CVPythonUniverse.
        theta_new = float(univ.MuonTheta())
        phi_new = float(univ.MuonPhi())

        if not all(
            math.isfinite(x)
            for x in (theta_mat, phi_mat, theta_new, phi_new)
        ):
            stats["skipped"] += 1
            return

        if not (
            0 <= theta_mat <= math.pi
            and 0 <= theta_new <= math.pi
        ):
            stats["skipped"] += 1
            return

        dtheta = math.degrees(theta_new - theta_mat)
        dphi = math.degrees(
            wrapped_delta_phi(phi_new, phi_mat)
        )

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
            "theta_new": math.degrees(theta_new),
            "phi_mat": math.degrees(phi_mat),
            "phi_new": math.degrees(phi_new),
            "dtheta": dtheta,
            "dphi": dphi,
            "p_det": components[:3],
            "p_beam": components[3:],
        }

        stats["outliers"].append(result)
        stats["outliers"].sort(
            key=lambda r: max(
                abs(r["dtheta"]),
                abs(r["dphi"]),
            ),
            reverse=True,
        )

        if len(stats["outliers"]) > max_print_outliers:
            stats["outliers"].pop()

    # Collect beam-angle comparisons for each universe.
    beam_stats = {}

    def compare_beam_angle(univ, counter):
        """
        Compare the varied muon direction with CV.

        Check:
          - theta variation
          - wrapped phi variation
          - momentum magnitude preservation
        """
        p_cv = cv_universe.MuonP3D_beam()
        p_var = univ.MuonP3D_beam()

        values = (
            float(p_cv.R()),
            float(p_var.R()),
            float(p_cv.Theta()),
            float(p_var.Theta()),
            float(p_cv.Phi()),
            float(p_var.Phi()),
        )

        if (
            not all(math.isfinite(v) for v in values)
            or p_cv.R() <= 0
            or p_var.R() <= 0
        ):
            return

        key = (
            str(univ.ShortName()),
            float(univ.GetSigma()),
            str(getattr(univ, "rotation", "")),
        )

        if key not in beam_stats:
            beam_stats[key] = {
                "n": 0,
                "max_dtheta": 0.0,
                "max_dphi": 0.0,
                "max_rel_dp": 0.0,
                "sum_dtheta2": 0.0,
                "sum_dphi2": 0.0,
            }

        stats = beam_stats[key]

        if stats["n"] >= max_debug_events:
            return

        dtheta = math.degrees(
            float(p_var.Theta() - p_cv.Theta())
        )

        dphi = math.degrees(
            wrapped_delta_phi(
                float(p_var.Phi()),
                float(p_cv.Phi()),
            )
        )

        rel_dp = abs(
            float(p_var.R() - p_cv.R()) / float(p_cv.R())
        )

        stats["n"] += 1

        stats["max_dtheta"] = max(
            stats["max_dtheta"], abs(dtheta)
        )
        stats["max_dphi"] = max(
            stats["max_dphi"], abs(dphi)
        )
        stats["max_rel_dp"] = max(
            stats["max_rel_dp"], rel_dp
        )

        stats["sum_dtheta2"] += dtheta * dtheta
        stats["sum_dphi2"] += dphi * dphi

    if debug_muon_angles:
        print("\n========== MUON ANGLE VALIDATION ==========")
        print("MAT: GetThetamu(), GetPhimu()")
        print("New: MuonTheta(), MuonPhi()")
        print("Forward and backward tracks included")
        print("Angles reported in degrees")
        print("===========================================\n")

    if debug_beam_angle:
        print("\n========== BEAM ANGLE VALIDATION ==========")
        print("Comparing beam-angle universes against CV")
        print("Checking angular shifts and |p| preservation")
        print("===========================================\n")

    # ============================================================
    # Standard reconstructed-event loop
    # ============================================================

    for counter in range(nEvents):

        if counter % 10000 == 0:
            print(counter)

            if setAlarm:
                signal.alarm(900)

        # Set CV to the same event before using it as a reference.
        # The normal universe loop will SetEntry again as usual.
        if debug_beam_angle:
            cv_universe.SetEntry(counter)

        for universe in chain.from_iterable(iter(universes.values())):

            universe.SetEntry(counter)
            universe.__dict__.pop(
                "_cached_new_true_l_over_e", None
            )
            universe.ResetWeight()

            if (
                debug_muon_angles
                and universe is cv_universe
            ):
                compare_muon_angles(universe, counter)

            if (
                debug_beam_angle
                and any(universe is u for u in beam_universes)
            ):
                compare_beam_angle(universe, counter)

            if (
                mc
                and AnalysisConfig.skip_2p2h
                and universe.mc_intType == 8
            ):
                continue

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
    # CV muon-angle comparison summary
    # ============================================================

    if debug_muon_angles:
        print("\n========== MUON ANGLE COMPARISON ==========")

        for direction in ("forward", "backward"):
            stats = angle_stats[direction]
            n = stats["n"]

            print("\n---", direction.upper(), "MUONS ---")
            print("Valid comparisons:", n)
            print("Invalid skipped:", stats["skipped"])

            if n == 0:
                continue

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

            print("\nLargest discrepancies:")

            for r in stats["outliers"]:
                print(
                    "\nEntry:", r["entry"],
                    "\n  theta MAT:", r["theta_mat"],
                    "\n  theta new:", r["theta_new"],
                    "\n  dtheta:", r["dtheta"],
                    "\n  phi MAT:", r["phi_mat"],
                    "\n  phi new:", r["phi_new"],
                    "\n  dphi:", r["dphi"],
                    "\n  p det:", r["p_det"],
                    "\n  p beam:", r["p_beam"],
                )

        print("\n===========================================\n")

    # ============================================================
    # Beam-angle systematic summary
    # ============================================================

    if debug_beam_angle:
        print("\n========== BEAM ANGLE SYSTEMATICS ==========")

        if not beam_stats:
            print("WARNING: No beam-angle universes were compared.")
            print("Check AnalysisConfig.exclude_universes.")

        for key, stats in beam_stats.items():
            n = stats["n"]

            if n == 0:
                continue

            print("\nUniverse:", key)
            print("Events:", n)

            print(
                "Delta theta RMS:",
                math.sqrt(stats["sum_dtheta2"] / n),
                "degrees",
            )
            print(
                "Delta phi RMS:",
                math.sqrt(stats["sum_dphi2"] / n),
                "degrees",
            )
            print(
                "Max |Delta theta|:",
                stats["max_dtheta"],
                "degrees",
            )
            print(
                "Max |Delta phi|:",
                stats["max_dphi"],
                "degrees",
            )
            print(
                "Max relative |p| difference:",
                stats["max_rel_dp"],
            )

        print("\n============================================\n")

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
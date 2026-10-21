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
    """The main code of event selection."""

    import array
    from config import PlotConfig

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

    for univ in chain.from_iterable(iter(universes.values())):
        univ.LoadTools(kin_cal, eventClassifier)

    Plots = preparePlots(universes, mc)
    nEvents = chainwrapper.GetEntries()
    debug_prints = 0

    # ============================================================
    # DEBUG: Low-Q2 pion CV effect
    #
    # Run with pure V1 baseline:
    #
    #   USE_MNVTUNE_V1 = True
    #   V1_ADD_DEUTERIUM = False
    #   V1_ADD_ZEXP = False
    #   V1_ADD_LOWQ2PI = False
    #   V1_ADD_COH = False
    #   V1_ADD_LE_COH = False
    #   V1_ADD_DIFFCOH = False
    #
    # LowQ2 is applied manually below, so it must NOT already be
    # included in base_weight.
    # ============================================================

    lowq2_stats = {
        "n_total": 0,
        "n_ccres": 0,
        "n_affected": 0,

        "sum_w_lowq2_ccres": 0.0,
        "sum_w_lowq2_affected": 0.0,

        "yield_base": 0.0,
        "yield_with_lowq2": 0.0,

        "yield_ccres_base": 0.0,
        "yield_ccres_with_lowq2": 0.0,
    }

    # ============================================================
    # EN4 debug histograms
    # ============================================================

    en4_edges = list(PlotConfig.NEUTRINO4_EE_BINNING)
    en4_binning = array.array("d", en4_edges)

    h_en4_base = ROOT.TH1D(
        "h_en4_base_debug",
        ";Energy (GeV);Weighted events",
        len(en4_edges) - 1,
        en4_binning
    )

    h_en4_lowq2 = ROOT.TH1D(
        "h_en4_lowq2_debug",
        ";Energy (GeV);Weighted events",
        len(en4_edges) - 1,
        en4_binning
    )

    h_en4_ccres_base = ROOT.TH1D(
        "h_en4_ccres_base_debug",
        ";Energy (GeV);Weighted CCRES events",
        len(en4_edges) - 1,
        en4_binning
    )

    h_en4_ccres_lowq2 = ROOT.TH1D(
        "h_en4_ccres_lowq2_debug",
        ";Energy (GeV);Weighted CCRES events",
        len(en4_edges) - 1,
        en4_binning
    )

    h_en4_ccres_raw = ROOT.TH1D(
        "h_en4_ccres_raw_debug",
        ";Energy (GeV);CCRES events",
        len(en4_edges) - 1,
        en4_binning
    )

    h_en4_total_raw = ROOT.TH1D(
        "h_en4_total_raw_debug",
        ";Energy (GeV);Selected events",
        len(en4_edges) - 1,
        en4_binning
    )

    # ============================================================
    # TRUE Q2 debug
    # ============================================================

    q2_edges = [
        0.0,
        0.05,
        0.10,
        0.20,
        0.30,
        0.40,
        0.50,
        0.70,
        1.00,
        1.50,
        2.00,
        3.00,
        5.00,
    ]

    q2_stats = []

    for _ in range(len(q2_edges) - 1):
        q2_stats.append({
            "n": 0,
            "sum_w": 0.0,
            "sum_base": 0.0,
            "sum_weighted": 0.0,
            "min_w": None,
            "max_w": None,
        })

    if mc:
        print(
            "[LOWQ2 DEBUG] channel =",
            SystematicsConfig.LowQ2PiWeightChannel
        )

    print("Total number of events RECO: ", {nEvents})

    if AnalysisConfig.testing and nEvents > 5000:
        nEvents = 5000

    print("plotRecoKin, mc ", mc)

    setAlarm = AnalysisConfig.grid

    # ============================================================
    # EVENT LOOP
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
            universe.__dict__.pop("_cached_new_true_l_over_e", None)
            universe.ResetWeight()

            if (
                mc
                and AnalysisConfig.skip_2p2h
                and universe.mc_intType == 8
            ):
                continue

            # Only update kin_cal and eventClassifier for non-vertical universes
            if not universe.IsVerticalOnly():
                kin_cal.CalculateKinematics(universe)
                eventClassifier.Classify(universe)

            # ====================================================
            # LOW-Q2 DEBUG
            #
            # Selected RECO CV MC events only.
            # ====================================================

            if (
                mc
                and universe.ShortName() == "cv"
                and eventClassifier.side_band is not None
            ):

                # Exact EN4 variable from PlotLibrary
                en4 = (
                    kin_cal.reco_E_lep
                    + kin_cal.reco_visE
                )

                base_weight = universe.GetWeight(False)
                w_lowq2 = universe.GetMyLowQ2PiWeight()

                # ------------------------------------------------
                # Inclusive selected sample
                # ------------------------------------------------

                lowq2_stats["n_total"] += 1

                lowq2_stats["yield_base"] += base_weight

                lowq2_stats["yield_with_lowq2"] += (
                    base_weight * w_lowq2
                )

                h_en4_total_raw.Fill(en4)

                h_en4_base.Fill(
                    en4,
                    base_weight
                )

                h_en4_lowq2.Fill(
                    en4,
                    base_weight * w_lowq2
                )

                # ------------------------------------------------
                # CC resonance component
                # ------------------------------------------------

                if (
                    universe.mc_current == 1
                    and universe.mc_intType == 2
                ):

                    lowq2_stats["n_ccres"] += 1
                    lowq2_stats["sum_w_lowq2_ccres"] += w_lowq2

                    lowq2_stats["yield_ccres_base"] += base_weight

                    lowq2_stats["yield_ccres_with_lowq2"] += (
                        base_weight * w_lowq2
                    )

                    h_en4_ccres_raw.Fill(en4)

                    h_en4_ccres_base.Fill(
                        en4,
                        base_weight
                    )

                    h_en4_ccres_lowq2.Fill(
                        en4,
                        base_weight * w_lowq2
                    )

                    if abs(w_lowq2 - 1.0) > 1.0e-6:

                        lowq2_stats["n_affected"] += 1

                        lowq2_stats[
                            "sum_w_lowq2_affected"
                        ] += w_lowq2

                    # --------------------------------------------
                    # TRUE Q2 diagnostic
                    # --------------------------------------------

                    q2 = universe.mc_Q2 * 1.0e-6

                    for iq in range(len(q2_edges) - 1):

                        qlo = q2_edges[iq]
                        qhi = q2_edges[iq + 1]

                        if qlo <= q2 < qhi:

                            s = q2_stats[iq]

                            s["n"] += 1
                            s["sum_w"] += w_lowq2
                            s["sum_base"] += base_weight
                            s["sum_weighted"] += (
                                base_weight * w_lowq2
                            )

                            if (
                                s["min_w"] is None
                                or w_lowq2 < s["min_w"]
                            ):
                                s["min_w"] = w_lowq2

                            if (
                                s["max_w"] is None
                                or w_lowq2 > s["max_w"]
                            ):
                                s["max_w"] = w_lowq2

                            break

            # ====================================================
            # NORMAL PLOT PROCESSING
            # ====================================================

            if (
                eventClassifier.side_band is not None
                or eventClassifier.is_true_signal
            ):

                for entry in Plots:
                    entry.Process(universe)

    # ============================================================
    # GLOBAL LOW-Q2 SUMMARY
    # ============================================================

    if mc:

        print("")
        print("============================================================")
        print("LOW-Q2 PION DEBUG: SELECTED CV EVENTS")
        print("============================================================")

        print(
            "channel                    =",
            SystematicsConfig.LowQ2PiWeightChannel
        )

        print(
            "selected events            =",
            lowq2_stats["n_total"]
        )

        print(
            "selected CCRES events      =",
            lowq2_stats["n_ccres"]
        )

        print(
            "events with LowQ2 w != 1   =",
            lowq2_stats["n_affected"]
        )

        if lowq2_stats["n_total"] > 0:

            print(
                "raw CCRES fraction         = {:.3f}%".format(
                    100.0
                    * lowq2_stats["n_ccres"]
                    / lowq2_stats["n_total"]
                )
            )

        if lowq2_stats["n_ccres"] > 0:

            print(
                "mean LowQ2 weight CCRES    = {:.6f}".format(
                    lowq2_stats["sum_w_lowq2_ccres"]
                    / lowq2_stats["n_ccres"]
                )
            )

        if lowq2_stats["n_affected"] > 0:

            print(
                "mean weight affected       = {:.6f}".format(
                    lowq2_stats["sum_w_lowq2_affected"]
                    / lowq2_stats["n_affected"]
                )
            )

        print("")
        print("----- weighted yield effect -----")

        print(
            "V1 baseline yield           = {:.6f}".format(
                lowq2_stats["yield_base"]
            )
        )

        print(
            "V1 + LowQ2 yield            = {:.6f}".format(
                lowq2_stats["yield_with_lowq2"]
            )
        )

        if lowq2_stats["yield_base"] != 0:

            print(
                "inclusive LowQ2 shift       = {:+.3f}%".format(
                    (
                        lowq2_stats["yield_with_lowq2"]
                        / lowq2_stats["yield_base"]
                        - 1.0
                    ) * 100.0
                )
            )

        print(
            "CCRES baseline yield        = {:.6f}".format(
                lowq2_stats["yield_ccres_base"]
            )
        )

        print(
            "CCRES + LowQ2 yield         = {:.6f}".format(
                lowq2_stats["yield_ccres_with_lowq2"]
            )
        )

        if lowq2_stats["yield_ccres_base"] != 0:

            print(
                "CCRES LowQ2 shift           = {:+.3f}%".format(
                    (
                        lowq2_stats["yield_ccres_with_lowq2"]
                        / lowq2_stats["yield_ccres_base"]
                        - 1.0
                    ) * 100.0
                )
            )

        if lowq2_stats["yield_base"] != 0:

            print(
                "weighted CCRES fraction     = {:.3f}%".format(
                    100.0
                    * lowq2_stats["yield_ccres_base"]
                    / lowq2_stats["yield_base"]
                )
            )

        # ========================================================
        # EN4 BIN-BY-BIN SUMMARY
        # ========================================================

        print("")
        print("============================================================")
        print("LOW-Q2 PION DEBUG: EN4 BIN-BY-BIN")
        print("============================================================")

        print(
            "{:<14s} {:>10s} {:>12s} {:>12s} "
            "{:>13s} {:>14s} {:>14s}".format(
                "EN4 bin",
                "raw N",
                "V1",
                "V1+LowQ2",
                "shift [%]",
                "CCRES frac [%]",
                "<w> CCRES"
            )
        )

        print("-" * 100)

        for i in range(1, h_en4_base.GetNbinsX() + 1):

            lo = h_en4_base.GetXaxis().GetBinLowEdge(i)
            hi = h_en4_base.GetXaxis().GetBinUpEdge(i)

            raw_n = h_en4_total_raw.GetBinContent(i)

            base = h_en4_base.GetBinContent(i)
            lowq2 = h_en4_lowq2.GetBinContent(i)

            ccres_base = h_en4_ccres_base.GetBinContent(i)
            ccres_lowq2 = h_en4_ccres_lowq2.GetBinContent(i)

            shift = (
                100.0 * (lowq2 / base - 1.0)
                if base != 0.0
                else 0.0
            )

            ccres_frac = (
                100.0 * ccres_base / base
                if base != 0.0
                else 0.0
            )

            mean_w_ccres = (
                ccres_lowq2 / ccres_base
                if ccres_base != 0.0
                else 1.0
            )

            print(
                "{:5.2f}-{:<6.2f} "
                "{:10.0f} "
                "{:12.3f} "
                "{:12.3f} "
                "{:13.3f} "
                "{:14.3f} "
                "{:14.4f}".format(
                    lo,
                    hi,
                    raw_n,
                    base,
                    lowq2,
                    shift,
                    ccres_frac,
                    mean_w_ccres
                )
            )

        # ========================================================
        # EN4 INTEGRATED CHECK
        # ========================================================

        en4_base_integral = h_en4_base.Integral(
            0,
            h_en4_base.GetNbinsX() + 1
        )

        en4_lowq2_integral = h_en4_lowq2.Integral(
            0,
            h_en4_lowq2.GetNbinsX() + 1
        )

        print("")
        print("----- EN4 integrated check -----")

        print(
            "EN4 V1 integral             = {:.6f}".format(
                en4_base_integral
            )
        )

        print(
            "EN4 V1 + LowQ2 integral     = {:.6f}".format(
                en4_lowq2_integral
            )
        )

        if en4_base_integral != 0.0:

            print(
                "EN4 LowQ2 shift             = {:+.3f}%".format(
                    (
                        en4_lowq2_integral
                        / en4_base_integral
                        - 1.0
                    ) * 100.0
                )
            )

        # ========================================================
        # TRUE Q2 SUMMARY
        # ========================================================

        print("")
        print("============================================================")
        print("LOW-Q2 PION DEBUG: TRUE Q2 FOR SELECTED CCRES")
        print("============================================================")

        print(
            "{:<14s} {:>8s} {:>12s} {:>14s} "
            "{:>12s} {:>12s}".format(
                "Q2 [GeV2]",
                "N",
                "<w> raw",
                "<w> weighted",
                "min w",
                "max w"
            )
        )

        print("-" * 80)

        for iq, s in enumerate(q2_stats):

            qlo = q2_edges[iq]
            qhi = q2_edges[iq + 1]

            if s["n"] == 0:
                continue

            mean_raw = (
                s["sum_w"]
                / s["n"]
            )

            mean_weighted = (
                s["sum_weighted"]
                / s["sum_base"]
                if s["sum_base"] != 0.0
                else 1.0
            )

            print(
                "{:5.2f}-{:<6.2f} "
                "{:8d} "
                "{:12.4f} "
                "{:14.4f} "
                "{:12.4f} "
                "{:12.4f}".format(
                    qlo,
                    qhi,
                    s["n"],
                    mean_raw,
                    mean_weighted,
                    s["min_w"],
                    s["max_w"]
                )
            )

        print("============================================================")
        print("")

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

        if st=="mc" and Truth:
            print("selecting truth")
            plotTruthKin(Utilities.fileChain(AnalysisConfig.playlist,"mc",AnalysisConfig.ntuple_tag,"Truth",AnalysisConfig.count[0],AnalysisConfig.count[1]),output_file)
        if Reco :
            print("selecting reco")
            plotRecoKin(st=="mc", Utilities.fileChain(AnalysisConfig.playlist,st,AnalysisConfig.ntuple_tag,None,AnalysisConfig.count[0],AnalysisConfig.count[1]), output_file)
        output_file.Close()
        print("selection is done for ", st, AnalysisConfig.playlist)
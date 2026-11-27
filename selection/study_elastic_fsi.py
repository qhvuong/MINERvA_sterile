#!/usr/bin/env python3
"""Run from the analysis repository, with its normal configuration activated.

Only this script and tools/ElasticFSIDiagnostics.py are added. No production
histograms are opened or overwritten. Diagnostic outputs use a new directory.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from array import array


def dump_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def dump_csv(path, rows):
    with path.open('w', newline='') as output:
        if rows:
            writer = csv.DictWriter(output, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--fsi-output', required=True)
    parser.add_argument('--fsi-top', type=int, default=10)
    parser.add_argument('--fsi-max-events', type=int, default=None)
    parser.add_argument('--fsi-verbose-events', type=int, default=6)
    opts, remaining = parser.parse_known_args()
    if opts.fsi_top < 1 or opts.fsi_verbose_events < 0 or (
        opts.fsi_max_events is not None and opts.fsi_max_events < 1
    ):
        parser.error('Invalid diagnostic event limit')
    # The production configuration parses sys.argv at import time.
    sys.argv = [sys.argv[0]] + remaining
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import ROOT
    from config.AnalysisConfig import AnalysisConfig
    from config import SystematicsConfig
    from config.CutConfig import SAMPLE_CUTS, KINEMATICS_CUTS
    from tools import SystematicsUniverse as SU, Utilities
    from tools.EventClassification import EventClassifier
    from tools.KinematicsCalculator import KinematicsCalculator
    from tools.PlotLibrary import PLOT_SETTINGS
    from tools.ElasticFSIDiagnostics import Accounting

    if AnalysisConfig.is_pc or SystematicsConfig.USE_SWAPPED:
        raise RuntimeError('Use the ordinary, unswapped MC sample for this study')
    if abs(int(SystematicsConfig.AnaNuPDG)) != 14:
        raise RuntimeError('Activate the CCnumu or CCnumubar analysis configuration')
    if not (SU.USE_MNVTUNE_V4P3P0 or SU.USE_MNVTUNE_V4P3P1):
        raise RuntimeError('Activate v4.3.0 or v4.3.1; other tune settings are incompatible')
    out = Path(opts.fsi_output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    status = dict(status='running', argv=remaining, playlist=AnalysisConfig.playlist,
                  count=list(AnalysisConfig.count), POT_normalized=False,
                  source_sha256=hashlib.sha256(Path(SU.__file__).read_bytes()).hexdigest(),
                  source=SU.__file__, active_tune='v4.3.1' if SU.USE_MNVTUNE_V4P3P1 else 'v4.3.0',
                  AnaNuPDG=int(SystematicsConfig.AnaNuPDG),
                  signal_cuts=list(SAMPLE_CUTS['Signal']),
                  kinematics_cuts=list(KINEMATICS_CUTS),
                  near_one_absolute_tolerance=1e-8)
    dump_json(out / 'audit.json', status)
    try:
        settings = PLOT_SETTINGS['Biased Neutrino Energy']
        if settings['name'] != 'EN4' or settings.get('cuts') or settings.get('weight_function'):
            raise RuntimeError('EN4 definition has changed; inspect custom cuts/weights first')
        edges = list(map(float, settings['binning'][0]))
        energy_getter = settings['value_getter'][0]
        status['EN4_edges_GeV'] = edges
        chain = Utilities.fileChain(AnalysisConfig.playlist, 'mc', AnalysisConfig.ntuple_tag,
                                    None, AnalysisConfig.count[0], AnalysisConfig.count[1])
        status['input_files'] = [str(x.GetTitle()) for x in chain.GetChain().GetListOfFiles()]
        status['tune_flags'] = {k: v for k, v in vars(SU).items()
                               if k.startswith(('USE_MNVTUNE_', 'V4P3P0_'))
                               and isinstance(v, (str, bool, int, float, type(None)))}
        dump_json(out / 'audit.json', status)
        universes = SU.GetAllSystematicsUniverses(chain, False, False, exclude=['all'])
        if set(universes) != {'cv'} or len(universes['cv']) != 1:
            raise RuntimeError('Expected exactly one CV universe')
        universe = universes['cv'][0]
        kin = KinematicsCalculator(correct_beam_angle=True, correct_MC_energy_scale=False,
                                   calc_true=True, is_pc=False)
        # Signal-only classification uses precisely the production Signal cuts.
        classifier = EventClassifier(classifiers=['Reco', 'Truth'], use_kin_cuts=True,
                                     use_sideband=[])
        universe.LoadTools(kin, classifier)
        # Make calculator type visible to PyROOT if it was not dictionary-exported.
        if not hasattr(ROOT.PlotUtils, 'weight_fsi'):
            if not ROOT.gInterpreter.Declare('#include "PlotUtils/weight_fsi.h"'):
                raise RuntimeError('Cannot expose PlotUtils::weight_fsi to PyROOT')
        independent = ROOT.PlotUtils.weight_fsi()
        independent.UseTrackingThreshold()
        vector_int = ('mc_er_ID', 'mc_er_status', 'mc_er_FD', 'mc_er_LD', 'mc_er_mother')
        vector_double = ('mc_er_Px', 'mc_er_Py', 'mc_er_Pz', 'mc_er_E')
        scalars = ('mc_incoming', 'mc_primaryLepton', 'mc_charm', 'mc_intType',
                   'mc_targetA', 'mc_targetZ', 'mc_resID', 'mc_er_nPart')

        def calculate(calculator, entry, verbose=False):
            args = [universe.GetInt(name) for name in scalars]
            args += [universe.GetVecInt(name) for name in vector_int]
            args += [universe.GetVecDouble(name) for name in vector_double]
            calculator.calcWeights(*(args + [entry, verbose]))
            return float(calculator.GetElasticWeight(1))

        ROOT.TH1.AddDirectory(False)
        histograms = {}
        def fill(category, energy, baseline, corrected):
            for label in ('TOTAL', category):
                for tag, weight in (('430', baseline), ('431', corrected)):
                    key = (label, tag)
                    if key not in histograms:
                        name = 'EN4_{}_{}'.format(label, tag)
                        h = ROOT.TH1D(name, ';Energy (GeV);Weighted events', len(edges)-1, array('d', edges))
                        h.Sumw2()
                        histograms[key] = h
                    histograms[key].Fill(energy, weight)

        accounting = Accounting(edges, opts.fsi_top)
        n = int(chain.GetEntries())
        if AnalysisConfig.testing:
            n = min(n, 2000)
        if opts.fsi_max_events is not None:
            n = min(n, opts.fsi_max_events)
        selected = 0
        max_fsi_difference = 0.
        max_nominal_difference = 0.
        for entry in range(n):
            if entry % 10000 == 0:
                print('[FSI] entry {}/{} selected {}'.format(entry, n, selected), flush=True)
            universe.SetEntry(entry)
            universe.ResetWeight()
            universe.__dict__.pop('_cached_new_true_l_over_e', None)
            if AnalysisConfig.skip_2p2h and universe.mc_intType == 8:
                continue
            kin.CalculateKinematics(universe)
            classifier.Classify(universe)
            if classifier.side_band != 'Signal':
                continue
            energy = float(energy_getter(universe))
            baseline = float(universe.GetStandardWeightV4p3p0())
            elastic = float(universe.GetFSIWeight(1))
            reference = calculate(independent, entry)
            corrected = baseline * elastic
            nominal = float(universe.GetWeight(False))
            expected = corrected if SU.USE_MNVTUNE_V4P3P1 else baseline
            if not all(math.isfinite(x) for x in (energy, baseline, elastic, reference, corrected, nominal)):
                raise RuntimeError('Nonfinite value at chain entry {}'.format(entry))
            max_fsi_difference = max(max_fsi_difference, abs(reference - elastic))
            max_nominal_difference = max(max_nominal_difference, abs(nominal - expected))
            if not math.isclose(reference, elastic, rel_tol=1e-10, abs_tol=1e-12):
                raise RuntimeError('FSI wrapper/calculator mismatch at entry {}: {} vs {}'.format(entry, elastic, reference))
            if not math.isclose(nominal, expected, rel_tol=1e-10, abs_tol=1e-12):
                raise RuntimeError('Nominal weight mismatch at entry {}'.format(entry))
            context = dict(chain_entry=entry, mc_incoming=int(universe.mc_incoming),
                           mc_intType=int(universe.mc_intType), mc_current=int(universe.mc_current),
                           true_Enu_GeV=float(universe.mc_incomingE)*1e-3)
            # Entry plus the exact input file list provides an unambiguous replay locator.
            accounting.add(str(classifier.truth_class), energy, baseline, elastic, context)
            fill(str(classifier.truth_class), energy, baseline, corrected)
            selected += 1

        rows = accounting.rows()
        totals = accounting.totals()
        dump_csv(out / 'weight_groups.csv', rows)
        dump_csv(out / 'integrated_yields.csv', totals)
        top = accounting.top_events()
        zeros = [event for group in accounting.zeros.values() for event in group]
        # Save complete GENIE records only for the bounded set of retained examples.
        examples = {e['chain_entry']: e for group in top for e in group['events']}
        examples.update({e['chain_entry']: e for e in zeros})
        for entry, event in sorted(examples.items()):
            universe.SetEntry(entry)
            event['genie_record'] = {name: int(universe.GetInt(name)) for name in scalars}
            event['genie_record'].update({name: list(map(int, universe.GetVecInt(name))) for name in vector_int})
            event['genie_record'].update({name: list(map(float, universe.GetVecDouble(name))) for name in vector_double})
        dump_json(out / 'top_events.json', dict(rankings=top, zero_examples=zeros))
        rootfile = ROOT.TFile.Open(str(out / 'EN4_comparison.root'), 'RECREATE')
        if not rootfile or rootfile.IsZombie():
            raise RuntimeError('Cannot create diagnostic ROOT file')
        for h in histograms.values():
            h.Write()
        rootfile.Close()
        status.update(status='complete', processed_entries=n, selected_events=selected,
                      skip_2p2h=bool(AnalysisConfig.skip_2p2h),
                      max_wrapper_calculator_difference=max_fsi_difference,
                      max_nominal_weight_difference=max_nominal_difference,
                      input_files=[str(x.GetTitle()) for x in chain.GetChain().GetListOfFiles()])
        dump_json(out / 'audit.json', status)
        print('\n[FSI] REGULAR EN4 BINS (raw MC weights; no POT scaling)', flush=True)
        for row in totals:
            if row['scope'] == 'regular':
                print(row, flush=True)
        # Show influential DIS examples first, plus zero-weight DIS examples.
        candidates = sorted(examples.values(), key=lambda e: (
            'DIS' in e['category'], e['region'] == 'regular', abs(e['delta'])), reverse=True)
        chosen = candidates[:opts.fsi_verbose_events]
        if opts.fsi_verbose_events >= 2:
            zero_dis = [e for e in zeros if 'DIS' in e['category'] and e['region'] == 'regular']
            if zero_dis and not any(e['elastic'] == 0 for e in chosen):
                chosen = chosen[:-1] + zero_dis[:1]
        for event in chosen:
            universe.SetEntry(event['chain_entry'])
            print('\n[FSI VERBOSE] {}'.format({k: v for k, v in event.items() if k != 'genie_record'}), flush=True)
            # Fresh object prevents an entry-cache hit during replay.
            replay = ROOT.PlotUtils.weight_fsi()
            replay.UseTrackingThreshold()
            calculate(replay, event['chain_entry'], verbose=True)
        print('\n[FSI] Outputs: {}'.format(out), flush=True)
    except Exception as exc:
        status.update(status='failed', error=str(exc))
        dump_json(out / 'audit.json', status)
        raise


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Before/after elastic FSI study, using one truth-defined cohort and common weights."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys


def json_out(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def csv_out(path, rows):
    with path.open('w', newline='') as f:
        if rows:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--fsi-output', required=True)
    parser.add_argument('--fsi-max-events', type=int)
    opts, remaining = parser.parse_known_args()
    if opts.fsi_max_events is not None and opts.fsi_max_events < 1:
        parser.error('--fsi-max-events must be positive')
    sys.argv = [sys.argv[0]] + remaining
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import ROOT
    from config.AnalysisConfig import AnalysisConfig as AC
    from config import SystematicsConfig as SC
    from config.SignalDef import IsFiducial
    from config.CutConfig import SAMPLE_CUTS, KINEMATICS_CUTS
    from tools import SystematicsUniverse as SU, Utilities
    from tools.EventClassification import EventClassifier
    from tools.KinematicsCalculator import KinematicsCalculator
    from tools.PlotLibrary import PLOT_SETTINGS
    from tools.ElasticFSISelection import (SelectionAccounting, generator_flux_function,
                                          NEUTRON_EDGES, TRUE_ENERGY_EDGES)
    if AC.is_pc or SC.USE_SWAPPED or abs(int(SC.AnaNuPDG)) != 14:
        raise RuntimeError('Use the ordinary unswapped CCnumu/CCnumubar MC config')
    if not (SU.USE_MNVTUNE_V4P3P0 or SU.USE_MNVTUNE_V4P3P1):
        raise RuntimeError('Activate v4.3.0 or v4.3.1')
    if AC.skip_2p2h:
        raise RuntimeError('Disable skip_2p2h for this study')
    out = Path(opts.fsi_output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    source = Path(SU.__file__).read_text()
    generator_weight, transformed_source = generator_flux_function(source, vars(SU))
    cuts = list(SAMPLE_CUTS['Signal'])
    cuts += ['Reco'+name for name in KINEMATICS_CUTS]
    stages = ['00_truth_fiducial'] + ['{:02d}_{}'.format(i+1, cut) for i, cut in enumerate(cuts)]
    stages += ['{:02d}_EN4_regular'.format(len(cuts)+1)]
    status = dict(status='running', diagnostic_version=2, argv=remaining,
                  playlist=AC.playlist, count=list(AC.count), AnaNuPDG=int(SC.AnaNuPDG),
                  cohort='right-sign CC, config.SignalDef.IsFiducial, no reconstructed cuts',
                  weight_basis='production V430 generator and flux; GEANT and MINOS factors removed',
                  POT_normalized=False, stages=stages, signal_cuts=cuts,
                  source_sha256=hashlib.sha256(source.encode()).hexdigest(),
                  derived_weight_ast_sha256=hashlib.sha256(transformed_source.encode()).hexdigest(),
                  neutron_split='minimum status-14 fate-1 neutron KE; candidate, not necessarily a contributor',
                  neutron_edges_MeV=list(NEUTRON_EDGES), true_Enu_edges_GeV=list(TRUE_ENERGY_EDGES),
                  limitation='denominator is events present in the reco ntuple, not a separate Truth tree')
    json_out(out/'audit.json', status)
    try:
        chain = Utilities.fileChain(AC.playlist, 'mc', AC.ntuple_tag, None, AC.count[0], AC.count[1])
        status['input_files'] = [str(f.GetTitle()) for f in chain.GetChain().GetListOfFiles()]
        status['tune_flags'] = {k:v for k,v in vars(SU).items()
                               if k.startswith(('USE_MNVTUNE_', 'V4P3P0_'))
                               and isinstance(v, (str, bool, int, float, type(None)))}
        universes = SU.GetAllSystematicsUniverses(chain, False, False, exclude=['all'])
        if set(universes) != {'cv'} or len(universes['cv']) != 1:
            raise RuntimeError('Expected only one CV universe')
        u = universes['cv'][0]
        kin = KinematicsCalculator(correct_beam_angle=True, correct_MC_energy_scale=False,
                                   calc_true=True, is_pc=False)
        classifier = EventClassifier(classifiers=['Reco','Truth'], use_kin_cuts=True, use_sideband=[])
        u.LoadTools(kin, classifier)
        if not hasattr(ROOT.PlotUtils, 'fate_code'):
            if not ROOT.gInterpreter.Declare('#include "PlotUtils/weight_fsi.h"'):
                raise RuntimeError('Cannot expose FSI classes')
        calc = ROOT.PlotUtils.weight_fsi()
        calc.UseTrackingThreshold()
        fates = ROOT.PlotUtils.fate_code()
        # Public calculator counters give the exact contribution multiplicities.
        getattr(calc, '__count_elastic_pion')
        getattr(calc, '__count_elastic_nucleon')
        scalars = ('mc_incoming','mc_primaryLepton','mc_charm','mc_intType',
                   'mc_targetA','mc_targetZ','mc_resID','mc_er_nPart')
        ivectors = ('mc_er_ID','mc_er_status','mc_er_FD','mc_er_LD','mc_er_mother')
        dvectors = ('mc_er_Px','mc_er_Py','mc_er_Pz','mc_er_E')
        en4settings = PLOT_SETTINGS['Biased Neutrino Energy']
        edges = list(map(float, en4settings['binning'][0]))
        status['EN4_edges_GeV'] = edges
        json_out(out/'audit.json', status)
        accounting = SelectionAccounting()
        production = {}
        max_fsi_diff = max_weight_diff = 0.
        processed = int(chain.GetEntries())
        if AC.testing:
            processed = min(processed, 2000)
        if opts.fsi_max_events is not None:
            processed = min(processed, opts.fsi_max_events)
        cohort = selected = regular = 0
        names = {1:'QE',2:'Delta',3:'DIS',4:'Coherent',8:'2p2h'}
        for entry in range(processed):
            if entry % 10000 == 0:
                print('[FSI selection] {}/{} cohort={} selected={}'.format(entry,processed,cohort,selected),flush=True)
            u.SetEntry(entry)
            u.ResetWeight()
            u.__dict__.pop('_cached_new_true_l_over_e', None)
            if u.mc_current != 1 or u.mc_incoming != int(SC.AnaNuPDG) or not IsFiducial(u):
                continue
            cohort += 1
            baseline = float(generator_weight(u))
            elastic = float(u.GetFSIWeight(1))
            integers = [u.GetInt(name) for name in scalars]
            vectors = [u.GetVecInt(name) for name in ivectors] + [u.GetVecDouble(name) for name in dvectors]
            calc.calcWeights(*(integers+vectors+[entry,False]))
            reference = float(calc.GetElasticWeight(1))
            if not all(math.isfinite(x) for x in (baseline,elastic,reference)):
                raise RuntimeError('Nonfinite weight at entry {}'.format(entry))
            max_fsi_diff = max(max_fsi_diff, abs(elastic-reference))
            if not math.isclose(elastic,reference,rel_tol=1e-10,abs_tol=1e-12):
                raise RuntimeError('FSI calculator mismatch at entry {}'.format(entry))
            npi = int(getattr(calc, '__count_elastic_pion'))
            nnuc = int(getattr(calc, '__count_elastic_nucleon'))
            # fate_code takes the same scalar list except mc_resID.
            fates.calcFates(*(integers[:6]+integers[7:]+vectors+[entry,False]))
            neutron_energies = [float(getattr(p,'__ke')) for p in fates.GetIntranukeParticles()
                                if int(getattr(p,'__pdg')) == 2112 and int(getattr(p,'__fate')) == 1]
            min_neutron = min(neutron_energies) if neutron_energies else None
            kin.CalculateKinematics(u)
            classifier.Classify(u)
            passed_stages = [stages[0]]
            passed = True
            for i, cut in enumerate(cuts):
                current = (not classifier.reco_cuts_passed[cut[7:]] if cut.startswith('Inverse')
                           else classifier.reco_cuts_passed[cut])
                passed = passed and bool(current)
                if passed:
                    passed_stages.append(stages[i+1])
            if passed != (classifier.side_band == 'Signal'):
                raise RuntimeError('Cut-flow/production selection mismatch at entry {}'.format(entry))
            category = names.get(int(u.mc_intType), 'Other')
            if int(u.mc_charm):
                category += '_charm'
            if passed:
                selected += 1
                en4 = float(en4settings['value_getter'][0](u))
                if not math.isfinite(en4):
                    raise RuntimeError('Nonfinite selected EN4 at entry {}'.format(entry))
                detector = 1.
                if SU.V4P3P0_ADD_GEANT_HADRON:
                    detector *= float(u.GetGeantHadronWeight())
                if SU.V4P3P0_ADD_MINOS_EFFICIENCY:
                    detector *= float(u.GetMyMinosEfficiencyWeight())
                full430 = float(u.GetStandardWeightV4p3p0())
                recomposed = baseline*detector
                if not math.isfinite(recomposed) or not math.isfinite(full430):
                    raise RuntimeError('Nonfinite detector weight')
                max_weight_diff = max(max_weight_diff,abs(full430-recomposed))
                if not math.isclose(full430,recomposed,rel_tol=1e-10,abs_tol=1e-12):
                    raise RuntimeError('Generator/detector factorization mismatch at entry {}'.format(entry))
                if edges[0] <= en4 < edges[-1]:
                    passed_stages.append(stages[-1])
                    regular += 1
                    key = str(classifier.truth_class)
                    row = production.setdefault(key,[0,0.,0.])
                    row[0] += 1
                    row[1] += full430
                    row[2] += full430*elastic
            accounting.add(category,passed_stages,baseline,elastic,npi,nnuc,min_neutron,
                           float(u.mc_incomingE)*1e-3)
        rows = accounting.rows()
        csv_out(out/'stage_totals.csv',[r for r in rows if r['split']=='all'])
        csv_out(out/'stage_splits.csv',[r for r in rows if r['split']!='all'])
        csv_out(out/'selection_efficiencies.csv',accounting.efficiencies(stages[0],stages[-1]))
        csv_out(out/'selected_production_weights.csv',[
            dict(category=c,events=v[0],sum_w430=v[1],sum_w431=v[2],ratio=v[2]/v[1] if v[1] else None)
            for c,v in sorted(production.items())])
        status.update(status='complete',processed_entries=processed,truth_cohort_events=cohort,
                      selected_cohort_events=selected,regular_selected_cohort_events=regular,
                      max_wrapper_calculator_difference=max_fsi_diff,
                      max_generator_detector_factorization_difference=max_weight_diff)
        json_out(out/'audit.json',status)
        print('\n[FSI selection] COMMON generator+flux weights',flush=True)
        if not cohort or not selected:
            print('[FSI selection] Empty cohort or selected sample: no before/after validation possible',flush=True)
        for r in rows:
            if r['split']=='all' and r['category'] in ('TOTAL','DIS','2p2h'):
                print(r,flush=True)
        print('[FSI selection] Outputs: {}'.format(out),flush=True)
    except Exception as exc:
        status.update(status='failed',error=str(exc))
        json_out(out/'audit.json',status)
        raise


if __name__ == '__main__':
    main()

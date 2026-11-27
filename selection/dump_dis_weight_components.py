#!/usr/bin/env python3
"""Run once per tune in separate processes; retain every selected muon DIS event."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys


def write_audit(path, audit):
    path.write_text(json.dumps(audit,indent=2,allow_nan=False)+'\n')


def main():
    p=argparse.ArgumentParser(add_help=False)
    p.add_argument('--pair-tune',required=True,choices=['v4p3p0','v4p3p1'])
    p.add_argument('--pair-output',required=True)
    p.add_argument('--pair-max-events',type=int)
    opts,remaining=p.parse_known_args()
    if opts.pair_max_events is not None and opts.pair_max_events<1:
        p.error('--pair-max-events must be positive')
    sys.argv=[sys.argv[0]]+remaining
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
    from config.AnalysisConfig import AnalysisConfig as AC
    from config import SystematicsConfig as SC
    from config.CutConfig import SAMPLE_CUTS,KINEMATICS_CUTS
    from tools import SystematicsUniverse as SU,Utilities
    from tools.KinematicsCalculator import KinematicsCalculator
    from tools.EventClassification import EventClassifier
    from tools.PlotLibrary import PLOT_SETTINGS
    from tools.FSIEventPair import WeightTrace,close
    if AC.is_pc or SC.USE_SWAPPED or abs(int(SC.AnaNuPDG))!=14 or AC.skip_2p2h:
        raise RuntimeError('Use ordinary unswapped CCnumu/CCnumubar MC; skip_2p2h must be off')
    # Configure the actual tune before universe creation/static GENIE setup.
    # Each command starts a new process; no source files are modified.
    SU.USE_MNVTUNE_V1=False
    SU.USE_MNVTUNE_V4P3P0=opts.pair_tune=='v4p3p0'
    SU.USE_MNVTUNE_V4P3P1=opts.pair_tune=='v4p3p1'
    SU.USE_MNVTUNE_V4P4P1=False
    source=Path(SU.__file__).read_text()
    trace=WeightTrace(source,vars(SU))
    out=Path(opts.pair_output).resolve()
    out.mkdir(parents=True,exist_ok=False)
    audit=dict(status='running',tune=opts.pair_tune,argv=remaining,
               playlist=AC.playlist,count=list(AC.count),AnaNuPDG=int(SC.AnaNuPDG),
               signal_cuts=list(SAMPLE_CUTS['Signal']),kinematics_cuts=list(KINEMATICS_CUTS),
               source_sha256=hashlib.sha256(source.encode()).hexdigest(),
               method_sha256=trace.method_sha256,component_labels=trace.labels,
               component_expressions=trace.expressions,
               shared_tune_flags={k:v for k,v in vars(SU).items() if k.startswith('V4P3P0_')
                                  and isinstance(v,(bool,int,float,str,type(None)))},
               sample='production Signal selection; CC muon-flavor DIS of either sign; no added truth-fiducial cuts',
               event_key='global chain_entry, valid only with identical ordered input files',
               POT_normalized=False)
    write_audit(out/'audit.json',audit)
    try:
        chain=Utilities.fileChain(AC.playlist,'mc',AC.ntuple_tag,None,AC.count[0],AC.count[1])
        audit['input_files']=[str(f.GetTitle()) for f in chain.GetChain().GetListOfFiles()]
        universes=SU.GetAllSystematicsUniverses(chain,False,False,exclude=['all'])
        if set(universes)!={'cv'} or len(universes['cv'])!=1:
            raise RuntimeError('Expected exactly one CV universe')
        u=universes['cv'][0]
        kin=KinematicsCalculator(correct_beam_angle=True,correct_MC_energy_scale=False,calc_true=True,is_pc=False)
        classifier=EventClassifier(classifiers=['Reco','Truth'],use_kin_cuts=True,use_sideband=[])
        u.LoadTools(kin,classifier)
        settings=PLOT_SETTINGS['Biased Neutrino Energy']
        if settings['name']!='EN4' or settings.get('cuts') or settings.get('weight_function'):
            raise RuntimeError('Inspect changed EN4 plot definition before running')
        edges=list(map(float,settings['binning'][0]))
        audit['EN4_edges_GeV']=edges
        n=int(chain.GetEntries())
        if AC.testing:n=min(n,2000)
        if opts.pair_max_events is not None:n=min(n,opts.pair_max_events)
        audit['processed_entries']=n
        write_audit(out/'audit.json',audit)
        fields=['chain_entry','truth_class','mc_incoming','mc_intType','true_Enu_GeV','EN4_GeV','region',
                'active_total','traced_baseline_total','fsi_evaluated','fsi_applied','reconstructed_total']
        for label in trace.labels:fields += [label,label+'_applied']
        selected=0
        max_difference=0.
        with (out/'events.csv').open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
            for entry in range(n):
                if entry%10000==0:
                    print('[DIS weights] tune={} entry={}/{} selected={}'.format(opts.pair_tune,entry,n,selected),flush=True)
                u.SetEntry(entry);u.ResetWeight()
                u.__dict__.pop('_cached_new_true_l_over_e',None)
                if u.mc_current!=1 or u.mc_intType!=3 or abs(int(u.mc_incoming))!=14:
                    continue
                kin.CalculateKinematics(u);classifier.Classify(u)
                if classifier.side_band!='Signal':continue
                # Actual active production weight is evaluated independently of tracing.
                active=float(u.GetWeight(False))
                baseline,components,applied=trace.evaluate(u)
                elastic=float(u.GetFSIWeight(1))
                fsi_applied=elastic if opts.pair_tune=='v4p3p1' else 1.
                reconstructed=baseline*fsi_applied
                energy=float(settings['value_getter'][0](u))
                if not all(math.isfinite(x) for x in (active,baseline,elastic,reconstructed,energy)):
                    raise RuntimeError('Nonfinite event at entry {}'.format(entry))
                if not close(active,reconstructed):
                    raise RuntimeError('Component/active-total mismatch at entry {}'.format(entry))
                max_difference=max(max_difference,abs(active-reconstructed))
                row=dict(chain_entry=entry,truth_class=str(classifier.truth_class),mc_incoming=int(u.mc_incoming),
                         mc_intType=int(u.mc_intType),true_Enu_GeV=float(u.mc_incomingE)*1e-3,EN4_GeV=energy,
                         region='underflow' if energy<edges[0] else 'overflow' if energy>=edges[-1] else 'regular',
                         active_total=active,traced_baseline_total=baseline,fsi_evaluated=elastic,
                         fsi_applied=fsi_applied,reconstructed_total=reconstructed)
                row.update(components)
                row.update({label+'_applied':int(label in applied) for label in trace.labels})
                writer.writerow(row);selected+=1
        audit.update(status='complete',selected_DIS_events=selected,max_component_total_difference=max_difference)
        write_audit(out/'audit.json',audit)
        print('[DIS weights] complete: {} events in {}'.format(selected,out),flush=True)
    except Exception as exc:
        audit.update(status='failed',error=str(exc));write_audit(out/'audit.json',audit)
        raise


if __name__=='__main__':main()

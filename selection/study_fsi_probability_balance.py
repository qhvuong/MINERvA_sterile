#!/usr/bin/env python3
"""Measure single-hadron compensation and its combination into DIS event weights."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys


def write_json(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def write_csv(path,rows):
    with path.open('w',newline='') as f:
        if rows:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main():
    parser=argparse.ArgumentParser(add_help=False)
    parser.add_argument('--balance-output',required=True)
    parser.add_argument('--balance-tree',choices=['reco','truth'],required=True)
    parser.add_argument('--balance-source',default='/exp/minerva/app/users/qvuong/MAT_AL9/MAT-MINERvA/weighters/weight_fsi.cxx')
    parser.add_argument('--balance-header',default=None)
    parser.add_argument('--balance-max-events',type=int)
    opts,remaining=parser.parse_known_args()
    if opts.balance_max_events is not None and opts.balance_max_events<1:parser.error('Invalid event limit')
    sys.argv=[sys.argv[0]]+remaining
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
    import ROOT
    from config.AnalysisConfig import AnalysisConfig as AC
    from config import SystematicsConfig as SC
    from config.SignalDef import IsFiducial
    from config.CutConfig import SAMPLE_CUTS,KINEMATICS_CUTS
    from tools import SystematicsUniverse as SU,Utilities
    from tools.EventClassification import EventClassifier
    from tools.KinematicsCalculator import KinematicsCalculator
    from tools.PlotLibrary import PLOT_SETTINGS
    from tools.FSIProbabilityBalance import diagnostic_cpp,generator_function,inspect_hadrons,BalanceAccumulator,KE_EDGES
    if AC.is_pc or SC.USE_SWAPPED or abs(int(SC.AnaNuPDG))!=14 or AC.skip_2p2h:
        raise RuntimeError('Use unswapped CCnumu/CCnumubar MC without skip_2p2h')
    if not (SU.USE_MNVTUNE_V4P3P0 or SU.USE_MNVTUNE_V4P3P1):raise RuntimeError('Use V430 or V431')
    out=Path(opts.balance_output).resolve();out.mkdir(parents=True,exist_ok=False)
    audit=dict(status='running',version=3,tree=opts.balance_tree,argv=remaining,
               playlist=AC.playlist,count=list(AC.count),AnaNuPDG=int(SC.AnaNuPDG),
               cohort='right-sign CC satisfying config.SignalDef.IsFiducial; all interaction types',
               baseline='unit and unchanged production V430 generator+flux without detector factors',
               source_path=opts.balance_source,
               selection_cuts=list(SAMPLE_CUTS['Signal']),kinematics_cuts=list(KINEMATICS_CUTS),
               SU_sha256=hashlib.sha256(Path(SU.__file__).read_bytes()).hexdigest(),
               KE_edges_MeV=list(KE_EDGES),POT_normalized=False,
               warning='alternative combinations and additive budgets are diagnostics, not proposed corrections')
    write_json(out/'audit.json',audit)
    try:
        sourcepath=Path(opts.balance_source).resolve()
        headerpath=Path(opts.balance_header).resolve() if opts.balance_header else sourcepath.with_suffix('.h')
        cpp=sourcepath.read_text();header=headerpath.read_text()
        audit.update(header_path=str(headerpath),cpp_sha256=hashlib.sha256(cpp.encode()).hexdigest(),
                     header_sha256=hashlib.sha256(header.encode()).hexdigest())
        chain=Utilities.fileChain(AC.playlist,'mc',AC.ntuple_tag,
                                  'Truth' if opts.balance_tree=='truth' else None,AC.count[0],AC.count[1])
        n=int(chain.GetEntries())
        if not n:raise RuntimeError('Input tree has no entries')
        raw=chain.GetChain();raw.LoadTree(0)
        # A missing generator record must never silently turn into a reco-only closure claim.
        for branch in ('mc_er_nPart','mc_er_ID','mc_er_status','mc_er_FD','mc_er_LD','mc_er_mother',
                       'mc_er_Px','mc_er_Py','mc_er_Pz','mc_er_E'):
            if not raw.GetBranch(branch):
                raise RuntimeError('{} tree lacks {}; need generator record to validate probability balance'.format(opts.balance_tree,branch))
        audit['input_files']=[str(f.GetTitle()) for f in raw.GetListOfFiles()]
        universes=SU.GetAllSystematicsUniverses(chain,False,False,exclude=['all'])
        if set(universes)!={'cv'} or len(universes['cv'])!=1:raise RuntimeError('Expected CV only')
        u=universes['cv'][0]
        derived=diagnostic_cpp(header,cpp,str(headerpath))
        (out/'calculator_clone.C').write_text(derived)
        if not ROOT.gInterpreter.Declare(derived):raise RuntimeError('Could not compile diagnostic calculator clone')
        calc=ROOT.PlotUtils.FSIBalanceCalculator();calc.UseTrackingThreshold()
        fate=ROOT.PlotUtils.fate_code()
        common_weight=generator_function(Path(SU.__file__).read_text(),vars(SU))
        is_reco=opts.balance_tree=='reco'
        if is_reco:
            kin=KinematicsCalculator(correct_beam_angle=True,correct_MC_energy_scale=False,calc_true=True,is_pc=False)
            classifier=EventClassifier(classifiers=['Reco','Truth'],use_kin_cuts=True,use_sideband=[])
            u.LoadTools(kin,classifier)
        en4=PLOT_SETTINGS['Biased Neutrino Energy'];edges=list(map(float,en4['binning'][0]))
        audit['EN4_edges_GeV']=edges
        if AC.testing:n=min(n,2000)
        if opts.balance_max_events is not None:n=min(n,opts.balance_max_events)
        audit['processed_entries']=n
        audit['tune_flags']={k:v for k,v in vars(SU).items() if k.startswith(('USE_MNVTUNE_','V4P3P0_'))
                             and isinstance(v,(str,bool,int,float,type(None)))}
        write_json(out/'audit.json',audit)
        scalars=('mc_incoming','mc_primaryLepton','mc_charm','mc_intType','mc_targetA','mc_targetZ','mc_resID','mc_er_nPart')
        ivectors=('mc_er_ID','mc_er_status','mc_er_FD','mc_er_LD','mc_er_mother')
        dvectors=('mc_er_Px','mc_er_Py','mc_er_Pz','mc_er_E')
        names={1:'QE',2:'Delta',3:'DIS',4:'Coherent',8:'2p2h'}
        accumulator=BalanceAccumulator();cohort=selected=0
        max_clone=max_combiner=max_factorization=0.
        event_fields=['chain_entry','category','targetA','targetZ','stages','true_Enu_GeV','generator_baseline','standard',
                      'pure_all','pure_retained','pure_retained_delta','without_0p70','without_1p03',
                      'extend_three_rule_to_ge4','single_hadron_balance','product_correlation','four_slot_cap',
                      'resonance_adjustment','event_combination_rules','ntot','preFSI_hadrons']
        hadron_fields=['chain_entry','category','targetA','stages','hadron_index','pdg','fate','ke_MeV',
                       'raw1','adjusted1','effective1','retained','eligible']
        with gzip.open(out/'events.csv.gz','wt',newline='') as ef,gzip.open(out/'hadrons.csv.gz','wt',newline='') as hf:
            ew=csv.DictWriter(ef,fieldnames=event_fields);ew.writeheader()
            hw=csv.DictWriter(hf,fieldnames=hadron_fields);hw.writeheader()
            for entry in range(n):
                if entry%10000==0:print('[FSI balance] {} entry {}/{} cohort {}'.format(opts.balance_tree,entry,n,cohort),flush=True)
                u.SetEntry(entry);u.ResetWeight()
                u.__dict__.pop('_cached_new_true_l_over_e',None)
                if u.mc_current!=1 or u.mc_incoming!=int(SC.AnaNuPDG) or not IsFiducial(u):continue
                cohort+=1
                integers=[int(u.GetInt(k)) for k in scalars]
                vectors=[u.GetVecInt(k) for k in ivectors]+[u.GetVecDouble(k) for k in dvectors]
                baseline=float(common_weight(u));standard=float(u.GetFSIWeight(1))
                calc.calcWeights(*(integers+vectors+[entry,False]))
                cloned=float(calc.GetElasticWeight(1))
                if not all(math.isfinite(x) for x in (baseline,standard,cloned)):raise RuntimeError('Nonfinite event weight at {}'.format(entry))
                max_clone=max(max_clone,abs(standard-cloned))
                if not math.isclose(standard,cloned,rel_tol=1e-10,abs_tol=1e-10):
                    raise RuntimeError('Installed library/source clone mismatch at {}'.format(entry))
                fate.calcFates(*(integers[:6]+integers[7:]+vectors+[entry,False]))
                hadrons=[]
                for h in fate.GetIntranukeParticles():
                    pdg=int(getattr(h,'__pdg'));f=int(getattr(h,'__fate'));ke=float(getattr(h,'__ke'))
                    weights=ROOT.PlotUtils.FSIBalanceRawWeights(calc,f,ke/1000.,integers[4],pdg)
                    hadrons.append(dict(pdg=pdg,fate=f,ke_MeV=ke,raw1=float(weights[0])))
                records,variants,budget,ntot=inspect_hadrons(hadrons,integers[3],integers[6],integers[2],integers[4])
                max_combiner=max(max_combiner,abs(standard-variants['standard']))
                if not math.isclose(standard,variants['standard'],rel_tol=1e-10,abs_tol=1e-10):
                    raise RuntimeError('Combination reconstruction mismatch at {}'.format(entry))
                cpp_ntot=int(getattr(calc,'__count_elastic_pion'))+int(getattr(calc,'__count_elastic_nucleon'))
                if cpp_ntot!=ntot:raise RuntimeError('Contributor count mismatch at {}'.format(entry))
                category=names.get(integers[3],'Other')+('_charm' if integers[2] else '')
                stages=['before']
                if is_reco:
                    kin.CalculateKinematics(u);classifier.Classify(u)
                    if classifier.side_band=='Signal':
                        selected+=1;stages.append('selected_all_EN4')
                        energy=float(en4['value_getter'][0](u))
                        if not math.isfinite(energy):raise RuntimeError('Nonfinite selected EN4')
                        if edges[0]<=energy<edges[-1]:stages.append('selected_regular_EN4')
                        # Check the new common baseline against the full production formula.
                        detector=1.
                        if SU.V4P3P0_ADD_GEANT_HADRON:detector*=u.GetGeantHadronWeight()
                        if SU.V4P3P0_ADD_MINOS_EFFICIENCY:detector*=u.GetMyMinosEfficiencyWeight()
                        full=float(u.GetStandardWeightV4p3p0())
                        max_factorization=max(max_factorization,abs(full-baseline*detector))
                        if not math.isclose(full,baseline*detector,rel_tol=1e-10,abs_tol=1e-10):raise RuntimeError('Baseline factorization failed')
                context=dict(chain_entry=entry,category=category,targetA=integers[4],targetZ=integers[5],
                             stages=';'.join(stages),true_Enu_GeV=float(u.mc_incomingE)*1e-3)
                ew.writerow(dict(context,generator_baseline=baseline,ntot=ntot,preFSI_hadrons=len(records),**variants,**budget))
                for i,h in enumerate(records):
                    hw.writerow(dict({k:context[k] for k in ('chain_entry','category','targetA','stages')},hadron_index=i,**h))
                accumulator.add(category,stages,baseline,variants,budget,records,ntot,integers[4],float(u.mc_incomingE)*1e-3)
        rows=accumulator.event_rows()
        write_csv(out/'event_balance.csv',rows)
        write_csv(out/'hadron_balance.csv',accumulator.hadron_rows())
        totals=[r for r in rows if r['split']=='all']
        write_json(out/'summary.json',dict(totals=totals))
        audit.update(status='complete',cohort_events=cohort,selected_events=selected,
                     max_source_clone_difference=max_clone,max_combination_difference=max_combiner,
                     max_baseline_factorization_difference=max_factorization,
                     normalization_validated=False)
        write_json(out/'audit.json',audit)
        for row in totals:
            if row['category']=='DIS':print('[DIS BALANCE]',row,flush=True)
        print('[FSI balance] Complete:',out,flush=True)
    except Exception as exc:
        audit.update(status='failed',error=str(exc));write_json(out/'audit.json',audit);raise


if __name__=='__main__':main()

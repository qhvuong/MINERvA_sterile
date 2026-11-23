#!/usr/bin/env python3
"""Read-only ME CV weight tracing using the actual selection and weight methods.

Place in studies/ and run in the normal GPVM MAT environment. No analysis files
are edited. v4.3.1 and its nested v4.3.0 weight are evaluated on the SAME events.
Only reconstructed Signal events are traced; this is not an efficiency study.
"""
import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import heapq
import json
import math
import os
from pathlib import Path
import sys
import numpy as np

METHODS = {
    'genie': 'GetGenieWeight', 'flux': 'GetFluxAndCVWeight',
    '2p2h': 'GetLowRecoil2p2hWeight', 'RPA': 'GetRPAWeight',
    'lowQ2pi': 'GetMyLowQ2PiWeight', 'MEcoherent': 'GetCOHPionWeight',
    'geant': 'GetGeantHadronWeight', 'MINOS': 'GetMyMinosEfficiencyWeight',
    'elasticFSI': 'GetFSIWeight',
}
COMPONENTS = ['genie','flux','2p2h','RPA','lowQ2pi','MEcoherent','diffractive','geant','MINOS','elasticFSI']


class Trace:
    def __init__(self,universe):
        self.u=universe
        self.reset()
        for component,method in METHODS.items():
            bound=getattr(universe,method)
            def wrapped(*args,_bound=bound,_component=component,**kwargs):
                value=float(_bound(*args,**kwargs))
                if not math.isfinite(value):
                    raise RuntimeError(f'Nonfinite {_component}: {value}')
                self.calls[_component].append(value)
                if _component=='elasticFSI':
                    self.fsi_arguments.append(list(args))
                return value
            setattr(universe,method,wrapped)
        bound430=universe.GetStandardWeightV4p3p0
        def nested430():
            value=float(bound430())
            self.base_weights.append(value)
            return value
        universe.GetStandardWeightV4p3p0=nested430

    def reset(self):
        self.calls=defaultdict(list)
        self.base_weights=[]
        self.fsi_arguments=[]

    def evaluate(self,SU):
        self.reset()
        self.u.ResetWeight()
        actual=float(self.u.GetWeight(False))
        if len(self.base_weights)!=1:
            raise RuntimeError('Expected one nested v4.3.0 evaluation; tracing failed or weight implementation changed')
        for key in ['genie','flux','RPA','elasticFSI']:
            if len(self.calls[key])!=1:
                raise RuntimeError(f'{key}: expected one captured call, got {self.calls[key]}')
        for key,values in self.calls.items():
            if len(values)>1:
                raise RuntimeError(f'Repeated component call: {key}: {values}; inspect before interpreting')
        if self.fsi_arguments != [[1]]:
            raise RuntimeError(f'Unexpected FSI arguments: {self.fsi_arguments}')
        w={key:(self.calls[key][0] if self.calls[key] else 1.) for key in COMPONENTS}
        # The only literal multiplicative correction in this v4.3.0 function.
        if SU.V4P3P0_ADD_DIFFRACTIVE and self.u.mc_intType==4 and self.u.mc_targetZ==1:
            w['diffractive']=1.4368
        base=self.base_weights[0]
        product430=math.prod(w[key] for key in COMPONENTS if key!='elasticFSI')
        product431=product430*w['elasticFSI']
        if not all(math.isfinite(x) for x in [actual,base,product430,product431]):
            raise RuntimeError('Nonfinite total CV weight')
        if not math.isclose(base,product430,rel_tol=1e-9,abs_tol=1e-12):
            raise RuntimeError(f'v4.3.0 product closure failed: {base} vs {product430}; calls={dict(self.calls)}')
        if not math.isclose(actual,product431,rel_tol=1e-9,abs_tol=1e-12):
            raise RuntimeError(f'v4.3.1 product closure failed: {actual} vs {product431}')
        return w,base,actual


def ratio(n,d):
    return float(n/d) if d!=0 else None


def write_csv(path,rows):
    if rows:
        with path.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def summarize(records,indices,label):
    rows=[]
    if not len(indices):
        return rows
    total=records['w431'][indices];base=records['w430'][indices]
    for j,name in enumerate(COMPONENTS):
        values=records['components'][indices,j]
        # Multiplication, not division: correct even for components equal to zero.
        without=np.prod(np.delete(records['components'][indices],j,axis=1),axis=1)
        rows.append(dict(group=label,component=name,n=len(indices),mean=float(values.mean()),
            median=float(np.median(values)),p95=float(np.quantile(values,.95)),p99=float(np.quantile(values,.99)),
            minimum=float(values.min()),maximum=float(values.max()),negative=int((values<0).sum()),
            zero=int((values==0).sum()),fraction_changed=float((np.abs(values-1)>1e-9).mean()),
            sum_w430=float(base.sum()),sum_w431=float(total.sum()),sum_without_component=float(without.sum()),
            full_over_without=ratio(total.sum(),without.sum()),
            share_of_full_yield_from_component_gt5=ratio(total[values>5].sum(),total.sum()),
            share_of_full_yield_from_component_gt10=ratio(total[values>10].sum(),total.sum())))
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path)
    p.add_argument('--playlist',required=True,choices=['1A','1B','1C','1D','1E','1F','1G','1L','1M','1N','1O','1P','5A','6A','6B','6C','6D','6E','6F','6G','6H','6I','6J'])
    p.add_argument('--selection',choices=['ryan','original'],default='ryan')
    p.add_argument('--count',type=int,nargs=2,default=[0,1],metavar=('FIRST_FILE','N_FILES'))
    p.add_argument('--max-events',type=int,default=-1,help='Maximum scanned reco entries, -1 means all in the file range')
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--no-plots',action='store_true',help='Save numerical results only; Matplotlib is not imported')
    a=p.parse_args()
    if a.count[0]<0 or a.count[1]<=0 or (a.max_events!=-1 and a.max_events<=0):
        p.error('Invalid count or max-events')
    root=(a.repo or Path(__file__).resolve().parents[1]).resolve()
    if not (root/'tools/SystematicsUniverse.py').is_file():
        p.error('Place in studies/ or provide --repo')
    out=a.output_dir.resolve()
    out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()):
        p.error('Use an empty output directory to avoid mixing diagnostic runs')
    cfg=root/'configs'/('CCNuMu' if a.playlist.startswith('1') else 'CCAntiNuMu')
    sys.path[:0]=[str(cfg),str(root)]
    os.environ['CCNUEROOT']=str(root);os.environ['CONFIGPATH']=str(cfg)
    os.environ['CCNuEConfig']=cfg.name
    os.chdir(root)
    sys.argv=[str(root/'selection/eventSelection.py'),'--playlist',f'me{a.playlist}_p8',
        '--ntuple_tag','MAD','--mc_only','--exclude-universes','all',
        '--selection_tag','ME_weight_diagnostic','--count',*map(str,a.count)]
    import ROOT
    ROOT.PyConfig.IgnoreCommandLineOptions=True
    ROOT.gROOT.SetBatch(True)
    from config import CutConfig,SystematicsConfig,PlotConfig
    if a.selection=='ryan':
        CutConfig.SAMPLE_CUTS['Signal']=[x for x in CutConfig.SAMPLE_CUTS['Signal'] if x not in ('DeadTime','Neutrino Helicity')]
        CutConfig.visE_RANGE=[0.,.7 if a.playlist.startswith('1') else .3]
    from config.AnalysisConfig import AnalysisConfig
    from tools import SystematicsUniverse as SU,Utilities
    from tools.KinematicsCalculator import KinematicsCalculator
    from tools.EventClassification import EventClassifier
    SU.USE_MNVTUNE_V1=False;SU.USE_MNVTUNE_V4P4P1=False
    SU.USE_MNVTUNE_V4P3P0=False;SU.USE_MNVTUNE_V4P3P1=True
    if SystematicsConfig.USE_NUE_CONSTRAINT or SystematicsConfig.USE_SWAPPED or SU.V4P3P0_ADD_LE_COH:
        raise RuntimeError('Diagnostic supports the unconstrained, non-swapped ME v4.3.x prescription only')
    info=dict(playlist=a.playlist,selection=a.selection,count=a.count,max_events=a.max_events,
        signal_cuts=CutConfig.SAMPLE_CUTS['Signal'],Eavail_range=CutConfig.visE_RANGE,
        Pt_range=CutConfig.RECO_PT_RANGE,PLOTUTILSROOT=os.environ.get('PLOTUTILSROOT'),
        tune='Trace actual v4.3.1 and nested v4.3.0; no systematic universes',status='STARTED',
        note='Component omission is a diagnostic variation, not a prescribed uncertainty or alternative tune.')
    source_paths=[root/'tools/SystematicsUniverse.py',root/'tools/KinematicsCalculator.py',root/'tools/EventClassification.py',cfg/'config/CutConfig.py']
    info['source_sha256']={str(path.relative_to(root)):hashlib.sha256(path.read_bytes()).hexdigest() for path in source_paths}
    info['flags']={k:v for k,v in vars(SU).items() if k.startswith('V4P3P0_') and isinstance(v,(str,bool,int,float,type(None)))}
    (out/'run.json').write_text(json.dumps(info,indent=2))
    print('[WEIGHT DIAGNOSTIC]',json.dumps(info,indent=2),flush=True)
    chain=Utilities.fileChain(AnalysisConfig.playlist,'mc',AnalysisConfig.ntuple_tag,None,*a.count)
    kin=KinematicsCalculator(correct_beam_angle=True,correct_MC_energy_scale=False,calc_true=True,is_pc=False)
    classifier=EventClassifier(classifiers=['Reco','Truth'],use_kin_cuts=True,use_sideband=[])
    universes=SU.GetAllSystematicsUniverses(chain,False,False,['all'])
    uv=[u for group in universes.values() for u in group]
    if len(uv)!=1:
        raise RuntimeError('Expected one CV universe')
    u=uv[0];u.LoadTools(kin,classifier)
    trace=Trace(u)
    edges=np.asarray(PlotConfig.NEUTRINO4_EE_BINNING,dtype=float)
    n_scan=chain.GetEntries()
    if a.max_events>0:
        n_scan=min(n_scan,a.max_events)
    values=[];w430=[];w431=[];energies=[];categories=[];targets=[];pdgs=[]
    top=[];max_closure=0.
    fields=['entry','category','pdg','intType','targetZ','true_Enu_GeV','EN4_GeV','Emu_GeV','Eavail_GeV',
            'w430','w431','two_p_two_h_evaluated','tuple_mc_cvweight_total',*COMPONENTS]
    stream=gzip.open(out/'selected_events.csv.gz','wt',newline='')
    writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
    try:
        for entry in range(n_scan):
            if entry%10000==0:
                print('[SCAN]',entry,'/',n_scan,'selected',len(values),flush=True)
            u.SetEntry(entry);u.ResetWeight()
            kin.CalculateKinematics(u);classifier.Classify(u)
            if classifier.side_band!='Signal':
                continue
            en4=float(kin.reco_E_lep+kin.reco_visE)
            if not math.isfinite(en4):
                raise RuntimeError(f'Nonfinite selected EN4 at entry {entry}')
            weights,base,full=trace.evaluate(SU)
            max_closure=max(max_closure,abs(full-math.prod(weights.values())))
            category=str(classifier.truth_class or 'Other')
            row=dict(entry=entry,category=category,pdg=int(u.mc_incoming),intType=int(u.mc_intType),
                targetZ=int(u.mc_targetZ),true_Enu_GeV=float(u.mc_incomingE)*1e-3,EN4_GeV=en4,
                Emu_GeV=float(kin.reco_E_lep),Eavail_GeV=float(kin.reco_visE),w430=base,w431=full,
                two_p_two_h_evaluated=bool(trace.calls['2p2h']),
                tuple_mc_cvweight_total=float(u.mc_cvweight_total),**weights)
            writer.writerow(row)
            values.append([weights[k] for k in COMPONENTS]);w430.append(base);w431.append(full)
            energies.append(en4);categories.append(category);targets.append(row['targetZ']);pdgs.append(row['pdg'])
            score=abs(full-base)
            item=(score,entry,row)
            if len(top)<50:heapq.heappush(top,item)
            elif score>top[0][0]:heapq.heapreplace(top,item)
    except BaseException:
        info.update(status='FAILED',selected_before_failure=len(values),failed_entry=entry)
        (out/'run.json').write_text(json.dumps(info,indent=2))
        raise
    finally:
        stream.close()
    if not values:
        raise RuntimeError('No selected events')
    records=dict(components=np.asarray(values),w430=np.asarray(w430),w431=np.asarray(w431))
    energy=np.asarray(energies);cats=np.asarray(categories);tz=np.asarray(targets);pdg=np.asarray(pdgs)
    groups={'All':np.arange(len(values))}
    groups.update({cat:np.flatnonzero(cats==cat) for cat in sorted(set(categories))})
    for cat in sorted(set(categories)):
        for z in sorted(set(tz[cats==cat])):
            groups[f'{cat}:targetZ={z}']=np.flatnonzero((cats==cat)&(tz==z))
    rows=[r for label,indices in groups.items() for r in summarize(records,indices,label)]
    write_csv(out/'component_summary.csv',rows)
    write_csv(out/'largest_FSI_yield_changes.csv',[x[2] for x in sorted(top,reverse=True)])
    binrows=[];saved={'edges':edges}
    for label,indices in groups.items():
        if ':targetZ' in label:continue
        mask=np.zeros(len(values),dtype=bool);mask[indices]=True
        old=np.histogram(energy[mask],edges,weights=records['w430'][mask])[0]
        new=np.histogram(energy[mask],edges,weights=records['w431'][mask])[0]
        saved[label+'_w430']=old;saved[label+'_w431']=new
        for j,name in enumerate(COMPONENTS):
            without=np.prod(np.delete(records['components'][mask],j,axis=1),axis=1)
            omitted=np.histogram(energy[mask],edges,weights=without)[0]
            saved[label+'_without_'+name]=omitted
            for b,(lo,hi) in enumerate(zip(edges[:-1],edges[1:])):
                binrows.append(dict(category=label,component=name,low=lo,high=hi,w430=old[b],w431=new[b],
                    without_component=omitted[b],full_over_without=ratio(new[b],omitted[b])))
    write_csv(out/'component_EN4.csv',binrows)
    np.savez_compressed(out/'spectra.npz',**saved)
    info.update(status='COMPLETED',entries_scanned=n_scan,selected=len(values),max_absolute_product_residual=max_closure,
        sum_w430=float(records['w430'].sum()),sum_w431=float(records['w431'].sum()),
        weighted_FSI_ratio=ratio(records['w431'].sum(),records['w430'].sum()),
        outside_EN4_range=int(((energy<edges[0])|(energy>=edges[-1])).sum()),
        negative_total_weights=int((records['w431']<0).sum()))
    (out/'run.json').write_text(json.dumps(info,indent=2,allow_nan=False))
    print('[RESULT]',json.dumps(info,indent=2),flush=True)
    if a.no_plots:
        print('Numerical outputs saved; plotting skipped:',out,flush=True)
        return
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,5,figsize=(17,7),constrained_layout=True)
    shown=[c for c in ['CCNuMuQE','CCNuMu2p2h','CCNuMuDelta','CCNuMuDIS','CCNuMuWrongSign'] if c in groups]
    for j,(ax,name) in enumerate(zip(axes.flat,COMPONENTS)):
        for cat in shown:
            x=records['components'][groups[cat],j]
            # Values outside the displayed 0..5 range remain in CSV summaries.
            ax.hist(x,bins=np.linspace(0,5,101),histtype='step',weights=np.ones(len(x))/len(x),label=cat.replace('CCNuMu',''))
        ax.axvline(1,color='gray',ls=':');ax.set(title=name,xlabel='Component weight',ylabel='Fraction / 0.05',yscale='log',xlim=(0,5))
    axes.flat[0].legend(fontsize=7)
    fig.suptitle('ME '+a.playlist+' — actual selected-event CV components; tails above 5 listed in CSV')
    fig.savefig(out/'component_distributions.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,5),constrained_layout=True)
    for label in ['All',*shown]:
        n=saved[label+'_w431'];d=saved[label+'_w430']
        ax.stairs(np.divide(n,d,out=np.full_like(n,np.nan),where=d!=0),edges,label=label.replace('CCNuMu',''))
    ax.axhline(1,color='gray',ls='--');ax.set(xlabel='Neutrino Energy Estimator (GeV)',ylabel='v4.3.1 / v4.3.0',title='ME '+a.playlist+' — same events, FSI effect',xlim=(.8,20))
    ax.legend();ax.grid(alpha=.2);fig.savefig(out/'FSI_effect_by_category.png',dpi=180);plt.close(fig)
    print('Outputs:',out,flush=True)


if __name__=='__main__':
    main()

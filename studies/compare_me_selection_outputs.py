#!/usr/bin/env python3
"""Compare the existing FullCuts/v4.3.1 and RyansCuts/v4.3.0,v4.3.1 ME outputs.

Read-only; no analysis-config imports. Run in the usual GPVM PyROOT environment.
Uses raw EN4 counts and each file's Meta/POT_Used. Excludes flow bins from yields.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np

CASES = {
    'full431': ('FullCuts', 'v4p3p1', 'Original cuts, v4.3.1'),
    'ryan431': ('RyansCuts', 'v4p3p1', 'Ryan cuts, v4.3.1'),
    'ryan430': ('RyansCuts', 'v4p3p0', 'Ryan cuts, v4.3.0'),
}
CATEGORIES = ['NCDiff','NuEElastic','NC','NonPhaseSpace','CCNuMuWrongSign',
              'CCNuMuQE','CCNuMuDelta','CCNuMuDIS','CCNuMu2p2h','CCNuMu','CC','Other']


def divide(n,d):
    return np.divide(n,d,out=np.full_like(np.asarray(n,dtype=float),np.nan),where=np.asarray(d)!=0)


def hist_arrays(h):
    if h.GetDimension() != 1:
        raise RuntimeError('Expected 1D histogram: '+h.GetName())
    n=h.GetNbinsX()
    e=np.array([h.GetBinLowEdge(i) for i in range(1,n+2)],dtype=float)
    y=np.array([h.GetBinContent(i) for i in range(1,n+1)],dtype=float)
    err=np.array([h.GetBinError(i) for i in range(1,n+1)],dtype=float)
    if not np.isfinite(y).all() or not np.isfinite(err).all() or not np.all(np.diff(e)>0):
        raise RuntimeError('Invalid histogram: '+h.GetName())
    return e,y,err


def read_file(ROOT,path,hname):
    f=ROOT.TFile.Open(str(path),'READ')
    if not f or f.IsZombie():
        raise RuntimeError('Cannot open '+str(path))
    try:
        if f.TestBit(ROOT.TFile.kRecovered):
            raise RuntimeError('Recovered/incomplete file: '+str(path))
        meta=f.Get('Meta')
        if not meta or not meta.GetBranch('POT_Used'):
            raise RuntimeError('Missing Meta/POT_Used: '+str(path))
        pots=[]
        for i in range(meta.GetEntries()):
            if meta.GetEntry(i)<=0:
                raise RuntimeError('Cannot read Meta entry')
            v=float(meta.GetLeaf('POT_Used').GetValue())
            if not math.isfinite(v) or v<0:
                raise RuntimeError('Invalid POT_Used')
            pots.append(v)
        pot=math.fsum(pots)
        if pot<=0:
            raise RuntimeError('Nonpositive total POT')
        h=f.Get(hname)
        if not h:
            raise RuntimeError(f'Missing {hname}; keys: {[k.GetName() for k in f.GetListOfKeys()]}')
        edges,y,err=hist_arrays(h)
        cats={}
        for category in CATEGORIES:
            ch=f.Get(hname+'_'+category)
            if ch:
                ce,cy,_=hist_arrays(ch)
                if not np.array_equal(ce,edges):
                    raise RuntimeError('Category binning differs: '+category)
                cats[category]=cy
        return dict(path=str(path),pot=pot,meta_entries=len(pots),edges=edges,y=y,err=err,
                    categories=cats,keys=[k.GetName() for k in f.GetListOfKeys()],
                    underflow=float(h.GetBinContent(0)),overflow=float(h.GetBinContent(h.GetNbinsX()+1)))
    finally:
        f.Close()


def safe_number(x):
    return float(x) if math.isfinite(x) else None


def write_csv(path,rows):
    if rows:
        with path.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base',type=Path,default=Path('/exp/minerva/data/users/qvuong'))
    p.add_argument('--output-dir',type=Path,default=Path('/exp/minerva/data/users/qvuong/ME_selection_output_comparison'))
    p.add_argument('--hist',default='EN4')
    a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
    import ROOT
    ROOT.PyConfig.IgnoreCommandLineOptions=True
    ROOT.gROOT.SetBatch(True)
    import PlotUtils  # noqa: F401
    data={};issues=[];metadata={};summary=[];binrows=[];categoryrows=[];effectrows=[]
    for beam,folder,channel in [('FHC','nu_mu','CCnumu'),('RHC','antinu_mu','CCnumubar')]:
        data[beam]={}
        for case,(cuts,tune,label) in CASES.items():
            tag=f'ME_{channel}_{cuts}_sameOffset_{tune}'
            pair={kind:read_file(ROOT,a.base/folder/f'kin_dist_{kind}me{beam}_{tag}_MAD.root',a.hist)
                  for kind in ('data','mc')}
            data[beam][case]=pair
            d,m=pair['data'],pair['mc']
            if not np.array_equal(d['edges'],m['edges']):
                raise RuntimeError(f'{beam} {case}: data/MC binning differs')
            scale=d['pot']/m['pot'];pred=m['y']*scale
            for kind,r in pair.items():
                metadata[f'{beam}_{case}_{kind}']={k:v for k,v in r.items() if k not in ('y','edges','err','categories')}
            category_total=sum(m['categories'].values(),np.zeros_like(m['y']))
            closure=bool(np.allclose(category_total,m['y'],rtol=1e-8,atol=1e-7))
            if not closure:
                issues.append(f'{beam} {case}: category sum differs from total; inspect keys before interpreting fractions')
            if 'CCNuMuQE' not in m['categories']:
                issues.append(f'{beam} {case}: missing QE histogram')
            e=d['edges']
            windows=[('all',np.ones_like(d['y'],dtype=bool)),('2to8',(e[:-1]>=2-1e-9)&(e[1:]<=8+1e-9))]
            for window,mask in windows:
                ds=float(d['y'][mask].sum());ms=float(pred[mask].sum());rawmc=float(m['y'][mask].sum())
                qe=m['categories'].get('CCNuMuQE',np.zeros_like(m['y']))
                summary.append(dict(beam=beam,case=case,window=window,data_POT=d['pot'],mc_POT=m['pot'],
                    data=ds,mc_at_data_POT=ms,data_over_mc=safe_number(ds/ms) if ms else None,
                    QE_fraction=safe_number(float(qe[mask].sum())/rawmc) if rawmc and 'CCNuMuQE' in m['categories'] else None,
                    category_closure=closure))
            for i,(lo,hi) in enumerate(zip(e[:-1],e[1:])):
                binrows.append(dict(beam=beam,case=case,low=lo,high=hi,data=d['y'][i],mc_at_data_POT=pred[i],
                    data_over_mc=safe_number(divide(d['y'],pred)[i]),
                    data_stat_error_over_mc=safe_number(divide(d['err'],pred)[i])))
            for category,y in m['categories'].items():
                for i,(lo,hi) in enumerate(zip(e[:-1],e[1:])):
                    categoryrows.append(dict(beam=beam,case=case,category=category,low=lo,high=hi,
                        mc_at_data_POT=y[i]*scale,fraction=safe_number(divide(y,m['y'])[i])))
            print(f'{beam} {case}: POT data={d["pot"]:.9g}, MC={m["pot"]:.9g}, '
                  f'D/MC={d["y"].sum()/pred.sum():.6f}, '
                  f'QE fraction={summary[-2]["QE_fraction"]}, category closure={closure}',flush=True)
        ref=data[beam]['full431']['data']['edges']
        if any(not np.array_equal(r[k]['edges'],ref) for r in data[beam].values() for k in ('data','mc')):
            raise RuntimeError(beam+': case binning differs')
        for kind in ('data','mc'):
            first=data[beam]['full431'][kind]
            for case in ('ryan431','ryan430'):
                current=data[beam][case][kind]
                if first['meta_entries']!=current['meta_entries'] or not math.isclose(first['pot'],current['pot'],rel_tol=1e-10):
                    issues.append(f'{beam} {kind}: {case} exposure differs from FullCuts; compare input coverage before interpreting selection changes')
        d0=data[beam]['ryan430']['data'];d1=data[beam]['ryan431']['data']
        same=(d0['meta_entries']==d1['meta_entries'] and math.isclose(d0['pot'],d1['pot'],rel_tol=1e-12)
              and np.array_equal(d0['y'],d1['y']) and np.array_equal(d0['err'],d1['err'])
              and d0['underflow']==d1['underflow'] and d0['overflow']==d1['overflow'])
        metadata[beam+'_Ryan_data_identical']=same
        if not same:
            issues.append(beam+': Ryan-cut data differ between tune-tagged runs; tune should not alter data')
        print(f'{beam}: Ryan-cut data identical between tune runs: {same}',flush=True)
        # Ratios here are per-POT yields; they do not mix changing data exposures.
        for effect,new,old in [('all_cut_changes','ryan431','full431'),('elastic_FSI_431_over_430','ryan431','ryan430')]:
            for kind in ('data','mc'):
                n=data[beam][new][kind];d=data[beam][old][kind]
                ratios=divide(n['y']/n['pot'],d['y']/d['pot'])
                for i in range(len(ratios)):
                    effectrows.append(dict(beam=beam,effect=effect,kind=kind,low=ref[i],high=ref[i+1],
                        new_over_old_per_POT=safe_number(ratios[i])))
    write_csv(a.output_dir/'summary.csv',summary)
    write_csv(a.output_dir/'bin_ratios.csv',binrows)
    write_csv(a.output_dir/'category_fractions.csv',categoryrows)
    write_csv(a.output_dir/'changes.csv',effectrows)
    report={'issues':issues,'files':metadata,
        'notes':['Raw CV counts; spectra divided by bin width only for display.',
                 'Ratio error bars show data statistical uncertainty only; no significance or chi-square calculated.',
                 'Variants reuse events and are correlated; ratio curves have no independent-sample errors.',
                 'FullCuts to RyansCuts combines DeadTime, helicity and FHC Eavail changes.',
                 'QE means the stored right-sign generator QE category, divided by all selected MC; wrong-sign QE is not separately resolved.',
                 'No legacy-flux-tagged outputs were supplied.']}
    (a.output_dir/'comparison.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors={'full431':'tab:blue','ryan431':'tab:orange','ryan430':'tab:green'}
    fig,axes=plt.subplots(2,2,figsize=(13,8),sharex='col',gridspec_kw={'height_ratios':[3,1.5]},constrained_layout=True)
    for j,beam in enumerate(('FHC','RHC')):
        ax,rat=axes[:,j]
        for case,(_,_,label) in CASES.items():
            d,m=data[beam][case]['data'],data[beam][case]['mc'];e=d['edges'];width=np.diff(e);x=(e[:-1]+e[1:])/2
            pred=m['y']*d['pot']/m['pot']
            ax.stairs(pred/width,e,color=colors[case],label=label+' MC')
            rat.errorbar(x,divide(d['y'],pred),yerr=divide(d['err'],pred),fmt='o',ms=3,color=colors[case],label=label)
            if case in ('full431','ryan431') or (case=='ryan430' and not metadata[beam+'_Ryan_data_identical']):
                ax.errorbar(x,d['y']/width,yerr=d['err']/width,fmt='o' if case=='full431' else 's',
                    ms=3,color=colors[case],mfc='none',label=label+' data')
        ax.set(title='ME '+beam,ylabel='Events / GeV');ax.legend(fontsize=8)
        rat.axhline(1,color='gray',ls='--');rat.set(xlabel='Neutrino Energy Estimator (GeV)',ylabel='Data / MC',xlim=(.8,20))
        rat.set_ylim(.6,1.4);rat.grid(alpha=.2)
    fig.suptitle('MINERvA Work In Progress\nMC normalized to each case’s data POT',fontsize=13)
    fig.savefig(a.output_dir/'spectra_and_data_mc.png',dpi=180);fig.savefig(a.output_dir/'spectra_and_data_mc.pdf');plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(12,7),sharex=True,constrained_layout=True)
    for j,beam in enumerate(('FHC','RHC')):
        for kind,style in [('data','--'),('mc','-')]:
            for i,(new,old) in enumerate([('ryan431','full431'),('ryan431','ryan430')]):
                n=data[beam][new][kind];d=data[beam][old][kind]
                axes[i,j].stairs(divide(n['y']/n['pot'],d['y']/d['pot']),n['edges'],ls=style,label=kind.upper())
        for i in range(2):
            axes[i,j].axhline(1,color='gray',ls=':');axes[i,j].grid(alpha=.2);axes[i,j].legend()
            axes[i,j].set(title='ME '+beam,xlim=(.8,20),ylabel='Ryan / original per POT' if i==0 else 'v4.3.1 / v4.3.0 per POT')
        axes[1,j].set_xlabel('Neutrino Energy Estimator (GeV)')
    fig.savefig(a.output_dir/'cut_and_FSI_changes.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,4),sharey=True,constrained_layout=True)
    for ax,beam in zip(axes,('FHC','RHC')):
        for case,(_,_,label) in CASES.items():
            m=data[beam][case]['mc']
            if 'CCNuMuQE' in m['categories']:
                ax.stairs(divide(m['categories']['CCNuMuQE'],m['y']),m['edges'],label=label,color=colors[case])
        ax.set(title='ME '+beam,xlabel='Neutrino Energy Estimator (GeV)',ylabel='Selected right-sign QE / total MC',xlim=(.8,20),ylim=(0,1))
        ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.savefig(a.output_dir/'QE_fractions.png',dpi=180);plt.close(fig)
    print('\nReports and plots:',a.output_dir)
    for issue in issues:
        print('[CHECK]',issue)


if __name__=='__main__':
    main()

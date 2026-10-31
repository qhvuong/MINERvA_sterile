#!/usr/bin/env python3
"""Read-only inspection of p8 LE5 truth-derived fluxes and audited event CV flux weights.
Requires PyROOT and matplotlib; does not re-run compute_flux or modify input ROOT files.
"""
import argparse
import csv
import math
import os
import ROOT
ROOT.gROOT.SetBatch(True)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DEFAULT_DIR = '/exp/minerva/data/users/qvuong/flux_studies/producedFluxes_p8_LE_newElectronResponse'
SPECIES = {'numu': 14, 'numubar': -14}

def histograms(f, prefix=''):
    result={}
    for key in f.GetListOfKeys():
        name=key.GetName()
        obj=key.ReadObj()
        path=prefix+name
        if obj.InheritsFrom('TDirectory'):
            result.update(histograms(obj,path+'/'))
        elif obj.InheritsFrom('TH1') and obj.GetDimension()==1:
            result[path]=obj
    return result

def extract(h):
    n=h.GetNbinsX()
    edges=np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1,n+2)])
    vals=np.array([h.GetBinContent(i) for i in range(1,n+1)])
    errs=np.array([h.GetBinError(i) for i in range(1,n+1)])
    return edges,vals,errs

def find_hist(hist, target):
    exact=[k for k in hist if k.split('/')[-1] == target]
    if len(exact)!=1:
        raise RuntimeError(f'Expected one {target}; found {exact}. Available: {list(hist)[:100]}')
    return hist[exact[0]]

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--dir',default=DEFAULT_DIR)
    p.add_argument('--audit',default='/exp/minerva/data/users/qvuong/flux_studies/rhc_cv_weight_audit.csv')
    p.add_argument('--outdir',default='/exp/minerva/data/users/qvuong/flux_studies/le5_flux_cv_checks')
    p.add_argument('--max-energy',type=float,default=20.)
    a=p.parse_args()
    os.makedirs(a.outdir,exist_ok=True)
    audit=[]
    if os.path.exists(a.audit):
        with open(a.audit,newline='') as f:
            audit=list(csv.DictReader(f))
        print(f'Audit rows loaded: {len(audit)} from {a.audit}')
    else:
        print(f'Audit CSV not found: {a.audit}; examining ROOT fluxes only')
    for species,pdg in SPECIES.items():
        fn=os.path.join(a.dir,f'LE5_{species}.root')
        f=ROOT.TFile.Open(fn,'READ')
        if not f or f.IsZombie():
            print(f'ERROR opening {fn}'); continue
        hs=histograms(f)
        print(f'\n=== {fn} ===')
        print('Histogram names:')
        for k in sorted(hs): print(' ',k)
        try:
            hu=find_hist(hs,'flux_E_unweighted')
            hc=find_hist(hs,'flux_E_cvweighted')
        except RuntimeError as e:
            print('Cannot compare flux:',e); f.Close(); continue
        edges,u,eu=extract(hu); ec,c,ec_err=extract(hc)
        if not np.array_equal(edges,ec):
            raise RuntimeError(f'Incompatible bin edges in {fn}')
        finite=np.isfinite(u)&np.isfinite(c)&(u>0)
        ratio=np.full_like(u,np.nan)
        ratio[finite]=c[finite]/u[finite]
        print('Bin [GeV]         unweighted flux         CV flux         CV/unweighted')
        for i in range(len(u)):
            if edges[i]>=a.max_energy: break
            print(f'{edges[i]:5.1f}-{edges[i+1]:5.1f}   {u[i]:18.8g}  {c[i]:18.8g}  {ratio[i]:13.6g}')
        if np.any(~np.isfinite(u)) or np.any(~np.isfinite(c)):
            print('WARNING: nonfinite flux bins')
        centers=(edges[:-1]+edges[1:])/2
        fig,ax=plt.subplots(figsize=(9,5))
        ax.step(edges[:-1],ratio,where='post',label='Produced flux CV/unweighted')
        ax.axhline(1,color='black',ls='--',lw=1)
        ax.set_xlim(0,a.max_energy)
        ax.set_xlabel(r'True $E_\nu$ (GeV)')
        ax.set_ylabel('CV-weighted / unweighted produced flux')
        ax.set_title(f'LE5 {species} (PDG {pdg})')
        ax.grid(alpha=.2)
        ax.legend()
        fig.tight_layout()
        out=os.path.join(a.outdir,f'LE5_{species}_produced_flux_ratio.png')
        fig.savefig(out,dpi=170)
        plt.close(fig)
        rows=[r for r in audit if int(r['pdg'])==pdg]
        if rows:
            print(f'\nAudited selected reconstructed events, PDG {pdg}: {len(rows)}')
            print('Etrue bin [GeV]    n audit    mean event GetFluxAndCVWeight     produced flux ratio')
            e=np.array([float(r['Etrue_GeV']) for r in rows])
            w=np.array([float(r['flux']) for r in rows])
            for i in range(len(u)):
                if edges[i]>=a.max_energy:break
                selection=(e>=edges[i])&(e<edges[i+1])&np.isfinite(w)
                if np.sum(selection):
                    print(f'{edges[i]:5.1f}-{edges[i+1]:5.1f}       {np.sum(selection):5d}       {np.mean(w[selection]):12.5f}            {ratio[i]:12.5f}')
            print('NOTE: event flux weight and produced CV/unweighted flux ratio are')
            print('NOT necessarily identical: the latter includes MC CV event weighting')
            print('used to reconstruct the flux; the former is a downstream reweighter.')
        f.Close()
        print('Saved',out)
    print('\nDone. No ROOT files modified.')

if __name__=='__main__': main()

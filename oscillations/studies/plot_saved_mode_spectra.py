#!/usr/bin/env python3
"""Compare saved fixed-K null/best spectra without ROOT or another fit.

python plot_saved_mode_spectra.py --six K6_DIR --full FULL_DIR --out PLOT_DIR
Both directories contain CASE_profiled_spectra.csv and optionally
sterile_refits.csv. Full mode counts are read from the files, not assumed.
Missing full-basis cases are explicitly marked, never filled with another K.
"""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

CASES = ('direct_all', 'ratio_all', 'ratios_dropped')
LABELS = ('Direct CC electron', 'All-sample ratios', 'Ratios excluded from profiling')
BLUE, ORANGE = '#167bb7', '#d97819'


def load(path):
    with path.open() as f:
        rows = list(csv.DictReader(f))
    out = {}
    for point in ('null', 'best'):
        selected = sorted((r for r in rows if r['point']==point),key=lambda r:int(r['bin']))
        if len(selected)!=54 or [int(r['bin']) for r in selected]!=list(range(1,55)):
            raise ValueError(f'{path}: expected all 54 scored bins at null and best')
        out[point] = {key:np.array([float(r[key]) for r in selected])
                      for key in ('data','mc','flux_shift','profiled_mc')}
        np.testing.assert_allclose(out[point]['mc']+out[point]['flux_shift'],
                                   out[point]['profiled_mc'],rtol=1e-10,atol=1e-8)
    np.testing.assert_allclose(out['null']['data'],out['best']['data'])
    ks={int(r['k']) for r in rows}
    if len(ks)!=1:raise ValueError(f'{path}: expected one K')
    out['k']=ks.pop()
    return out


def decorate(ax, case, labels=False):
    boundaries=(6.5,18.5,30.5,42.5)
    for x in boundaries:ax.axvline(x,color='#b8bec5',lw=.8,ls=':')
    if labels:
        titles=['Elastic','FHC CC e','RHC CC e','FHC CC mu','RHC CC mu']
        if case!='direct_all':titles[1:3]=['FHC mu/e','RHC mu/e']
        for x,title in zip((3.5,12.5,24.5,36.5,48.5),titles):
            ax.text(x,1.025,title,transform=ax.get_xaxis_transform(),ha='center',fontsize=9)
    ax.set_xlim(.5,54.5);ax.grid(axis='y',alpha=.16)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--six',type=Path,required=True)
    p.add_argument('--full',type=Path,required=True)
    p.add_argument('--out',type=Path,default=Path('mode_spectrum_comparison'))
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    scans=[]
    for case in CASES:
        six=load(args.six/f'{case}_profiled_spectra.csv')
        if six['k']!=6:raise ValueError('The --six directory must contain K=6 spectra')
        path=args.full/f'{case}_profiled_spectra.csv'
        full=load(path) if path.exists() else None
        if full:
            expected=30 if case=='ratios_dropped' else 54
            if full['k']!=expected:raise ValueError(f'{case}: full basis must have K={expected}')
            for key in ('data','mc'):
                np.testing.assert_allclose(six['null'][key],full['null'][key],rtol=1e-10,atol=1e-8)
        scans.append((six,full))
    fig,axs=plt.subplots(3,2,figsize=(16,10),sharex=True,layout='constrained')
    for case,label,(six,full),row in zip(CASES,LABELS,scans,axs):
        denom=six['null']['mc'];x=np.arange(1,55)
        for point,ax in zip(('null','best'),row):
            ax.axhline(1,color='#969da5',lw=1)
            ax.plot(x,six[point]['data']/denom,'o',color='black',ms=3,label='Data (no error bars in saved CSV)')
            ax.step(x,six[point]['profiled_mc']/denom,where='mid',color=BLUE,lw=1.8,label='K=6')
            if full:
                ax.step(x,full[point]['profiled_mc']/denom,where='mid',color=ORANGE,lw=1.8,label=f"Full basis, K={full['k']}")
            else:
                ax.text(.5,.06,'Full-basis K=30 fit not yet available',transform=ax.transAxes,
                        color='#a14f15',ha='center',fontsize=10,bbox=dict(fc='white',ec='none',alpha=.85))
            ax.set_title(label+(' — null' if point=='null' else ' — sterile best fit'),fontsize=12,pad=25)
            ax.set_ylabel('Ratio to raw null')
            decorate(ax,case,labels=True)
            if case=='direct_all' and point=='best':ax.legend(fontsize=9,loc='lower right')
    for ax in axs[-1]:ax.set_xlabel('Stitched bin')
    fig.suptitle('Normalized prior: K=6 versus the full available basis\n'
                 'Each hypothesis has its own fitted flux correction; current updated-statistics / EL inputs',fontsize=16)
    fig.savefig(args.out/'all_cases_mode_comparison.png',dpi=180)
    plt.close(fig)
    for case,label,(six,full) in zip(CASES,LABELS,scans):
        if full is None:continue
        fig,axs=plt.subplots(2,2,figsize=(16,7.5),sharex=True,
                             gridspec_kw={'height_ratios':[2.5,1]},layout='constrained')
        x=np.arange(1,55);denom=six['null']['mc']
        for j,point in enumerate(('null','best')):
            ax=axs[0,j]
            ax.step(x,denom,where='mid',color='#969da5',ls=':',label='Raw null')
            ax.plot(x,six[point]['data'],'o',color='black',ms=3,label='Data')
            ax.step(x,six[point]['profiled_mc'],where='mid',color=BLUE,lw=1.8,label='K=6')
            ax.step(x,full[point]['profiled_mc'],where='mid',color=ORANGE,lw=1.8,label=f"K={full['k']}")
            if all(np.all(a>0) for a in (six[point]['data'],denom,six[point]['profiled_mc'],full[point]['profiled_mc'])):ax.set_yscale('log')
            ax.set_ylabel('Stitched-bin content');ax.set_title('Profiled null' if point=='null' else 'Profiled sterile best fit',pad=25)
            decorate(ax,case,True)
            if j==1:ax.legend(fontsize=9,loc='best')
            ax=axs[1,j];ax.axhline(1,color='#969da5',lw=1)
            ax.plot(x,six[point]['data']/denom,'o',color='black',ms=3)
            ax.step(x,six[point]['profiled_mc']/denom,where='mid',color=BLUE,lw=1.6)
            ax.step(x,full[point]['profiled_mc']/denom,where='mid',color=ORANGE,lw=1.6)
            decorate(ax,case);ax.set(ylabel='Ratio to raw null',xlabel='Stitched bin')
        fig.suptitle(label+f" — K=6 versus K={full['k']}\nSaved central spectra; statistical errors were not included in these CSVs",fontsize=15)
        fig.savefig(args.out/f'{case}_mode_comparison.png',dpi=180)
        plt.close(fig)
    print(f'Wrote plots to {args.out}; missing full cases are marked explicitly.')


if __name__=='__main__':main()

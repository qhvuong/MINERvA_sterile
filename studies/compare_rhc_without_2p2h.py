#!/usr/bin/env python3
"""Compare stored RHC MC with and without its CCNuMu2p2h truth category.

This subtracts the 2p2h truth category ALREADY PRESENT in each ROOT MC file.
It does not subtract separately merged 2p2h productions and does not disable
2p2h CV weights on non-MEC events. No new event selection is needed.
"""
import argparse
from pathlib import Path
import numpy as np
import ROOT
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT.gROOT.SetBatch(True)

BASE = Path('/exp/minerva/data/users/qvuong/antinu_mu')
TAGS = {
    'No helicity': 'ValidateCuts_CCnumubar_Inclusive_noLeptonAngle_noHelicity_v1',
    'With helicity': 'ValidateCuts_CCnumubar_Inclusive_noLeptonAngle_v1',
}
SLICES = [(0,5),(5,10),(10,20),(20,40)]

def load(path, name, dim):
    f=ROOT.TFile.Open(str(path))
    if not f or f.IsZombie(): raise RuntimeError(f'Could not open {path}')
    h=f.Get(name)
    if not h or h.GetDimension()!=dim:
        f.Close(); raise RuntimeError(f'No {dim}D {name} in {path}')
    x=np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1,h.GetNbinsX()+2)],float)
    if dim==1:
        a=np.array([h.GetBinContent(i) for i in range(1,h.GetNbinsX()+1)],float)
    else:
        y=np.array([h.GetYaxis().GetBinLowEdge(j) for j in range(1,h.GetNbinsY()+2)],float)
        a=np.array([[h.GetBinContent(i,j) for j in range(1,h.GetNbinsY()+1)] for i in range(1,h.GetNbinsX()+1)],float)
    f.Close()
    if not np.isfinite(a).all(): raise RuntimeError(f'Nonfinite bins: {path}, {name}')
    return (x,a) if dim==1 else (x,y,a)

def ratio(d,mc,scale):
    den=mc*scale
    return np.divide(d,den,out=np.full(d.shape,np.nan,float),where=den>0)

def lines_plot(x,d,nominal,without,scale,out,title):
    center=(x[:-1]+x[1:])/2
    half=(x[1:]-x[:-1])/2
    rnom=ratio(d,nominal,scale); rno=ratio(d,without,scale)
    fig,ax=plt.subplots(figsize=(9,5.5))
    ax.errorbar(center,rnom,xerr=half,fmt='o-',label='Nominal MC',markersize=4)
    ax.errorbar(center,rno,xerr=half,fmt='s--',label='MC minus CCNuMu2p2h',markersize=4)
    ax.axhline(1,c='k',lw=1,ls=':')
    ax.set(xlabel='Reconstructed muon energy (GeV)',ylabel='Data / scaled MC',title=title)
    ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out,dpi=180);plt.close(fig)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--mc-scale',type=float,default=0.108,help='MC-to-data multiplier; verify POT provenance')
    p.add_argument('--outdir',default='/exp/minerva/data/users/qvuong/studies/rhc_without_2p2h')
    a=p.parse_args()
    out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
    for label,tag in TAGS.items():
        path=lambda typ:BASE/f'kin_dist_{typ}le5_{tag}_MAD.root'
        x,d=load(path('data'),'Eel',1)
        xm,m=load(path('mc'),'Eel',1)
        xp,mec=load(path('mc'),'Eel_CCNuMu2p2h',1)
        if not (np.array_equal(x,xm) and np.array_equal(x,xp)): raise RuntimeError('1D binning differs')
        without=m-mec
        print('\n===',label,'===')
        print(f'Eel: Data {d.sum():.1f}; MC {m.sum():.3f}; MEC {mec.sum():.3f} ({mec.sum()/m.sum()*100:.2f}% of MC); non-MEC {without.sum():.3f}')
        print(f'D/MC nominal: {d.sum()/(a.mc_scale*m.sum()):.4f}; minus MEC: {d.sum()/(a.mc_scale*without.sum()):.4f}')
        print('E_low E_high   MEC/MC%  D/MC nominal  D/MC minus-MEC')
        rn=ratio(d,m,a.mc_scale);ro=ratio(d,without,a.mc_scale)
        for i in range(len(d)):
            frac=100*mec[i]/m[i] if m[i]>0 else float('nan')
            print(f'{x[i]:5.2f} {x[i+1]:5.2f} {frac:9.2f} {rn[i]:13.4f} {ro[i]:17.4f}')
        slug='no_helicity' if label=='No helicity' else 'with_helicity'
        lines_plot(x,d,m,without,a.mc_scale,out/f'energy_nominal_vs_no_2p2h_{slug}.png',f'RHC {label}: effect of removing stored 2p2h category')
        xx,yy,dd=load(path('data'),'Elepton_LepAngle',2)
        xm,ym,mm=load(path('mc'),'Elepton_LepAngle',2)
        xp,yp,mecc=load(path('mc'),'Elepton_LepAngle_CCNuMu2p2h',2)
        if not all((np.array_equal(xx,q) for q in (xm,xp))) or not all((np.array_equal(yy,q) for q in (ym,yp))): raise RuntimeError('2D binning differs')
        if np.any(mm-mecc < -1e-6): print('NOTE: some negative non-MEC weighted bins; inspect MC weight signs')
        print('2D visible: MC MEC fraction %.2f%%'% (100*mecc.sum()/mm.sum()))
        fig,axs=plt.subplots(2,2,figsize=(12,8),sharex=True,sharey=True)
        for ax,(lo,hi) in zip(axs.flat,SLICES):
            sel=np.where((yy[:-1]>=lo-1e-9)&(yy[1:]<=hi+1e-9))[0]
            ds=dd[:,sel].sum(axis=1); ms=mm[:,sel].sum(axis=1); mec_s=mecc[:,sel].sum(axis=1)
            w=ms-mec_s
            r1=ratio(ds,ms,a.mc_scale);r2=ratio(ds,w,a.mc_scale)
            c=(xx[:-1]+xx[1:])/2
            ax.plot(c,r1,'o-',markersize=3,label='Nominal')
            ax.plot(c,r2,'s--',markersize=3,label='Minus 2p2h')
            ax.axhline(1,c='k',ls=':',lw=1);ax.grid(alpha=.2)
            ax.set_title(f'{lo}–{hi}°: MEC fraction={100*mec_s.sum()/ms.sum():.1f}%',fontsize=10)
            print(f'  Angle {lo:2g}–{hi:2g} deg: D/MC nominal={ds.sum()/(a.mc_scale*ms.sum()):.4f}, minus MEC={ds.sum()/(a.mc_scale*w.sum()):.4f}, MEC fraction={100*mec_s.sum()/ms.sum():.2f}%')
        axs.flat[0].legend(fontsize=8)
        fig.supxlabel('Reconstructed muon energy (GeV)');fig.supylabel('Data / scaled MC')
        fig.suptitle(f'RHC {label}: MC 2p2h subtraction by muon angle')
        fig.tight_layout();fig.savefig(out/f'energy_by_angle_nominal_vs_no_2p2h_{slug}.png',dpi=180);plt.close(fig)
    print('\nSaved plots:',out)
    print('CAUTION: This is category subtraction for diagnostics only; assess full MC composition and separate 2p2h production before drawing physics conclusions.')

if __name__=='__main__': main()

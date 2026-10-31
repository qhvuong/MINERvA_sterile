#!/usr/bin/env python3
"""RHC MINOS-matched 2D E_mu versus theta diagnostic from existing ROOT files.

Example:
 python studies/compare_rhc_energy_angle_2d.py --mc-scale 0.108

Scale is an MC-to-data exposure multiplier; confirm its provenance before physical interpretation.
"""
import argparse
from pathlib import Path
import numpy as np
import ROOT
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

ROOT.gROOT.SetBatch(True)
BASE = Path('/exp/minerva/data/users/qvuong/antinu_mu')
TAGS = {
    'No helicity': 'ValidateCuts_CCnumubar_Inclusive_noLeptonAngle_noHelicity_v1',
    'With helicity': 'ValidateCuts_CCnumubar_Inclusive_noLeptonAngle_v1',
}
GROUPS = [(0,5), (5,10), (10,20), (20,40)]

def read_hist(path, name):
    f=ROOT.TFile.Open(str(path))
    if not f or f.IsZombie():
        raise RuntimeError(f'Cannot open: {path}')
    h=f.Get(name)
    if not h or h.GetDimension()!=2:
        raise RuntimeError(f'Missing 2D histogram {name} in {path}')
    nx,ny=h.GetNbinsX(),h.GetNbinsY()
    x=np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1,nx+2)], dtype=float)
    y=np.array([h.GetYaxis().GetBinLowEdge(j) for j in range(1,ny+2)], dtype=float)
    z=np.array([[h.GetBinContent(i,j) for j in range(1,ny+1)] for i in range(1,nx+1)], dtype=float)
    e2=np.array([[h.GetBinError(i,j)**2 for j in range(1,ny+1)] for i in range(1,nx+1)], dtype=float)
    f.Close()
    if not np.isfinite(z).all() or not np.isfinite(e2).all():
        raise RuntimeError(f'Nonfinite bins in {path}')
    return x,y,z,e2

def ratio_plot(x,y,data,mc,mc_scale,out,title,min_mc):
    expected=mc*mc_scale
    ratio=np.divide(data,expected,out=np.full_like(data,np.nan),where=expected>0)
    ratio=np.where((expected>=min_mc)&(data>=3),ratio,np.nan)
    fig,ax=plt.subplots(figsize=(10,6.5))
    im=ax.pcolormesh(x,y,ratio.T,shading='flat',cmap='RdBu_r',norm=TwoSlopeNorm(vmin=0.6,vcenter=1,vmax=1.4))
    cb=fig.colorbar(im,ax=ax,label='Data / scaled MC')
    ax.set_xlabel('Reconstructed muon energy (GeV)')
    ax.set_ylabel('Reconstructed muon angle (deg)')
    ax.set_title(title)
    fig.tight_layout();fig.savefig(out,dpi=190);plt.close(fig)

def projection_by_angle(x,y,d,m,vd,vm,mc_scale,out,title):
    fig,axs=plt.subplots(2,2,figsize=(12,8),sharex=True,sharey=True)
    for ax,(lo,hi) in zip(axs.flat,GROUPS):
        selected=np.where((y[:-1]>=lo-1e-8)&(y[1:]<=hi+1e-8))[0]
        dproj=d[:,selected].sum(axis=1); mproj=m[:,selected].sum(axis=1)*mc_scale
        vdproj=vd[:,selected].sum(axis=1); vmproj=vm[:,selected].sum(axis=1)*(mc_scale**2)
        valid=(mproj>0)&(dproj>=3)
        r=np.divide(dproj,mproj,out=np.full_like(dproj,np.nan),where=mproj>0)
        er=np.full_like(r,np.nan)
        er[valid]=np.sqrt(vdproj[valid]/mproj[valid]**2 + dproj[valid]**2*vmproj[valid]/mproj[valid]**4)
        centers=(x[:-1]+x[1:])/2; halfwidth=(x[1:]-x[:-1])/2
        ax.errorbar(centers[valid],r[valid],xerr=halfwidth[valid],yerr=er[valid],fmt='o',markersize=3,capsize=2)
        ax.axhline(1,c='black',ls='--',lw=1)
        ax.set_title(f'{lo}°–{hi}° (Data={dproj.sum():.0f}, scaled MC={mproj.sum():.0f})',fontsize=10)
        ax.grid(alpha=.2)
        print(f'  theta {lo:2g}-{hi:2g} deg: data={dproj.sum():.0f}, scaled MC={mproj.sum():.1f}, D/MC={dproj.sum()/mproj.sum() if mproj.sum()>0 else float("nan"):.4f}')
    fig.supxlabel('Reconstructed muon energy (GeV)');fig.supylabel('Data / scaled MC');fig.suptitle(title)
    fig.tight_layout();fig.savefig(out,dpi=190);plt.close(fig)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--mc-scale',type=float,default=None,help='MC-to-data scaling (e.g. POT_data / POT_MC); must verify')
    p.add_argument('--min-scaled-mc',type=float,default=10,help='Minimum scaled MC per 2D cell shown')
    p.add_argument('--outdir',default='/exp/minerva/data/users/qvuong/studies/rhc_energy_angle_2d')
    a=p.parse_args()
    out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
    pairs={}
    for label,tag in TAGS.items():
        p_data=BASE/f'kin_dist_datale5_{tag}_MAD.root'
        p_mc=BASE/f'kin_dist_mcle5_{tag}_MAD.root'
        x,y,d,vd=read_hist(p_data,'Elepton_LepAngle')
        xm,ym,m,vm=read_hist(p_mc,'Elepton_LepAngle')
        if not (np.array_equal(x,xm) and np.array_equal(y,ym)):
            raise RuntimeError(f'Data/MC binning mismatch: {label}')
        pairs[label]=(x,y,d,vd,m,vm)
        print(f'\n{label}: 2D visible data={d.sum():.0f}, MC stored={m.sum():.3f}; angle range={y[0]:g}–{y[-1]:g} deg')
        if a.mc_scale is not None:
            print(f'  scaled MC={m.sum()*a.mc_scale:.2f}; D/MC={d.sum()/(m.sum()*a.mc_scale):.5f}')
            suffix='no_helicity' if label=='No helicity' else 'with_helicity'
            ratio_plot(x,y,d,m,a.mc_scale,out/f'ratio_heatmap_{suffix}.png',f'RHC {label}: 2D data/MC',a.min_scaled_mc)
            print('Angular slice projection:')
            projection_by_angle(x,y,d,m,vd,vm,a.mc_scale,out/f'energy_ratios_by_angle_{suffix}.png',f'RHC {label}: D/MC vs Eμ by angle')
    x,y,d0,vd0,m0,vm0=pairs['No helicity']
    x1,y1,d1,vd1,m1,vm1=pairs['With helicity']
    if not (np.array_equal(x,x1) and np.array_equal(y,y1)):
        raise RuntimeError('Different 2D binnings before vs after helicity')
    if np.any(d1>d0+1e-6):print('WARNING: data after > before in at least one cell; verify selections')
    # Exposure-independent helicity survival double ratio in each cell.
    eff_d=np.divide(d1,d0,out=np.full_like(d0,np.nan),where=d0>0)
    eff_m=np.divide(m1,m0,out=np.full_like(m0,np.nan),where=m0>0)
    double=np.divide(eff_d,eff_m,out=np.full_like(eff_d,np.nan),where=eff_m>0)
    double=np.where((d0>=10)&(m0>=50),double,np.nan)
    fig,ax=plt.subplots(figsize=(10,6.5))
    im=ax.pcolormesh(x,y,double.T,cmap='RdBu_r',norm=TwoSlopeNorm(vmin=.6,vcenter=1,vmax=1.4),shading='flat')
    fig.colorbar(im,ax=ax,label='(Data after / before) / (MC after / before)')
    ax.set(xlabel='Reconstructed muon energy (GeV)',ylabel='Reconstructed muon angle (deg)',title='RHC helicity double survival ratio (POT-independent)')
    fig.tight_layout();fig.savefig(out/'helicity_double_ratio_2d.png',dpi=190);plt.close(fig)
    print(f'\nSaved diagnostics in {out.resolve()}')
    if a.mc_scale is None:print('No MC scale supplied: produced only exposure-independent helicity 2D double ratio.')

if __name__=='__main__':main()

#!/usr/bin/env python3
"""Compare RHC data/MC before and after helicity, with MINOS match in both.

Usage:
  python studies/compare_rhc_data_mc_helicity.py
  python studies/compare_rhc_data_mc_helicity.py --mc-scale 0.108

--mc-scale is the SINGLE external factor multiplying stored MC to data exposure.
Do not use the example value unless verified from production POT and normalization.
Without --mc-scale, only shape and normalization-independent quantities are plotted.
"""
import argparse
import os
import numpy as np
import ROOT
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT.gROOT.SetBatch(True)
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--directory', default='/exp/minerva/data/users/qvuong/antinu_mu')
p.add_argument('--outdir', default='studies/rhc_data_mc_helicity')
p.add_argument('--mc-scale', type=float, default=None,
               help='Verified multiplier for stored MC to match data POT/exposure; NOT inferred from D/MC')
a = p.parse_args()
os.makedirs(a.outdir, exist_ok=True)
if a.mc_scale is not None and (not np.isfinite(a.mc_scale) or a.mc_scale <= 0):
    p.error('--mc-scale must be finite and positive')

TAGS = {
    'Before helicity': 'ValidateCuts_CCnumubar_FV_InverseHasNoBackExitingTracks_DeadTime_LeptonAngle_v1',
    'After helicity': 'ValidateCuts_CCnumubar_Inclusive_v1',
}
RIGHT = ['CCNuMuQE', 'CCNuMuDelta', 'CCNuMuDIS', 'CCNuMu2p2h', 'CCNuMu']
WRONG = ['CCNuMuWrongSign']
OTHER = ['CC', 'NC', 'NuEElastic', 'NCDiff', 'Other']

def hist(path, names):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError('Cannot open ' + path)
    out = {}
    try:
        for name in names:
            h = f.Get(name)
            if not h:
                raise RuntimeError('Missing histogram %s in %s' % (name,path))
            edges = np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1,h.GetNbinsX()+2)])
            vals = np.array([h.GetBinContent(i) for i in range(1,h.GetNbinsX()+1)])
            errs = np.array([h.GetBinError(i) for i in range(1,h.GetNbinsX()+1)])
            if not np.isfinite(vals).all() or not np.isfinite(errs).all():
                raise RuntimeError('Nonfinite values/errors in %s of %s' % (name,path))
            out[name] = (edges,vals,errs)
    finally:
        f.Close()
    return out

result={}
for label, tag in TAGS.items():
    datafile = os.path.join(a.directory, 'kin_dist_datale5_'+tag+'_MAD.root')
    mcfile = os.path.join(a.directory, 'kin_dist_mcle5_'+tag+'_MAD.root')
    d = hist(datafile,['Eel'])['Eel']
    mc = hist(mcfile,['Eel']+['Eel_'+c for c in RIGHT+WRONG+OTHER])
    edges, dv, de = d
    _, mv, me = mc['Eel']
    for key,(e,_,_) in mc.items():
        if not np.allclose(edges,e,atol=1e-9,rtol=0):
            raise RuntimeError('Binning mismatch in '+key)
    comps = {k:sum(mc['Eel_'+c][1] for c in cats) for k,cats in
             [('right',RIGHT),('wrong',WRONG),('other',OTHER)]}
    closure = mv - sum(comps.values())
    if not np.allclose(mv, sum(comps.values()), atol=1e-6, rtol=1e-9):
        raise RuntimeError('MC truth categories do not close for '+label)
    result[label] = dict(edges=edges, data=dv, dataerr=de, mc=mv, mcerr=me, **comps)
    print('\n%s\nData = %.3f, MC stored = %.3f, raw D/MC = %.6f' %
          (label,dv.sum(),mv.sum(),dv.sum()/mv.sum()))
    print('Truth components: right %.2f, wrong %.2f, other %.2f' %
          tuple(comps[k].sum() for k in ('right','wrong','other')))

before,after = result.values()
if not np.allclose(before['edges'],after['edges'],atol=1e-9,rtol=0):
    raise RuntimeError('Selections have different Eel binning')
if np.any(after['data'] > before['data'] + 1e-8):
    raise RuntimeError('Data selection not nested bin by bin')
if np.any(after['mc'] > before['mc'] + 1e-8):
    print('WARNING: MC after larger than before in at least one bin. Check nested productions.')
edges = before['edges']; centers=(edges[:-1]+edges[1:])/2
xerr=np.diff(edges)/2

def divide(x,y):
    return np.divide(x,y,out=np.full_like(x,np.nan,dtype=float),where=y>0)

r_before=divide(before['data'],before['mc'])
r_after=divide(after['data'],after['mc'])
double=divide(r_after,r_before)
# Use existing sumw2 errors for the two data/MC ratios; no covariance assumed here.
# Double-ratio error intentionally omitted: nested samples are correlated.

def ratioerr(res):
    d,m,de,me=(res[k] for k in ['data','mc','dataerr','mcerr'])
    err=np.full_like(d,np.nan,dtype=float)
    valid=(d>0)&(m>0)
    err[valid] = (d[valid]/m[valid])*np.sqrt((de[valid]/d[valid])**2+(me[valid]/m[valid])**2)
    return err

print('\n E_low E_high    D/MC before  D/MC after  after/before   right% after  wrong% after')
for i in range(len(centers)):
    print('%6.2f %6.2f  %12.5f %11.5f %12.5f %12.2f %12.2f' %
          (edges[i],edges[i+1],r_before[i],r_after[i],double[i],
           100*divide(after['right'],after['mc'])[i],
           100*divide(after['wrong'],after['mc'])[i]))
print('\nIntegrated helicity double survival ratio:',
      (after['data'].sum()/before['data'].sum())/(after['mc'].sum()/before['mc'].sum()))

scale = 1.0 if a.mc_scale is None else a.mc_scale
fig,ax=plt.subplots(figsize=(10,5.5))
for label,res,marker in [('Before helicity',before,'o'),('After helicity',after,'s')]:
    r=divide(res['data'],res['mc'])*scale**-1
    ax.errorbar(centers,r,yerr=ratioerr(res)/scale,xerr=xerr,fmt=marker,capsize=2,label=label)
ax.set_xlabel('Reconstructed muon energy (GeV)')
ax.set_ylabel('Data / MC'+(' (exposure normalized)' if a.mc_scale is not None else ' (unscaled MC)'))
ax.set_title('RHC: MINOS-matched data/MC before vs after helicity')
ax.grid(alpha=.25);ax.legend();fig.tight_layout()
fig.savefig(os.path.join(a.outdir,'data_mc_before_after.png'),dpi=180);plt.close(fig)

fig,ax=plt.subplots(figsize=(10,5.5))
ax.plot(centers,double,'o-',label='(D/MC) after / before')
ax.axhline(1,ls='--',color='black',alpha=.7)
ax.set_xlabel('Reconstructed muon energy (GeV)')
ax.set_ylabel('Helicity double survival ratio')
ax.set_title('RHC: helicity effect (independent of global POT scale)')
ax.grid(alpha=.25);ax.legend();fig.tight_layout()
fig.savefig(os.path.join(a.outdir,'helicity_double_ratio.png'),dpi=180);plt.close(fig)

# Shape-only plot to compare energy dependence without invoking POT scaling.
fig,ax=plt.subplots(figsize=(10,5.5))
for label,res,marker in [('Before helicity',before,'o'),('After helicity',after,'s')]:
    dr=divide(res['data'],res['data'].sum()*np.ones_like(res['data']))
    mr=divide(res['mc'],res['mc'].sum()*np.ones_like(res['mc']))
    ax.plot(centers,divide(dr,mr),marker+'-',label=label)
ax.axhline(1,ls='--',color='black',alpha=.7)
ax.set_xlabel('Reconstructed muon energy (GeV)')
ax.set_ylabel('(Data fraction in bin) / (MC fraction in bin)')
ax.set_title('RHC: shape-only comparison')
ax.grid(alpha=.25);ax.legend();fig.tight_layout()
fig.savefig(os.path.join(a.outdir,'shape_only_data_mc.png'),dpi=180);plt.close(fig)

print('\nPlots saved to',a.outdir)
if a.mc_scale is None:
    print('IMPORTANT: --mc-scale not set: absolute data/MC ratios are NOT POT normalized.')
    print('Supply verified MC-to-data exposure multiplier to interpret absolute deficit.')

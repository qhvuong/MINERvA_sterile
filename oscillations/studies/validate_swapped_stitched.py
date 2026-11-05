#!/usr/bin/env python3
"""Read-only saved-component and production-oscillation checks; no fit/profiling.
Run from oscillations/ with the usual analysis environment.
"""
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np


def vec(h):
    return np.array([h.GetBinContent(i) for i in range(1, h.GetNbinsX()+1)])


def averages(h, dm2):
    # Stitched templates: X = reconstructed bin, Y = true L/E.
    edges = np.array([h.GetYaxis().GetBinLowEdge(j) for j in range(1,h.GetNbinsY()+2)])
    widths = np.diff(edges)
    if not np.all(widths > 0):
        raise ValueError('Invalid L/E edges')
    mid = (edges[1:]+edges[:-1])/2
    k = 1.27*dm2
    sin2 = .5*(1-np.cos(2*k*mid)*np.sinc(k*widths/np.pi))
    contents = np.array([[h.GetBinContent(i,j) for j in range(1,h.GetNbinsY()+1)]
                         for i in range(1,h.GetNbinsX()+1)])
    norms = contents.sum(axis=1)
    return np.divide(contents@sin2,norms,out=np.zeros_like(norms),where=norms!=0), contents


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tag',default='predNueel_noRatio_p8_v4p3p0')
    ap.add_argument('--output',type=Path,default=Path('swapped_stitched_validation'))
    ap.add_argument('--expect-fhc-swap',type=float)
    ap.add_argument('--expect-rhc-swap',type=float)
    args=ap.parse_args()
    osc=Path(__file__).resolve().parents[1]
    sys.path.insert(0,str(osc.parent))
    # Project config parsers must not consume diagnostic arguments.
    sys.argv=[sys.argv[0]]
    import ROOT
    import PlotUtils
    from tools.StitchedHistogram import StitchedHistogram, minBinCont
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    args.output.mkdir(parents=True,exist_ok=True)
    cfg=json.loads((osc/('HIST_CONFIG_'+args.tag+'.json')).read_text())
    path=osc/'rootfiles'/('NuE_stitched_hists_'+args.tag+'.root')
    f=ROOT.TFile.Open(str(path),'READ')
    if not f or f.IsZombie():raise RuntimeError('Cannot open '+str(path))
    def get(name):
        obj=f.Get(name)
        if not obj:raise RuntimeError('Missing '+name)
        return obj
    # Read only the objects needed by the production OscillateHistogram method.
    # Keep f open; avoid Load(), covariance construction and TFile.Close().
    h=StitchedHistogram('sample')
    for attr,key in [('mc_hist','mc'),('nue_hist','nue'),('numu_hist','numu'),
                     ('swap_hist','swap'),('elastic_id','nueelID'),('ratio_id','ratioID'),
                     ('nue_template','nue_template'),('numu_template','numu_template'),
                     ('swap_template','swap_template')]:
        setattr(h,attr,get('sample_'+key))
    E,M,S=map(vec,[h.nue_hist,h.numu_hist,h.swap_hist])
    rid=vec(h.ratio_id)>.5
    failures=[]
    def check(label,actual,expected):
        a,b=np.asarray(actual),np.asarray(expected)
        ok=a.shape==b.shape and np.all(np.isfinite(a)) and np.allclose(a,b,rtol=1e-8,atol=1e-7)
        print(('PASS ' if ok else 'FAIL ')+label)
        if not ok:failures.append(label)
    coverage=[]
    swap_norms=averages(h.swap_template,0)[1].sum(axis=1)
    for sample,limits in cfg.items():
        ix=np.arange(limits['start'],limits['end']+1)
        local=get('swap_'+sample)
        if sample in ('fhc_nue_selection','rhc_nue_selection'):
            expected_integral=(args.expect_fhc_swap if sample.startswith('fhc')
                               else args.expect_rhc_swap)
            if expected_integral is not None:
                check(sample+' saved full swap integral vs stitch log',
                      local.Integral(),expected_integral)
        mc=get('mc_'+sample)
        selected=[i for i in range(1,mc.GetNbinsX()+1) if mc.GetBinContent(i)>minBinCont]
        check(sample+' local -> stitched swap mapping',S[ix],
              [local.GetBinContent(i) for i in selected])
        band=local.GetVertErrorBand('Flux')
        global_band=h.swap_hist.GetVertErrorBand('Flux')
        check(sample+' Flux universe count',band.GetNHists(),1000)
        check(sample+' global Flux universe count',global_band.GetNHists(),1000)
        if band.GetNHists()==global_band.GetNHists():
            a=np.array([[global_band.GetHist(int(u)).GetBinContent(int(i)+1) for i in ix]
                        for u in range(band.GetNHists())])
            b=np.array([[band.GetHist(int(u)).GetBinContent(int(i)) for i in selected]
                        for u in range(band.GetNHists())])
            check(sample+' all saved swap Flux-universe bins',a,b)
        if ('nue_selection' in sample or 'ratio' in sample):
            ratio=S[ix].sum()/E[ix].sum()
            print('  {}: intrinsic={:.8g} swap={:.8g} swap/intrinsic={:.5g}'.format(
                sample,E[ix].sum(),S[ix].sum(),ratio))
            # Informational ratio: do not require an arbitrary expected magnitude.
            for i in ix:
                if S[i]>0 and swap_norms[i]<=0:failures.append('empty swap L/E bin '+str(i+1))
                coverage.append(dict(sample=sample,global_bin=i+1,swap=S[i],
                                     template_in_range=swap_norms[i],
                                     fraction_missing=1-swap_norms[i]/S[i] if S[i] else 0))
    rows=[]
    points={'null':(0,0,0,0),'mu_only':(2,0,.1,0),'e_only':(2,.03,0,0),
            'mixed_2':(2,.03,.1,0),'mixed_39':(39.2,.002,.12,0)}
    for label,p in points.items():
        dm,ue,um,ut=p
        se,_=averages(h.nue_template,dm)
        sm,_=averages(h.numu_template,dm)
        ss,_=averages(h.swap_template,dm)
        appearance=S*4*ue*um*ss
        electron=E*(1-4*ue*(1-ue)*se)+appearance
        muon=M*(1-4*um*(1-um)*sm)
        expected=np.where(rid,np.divide(muon,electron,out=np.zeros_like(muon),where=electron!=0),muon+electron)
        h.OscillateHistogram(*p)
        actual=vec(h.GetOscillatedHistogram())
        check(label+' production vs independent oscillation',actual,expected)
        if label=='null':check('null closure to stored MC',actual,vec(h.mc_hist))
        saved=h.swap_hist
        try:
            h.swap_hist=saved.Clone();h.swap_hist.Reset()
            h.OscillateHistogram(*p)
            without=vec(h.GetOscillatedHistogram())
        finally:h.swap_hist=saved
        if ue==0 or um==0:check(label+' swap-off invariance',actual,without)
        elif not np.any(np.abs(actual-without)>1e-7):failures.append(label+' appearance has no effect')
        for sample,lim in cfg.items():
            ix=np.arange(lim['start'],lim['end']+1)
            if 'nue_selection' in sample or 'ratio' in sample:
                print('  {} {}: intrinsic loss={:.7g}, appearance gain={:.7g}, net electron change={:.7g}'.format(
                    label,sample,np.sum(E[ix]*4*ue*(1-ue)*se[ix]),appearance[ix].sum(),(electron[ix]-E[ix]).sum()))
            for i in ix:
                rows.append(dict(point=label,sample=sample,global_bin=i+1,nominal_e=E[i],swap=S[i],
                                 avg_sin2_swap=ss[i],appearance_gain=appearance[i],
                                 electron_prediction=electron[i],prediction=actual[i],
                                 prediction_without_swap=without[i]))
    for name,data in [('oscillation_bins.csv',rows),('template_coverage.csv',coverage)]:
        if data:
            with open(args.output/name,'w',newline='') as stream:
                w=csv.DictWriter(stream,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
    print('\nSUMMARY: {} failures'.format(len(failures)))
    for failure in failures:print('FAIL:',failure)
    print('Template coverage differences are reported, not automatically treated as failures.')
    print('Universe mapping checks verify saved bookkeeping, not event-level parent-flux weighting.')
    return int(bool(failures))

if __name__=='__main__':sys.exit(main())

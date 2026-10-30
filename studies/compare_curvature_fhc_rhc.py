#!/usr/bin/env python3
"""Compare helicity-only vs helicity+curvature LE CC muon selections.

Run in MINERvA PyROOT environment. Input paths default to the Oct 2026 productions.
Nothing is written to the input ROOT files.
"""
import argparse
import csv
import math
import os
import ROOT

ROOT.gROOT.SetBatch(True)
TAG = {
    'FHC': ('/exp/minerva/data/users/qvuong/nu_mu', 'leFHC', 'ValidateCuts_CCnumu_Inclusive_v1', 'ValidateCuts_CCnumu_Inclusive_Curvature_v1'),
    'RHC': ('/exp/minerva/data/users/qvuong/antinu_mu', 'le5', 'ValidateCuts_CCnumubar_Inclusive_v1', 'ValidateCuts_CCnumubar_Inclusive_Curvature_v1'),
}

def paths(mode):
    directory, playlist, old, new = TAG[mode]
    return { (cut, typ): os.path.join(directory, f'kin_dist_{typ}{playlist}_{tag}_MAD.root')
             for cut, tag in [('H',old), ('HC',new)] for typ in ('data','mc')}

def open_file(path):
    f = ROOT.TFile.Open(path, 'READ')
    if not f or f.IsZombie():
        raise RuntimeError(f'Cannot open ROOT file: {path}')
    return f

def pot(f):
    t=f.Get('Meta')
    if not t or not t.InheritsFrom('TTree'):
        raise RuntimeError(f'Missing Meta tree: {f.GetName()}')
    # Match production Utilities.getPOTFromFile, rather than assuming one entry.
    ch=ROOT.TChain('Meta')
    ch.Add(f.GetName())
    return float(ROOT.PlotUtils.POTCounter().getPOTfromTChain(ch))

def hist(f,name):
    obj=f.Get(name)
    if not obj or not obj.InheritsFrom('TH1'):
        raise RuntimeError(f'Missing TH1 {name!r} in {f.GetName()}')
    if obj.GetDimension()!=1:
        raise RuntimeError(f'{name!r} is {obj.GetDimension()}D; this script expects 1D')
    return obj

def edges_match(a,b):
    if a.GetNbinsX()!=b.GetNbinsX(): return False
    xa,xb=a.GetXaxis(),b.GetXaxis()
    return all(math.isclose(xa.GetBinLowEdge(i),xb.GetBinLowEdge(i),rel_tol=1e-10,abs_tol=1e-10)
               for i in range(1,a.GetNbinsX()+2))

def ratio(a,b):
    return a/b if b!=0 else math.nan

def fmt(v, digits=4):
    return f'{v:.{digits}f}' if math.isfinite(v) else 'n/a'

def do_mode(mode,key,rows):
    p=paths(mode)
    files={k:open_file(v) for k,v in p.items()}
    try:
        h={k:hist(f,key) for k,f in files.items()}
        baseline=h['H','data']
        for key2,v in h.items():
            if not edges_match(baseline,v): raise RuntimeError(f'Incompatible bins for {mode} {key2}')
        po={k:pot(f) for k,f in files.items()}
        for k,v in po.items():
            if not math.isfinite(v) or v<=0: raise RuntimeError(f'Nonpositive POT {mode}/{k}: {v}')
        print(f'\n===== {mode} / {key} =====')
        for (cut,typ),path in p.items():
            print(f'  {cut:2} {typ:4} POT={po[cut,typ]:.9g}  {path}')
        print('  H = helicity only; HC = helicity + curvature')
        print('  Yields in both selections scaled to the SAME baseline (H) data POT')
        print('  Excludes underflow/overflow. No bin-width scaling.')
        print('  Caution: survival fractions are weighted for MC and raw for data.')
        ref_pot=po['H','data']
        scales={k:ref_pot/v for k,v in po.items()}
        sumy={k:0.0 for k in h}
        print(f'{"E range":>15} {"D_H":>9} {"D_HC":>9} {"MC_H":>10} {"MC_HC":>10} {"D/MC H":>9} {"D/MC HC":>9} {"epsD":>8} {"epsMC":>8} {"epsD/epsMC":>11}')
        for i in range(1,baseline.GetNbinsX()+1):
            y={k:float(v.GetBinContent(i))*scales[k] for k,v in h.items()}
            for k,v in y.items(): sumy[k]+=v
            dh,dc,mh,mc=(y['H','data'],y['HC','data'],y['H','mc'],y['HC','mc'])
            ed=ratio(dc,dh); em=ratio(mc,mh)
            fields=dict(mode=mode,variable=key,lo=baseline.GetXaxis().GetBinLowEdge(i),hi=baseline.GetXaxis().GetBinUpEdge(i),
                        data_h=dh,data_hc=dc,mc_h=mh,mc_hc=mc,
                        dmc_h=ratio(dh,mh),dmc_hc=ratio(dc,mc),eps_data=ed,eps_mc=em,
                        double_ratio=ratio(ed,em),delta_dmc=ratio(dc,mc)-ratio(dh,mh))
            rows.append(fields)
            print(f'{fields["lo"]:6g}-{fields["hi"]:<8g} {dh:9.1f} {dc:9.1f} {mh:10.1f} {mc:10.1f} {fmt(fields["dmc_h"]):>9} {fmt(fields["dmc_hc"]):>9} {fmt(ed):>8} {fmt(em):>8} {fmt(fields["double_ratio"]):>11}')
        dh,dc,mh,mc=(sumy['H','data'],sumy['HC','data'],sumy['H','mc'],sumy['HC','mc'])
        ed=ratio(dc,dh);em=ratio(mc,mh)
        print(f'INTEGRATED: D_H={dh:.2f} D_HC={dc:.2f} MC_H={mh:.2f} MC_HC={mc:.2f}')
        print(f'            D/MC H={fmt(ratio(dh,mh))}, D/MC HC={fmt(ratio(dc,mc))}, data survival={fmt(ed)}, MC survival={fmt(em)}, double ratio={fmt(ratio(ed,em))}')
        if dc > dh*1.00001 or mc > mh*1.00001:
            print('WARNING: post-cut exposure-normalized yield exceeds pre-cut; check that production settings and file inputs match.')
        if not math.isclose(po['H','data'],po['HC','data'],rel_tol=5e-4):
            print('WARNING: data POT differs between productions; compare playlists and processing scope.')
        if not math.isclose(po['H','mc'],po['HC','mc'],rel_tol=5e-4):
            print('WARNING: MC POT differs between productions; compare playlists and processing scope.')
    finally:
        for f in files.values(): f.Close()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--key',default='Eel',help='1D histogram name, default Eel')
    ap.add_argument('--modes',nargs='+',choices=['FHC','RHC'],default=['FHC','RHC'])
    ap.add_argument('--csv',default=None)
    a=ap.parse_args()
    rows=[]
    for mode in a.modes: do_mode(mode,a.key,rows)
    if a.csv:
        os.makedirs(os.path.dirname(os.path.abspath(a.csv)),exist_ok=True)
        with open(a.csv,'w',newline='') as f:
            out=csv.DictWriter(f,fieldnames=list(rows[0]));out.writeheader();out.writerows(rows)
        print('\nCSV:',a.csv)
    print('\nInterpretation: double ratio < 1 means curvature lowers data/MC in that energy bin; > 1 raises it.')
    print('This tests incremental cut acceptance, not the correctness of the underlying MINOS charge/curvature calibration.')

if __name__=='__main__':main()

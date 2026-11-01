#!/usr/bin/env python3
"""Read-only LE5 standard-only benchmark against Maria's MINERvA EventSelection ROOT file.

Requirements: FNAL GPVM PyROOT/PlotUtils environment; no numpy or pandas.
Never includes the dedicated 2p2h production. The x and y axes are transposed.
"""
import argparse
import csv
import math
from pathlib import Path
import ROOT

ROOT.gROOT.SetBatch(True)
BASE = Path('/exp/minerva/data/users/qvuong/antinu_mu')
MARIA = Path('/exp/minerva/data/users/mmehmood/2026-10-24-NoNuEConstraint/MnvTunev1_minerva5_noNuE_Nov1/Hists_EventSelection_MnvTunev1_minerva5_noNuE_Nov1_minerva5_nosys_t99_z99_AntiNu.root')
TAG = 'ValidateCuts_CCnumubar_Inclusive_Curvature_v1'

COMPONENTS = {
    'QE': ('reco_QE_pZmu_pTmu', 'Lepton_Pt_Ppar_CCNuMuQE'),
    'RES': ('reco_resonant_pZmu_pTmu', 'Lepton_Pt_Ppar_CCNuMuDelta'),
    # Deliberately do not force a DIS match: Maria splits true/soft DIS,
    # and the truth definitions must be checked before this is valid.
}

def edges(h, axis):
    a = h.GetXaxis() if axis == 'x' else h.GetYaxis()
    return [float(a.GetBinLowEdge(k)) for k in range(1, a.GetNbins() + 2)]

def matching(m, u, tol=1e-8):
    ex, ey = edges(m, 'x'), edges(m, 'y')
    ux, uy = edges(u, 'x'), edges(u, 'y')
    return len(ex)==len(uy) and len(ey)==len(ux) and all(abs(x-y)<=tol for x,y in zip(ex,uy)) and all(abs(x-y)<=tol for x,y in zip(ey,ux))

def integral(h):
    return float(h.Integral(1,h.GetNbinsX(),1,h.GetNbinsY()))

def compare(label, m, u, rows, min_mc):
    if not m or not u:
        print('[SKIP] %s: a histogram is missing' % label)
        return
    if not (m.InheritsFrom('TH2') and u.InheritsFrom('TH2')):
        print('[SKIP] %s: not TH2 objects' % label)
        return
    if not matching(m,u):
        print('[FAIL] %s: bin edges do not match after swapping axes' % label)
        print(' Maria X/Y:', edges(m,'x'),edges(m,'y'))
        print(' Yours X/Y:', edges(u,'x'),edges(u,'y'))
        return
    im, iu = integral(m), integral(u)
    print('\n--- %s ---' % label)
    print('Maria visible=%.8f  yours visible=%.8f  yours/Maria=%.9f  difference=%.6f' %
          (im,iu,iu/im if im else float('nan'),iu-im))
    significant = []
    sumweighted = 0.
    weights = 0.
    maxall = None
    exact = 0
    bothzero = 0
    for i in range(1,m.GetNbinsX()+1):
        for j in range(1,m.GetNbinsY()+1):
            mv = float(m.GetBinContent(i,j))
            uv = float(u.GetBinContent(j,i))
            delta = uv - mv
            rel = delta/mv if mv != 0 else float('nan')
            if abs(delta)<1e-7:
                exact += 1
            if mv==0 and uv==0:
                bothzero += 1
            pzlo,pzhi = edges(m,'x')[i-1:i+1]
            ptlo,pthi = edges(m,'y')[j-1:j+1]
            rows.append(dict(group=label,pz_lo=pzlo,pz_hi=pzhi,pt_lo=ptlo,pt_hi=pthi,
                             maria=mv,ours=uv,delta=delta,relative=rel))
            if mv >= min_mc:
                significant.append((abs(rel),rel,pzlo,pzhi,ptlo,pthi,mv,uv))
                sumweighted += mv*rel*rel
                weights += mv
                if maxall is None or abs(rel)>abs(maxall[1]):
                    maxall = (pzlo,rel)
    print('2D bins=%d  exact (abs diff<1e-7)=%d  both zero=%d' % (m.GetNbinsX()*m.GetNbinsY(),exact,bothzero))
    print('Bins with Maria >= %.1f: %d; weighted RMS relative diff=%.4f%%' %
          (min_mc,len(significant),100*math.sqrt(sumweighted/weights) if weights else float('nan')))
    print('Top 12 per-bin fractional differences (Maria >= %.1f):' % min_mc)
    for _,rr,pl,ph,tl,th,mv,uv in sorted(significant,reverse=True)[:12]:
        print(' ppar [%5.2f,%5.2f]  pt [%5.3f,%5.3f]  Maria=%10.3f ours=%10.3f diff=%+7.2f%%' %
              (pl,ph,tl,th,mv,uv,100*rr))

def get2(f,key):
    obj=f.Get(key)
    if not obj:
        print('[MISSING]',f.GetName(), key)
    return obj

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--maria',default=str(MARIA))
    ap.add_argument('--our-data',default=str(BASE / ('kin_dist_datale5_'+TAG+'_MAD.root')))
    ap.add_argument('--our-mc',default=str(BASE / ('kin_dist_mcle5_'+TAG+'_MAD_no2p2h.root')))
    ap.add_argument('--out',default='maria_le5_2d_bin_comparison.csv')
    ap.add_argument('--min-mc',type=float,default=100.,help='Maria MC events threshold for largest relative differences')
    args=ap.parse_args()
    for p in (args.maria,args.our_data,args.our_mc):
        if not Path(p).is_file():
            ap.error('Missing input ROOT file: '+p)
    fm=ROOT.TFile.Open(args.maria)
    fd=ROOT.TFile.Open(args.our_data)
    fc=ROOT.TFile.Open(args.our_mc)
    for f in (fm,fd,fc):
        if not f or f.IsZombie():
            raise RuntimeError('Cannot open an input ROOT file')
    print('Maria file:',args.maria)
    print('Our DATA:',args.our_data)
    print('Our standard-only MC:',args.our_mc)
    print('NOTE: comparing RAW histogram yields; if MC POT differs, normalize before rate conclusions.')
    pars={}
    for name in ('DataPOT','MCPOT'):
        obj=fm.Get(name)
        if obj:
            pars[name]=float(obj.GetVal())
            print('Maria %s = %.12g' % (name,pars[name]))
    if len(pars)==2:
        print('Maria data/MC POT = %.12g' % (pars['DataPOT']/pars['MCPOT']))
    rows=[]
    compare('DATA',get2(fm,'h_data_pZmu_pTmu'),get2(fd,'Lepton_Pt_Ppar'),rows,20.)
    compare('MC_STANDARD',get2(fm,'h_mc_pZmu_pTmu'),get2(fc,'Lepton_Pt_Ppar'),rows,args.min_mc)
    for label,(mk,uk) in COMPONENTS.items():
        compare('MC_'+label,get2(fm,mk),get2(fc,uk),rows,args.min_mc)
    # Print a category closure test; do not assume same truth definitions.
    cats=['reco_QE','reco_2p2h','reco_resonant','reco_trueDIS','reco_softDIS','reco_other']
    tot=get2(fm,'h_mc_pZmu_pTmu')
    if tot:
        vals=[]
        for cat in cats:
            h=fm.Get(cat+'_pZmu_pTmu')
            vals.append((cat, integral(h) if h else float('nan')))
        print('\nMaria MC category totals (visible 2D):',vals)
        print('Maria category sum - h_mc = %.8f' % (sum(v for _,v in vals)-integral(tot)))
    with open(args.out,'w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=['group','pz_lo','pz_hi','pt_lo','pt_hi','maria','ours','delta','relative'])
        writer.writeheader()
        writer.writerows(rows)
    print('\nWrote',args.out)
    print('Do NOT compare raw total MC if production POT differs. No dedicated MEC included here.')
    # Keep ROOT ownership alive until script exit; avoid explicitly Close() with PyROOT custom objects.

if __name__=='__main__':
    main()

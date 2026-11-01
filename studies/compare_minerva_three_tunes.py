#!/usr/bin/env python3
"""Read-only three-tune MINERvA ROOT comparison. Requires PyROOT on GPVM.
Never manually Close() ROOT files (ROOT/cppyy ownership workaround).
"""
import argparse, csv, math, os
DEFAULTS={
 'FHC':('/exp/minerva/data/users/qvuong/nu_mu','leFHC','ValidateCuts_CCnumu_MyCuts'),
 'RHC':('/exp/minerva/data/users/qvuong/antinu_mu','le5','ValidateCuts_CCnumubar_MyCuts')}
TUNES=('v1','v4p3p0','v4p4p1')
KEEP=[]
def op(ROOT,p):
 f=ROOT.TFile.Open(p,'READ')
 if not f or f.IsZombie(): raise RuntimeError('Cannot open '+p)
 KEEP.append(f);return f
def get(f,k,required=True):
 h=f.Get(k)
 if not h or not h.InheritsFrom('TH1') or h.GetDimension()!=1:
  if required: raise RuntimeError('Missing histogram {} in {}'.format(k,f.GetName()))
  return None
 return h
def edges(h):return [h.GetXaxis().GetBinLowEdge(i) for i in range(1,h.GetNbinsX()+2)]
def vals(h):return [float(h.GetBinContent(i)) for i in range(1,h.GetNbinsX()+1)]
def safe(a,b):return a/b if b else float('nan')
def compatible(a,b):
 x,y=edges(a),edges(b)
 return len(x)==len(y) and all(abs(i-j)<=1.e-8*max(1.,abs(i),abs(j)) for i,j in zip(x,y))
def report_hist(label,hd,mc,sc,w,beam,cat='inclusive'):
 ref=hd[TUNES[0]]
 for t in TUNES:
  if not compatible(ref,hd[t]) or not compatible(ref,mc[t]):
   raise RuntimeError('Binning mismatch for {} {}'.format(beam,label))
 for t in TUNES[1:]:
  deviation=max((abs(a-b) for a,b in zip(vals(ref),vals(hd[t]))),default=0)
  if deviation>1e-8: print('WARNING: DATA changed {} {} max bin difference {}'.format(beam,t,deviation))
 for i,(lo,hi) in enumerate(zip(edges(ref)[:-1],edges(ref)[1:]),1):
  d=ref.GetBinContent(i); m={t:mc[t].GetBinContent(i) for t in TUNES}
  w.writerow(dict(beam=beam,category=cat,hist=label,low_GeV=lo,high_GeV=hi,data=d,
   **{'mc_'+t:m[t] for t in TUNES},
   v430_over_v1=safe(m['v4p3p0'],m['v1']),v441_over_v1=safe(m['v4p4p1'],m['v1']),
   v441_over_v430=safe(m['v4p4p1'],m['v4p3p0']),
   **{'data_mc_'+t:safe(d,m[t]*sc[t]) for t in TUNES}))
  if cat=='inclusive':
   print('{:>9.2f}–{:<6.2f} {:>11.1f} {:>11.1f} {:>11.1f} {:>8.4f} {:>8.4f} {:>8.4f}  {:>7.3f} {:>7.3f} {:>7.3f}'.format(lo,hi,*(m[t] for t in TUNES),safe(m['v4p3p0'],m['v1']),safe(m['v4p4p1'],m['v1']),safe(m['v4p4p1'],m['v4p3p0']),*(safe(d,m[t]*sc[t]) for t in TUNES)))
 totals={t:sum(vals(mc[t])) for t in TUNES}
 return totals

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--beam',choices=['both','FHC','RHC'],default='both')
 p.add_argument('--hist',default='Eel')
 p.add_argument('--categories',action='store_true')
 p.add_argument('--out',default='three_tune_comparison.csv')
 p.add_argument('--scale-fhc',type=float,default=1.,help='MC to data POT scale for FHC (default raw=1)')
 p.add_argument('--scale-rhc',type=float,default=1.,help='MC to data POT scale for RHC (default raw=1)')
 p.add_argument('--all-components',action='store_true',help='Print bin-by-bin category ratios too, not only totals')
 a=p.parse_args()
 import ROOT
 ROOT.gROOT.SetBatch(True)
 fields=['beam','category','hist','low_GeV','high_GeV','data']+['mc_'+x for x in TUNES]+['v430_over_v1','v441_over_v1','v441_over_v430']+['data_mc_'+x for x in TUNES]
 with open(a.out,'w',newline='') as output:
  w=csv.DictWriter(output,fieldnames=fields);w.writeheader()
  for beam in (('FHC','RHC') if a.beam=='both' else (a.beam,)):
   folder,playlist,stem=DEFAULTS[beam]
   files={(t,k):op(ROOT,os.path.join(folder,'kin_dist_{}{}_{}_{}_MAD.root'.format(k,playlist,stem,t))) for t in TUNES for k in ('mc','data')}
   hd={t:get(files[t,'data'],a.hist) for t in TUNES}
   mc={t:get(files[t,'mc'],a.hist) for t in TUNES}
   scale=a.scale_fhc if beam=='FHC' else a.scale_rhc
   scales={t:scale for t in TUNES}
   print('\n'+ '='*112+'\n{} {} | MC-to-data scale {}{}'.format(beam,a.hist,scale,' (RAW: D/MC not physical)' if scale==1 else ''))
   print('       E mu             V1         v4.3.0       v4.4.1  430/V1    441/V1    441/430   D/MC V1   D/MC 430  D/MC 441')
   tot=report_hist(a.hist,hd,mc,scales,w,beam)
   print('TOTAL MC: '+ ' | '.join('{}={:.2f}'.format(t,tot[t]) for t in TUNES))
   print('TOTAL ratios 430/V1={:.5f} 441/V1={:.5f} 441/430={:.5f}'.format(safe(tot['v4p3p0'],tot['v1']),safe(tot['v4p4p1'],tot['v1']),safe(tot['v4p4p1'],tot['v4p3p0'])))
   if a.categories:
    keys=[set(k.GetName() for k in files[t,'mc'].GetListOfKeys()) for t in TUNES]
    cats=sorted(k for k in set.intersection(*keys) if k.startswith(a.hist+'_') and '_CV_WithErr' not in k and all(get(files[t,'mc'],k,False) for t in TUNES))
    print('\nCategory totals (raw MC):')
    print('{:<32} {:>12} {:>12} {:>12} {:>9} {:>9} {:>9}'.format('category','V1','v4.3.0','v4.4.1','430/V1','441/V1','441/430'))
    for cat in cats:
     cm={t:get(files[t,'mc'],cat) for t in TUNES}
     if any(not compatible(mc['v1'],h) for h in cm.values()):print('SKIP category binning mismatch',cat);continue
     sums={t:sum(vals(cm[t])) for t in TUNES}
     print('{:<32} {:>12.1f} {:>12.1f} {:>12.1f} {:>9.4f} {:>9.4f} {:>9.4f}'.format(cat.replace(a.hist+'_',''),*(sums[t] for t in TUNES),safe(sums['v4p3p0'],sums['v1']),safe(sums['v4p4p1'],sums['v1']),safe(sums['v4p4p1'],sums['v4p3p0'])))
     # per-category bin data/MC is NOT meaningful, so CSV fields are left blank for data ratios
     for i,(lo,hi) in enumerate(zip(edges(cm['v1'])[:-1],edges(cm['v1'])[1:]),1):
      v={t:float(cm[t].GetBinContent(i)) for t in TUNES}
      w.writerow(dict(beam=beam,category=cat,hist=a.hist,low_GeV=lo,high_GeV=hi,data='',
       **{'mc_'+t:v[t] for t in TUNES},v430_over_v1=safe(v['v4p3p0'],v['v1']),v441_over_v1=safe(v['v4p4p1'],v['v1']),v441_over_v430=safe(v['v4p4p1'],v['v4p3p0']),
       **{'data_mc_'+t:'' for t in TUNES}))
     if a.all_components:
      for i,(lo,hi) in enumerate(zip(edges(cm['v1'])[:-1],edges(cm['v1'])[1:]),1):
       v=[float(cm[t].GetBinContent(i)) for t in TUNES]
       print('  {:>4g}–{:<4g} yields {} ratios 430/V1={:.4f} 441/V1={:.4f} 441/430={:.4f}'.format(lo,hi,','.join('{:.3g}'.format(x) for x in v),safe(v[1],v[0]),safe(v[2],v[0]),safe(v[2],v[1])))
 print('\nWrote',os.path.abspath(a.out))
if __name__=='__main__':main()

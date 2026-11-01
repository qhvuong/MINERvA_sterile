#!/usr/bin/env python3
"""Read-only inventory and Maria LE1 cross-production comparison (PyROOT)."""
import ROOT, os, math
ROOT.gROOT.SetBatch(True)
BASE='/exp/minerva/data/users/mmehmood/2026-10-24-NoNuEConstraint'
SETS={
 'Maria LE5': ('MnvTunev1_minerva5_noNuE_Nov1','minerva5','AntiNu'),
 'Maria LE1 Oct23': ('MnvTunev1_minerva1_NoNuE_v5_2026-10-23','minerva1','Nu'),
 'Maria LE1 Oct27': ('MnvTunev1_minerva1_NoNuE_v5_2026-10-27','minerva1','Nu'),
}
USER={
 'Our LE5 DATA':'/exp/minerva/data/users/qvuong/antinu_mu/kin_dist_datale5_ValidateCuts_CCnumubar_Inclusive_Curvature_v1_MAD.root',
 'Our LE5 standard MC':'/exp/minerva/data/users/qvuong/antinu_mu/kin_dist_mcle5_ValidateCuts_CCnumubar_Inclusive_Curvature_v1_MAD_no2p2h.root',
 'Our LE5 combined MC':'/exp/minerva/data/users/qvuong/antinu_mu/kin_dist_mcle5_ValidateCuts_CCnumubar_Inclusive_Curvature_v1_MAD.root',
 'Our LE1 DATA':'/exp/minerva/data/users/qvuong/nu_mu/kin_dist_datale1_Maria_LE1_p6_CV_MAD.root',
 'Our LE1 MC':'/exp/minerva/data/users/qvuong/nu_mu/kin_dist_mcle1_Maria_LE1_p6_CV_MAD.root',
}

def op(path):
 if not os.path.isfile(path): print('MISSING',path); return None
 f=ROOT.TFile.Open(path)
 if not f or f.IsZombie(): print('BAD ROOT',path); return None
 return f

def integral(h):
 if h.InheritsFrom('TH2'): return h.Integral(1,h.GetNbinsX(),1,h.GetNbinsY())
 return h.Integral(1,h.GetNbinsX())

def edges(axis):
 return [round(float(axis.GetBinLowEdge(i)),6) for i in range(1,axis.GetNbins()+2)]

def info(f,label,is_maria):
 print('\n'+'='*85+'\n'+label+'\n'+f.GetName())
 for key in ('DataPOT','MCPOT','DataRecoCuts','MCReco_RecoCuts','DataPassThetaCut','MCRecoPassThetaCut'):
  o=f.Get(key)
  if o and hasattr(o,'GetVal'): print(' ',key,'=',o.GetVal())
 for k in f.GetListOfKeys():
  o=f.Get(k.GetName()); name=k.GetName()
  if not o or not o.InheritsFrom('TH1'):continue
  wanted= name in ('h_mc_pZmu','h_mc_pTmu','h_data_pZmu','h_data_pTmu','h_mc_pZmu_pTmu','h_data_pZmu_pTmu','reco_2p2h_pZmu_pTmu') if is_maria else any(t in name.lower() for t in ('pz','pt','parallel','transverse','pot','eel'))
  if not wanted:continue
  shape=f'{o.GetNbinsX()}x{o.GetNbinsY()}' if o.InheritsFrom('TH2') else str(o.GetNbinsX())
  print(f'  {name:46s} {o.ClassName():18s} bins={shape:8s} integral={integral(o):.8g}')
  if (is_maria and name in ('h_mc_pZmu_pTmu','h_data_pZmu_pTmu')) or (not is_maria and o.InheritsFrom('TH2') and any(t in name.lower() for t in ('pz','pt','parallel','transverse'))):
   print('      X edges:',edges(o.GetXaxis()))
   print('      Y edges:',edges(o.GetYaxis()))
 return

def le1comparison(a,b):
 print('\n'+'='*85+'\nMARIA LE1 Oct23 vs Oct27 BIN COMPARISON')
 for name in ('h_data_pZmu','h_mc_pZmu','h_data_pTmu','h_mc_pTmu','h_data_pZmu_pTmu','h_mc_pZmu_pTmu'):
  x=a.Get(name);y=b.Get(name)
  if not x or not y:continue
  shape_eq=edges(x.GetXaxis())==edges(y.GetXaxis()) and (not x.InheritsFrom('TH2') or edges(x.GetYaxis())==edges(y.GetYaxis()))
  va,vb=integral(x),integral(y)
  print(f'{name:25s} oct23={va:.8g} oct27={vb:.8g} ratio={vb/va if va else float("nan"):.7g} bins_match={shape_eq}')
  if shape_eq:
   diffs=[]
   for ix in range(1,x.GetNbinsX()+1):
    for iy in range(1,(x.GetNbinsY()+1) if x.InheritsFrom('TH2') else 2):
     v1=x.GetBinContent(ix,iy) if x.InheritsFrom('TH2') else x.GetBinContent(ix)
     v2=y.GetBinContent(ix,iy) if y.InheritsFrom('TH2') else y.GetBinContent(ix)
     diffs.append(abs(v1-v2))
   print(f'  max absolute bin difference = {max(diffs):.7g}')

opened={}
for label,(folder,playlist,beam) in SETS.items():
 path=os.path.join(BASE,folder,f'Hists_EventSelection_{folder}_{playlist}_nosys_t99_z99_{beam}.root')
 f=op(path)
 if f: info(f,label,True);opened[label]=f
if 'Maria LE1 Oct23' in opened and 'Maria LE1 Oct27' in opened:le1comparison(opened['Maria LE1 Oct23'],opened['Maria LE1 Oct27'])
for label,path in USER.items():
 f=op(path)
 if f: info(f,label,False);f.Close()
for f in opened.values():f.Close()

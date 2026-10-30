#!/usr/bin/env python3
"""Compare same-sample ROOT histograms with/without MINOS CV factor."""
import argparse
import ROOT
ROOT.gROOT.SetBatch(True)
p=argparse.ArgumentParser()
p.add_argument('--on',required=True)
p.add_argument('--off',required=True)
p.add_argument('--keys',nargs='+',default=['Eel','Eel_CCNuMu2p2h','Elepton_Eavail'])
a=p.parse_args()
f1=ROOT.TFile.Open(a.on)
f0=ROOT.TFile.Open(a.off)
for tag,f in [('on',f1),('off',f0)]:
    if not f or f.IsZombie(): raise RuntimeError('Cannot open %s file' % tag)
print('File keys: on=%d off=%d' % (f1.GetListOfKeys().GetSize(),f0.GetListOfKeys().GetSize()))
for key in a.keys:
    h1=f1.Get(key); h0=f0.Get(key)
    if not h1 or not h0:
        print('MISSING',key,'on=',bool(h1),'off=',bool(h0));continue
    if not h1.InheritsFrom('TH1') or not h0.InheritsFrom('TH1'):
        print('SKIP non-TH1',key);continue
    if h1.GetNcells()!=h0.GetNcells():
        print('INCOMPATIBLE',key);continue
    s1=sum(h1.GetBinContent(i) for i in range(h1.GetNcells()))
    s0=sum(h0.GetBinContent(i) for i in range(h0.GetNcells()))
    diffs=[(i,h1.GetBinContent(i)-h0.GetBinContent(i)) for i in range(h1.GetNcells())]
    print('%-30s on=%12.5f off=%12.5f on/off=%9.5f maxabsdiff=%g changedbins=%d' %
        (key,s1,s0,s1/s0 if s0 else float('nan'), max(abs(v) for _,v in diffs),sum(v!=0 for _,v in diffs)))
    if h1.GetDimension()==1:
        for i in range(1,h1.GetNbinsX()+1):
            u=h1.GetBinContent(i);v=h0.GetBinContent(i)
            print('  bin [%6g, %6g): on=%10.3f off=%10.3f ratio=%7.4f' %
                  (h1.GetXaxis().GetBinLowEdge(i),h1.GetXaxis().GetBinUpEdge(i),u,v,u/v if v else float('nan')))
f1.Close();f0.Close()

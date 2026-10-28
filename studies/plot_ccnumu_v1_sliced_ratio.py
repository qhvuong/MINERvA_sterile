#!/usr/bin/env python3
"""CCnumu v1: data/MC ratios in pT and p_parallel slices (ROOT/PyROOT).
Run on the MINERvA GPVM with the usual ROOT/PlotUtils environment.
"""
import argparse
import math
import os
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetEndErrorSize(4)

DATA = '/exp/minerva/data/users/qvuong/nu_mu/kin_dist_dataleFHC_ValidateDataMC_CCnumu_myInclusive_v1_MAD.root'
MC = '/exp/minerva/data/users/qvuong/nu_mu/kin_dist_mcleFHC_ValidateDataMC_CCnumu_myInclusive_v1_MAD.root'


def open_file(path):
    f = ROOT.TFile.Open(path, 'READ')
    if not f or f.IsZombie():
        raise RuntimeError('Cannot open ROOT file: ' + path)
    return f


def print_2d_keys(f):
    result = []
    for key in f.GetListOfKeys():
        obj = key.ReadObj()
        if obj.InheritsFrom('TH2'):
            result.append(key.GetName())
    return sorted(result)


def choose_pair(fd, fm, requested):
    kd = set(print_2d_keys(fd))
    km = set(print_2d_keys(fm))
    shared = sorted(kd & km)
    if requested:
        if requested not in shared:
            raise RuntimeError('Histogram %r not shared by both files. Available shared TH2 keys:\n%s' % (requested, '\n'.join(shared)))
        return requested
    candidates = [k for k in shared if ('Ppar_Pt' in k or 'Pt_Ppar' in k or 'Ppar_Pt' in k.replace('p','P'))]
    candidates = [k for k in candidates if not any(k.endswith('_'+s) for s in ('CCNuMuQE','CCNuMuRES','CCNuMuDIS','Other','CCPi','NC','CCNuMu2p2h'))]
    exact = [k for k in candidates if k in ('Lepton_Ppar_Pt', 'Lepton_Pt_Ppar')]
    if len(exact) == 1:
        return exact[0]
    if len(candidates) == 1:
        return candidates[0]
    raise RuntimeError('Cannot uniquely select total 2D histogram. Pass --hist KEY.\nShared 2D keys:\n'+'\n'.join(shared))


def central_th2(obj, name):
    if hasattr(obj, 'GetCVHistoWithStatError'):
        h = obj.GetCVHistoWithStatError()
    else:
        h = obj.Clone(name+'_raw')
    h.SetName(name)
    h.SetDirectory(0)
    h.Sumw2()
    return h


def pot(f):
    meta = f.Get('Meta')
    if not meta:
        return None
    branch = meta.GetBranch('POT_Used')
    if not branch:
        return None
    total = 0.
    for ent in meta:
        total += float(ent.POT_Used)
    return total if math.isfinite(total) and total > 0 else None


def axis_edges(axis):
    return [axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins()+2)]


def hist_projection(h, along_x, bin_number, name):
    if along_x:
        out = h.ProjectionX(name, bin_number, bin_number, 'e')
    else:
        out = h.ProjectionY(name, bin_number, bin_number, 'e')
    out.SetDirectory(0)
    return out


def draw_pages(data, mc, axis_x, slice_by_x, outdir, label, ratio_min, ratio_max):
    """slice_by_x=True -> project Y in each X bin; False -> project X in each Y bin."""
    slice_axis = data.GetXaxis() if slice_by_x else data.GetYaxis()
    nb = slice_axis.GetNbins()
    nc = 3
    nr = math.ceil(nb/nc)
    c = ROOT.TCanvas('c_'+label, label, 460*nc, 330*nr)
    c.Divide(nc, nr, .008, .015)
    keep = []
    for i in range(1, nb+1):
        c.cd(i)
        pad = ROOT.gPad
        pad.SetTopMargin(.13)
        pad.SetBottomMargin(.17)
        pad.SetLeftMargin(.17)
        pad.SetRightMargin(.04)
        d = hist_projection(data, not slice_by_x, i, 'data_%s_%d'%(label,i))
        m = hist_projection(mc, not slice_by_x, i, 'mc_%s_%d'%(label,i))
        ratio = d.Clone('ratio_%s_%d'%(label,i))
        ratio.SetDirectory(0)
        ratio.Divide(m)  # ROOT statistical propagation; zero denominator bins give zero
        ratio.SetTitle('')
        ratio.SetMinimum(ratio_min)
        ratio.SetMaximum(ratio_max)
        ratio.GetYaxis().SetTitle('Data / MC')
        ratio.GetYaxis().SetTitleSize(.048)
        ratio.GetYaxis().SetLabelSize(.045)
        ratio.GetYaxis().SetTitleOffset(1.3)
        ratio.GetXaxis().SetTitle('p_{T}^{#mu} (GeV/c)' if slice_by_x == axis_x else 'p_{||}^{#mu} (GeV/c)')
        ratio.GetXaxis().SetTitleSize(.048)
        ratio.GetXaxis().SetLabelSize(.045)
        ratio.GetXaxis().SetTitleOffset(1.18)
        ratio.SetMarkerStyle(20)
        ratio.SetMarkerSize(.75)
        ratio.SetLineColor(ROOT.kBlack)
        ratio.Draw('E1')
        line = ROOT.TLine(ratio.GetXaxis().GetXmin(),1.,ratio.GetXaxis().GetXmax(),1.)
        line.SetLineStyle(2)
        line.SetLineColor(ROOT.kRed+1)
        line.Draw('SAME')
        lo,hi = slice_axis.GetBinLowEdge(i),slice_axis.GetBinUpEdge(i)
        tex = ROOT.TLatex()
        tex.SetNDC(True)
        tex.SetTextFont(42)
        tex.SetTextSize(.055)
        sym = 'p_{||}^{#mu}' if slice_by_x == axis_x else 'p_{T}^{#mu}'
        tex.DrawLatex(.19,.90, '%g < %s < %g'%(lo,sym,hi))
        keep.extend([d,m,ratio,line,tex])
    c.cd()
    base = os.path.join(outdir, label)
    c.SaveAs(base+'.pdf')
    c.SaveAs(base+'.png')
    return base


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--data', default=DATA)
    p.add_argument('--mc', default=MC)
    p.add_argument('--hist', help='Exact shared total TH2 histogram key, e.g. Lepton_Ppar_Pt')
    p.add_argument('--outdir', default='ccnumu_v1_sliced_ratios')
    p.add_argument('--ratio-min', type=float, default=0.4)
    p.add_argument('--ratio-max', type=float, default=1.6)
    p.add_argument('--scale', type=float, help='Override MC scale (otherwise data POT / MC POT)')
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    fd,fm = open_file(args.data),open_file(args.mc)
    key = choose_pair(fd,fm,args.hist)
    rawd,rawm = fd.Get(key),fm.Get(key)
    d,m = central_th2(rawd,'d2'),central_th2(rawm,'m2')
    dx,dy = axis_edges(d.GetXaxis()),axis_edges(d.GetYaxis())
    if dx != axis_edges(m.GetXaxis()) or dy != axis_edges(m.GetYaxis()):
        raise RuntimeError('Data and MC histogram binning differs.')
    # Infer orientation from key; don't silently assume axes.
    if 'Ppar_Pt' in key:
        x_is_ppar = True
    elif 'Pt_Ppar' in key:
        x_is_ppar = False
    else:
        raise RuntimeError('Unknown histogram orientation in key %r; use a key containing Ppar_Pt or Pt_Ppar.'%key)
    dp,mp = pot(fd),pot(fm)
    if args.scale is None:
        if dp is None or mp is None:
            raise RuntimeError('Missing Meta/POT_Used. Supply --scale or inspect the files.')
        scale = dp/mp
    else:
        scale = args.scale
    m.Scale(scale)
    print('Histogram:',key,' x=p_parallel:',x_is_ppar)
    print('Data POT:',dp,' MC POT:',mp,' MC scaling:',scale)
    print('X edges:',dx)
    print('Y edges:',dy)
    print('Data events (weighted integral):',d.Integral(),'MC (POT-scaled):',m.Integral())
    # For pT in p_parallel slices: slice X if x is p_parallel.
    a = draw_pages(d,m,x_is_ppar, x_is_ppar,args.outdir,'pt_in_ppar_slices',args.ratio_min,args.ratio_max)
    # For p_parallel in pT slices: slice X if x is pT.
    b = draw_pages(d,m,x_is_ppar,not x_is_ppar,args.outdir,'ppar_in_pt_slices',args.ratio_min,args.ratio_max)
    print('Saved:',a+'.pdf/png',b+'.pdf/png')
    fd.Close();fm.Close()

if __name__ == '__main__':
    main()

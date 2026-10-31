#!/usr/bin/env python3
"""Audit LE5 standard + dedicated MEC merging, without modifying ROOT files.
Run in GPVM environment with PyROOT. Outputs a CSV and human-readable tables.
"""
import argparse, csv, math, os
import ROOT
ROOT.gROOT.SetBatch(True)

DEFAULT_TAG = 'ValidateCuts_CCnumubar_Inclusive_noLeptonAngle_v1'
DEFAULT_DIR = '/exp/minerva/data/users/qvuong/antinu_mu'
POT_STD = 9.513217104576512e20
POT_SPECIAL = 2.88863442755131e21


def read_hist(path, name):
    f = ROOT.TFile.Open(path, 'READ')
    if not f or f.IsZombie():
        raise RuntimeError('Cannot open ROOT file: ' + path)
    h = f.Get(name)
    if not h or not h.InheritsFrom('TH1') or h.GetDimension() != 1:
        names = [k.GetName() for k in f.GetListOfKeys() if k.GetName().startswith('Eel')]
        raise RuntimeError('Histogram %r missing from %s; Eel histograms: %s' % (name, path, names))
    # Extract values while file is open, no ownership or lifetime ambiguity.
    n = h.GetNbinsX()
    edges = [float(h.GetXaxis().GetBinLowEdge(i)) for i in range(1, n+2)]
    values = [float(h.GetBinContent(i)) for i in range(1, n+1)]
    errs = [float(h.GetBinError(i)) for i in range(1, n+1)]
    f.Close()
    return edges, values, errs


def same_edges(a, b):
    return len(a) == len(b) and all(math.isclose(x, y, abs_tol=1e-9) for x, y in zip(a,b))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dir', default=DEFAULT_DIR)
    ap.add_argument('--tag', default=DEFAULT_TAG)
    ap.add_argument('--hist', default='Eel', help='e.g. Eel or Eel_CCNuMu2p2h')
    ap.add_argument('--mc-scale', type=float, default=0.108, help='Data POT / standard MC POT; only used for optional data comparison')
    ap.add_argument('--out', default='/exp/minerva/data/users/qvuong/studies/le5_mec_merge_check.csv')
    args = ap.parse_args()
    base = 'kin_dist_mcle5_' + args.tag + '_MAD'
    paths = {
        'standard': os.path.join(args.dir, base + '_no2p2h.root'),
        'special': os.path.join(args.dir, base + '_only2p2h.root'),
        'combined': os.path.join(args.dir, base + '.root'),
    }
    for label, path in paths.items():
        if not os.path.isfile(path):
            raise SystemExit('Missing %s file: %s' % (label,path))
    arrays = {k: read_hist(v,args.hist) for k,v in paths.items()}
    edges = arrays['standard'][0]
    for k in arrays:
        if not same_edges(edges, arrays[k][0]):
            raise RuntimeError('Histogram binning differs for '+k)
    s=arrays['standard'][1]; p=arrays['special'][1]; c=arrays['combined'][1]
    scale=POT_STD/POT_SPECIAL
    tol=1e-5
    print('Histogram:',args.hist)
    print('Dedicated special POT scale: %.12g' % scale)
    print('NOTE: standard-only file is named no2p2h but may contain MEC events.\n')
    print('%13s %13s %13s %13s %13s %13s %11s' % ('Bin GeV','Standard','Special raw','Special scaled','Combined','Delta merge','Special %'))
    rows=[]
    for i,(lo,hi) in enumerate(zip(edges[:-1],edges[1:])):
        expected=s[i]+scale*p[i]
        residual=c[i]-expected
        fraction=scale*p[i]/c[i] if c[i] else float('nan')
        rows.append(dict(lo=lo,hi=hi,standard=s[i],special_raw=p[i],special_scaled=scale*p[i],combined=c[i],expected=expected,merge_residual=residual,special_fraction=fraction))
        print('%5.1f-%-6.1f %13.5g %13.5g %13.5g %13.5g %13.4g %10.2f%%' % (lo,hi,s[i],p[i],scale*p[i],c[i],residual,100*fraction))
    total=[sum(s),sum(p),sum(c)]
    predicted=total[0]+scale*total[1]
    print('\nIntegrals (visible bins): standard %.9g, special raw %.9g, special scaled %.9g, combined %.9g' % (total[0],total[1],scale*total[1],total[2]))
    print('Combined - (standard + scaled special) = %.9g (relative %.4g)' % (total[2]-predicted,(total[2]-predicted)/max(abs(total[2]),1e-12)))
    max_rel=max(abs(r['merge_residual'])/max(abs(r['combined']),1e-12) for r in rows)
    print('Largest relative bin mismatch: %.5g' % max_rel)
    if max_rel>tol:
        print('WARNING: Bin closure exceeds tolerance; inspect merge code and POT scaling before interpreting physics.')
    else:
        print('PASS: Bin-by-bin merge closes at relative tolerance %.0g' % tol)
    if args.hist=='Eel':
        data_path=os.path.join(args.dir,'kin_dist_datale5_'+args.tag+'_MAD.root')
        if os.path.isfile(data_path):
            de,d,_=read_hist(data_path,args.hist)
            if same_edges(edges,de):
                print('\nData comparison using MC normalization %.8g:' % args.mc_scale)
                print('%13s %14s %14s %14s' % ('Bin GeV','Data','D/MC full','D/MC standard'))
                for i,(lo,hi) in enumerate(zip(edges[:-1],edges[1:])):
                    print('%5.1f-%-6.1f %14.5g %14.4f %14.4f' % (lo,hi,d[i],d[i]/(args.mc_scale*c[i]) if c[i]>0 else float('nan'),d[i]/(args.mc_scale*s[i]) if s[i]>0 else float('nan')))
            else:
                print('Data bins do not match MC, omitting D/MC.')
    with open(args.out,'w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    print('\nWrote:',args.out)

if __name__=='__main__': main()

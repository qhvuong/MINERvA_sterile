#!/usr/bin/env python3
"""Audit combined MINERvA LE5 MC against nominal and special 2p2h ROOT outputs.

Run on a MINERvA GPVM in the normal ROOT/PlotUtils environment.
Read-only: does not modify any input ROOT file.
"""
import argparse
import math
import os
import sys
from collections import Counter

import ROOT
ROOT.gROOT.SetBatch(True)

DEFAULT_DIR = '/exp/minerva/data/users/qvuong/antinu_mu'
DEFAULT_TAG = 'ValidateCuts_CCnumubar_Inclusive_Curvature_v1'
DEFAULT_RATIO = 9.513217104576512e20 / 2.88863442755131e21


def open_root(path):
    f = ROOT.TFile.Open(path, 'READ')
    if not f or f.IsZombie():
        raise RuntimeError('Cannot open ROOT file: ' + path)
    return f


def get_hists(f):
    """Return key names for histograms (1D/2D/3D); ignore TTrees and metadata."""
    result = {}
    for key in f.GetListOfKeys():
        name = key.GetName()
        obj = f.Get(name)
        if obj and obj.InheritsFrom('TH1'):
            result[name] = obj
    return result


def cells(h):
    # ROOT's global bin indices include flow bins, including TH2/TH3 corners.
    return range(h.GetNcells())


def integral(h, flow=False):
    if not flow:
        if h.InheritsFrom('TH3'):
            return h.Integral(1, h.GetNbinsX(), 1, h.GetNbinsY(), 1, h.GetNbinsZ())
        if h.InheritsFrom('TH2'):
            return h.Integral(1, h.GetNbinsX(), 1, h.GetNbinsY())
        return h.Integral(1, h.GetNbinsX())
    return sum(h.GetBinContent(i) for i in cells(h))


def geometry(h):
    return (h.GetDimension(), h.GetNbinsX(), h.GetNbinsY(), h.GetNbinsZ(), h.GetNcells())


def bin_label(h, b):
    if h.GetDimension() == 1:
        if b == 0:
            return 'underflow'
        if b == h.GetNbinsX() + 1:
            return 'overflow'
        ax = h.GetXaxis()
        return '[{:.4g}, {:.4g})'.format(ax.GetBinLowEdge(b), ax.GetBinUpEdge(b))
    return 'global bin {}'.format(b)


def is_analysis_hist(name):
    low = name.lower()
    return not any(s in low for s in ('cutflow', 'cut_flow', 'pot', 'metadata'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dir', default=DEFAULT_DIR)
    p.add_argument('--tag', default=DEFAULT_TAG)
    p.add_argument('--playlist', default='le5')
    p.add_argument('--pot-nominal', type=float, default=9.513217104576512e20)
    p.add_argument('--pot-special', type=float, default=2.88863442755131e21)
    p.add_argument('--scale', type=float, default=None, help='Override special POT scaling')
    p.add_argument('--filter', default='', help='Substring to restrict printed histogram diagnostics')
    p.add_argument('--top', type=int, default=25)
    p.add_argument('--list-keys', action='store_true', help='Print all histogram names matching --filter')
    p.add_argument('--csv', default='audit_ccnumubar_2p2h.csv')
    args = p.parse_args()
    stem = os.path.join(args.dir, 'kin_dist_mc{}_{}_MAD'.format(args.playlist, args.tag))
    paths = {'nominal': stem+'_no2p2h.root', 'special': stem+'_only2p2h.root', 'combined': stem+'.root'}
    for k, path in paths.items():
        print('{}: {} ({} bytes)'.format(k, path, os.path.getsize(path) if os.path.exists(path) else 'MISSING'))
    files = {k: open_root(v) for k, v in paths.items()}
    hists = {k: get_hists(f) for k, f in files.items()}
    scale = args.scale if args.scale is not None else args.pot_nominal/args.pot_special
    print('\nSpecial scale = {:.12g}; nominal POT = {:.9e}; special POT = {:.9e}'.format(scale,args.pot_nominal,args.pot_special))
    for k, h in hists.items():
        print('{}: {} histogram keys'.format(k, len(h)))
    names = {k: set(h) for k,h in hists.items()}
    shared = sorted(set.intersection(*names.values()))
    print('Common histogram keys: {}'.format(len(shared)))
    for k in paths:
        only = sorted(names[k] - set.union(*(names[other] for other in paths if other != k)))
        print('Keys exclusive to {} ({}): {}'.format(k, len(only), only[:10]))
    print('Keys in nominal but absent from combined:', sorted(names['nominal']-names['combined'])[:20])
    print('Keys in special but absent from combined:', sorted(names['special']-names['combined'])[:20])

    # Audit component addition for all common histograms, including flow bins.
    bad, checked, maximum = [], 0, (0.0, '')
    for name in shared:
        a, b, c = (hists[k][name] for k in ('nominal','special','combined'))
        if not (geometry(a)==geometry(b)==geometry(c)):
            bad.append((name, 'geometry mismatch'))
            continue
        checked += 1
        max_err = max(abs(c.GetBinContent(i) - a.GetBinContent(i) - scale*b.GetBinContent(i)) for i in cells(c))
        size = max(1., max(abs(c.GetBinContent(i)) for i in cells(c)))
        if max_err > 1e-6*size:
            bad.append((name, 'max bin residual {:.6g} (relative {:.3g})'.format(max_err,max_err/size)))
        if max_err > maximum[0]: maximum = (max_err, name)
    print('\n=== ADDITION CLOSURE: combined = nominal + scale * special ===')
    print('Compatible/common histograms checked: {}; failing: {}'.format(checked,len(bad)))
    print('Largest absolute bin discrepancy: {:.6g} ({})'.format(*maximum))
    for item in bad[:args.top]: print('  FAIL:', *item)
    if len(bad)>args.top: print('  ... {} additional failures'.format(len(bad)-args.top))

    # Report all analysis histograms with nonzero special yield, including all shared bins.
    rows=[]
    for name in shared:
        if args.filter.lower() not in name.lower(): continue
        a,b,c=(hists[k][name] for k in ('nominal','special','combined'))
        if geometry(a)!=geometry(b) or geometry(a)!=geometry(c): continue
        nom=integral(a)
        spec=integral(b)*scale
        tot=integral(c)
        frac=spec/tot if abs(tot)>1e-12 else float('nan')
        excess=tot-nom
        rows.append((name, nom, spec, tot, frac, excess, a.GetDimension()))
    rows.sort(key=lambda r: -abs(r[2]))
    if args.list_keys:
        print('\n=== MATCHED HISTOGRAM KEYS ===')
        for row in sorted(rows):
            print(row[0])
        print('Matched {} histograms for filter {!r}'.format(len(rows), args.filter))
    print('\n=== SPECIAL-SAMPLE CONTRIBUTION (hist integrals, excluding flow) ===')
    print('NOTE: Histograms can include overlapping samples, universes, cut-flow, and truth categories; do not sum across names.')
    print('{:<75} {:>11} {:>11} {:>11} {:>10}'.format('Histogram', 'Nominal', 'Special*POT', 'Combined', 'Frac'))
    for name,nom,spec,tot,frac,excess,dim in rows[:args.top]:
        print('{:<75} {:>11.5g} {:>11.5g} {:>11.5g} {:>10.3%}'.format(name[:75],nom,spec,tot,frac))

    import csv
    with open(args.csv, 'w', newline='') as fout:
        w=csv.writer(fout)
        w.writerow(['histogram','dimension','nominal','special_scaled','combined','special_fraction','combined_minus_nominal'])
        for name,nom,spec,tot,frac,excess,dim in rows: w.writerow([name,dim,nom,spec,tot,frac,excess])
    print('\nWrote:', os.path.abspath(args.csv))
    print('\n=== MOST SIGNIFICANT SPECIAL FRACTION AMONG 1D BINS ===')
    candidates=[]
    for name,_,_,_,_,_,dim in rows:
        if dim!=1 or not is_analysis_hist(name): continue
        a,b,c=(hists[k][name] for k in ('nominal','special','combined'))
        for i in range(1,c.GetNbinsX()+1):
            total=c.GetBinContent(i)
            sp=scale*b.GetBinContent(i)
            if total>0 and sp>0:
                candidates.append((sp/total, name, bin_label(c,i), total, sp))
    candidates.sort(reverse=True)
    for frac,name,binrng,total,sp in candidates[:args.top]:
        print('  {:>7.1%} {:<67} bin {:<20} combined={:>10.4g} special={:>10.4g}'.format(frac,name[:67],binrng,total,sp))
    print('\nInterpretation: Closure checks combining arithmetic, NOT whether the special sample was produced/weighted correctly.')
    print('A high special fraction alone does NOT demonstrate a bug; inspect MC truth categories, weight logic, and data/MC as well.')
    for f in files.values(): f.Close()


if __name__=='__main__':
    main()

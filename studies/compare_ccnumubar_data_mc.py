#!/usr/bin/env python3
"""CC anti-numu LE5 data vs nominal + POT-scaled special 2p2h audit.

Run from an initialized MINERvA PyROOT environment. Read-only analysis.
"""
import argparse
import csv
import math
import os
import sys

import ROOT
ROOT.gROOT.SetBatch(True)

BASE = '/exp/minerva/data/users/qvuong/antinu_mu'
TAG = 'ValidateCuts_CCnumubar_Inclusive_Curvature_v1'
STEM = f'kin_dist_mcle5_{TAG}_MAD'
DEFAULTS = {
    'data': os.path.join(BASE, f'kin_dist_datale5_{TAG}_MAD.root'),
    'nominal': os.path.join(BASE, STEM + '_no2p2h.root'),
    'special': os.path.join(BASE, STEM + '_only2p2h.root'),
    'combined': os.path.join(BASE, STEM + '.root'),
}


def open_file(path):
    f = ROOT.TFile.Open(path, 'READ')
    if not f or f.IsZombie():
        raise RuntimeError('Cannot open ROOT file: ' + path)
    return f


def pot_from_file(path):
    # Match tools.Utilities.getPOTFromFile and the existing plotting workflow.
    from tools import Utilities
    pot = Utilities.getPOTFromFile(path)
    if pot is None or not math.isfinite(float(pot)) or pot <= 0:
        raise RuntimeError('Invalid POT for ' + path + ': ' + str(pot))
    return float(pot)


def get_hist(file, key):
    obj = file.Get(key)
    if not obj:
        raise RuntimeError('Missing histogram: ' + key + ' in ' + file.GetName())
    if not obj.InheritsFrom('TH1') or obj.GetDimension() != 1:
        raise RuntimeError(key + ' is not a 1D TH1')
    return obj


def matched_bins(hists):
    ref = hists[0]
    for h in hists[1:]:
        if h.GetNbinsX() != ref.GetNbinsX():
            raise RuntimeError('Different number of bins')
        for i in range(1, ref.GetNbinsX() + 2):
            if not math.isclose(ref.GetXaxis().GetBinLowEdge(i), h.GetXaxis().GetBinLowEdge(i), rel_tol=1e-9, abs_tol=1e-9):
                raise RuntimeError('Mismatched bin edges')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--key', default='Eel', help='Exact histogram key, e.g. Eel or Lepton_Ppar')
    ap.add_argument('--csv', default='audit_ccnumubar_Eel_data_mc.csv')
    ap.add_argument('--data', default=DEFAULTS['data'])
    ap.add_argument('--nominal', default=DEFAULTS['nominal'])
    ap.add_argument('--special', default=DEFAULTS['special'])
    ap.add_argument('--combined', default=DEFAULTS['combined'])
    ap.add_argument('--data-pot', type=float, help='Override data POT if required')
    ap.add_argument('--mc-pot', type=float, help='Override nominal/combined MC POT if required')
    ap.add_argument('--special-pot', type=float, help='Override special MC POT if required')
    args = ap.parse_args()

    paths = {k: getattr(args, k) for k in DEFAULTS}
    files = {k: open_file(path) for k, path in paths.items()}
    pots = {
        'data': args.data_pot if args.data_pot is not None else pot_from_file(paths['data']),
        'nominal': args.mc_pot if args.mc_pot is not None else pot_from_file(paths['nominal']),
        'special': args.special_pot if args.special_pot is not None else pot_from_file(paths['special']),
    }
    scale_special_to_nominal = pots['nominal'] / pots['special']
    scale_mc_to_data = pots['data'] / pots['nominal']
    print('POT data={:.12g} nominal={:.12g} special={:.12g}'.format(pots['data'], pots['nominal'], pots['special']))
    print('special -> nominal {:.12g}; nominal -> data {:.12g}'.format(scale_special_to_nominal, scale_mc_to_data))
    print('NOTE: no scaling by bin width; ratios are invariant for matched bins.\n')

    hs = {k: get_hist(f, args.key) for k, f in files.items()}
    matched_bins(list(hs.values()))
    rows = []
    totals = {k: 0.0 for k in ['data', 'nominal', 'special', 'combined']}
    max_closure = 0.0
    for i in range(1, hs['data'].GetNbinsX() + 1):
        lo = hs['data'].GetXaxis().GetBinLowEdge(i)
        hi = hs['data'].GetXaxis().GetBinUpEdge(i)
        d = hs['data'].GetBinContent(i)
        n = hs['nominal'].GetBinContent(i) * scale_mc_to_data
        s = hs['special'].GetBinContent(i) * scale_special_to_nominal * scale_mc_to_data
        m = hs['combined'].GetBinContent(i) * scale_mc_to_data
        closure = m - n - s
        max_closure = max(max_closure, abs(closure))
        for k, v in [('data', d), ('nominal', n), ('special', s), ('combined', m)]:
            totals[k] += v
        deficit = m - d
        rows.append({
            'bin': i, 'low': lo, 'high': hi, 'data': d,
            'nominal_without_2p2h': n, 'special_2p2h': s, 'total_mc': m,
            'data_over_total_mc': d / m if m > 0 else float('nan'),
            'data_over_nominal': d / n if n > 0 else float('nan'),
            'special_fraction_of_mc': s / m if m > 0 else float('nan'),
            'mc_minus_data': deficit,
            'excess_over_special': deficit / s if s > 0 else float('nan'),
            'closure_residual': closure,
        })

    print('Histogram:', args.key)
    print('All yields normalized to DATA POT; flow bins excluded.')
    print('{:>12} {:>12} {:>12} {:>12} {:>12} {:>10} {:>10} {:>10} {:>10}'.format(
        'bin', 'Data', 'Nominal', '2p2h', 'MC total', 'D/MC', 'D/Nom', '2p2h/MC', 'Excess/2p2h'))
    for r in rows:
        print('{low:5.3g}-{high:<6.3g} {data:12.1f} {nominal_without_2p2h:12.1f} {special_2p2h:12.1f} '
              '{total_mc:12.1f} {data_over_total_mc:10.3f} {data_over_nominal:10.3f} '
              '{special_fraction_of_mc:10.1%} {excess_over_special:10.2f}'.format(**r))

    print('\nINTEGRATED (within histogram x range)')
    print('Data={data:.3f} Nominal={nominal:.3f} Special={special:.3f} MC={combined:.3f}'.format(**totals))
    if totals['combined'] > 0:
        print('D/MC={:.4f}; special/MC={:.2%}'.format(totals['data']/totals['combined'], totals['special']/totals['combined']))
    if totals['special'] > 0:
        print('(MC - Data)/special={:.4f}'.format((totals['combined']-totals['data'])/totals['special']))
    print('Max |MC - nominal - special| after scaling={:.6g}'.format(max_closure))
    print('\nInterpretation: positive Excess/2p2h = MC > data; >1 means removing all 2p2h would not close the gap.')
    print('Nominal excludes 2p2h as previously verified; this script does not adjust model or fit data.')
    outdir = os.path.dirname(os.path.abspath(args.csv))
    os.makedirs(outdir, exist_ok=True)
    with open(args.csv, 'w', newline='') as fout:
        w = csv.DictWriter(fout, fieldnames=list(rows[0]) if rows else ['bin'])
        w.writeheader()
        w.writerows(rows)
    print('CSV:', args.csv)
    for f in files.values():
        f.Close()


if __name__ == '__main__':
    main()

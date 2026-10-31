#!/usr/bin/env python3
"""Compare MAT GetFluxAndCVWeight values in an audit CSV with installed LE5 flux ratios.

Read-only. Run in the GPVM environment with PyROOT. No pandas needed.
"""
import argparse
import csv
import math
import os
from collections import defaultdict

import ROOT
ROOT.gROOT.SetBatch(True)

DEFAULT_DIR = 'custom_plotutils/data/flux'


# Keep ROOT files alive while their histograms are used. Four small files only.
# Avoid Clone()/SetDirectory()/Close() lifetime issues with PyROOT/cppyy.
_OPEN_FILES = []


def load_hist(path, name):
    f = ROOT.TFile.Open(path, 'READ')
    if not f or f.IsZombie():
        raise RuntimeError('Cannot open {}'.format(path))
    _OPEN_FILES.append(f)
    h = f.Get(name)
    if not h or not h.InheritsFrom('TH1'):
        raise RuntimeError('Missing histogram {} in {}'.format(name, path))
    return h


def get_pair(directory, pdg):
    tag = 'pdg{}-minerva5.root'.format(pdg)
    cv = load_hist(os.path.join(directory, 'flux-gen2thin-' + tag), 'flux_E_cvweighted')
    un = load_hist(os.path.join(directory, 'flux-g4numiv5-' + tag), 'flux_E_unweighted')
    if cv.GetNbinsX() != un.GetNbinsX():
        raise RuntimeError('Binning size mismatch for PDG {}'.format(pdg))
    for b in range(1, cv.GetNbinsX() + 2):
        a, c = cv.GetXaxis().GetBinLowEdge(b), un.GetXaxis().GetBinLowEdge(b)
        if not math.isclose(a, c, rel_tol=1e-12, abs_tol=1e-10):
            raise RuntimeError('Bin edges differ for PDG {}'.format(pdg))
    return cv, un


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--audit', default='rhc_cv_weight_audit.csv')
    ap.add_argument('--flux-dir', default=DEFAULT_DIR)
    ap.add_argument('--out', default='flux_weight_closure.csv')
    ap.add_argument('--print-events', type=int, default=12)
    args = ap.parse_args()
    if not os.path.isfile(args.audit):
        ap.error('Audit CSV not found: {}'.format(args.audit))
    pairs = {pdg: get_pair(args.flux_dir, pdg) for pdg in [-14, 14]}
    fields = ['entry', 'pdg', 'Etrue_GeV', 'root_bin', 'mat_flux_weight',
              'hist_cv', 'hist_unweighted', 'hist_ratio', 'mat_over_hist_ratio',
              'difference', 'status']
    stats = defaultdict(list)
    counts = defaultdict(int)
    printed = 0
    with open(args.audit, newline='') as src, open(args.out, 'w', newline='') as dst:
        writer = csv.DictWriter(dst, fieldnames=fields)
        writer.writeheader()
        for event in csv.DictReader(src):
            pdg = int(event['pdg'])
            if pdg not in pairs:
                counts['other_species'] += 1
                continue
            e = float(event['Etrue_GeV'])
            w = float(event['flux'])
            hc, hu = pairs[pdg]
            axis = hc.GetXaxis()
            ibin = axis.FindFixBin(e)
            cv = float(hc.GetBinContent(ibin))
            un = float(hu.GetBinContent(ibin))
            in_range = 1 <= ibin <= hc.GetNbinsX()
            valid = in_range and math.isfinite(un) and un > 0 and math.isfinite(cv) and math.isfinite(w)
            ratio = cv / un if valid else float('nan')
            scale = w / ratio if valid and ratio != 0 else float('nan')
            diff = w - ratio if valid else float('nan')
            status = 'OK' if valid and ratio > 0 else ('OUT_OF_RANGE' if not in_range else 'INVALID_RATIO')
            counts[(pdg, status)] += 1
            if status == 'OK':
                stats[pdg].append((e, w, ratio, scale, diff))
            writer.writerow(dict(entry=event['entry'], pdg=pdg, Etrue_GeV=e,
                                 root_bin=ibin, mat_flux_weight=w, hist_cv=cv,
                                 hist_unweighted=un, hist_ratio=ratio,
                                 mat_over_hist_ratio=scale, difference=diff, status=status))
            if printed < args.print_events:
                print('entry {:>6} PDG {:>3} Enu={:7.3f} GeV MAT={:9.5f} '
                      'CV/un={:9.5f} MAT/ratio={:9.5f} {}'.format(
                          event['entry'], pdg, e, w, ratio, scale, status))
                printed += 1
    print('\nComparison by species:')
    for pdg in [-14, 14]:
        rows = stats[pdg]
        print('\nPDG {}: {} valid rows'.format(pdg, len(rows)))
        if not rows:
            continue
        for name, idx in [('MAT weight', 1), ('CV/unweighted', 2),
                          ('MAT divided by ratio', 3), ('MAT minus ratio', 4)]:
            vals = [r[idx] for r in rows]
            print('  {:22s} mean={:.6g}  min={:.6g}  max={:.6g}'.format(
                  name, sum(vals)/len(vals), min(vals), max(vals)))
        # Compare binned averages in true neutrino energy, not reco muon energy.
        for lo, hi in [(0, 2), (2, 3), (3, 5), (5, 8), (8, 12), (12, 20), (20, 100)]:
            subset = [r for r in rows if lo <= r[0] < hi]
            if subset:
                wm = sum(r[1] for r in subset) / len(subset)
                wr = sum(r[2] for r in subset) / len(subset)
                print('  Enu [{:>2},{:>3}) n={:4} MAT={:8.4f} ratio={:8.4f} delta={:+8.4f}'.format(
                      lo, hi, len(subset), wm, wr, wm-wr))
    print('\nSkipped/flagged:', dict(counts))
    print('Wrote:', os.path.abspath(args.out))
    print('\nNote: Matching is a hypothesis, not a guaranteed MAT convention. '
          'Differences can reflect alternative flux inputs or interpolation/normalization.')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Compare LE/ME MHRW standard-selection ROOT outputs (CV bin contents)."""
import argparse
import csv
import math
from pathlib import Path
import ROOT

ROOT.gROOT.SetBatch(True)

DEFAULT_DIR = Path('/exp/minerva/data/users/qvuong/nu_mu')

def read_hists(path):
    f = ROOT.TFile.Open(str(path), 'READ')
    if not f or f.IsZombie():
        raise RuntimeError(f'Cannot open {path}')
    out = {}
    def walk(directory, prefix=''):
        for key in directory.GetListOfKeys():
            name = key.GetName()
            obj = key.ReadObj()
            full = f'{prefix}{name}'
            if obj.InheritsFrom('TDirectory'):
                walk(obj, full + '/')
            elif obj.InheritsFrom('TH1'):
                # The central histogram of an MnvH1D/H2D is accessible through GetBinContent.
                # Do NOT sum error-band universes or include them as independent hists.
                d = obj.GetDimension()
                if d > 2: continue
                content = [float(obj.GetBinContent(i)) for i in range(obj.GetNcells())]
                out[full] = (d, content, obj.ClassName())
    walk(f)
    f.Close()
    return out

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--le', type=Path, default=DEFAULT_DIR/'kin_disttest_mcle1_p8_MHRW_LE_MAD_0.root')
    ap.add_argument('--me', type=Path, default=DEFAULT_DIR/'kin_disttest_mcle1_p8_MHRW_ME_MAD_0.root')
    ap.add_argument('--csv', type=Path, default=Path('studies/mhrw_selection_comparison.csv'))
    ap.add_argument('--top', type=int, default=25)
    args = ap.parse_args()
    le, me = read_hists(args.le), read_hists(args.me)
    common = sorted(le.keys() & me.keys())
    print(f'LE histograms: {len(le)} | ME histograms: {len(me)} | common: {len(common)}')
    if not common:
        print('LE names:', list(le)[:20]); print('ME names:', list(me)[:20]); return
    rows=[]
    for name in common:
        dl,vl,cl = le[name]; dm,vm,cm = me[name]
        if dl != dm or len(vl) != len(vm):
            print('SKIPPING mismatch:', name); continue
        # Cell 0 is underflow (in 1D), so compute all-cell sum and nonzero bins.
        # For TH2 the ROOT global-cell layout includes flow bins as well.
        sl,sm=sum(vl),sum(vm)
        delta = 100*(sm/sl-1) if sl else float('nan')
        sum_abs = sum(abs(a-b) for a,b in zip(vl,vm))
        nonzero = sum(abs(a-b)>1e-10 for a,b in zip(vl,vm))
        max_rel = max((abs(b/a-1)*100 for a,b in zip(vl,vm) if a != 0),default=float('nan'))
        rows.append({'histogram':name,'class':cl,'dimension':dl,'LE_sum_all_cells':sl,'ME_sum_all_cells':sm,'ME_over_LE_minus_1_pct':delta,'sum_abs_bin_diff':sum_abs,'n_changed_bins':nonzero,'max_bin_pct_for_nonzero_LE':max_rel})
    args.csv.parent.mkdir(parents=True,exist_ok=True)
    with args.csv.open('w',newline='') as fh:
        writer=csv.DictWriter(fh,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
    changed=[r for r in rows if r['n_changed_bins']]
    print(f'Histograms with changed CV bin contents: {len(changed)} / {len(rows)}')
    print('\nLargest fractional changes in summed CV contents (include flow bins):')
    rank=sorted((r for r in changed if math.isfinite(r['ME_over_LE_minus_1_pct'])),key=lambda r:abs(r['ME_over_LE_minus_1_pct']),reverse=True)
    for r in rank[:args.top]:
        print(f" {r['ME_over_LE_minus_1_pct']:+9.4f}% | LE={r['LE_sum_all_cells']:.8g} ME={r['ME_sum_all_cells']:.8g} | changed bins={r['n_changed_bins']:3d} | {r['histogram']}")
    if not changed: print('No changed central histogram bins were found.')
    print(f'\nFull comparison: {args.csv}')
    print('Note: this compares saved central histograms, not systematic universes or errors.')

if __name__ == '__main__': main()

#!/usr/bin/env python3
"""Read-only inventory and POT audit of MINERvA FHC CC-nue v4.3.1 files.

Run inside an FNAL environment with PyROOT available. The audit deliberately does
not guess histogram names or claim that energy-binned dE/dx shapes exist.
"""
import argparse
import csv
import os
import re
from pathlib import Path

DEFAULT_DATA = '/exp/minerva/data/users/qvuong/nu_e/kin_dist_dataleFHC_CCnue_FullCuts_v4p3p1_MAD.root'
DEFAULT_MC = '/exp/minerva/data/users/qvuong/nu_e/kin_dist_mcleFHC_CCnue_FullCuts_v4p3p1_MAD.root'
DEFAULT_TUNED = '/exp/minerva/data/users/qvuong/nu_e/bkgfit_leFHC_N4_tune_CCnue_FullCuts_v4p3p1_MAD.root'


def metadata_pot(fin):
    meta = fin.Get('Meta')
    if not meta or not meta.InheritsFrom('TTree'):
        return None, 'Meta TTree missing'
    branches = {b.GetName() for b in meta.GetListOfBranches()}
    for key in ('POT_Used', 'POTUsed', 'used_POT', 'UsedPOT', 'POT', 'pot', 'POT_Total', 'total_pot'):
        if key in branches:
            val = 0.0
            for i in range(meta.GetEntries()):
                meta.GetEntry(i)
                val += float(getattr(meta, key))
            return val, key
    return None, 'No recognized POT branch; available: ' + ', '.join(sorted(branches))


def walk(directory, prefix=''):
    for k in directory.GetListOfKeys():
        name = k.GetName()
        path = f'{prefix}/{name}' if prefix else name
        obj = k.ReadObj()
        if obj.InheritsFrom('TDirectory'):
            yield from walk(obj, path)
        elif obj.InheritsFrom('TH1'):
            yield path, obj


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data', default=DEFAULT_DATA)
    p.add_argument('--mc', default=DEFAULT_MC)
    p.add_argument('--tuned', default=DEFAULT_TUNED)
    p.add_argument('--out', default='fhc_ccnue_dedx_audit')
    p.add_argument('--dedx-cut', type=float, default=2.64)
    args = p.parse_args()

    import ROOT
    ROOT.gROOT.SetBatch(True)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    sources = {'data': args.data, 'mc': args.mc, 'tuned': args.tuned}
    handles = {}
    pots = {}
    rows = []
    summary = []
    for label, path in sources.items():
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            raise RuntimeError(f'Cannot open {label} ROOT file: {path}')
        handles[label] = f
        pot, branch = metadata_pot(f)
        pots[label] = pot
        print(f'{label}: {path}\n  POT={pot}, source={branch}')
        summary.append({'file': label, 'path': path, 'pot': pot, 'pot_source': branch})

        for objpath, h in walk(f):
            is_2d = h.InheritsFrom('TH2')
            dims = 3 if h.InheritsFrom('TH3') else (2 if is_2d else 1)
            title = str(h.GetTitle())
            s = (objpath + ' ' + title).lower()
            relevant = bool(re.search(r'dedx|de/dx|d[e]?dx|energy|en4|sideband|signal|eavail|lepton', s))
            rows.append({
                'file': label, 'hist_path': objpath, 'class': h.ClassName(),
                'dimensions': dims, 'nbins_x': h.GetNbinsX(),
                'nbins_y': h.GetNbinsY() if dims >= 2 else '',
                'title': title, 'relevant': int(relevant),
                'integral': float(h.Integral()),
            })
    if not pots['data'] or not pots['mc']:
        print('WARNING: Cannot determine data/MC POT ratio; no scaled comparison will be made.')
    else:
        print(f'\nMC-to-data POT scale: {pots["data"] / pots["mc"]:.12g}')

    with (out / 'histogram_inventory.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]) if rows else
                           ['file','hist_path','class','dimensions','nbins_x','nbins_y','title','relevant','integral'])
        w.writeheader()
        w.writerows(rows)
    with (out / 'pot_summary.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['file','path','pot','pot_source'])
        w.writeheader()
        w.writerows(summary)

    # Compare histograms only when data and MC have IDENTICAL object paths and bin edges.
    data_hists = dict(walk(handles['data']))
    mc_hists = dict(walk(handles['mc']))
    yields = []
    scale = pots['data']/pots['mc'] if pots['data'] and pots['mc'] else None
    for path in sorted(set(data_hists) & set(mc_hists)):
        hd, hm = data_hists[path], mc_hists[path]
        if hd.InheritsFrom('TH2') or hm.InheritsFrom('TH2'):
            continue
        if hd.GetNbinsX() != hm.GetNbinsX():
            continue
        n = hd.GetNbinsX()
        if any(abs(hd.GetXaxis().GetBinLowEdge(i)-hm.GetXaxis().GetBinLowEdge(i)) > 1e-8 for i in range(1,n+2)):
            continue
        name = (path + ' ' + hd.GetTitle()).lower()
        if not any(k in name for k in ('energy', 'en4', 'sideband', 'dedx', 'de/dx')):
            continue
        for lo, hi in ((0.8,2.0), (2.0,2.5)):
            # Only integrate whole histogram bins entirely inside the interval.
            bins = [i for i in range(1,n+1) if hd.GetXaxis().GetBinLowEdge(i) >= lo-1e-8 and hd.GetXaxis().GetBinUpEdge(i) <= hi+1e-8]
            if not bins:
                continue
            yields.append({'hist_path': path, 'energy_lo': lo, 'energy_hi': hi,
                           'num_bins': len(bins), 'data_integral': sum(hd.GetBinContent(i) for i in bins),
                           'mc_raw_integral': sum(hm.GetBinContent(i) for i in bins),
                           'mc_at_data_pot': sum(hm.GetBinContent(i) for i in bins)*scale if scale is not None else '',
                           'note': 'Only valid as energy-binned yields if X axis is energy; inspect inventory/title.'})
    with (out / 'candidate_energy_yields.csv').open('w', newline='') as f:
        cols = ['hist_path','energy_lo','energy_hi','num_bins','data_integral','mc_raw_integral','mc_at_data_pot','note']
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(yields)

    print(f'\nInventory: {out / "histogram_inventory.csv"} ({len(rows)} histograms)')
    print(f'Candidate energy comparisons: {out / "candidate_energy_yields.csv"} ({len(yields)} rows)')
    print('NOTE: Candidate comparisons require axis/name validation. Do not interpret dE/dx-sideband yields from energy-only histograms.')
    print(f'Front dE/dx cut for reference: {args.dedx_cut} MeV/cm')
    for f in handles.values():
        f.Close()


if __name__ == '__main__':
    main()

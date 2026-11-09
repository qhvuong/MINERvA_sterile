#!/usr/bin/env python3
"""Read-only CV diagnostic for MINERvA LE CC newBinning selections.

Requires PyROOT on the GPVM. No systematic-universe or background-fit inputs.
Outputs: inventory.csv, pot_summary.csv, cv_bin_yields.csv,
         mc_components.csv, region_summary.csv, binning_comparison.csv,
         ratio_diagnostics.csv, diagnostic_plots.pdf, summary.txt.

The displayed CC ratios are *raw selected-event* diagnostics, NOT the
background-subtracted ratios in the stitched oscillation fit.
"""
import argparse
import csv
import math
import os
from collections import defaultdict
from pathlib import Path

import ROOT
ROOT.gROOT.SetBatch(True)

BASE = Path('/exp/minerva/data/users/qvuong')
SAMPLES = {
    'fhc_ccnue':   ('nu_e',     'leFHC', 'CCnue'),
    'fhc_ccnumu':  ('nu_mu',    'leFHC', 'CCnumu'),
    'rhc_ccnuebar':('antinu_e', 'le5',   'CCnuebar'),
    'rhc_ccnumubar':('antinu_mu','le5', 'CCnumubar'),
}
COMPONENTS_SIGNAL = {'CCNuEQE', 'CCNuEDelta', 'CCNuEDIS',
                     'CCNuE2p2h', 'CCNuE', 'CCNuEWrongSign'}
NOMINAL_CC_EDGES = [0.8, 2., 2.5, 3., 3.8, 4.6, 5.5, 6.5,
                    8., 10., 12.5, 16., 20.]


def output_csv(path, rows, fields):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def read_pot(f):
    t = f.Get('Meta')
    if not t:
        return None, 'no Meta tree'
    branches = {b.GetName() for b in t.GetListOfBranches()}
    for candidate in ('POT_Used','POTUsed','used_POT','UsedPOT','POT','pot','POT_Total','total_pot'):
        if candidate in branches:
            total = 0.
            for entry in t:
                total += float(getattr(entry, candidate))
            return total, candidate
    return None, 'no recognized POT branch: ' + ','.join(sorted(branches))


def is_hist(obj):
    return bool(obj) and obj.InheritsFrom('TH1') and obj.GetDimension() == 1


def get_hist(f, key):
    h = f.Get(key)
    return h if is_hist(h) else None


def axis_edges(h):
    ax = h.GetXaxis()
    return [float(ax.GetBinLowEdge(i)) for i in range(1, h.GetNbinsX()+1)] + [float(ax.GetBinUpEdge(h.GetNbinsX()))]


def aligned(a, b):
    return len(a) == len(b) and all(math.isclose(x, y, rel_tol=1e-8, abs_tol=1e-7) for x,y in zip(a,b))


def integral_on_range(h, lo, hi):
    """Exact integration only when requested boundaries coincide with existing edges."""
    if h is None:
        return None
    e = axis_edges(h)
    matches_lo = [i for i,x in enumerate(e) if math.isclose(x,lo,abs_tol=1e-6,rel_tol=1e-8)]
    matches_hi = [i for i,x in enumerate(e) if math.isclose(x,hi,abs_tol=1e-6,rel_tol=1e-8)]
    if not matches_lo or not matches_hi or matches_hi[0] <= matches_lo[0]:
        return None
    return sum(float(h.GetBinContent(i)) for i in range(matches_lo[0]+1, matches_hi[0]+1))


def inventory(f, sample, kind):
    rows=[]
    for key in f.GetListOfKeys():
        name=key.GetName()
        obj=key.ReadObj()
        if not obj.InheritsFrom('TH1'):
            continue
        rows.append(dict(sample=sample,kind=kind,histogram=name,
                         dimension=int(obj.GetDimension()), bins_x=int(obj.GetNbinsX()),
                         xmin=float(obj.GetXaxis().GetXmin()), xmax=float(obj.GetXaxis().GetXmax())))
    return rows


def components(f, region):
    """Restrict to direct KEY entries named region_<category>, avoiding ROOT nested paths."""
    names={k.GetName() for k in f.GetListOfKeys()}
    return {name[len(region)+1:]:get_hist(f,name)
            for name in sorted(names) if name.startswith(region+'_') and
            get_hist(f,name) is not None}


def compatible_pairs(files):
    return [(x,y) for x,y in [('fhc_ccnumu','fhc_ccnue'),
                              ('rhc_ccnumubar','rhc_ccnuebar')]
            if x in files and y in files]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base', type=Path, default=BASE)
    ap.add_argument('--out', type=Path, default=Path('cv_new_binning_diagnostics'))
    ap.add_argument('--tag', default='ValidateSideband')
    ap.add_argument('--suppress-plots', action='store_true')
    args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    opened={}
    pot_rows=[]; inventory_rows=[]; bins=[]; components_rows=[]; summaries=[]; comparison=[]; ratios=[]
    notes=[]
    for sample,(folder,playlist,selection) in SAMPLES.items():
        suffix=f'{playlist}_{args.tag}_{selection}_newBinning_MAD.root'
        data_path=args.base/folder/f'kin_dist_data{suffix}'
        mc_path=args.base/folder/f'kin_dist_mc{suffix}'
        if not data_path.is_file() or not mc_path.is_file():
            notes.append(f'MISSING {sample}: data={data_path.exists()}, mc={mc_path.exists()}')
            continue
        df=ROOT.TFile.Open(str(data_path))
        mf=ROOT.TFile.Open(str(mc_path))
        if not df or df.IsZombie() or not mf or mf.IsZombie():
            notes.append(f'FAILED TO OPEN {sample}')
            continue
        pd,sd=read_pot(df); pm,sm=read_pot(mf)
        if pd is None or pm is None or pd<=0 or pm<=0:
            notes.append(f'INVALID POT {sample}: data={pd} ({sd}), mc={pm} ({sm}); no yields computed')
            df.Close();mf.Close();continue
        scale=pd/pm
        pot_rows.append(dict(sample=sample,data_file=str(data_path),mc_file=str(mc_path),
                             data_pot=pd,mc_pot=pm,mc_to_data_scale=scale,
                             data_pot_branch=sd,mc_pot_branch=sm))
        inventory_rows.extend(inventory(df,sample,'data'))
        inventory_rows.extend(inventory(mf,sample,'mc'))
        opened[sample]=(df,mf,scale)
        regions=['EN4']
        if 'nue' in sample:
            regions+=['EN4_dEdX']
        else:
            if get_hist(df,'EN4_dEdX') and get_hist(mf,'EN4_dEdX'):
                regions+=['EN4_dEdX']
        for region in regions:
            hd=get_hist(df,region); hm=get_hist(mf,region)
            if hd is None or hm is None:
                notes.append(f'NO PAIRED {region} histograms for {sample}')
                continue
            e=axis_edges(hd)
            if not aligned(e,axis_edges(hm)):
                notes.append(f'MISMATCHED data/MC edges for {sample} {region}: {e} / {axis_edges(hm)}')
                continue
            parts=components(mf,region)
            parts={k:v for k,v in parts.items() if aligned(e,axis_edges(v))}
            n=hd.GetNbinsX()
            total_data=sum(float(hd.GetBinContent(i)) for i in range(1,n+1))
            total_mc=sum(float(hm.GetBinContent(i))*scale for i in range(1,n+1))
            summaries.append(dict(sample=sample,region=region,nbins=n,edges=';'.join(f'{x:g}' for x in e),
                                  data=total_data,mc=total_mc,data_over_mc=total_data/total_mc if total_mc else '',
                                  diff=total_data-total_mc,scale=scale))
            for i in range(1,n+1):
                lo,hi=e[i-1:i+1]
                d=float(hd.GetBinContent(i));m=float(hm.GetBinContent(i))*scale
                stat_err=float(hd.GetBinError(i))
                compvals={k:float(v.GetBinContent(i))*scale for k,v in parts.items()}
                comp_sum=sum(compvals.values())
                sig=sum(v for k,v in compvals.items() if k in COMPONENTS_SIGNAL)
                pi0=sum(v for k,v in compvals.items() if 'Pi0' in k)
                row=dict(sample=sample,region=region,bin=i,lo=lo,hi=hi,data=d,mc=m,
                         diff=d-m,ratio=d/m if m>0 else '',stat_err=stat_err,
                         raw_pull=(d-m)/stat_err if stat_err>0 else '',
                         cc_nue_signal_mc=sig if parts else '',
                         pi0_mc=pi0 if parts else '',
                         signal_purity=sig/m if m>0 and parts else '',
                         pi0_fraction=pi0/m if m>0 and parts else '',
                         categories_sum=comp_sum if parts else '',
                         sum_minus_total=comp_sum-m if parts else '',
                         sample_key=f'{sample}:{region}:{lo:.5f}:{hi:.5f}')
                bins.append(row)
                if parts and abs(comp_sum-m)>1e-3*max(1,abs(m)):
                    notes.append(f'COMPONENT SUM discrepancy {sample} {region} bin {i}: '+
                                 f'components={comp_sum:.3f}, total={m:.3f}; inspect categories')
                for cat,val in compvals.items():
                    components_rows.append(dict(sample=sample,region=region,bin=i,lo=lo,hi=hi,
                                                category=cat,mc_scaled=val,fraction=val/m if m else ''))
            # Compare only windows with EXACT existing source boundaries.
            for lo,hi in [(0.8,1.),(1.,1.5),(1.,2.),(1.5,2.),(2.,2.5),
                          (0.8,2.),(0.8,2.5),(1.,2.5)]:
                d=integral_on_range(hd,lo,hi)
                m=integral_on_range(hm,lo,hi)
                if d is None or m is None:
                    comparison.append(dict(sample=sample,region=region,lo=lo,hi=hi,
                                           status='NOT EXACTLY AVAILABLE (no fractional-bin approximation)'))
                else:
                    m*=scale
                    comparison.append(dict(sample=sample,region=region,lo=lo,hi=hi,
                                           status='OK',data=d,mc=m,data_over_mc=d/m if m else '',diff=d-m))
        print(f'{sample}: POT scale={scale:.9g}; main EN4 bins='+
              str(getattr(get_hist(df,'EN4'),'GetNbinsX',lambda:'MISSING')()))

    # Raw selected-event ratios should be interpreted only as diagnostics.
    for num,den in compatible_pairs(opened):
        dn,mn,sn=opened[num]
        dd,md,sd=opened[den]
        num_data=get_hist(dn,'EN4'); den_data=get_hist(dd,'EN4')
        num_mc=get_hist(mn,'EN4'); den_mc=get_hist(md,'EN4')
        if not all([num_data,den_data,num_mc,den_mc]):
            notes.append(f'NO EN4 for diagnostic ratio {num}/{den}')
            continue
        e=axis_edges(num_data)
        if not all(aligned(e,axis_edges(h)) for h in [den_data,num_mc,den_mc]):
            notes.append(f'MISMATCHED EN4 bin edges; cannot calculate {num}/{den} selected-event ratios')
            continue
        # Even if equal binning, the FHC/RHC within-mode POT exposures must also match.
        npot=next(p['data_pot'] for p in pot_rows if p['sample']==num)
        dpot=next(p['data_pot'] for p in pot_rows if p['sample']==den)
        if not math.isclose(npot,dpot,rel_tol=5e-3):
            notes.append(f'DIFFERENT DATA POT {num}/{den}: {npot} vs {dpot}; '+
                         'raw ratios are NOT directly exposure-matched')
        for i in range(1,num_data.GetNbinsX()+1):
            nd=float(num_data.GetBinContent(i));ed=float(den_data.GetBinContent(i))
            nm=float(num_mc.GetBinContent(i))*sn;em=float(den_mc.GetBinContent(i))*sd
            ratios.append(dict(mode=num.split('_')[0],bin=i,lo=e[i-1],hi=e[i],
                               num_data=nd,den_data=ed,num_mc=nm,den_mc=em,
                               raw_data_ratio=nd/ed if ed else '',raw_mc_ratio=nm/em if em else '',
                               ratio_data_over_mc=(nd/ed)/(nm/em) if ed and em and nm else '',
                               num_data_pot=npot,den_data_pot=dpot,
                               note='RAW SELECTED EVENT RATIO ONLY; not background-subtracted'))

    output_csv(args.out/'pot_summary.csv',pot_rows,
               ['sample','data_file','mc_file','data_pot','mc_pot','mc_to_data_scale','data_pot_branch','mc_pot_branch'])
    output_csv(args.out/'histogram_inventory.csv',inventory_rows,
               ['sample','kind','histogram','dimension','bins_x','xmin','xmax'])
    output_csv(args.out/'cv_bin_yields.csv',bins,
               ['sample','region','bin','lo','hi','data','mc','diff','ratio','stat_err','raw_pull',
                'cc_nue_signal_mc','pi0_mc','signal_purity','pi0_fraction','categories_sum',
                'sum_minus_total','sample_key'])
    output_csv(args.out/'mc_components.csv',components_rows,
               ['sample','region','bin','lo','hi','category','mc_scaled','fraction'])
    output_csv(args.out/'region_summary.csv',summaries,
               ['sample','region','nbins','edges','data','mc','data_over_mc','diff','scale'])
    output_csv(args.out/'binning_comparison.csv',comparison,
               ['sample','region','lo','hi','status','data','mc','data_over_mc','diff'])
    output_csv(args.out/'ratio_diagnostics.csv',ratios,
               ['mode','bin','lo','hi','num_data','den_data','num_mc','den_mc',
                'raw_data_ratio','raw_mc_ratio','ratio_data_over_mc',
                'num_data_pot','den_data_pot','note'])

    with open(args.out/'summary.txt','w') as f:
        f.write('MINERvA LE newBinning CV audit\n')
        f.write('Raw selected-event ratio is NOT background-subtracted.\n')
        f.write('raw_pull uses data statistical errors ONLY (NOT a full-systematics significance).\n')
        f.write('Bin-window integrations are made ONLY for exactly matching energy edges.\n\n')
        for r in summaries:
            f.write(f"{r['sample']:17s} {r['region']:10s}: data={r['data']:.2f}, "
                    f"MC={r['mc']:.2f}, D/MC={r['data_over_mc']:.3f}, "
                    f"nbins={r['nbins']}, edges={r['edges']}\n")
        f.write('\nWarnings/notes:\n')
        for note in notes:f.write(' - '+note+'\n')

    if not args.suppress_plots:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_pdf import PdfPages
            with PdfPages(args.out/'diagnostic_plots.pdf') as pdf:
                for sample in SAMPLES:
                    for region in ['EN4','EN4_dEdX']:
                        subset=[r for r in bins if r['sample']==sample and r['region']==region]
                        if not subset:continue
                        x=[(r['lo']+r['hi'])/2 for r in subset]
                        widths=[(r['hi']-r['lo']) for r in subset]
                        fig,ax=plt.subplots(figsize=(10,5))
                        ax.axhline(1,linestyle='--',color='gray',linewidth=1)
                        ax.errorbar(x,[r['ratio'] if r['ratio'] != '' else float('nan') for r in subset],
                                    xerr=[w/2 for w in widths],
                                    yerr=[r['stat_err']/r['mc'] if r['mc'] else 0 for r in subset],
                                    fmt='o',capsize=2)
                        ax.set_xlabel(r'Reconstructed $E_\nu$ (GeV)')
                        ax.set_ylabel('Data / prefit MC (statistical bars only)')
                        ax.set_title(f'{sample}  {region}   CV newBinning')
                        ax.set_xlim(min(r['lo'] for r in subset), max(r['hi'] for r in subset))
                        ax.grid(alpha=.25)
                        fig.tight_layout();pdf.savefig(fig);plt.close(fig)
                for mode in ['fhc','rhc']:
                    subset=[r for r in ratios if r['mode']==mode]
                    if not subset:continue
                    fig,ax=plt.subplots(figsize=(10,5))
                    x=[(r['lo']+r['hi'])/2 for r in subset]
                    ax.axhline(1,linestyle='--',color='gray',linewidth=1)
                    ax.plot(x,[r['ratio_data_over_mc'] if r['ratio_data_over_mc']!='' else float('nan') for r in subset],marker='o')
                    ax.set_xlabel(r'Reconstructed $E_\nu$ (GeV)')
                    ax.set_ylabel('(raw data mu/e) / (raw MC mu/e)')
                    ax.set_title(f'{mode.upper()} raw selected-event ratio check (NOT oscillation-fit ratio)')
                    ax.grid(alpha=.25);fig.tight_layout();pdf.savefig(fig);plt.close(fig)
        except ImportError:
            notes.append('Matplotlib missing; no PDF created. CSVs and summary are valid.')
            print('Matplotlib unavailable: skipped PDF plots.')
    for df,mf,_ in opened.values():df.Close();mf.Close()
    print('Wrote diagnostics to',args.out)
    print('WARNING: CC raw ratios use total selected events, not background-subtracted CCnue.')
    for n in notes:print('NOTE:',n)


if __name__=='__main__':
    main()

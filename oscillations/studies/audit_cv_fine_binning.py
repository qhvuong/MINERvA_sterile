#!/usr/bin/env python3
"""Read-only MINERvA LE fine-binning CV audit for four CC selections.

Scientific intent: diagnose the data/MC residuals, NOT optimize bins, cut
away discrepancies, calculate sterile significance, or decide which threshold
is preferable from observed chi2. Compare 0.8 and 1.0 GeV thresholds from
one common fine-binned CV production. The raw mu/e ratios are NOT the
background-subtracted ratios entering the sterile oscillation fit.

Requires PyROOT; matplotlib optional. No systematic universes or bkgfit.
"""
import argparse
import csv
import math
from pathlib import Path

import ROOT
ROOT.gROOT.SetBatch(True)

BASE = Path('/exp/minerva/data/users/qvuong')
SAMPLES = {
    'fhc_ccnue': ('nu_e', 'leFHC', 'CCnue'),
    'fhc_ccnumu': ('nu_mu', 'leFHC', 'CCnumu'),
    'rhc_ccnuebar': ('antinu_e', 'le5', 'CCnuebar'),
    'rhc_ccnumubar': ('antinu_mu', 'le5', 'CCnumubar'),
}
FINE_EDGES = [0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.8, 4.6, 5.5,
              6.5, 8.0, 10.0, 12.5, 16.0, 20.0]
OLD_EDGES = [0.8, 2.0, 2.5, 3.0, 3.8, 4.6, 5.5,
             6.5, 8.0, 10.0, 12.5, 16.0, 20.0]
TEST_EDGES = [1.0, 2.0, 2.5, 3.0, 3.8, 4.6, 5.5,
              6.5, 8.0, 10.0, 12.5, 16.0, 20.0]
LOW_WINDOWS = [(0.8, 1.0), (1.0, 1.5), (1.5, 2.0),
               (2.0, 2.5), (0.8, 2.0), (1.0, 2.0),
               (0.8, 2.5), (1.0, 2.5)]
# Only *exclusive, top-level* category names. The PDG subcategories, sideband
# subcategories, and histogram TOTALS must never be included in a sum.
E_CATEGORIES = {
    'CCNuEQE', 'CCNuEDelta', 'CCNuEDIS', 'CCNuE2p2h', 'CCNuE',
    'CCNuEWrongSign', 'NCPi0', 'NCCohPi0', 'CCPi0', 'NCPi', 'CCPi',
    'CCOther', 'NCOther', 'NCDiff', 'NuEElastic', 'Other',
}
MU_CATEGORIES = {
    'CCNuMu', 'CCNuMu2p2h', 'CCNuMuDIS', 'CCNuMuDelta',
    'CCNuMuQE', 'CCNuMuWrongSign', 'CC', 'NC', 'NCDiff',
    'NuEElastic', 'Other',
}
E_SIGNAL_CATEGORIES = {
    'CCNuEQE', 'CCNuEDelta', 'CCNuEDIS', 'CCNuE2p2h',
    'CCNuE', 'CCNuEWrongSign',
}
PI0_CATEGORIES = {'NCPi0', 'NCCohPi0', 'CCPi0'}


def near(x, y):
    return math.isclose(float(x), float(y), rel_tol=1e-9, abs_tol=1e-6)


def edges(h):
    ax = h.GetXaxis()
    return [float(ax.GetBinLowEdge(i)) for i in range(1, h.GetNbinsX()+1)] + [
        float(ax.GetBinUpEdge(h.GetNbinsX()))]


def matching(a, b):
    return len(a) == len(b) and all(near(x, y) for x, y in zip(a, b))


def index_of(e, value):
    return next((i for i, x in enumerate(e) if near(x, value)), None)


def is_1d_hist(h):
    return bool(h) and h.InheritsFrom('TH1') and h.GetDimension() == 1


def hist(f, name):
    h = f.Get(name)
    return h if is_1d_hist(h) else None


def pot(f):
    meta = f.Get('Meta')
    if not meta:
        return None, 'MISSING Meta'
    present = {b.GetName() for b in meta.GetListOfBranches()}
    for key in ('POT_Used', 'POTUsed', 'used_POT', 'UsedPOT',
                'POT', 'pot', 'POT_Total', 'total_pot'):
        if key in present:
            return sum(float(getattr(entry, key)) for entry in meta), key
    return None, 'MISSING POT branch (' + ','.join(sorted(present)) + ')'


def yield_range(h, lo, hi, scale=1.0):
    """Return (yield, statistical error), only for exactly aligned edges."""
    e = edges(h)
    i, j = index_of(e, lo), index_of(e, hi)
    if i is None or j is None or i >= j:
        return None
    val = sum(float(h.GetBinContent(b)) for b in range(i+1, j+1)) * scale
    err2 = sum(float(h.GetBinError(b))**2 for b in range(i+1, j+1))
    return val, math.sqrt(err2) * scale


def csvout(path, rows, fields):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def open_root(path):
    f = ROOT.TFile.Open(str(path))
    return f if f and not f.IsZombie() else None


def select_inputs(base, folder, playlist, tag, channel, binning_tag, notes):
    root_dir = base / folder
    if binning_tag.lower() != 'auto':
        stem = f'{playlist}_{tag}_{channel}_{binning_tag}_MAD.root'
        dp = root_dir / ('kin_dist_data' + stem)
        mp = root_dir / ('kin_dist_mc' + stem)
        return dp, mp
    # Discover matching productions using actual fine-binned EN4 edges,
    # never by arbitrary newest file or name alone.
    pattern = f'kin_dist_data{playlist}_{tag}_{channel}_*_MAD.root'
    candidates = []
    for dp in sorted(root_dir.glob(pattern)):
        suffix = dp.name[len('kin_dist_data'):]
        mp = root_dir / ('kin_dist_mc' + suffix)
        if not mp.is_file():
            continue
        df, mf = open_root(dp), open_root(mp)
        ok = bool(df and mf and hist(df, 'EN4') and hist(mf, 'EN4') and
                  matching(edges(hist(df, 'EN4')), FINE_EDGES) and
                  matching(edges(hist(mf, 'EN4')), FINE_EDGES))
        if df: df.Close()
        if mf: mf.Close()
        if ok:
            candidates.append((dp, mp))
    if len(candidates) != 1:
        notes.append(f'{folder}/{channel}: expected ONE fine-binned production, '
                     f'found {len(candidates)}; specify --binning-tag to disambiguate. '
                     f'Candidates: {[str(p[0]) for p in candidates]}')
        return None, None
    return candidates[0]


def category_hists(f, region, is_electron, e, notes, sample):
    allowed = E_CATEGORIES if is_electron else MU_CATEGORIES
    out = {}
    for name in sorted(allowed):
        h = hist(f, f'{region}_{name}')
        if h is not None:
            if matching(edges(h), e):
                out[name] = h
            else:
                notes.append(f'{sample}/{region}: mismatched component axes: {name}')
    if not out:
        notes.append(f'{sample}/{region}: NO recognized exclusive MC categories')
    return out


def summarize_ratio(mode, lo, hi, opened, notes):
    num = f'{mode}_ccnumu' if mode == 'fhc' else 'rhc_ccnumubar'
    den = f'{mode}_ccnue' if mode == 'fhc' else 'rhc_ccnuebar'
    if num not in opened or den not in opened:
        return None
    dn, mn, sn, pnum = opened[num]
    de, me, se, pden = opened[den]
    if not near(pnum, pden):
        notes.append(f'WARNING: {num}/{den} different data POT: {pnum:.6g}/{pden:.6g}')
    nd = yield_range(hist(dn,'EN4'),lo,hi)
    ed = yield_range(hist(de,'EN4'),lo,hi)
    nm = yield_range(hist(mn,'EN4'),lo,hi,sn)
    em = yield_range(hist(me,'EN4'),lo,hi,se)
    if any(x is None for x in [nd,ed,nm,em]):
        return None
    dn_y, _ = nd; de_y, _ = ed; mn_y, _ = nm; me_y, _ = em
    data_ratio = dn_y/de_y if de_y > 0 else None
    mc_ratio = mn_y/me_y if me_y > 0 else None
    double_ratio = data_ratio/mc_ratio if mc_ratio and data_ratio is not None else None
    # Approximate raw-count statistical uncertainty, not a full covariance.
    double_ratio_err = (double_ratio * math.sqrt(1/dn_y + 1/de_y)
                        if double_ratio is not None and dn_y > 0 and de_y > 0 else None)
    return dict(mode=mode, lo=lo, hi=hi, num_data=dn_y, den_data=de_y,
                num_mc=mn_y, den_mc=me_y, raw_data_ratio=data_ratio,
                raw_mc_ratio=mc_ratio, ratio_data_over_mc=double_ratio,
                raw_ratio_stat_err=double_ratio_err,
                note='RAW SELECTED EVENTS, not bkg-subtracted or sterile fit')



def fmt_number(value, spec=".3f"):
    if value is None:
        return "undefined"
    return format(value, spec)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base', type=Path, default=BASE)
    ap.add_argument('--out', type=Path, default=BASE/'surfaces'/'cv_fine_binning_diagnostics')
    ap.add_argument('--tag', default='ValidateSideband')
    ap.add_argument('--binning-tag', default='auto',
                    help='Production suffix, e.g. fineBinning; auto searches by exact bin edges')
    ap.add_argument('--suppress-plots', action='store_true')
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    notes, opened = [], {}
    inventories, pots, bins, comps, windows, ratios = [], [], [], [], [], []
    region_summaries = []; rebinned = []
    try:
        for sample, (folder, playlist, channel) in SAMPLES.items():
            dp, mp = select_inputs(args.base, folder, playlist, args.tag,
                                   channel, args.binning_tag, notes)
            if dp is None:
                continue
            if not dp.is_file() or not mp.is_file():
                notes.append(f'{sample}: MISSING input data={dp} mc={mp}')
                continue
            df, mf = open_root(dp), open_root(mp)
            if not df or not mf:
                notes.append(f'{sample}: ROOT OPEN FAILED data={dp} mc={mp}')
                if df: df.Close()
                if mf: mf.Close()
                continue
            pd, pdk = pot(df); pm, pmk = pot(mf)
            if not pd or not pm or pd <= 0 or pm <= 0:
                notes.append(f'{sample}: INVALID POT data={pd} ({pdk}) mc={pm} ({pmk})')
                df.Close(); mf.Close()
                continue
            hd, hm = hist(df, 'EN4'), hist(mf, 'EN4')
            if not hd or not hm:
                notes.append(f'{sample}: MISSING EN4 histograms')
                df.Close(); mf.Close()
                continue
            if not matching(edges(hd), FINE_EDGES) or not matching(edges(hm), FINE_EDGES):
                notes.append(f'{sample}: WRONG BINS; expected {FINE_EDGES}; '
                             f'data={edges(hd)} mc={edges(hm)}; SKIPPING')
                df.Close(); mf.Close()
                continue
            scale = pd/pm
            opened[sample] = (df, mf, scale, pd)
            pots.append(dict(sample=sample, data_file=str(dp), mc_file=str(mp),
                             data_pot=pd, mc_pot=pm, mc_to_data_scale=scale,
                             data_pot_branch=pdk, mc_pot_branch=pmk))
            for kind, f in [('data', df), ('mc', mf)]:
                for key in f.GetListOfKeys():
                    o = key.ReadObj()
                    if o.InheritsFrom('TH1'):
                        inventories.append(dict(sample=sample, kind=kind,
                                                histogram=key.GetName(), dimension=o.GetDimension(),
                                                bins_x=o.GetNbinsX(), xmin=o.GetXaxis().GetXmin(),
                                                xmax=o.GetXaxis().GetXmax()))
            is_electron = sample.endswith('ccnue') or sample.endswith('ccnuebar')
            regions = ['EN4'] + (['EN4_dEdX'] if is_electron else [])
            for region in regions:
                hd, hm = hist(df,region), hist(mf,region)
                if not hd or not hm:
                    notes.append(f'{sample}/{region}: missing paired histograms')
                    continue
                e = edges(hd)
                if not matching(e, FINE_EDGES) or not matching(e, edges(hm)):
                    notes.append(f'{sample}/{region}: wrong bin edges; SKIPPING')
                    continue
                parts = category_hists(mf, region, is_electron, e, notes, sample)
                region_d = yield_range(hd,e[0],e[-1])[0]
                region_m = yield_range(hm,e[0],e[-1],scale)[0]
                region_summaries.append(dict(sample=sample,region=region,
                                             data=region_d,mc=region_m,
                                             data_over_mc=region_d/region_m if region_m else None,
                                             bins=hd.GetNbinsX(), edges=';'.join(map(str,e))))
                for i,(lo,hi) in enumerate(zip(e[:-1],e[1:]),1):
                    d,de = yield_range(hd,lo,hi)
                    m,me = yield_range(hm,lo,hi,scale)
                    vals = {name: yield_range(h,lo,hi,scale)[0] for name,h in parts.items()}
                    cs = sum(vals.values()) if vals else None
                    valid = bool(vals) and abs(cs-m)<=max(1e-4, 1e-3*max(abs(m),1.0))
                    if vals and not valid:
                        notes.append(f'{sample}/{region} {lo:g}-{hi:g}: exclusive component '
                                     f'sum {cs:.6g} != total {m:.6g}; purity INVALID')
                    sig = sum(vals.get(cat,0.) for cat in E_SIGNAL_CATEGORIES) if valid and is_electron else None
                    pi0 = sum(vals.get(cat,0.) for cat in PI0_CATEGORIES) if valid and is_electron else None
                    bins.append(dict(sample=sample,region=region,bin=i,lo=lo,hi=hi,
                                     data=d,mc=m,data_stat_err=de,mc_stat_err=me,
                                     residual=d-m,data_over_mc=d/m if m>0 else None,
                                     raw_stat_residual=(d-m)/math.hypot(de,me) if de or me else None,
                                     components_sum=cs,components_valid=valid,
                                     signal_mc=sig,signal_fraction=sig/m if sig is not None and m>0 else None,
                                     pi0_mc=pi0,pi0_fraction=pi0/m if pi0 is not None and m>0 else None))
                    for name,val in vals.items():
                        comps.append(dict(sample=sample,region=region,bin=i,lo=lo,hi=hi,
                                          category=name,mc_at_data_pot=val,
                                          fraction=val/m if m else None))
                # Rebuild both existing analysis definitions from the same CV
                # production: no interpolation, and no refitting is implied.
                for label, requested in [('old_0p8', OLD_EDGES),
                                         ('test_1p0', TEST_EDGES)]:
                    for j,(blo,bhi) in enumerate(zip(requested[:-1], requested[1:]),1):
                        dy = yield_range(hd,blo,bhi)
                        my = yield_range(hm,blo,bhi,scale)
                        if dy is not None and my is not None:
                            rebinned.append(dict(sample=sample,region=region,
                                                 binning=label,bin=j,lo=blo,hi=bhi,
                                                 data=dy[0],mc=my[0],
                                                 data_stat_err=dy[1],mc_stat_err=my[1],
                                                 residual=dy[0]-my[0],
                                                 data_over_mc=dy[0]/my[0] if my[0]>0 else None,
                                                 source='SAME FINE CV PRODUCTION'))
                for lo,hi in LOW_WINDOWS:
                    d = yield_range(hd,lo,hi)
                    m = yield_range(hm,lo,hi,scale)
                    if d is None or m is None:
                        notes.append(f'{sample}/{region}: cannot integrate exact window {lo}-{hi}')
                        continue
                    windows.append(dict(sample=sample,region=region,lo=lo,hi=hi,
                                        data=d[0],mc=m[0],data_stat_err=d[1],mc_stat_err=m[1],
                                        residual=d[0]-m[0],data_over_mc=d[0]/m[0] if m[0]>0 else None,
                                        from_same_fine_production=True))
            print(f'{sample}: {dp.name}, POT scale={scale:.9g}, fine bins verified')

        # All available intervals, plus old/test analysis bins reconstructed exactly
        # from the SAME finer-binned histograms.
        intervals = sorted(set(LOW_WINDOWS +
                               list(zip(OLD_EDGES[:-1],OLD_EDGES[1:])) +
                               list(zip(TEST_EDGES[:-1],TEST_EDGES[1:]))))
        for mode in ['fhc','rhc']:
            for lo,hi in intervals:
                r = summarize_ratio(mode,lo,hi,opened,notes)
                if r is not None:
                    r['window_kind'] = ('old_and_test' if (lo,hi) in zip(OLD_EDGES[:-1],OLD_EDGES[1:])
                                        and (lo,hi) in zip(TEST_EDGES[:-1],TEST_EDGES[1:]) else
                                        'old' if (lo,hi) in zip(OLD_EDGES[:-1],OLD_EDGES[1:]) else
                                        'test' if (lo,hi) in zip(TEST_EDGES[:-1],TEST_EDGES[1:]) else 'diagnostic')
                    ratios.append(r)

        csvout(args.out/'pot_summary.csv',pots,
               ['sample','data_file','mc_file','data_pot','mc_pot','mc_to_data_scale',
                'data_pot_branch','mc_pot_branch'])
        csvout(args.out/'histogram_inventory.csv',inventories,
               ['sample','kind','histogram','dimension','bins_x','xmin','xmax'])
        csvout(args.out/'cv_bin_yields.csv',bins,
               ['sample','region','bin','lo','hi','data','mc','data_stat_err','mc_stat_err',
                'residual','data_over_mc','raw_stat_residual','components_sum',
                'components_valid','signal_mc','signal_fraction','pi0_mc','pi0_fraction'])
        csvout(args.out/'mc_components.csv',comps,
               ['sample','region','bin','lo','hi','category','mc_at_data_pot','fraction'])
        csvout(args.out/'low_energy_windows.csv',windows,
               ['sample','region','lo','hi','data','mc','data_stat_err','mc_stat_err',
                'residual','data_over_mc','from_same_fine_production'])
        csvout(args.out/'old_test_rebinned.csv',rebinned,
               ['sample','region','binning','bin','lo','hi','data','mc',
                'data_stat_err','mc_stat_err','residual','data_over_mc','source'])
        csvout(args.out/'region_summary.csv',region_summaries,
               ['sample','region','data','mc','data_over_mc','bins','edges'])
        csvout(args.out/'ratio_diagnostics.csv',ratios,
               ['mode','lo','hi','window_kind','num_data','den_data','num_mc','den_mc',
                'raw_data_ratio','raw_mc_ratio','ratio_data_over_mc','raw_ratio_stat_err','note'])
        with open(args.out/'summary.txt','w') as f:
            f.write('MINERvA LE fine-binned CV audit (read-only)\n')
            f.write('SCIENCE POLICY: Do not choose binning to eliminate a data/MC residual; '
                    'investigate physics, reconstruction and modeling causes first.\n')
            f.write('All energy thresholds are integrated from the SAME fine CV sample.\n')
            f.write('Raw selected-event mu/e ratios are NOT the background-subtracted sterile-fit ratios.\n')
            f.write('Stat errors are data/MC statistical errors ONLY, NOT full systematic uncertainties.\n')
            f.write('No background fit or sterile interpretation is performed.\n\n')
            f.write('Fine bin edges: '+str(FINE_EDGES)+'\n\n')
            for r in region_summaries:
                f.write(f"{r['sample']:17s} {r['region']:10s}: data={fmt_number(r['data'], '.3f')}, "
                        f"MC={fmt_number(r['mc'], '.3f')}, Data/MC={fmt_number(r['data_over_mc'], '.3f')}\n")
            f.write('\nLOW-ENERGY RAW DATA/MC (do not rank/select configurations):\n')
            for r in windows:
                if r['region']=='EN4' and (r['lo'],r['hi']) in [(0.8,1.0),(1.0,1.5),(1.5,2.0),(2.0,2.5)]:
                    f.write(f"{r['sample']:17s} [{r['lo']:3g},{r['hi']:3g}] "
                            f"data={r['data']:9.2f} MC={r['mc']:9.2f} "
                            f"D/MC={fmt_number(r['data_over_mc'], '.3f')}\n")
            f.write('\nWarnings and notes:\n')
            for n in notes: f.write(' - '+n+'\n')
        if not args.suppress_plots:
            try:
                import matplotlib
                matplotlib.use('Agg')
                import matplotlib.pyplot as plt
                from matplotlib.backends.backend_pdf import PdfPages
                with PdfPages(args.out/'diagnostic_plots.pdf') as pdf:
                    for sample in SAMPLES:
                        for region in ('EN4','EN4_dEdX'):
                            sub=[x for x in bins if x['sample']==sample and x['region']==region]
                            if not sub: continue
                            x=[(r['lo']+r['hi'])/2 for r in sub]
                            xe=[(r['hi']-r['lo'])/2 for r in sub]
                            fig,ax=plt.subplots(figsize=(10,5))
                            ax.axhline(1,linestyle='--',linewidth=1,alpha=.6)
                            ax.errorbar(x,[float('nan') if r['data_over_mc'] is None else r['data_over_mc'] for r in sub],xerr=xe,
                                        yerr=[r['data_stat_err']/r['mc'] if r['mc'] else 0 for r in sub],
                                        fmt='o',capsize=3)
                            ax.set_xlabel(r'Reconstructed $E_\nu$ (GeV)')
                            ax.set_ylabel('Data / prefit MC (statistical bars only)')
                            ax.set_title(f'{sample}: {region} fine CV')
                            ax.grid(alpha=.2); fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
                            # Close look at near-threshold data and MC, without hiding residuals.
                            low=[r for r in sub if r['hi']<=2.5+1e-7]
                            if low:
                                fig,ax=plt.subplots(figsize=(8,5))
                                lx=[(r['lo']+r['hi'])/2 for r in low]
                                lw=[r['hi']-r['lo'] for r in low]
                                ax.errorbar(lx,[r['data'] for r in low],
                                            yerr=[r['data_stat_err'] for r in low],
                                            xerr=[w/2 for w in lw],fmt='o',capsize=3,label='Data')
                                ax.step([r['lo'] for r in low]+[low[-1]['hi']],
                                        [r['mc'] for r in low]+[low[-1]['mc']],
                                        where='post',label='CV MC at data POT')
                                ax.set_xlabel(r'Reconstructed $E_\nu$ (GeV)')
                                ax.set_ylabel('Events / bin')
                                ax.set_title(f'{sample} {region}: low-energy counts')
                                ax.legend();ax.grid(alpha=.2)
                                fig.tight_layout();pdf.savefig(fig);plt.close(fig)
                    for mode in ['fhc','rhc']:
                        sub=[r for r in ratios if r['mode']==mode
                             and (r['lo'],r['hi']) in [(0.8,1),(1,1.5),(1.5,2),(2,2.5)]]
                        if not sub:continue
                        fig,ax=plt.subplots(figsize=(8,5))
                        ax.axhline(1,linestyle='--',linewidth=1,alpha=.6)
                        ax.errorbar([(r['lo']+r['hi'])/2 for r in sub],
                                    [float('nan') if r['ratio_data_over_mc'] is None else r['ratio_data_over_mc'] for r in sub],
                                    xerr=[(r['hi']-r['lo'])/2 for r in sub],
                                    yerr=[r['raw_ratio_stat_err'] if r['raw_ratio_stat_err'] is not None else 0.0 for r in sub],
                                    fmt='o',capsize=3)
                        ax.set_xlabel(r'Reconstructed $E_\nu$ (GeV)')
                        ax.set_ylabel('(raw data mu/e)/(raw MC mu/e)')
                        ax.set_title(f'{mode.upper()}: RAW selected-event ratios, not sterile fit')
                        ax.grid(alpha=.2);fig.tight_layout();pdf.savefig(fig);plt.close(fig)
            except ImportError:
                notes.append('Matplotlib unavailable; CSV and text outputs still produced')
        print('Wrote outputs to:',args.out)
        print('Validated sample count:',len(opened),'of',len(SAMPLES))
        for n in notes: print('NOTE:',n)
    finally:
        for df,mf,_,_ in opened.values():
            df.Close();mf.Close()


if __name__=='__main__':
    main()

#!/usr/bin/env python3
"""RHC MINOS-matched lepton-angle data/MC diagnostic; no new grid production.

Run: python studies/compare_rhc_lepton_angle.py
Optional: --mc-scale <verified MC->data exposure multiplier>
          --hist Lep_Angle  (or exact histogram name)
"""
import argparse
import os
import re
import numpy as np
import ROOT
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT.gROOT.SetBatch(True)

BEFORE = 'ValidateCuts_CCnumubar_FV_InverseHasNoBackExitingTracks_DeadTime_LeptonAngle_v1'
AFTER = 'ValidateCuts_CCnumubar_Inclusive_v1'


def open_file(base, sample, tag):
    path = os.path.join(base, f'kin_dist_{sample}le5_{tag}_MAD.root')
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError(f'Cannot open {path}')
    return f


def hist_candidates(f):
    keys = []
    for key in f.GetListOfKeys():
        name = key.GetName()
        obj = f.Get(name)
        if obj and obj.InheritsFrom('TH1') and not obj.InheritsFrom('TH2'):
            keys.append(name)
    return keys


def choose_hist(files, requested=None):
    common = set.intersection(*(set(hist_candidates(f)) for f in files))
    if requested:
        if requested in common:
            return requested
        matching = [n for n in common if n.lower() == requested.lower()]
        if len(matching) == 1:
            return matching[0]
        raise RuntimeError(f'Histogram {requested!r} not shared by all files. Angle-like common keys: {sorted(n for n in common if "angle" in n.lower() or "theta" in n.lower())}')
    # Prefer the named reco lepton angle, not truth/reco residual or category histograms.
    exact = ['Lep_Angle', 'Lepton_Angle', 'LeptonAngle', 'reco_theta_lep']
    for key in exact:
        if key in common:
            return key
    matches = [n for n in common if re.search(r'(lep.*angle|lepton.*angle)', n, re.I)
               and not re.search(r'(truth|true|resid|response|migration|_CC|_NC|_Other)', n, re.I)]
    if len(matches) == 1:
        return matches[0]
    raise RuntimeError('Unable to select a unique reco lepton-angle histogram. '
                       f'Candidates: {sorted(matches)}; all angle-like keys: '
                       f'{sorted(n for n in common if "angle" in n.lower() or "theta" in n.lower())}')


def read_hist(f, key):
    h = f.Get(key)
    if not h or not h.InheritsFrom('TH1') or h.InheritsFrom('TH2'):
        raise RuntimeError(f'Missing TH1 {key}')
    n = h.GetNbinsX()
    edges = np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1, n+2)], dtype=float)
    contents = np.array([h.GetBinContent(i) for i in range(1, n+1)], dtype=float)
    errors = np.array([h.GetBinError(i) for i in range(1, n+1)], dtype=float)
    if not np.all(np.isfinite(contents)) or not np.all(np.isfinite(errors)):
        raise RuntimeError(f'Nonfinite contents/errors in {key}')
    return edges, contents, errors


def safe_div(a, b):
    result = np.full(np.broadcast_shapes(np.shape(a), np.shape(b)), np.nan)
    np.divide(a, b, out=result, where=np.asarray(b) > 0)
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--directory', default='/exp/minerva/data/users/qvuong/antinu_mu')
    p.add_argument('--outdir', default='studies/rhc_lepton_angle')
    p.add_argument('--hist', default=None, help='Explicit ROOT TH1 name if auto-detection fails')
    p.add_argument('--mc-scale', type=float, default=None, help='Verified MC-to-data exposure multiplier; not fitted to the data')
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    files = {('before','data'): open_file(a.directory, 'data', BEFORE),
             ('before','mc'): open_file(a.directory, 'mc', BEFORE),
             ('after','data'): open_file(a.directory, 'data', AFTER),
             ('after','mc'): open_file(a.directory, 'mc', AFTER)}
    try:
        key = choose_hist(list(files.values()), a.hist)
        print(f'Using histogram: {key}')
        arrays = {k: read_hist(f, key) for k, f in files.items()}
        edges = arrays[('before','data')][0]
        if any(not np.allclose(edges, item[0]) for item in arrays.values()):
            raise RuntimeError('Inconsistent binning between productions')
        x = (edges[1:] + edges[:-1]) / 2
        widths = edges[1:] - edges[:-1]
        summary = {}
        for phase in ['before','after']:
            d, de = arrays[(phase,'data')][1:]
            m, me = arrays[(phase,'mc')][1:]
            raw = safe_div(d,m)
            D, M = np.sum(d), np.sum(m)
            shape = raw/(D/M) if D > 0 and M > 0 else raw * np.nan
            summary[phase] = {'d':d,'de':de,'m':m,'me':me,'raw':raw,'shape':shape,'D':D,'M':M}
            print(f'\n{phase.upper()}: visible data={D:.3f}, MC stored={M:.3f}, raw D/MC={D/M:.6f}')
            print(f'  under/overflow omitted from visible-bin calculations')
        before, after = summary['before'],summary['after']
        double = safe_div(after['raw'],before['raw'])
        print('\n Angle_low Angle_high     raw_before    raw_after   after/before   shape_before   shape_after')
        for i in range(len(x)):
            print(f'{edges[i]:10.3f} {edges[i+1]:10.3f} {before["raw"][i]:13.6f} '
                  f'{after["raw"][i]:12.6f} {double[i]:14.6f} '
                  f'{before["shape"][i]:14.6f} {after["shape"][i]:13.6f}')
        # Plot 1: data and MC normalized to unit area separately: distribution shape comparison.
        fig, axs = plt.subplots(1, 2, figsize=(12, 4.8), sharex=True)
        for ax, phase in zip(axs, ['before','after']):
            q = summary[phase]
            ax.stairs(q['m']/q['M']/widths, edges, label='MC shape')
            ax.errorbar(x, q['d']/q['D']/widths,
                        yerr=q['de']/q['D']/widths,
                        fmt='o', ms=3, capsize=1, label='Data shape')
            ax.set_title('MINOS matched, ' + ('no helicity' if phase=='before' else '+ helicity'))
            ax.set_xlabel('Reconstructed lepton angle (histogram units)')
            ax.grid(alpha=.25)
        axs[0].set_ylabel('Unit-area density')
        axs[1].legend()
        fig.tight_layout()
        fig.savefig(os.path.join(a.outdir, 'angle_shape_overlay.png'), dpi=180)
        plt.close(fig)
        # Plot 2: normalization independent shape-ratio data/MC in both selections.
        fig, ax = plt.subplots(figsize=(9,5.4))
        for phase, marker in [('before','o'), ('after','s')]:
            q=summary[phase]
            ax.plot(x,q['shape'],marker=marker,label=phase+' helicity' if phase=='after' else 'before helicity')
        ax.axhline(1,color='black',ls='--',lw=1)
        ax.set_xlabel('Reconstructed lepton angle (histogram units)')
        ax.set_ylabel('(Data/MC in bin) / (Data/MC integrated)')
        ax.legend();ax.grid(alpha=.25);fig.tight_layout()
        fig.savefig(os.path.join(a.outdir, 'angle_shape_data_mc_ratio.png'), dpi=180)
        plt.close(fig)
        # Plot 3: effect of helicity on D/MC; POT cancels if exposure convention identical.
        fig, ax=plt.subplots(figsize=(9,5.4))
        ax.plot(x,double,'o-')
        ax.axhline(1,color='black',ls='--',lw=1)
        ax.set_xlabel('Reconstructed lepton angle (histogram units)')
        ax.set_ylabel('(Data/MC after helicity) / (Data/MC before helicity)')
        ax.grid(alpha=.25);fig.tight_layout()
        fig.savefig(os.path.join(a.outdir,'angle_helicity_double_ratio.png'),dpi=180)
        plt.close(fig)
        if a.mc_scale is not None:
            if a.mc_scale <= 0:
                raise ValueError('--mc-scale must be positive')
            fig,ax=plt.subplots(figsize=(9,5.4))
            for phase,marker in [('before','o'),('after','s')]:
                q=summary[phase]
                ax.plot(x,q['raw']/a.mc_scale,marker=marker,label=phase)
            ax.axhline(1,color='black',ls='--',lw=1)
            ax.set_xlabel('Reconstructed lepton angle (histogram units)')
            ax.set_ylabel('POT-normalized data / MC')
            ax.legend();ax.grid(alpha=.25);fig.tight_layout()
            fig.savefig(os.path.join(a.outdir,'angle_pot_normalized_data_mc_ratio.png'),dpi=180)
            plt.close(fig)
        print(f'\nSaved figures to {a.outdir}')
        if a.mc_scale is None:
            print('No POT multiplier supplied: shape comparisons and helicity double ratio are valid; absolute D/MC remains unnormalized.')
        print('Check selection archives have identical MC weights/POT conventions before attributing differences to helicity.')
    finally:
        for f in files.values():
            f.Close()

if __name__ == '__main__':
    main()

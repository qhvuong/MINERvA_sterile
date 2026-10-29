#!/usr/bin/env python3
"""Figure-2-style pT panels of Eavail for a chosen p_parallel slice.

Input histogram: Eavail_Pt_Ppar_3D (X=Eavail, Y=pT, Z=p_parallel).
Plots differential densities, divided by dEavail * dpT * dp_parallel,
and multiplied by 1e-5 for display. Components are unstacked.
"""
import argparse
from pathlib import Path
import numpy as np
import ROOT
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)
BASE = '/exp/minerva/data/users/qvuong/nu_mu/'
TAG = 'ValidateDataMC_CCnumu_DansCuts_v4p4p1_LE_MHRW'
DEFAULT_DATA = BASE + f'kin_dist_dataleFHC_{TAG}_MAD.root'
DEFAULT_MC = BASE + f'kin_dist_mcleFHC_{TAG}_MAD.root'
STEM = 'Eavail_Pt_Ppar_3D'
SCALE = 1e-5
SIGNAL = [('QE', 'CCNuMuQE', '#27823d'),
          ('Resonant', 'CCNuMuDelta', '#3979c4'),
          ('2p2h', 'CCNuMu2p2h', '#df9229')]
BACKGROUND = ['CCNuMuDIS', 'CCNuMu', 'CC', 'Other', 'NuEElastic',
              'NCDiff', 'CCNuMuWrongSign', 'NC']


def get_th3(f, name):
    h = f.Get(name)
    if not h or not h.InheritsFrom('TH3'):
        raise RuntimeError(f'Missing TH3: {name} in {f.GetName()}')
    return h


def axis_edges(axis):
    return np.array([axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins() + 2)], dtype=float)


def density_slice(h, pt_bin, ppar_bin):
    ex = axis_edges(h.GetXaxis())
    volume = (np.diff(ex) * h.GetYaxis().GetBinWidth(pt_bin)
              * h.GetZaxis().GetBinWidth(ppar_bin))
    vals = np.array([h.GetBinContent(i, pt_bin, ppar_bin)
                     for i in range(1, h.GetNbinsX() + 1)], dtype=float)
    errs = np.array([h.GetBinError(i, pt_bin, ppar_bin)
                     for i in range(1, h.GetNbinsX() + 1)], dtype=float)
    return vals / volume, errs / volume


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data-file', default=DEFAULT_DATA)
    ap.add_argument('--mc-file', default=DEFAULT_MC)
    ap.add_argument('--mc-scale', type=float, required=True,
                    help='Data POT divided by MC POT')
    ap.add_argument('--ppar-bin', type=int, default=3, help='1-based p_parallel bin, default 3 = 2.5 to 3.0 GeV/c')
    ap.add_argument('--out', default='studies/dans_figure2_reco.png')
    ap.add_argument('--max-eavail', type=float, default=None)
    ap.add_argument('--panel-scales', type=float, nargs=9, default=[1.0]*9,
                    metavar='FACTOR', help='Nine optional cosmetic multipliers, one per pT panel')
    args = ap.parse_args()
    if args.mc_scale <= 0 or any(x <= 0 for x in args.panel_scales):
        ap.error('--mc-scale and all --panel-scales must be positive')
    data_file = ROOT.TFile.Open(args.data_file, 'READ')
    mc_file = ROOT.TFile.Open(args.mc_file, 'READ')
    if not data_file or data_file.IsZombie() or not mc_file or mc_file.IsZombie():
        raise RuntimeError('Cannot open ROOT data/MC files')
    try:
        data = get_th3(data_file, STEM)
        total = get_th3(mc_file, STEM)
        signals = [(label, get_th3(mc_file, STEM + '_' + suffix), color)
                   for label, suffix, color in SIGNAL]
        backgrounds = [get_th3(mc_file, STEM + '_' + suffix) for suffix in BACKGROUND]
        e_edges = axis_edges(data.GetXaxis())
        pt_edges = axis_edges(data.GetYaxis())
        par_edges = axis_edges(data.GetZaxis())
        if len(pt_edges) - 1 != 9:
            raise RuntimeError('4-4-1 layout expects nine pT bins')
        if not 1 <= args.ppar_bin < len(par_edges):
            ap.error(f'--ppar-bin must be 1 to {len(par_edges)-1}')
        for h in [total] + [h for _, h, _ in signals] + backgrounds:
            if any(not np.allclose(ref, axis_edges(axis), rtol=0, atol=1e-10)
                   for ref, axis in ((e_edges, h.GetXaxis()),
                                     (pt_edges, h.GetYaxis()),
                                     (par_edges, h.GetZaxis()))):
                raise RuntimeError(f'Histogram binning differs: {h.GetName()}')
        fig = plt.figure(figsize=(16, 11.5))
        gs = fig.add_gridspec(3, 4, left=0.075, right=0.98, bottom=0.085,
                              top=0.89, hspace=0.48, wspace=0.28)
        axes = [fig.add_subplot(gs[r, c]) for r, c in
                [(0,0),(0,1),(0,2),(0,3),(1,0),(1,1),(1,2),(1,3),(2,0)]]
        e_centers = 0.5 * (e_edges[:-1] + e_edges[1:])
        max_mismatch = 0
        for j, ax in enumerate(axes, start=1):
            factor = args.panel_scales[j-1]
            d, de = density_slice(data, j, args.ppar_bin)
            t, _ = density_slice(total, j, args.ppar_bin)
            d = d * SCALE * factor
            de = de * SCALE * factor
            t = t * args.mc_scale * SCALE * factor
            components = []
            for label, h, color in signals:
                v, _ = density_slice(h, j, args.ppar_bin)
                components.append((label, v * args.mc_scale * SCALE * factor, color))
            bg = sum((density_slice(h, j, args.ppar_bin)[0] for h in backgrounds),
                     np.zeros(len(e_centers))) * args.mc_scale * SCALE * factor
            predicted = bg + sum((v for _, v, _ in components), np.zeros_like(bg))
            max_mismatch = max(max_mismatch, float(np.max(np.abs(predicted - t))))
            if not np.allclose(predicted, t, rtol=1e-4, atol=1e-7):
                print(f'[WARN] pT panel {j}: signal + background differs from inclusive MC')
            ax.stairs(t, e_edges, color='#bb2a27', linewidth=1.8)
            for _, v, color in components:
                ax.stairs(v, e_edges, color=color, linewidth=1.5)
            ax.stairs(bg, e_edges, color='#a74491', linewidth=1.5)
            ax.errorbar(e_centers, d, yerr=de, color='black', fmt='o',
                        ms=3.1, capsize=0, linewidth=1)
            ax.set_title(rf'${pt_edges[j-1]:g} < p_T < {pt_edges[j]:g}$ GeV/c', fontsize=10)
            if factor != 1:
                ax.text(0.97, 0.94, rf'$\times {factor:g}$',
                        transform=ax.transAxes, va='top', ha='right', fontsize=10)
            ax.set_xlabel(r'$E_{\mathrm{avail}}$ (GeV)', fontsize=10)
            ax.set_ylabel(r'$10^{-5}$ events / (GeV$^3$/$c^2$)', fontsize=9)
            ax.set_ylim(bottom=0)
            ax.grid(alpha=0.15, axis='y')
            if args.max_eavail is not None:
                ax.set_xlim(e_edges[0], args.max_eavail)
            volume = (np.diff(e_edges) * (pt_edges[j]-pt_edges[j-1])
                      * (par_edges[args.ppar_bin]-par_edges[args.ppar_bin-1]))
            # Undo cosmetic scaling to print actual weighted event totals.
            nev_data = np.sum((d / (SCALE*factor)) * volume)
            nev_mc = np.sum((t / (SCALE*factor)) * volume)
            print(f'pT {pt_edges[j-1]:g}-{pt_edges[j]:g}: '
                  f'data={nev_data:.6g}, POT-scaled MC={nev_mc:.6g}, '
                  f'D/MC={nev_data/nev_mc if nev_mc else float("nan"):.4g}')
        legend = [Line2D([0], [0], color='black', marker='o', linestyle='None', label='Data'),
                  Line2D([0], [0], color='#bb2a27', lw=2, label='Total MC')]
        legend += [Line2D([0], [0], color=color, lw=2, label=label)
                   for label, _, color in SIGNAL]
        legend.append(Line2D([0], [0], color='#a74491', lw=2, label='Background'))
        fig.legend(handles=legend, loc='center', bbox_to_anchor=(0.64, 0.17),
                   fontsize=12, frameon=False, ncol=2)
        fig.suptitle('Reconstructed available energy in muon $p_T$ bins\n'
                     + rf'${par_edges[args.ppar_bin-1]:g} < p_{{\parallel}} < '
                     + rf'{par_edges[args.ppar_bin]:g}$ GeV/c', fontsize=15)
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=180, bbox_inches='tight')
        fig.savefig(out.with_suffix('.pdf'), bbox_inches='tight')
        plt.close(fig)
        print('Saved:', out, 'and', out.with_suffix('.pdf'))
        print('Maximum displayed category-sum discrepancy:', max_mismatch)
    finally:
        data_file.Close()
        mc_file.Close()


if __name__ == '__main__':
    main()

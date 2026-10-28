#!/usr/bin/env python3
"""Figure-2-style reconstructed Eavail panels in reconstructed muon pT bins.

Assumes the desired p_parallel cut was applied before histogram filling.
MC exposure scaling is passed explicitly with --mc-scale.
"""
import argparse
from pathlib import Path

import ROOT
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

BASE = '/exp/minerva/data/users/qvuong/nu_mu/'
STEM = 'Lepton_Pt_vs_Eavail'
DEFAULT_DATA = BASE + 'kin_dist_dataleFHC_ValidateDataMC_CCnumu_DansCuts_v4p4p1_LE_MHRW_MAD.root'
DEFAULT_MC = BASE + 'kin_dist_mcleFHC_ValidateDataMC_CCnumu_DansCuts_v4p4p1_LE_MHRW_MAD.root'

# These are mutually exclusive classifications in this analysis.
CATEGORIES = [
    ('QE', 'CCNuMuQE'),
    ('2p2h', 'CCNuMu2p2h'),
    ('RES / Delta', 'CCNuMuDelta'),
    ('DIS', 'CCNuMuDIS'),
    ('Other CCNuMu', 'CCNuMu'),
    ('Other CC', 'CC'),
    ('Other', 'Other'),
    ('Nu-e elastic', 'NuEElastic'),
    ('NC diffractive', 'NCDiff'),
    ('Wrong-sign CCNuMu', 'CCNuMuWrongSign'),
    ('NC', 'NC'),
]


def get_hist(root_file, name):
    hist = root_file.Get(name)
    if not hist or not hist.InheritsFrom('TH2'):
        raise RuntimeError(f'Missing TH2 histogram: {name}')
    return hist


def edges(axis):
    return np.array([axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins() + 2)], dtype=float)


def pt_slice(hist, ix):
    return (np.array([hist.GetBinContent(ix, iy) for iy in range(1, hist.GetNbinsY() + 1)], dtype=float),
            np.array([hist.GetBinError(ix, iy) for iy in range(1, hist.GetNbinsY() + 1)], dtype=float))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-file', default=DEFAULT_DATA)
    parser.add_argument('--mc-file', default=DEFAULT_MC)
    parser.add_argument('--mc-scale', type=float, required=True,
                        help='Data POT / MC POT, or 1 if MC histograms already POT-scaled')
    parser.add_argument('--out', default='studies/dans_figure2_reco.png')
    parser.add_argument('--max-eavail', type=float, default=None)
    parser.add_argument('--columns', type=int, default=3)
    parser.add_argument('--ppar-label', default=r'$2.5 < p_{\parallel} < 3.0$ GeV/c')
    args = parser.parse_args()
    if args.mc_scale <= 0 or args.columns < 1:
        parser.error('--mc-scale and --columns must be positive')

    data_file = ROOT.TFile.Open(args.data_file, 'READ')
    mc_file = ROOT.TFile.Open(args.mc_file, 'READ')
    if not data_file or data_file.IsZombie() or not mc_file or mc_file.IsZombie():
        raise RuntimeError('Cannot open data or MC ROOT file')

    data = get_hist(data_file, STEM)
    total = get_hist(mc_file, STEM)
    components = [(label, get_hist(mc_file, STEM + '_' + suffix)) for label, suffix in CATEGORIES]
    pt_edges, e_edges = edges(data.GetXaxis()), edges(data.GetYaxis())
    for label, hist in [('Total MC', total)] + components:
        if (not np.array_equal(pt_edges, edges(hist.GetXaxis())) or
                not np.array_equal(e_edges, edges(hist.GetYaxis()))):
            raise RuntimeError(f'Histogram binning mismatch: {label}')

    npt = len(pt_edges) - 1
    ncols = min(args.columns, npt)
    nrows = (npt + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.1 * ncols, 3.8 * nrows), squeeze=False)
    e_centers = (e_edges[:-1] + e_edges[1:]) / 2
    widths = np.diff(e_edges)
    colors = plt.get_cmap('tab20').colors[:len(components)]
    max_abs_diff = 0.0

    for ix in range(1, npt + 1):
        ax = axes[(ix - 1) // ncols, (ix - 1) % ncols]
        observed, observed_err = pt_slice(data, ix)
        nominal, nominal_err = pt_slice(total, ix)
        nominal *= args.mc_scale
        nominal_err *= args.mc_scale
        slices = [(label, pt_slice(hist, ix)[0] * args.mc_scale) for label, hist in components]
        stacked = np.sum([v for _, v in slices], axis=0)
        diff = stacked - nominal
        max_abs_diff = max(max_abs_diff, float(np.max(np.abs(diff))))
        bad = np.abs(diff) > (1e-4 * np.maximum(1.0, np.abs(nominal)))
        if np.any(bad):
            print(f'[WARNING] pT bin {ix}: category sum != total MC in {np.count_nonzero(bad)} Eavail bins; '
                  f'max |sum-total|={np.max(np.abs(diff)):.5g}')

        bottom = np.zeros_like(nominal)
        for j, (label, values) in enumerate(slices):
            ax.bar(e_edges[:-1], values, widths, bottom=bottom, align='edge',
                   color=colors[j], linewidth=0, label=label)
            bottom += values

        ax.stairs(nominal, e_edges, color='black', linewidth=1.3, label='Total MC')
        ax.errorbar(e_centers, observed, yerr=observed_err, fmt='ko', ms=3,
                    capsize=0, label='Data')
        ax.set_title(f'{pt_edges[ix-1]:g} < $p_T$ < {pt_edges[ix]:g} GeV/c')
        ax.set_xlabel(r'$E_{\mathrm{avail}}$ (GeV)')
        ax.set_ylabel('Events / Eavail bin')
        ax.set_ylim(bottom=0)
        ax.grid(alpha=0.17, axis='y')
        if args.max_eavail is not None:
            ax.set_xlim(e_edges[0], args.max_eavail)
        ratio = np.sum(observed) / np.sum(nominal) if np.sum(nominal) > 0 else np.nan
        print(f'pT [{pt_edges[ix-1]:g}, {pt_edges[ix]:g}]: '
              f'data={np.sum(observed):.4g}, MC={np.sum(nominal):.4g}, '
              f'stack={np.sum(stacked):.4g}, D/MC={ratio:.4g}')

    for ax in axes.flat[npt:]:
        ax.set_visible(False)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.985),
               ncol=5, fontsize=9)
    fig.suptitle('Reconstructed available energy in muon $p_T$ bins\n' + args.ppar_label,
                 fontsize=15, y=1.075)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180, bbox_inches='tight')
    fig.savefig(out.with_suffix('.pdf'), bbox_inches='tight')
    plt.close(fig)
    print('Saved:', out)
    print('Max absolute category-sum vs total-MC difference:', max_abs_diff)
    print('Verify Ppar selection stage and MC POT scaling before physics comparisons.')
    data_file.Close()
    mc_file.Close()


if __name__ == '__main__':
    main()

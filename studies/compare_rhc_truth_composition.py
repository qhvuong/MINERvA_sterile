#!/usr/bin/env python3
"""Compare RHC truth-category composition with/without reconstructed helicity.

Uses existing CV MnvH1D histograms. Does not apply POT scaling: composition
fractions and same-sample helicity efficiencies do not require it.
"""
import os
import numpy as np
import ROOT
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT.gROOT.SetBatch(True)

DIRECTORY = '/exp/minerva/data/users/qvuong/antinu_mu'
TAGS = {
    'Before helicity': 'ValidateCuts_CCnumubar_FV_InverseHasNoBackExitingTracks_DeadTime_LeptonAngle_v1',
    'After helicity': 'ValidateCuts_CCnumubar_Inclusive_v1',
}
OUTDIR = 'studies/rhc_truth_composition'
os.makedirs(OUTDIR, exist_ok=True)

# First-match, mutually exclusive classification from CCAntiNuMu/SignalDef.py.
RIGHT = ['CCNuMuQE', 'CCNuMuDelta', 'CCNuMuDIS', 'CCNuMu2p2h', 'CCNuMu']
WRONG = ['CCNuMuWrongSign']
OTHER = ['CC', 'NC', 'NuEElastic', 'NCDiff', 'Other']
CATEGORIES = WRONG + RIGHT + OTHER

def read_hist(f, name):
    h = f.Get(name)
    if not h:
        raise RuntimeError(f'Missing {name} in {f.GetName()}')
    edges = np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1, h.GetNbinsX()+2)])
    vals = np.array([h.GetBinContent(i) for i in range(1, h.GetNbinsX()+1)], dtype=float)
    if not np.isfinite(vals).all():
        raise RuntimeError(f'Nonfinite CV bins in {name} of {f.GetName()}')
    return edges, vals

def load(tag):
    path = os.path.join(DIRECTORY, f'kin_dist_mcle5_{tag}_MAD.root')
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError(f'Cannot open {path}')
    try:
        edges, inclusive = read_hist(f, 'Eel')
        cats = {}
        for name in CATEGORIES:
            e, v = read_hist(f, 'Eel_' + name)
            if not np.allclose(e, edges, rtol=0, atol=1e-9):
                raise RuntimeError(f'Binning mismatch in {name}: {path}')
            cats[name] = v
    finally:
        f.Close()
    total = sum(cats.values())
    residual = inclusive - total
    print(f'\n=== {tag} ===')
    print('Visible inclusive MC:', inclusive.sum())
    print('Visible sum of categories:', total.sum())
    print('Closure difference:', residual.sum())
    print('Max absolute bin difference:', np.max(np.abs(residual)))
    for i in np.where(np.abs(residual) > np.maximum(1e-6, 1e-6*np.abs(inclusive)))[0]:
        print(f'  CLOSURE WARNING {edges[i]:.2f}-{edges[i+1]:.2f}: inclusive={inclusive[i]:.7g} categories={total[i]:.7g} residual={residual[i]:.7g}')
    right = sum(cats[x] for x in RIGHT)
    wrong = sum(cats[x] for x in WRONG)
    other = sum(cats[x] for x in OTHER)
    print(f'Right sign: {right.sum():.5f} | wrong sign: {wrong.sum():.5f} | other: {other.sum():.5f}')
    denom = inclusive.sum()
    if denom:
        print(f'Fractions: right={right.sum()/denom:.4%}, wrong={wrong.sum()/denom:.4%}, other={other.sum()/denom:.4%}')
    return {'edges':edges,'all':inclusive,'cats':cats,'right':right,'wrong':wrong,'other':other}

results = {label: load(tag) for label,tag in TAGS.items()}
before = results['Before helicity']
after = results['After helicity']
if not np.allclose(before['edges'], after['edges']):
    raise RuntimeError('Before and after histograms have different binning')
edges = before['edges']
centers = (edges[:-1] + edges[1:])/2
widths = np.diff(edges)

print('\n=== PER-BIN COMPOSITION (weighted MC, visible bins) ===')
print('E_low E_high   before_right% before_wrong% after_right% after_wrong%  helicity_eff_right  helicity_eff_wrong')
for i in range(len(centers)):
    def frac(part, whole):
        return part/whole if whole > 0 else float('nan')
    print(f'{edges[i]:5.2f} {edges[i+1]:6.2f}'
          f'   {100*frac(before["right"][i],before["all"][i]):10.2f}'
          f' {100*frac(before["wrong"][i],before["all"][i]):11.2f}'
          f' {100*frac(after["right"][i],after["all"][i]):11.2f}'
          f' {100*frac(after["wrong"][i],after["all"][i]):11.2f}'
          f' {frac(after["right"][i],before["right"][i]):19.3f}'
          f' {frac(after["wrong"][i],before["wrong"][i]):19.3f}')

# Plot stacked individual categories (each selection on separate axes)
fig, axes = plt.subplots(2, 1, figsize=(11, 10), sharex=True, constrained_layout=True)
for ax, (label, res) in zip(axes, results.items()):
    bottom = np.zeros_like(res['all'])
    for name in CATEGORIES:
        values = res['cats'][name]
        ax.bar(edges[:-1], values, widths, bottom=bottom, align='edge', label=name, linewidth=0)
        bottom += values
    ax.stairs(res['all'], edges, color='black', linewidth=1.4, label='Inclusive Eel')
    ax.set_ylabel('Weighted MC yield / bin')
    ax.set_title(label)
    ax.legend(ncol=4, fontsize=8)
axes[-1].set_xlabel('Reconstructed lepton energy (GeV)')
fig.savefig(os.path.join(OUTDIR, 'stacked_truth_categories.png'), dpi=180)
plt.close(fig)

# A composition fraction figure makes comparison between selections easier.
fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True, constrained_layout=True)
for ax, (label, res) in zip(axes, results.items()):
    if np.any(res['all'] <= 0):
        print('WARNING: nonpositive inclusive MC bin; omitting percentages there')
    with np.errstate(invalid='ignore', divide='ignore'):
        y_right = np.where(res['all'] > 0, 100*res['right']/res['all'], 0)
        y_wrong = np.where(res['all'] > 0, 100*res['wrong']/res['all'], 0)
        y_other = np.where(res['all'] > 0, 100*res['other']/res['all'], 0)
    ax.bar(edges[:-1], y_right, widths, align='edge', label='Right-sign anti-numu')
    ax.bar(edges[:-1], y_wrong, widths, bottom=y_right, align='edge', label='Wrong-sign numu')
    ax.bar(edges[:-1], y_other, widths, bottom=y_right+y_wrong, align='edge', label='Other')
    ax.set_ylim(0, 105)
    ax.set_ylabel('MC composition (%)')
    ax.set_title(label)
    ax.legend(loc='lower right')
axes[-1].set_xlabel('Reconstructed lepton energy (GeV)')
fig.savefig(os.path.join(OUTDIR, 'truth_flavor_fractions.png'), dpi=180)
plt.close(fig)

# Right-sign/wrong-sign survival relative to their own pre-helicity yields.
fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
for key, label in [('right', 'Right-sign anti-numu'), ('wrong', 'Wrong-sign numu'), ('other', 'Other')]:
    with np.errstate(divide='ignore', invalid='ignore'):
        rate = np.where(before[key] > 0, after[key]/before[key], np.nan)
    ax.stairs(rate, edges, label=label)
ax.set_xlabel('Reconstructed lepton energy (GeV)')
ax.set_ylabel('MC survival fraction (after / before)')
ax.set_title('Helicity cut survival by true-flavor component')
ax.legend()
fig.savefig(os.path.join(OUTDIR, 'helicity_survival_by_truth_group.png'), dpi=180)
plt.close(fig)

print('\nSaved plots under:', OUTDIR)
print('NOTE: Composition uses weighted MC as stored; no POT scaling or separately produced 2p2h file is added.')

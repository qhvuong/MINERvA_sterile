
#!/usr/bin/env python3
"""Read-only Eavail diagnostic: CV tune changes in electron-kinematic slices.
Run on MINERvA GPVM with PyROOT, NumPy, Matplotlib. No event selection reruns.
"""
import argparse
import csv
from pathlib import Path

import ROOT  # Load ROOT first to avoid FNAL libstdc++ conflicts.
ROOT.gROOT.SetBatch(True)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = Path('/exp/minerva/data/users/qvuong')
TUNES = {
    'V1': 'ValidateTune_{sample}_noEavail_v1',
    'V431': 'ValidateTune_{sample}_noEavail_v4p3p1',
    'V431_noMENU1PI': 'ValidateTune_{sample}_noEavail_v4p3p1_noMENU1PI',
    'V441': 'ValidateTune_{sample}_noEavail_v4p4p1',
}
# x: electron variable, y: available energy (as defined in PlotLibrary)
VARIABLES = {
    'Elepton_Eavail': 'Electron energy',
    'Lepton_Pt_vs_Eavail': 'Electron pT',
    'Lepton_Ppar_vs_Eavail': 'Electron p-parallel',
}
CATEGORIES = ('CCNuE', 'CCNuEQE', 'CCNuE2p2h', 'CCNuEDelta',
              'CCNuEDIS', 'CCPi', 'CCPi0', 'NCPi', 'NCPi0',
              'CCNuEWrongSign', 'NCCohPi0')
KEEP = []


def read_file(p):
    f = ROOT.TFile.Open(str(p), 'READ')
    if not f or f.IsZombie():
        raise RuntimeError('Cannot read ' + str(p))
    KEEP.append(f)  # Keep ROOT files alive while reading objects.
    meta = f.Get('Meta')
    if not meta or not meta.GetBranch('POT_Used'):
        raise RuntimeError('Missing Meta.POT_Used: ' + str(p))
    pot = sum(float(entry.POT_Used) for entry in meta)
    if pot <= 0:
        raise RuntimeError('Bad POT: ' + str(p))
    return f, pot


def read_h2(f, key):
    h = f.Get(key)
    if not h or not h.InheritsFrom('TH2'):
        return None
    x = np.array([h.GetXaxis().GetBinLowEdge(i)
                  for i in range(1, h.GetNbinsX() + 2)])
    y = np.array([h.GetYaxis().GetBinLowEdge(i)
                  for i in range(1, h.GetNbinsY() + 2)])
    z = np.array([[h.GetBinContent(i, j)
                   for j in range(1, h.GetNbinsY() + 1)]
                  for i in range(1, h.GetNbinsX() + 1)], dtype=float)
    return x, y, z


def same_bins(a, b):
    return len(a) == len(b) and np.allclose(a, b, atol=1e-10, rtol=1e-9)


def div(a, b):
    return np.divide(a, b, out=np.full_like(a, np.nan, dtype=float), where=b != 0)


def write_csv(path, rows, fields):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def make_plots(out, stem, title, xlabel, x, y, data, mc):
    # Data/MC panels, plotted versus Eavail for each electron-variable slice.
    n = len(x) - 1
    fig, axes = plt.subplots(n, 1, figsize=(9, max(3.4, 2.25 * n)), squeeze=False)
    yc = (y[:-1] + y[1:]) / 2
    for i, ax in enumerate(axes[:, 0]):
        for tune, values in mc.items():
            ax.plot(yc, div(data[i], values[i]), '.-', markersize=4, label=tune)
        ax.axhline(1, color='0.45', linestyle='--', linewidth=1)
        ax.set_ylabel('Data / MC')
        ax.set_title(f'{xlabel}: {x[i]:.3g}–{x[i+1]:.3g} GeV', fontsize=10)
        ax.grid(alpha=0.25)
        if i == 0:
            ax.legend(fontsize=8, ncol=2)
    axes[-1, 0].set_xlabel('Available energy (GeV)')
    fig.suptitle(title, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    fig.savefig(out / f'{stem}_sliced_datamc.pdf')
    plt.close(fig)

    # Heatmaps: V441-V1 increase, and residual Data-V441.
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    fields = [mc['V441'] - mc['V1'], data - mc['V441']]
    headings = ['V4.4.1 - V1 (MC events)', 'Data - V4.4.1 (events)']
    for ax, values, heading in zip(axes, fields, headings):
        mesh = ax.pcolormesh(x, y, values.T, shading='auto', cmap='coolwarm')
        fig.colorbar(mesh, ax=ax)
        ax.set_xlabel(xlabel + ' (GeV)')
        ax.set_ylabel('Available energy (GeV)')
        ax.set_title(heading)
    fig.suptitle(title)
    fig.savefig(out / f'{stem}_heatmaps.png', dpi=160)
    plt.close(fig)


def run_sample(base, out, sample, playlist, summary):
    folder = 'nu_e' if sample == 'CCnue' else 'antinu_e'
    mcs = {}
    data = None
    for name, template in TUNES.items():
        tag = template.format(sample=sample)
        mc_path = base / folder / f'kin_dist_mc{playlist}_{tag}_MAD.root'
        data_path = base / folder / f'kin_dist_data{playlist}_{tag}_MAD.root'
        if not mc_path.exists():
            summary.append('MISSING FILE ' + str(mc_path))
            return
        fm, pm = read_file(mc_path)
        fd, pd = read_file(data_path)
        mcs[name] = (fm, pd / pm)
        if data is None:
            data = fd
        summary.append(f'{sample} {name}: MC POT={pm:.6g}; data POT={pd:.6g}')
    for key, xlabel in VARIABLES.items():
        h = read_h2(data, key)
        if h is None:
            summary.append(f'{sample} {key}: missing DATA histogram; skip')
            continue
        x, y, d = h
        pred = {}
        for name, (f, scale) in mcs.items():
            v = read_h2(f, key)
            if v is None or not same_bins(v[0], x) or not same_bins(v[1], y):
                summary.append(f'{sample} {key} {name}: missing or mismatched MC; skip variable')
                pred = {}
                break
            pred[name] = v[2] * scale
        if len(pred) != 4:
            continue
        dest = out / sample / key
        dest.mkdir(parents=True, exist_ok=True)
        title = f'{sample}  {key}'
        make_plots(dest, key, title, xlabel, x, y, d, pred)
        rows = []
        for i in range(len(x) - 1):
            for j in range(len(y) - 1):
                row = dict(x_low=x[i], x_high=x[i+1], eavail_low=y[j],
                           eavail_high=y[j+1], data=d[i,j],
                           **{name: pred[name][i,j] for name in TUNES})
                row['V441_minus_V1'] = pred['V441'][i,j] - pred['V1'][i,j]
                row['data_minus_V441'] = d[i,j] - pred['V441'][i,j]
                row['data_over_V441'] = float(div(np.array([d[i,j]]), np.array([pred['V441'][i,j]]))[0])
                rows.append(row)
        write_csv(dest / 'bins_2d.csv', rows, list(rows[0]))
        summary.append(f'\n{sample} | {key}: data={d.sum():.1f}; '
                       + '; '.join(f'{name}={pred[name].sum():.1f}' for name in TUNES))
        summary.append('  NOTE: 2D integral is limited to the saved x AND Eavail ranges.')
        slices = []
        for i in range(len(x) - 1):
            dd = float(d[i].sum())
            a = float(pred['V1'][i].sum())
            b = float(pred['V441'][i].sum())
            slices.append(dict(x_low=x[i], x_high=x[i+1], data=dd,
                               V1=a, V441=b, delta=b-a,
                               D_over_V1=dd/a if a else float('nan'),
                               D_over_V441=dd/b if b else float('nan')))
        write_csv(dest / 'x_slices.csv', slices, list(slices[0]))
        for s in slices:
            summary.append('  x=[{x_low:.3g},{x_high:.3g}] data={data:.0f} '
                           'V1={V1:.1f} V441={V441:.1f} increase={delta:+.1f} '
                           'D/V441={D_over_V441:.3f}'.format(**s))
        biggest = sorted(rows, key=lambda r: r['V441_minus_V1'], reverse=True)[:8]
        summary.append('  Largest positive MC increases (x bin, Eavail bin, delta, Data-V441):')
        for r in biggest:
            summary.append(f"    x {r['x_low']:.3g}-{r['x_high']:.3g}; "
                           f"Eavail {r['eavail_low']:.3g}-{r['eavail_high']:.3g}; "
                           f"delta={r['V441_minus_V1']:+.1f}; "
                           f"residual={r['data_minus_V441']:+.1f}")
        # Interaction-category TH2, when available.
        category_rows = []
        for category in CATEGORIES:
            name = f'{key}_{category}'
            h1 = read_h2(mcs['V1'][0], name)
            h4 = read_h2(mcs['V441'][0], name)
            if h1 is None or h4 is None:
                continue
            if not (same_bins(h1[0], x) and same_bins(h1[1], y)
                    and same_bins(h4[0], x) and same_bins(h4[1], y)):
                continue
            change = h4[2] * mcs['V441'][1] - h1[2] * mcs['V1'][1]
            category_rows.append(dict(category=category, V441_minus_V1=change.sum()))
        if category_rows:
            write_csv(dest / 'category_changes.csv', category_rows,
                      ['category', 'V441_minus_V1'])
            summary.append('  2D interaction categories (V441-V1): ' + ', '.join(
                f"{z['category']}={z['V441_minus_V1']:+.1f}" for z in
                sorted(category_rows, key=lambda r: abs(r['V441_minus_V1']), reverse=True)))
        else:
            summary.append('  No compatible 2D component histograms; use existing 1D component CSVs.')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--base', type=Path, default=BASE)
    p.add_argument('--outdir', type=Path,
                   default=BASE / 'studies' / 'eavail_physics_diagnostics')
    p.add_argument('--fhc-playlist', choices=['le1','leFHC'], default='le1')
    args = p.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    summary = ["Eavail tune diagnostic - POT-normalized, CV-only, saved histogram ranges"]
    run_sample(args.base, args.outdir, 'CCnue', args.fhc_playlist, summary)
    run_sample(args.base, args.outdir, 'CCnuebar', 'le5', summary)
    report = '\n'.join(summary)
    (args.outdir / 'SUMMARY.txt').write_text(report + '\n')
    print(report)
    print('\nREPORT:', args.outdir / 'SUMMARY.txt')


if __name__ == '__main__':
    main()

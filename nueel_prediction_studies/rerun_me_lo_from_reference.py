#!/usr/bin/env python3
"""Recompute LO using exposure and flux inputs recorded in an existing prediction."""
import argparse
import csv
import json
from pathlib import Path
import shlex
import subprocess
import sys


def records(meta):
    if 'exposures' in meta:
        return [meta]
    return [r for item in meta.get('inputs', [])
            for r in records(item.get('metadata', item))]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--sin2thetaW', type=float, default=0.23129)
    a = p.parse_args()
    if a.reference.resolve() == a.output.resolve() or a.output.exists():
        p.error('Choose a new output path; existing predictions will not be overwritten.')
    import ROOT
    import PlotUtils.LoadPlotUtilsLib
    import numpy as np
    ROOT.gROOT.SetBatch(True)
    ref = ROOT.TFile.Open(str(a.reference), 'READ')
    if not ref or ref.IsZombie() or ref.TestBit(ROOT.TFile.kRecovered):
        raise RuntimeError('Invalid reference file: ' + str(a.reference))
    obj = ref.Get('prediction_metadata_json')
    if not obj:
        raise RuntimeError('Missing reference metadata')
    meta = json.loads(str(obj.GetString().Data()))
    recs = records(meta)
    if not recs:
        raise RuntimeError('Cannot find original prediction inputs in reference metadata')
    beam = recs[0]['beam']
    if any(r['beam'] != beam or r.get('response_json') for r in recs):
        raise RuntimeError('Expected one beam and the existing true-energy prediction')
    shared = ('target_electrons', 'input_flux_area_unit', 'true_edges_GeV')
    for key in shared:
        if any(r[key] != recs[0][key] for r in recs):
            raise RuntimeError('Inconsistent reference setting: ' + key)
    exposures = {}
    sources = []
    for r in recs:
        for playlist, pot in r['exposures'].items():
            if playlist in exposures:
                raise RuntimeError('Duplicate playlist: ' + playlist)
            exposures[playlist] = pot
        sources.extend(r['inputs'])
    parents = {Path(s['path']).parent for s in sources}
    if len(parents) != 1:
        raise RuntimeError('Flux inputs must share a directory for this runner')
    flux_dir = parents.pop()
    patterns = ('{playlist}_{flavor}.root',
                'flux-gen2thin-pdg{pdg}-minervame{tag}.root')
    pdgs = {'nue': 12, 'nuebar': -12, 'numu': 14, 'numubar': -14}
    pattern = next((pat for pat in patterns if all(
        flux_dir / pat.format(playlist=s['playlist'], tag=s['playlist'][2:],
                              flavor=s['flavor'], pdg=pdgs[s['flavor']])
        == Path(s['path']) for s in sources)), None)
    if pattern is None:
        raise RuntimeError('Unrecognized recorded flux filenames; no inputs guessed')
    for s in sources:
        if not Path(s['path']).is_file():
            raise FileNotFoundError(s['path'])
    expected = {(pl, fl) for pl in exposures for fl in pdgs}
    actual = [(s['playlist'], s['flavor']) for s in sources]
    if set(actual) != expected or len(actual) != len(expected):
        raise RuntimeError('Reference does not contain exactly four species per playlist')
    script = Path(__file__).with_name('make_me_nue_elastic_prediction.py')
    cmd = [sys.executable, '-u', str(script), '--beam', beam, '--xsec', 'lo',
           '--sin2thetaW', str(a.sin2thetaW), '--flux-dir', str(flux_dir),
           '--filename-pattern', pattern,
           '--target-electrons', str(recs[0]['target_electrons']),
           '--flux-area-unit', recs[0]['input_flux_area_unit'],
           '--hist-name', recs[0].get('flux_hist_name', 'flux_E_cvweighted'),
           '--flux-band-name', recs[0].get('flux_band_name', 'Flux'),
           '--n-universes', str(recs[0].get('n_universes', 1000)),
           '--enu-order', str(recs[0].get('enu_order', 8)),
           '--electron-order', str(recs[0].get('electron_order', 32)),
           '--check-quadrature', '--output', str(a.output)]
    for pl, pot in exposures.items():
        cmd.extend(['--playlist-pot', f'{pl}={pot}'])
    print('Reference:', a.reference, flush=True)
    print('Target electrons:', recs[0]['target_electrons'], flush=True)
    print('Exposures:', exposures, flush=True)
    print('Running:', shlex.join(cmd), flush=True)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with a.output.with_suffix('.log').open('w') as log:
        log.write(shlex.join(cmd) + '\n'); log.flush()
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in proc.stdout:
            print(line, end='', flush=True); log.write(line); log.flush()
        code = proc.wait()
    if code:
        raise RuntimeError(f'Prediction exited with {code}; inspect {a.output.with_suffix(".log")}')
    out = ROOT.TFile.Open(str(a.output), 'READ')
    if not out or out.IsZombie() or out.TestBit(ROOT.TFile.kRecovered):
        raise RuntimeError('Invalid LO output')
    values = {}
    rows = []
    edges = None
    for flavor in ('nue', 'nuebar', 'numu', 'numubar', 'total'):
        old = ref.Get('h_nue_elastic_' + flavor)
        new = out.Get('h_nue_elastic_' + flavor)
        if not old or not new or old.GetNbinsX() != new.GetNbinsX():
            raise RuntimeError('Missing/mismatched flavor histogram: ' + flavor)
        n = old.GetNbinsX()
        e = np.array([old.GetBinLowEdge(i) for i in range(1, n+2)])
        np.testing.assert_allclose(e, [new.GetBinLowEdge(i) for i in range(1, n+2)], rtol=0, atol=1e-10)
        np.testing.assert_allclose(e, recs[0]['true_edges_GeV'], rtol=0, atol=1e-10)
        edges = e
        v = np.array([new.GetBinContent(i) for i in range(1, n+1)])
        b = np.array([old.GetBinContent(i) for i in range(1, n+1)])
        values[flavor] = (v, b)
        print(f'{flavor:8s} LO={v.sum():.8f}  reference={b.sum():.8f}', flush=True)
        for i in range(n):
            rows.append([flavor, e[i], e[i+1], v[i], b[i], b[i]/v[i] if v[i] else ''])
    for index in (0, 1):
        np.testing.assert_allclose(values['total'][index],
            sum(values[f][index] for f in pdgs), rtol=1e-10, atol=1e-9)
    with a.output.with_name(a.output.stem+'_LO_vs_reference.csv').open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['flavor', 'E_low_GeV', 'E_high_GeV', 'LO_events', 'reference_events', 'reference_over_LO'])
        w.writerows(rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter

    display_edges = np.array([.8, 2, 3, 5, 7, 9, 20.])
    np.testing.assert_allclose(edges, [.8, 2, 3, 5, 7, 9, 100.])

    lo, baseline = values['total']
    scale = 2 / np.diff(display_edges)

    fig, (ax, ratio) = plt.subplots(
        2, 1, figsize=(7, 7), sharex=True,
        gridspec_kw={'height_ratios': [3, 1]},
        constrained_layout=True,
    )

    ax.stairs(
        lo * scale, display_edges,
        label='Tree level (all four species)', color='tab:blue',
    )
    ax.stairs(
        baseline * scale, display_edges,
        label='Existing Tomalak-based prediction', color='tab:red',
    )
    ax.set(title='ME ' + beam, ylabel='Events / 2 GeV')
    ax.legend(fontsize=9)
    ax.text(
        .02, .02,
        'Full 9–100 GeV yield displayed in 9–20 GeV bin',
        transform=ax.transAxes, fontsize=8,
    )

    ratio.stairs(
        np.divide(
            baseline, lo,
            out=np.full_like(lo, np.nan), where=lo > 0,
        ),
        display_edges,
        baseline=None,
        color='tab:blue',
    )
    ratio.axhline(1, color='gray', linestyle='--')
    ratio.set(
        xlabel='True total electron energy (GeV)',
        ylabel='Reference / LO',
        xlim=(0, 20),
        ylim=(0.95, 1.05),
    )
    ratio.set_yticks([0.95, 0.975, 1.00, 1.025, 1.05])
    ratio.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ratio.grid(axis='y', alpha=.25)

    fig.savefig(
        a.output.with_name(a.output.stem + '_LO_vs_reference.png'),
        dpi=180,
    )
    plt.close(fig)
    out.Close()
    ref.Close()


if __name__ == '__main__':
    main()

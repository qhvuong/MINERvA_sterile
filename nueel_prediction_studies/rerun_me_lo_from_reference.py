#!/usr/bin/env python3
"""Recompute LO using exposure and flux inputs recorded in an existing prediction."""
import argparse
import csv
import json
import math
from pathlib import Path
import shlex
import subprocess
import sys


def records(meta):
    if 'exposures' in meta:
        return [meta]
    return [r for item in meta.get('inputs', [])
            for r in records(item.get('metadata', item))]


def run_prediction(cmd, output):
    print('Running:', shlex.join(cmd), flush=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.with_suffix('.log').open('w') as log:
        log.write(shlex.join(cmd) + '\n'); log.flush()
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in proc.stdout:
            print(line, end='', flush=True); log.write(line); log.flush()
        code = proc.wait()
    if code:
        raise RuntimeError(f'Prediction exited with {code}; inspect {output.with_suffix(".log")}')


def read_prediction(ROOT, np, path, band_name, nuniv):
    f = ROOT.TFile.Open(str(path), 'READ')
    if not f or f.IsZombie() or f.TestBit(ROOT.TFile.kRecovered):
        raise RuntimeError('Invalid prediction: ' + str(path))
    result = {}
    try:
        obj = f.Get('prediction_metadata_json')
        if not obj:
            raise RuntimeError('Missing metadata: ' + str(path))
        result['metadata'] = json.loads(str(obj.GetString().Data()))
        for flavor in ('nue','nuebar','numu','numubar','total'):
            h = f.Get('h_nue_elastic_' + flavor)
            if not h or not hasattr(h, 'HasVertErrorBand') or not h.HasVertErrorBand(band_name):
                raise RuntimeError(f'{path}: missing {flavor} histogram/Flux band')
            band = h.GetVertErrorBand(band_name)
            if band.GetNHists() != nuniv:
                raise RuntimeError(f'{path}: wrong universe count for {flavor}')
            n = h.GetNbinsX()
            cv = np.array([h.GetBinContent(i) for i in range(n+2)])
            universes = np.array([[band.GetHist(u).GetBinContent(i) for i in range(n+2)] for u in range(nuniv)])
            if not np.isfinite(cv).all() or not np.isfinite(universes).all():
                raise RuntimeError(f'{path}: nonfinite {flavor} prediction')
            np.testing.assert_allclose(cv,[band.GetBinContent(i) for i in range(n+2)],rtol=1e-10,atol=1e-10)
            result[flavor] = dict(edges=np.array([h.GetBinLowEdge(i) for i in range(1,n+2)]),
                                  cv=cv,universes=universes)
        for key in ('cv','universes'):
            np.testing.assert_allclose(result['total'][key],
                sum(result[f][key] for f in ('nue','nuebar','numu','numubar')),rtol=1e-10,atol=1e-9)
    finally:
        f.Close()
    return result


def verify_flux_sum(np, combined, singles):
    checks = {}
    for flavor in ('nue','nuebar','numu','numubar','total'):
        for single in singles:
            np.testing.assert_allclose(single[flavor]['edges'],combined[flavor]['edges'],rtol=0,atol=1e-10)
        checks[flavor] = {}
        for key in ('cv','universes'):
            summed = sum(single[flavor][key] for single in singles)
            np.testing.assert_allclose(summed,combined[flavor][key],rtol=1e-10,atol=1e-8)
            checks[flavor][key+'_max_absolute_difference'] = float(np.max(np.abs(summed-combined[flavor][key])))
    return checks


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--sin2thetaW', type=float, default=0.23129)
    p.add_argument('--single-flux-dir', type=Path,
                   help='Default: sibling single_flux directory beside combined_ME, or OUTPUT parent/single_flux')
    p.add_argument('--reuse-combined', action='store_true',
                   help='Keep an existing LO combined output; produce single-flux files and verify their sum against it')
    p.add_argument('--no-plots', action='store_true',help='Skip Matplotlib; ROOT and CSV outputs are still produced')
    a = p.parse_args()
    if a.reference.resolve() == a.output.resolve():
        p.error('Reference and LO output must be different files.')
    if a.output.exists() and not a.reuse_combined:
        p.error('LO output exists; use --reuse-combined to retain it, or choose a new output.')
    if a.reuse_combined and not a.output.is_file():
        p.error('--reuse-combined requires an existing LO output file.')
    if not math.isfinite(a.sin2thetaW) or not 0 < a.sin2thetaW < 1:
        p.error('Require 0 < sin2thetaW < 1.')
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
    single_dir = a.single_flux_dir or (
        a.output.parent.parent/'single_flux' if a.output.parent.name=='combined_ME'
        else a.output.parent/'single_flux')
    single_paths = {pl:single_dir/f'me_nue_elastic_{pl}.root' for pl in exposures}
    for path in single_paths.values():
        if path.resolve() in (a.output.resolve(),a.reference.resolve()):
            p.error('Single-flux output collides with reference/combined output: ' + str(path))
        for target in (path,path.with_suffix('.json'),path.with_suffix('.csv'),path.with_suffix('.log'),
                       path.with_name(path.stem+'_paper_covariance_12.csv')):
            if target.exists():
                p.error('Single-flux output already exists; choose a new --single-flux-dir: ' + str(target))
    band_name = recs[0].get('flux_band_name','Flux')
    nuniv = int(recs[0].get('n_universes',1000))
    if a.reuse_combined:
        existing = read_prediction(ROOT,np,a.output,band_name,nuniv)
        saved = existing['metadata']
        if saved.get('xsec_model') != 'lo' or saved.get('beam') != beam:
            raise RuntimeError('Existing combined output is not the expected LO beam prediction')
        if saved.get('exposures') != exposures or saved.get('sin2thetaW') != a.sin2thetaW:
            raise RuntimeError('Existing combined exposure/weak-angle settings differ')
        for key in shared:
            if saved.get(key) != recs[0][key]:
                raise RuntimeError('Existing combined setting differs: ' + key)
        recorded = {(x['playlist'],x['flavor'],str(Path(x['path']).resolve())) for x in saved.get('inputs',[])}
        expected_inputs = {(x['playlist'],x['flavor'],str(Path(x['path']).resolve())) for x in sources}
        if recorded != expected_inputs:
            raise RuntimeError('Existing combined flux source paths differ from reference')
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
           '--check-quadrature']
    print('Reference:', a.reference, flush=True)
    print('Target electrons:', recs[0]['target_electrons'], flush=True)
    print('Exposures:', exposures, flush=True)
    print('Single-flux outputs:',single_dir,flush=True)
    for pl,pot in exposures.items():
        single_cmd = cmd + ['--playlist',pl,'--pot',str(pot),'--output',str(single_paths[pl])]
        run_prediction(single_cmd,single_paths[pl])
    if not a.reuse_combined:
        combined_cmd = cmd + ['--output',str(a.output)]
        for pl,pot in exposures.items():
            combined_cmd.extend(['--playlist-pot',f'{pl}={pot}'])
        run_prediction(combined_cmd,a.output)
    combined = read_prediction(ROOT,np,a.output,band_name,nuniv)
    singles = [read_prediction(ROOT,np,path,band_name,nuniv) for path in single_paths.values()]
    checks = verify_flux_sum(np,combined,singles)
    closure_path = a.output.with_name(a.output.stem+'_single_flux_closure.json')
    closure_path.write_text(json.dumps(dict(status='PASS',combined=str(a.output),
        single_flux_outputs={pl:str(path) for pl,path in single_paths.items()},
        exposures=exposures,n_universes=nuniv,checks=checks),indent=2))
    print('PASS: single-flux CVs and all universes sum to the combined prediction for every flavor.',flush=True)
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
    if a.no_plots:
        out.Close();ref.Close()
        print('Numerical outputs complete; plotting skipped.',flush=True)
        return
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

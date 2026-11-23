#!/usr/bin/env python3
"""Read-only installed ME flux / CV-weight audit. Run in the GPVM PyROOT setup.

Writes reports only. Does not change ROOT inputs, fluxes, or selection settings.
The runtime factory check uses the installed MAT library, not a replacement
implementation. Its FRW filename messages are retained in per-playlist logs.
This does not inspect the historical grid tarball or reweight selected events.
"""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import sys

BASE = Path('/exp/minerva/app/users/qvuong/MAT_AL9')
SPECIES = {14: 'numu', -14: 'numubar', 12: 'nue', -12: 'nuebar'}
GROUPS = {
    '1D': ['1A', '1B', '1C', '1D', '1E', '1F'],
    '1M': ['1G', '1L', '1M'], '1N': ['1N', '1O', '1P'],
    '5A': ['5A'],
    '6A': ['6A', '6B', '6C', '6D', '6E', '6F', '6G', '6H', '6I', '6J'],
}
ENERGIES = [round(0.8 + 0.1*i, 8) for i in range(193)]


def root_setup():
    import ROOT
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(True)
    import PlotUtils  # noqa: F401; load the user's normal dictionaries
    ROOT.TH1.AddDirectory(False)
    return ROOT


def read_hist(ROOT, path, name):
    f = ROOT.TFile.Open(str(path), 'READ')
    if not f or f.IsZombie():
        raise RuntimeError('Cannot open ' + str(path))
    try:
        if f.TestBit(ROOT.TFile.kRecovered):
            raise RuntimeError('Recovered file: ' + str(path))
        h = f.Get(name)
        if not h:
            raise RuntimeError(f'Missing {name} in {path}')
        h = h.Clone()
        h.SetDirectory(0)
        return h
    finally:
        f.Close()


def edges(h):
    return [float(h.GetBinLowEdge(i)) for i in range(1, h.GetNbinsX()+2)]


def equal_cv(a, b):
    return edges(a) == edges(b) and all(
        math.isclose(a.GetBinContent(i), b.GetBinContent(i), rel_tol=1e-10, abs_tol=0)
        for i in range(a.GetNbinsX()+2))


def ratio(n, d):
    return n/d if math.isfinite(n) and math.isfinite(d) and d != 0 else None


def weight(n, d, e):
    # MAT GetFluxCVWeight: interpolate each flux, then divide. It returns
    # unity for zero numerator/denominator; expose such bins separately.
    nv, dv = float(n.Interpolate(e)), float(d.Interpolate(e))
    return (1.0 if nv == 0 or dv == 0 else ratio(nv, dv)), nv, dv


def runtime_child(playlist, out):
    ROOT = root_setup()
    frw = ROOT.PlotUtils.flux_reweighter('minervame'+playlist, 14, False, 1000)
    result = {'playlist': playlist, 'libraries': str(ROOT.gSystem.GetLibraries()), 'rows': []}
    # A muon factory holds both signs. This is the runtime path relevant to
    # CCnumu(bar); the direct file check below includes all four species.
    for pdg in (14, -14):
        num, den = frw.GetFluxReweighted(pdg), frw.GetFluxGenerated(pdg)
        result.setdefault('histograms', {})[str(pdg)] = {
            'weighted_edges': edges(num), 'generated_edges': edges(den)}
        for e in ENERGIES:
            result['rows'].append({'pdg': pdg, 'energy': e,
                'weight': float(frw.GetFluxCVWeight(e, pdg)),
                'weighted': float(num.Interpolate(e)),
                'generated': float(den.Interpolate(e))})
    Path(out).write_text(json.dumps(result, indent=2, allow_nan=False))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--custom', type=Path, default=BASE/'CC-NuE-XSec/custom_plotutils/data/flux')
    p.add_argument('--legacy', type=Path, default=BASE/'opt/lib/data/flux')
    p.add_argument('--produced', type=Path, default=Path('/exp/minerva/data/users/qvuong/flux_studies/producedFluxes_p8_ME_offsetBeamFocus'))
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--skip-runtime', action='store_true')
    p.add_argument('--runtime-child', nargs=2, metavar=('PLAYLIST', 'JSON'))
    a = p.parse_args()
    if a.runtime_child:
        runtime_child(*a.runtime_child)
        return
    a.output_dir.mkdir(parents=True, exist_ok=True)
    ROOT = root_setup()
    report = {'custom': str(a.custom), 'legacy': str(a.legacy),
        'produced': str(a.produced), 'PLOTUTILSROOT': os.environ.get('PLOTUTILSROOT'),
        'libraries': str(ROOT.gSystem.GetLibraries()), 'files': [], 'issues': [],
        'scope': 'Current installed files/library; historical grid payload not verified.'}
    rows, installed = [], {}
    for tag in GROUPS:
        for pdg, species in SPECIES.items():
            try:
                hs = {}
                for label, directory in [('custom', a.custom), ('legacy', a.legacy)]:
                    for version, name in [('gen2thin', 'flux_E_cvweighted'), ('g4numiv6', 'flux_E_unweighted')]:
                        path = directory/f'flux-{version}-pdg{pdg}-minervame{tag}.root'
                        h = read_hist(ROOT, path, name)
                        hs[label, version] = h
                        n_univ = int(h.GetVertErrorBand('Flux').GetNHists()) if h.HasVertErrorBand('Flux') else 0
                        entry = {'path': str(path), 'histogram': name, 'flux_universes': n_univ}
                        if label == 'custom':
                            installed[tag, pdg, version] = h
                            if n_univ != 1000:
                                report['issues'].append(f'{path}: Flux universes={n_univ}')
                            for i in range(h.GetNbinsX()+2):
                                if not math.isfinite(h.GetBinContent(i)) or h.GetBinContent(i) < 0:
                                    report['issues'].append(f'{path}: invalid CV bin {i}')
                            source = a.produced/f'ME{tag}_{species}.root'
                            try:
                                sh = read_hist(ROOT, source, name)
                                entry['source_cv_matches'] = equal_cv(h, sh)
                                if not entry['source_cv_matches']:
                                    report['issues'].append(f'{path}: CV differs from {source}')
                            except Exception as error:
                                entry['source_error'] = str(error)
                                report['issues'].append(str(error))
                        report['files'].append(entry)
                for label in ('custom', 'legacy'):
                    if edges(hs[label, 'gen2thin']) != edges(hs[label, 'g4numiv6']):
                        report['issues'].append(f'{tag} {pdg}: {label} numerator/denominator binning differs')
                for e in ENERGIES:
                    cw, cn, cd = weight(hs['custom', 'gen2thin'], hs['custom', 'g4numiv6'], e)
                    lw, ln, ld = weight(hs['legacy', 'gen2thin'], hs['legacy', 'g4numiv6'], e)
                    rows.append({'flux': tag, 'pdg': pdg, 'energy': e,
                        'custom_weight': cw, 'legacy_weight': lw,
                        'weight_custom_over_legacy': ratio(cw, lw) if cw is not None and lw is not None else None,
                        'weighted_custom_over_legacy': ratio(cn, ln),
                        'generated_custom_over_legacy': ratio(cd, ld),
                        'zero_flux_fallback': cn == 0 or cd == 0 or ln == 0 or ld == 0})
                right_sign = (pdg == 14 if tag.startswith('1') else pdg == -14)
                if right_sign:
                    peak = [r['weight_custom_over_legacy'] for r in rows
                            if r['flux'] == tag and r['pdg'] == pdg and 2 <= r['energy'] <= 8
                            and r['weight_custom_over_legacy'] is not None and not r['zero_flux_fallback']]
                    print(f'{tag} PDG {pdg}: CV weight custom/legacy range at 2–8 GeV: '
                          + (f'{min(peak):.6f} to {max(peak):.6f}' if peak else 'UNAVAILABLE'), flush=True)
            except Exception as error:
                report['issues'].append(f'{tag} PDG {pdg}: {error}')
    if rows:
        with (a.output_dir/'cv_weight_comparison.csv').open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    if not a.skip_runtime:
        for representative, playlists in GROUPS.items():
            for playlist in playlists:
                output = (a.output_dir/f'runtime_{playlist}.json').resolve()
                log = a.output_dir/f'runtime_{playlist}.log'
                command = [sys.executable, '-u', str(Path(__file__).resolve()),
                           '--output-dir', str(a.output_dir), '--runtime-child', playlist, str(output)]
                with log.open('w') as stream:
                    try:
                        run = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, timeout=180)
                        if run.returncode != 0:
                            raise RuntimeError(f'child exited {run.returncode}; see {log}')
                        result = json.loads(output.read_text())
                        bad = 0
                        for r in result['rows']:
                            n = installed[representative, r['pdg'], 'gen2thin']
                            d = installed[representative, r['pdg'], 'g4numiv6']
                            expected, nv, dv = weight(n, d, r['energy'])
                            # Check actual runtime histograms as well as their ratio:
                            # a common wrong normalization can cancel in the ratio.
                            if expected is None or not all(math.isclose(x, y, rel_tol=1e-8, abs_tol=0)
                                for x, y in [(r['weight'], expected), (r['weighted'], nv), (r['generated'], dv)]):
                                bad += 1
                        print(f'Runtime {playlist} -> {representative}: {bad} mismatches', flush=True)
                        if bad:
                            report['issues'].append(f'Runtime {playlist}: {bad} mismatches; see {output} and {log}')
                    except Exception as error:
                        report['issues'].append(f'Runtime {playlist}: {error}')
    else:
        report['runtime'] = 'NOT CHECKED (--skip-runtime)'
    report['status'] = 'CHECK' if report['issues'] else 'CHECKS_MATCH'
    (a.output_dir/'audit.json').write_text(json.dumps(report, indent=2, allow_nan=False))
    print('\nStatus:', report['status'])
    for issue in report['issues']:
        print('[CHECK]', issue)
    print('Reports:', a.output_dir)


if __name__ == '__main__':
    main()

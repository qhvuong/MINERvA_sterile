#!/usr/bin/env python3
"""Compare ME nu-e predictions made with produced and legacy fluxes.

Accepts combined beam ROOT files or the five exposure-normalized playlist
outputs. Sums matching universes within each flux set; never pairs universes
between different flux productions. Does not modify prediction ROOT files.
"""
import argparse
import csv
import json
from pathlib import Path

GROUPS = {'FHC': ('ME1D', 'ME1M', 'ME1N'), 'RHC': ('ME5A', 'ME6A')}


def metadata_records(meta):
    """Flatten combination metadata to its per-prediction metadata records."""
    children = [v['metadata'] for v in meta.get('inputs', [])
                if isinstance(v, dict) and isinstance(v.get('metadata'), dict)]
    if children:
        return [record for child in children for record in metadata_records(child)]
    return [meta]


def load_prediction(ROOT, np, directory, beam, reco):
    suffix = '_reco' if reco else ''
    combined = directory / ('me_nue_elastic_' + beam + suffix + '.root')
    paths = [combined] if combined.is_file() else [
        directory / ('me_nue_elastic_' + p + '.root') for p in GROUPS[beam]]
    name = 'h_nue_elastic_' + ('reco_' if reco else '') + 'total'
    total = universes = edges = None
    records = []
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
        f = ROOT.TFile.Open(str(path), 'READ')
        if not f or f.IsZombie() or f.TestBit(ROOT.TFile.kRecovered):
            raise RuntimeError('Invalid ROOT input: ' + str(path))
        try:
            h = f.Get(name)
            if not h or not h.HasVertErrorBand('Flux'):
                raise RuntimeError('Missing ' + name + ' or Flux band: ' + str(path))
            band = h.GetVertErrorBand('Flux')
            n, nu = h.GetNbinsX(), band.GetNHists()
            if nu < 2:
                raise ValueError('Need at least two Flux universes: ' + str(path))
            e = np.array([h.GetBinLowEdge(i) for i in range(1, n+2)])
            cv = np.array([h.GetBinContent(i) for i in range(1, n+1)])
            u = np.array([[band.GetHist(j).GetBinContent(i) for i in range(1, n+1)]
                          for j in range(nu)])
            if not all(np.isfinite(v).all() for v in (e, cv, u)):
                raise ValueError('Nonfinite prediction: ' + str(path))
            if np.any(np.diff(e) <= 0):
                raise ValueError('Invalid bin edges: ' + str(path))
            if total is None:
                total, universes, edges = cv, u, e
            else:
                np.testing.assert_allclose(e, edges, rtol=0, atol=1e-10)
                if u.shape != universes.shape:
                    raise ValueError('Universe counts differ within ' + str(directory))
                total += cv
                universes += u
            obj = f.Get('prediction_metadata_json')
            if not obj:
                raise ValueError('Missing prediction metadata: ' + str(path))
            records += metadata_records(json.loads(str(obj.GetString().Data())))
        finally:
            f.Close()
    exposure = {}
    settings = {}
    for record in records:
        if record.get('beam') != beam:
            raise ValueError('Metadata beam mismatch in ' + str(directory))
        for p, pot in record.get('exposures', {}).items():
            if p in exposure:
                raise ValueError('Duplicate exposure for ' + p)
            exposure[p] = float(pot)
        for key in ('target_electrons', 'xsec_total_scales', 'observable', 'response_json'):
            if key in record:
                if key in settings and settings[key] != record[key]:
                    raise ValueError('Inconsistent metadata setting: ' + key)
                settings[key] = record[key]
    if not exposure:
        raise ValueError('No exposure metadata in ' + str(directory))
    return dict(cv=total, universes=universes, edges=edges,
                sigma=np.std(universes, axis=0, ddof=0),
                exposures=exposure, settings=settings, paths=list(map(str, paths)))


def verify_comparison(np, produced, legacy):
    np.testing.assert_allclose(produced['edges'], legacy['edges'], rtol=0, atol=1e-10)
    if produced['exposures'].keys() != legacy['exposures'].keys():
        raise ValueError('Different playlist exposures between flux productions')
    for p in produced['exposures']:
        np.testing.assert_allclose(produced['exposures'][p], legacy['exposures'][p], rtol=1e-10)
    if produced['settings'] != legacy['settings']:
        raise ValueError('Prediction settings differ; compare metadata before interpreting a flux-only change')
    if np.any(legacy['cv'] <= 0):
        raise ValueError('Legacy CV has nonpositive bins; cannot construct the ratio')


def draw(np, ax, ratio, beam, produced, legacy, data=None, data_err=None, reco=False):
    np.testing.assert_allclose(produced['edges'], [.8, 2, 3, 5, 7, 9, 100], atol=1e-9)
    edges = np.array([.8, 2, 3, 5, 7, 9, 20.])
    scale = 2 / np.diff(edges)
    colors = ('tab:blue', 'tab:orange')
    for pred, color, label in zip((produced, legacy), colors, ('Produced flux', 'Legacy flux')):
        cv, sigma = pred['cv'], pred['sigma']
        ax.stairs((cv-sigma)*scale, edges, baseline=(cv+sigma)*scale,
                  fill=True, color=color, alpha=.22, linewidth=0,
                  label=label + r': $\pm1\sigma$ flux')
        ax.stairs(cv*scale, edges, color=color, linewidth=1.8, label=label + ': CV')
    if data is not None:
        ax.errorbar((edges[:-1]+edges[1:])/2, data*scale,
                    xerr=np.diff(edges)/2, yerr=data_err*scale,
                    fmt='ko', markersize=4, linewidth=1, label='Published data')
    ax.set(title='ME ' + beam, ylabel='Events / 2 GeV')
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=9, loc='upper right')
    ax.text(.97, .69, 'Response-folded' if reco else 'True-energy prediction (diagnostic)',
            transform=ax.transAxes, ha='right', fontsize=8)
    difference = 100. * (produced['cv'] / legacy['cv'] - 1.)
    ratio.stairs(difference, edges, color=colors[0], linewidth=1.8)
    limit = max(.1, float(np.max(np.abs(difference))) * 1.3)
    ratio.set(xlabel='Electron energy (GeV)', ylabel='CV difference (%)',
              xlim=(0,20), ylim=(-limit, limit))
    ratio.set_xticks(np.arange(0,21,2))
    ratio.axhline(0, color='0.4', lw=.7, linestyle='--')
    ratio.text(.98, .04, r'Last bin: all $E_e>9$ GeV',
               transform=ratio.transAxes, ha='right', fontsize=8)
    for panel in (ax,ratio):
        panel.minorticks_on()
        panel.tick_params(which='both', direction='in', top=True, right=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--produced-dir', type=Path, required=True)
    ap.add_argument('--legacy-dir', type=Path, default=Path('legacy_ME'))
    ap.add_argument('--output-dir', type=Path, default=Path('compare_ME_flux_predictions'))
    ap.add_argument('--paper-inputs', type=Path, default=Path('me_nue_elastic_paper_inputs.json'))
    ap.add_argument('--no-data', action='store_true')
    ap.add_argument('--reco', action='store_true')
    args = ap.parse_args()
    import ROOT
    import PlotUtils.LoadPlotUtilsLib
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    paper = None if args.no_data else json.loads(args.paper_inputs.read_text())
    if paper:
        data = np.asarray(paper['data_13'], float)[:12]
        errors = np.sqrt(np.diag(np.asarray(paper['covariance_13'], float))[:12])
    fig, axes = plt.subplots(2,2,figsize=(13,8),sharex=True,
        gridspec_kw={'height_ratios':[3,1]},constrained_layout=True)
    summary = {'flux_band': 'population standard deviation about the universe mean; drawn around CV',
               'lower_panel': '100 * (produced CV / legacy CV - 1); flux bands shown only in upper panels',
               'beams': {}}
    for col, beam in enumerate(GROUPS):
        produced = load_prediction(ROOT,np,args.produced_dir,beam,args.reco)
        legacy = load_prediction(ROOT,np,args.legacy_dir,beam,args.reco)
        verify_comparison(np,produced,legacy)
        observed = None if paper is None else data[6*col:6*(col+1)]
        err = None if paper is None else errors[6*col:6*(col+1)]
        draw(np,axes[0,col],axes[1,col],beam,produced,legacy,observed,err,args.reco)
        single, panels = plt.subplots(2,1,figsize=(7,7),sharex=True,
            gridspec_kw={'height_ratios':[3,1]},constrained_layout=True)
        draw(np,panels[0],panels[1],beam,produced,legacy,observed,err,args.reco)
        stem = args.output_dir / ('me_nueel_produced_vs_legacy_' + beam)
        for ext in ('png','pdf'):
            single.savefig(str(stem)+'.'+ext,dpi=180)
        plt.close(single)
        with Path(str(stem)+'.csv').open('w',newline='') as f:
            writer=csv.writer(f)
            writer.writerow(['E_low_GeV','E_high_GeV','produced_CV','produced_flux_sigma',
                             'legacy_CV','legacy_flux_sigma','produced_over_legacy_CV'])
            writer.writerows(zip(produced['edges'][:-1],produced['edges'][1:],produced['cv'],
                produced['sigma'],legacy['cv'],legacy['sigma'],produced['cv']/legacy['cv']))
        info = dict(produced_total=float(produced['cv'].sum()),legacy_total=float(legacy['cv'].sum()),
            total_ratio=float(produced['cv'].sum()/legacy['cv'].sum()),
            produced_total_flux_sigma=float(np.std(produced['universes'].sum(axis=1))),
            legacy_total_flux_sigma=float(np.std(legacy['universes'].sum(axis=1))),
            produced_universes=len(produced['universes']),legacy_universes=len(legacy['universes']),
            exposures=produced['exposures'],settings=produced['settings'],
            produced_inputs=produced['paths'],legacy_inputs=legacy['paths'])
        summary['beams'][beam]=info
        print(f"{beam}: produced={info['produced_total']:.6f}, legacy={info['legacy_total']:.6f}, ratio={info['total_ratio']:.6f}")
    fig.suptitle('ME neutrino–electron prediction: produced versus legacy flux')
    stem = args.output_dir / 'me_nueel_produced_vs_legacy_combined'
    for ext in ('png','pdf'):
        fig.savefig(str(stem)+'.'+ext,dpi=180)
    plt.close(fig)
    (args.output_dir/'comparison_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('Saved plots, bin counts, flux uncertainties, and summary in',args.output_dir)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Sum exposure-normalized ME predictions, preserving paired Flux universes.

Run in nueel_prediction_studies. Default: diagnostic truth/reco overlay.
Use --reco only when every input contains response-folded histograms.
"""
import argparse
from array import array
import json
from pathlib import Path

GROUPS = {'FHC': ('ME1D', 'ME1M', 'ME1N'), 'RHC': ('ME5A', 'ME6A')}
FLAVORS = ('total', 'nue', 'nuebar', 'numu', 'numubar')


def draw_paper_style(ax, ratio, beam, pred, std, observed, err, reco=False):
    """Figure-6 display: full >9 GeV count represented by the 9--20 bin."""
    import numpy as np
    edges = np.array([.8, 2., 3., 5., 7., 9., 20.])
    widths = np.diff(edges)
    centers = (edges[:-1] + edges[1:]) / 2
    scale = 2. / widths
    y, sy = pred * scale, std * scale
    if np.any(pred <= 0):
        raise ValueError('Cannot form data/prediction with nonpositive prediction bins')
    ax.stairs(y, edges, color='red', linewidth=1.8, label='Prediction')
    ax.errorbar(centers, y, xerr=widths/2, yerr=sy, fmt='none',
                ecolor='red', elinewidth=1.2)
    ax.errorbar(centers, observed*scale, xerr=widths/2, yerr=err*scale,
                fmt='ko', markersize=4, elinewidth=1.2, label='Data')
    ratio.errorbar(centers, observed/pred, xerr=widths/2, yerr=err/pred,
                   fmt='ko', markersize=4, elinewidth=1.2)
    # Separate prediction uncertainty avoids assuming independence from the
    # measurement's small beam-related correction uncertainties.
    ratio.fill_between(edges, np.r_[1-std/pred, (1-std/pred)[-1]],
        np.r_[1+std/pred, (1+std/pred)[-1]], step='post', color='red', alpha=.14)
    ratio.axhline(1, color='red', linewidth=1)
    ratio.set(xlabel='Electron Energy (GeV)', ylabel='Data / Prediction',
              xlim=(0, 20), ylim=(0, 1.8))
    ratio.set_xticks(np.arange(0, 21, 2))
    ax.set(ylabel='N Events / 2.0 GeV', title='ME ' + beam, ylim=(0, None))
    ax.set_ylim(0, 1.25*max(np.max(y+sy), np.max((observed+err)*scale)))
    ax.legend(frameon=False, fontsize=10)
    ax.text(.97, .75, 'Response-folded prediction' if reco else 'Truth prediction: diagnostic overlay',
            transform=ax.transAxes, ha='right', fontsize=8)
    for panel in (ax, ratio):
        panel.minorticks_on()
        panel.tick_params(which='both', direction='in', top=True, right=True)
        panel.tick_params(which='major', length=6)
        panel.tick_params(which='minor', length=3)
    ratio.text(.98, .04, 'Last bin: all $E_e>9$ GeV', transform=ratio.transAxes,
               ha='right', fontsize=8)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input-dir', type=Path, default=Path('.'))
    ap.add_argument('--output-dir', type=Path, default=Path('combined_ME'))
    ap.add_argument('--paper-inputs', type=Path, default=Path('me_nue_elastic_paper_inputs.json'))
    ap.add_argument('--reco', action='store_true', help='Use existing h_nue_elastic_reco_* inputs')
    args = ap.parse_args()
    import numpy as np
    import ROOT
    import PlotUtils.LoadPlotUtilsLib
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    paper = json.loads(args.paper_inputs.read_text())
    data = np.asarray(paper['data_13'], float)[:12]
    cov = np.asarray(paper['covariance_13'], float)[:12, :12]
    if data.shape != (12,) or cov.shape != (12, 12):
        raise ValueError('Expected twelve published bins and their joint covariance')
    np.testing.assert_allclose(cov, cov.T, atol=1e-10)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    prefix = 'h_nue_elastic_reco_' if args.reco else 'h_nue_elastic_'
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True,
                             gridspec_kw={'height_ratios': [3, 1]}, layout='constrained')
    summary = {}
    for col, (beam, groups) in enumerate(GROUPS.items()):
        sums = {}
        edges = None
        metadata = []
        for group in groups:
            path = args.input_dir / ('me_nue_elastic_' + group + '.root')
            f = ROOT.TFile.Open(str(path), 'READ')
            if not f or f.IsZombie() or f.TestBit(ROOT.TFile.kRecovered):
                raise RuntimeError('Invalid input: ' + str(path))
            try:
                obj = f.Get('prediction_metadata_json')
                meta = json.loads(str(obj.GetString().Data())) if obj else {}
                if meta.get('beam', beam) != beam:
                    raise ValueError('Wrong beam in ' + str(path))
                metadata.append({'file': str(path), 'metadata': meta})
                for flavor in FLAVORS:
                    name = prefix + flavor
                    h = f.Get(name)
                    if not h or not h.HasVertErrorBand('Flux'):
                        raise RuntimeError('Missing histogram/Flux band: ' + str(path) + ':' + name)
                    n = h.GetNbinsX()
                    e = np.array([h.GetBinLowEdge(i) for i in range(1, n + 2)])
                    if edges is None:
                        edges = e
                    np.testing.assert_allclose(e, edges, rtol=0, atol=1e-10)
                    band = h.GetVertErrorBand('Flux')
                    if band.GetNHists() != 1000:
                        raise ValueError('Expected 1000 Flux universes: ' + str(path))
                    cv = np.array([h.GetBinContent(i) for i in range(n + 2)])
                    u = np.array([[band.GetHist(j).GetBinContent(i) for i in range(n + 2)]
                                  for j in range(1000)])
                    if not np.isfinite(cv).all() or not np.isfinite(u).all():
                        raise ValueError('Nonfinite prediction: ' + str(path))
                    if flavor not in sums:
                        sums[flavor] = [cv.copy(), u.copy()]
                    else:
                        sums[flavor][0] += cv
                        sums[flavor][1] += u
            finally:
                f.Close()
        np.testing.assert_allclose(edges, paper['numerical_edges_GeV'], rtol=0, atol=1e-10)
        for component in (0, 1):
            np.testing.assert_allclose(sums['total'][component],
                sum(sums[k][component] for k in FLAVORS[1:]), rtol=1e-10, atol=1e-9)
        out = args.output_dir / ('me_nue_elastic_' + beam + ('_reco' if args.reco else '') + '.root')
        fout = ROOT.TFile(str(out), 'RECREATE')
        if not fout or fout.IsZombie():
            raise RuntimeError('Cannot create ' + str(out))
        keep = []
        for flavor, (cv, u) in sums.items():
            name = prefix + flavor
            h = ROOT.PlotUtils.MnvH1D(name, name, len(edges)-1, array('d', edges))
            h.SetDirectory(0)
            h.AddVertErrorBand('Flux', 1000)
            band = h.GetVertErrorBand('Flux')
            for i in range(len(cv)):
                h.SetBinContent(i, float(cv[i])); h.SetBinError(i, 0.)
                band.SetBinContent(i, float(cv[i])); band.SetBinError(i, 0.)
                for j in range(1000):
                    hu = band.GetHist(j)
                    hu.SetBinContent(i, float(u[j, i])); hu.SetBinError(i, 0.)
            fout.cd(); h.Write(); keep.append(h)
        # Published data/covariance are copied once, never summed over playlists.
        for b, start in [('FHC', 0), ('RHC', 6)]:
            h = ROOT.TH1D('paper_data_' + b, 'Published measurement', 6, array('d', edges))
            h.SetDirectory(0)
            for i in range(6):
                h.SetBinContent(i+1, float(data[start+i]))
                h.SetBinError(i+1, float(np.sqrt(cov[start+i, start+i])))
            fout.cd(); h.Write(); keep.append(h)
        for name in ('paper_measurement_covariance_12', 'cov_nonflux_ME_FHC_RHC'):
            mat = ROOT.TMatrixD(12, 12)
            for i in range(12):
                for j in range(12):
                    mat[i][j] = float(cov[i, j])
            fout.cd(); mat.Write(name)
        meta = {'beam': beam, 'inputs': metadata, 'operation': 'sum; no additional exposure scaling',
                'observable': 'reco' if args.reco else 'true electron energy',
                'universe_pairing': 'same index summed across playlists and flavors',
                'paper_inputs': str(args.paper_inputs), 'flux_std_reference': 'ensemble mean',
                'note': 'Truth/reco overlay is diagnostic; no response applied by this script.'}
        fout.cd(); ROOT.TObjString(json.dumps(meta)).Write('prediction_metadata_json')
        fout.Close()
        cv, u = sums['total']
        pred = cv[1:-1]
        std = np.std(u[:, 1:-1], axis=0, ddof=0)
        observed = data[col*6:(col+1)*6]
        err = np.sqrt(np.diag(cov)[col*6:(col+1)*6])
        ax, ratio = axes[0, col], axes[1, col]
        draw_paper_style(ax, ratio, beam, pred, std, observed, err, args.reco)
        single, panels = plt.subplots(2, 1, figsize=(7, 7), sharex=True,
            gridspec_kw={'height_ratios': [3, 1], 'hspace': .03})
        single.subplots_adjust(left=.15, right=.97, bottom=.11, top=.94)
        draw_paper_style(panels[0], panels[1], beam, pred, std, observed, err, args.reco)
        single_stem = args.output_dir / ('me_nue_elastic_' + beam + '_figure6_style')
        single.savefig(str(single_stem)+'.png', dpi=180)
        single.savefig(str(single_stem)+'.pdf')
        plt.close(single)
        summary[beam] = {'prediction_total': float(pred.sum()), 'data_total': float(observed.sum()),
                         'prediction_bins': pred.tolist(), 'data_bins': observed.tolist()}
        print(f'{beam}: prediction={pred.sum():.6f}, data={observed.sum():.2f}; wrote {out}')
    fig.suptitle('ME neutrino–electron scattering: combined exposures' +
                 (' (response-folded)' if args.reco else '\nDiagnostic: true-energy prediction versus efficiency-corrected reconstructed-energy data'), fontsize=12)
    stem = args.output_dir / ('me_nue_elastic_combined' + ('_reco' if args.reco else '_truth_diagnostic'))
    fig.savefig(str(stem)+'.png', dpi=180)
    fig.savefig(str(stem)+'.pdf')
    plt.close(fig)
    stem.with_suffix('.json').write_text(json.dumps(summary, indent=2))
    print('Saved comparison:', str(stem)+'.png', 'and .pdf')


if __name__ == '__main__':
    main()

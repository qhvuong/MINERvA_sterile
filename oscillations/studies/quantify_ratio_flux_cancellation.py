#!/usr/bin/env python3
"""Measure common-flux cancellation using paired MnvH1D flux universes.

Run from oscillations/ in the usual ROOT/PlotUtils environment. No fit runs.
Centered ensemble moments use 1/N, matching the normalized-prior convention.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np


def covariance(x, y=None):
    """Population ensemble covariance, centered separately in each bin."""
    x = np.asarray(x, dtype=float)
    y = x if y is None else np.asarray(y, dtype=float)
    return (x-x.mean(axis=0)).T @ (y-y.mean(axis=0)) / len(x)


def cancellation(mu_cv, e_cv, mu_u, e_u):
    """Compare paired ratios with an exactly independent ensemble baseline.

    The independent baseline pairs each numerator universe with every
    denominator universe with equal weight. Its moments factorize, so it
    needs no random shuffle and no N-by-N array. This isolates the effect
    of pairing/correlation even for the nonlinear M/E ratio.
    """
    mu_cv, e_cv, mu_u, e_u = map(np.asarray, (mu_cv, e_cv, mu_u, e_u))
    if mu_u.shape != e_u.shape or mu_u.ndim != 2 or mu_u.shape[1] != len(mu_cv) or len(e_cv) != len(mu_cv):
        raise ValueError('Muon/electron arrays or CV bin counts do not match')
    if not all(np.isfinite(a).all() for a in (mu_cv, e_cv, mu_u, e_u)):
        raise ValueError('Non-finite CV or universe content')
    if np.any(mu_cv == 0) or np.any(e_cv == 0) or np.any(e_u == 0):
        raise ValueError('Undefined fractional shift or ratio: zero CV or universe denominator')
    dm = mu_u / mu_cv - 1.
    de = e_u / e_cv - 1.
    r_cv = mu_cv/e_cv
    dr = (mu_u/e_u)/r_cv - 1.
    v_m, v_e, cross = covariance(dm), covariance(de), covariance(dm, de)
    v_uncorr = v_m + v_e
    v_linear = covariance(dm-de)
    np.testing.assert_allclose(v_linear, v_uncorr-cross-cross.T, rtol=1e-10, atol=1e-13)
    var_m, var_e, cov_me = np.diag(v_m), np.diag(v_e), np.diag(cross)
    var_uncorr, var_linear, var_exact = np.diag(v_uncorr), np.diag(v_linear), np.diag(covariance(dr))
    # Var(XY) = Var(X) E[Y^2] + E[X]^2 Var(Y), when independent.
    # X=M/M_CV, Y=E_CV/E. Compute centered variances for stability.
    x, y = mu_u/mu_cv, e_cv/e_u
    vx = np.mean((x-x.mean(axis=0))**2, axis=0)
    vy = np.mean((y-y.mean(axis=0))**2, axis=0)
    var_ind_exact = vx*np.mean(y*y, axis=0) + np.mean(x, axis=0)**2*vy

    def divide(a, b):
        return np.divide(a, b, out=np.full_like(a, np.nan), where=b > 0)

    values = dict(
        sigma_mu_fraction=np.sqrt(var_m), sigma_e_fraction=np.sqrt(var_e),
        rho_mu_e=divide(cov_me, np.sqrt(var_m*var_e)),
        sigma_ratio_uncorrelated_linear_fraction=np.sqrt(var_uncorr),
        sigma_ratio_paired_linear_fraction=np.sqrt(np.maximum(var_linear, 0.)),
        covariance_removed_linear_fraction_sq=2*cov_me,
        variance_cancelled_linear_pct=100*divide(2*cov_me, var_uncorr),
        sigma_ratio_independent_exact_fraction=np.sqrt(var_ind_exact),
        sigma_ratio_paired_exact_fraction=np.sqrt(var_exact),
        variance_cancelled_exact_pct=100*divide(var_ind_exact-var_exact, var_ind_exact),
        sigma_reduction_exact_pct=100*(1-np.sqrt(divide(var_exact, var_ind_exact))),
        ratio_mean_shift_fraction=dr.mean(axis=0),
        ratio_rms_about_cv_fraction=np.sqrt(np.mean(dr*dr, axis=0)),
        nonlinear_difference_rms_fraction=np.sqrt(np.mean((dr-(dm-de))**2, axis=0)),
    )
    arrays = dict(mu_fractional_shifts=dm, e_fractional_shifts=de,
                  ratio_fractional_shifts=dr, covariance_mu_fraction=v_m,
                  covariance_e_fraction=v_e, covariance_mu_e_fraction=cross,
                  covariance_ratio_linear_fraction=v_linear,
                  covariance_ratio_exact_fraction=covariance(dr))
    return values, arrays


def histogram_arrays(hist, band_name, expected_n):
    if not hist or not hasattr(hist, 'HasVertErrorBand'):
        raise RuntimeError('Missing MnvH1D histogram')
    if not hist.HasVertErrorBand(band_name):
        raise RuntimeError(f'{hist.GetName()}: missing {band_name} error band')
    band = hist.GetVertErrorBand(band_name)
    n, bins = band.GetNHists(), hist.GetNbinsX()
    if n != expected_n:
        raise RuntimeError(f'{hist.GetName()}: expected {expected_n} universes, found {n}')
    cv = np.array([hist.GetBinContent(i) for i in range(1, bins+1)])
    u = np.array([[band.GetHist(j).GetBinContent(i) for i in range(1, bins+1)] for j in range(n)])
    edges = np.array([hist.GetXaxis().GetBinLowEdge(i) for i in range(1, bins+2)])
    return cv, u, edges


def assert_match(label, a, b, rtol, atol):
    try:
        np.testing.assert_allclose(a, b, rtol=rtol, atol=atol)
    except AssertionError as exc:
        raise RuntimeError(f'{label}: closure failed; verify that files contain the same selections, POT normalization and universe ordering. {exc}') from exc
    return float(np.max(np.abs(a-b)))


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def plot(rows, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for sample in sorted(set(r['sample'] for r in rows)):
        fig, axs = plt.subplots(2, 2, figsize=(13, 8), sharex='col', constrained_layout=True)
        for col, beam in enumerate(('fhc', 'rhc')):
            rr = [r for r in rows if r['sample'] == sample and r['beam'] == beam]
            if not rr:
                axs[0,col].set_visible(False); axs[1,col].set_visible(False); continue
            x = np.arange(1, len(rr)+1)
            for key, label, color, ls in (
                ('sigma_mu_fraction','Muon flux uncertainty','#377eb8','-'),
                ('sigma_e_fraction','Electron flux uncertainty','#ff7f00','-'),
                ('sigma_ratio_independent_exact_fraction','Ratio: independent flux ensembles','#777777','--'),
                ('sigma_ratio_paired_exact_fraction','Ratio: paired flux universes','#228833','-')):
                axs[0,col].plot(x, [100*r[key] for r in rr], marker='o', ms=3, color=color, ls=ls, label=label)
            axs[0,col].set_title(f'{beam.upper()} — {sample.upper()}'); axs[0,col].set_ylabel('Fractional flux standard deviation (%)')
            axs[0,col].legend(fontsize=8); axs[0,col].grid(alpha=.2)
            axs[1,col].plot(x, [r['variance_cancelled_exact_pct'] for r in rr], 'o-', label='Exact ratio')
            axs[1,col].plot(x, [r['variance_cancelled_linear_pct'] for r in rr], 's--', label='First-order covariance formula')
            axs[1,col].axhline(0, color='k', lw=.7); axs[1,col].set_ylabel('Ratio variance cancelled (%)')
            axs[1,col].set_xlabel('Neutrino energy bin (GeV edges below)')
            axs[1,col].set_xticks(x, [f"{r['energy_low_GeV']:g}–{r['energy_high_GeV']:g}" for r in rr], rotation=55, fontsize=8)
            axs[1,col].grid(alpha=.2); axs[1,col].legend(fontsize=8)
        fig.suptitle('Cancellation of common flux uncertainty; centered 1/N ensemble moments')
        fig.savefig(out/f'flux_cancellation_{sample}.png', dpi=180)
        fig.savefig(out/f'flux_cancellation_{sample}.pdf'); plt.close(fig)


def plot_simple(rows, out):
    """Old response-study style: three uncertainty curves, separate correlation."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for sample in sorted(set(r['sample'] for r in rows)):
        for kind in ('uncertainties', 'correlation'):
            fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.5), sharey=True,
                                     constrained_layout=True)
            for ax, beam in zip(axes, ('fhc', 'rhc')):
                rr = [r for r in rows if r['sample'] == sample and r['beam'] == beam]
                if not rr:
                    ax.set_visible(False); continue
                x = np.arange(1, len(rr)+1)
                if kind == 'uncertainties':
                    for key, label, color, marker in (
                        ('sigma_mu_fraction', 'CC muon sample', '#377eb8', 'o'),
                        ('sigma_e_fraction', 'CC electron sample', '#ff7f00', 's'),
                        ('sigma_ratio_paired_exact_fraction', r'CC $\mu/e$ ratio', '#228833', '^')):
                        ax.plot(x, [100*r[key] for r in rr], marker=marker,
                                color=color, lw=2, ms=4, label=label)
                    ax.set_ylim(bottom=0)
                    ax.set_ylabel('Fractional flux uncertainty (%)')
                    ax.legend(fontsize=10)
                else:
                    ax.plot(x, [r['rho_mu_e'] for r in rr], 'o-', color='#574b90', lw=2, ms=4)
                    ax.set_ylim(0, 1.05)
                    # Keep anticorrelation visible for other future inputs.
                    if any(r['rho_mu_e'] < 0 for r in rr): ax.set_ylim(-1.05, 1.05)
                    ax.set_ylabel(r'Muon/electron flux correlation $\rho_{\mu e}$')
                ax.set_title(beam.upper() + (' (antineutrino beam)' if beam == 'rhc' else ' (neutrino beam)'))
                ax.set_xticks(x)
                ax.set_xticklabels([f"{r['energy_low_GeV']:g}–{r['energy_high_GeV']:g}" for r in rr], rotation=55, fontsize=8)
                ax.set_xlabel('Reconstructed neutrino energy bin (GeV)')
                ax.grid(alpha=.25)
            if kind == 'uncertainties':
                ymax = max(100*r[key] for r in rows if r['sample'] == sample
                           for key in ('sigma_mu_fraction', 'sigma_e_fraction', 'sigma_ratio_paired_exact_fraction'))
                axes[0].set_ylim(0, max(1., 1.1*ymax))
            title = ('Flux uncertainty before and after forming the ratio' if kind == 'uncertainties'
                     else 'Muon/electron flux correlations')
            fig.suptitle(title + f' — {sample.upper()}')
            for ext in ('png', 'pdf'):
                fig.savefig(out/f'flux_cancellation_{kind}_{sample}.{ext}', dpi=180)
            plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--from-csv', type=Path, help='Redraw an existing flux_cancellation.csv without ROOT or recomputing uncertainties')
    p.add_argument('--direct-tag', default='predNueel_noRatio_p8_v4p3p1_updatedStat_newEL')
    p.add_argument('--ratio-tag', default='predNueel_p8_v4p3p1_updatedStat_newEL')
    p.add_argument('--root-dir', type=Path, help='Default: oscillations/rootfiles next to this studies directory')
    p.add_argument('--beam', choices=('fhc','rhc','both'), default='both')
    p.add_argument('--samples', choices=('mc','data'), nargs='+', default=['mc'],
                   help='MC by default; optional data checks background-subtracted measured spectra flux bands')
    p.add_argument('--band', default='Flux')
    p.add_argument('--expected-universes', type=int, default=1000)
    p.add_argument('--rtol', type=float, default=1e-7)
    p.add_argument('--atol', type=float, default=1e-9)
    p.add_argument('--out', type=Path, default=Path('flux_cancellation_updatedStat_newEL'))
    p.add_argument('--no-plots', action='store_true')
    args = p.parse_args()
    if args.expected_universes < 2: p.error('Need at least two universes')
    if args.from_csv:
        with args.from_csv.open(newline='') as f:
            rows = [{k:(v if k in ('sample', 'beam') else float(v)) for k,v in row.items()}
                    for row in csv.DictReader(f)]
        if not rows: p.error('Saved cancellation table is empty')
        args.out.mkdir(parents=True, exist_ok=True)
        if not args.no_plots:
            plot(rows, args.out)
            plot_simple(rows, args.out)
        print(f'Redrawn saved cancellation results: {args.out.resolve()}')
        return
    # Import ROOT only after parsing: --help and numerical functions work without ROOT.
    import ROOT
    import PlotUtils  # registers MnvH1D dictionaries
    ROOT.gROOT.SetBatch(True)
    root_dir = args.root_dir or Path(__file__).resolve().parents[1]/'rootfiles'
    paths = {name: root_dir/f'NuE_stitched_hists_{tag}.root' for name,tag in [('direct',args.direct_tag),('ratio',args.ratio_tag)]}
    files = {}
    try:
        for key,path in paths.items():
            f = ROOT.TFile.Open(str(path), 'READ')
            if not f or f.IsZombie(): raise RuntimeError(f'Cannot open {path}')
            files[key] = f
        rows, archive, closure = [], {}, []
        beams = ('fhc','rhc') if args.beam == 'both' else (args.beam,)
        for sample in dict.fromkeys(args.samples):
            for beam in beams:
                def read(source, name):
                    h = files[source].Get(name)
                    if not h: raise RuntimeError(f'{paths[source]}: missing {name}')
                    return histogram_arrays(h, args.band, args.expected_universes)
                m, mu, edges = read('direct',f'{sample}_{beam}_numu_selection')
                e, eu, ee = read('direct',f'{sample}_{beam}_nue_selection')
                r, ru, re = read('ratio',f'{sample}_{beam}_ratio')
                if not (np.array_equal(edges,ee) and np.array_equal(edges,re)):
                    raise RuntimeError(f'{sample} {beam}: physical bin edges differ across samples')
                vals, arrays = cancellation(m,e,mu,eu)
                check = dict(sample=sample,beam=beam,
                    ratio_cv_max_abs_difference=assert_match('Stored ratio CV',r,m/e,args.rtol,args.atol),
                    ratio_universes_max_abs_difference=assert_match('Stored ratio universes',ru,mu/eu,args.rtol,args.atol))
                # The muon yield is reused by the ratio fit; verify it too.
                rm,rmu,rme = read('ratio',f'{sample}_{beam}_numu_selection')
                if not np.array_equal(edges,rme): raise RuntimeError('Muon bin-edge mismatch')
                check['muon_cv_max_abs_difference'] = assert_match('Muon CV across configurations',rm,m,args.rtol,args.atol)
                check['muon_universes_max_abs_difference'] = assert_match('Muon universes across configurations',rmu,mu,args.rtol,args.atol)
                closure.append(check)
                for i in range(len(m)):
                    rows.append(dict(sample=sample,beam=beam,bin=i+1,energy_low_GeV=edges[i],energy_high_GeV=edges[i+1],
                                     universes=len(mu),**{k:float(v[i]) for k,v in vals.items()}))
                for key,value in arrays.items(): archive[f'{sample}_{beam}__{key}'] = value
                archive[f'{sample}_{beam}__energy_edges_GeV'] = edges
                archive[f'{sample}_{beam}__mu_cv'] = m
                archive[f'{sample}_{beam}__e_cv'] = e
                archive[f'{sample}_{beam}__ratio_cv'] = r
        args.out.mkdir(parents=True,exist_ok=True)
        write_csv(args.out/'flux_cancellation.csv',rows)
        write_csv(args.out/'closure_checks.csv',closure)
        np.savez_compressed(args.out/'flux_cancellation_arrays.npz',**archive)
        metadata = dict(inputs={k:str(v.resolve()) for k,v in paths.items()},band=args.band,
                        input_file_stats={k:dict(size_bytes=v.stat().st_size,mtime_ns=v.stat().st_mtime_ns) for k,v in paths.items()},
                        expected_universes=args.expected_universes,ensemble_normalization='1/N, mean-centered',
                        interpretation='Positive cancellation removes variance; negative means amplification. Exact independent baseline uses all N x N numerator/denominator pairings. MC cancellation is not the complete data-minus-MC residual prior or a measurement of additional profiling information.',
                        universe_alignment='Same index across matched inputs is assumed. Ratio-universe closure is required; index semantics cannot be inferred from histograms alone.',
                        rtol=args.rtol,atol=args.atol)
        (args.out/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
        print('\nsample beam bin  mu_sigma%  e_sigma%    rho  ratio_indep%  ratio_paired%  variance_cancelled%  sigma_reduction%')
        for r in rows:
            print(f"{r['sample']:6s} {r['beam']:4s} {r['bin']:3d} {100*r['sigma_mu_fraction']:10.3f} {100*r['sigma_e_fraction']:9.3f} {r['rho_mu_e']:6.3f} {100*r['sigma_ratio_independent_exact_fraction']:13.3f} {100*r['sigma_ratio_paired_exact_fraction']:14.3f} {r['variance_cancelled_exact_pct']:20.2f} {r['sigma_reduction_exact_pct']:17.2f}")
        print(f'\nStored ratios and reused muons close to matched direct inputs. Saved: {args.out.resolve()}')
        if not args.no_plots:
            plot(rows,args.out)
            plot_simple(rows,args.out)
    finally:
        for f in files.values(): f.Close()


if __name__ == '__main__':
    main()
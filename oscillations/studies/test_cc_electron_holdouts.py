#!/usr/bin/env python3
"""Remove both CC-electron/ratio blocks; retain elastic and both CC muon blocks.

Install under oscillations/studies/. Uses OscillationFitter.DoFit unchanged.
All sterile parameters remain free. No fixed mass or extra rate parameters.
The normalized K-mode basis is built from the retained 30 bins. The omitted
electron/ratio blocks enter neither basis construction nor the pull solve
or minimized chi-square. One source configuration suffices for the core fit.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import sys

TAGS = {'direct_all': 'predNueel_noRatio_p8_v4p3p1',
        'ratio_all': 'predNueel_p8_v4p3p1'}
CORE = ('fhc_elastic', 'fhc_numu_selection', 'rhc_numu_selection')
ELECTRON = {'direct_all': ('fhc_nue_selection', 'rhc_nue_selection'),
            'ratio_all': ('fhc_ratio', 'rhc_ratio')}


class TrainingProfile:
    """Explicit unit-penalty pull solve using ONLY retained observations."""
    def __init__(self, C, modes, keep):
        from scipy.linalg import cho_factor, cho_solve, solve_triangular
        self.solve = solve_triangular
        self.cho_solve = cho_solve
        self.keep = np.asarray(keep, dtype=int)
        self.modes = modes
        self.L = np.linalg.cholesky(C[np.ix_(self.keep, self.keep)])
        self.X = self.solve(self.L, modes[self.keep], lower=True)
        self.factor = cho_factor(np.eye(modes.shape[1]) + self.X.T @ self.X)

    def evaluate(self, residual):
        y = self.solve(self.L, residual[self.keep], lower=True, check_finite=False)
        pulls = self.cho_solve(self.factor, self.X.T @ y, check_finite=False)
        after = y - self.X @ pulls
        q_resid, penalty = float(after @ after), float(pulls @ pulls)
        return q_resid + penalty, q_resid, penalty, self.modes @ pulls


def vector(hist):
    return np.array([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX()+1)])


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    global np
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--k', type=int, default=6)
    parser.add_argument('--cases', nargs='+', choices=TAGS, default=['direct_all'],
                        help='Input source; direct_all alone suffices. A second source is an optional cross-check.')
    parser.add_argument('--tag-suffix', default='_updatedStat_newEL')
    parser.add_argument('--seed', type=int, default=431)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    if args.k < 1 or args.seed < 0:
        parser.error('--k must be positive and --seed nonnegative')
    out = args.out or Path(f'cc_electron_holdouts_K{args.k}')
    out.mkdir(parents=True, exist_ok=True)
    osc = Path(__file__).resolve().parents[1]
    sys.path[:0] = [str(osc), str(osc.parent)]
    os.environ.setdefault('CCNUEROOT', str(osc.parent))
    sys.argv = [sys.argv[0]]
    # ROOT must initialize before NumPy/SciPy in the GPVM environment.
    import ROOT
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    import numpy as np
    from scipy.linalg import cho_factor, cho_solve
    from fitData_normalizedFluxPrior import NormalizedFluxPrior, RetainedModeFluxPrior
    from tools.Fitters import OscillationFitter
    from tools.StitchedHistogram import StitchedHistogram

    rows = []
    for index, case in enumerate(dict.fromkeys(args.cases)):
        np.random.seed(args.seed)
        tag = TAGS[case] + args.tag_suffix
        cfg = osc / f'HIST_CONFIG_{tag}.json'
        rootfile = osc / 'rootfiles' / f'NuE_stitched_hists_{tag}.root'
        if not cfg.is_file() or not rootfile.is_file():
            raise FileNotFoundError(f'Need {cfg} and {rootfile}')
        groups = json.loads(cfg.read_text())
        if set(groups) != set(CORE + ELECTRON[case]):
            raise ValueError(f'{case}: unexpected sample blocks {list(groups)}')
        blocks = {name: np.arange(info['start'], info['end']+1)
                  for name, info in groups.items()}
        keep = np.sort(np.concatenate([blocks[name] for name in CORE]))
        held = np.sort(np.concatenate([blocks[name] for name in ELECTRON[case]]))
        np.testing.assert_array_equal(np.sort(np.r_[keep, held]), np.arange(54))
        if len(keep) != 30 or len(held) != 24:
            raise ValueError('Expected 30 retained bins and 24 omitted bins')

        h = StitchedHistogram('sample', include_shared_ratio_stats=True)
        h.SetHistConfig(str(cfg))
        h.Load(str(rootfile))
        dh, mh = h.GetDataHistogram(), h.GetMCHistogram()
        d, m = vector(dh), vector(mh)
        C = np.asarray(h.GetCovarianceMatrix(True), dtype=float)
        A = np.asarray(h.GetAMatrix(), dtype=float)
        band = dh.GetVertErrorBand('Flux')
        if A.shape != (1000, 54) or not band or band.GetNHists() != 1000:
            raise ValueError('Expected aligned 1000-universe flux bands and 54 bins')
        A = np.array([vector(band.GetHist(u))-d for u in range(1000)]) - A
        full = NormalizedFluxPrior(C, A)
        full.assert_flux_covariance(h.GetCovarianceMatrix(False)-C)
        # Build the modes using ONLY retained bins. Identical retained inputs
        # then define the same prior regardless of their parent configuration.
        basis = RetainedModeFluxPrior(C, A, k=args.k, profile_bins=keep)
        F = basis.modes
        prior = TrainingProfile(C, F, keep)
        V = C + F @ F.T
        retained_factor = cho_factor(V[np.ix_(keep, keep)])
        cross = V[np.ix_(held, keep)]
        cond = V[np.ix_(held, held)] - cross @ cho_solve(retained_factor, cross.T)
        cond_factor = cho_factor(cond)

        def prediction(point):
            h.OscillateHistogram(*point, False, False)
            return vector(h.GetOscillatedHistogram())

        def evaluate(pred):
            r = d-pred
            q, residual, penalty, shift = prior.evaluate(r)
            alpha = cho_solve(retained_factor, r[keep])
            np.testing.assert_allclose(q, r[keep] @ alpha, rtol=1e-8, atol=1e-7)
            cond_shift = cross @ alpha
            leftover = r[held]-cond_shift
            conditional = float(leftover @ cho_solve(cond_factor, leftover))
            return q, residual, penalty, shift, conditional, cond_shift

        np.testing.assert_allclose(prediction([0., 0., 0., 0.]), m, rtol=1e-7, atol=1e-8)
        null = evaluate(m)

        class HeldOutFitter(OscillationFitter):
            def CalChi2(self, x):
                pred = prediction([100.*x[0], *x[1:]])
                return prior.evaluate(d-pred)[0]

        print(f'INPUT {case}: {rootfile}', flush=True)
        print(f'K={args.k}; retained-bin basis; fit/score ONLY 30 bins: {CORE}; '
              f'remove 24 bins: {ELECTRON[case]}', flush=True)
        fitter = HeldOutFitter(h, exclude='', marginalize_flux=False, hist_config=str(cfg))
        score, params = fitter.DoFit()
        point = [float(params[name]) for name in ('m', 'ue4', 'umu4', 'utau4')]
        best_pred = prediction(point)
        best = evaluate(best_pred)
        np.testing.assert_allclose(score, best[0], rtol=1e-8, atol=1e-7)
        # The null is always an allowed candidate, including a flat likelihood.
        if null[0] < best[0]:
            point, best_pred, best = [0., 0., 0., 0.], m.copy(), null
        result = dict(case=case, k=args.k, hist_config_tag=tag,
                      flux_source='mean_centered_residual', basis_scope='retained_30_bins',
                      retained_samples='+'.join(CORE), removed_samples='+'.join(ELECTRON[case]),
                      n_profile_bins=30, n_score_bins=30, n_held=24,
                      chi2_null=null[0], chi2_best=best[0], delta_chi2=null[0]-best[0],
                      residual_null=null[1], penalty_null=null[2],
                      residual_best=best[1], penalty_best=best[2],
                      dm2=point[0], ue4_sq=point[1], umu4_sq=point[2], utau4_sq=point[3],
                      held_conditional_null=null[4], held_conditional_best=best[4],
                      seed=args.seed)
        rows.append(result)
        write_csv(out/'holdout_summary.csv', rows)
        (out/f'{case}_K{args.k}_fit.json').write_text(json.dumps(result, indent=2)+'\n')
        spectra = []
        for i in range(len(d)):
            j = int(np.flatnonzero(held == i)[0]) if i in held else None
            spectra.append(dict(bin=i+1, used_in_fit=bool(i in keep), data=float(d[i]),
                                null_cv=float(m[i]), best_cv=float(best_pred[i]),
                                null_flux_shift=float(null[3][i]), best_flux_shift=float(best[3][i]),
                                null_conditional_prediction=float(m[i]+null[5][j]) if j is not None else None,
                                best_conditional_prediction=float(best_pred[i]+best[5][j]) if j is not None else None))
        write_csv(out/f'{case}_K{args.k}_spectra.csv', spectra)
        print(f'{case}: chi2_null={null[0]:.6f}; chi2_best={best[0]:.6f}; '
              f'delta_chi2={null[0]-best[0]:.6f}; dm2={point[0]:.6f}; '
              f'Ue4^2={point[1]:.6g}; Umu4^2={point[2]:.6g}; Utau4^2={point[3]:.6g}', flush=True)
    print(f'Results: {out.resolve()}', flush=True)


if __name__ == '__main__':
    main()

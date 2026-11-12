#!/usr/bin/env python3
"""Check the uncertainty implied by raw-universe quadratic profiling.

Synthetic mode needs NumPy and SciPy. With --current-root, install this file in
oscillations/studies/ and use the usual GPVM ROOT/PlotUtils environment.
No sterile fit runs, no input ROOT file is changed, and no K truncation is used.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys



def runtime_state():
    """No scientific imports: useful even if ROOT aborts at initialization."""
    libs=[]
    path=Path('/proc/self/maps')
    if path.exists():
        libs=sorted({line.split()[-1] for line in path.read_text().splitlines()
                     if 'libstdc++' in line or 'libgcc_s' in line})
    return dict(python_executable=sys.executable,python_version=sys.version,
                rootsys=os.environ.get('ROOTSYS'),loaded_cxx_libraries=libs)


def initialize_runtime(args):
    global np, cho_factor, cho_solve, solve_triangular
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/'runtime_before_imports.json').write_text(json.dumps(runtime_state(),indent=2)+'\n')
    version=None
    if args.current_root:
        print('Loading ROOT before NumPy/SciPy...',flush=True)
        sys.argv=[sys.argv[0]]
        import ROOT
        ROOT.PyConfig.IgnoreCommandLineOptions=True
        ROOT.gROOT.SetBatch(True)
        version=ROOT.gROOT.GetVersion()
        print(f'ROOT {version} initialized.',flush=True)
    import numpy as np
    import scipy
    from scipy.linalg import cho_factor, cho_solve, solve_triangular
    info=runtime_state()
    info.update(root_version=version,numpy_version=np.__version__,scipy_version=scipy.__version__)
    (args.out/'runtime_after_imports.json').write_text(json.dumps(info,indent=2)+'\n')
    return info


def latex_escape(s):
    # Case labels and test names only; numeric cells are separately formatted.
    return str(s).replace('\\',r'\textbackslash{}').replace('_',r'\_').replace('%',r'\%').replace('&',r'\&').replace('#',r'\#')


def write_tex_report(out, scores, checks, scales, metadata):
    actual=metadata['mode']=='production ROOT null inputs'
    label='Actual stitched-input tests' if actual else 'Reproducible synthetic tests'
    lines=[r'\section{'+label+r': generated numerical report}',
           'This section is generated directly from the accompanying test output. '
           'It concerns a fixed linear response, unrestricted coefficients, and identical pull-solve and scoring bins.',
           r'\par\medskip\noindent\textbf{Run status:} '+latex_escape(metadata['run_status'])+'.',
           r'\par\noindent Assertion tolerance: relative '+f"{metadata['rtol']:g}"+
           '; absolute tolerance for general algebra checks is the relative tolerance times 0.01. '
           'The legacy-shift comparison uses an absolute tolerance of $10^{-6}$.']
    passed=[r for r in checks if r['status']=='PASS']
    lines += [r'\par\noindent Passed numerical assertions: '+str(len(passed))+'.']
    if passed:
        mx=max(float(r['relative_difference']) for r in passed)
        lines += [r'Largest reported relative difference among passed assertions: '+f'{mx:.4g}'+'.']
    if metadata.get('error'):
        lines += [r'\par\noindent\textbf{Failure:} '+latex_escape(metadata['error'])+'. No successful overall result is claimed.']
    methods=['raw_lambda1','normalized_lambda1','duplicated_raw_lambda1','duplicated_normalized_lambda1']
    if scores:
        lines += [r'\subsection{Profile scores before and after duplicating every universe}',
                  r'\begin{center}\small',r'\begin{tabular}{lrrrr}',r'\toprule',
                  r'Case / response / centering & Raw & Normalized & Raw dup. & Norm. dup.\\',r'\midrule']
        for case in sorted(set(r['case'] for r in scores)):
            vals={r['method']:r['score'] for r in scores if r['case']==case}
            if not all(k in vals for k in methods):continue
            name=case.replace('direct_all','Direct').replace('ratio_all','Ratios').replace('synthetic_','Toy / ').replace('rank_deficient','rank deficient').replace('_',' / ')
            lines.append(latex_escape(name)+' & '+' & '.join(f'{vals[k]:.6f}' for k in methods)+r'\\')
        lines += [r'\bottomrule\end{tabular}\end{center}',
                  'Duplication leaves the supplied distribution of universe shifts unchanged. '
                  'Agreement of the two normalized scores tests representation invariance. '
                  r'A change in the raw score at fixed unit penalty is compared separately against the predicted $\lambda/2$ identity.',
                  r'\subsection{Subset representativeness is measured separately}',
                  r'\begin{center}\small\begin{tabular}{lr}\toprule',
                  r'Case & $\|V_{\rm subset}-V_{\rm ens}\|_F/\|V_{\rm ens}\|_F$\\\midrule']
        for case in sorted(set(r['case'] for r in scores)):
            row=next(r for r in scores if r['case']==case)
            lines.append(latex_escape(case)+f" & {row['subset_covariance_relative_difference']:.6g}"+r'\\')
        lines += [r'\bottomrule\end{tabular}\end{center}',
                  'This is a descriptive Frobenius-norm difference in the original observable coordinates, not a coverage criterion or proof that the subset is adequate.']
        lines += [r'\subsection{Duplication test on the selected raw subset}',
                  r'This repeats the representation test on the selected $N_{\rm prof}$ rows themselves. Both copies are allowed in the duplicated solve. The normalized columns use the subset\textquotesingle s own second moment ($1/N_{\rm prof}$), solely for this invariance test; they do not replace the full-ensemble reference covariance.',
                  r'\begin{center}\small\begin{tabular}{lrrrr}\toprule',
                  r'Case & Raw & Raw dup. & Subset norm. & Norm. dup.\\\midrule']
        for case in sorted(set(r['case'] for r in scores)):
            vals={r['method']:r['score'] for r in scores if r['case']==case}
            keys=['subset_raw_lambda1','subset_duplicated_raw_lambda1','subset_normalized_lambda1','subset_duplicated_normalized_lambda1']
            if all(k in vals for k in keys):
                name=case.replace('direct_all','Direct').replace('ratio_all','Ratios').replace('synthetic_','Toy / ').replace('rank_deficient','rank deficient').replace('_',' / ')
                lines.append(latex_escape(name)+' & '+' & '.join(f'{vals[k]:.6f}' for k in keys)+r'\\')
        lines += [r'\bottomrule\end{tabular}\end{center}']
    if scales:
        lines += [r'\subsection{Direct measurement of the implied uncertainty scale}',
                  r'The table reports the minimum and maximum over nonzero-variance bins of $\sqrt{(V_{\rm implied})_{ii}/(V_{\rm ens})_{ii}}$ at raw $\lambda=1$. The full-ensemble column must equal $\sqrt{N}$. The subset column is measured, not assumed to equal $\sqrt{N_{\rm prof}}$.',
                  r'\begin{center}\small\begin{tabular}{lrr}\toprule',
                  r'Case & Full raw: range & Selected raw subset: range\\\midrule']
        for case in sorted(set(row['case'] for row in scales)):
            cells=[]
            for basis in ['full','subset']:
                vals=[row['sigma_ratio_to_full_ensemble'] for row in scales
                      if row['case']==case and row['basis']==basis and row['lambda_value']==1.
                      and np.isfinite(row['sigma_ratio_to_full_ensemble'])]
                cells.append(f'{min(vals):.4f}--{max(vals):.4f}' if vals else 'undefined')
            lines.append(latex_escape(case)+' & '+' & '.join(cells)+r'\\')
        lines += [r'\bottomrule\end{tabular}\end{center}']
    if passed:
        lines += [r'\subsection{Numerical assertion summary}',
                  r'\begin{center}\small\begin{tabular}{lrr}\toprule',
                  r'Case & Passed checks & Maximum relative difference\\\midrule']
        for case in sorted(set(row['case'] for row in passed)):
            rows=[row for row in passed if row['case']==case]
            lines.append(latex_escape(case)+f" & {len(rows)} & {max(float(row['relative_difference']) for row in rows):.4g}"+r'\\')
        lines += [r'\bottomrule\end{tabular}\end{center}']
    diagnostics=[r for r in checks if r['status']=='DIAGNOSTIC_ONLY']
    if diagnostics:
        lines += [r'\subsection{ROOT covariance convention diagnostic}']
        for row in diagnostics:
            lines += [latex_escape(row['case'])+': relative difference in non-flux-whitened Frobenius norm = '+f"{row['relative_difference']:.6g}"+r'.\par']
        lines += ['This diagnostic has no automatic pass/fail threshold: a different centering or robust/IQR estimator must be compared on its own definition.']
    lines += [r'\subsection{Scope and reproducibility}',
              'Individual checks are recorded in '+r'\texttt{assertions.csv}'+
              '; score decompositions in '+r'\texttt{scores.csv}'+
              '; per-bin implied variance and standard-deviation ratios in '+r'\texttt{implied\_scale\_by\_bin.csv}'+
              '. Runtime versions, loaded C++ libraries, tolerances and run status accompany the report. Production runs also save input paths, SHA-256 hashes and matrix snapshots. '
              'These tests establish or falsify identities for the specified objective. They do not establish an optimal elbow, complete physical flux uncertainty, or oscillation coverage.']
    (out/'test_report.tex').write_text('\n\n'.join(lines)+'\n')


def profile(A, C, r, lam=1.):
    """SVD solve for the UNNORMALIZED A supplied here; no hidden rescaling."""
    L = np.linalg.cholesky(C)
    X = solve_triangular(L, A.T, lower=True)
    y = solve_triangular(L, r, lower=True)
    U, s, Vt = np.linalg.svd(X, full_matrices=False)
    a = Vt.T @ ((s/(s*s+lam))*(U.T@y))
    shift = A.T@a
    after = solve_triangular(L, r-shift, lower=True)
    penalty = float(lam*(a@a))
    return dict(score=float(after@after)+penalty, residual=float(after@after),
                penalty=penalty, shift=shift, a=a)


def covscore(C, V, r):
    return float(r@cho_solve(cho_factor(C+V),r))


def relative(a,b):
    return float(np.linalg.norm(np.asarray(a)-np.asarray(b))/max(np.linalg.norm(b),1e-30))


def write_csv(path, rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def audit(label, A, C, r, nprof, assertions, scores, scales, rtol):
    N,n=A.shape
    def check(name, actual, expected):
        err=relative(actual,expected)
        np.testing.assert_allclose(actual,expected,rtol=rtol,atol=rtol*1e-2,
                                   err_msg=f'{label}: {name}')
        assertions.append(dict(case=label,test=name,relative_difference=err,status='PASS'))
    V=A.T@A/N
    raw=profile(A,C,r)
    norm=profile(A/np.sqrt(N),C,r)
    check('raw profile = covariance C + A.T A',raw['score'],covscore(C,A.T@A,r))
    check('normalized profile = covariance C + A.T A/N',norm['score'],covscore(C,V,r))
    matched=profile(A,C,r,lam=N)
    check('raw lambda=N = normalized lambda=1 score',matched['score'],norm['score'])
    check('raw lambda=N = normalized lambda=1 shift',matched['shift'],norm['shift'])
    check('raw lambda=N = normalized lambda=1 penalty',matched['penalty'],norm['penalty'])
    # Duplicate every row; centered and CV-referenced arrays both retain their convention.
    duplicated=np.concatenate([A,A])
    dupraw=profile(duplicated,C,r)
    dupnorm=profile(duplicated/np.sqrt(2*N),C,r)
    check('duplicated raw covariance doubles',duplicated.T@duplicated,2*(A.T@A))
    check('duplicated normalized covariance unchanged',duplicated.T@duplicated/(2*N),V)
    check('duplicated normalized profile score unchanged',dupnorm['score'],norm['score'])
    check('duplicated normalized profile shift unchanged',dupnorm['shift'],norm['shift'])
    check('duplicated raw score = original raw lambda=0.5',dupraw['score'],profile(A,C,r,lam=.5)['score'])
    # An explicit same-shift construction, independent of what the optimizer selects.
    split=np.r_[raw['a']/2,raw['a']/2]
    check('split duplicate coefficients preserve physical shift',duplicated.T@split,raw['shift'])
    check('split duplicate coefficients halve unit penalty',split@split,(raw['a']@raw['a'])/2)
    # Independent coefficient-space normal-equation solve on a manageable subset.
    M=min(nprof,N);As=A[:M]
    inv=cho_solve(cho_factor(C),np.eye(n))
    aa=cho_solve(cho_factor(As@inv@As.T+np.eye(M)),As@inv@r)
    ss=profile(As,C,r)
    check('SVD vs coefficient normal equations: subset shift',ss['shift'],As.T@aa)
    check('SVD vs coefficient normal equations: subset coefficients',ss['a'],aa)
    check('subset raw profile = its own implied covariance',ss['score'],covscore(C,As.T@As,r))
    subset_dups=np.concatenate([As,As])
    subset_dupraw=profile(subset_dups,C,r)
    subset_norm=profile(As/np.sqrt(M),C,r)
    subset_dupnorm=profile(subset_dups/np.sqrt(2*M),C,r)
    check('selected raw subset duplication = original subset lambda=0.5',subset_dupraw['score'],profile(As,C,r,lam=.5)['score'])
    check('selected normalized subset duplication preserves score',subset_dupnorm['score'],subset_norm['score'])
    # Do NOT assume a subset is representative. Measure that separately.
    subset_relative=relative(As.T@As/M,V)
    for name,result in [('raw_lambda1',raw),('normalized_lambda1',norm),
                        ('raw_lambdaN',matched),('duplicated_raw_lambda1',dupraw),
                        ('duplicated_normalized_lambda1',dupnorm),('subset_raw_lambda1',ss),
                        ('subset_duplicated_raw_lambda1',subset_dupraw),('subset_normalized_lambda1',subset_norm),
                        ('subset_duplicated_normalized_lambda1',subset_dupnorm)]:
        scores.append(dict(case=label,method=name,original_N=N,n_bins=n,nprof=M,
                           score=result['score'],residual=result['residual'],penalty=result['penalty'],
                           subset_covariance_relative_difference=subset_relative))
    for lam in [1.,12.,float(M),float(N)]:
        for name,arr in [('full',A),('subset',As)]:
            own_N=len(arr);implied=arr.T@arr/lam
            for i in range(n):
                ref=V[i,i]
                scales.append(dict(case=label,basis=name,lambda_value=lam,n_rows=own_N,bin=i+1,
                    reference_variance=ref,implied_variance=implied[i,i],
                    variance_ratio_to_full_ensemble=implied[i,i]/ref if ref>0 else float('nan'),
                    sigma_ratio_to_full_ensemble=np.sqrt(implied[i,i]/ref) if ref>0 else float('nan')))
    print(f'{label}: PASS; raw lambda1 Q={raw["score"]:.8g}; normalized Q={norm["score"]:.8g}; '
          f'duplicated raw Q={dupraw["score"]:.8g}; duplicated normalized Q={dupnorm["score"]:.8g}',flush=True)


def root_inputs(args, checks):
    osc=Path(__file__).resolve().parents[1]
    sys.path[:0]=[str(osc),str(osc.parent)]
    os.environ.setdefault('CCNUEROOT',str(osc.parent))
    sys.argv=[sys.argv[0]]
    import ROOT
    ROOT.gROOT.SetBatch(True)
    from tools.StitchedHistogram import StitchedHistogram
    def vec(h):return np.array([h.GetBinContent(i) for i in range(1,h.GetNbinsX()+1)])
    bases={'direct_all':'predNueel_noRatio_p8_v4p3p1','ratio_all':'predNueel_p8_v4p3p1'}
    for case in args.cases:
        tag=bases[case]+args.tag_suffix
        cfg=osc/f'HIST_CONFIG_{tag}.json'
        path=osc/'rootfiles'/f'NuE_stitched_hists_{tag}.root'
        if not cfg.is_file() or not path.is_file():raise FileNotFoundError(f'Need {cfg} and {path}')
        print(f'INPUT {case}: {path}',flush=True)
        h=StitchedHistogram('sample',include_shared_ratio_stats=True)
        h.SetHistConfig(str(cfg));h.Load(str(path))
        dh,mh=h.GetDataHistogram(),h.GetMCHistogram()
        d,m=vec(dh),vec(mh)
        C=np.asarray(h.GetCovarianceMatrix(True),dtype=float)
        A=np.asarray(h.GetAMatrix(),dtype=float)
        if A.shape!=(1000,len(d)):raise ValueError(f'Expected 1000 universes, got {A.shape}')
        if args.check_legacy:
            from tools.Fitters import FluxFitter
            ff=FluxFitter(h,exclude='',lam=1.,useOsc=False,
                          profile_n_universes=min(args.nprof,len(A)),hist_config=str(cfg))
            _,actual=ff.GetProfileAMatrixAndUniverses()
            expected=A[:min(args.nprof,len(A))]
            np.testing.assert_allclose(actual,expected,rtol=args.rtol,atol=1e-9)
            q=profile(actual,C,d-m)
            np.testing.assert_allclose(actual.T@ff.GetFluxSolution(),q['shift'],rtol=args.rtol,atol=1e-6)
            checks.append(dict(case=case,test='legacy FluxFitter null raw-MC subset shift agrees',
                               relative_difference=relative(actual.T@ff.GetFluxSolution(),q['shift']),status='PASS'))
        band=dh.GetVertErrorBand('Flux')
        if not band or band.GetNHists()!=len(A):raise ValueError('Missing/alignment mismatch in data Flux band')
        residualA=np.array([vec(band.GetHist(u))-d for u in range(len(A))])-A
        target=np.asarray(h.GetCovarianceMatrix(False),dtype=float)-C
        np.savez_compressed(args.out/f'{case}_input_arrays.npz',A_mc_cv=A,A_residual_cv=residualA,C_nonflux=C,residual=d-m,root_flux_covariance=target)
        def sha256(p):
            digest=hashlib.sha256()
            with p.open('rb') as stream:
                for block in iter(lambda:stream.read(1024*1024),b''):
                    digest.update(block)
            return digest.hexdigest()
        info=dict(root_path=str(path.resolve()),config_path=str(cfg.resolve()),root_size_bytes=path.stat().st_size,root_mtime_ns=path.stat().st_mtime_ns,
                  root_sha256=sha256(path),config_sha256=sha256(cfg),script_sha256=sha256(Path(__file__)))
        (args.out/f'{case}_input_provenance.json').write_text(json.dumps(info,indent=2)+'\n')
        centered=residualA-residualA.mean(axis=0)
        L=np.linalg.cholesky(C)
        def white(V):return solve_triangular(L,solve_triangular(L,V,lower=True).T,lower=True).T
        delta=white(centered.T@centered/len(A)-target)
        # Diagnostic only: ROOT may use a different covariance estimator/convention.
        checks.append(dict(case=case,test='ROOT Flux covariance vs centered residual moment (diagnostic)',
            relative_difference=float(np.linalg.norm(delta)/max(np.linalg.norm(white(target)),1e-30)),
            status='DIAGNOSTIC_ONLY'))
        for source,aa in [('mc',A),('residual',residualA)]:
            for centering in ['cv','mean']:
                used=aa if centering=='cv' else aa-aa.mean(axis=0)
                yield f'{case}_{source}_{centering}',used,C,d-m


def plots(out,scores):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for case in sorted(set(r['case'] for r in scores)):
        rr=[r for r in scores if r['case']==case and r['method'] in
            ('raw_lambda1','normalized_lambda1','duplicated_raw_lambda1','duplicated_normalized_lambda1')]
        fig,ax=plt.subplots(figsize=(9,4.8),constrained_layout=True)
        ax.bar(np.arange(len(rr)),[r['score'] for r in rr],color=['#b36521','#087f8c','#d6a775','#6db1b9'])
        ax.set_xticks(np.arange(len(rr)));ax.set_xticklabels(['Raw','Normalized','Raw, duplicated','Normalized, duplicated'],fontsize=9)
        ax.set_ylabel('Minimized penalized score');ax.set_title(case.replace('_',' '))
        for i,r in enumerate(rr):ax.annotate(f'{r["score"]:.4g}',(i,r['score']),xytext=(0,4),textcoords='offset points',ha='center')
        ax.set_ylim(0,max(r['score'] for r in rr)*1.15);ax.grid(axis='y',alpha=.2)
        fig.savefig(out/f'{case}_duplication.png',dpi=160);plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--current-root',action='store_true',help='Read the actual stitched ROOT files; otherwise run reproducible synthetic tests')
    p.add_argument('--cases',nargs='+',choices=['direct_all','ratio_all'],default=['direct_all','ratio_all'])
    p.add_argument('--tag-suffix',default='_updatedStat_newEL')
    p.add_argument('--nprof',type=int,default=100,help='Subset size for raw profiling and independent normal-equation checks')
    p.add_argument('--check-legacy',action='store_true',help='Also check the actual legacy FluxFitter at the null (ROOT mode only)')
    p.add_argument('--out',type=Path,default=Path('flux_prior_scale_tests'))
    p.add_argument('--rtol',type=float,default=1e-6)
    p.add_argument('--no-plots',action='store_true')
    p.add_argument('--runtime-check-only',action='store_true',help='Initialize ROOT first (with --current-root), then NumPy/SciPy; write runtime diagnostics and exit')
    args=p.parse_args()
    if args.nprof<1:p.error('--nprof must be positive')
    if args.check_legacy and not args.current_root:p.error('--check-legacy requires --current-root')
    runtime=initialize_runtime(args)
    if args.runtime_check_only:
        print(f'Runtime initialized; diagnostics saved in {args.out.resolve()}')
        return
    assertions,scores,scales=[],[],[]
    if args.current_root:
        inputs=root_inputs(args,assertions)
    else:
        rng=np.random.default_rng(20261111);n,N=8,120
        H=rng.normal(size=(n,n));C=H@H.T+np.eye(n)
        A=rng.normal(size=(N,n))@np.diag(np.linspace(.3,2,n))
        r=rng.normal(size=n)*4
        Arank=rng.normal(size=(N,3))@rng.normal(size=(3,n))
        inputs=[('synthetic_cv',A,C,r),('synthetic_mean',A-A.mean(axis=0),C,r),
                ('synthetic_rank_deficient',Arank,C,r)]
    failure=None
    label='input loading'
    try:
        for label,A,C,r in inputs:
            audit(label,A,C,r,args.nprof,assertions,scores,scales,args.rtol)
    except Exception as exc:
        failure=exc
        assertions.append(dict(case=label,test=str(exc),relative_difference=float('nan'),status='FAIL'))
    args.out.mkdir(parents=True,exist_ok=True)
    write_csv(args.out/'assertions.csv',assertions)
    if scores:write_csv(args.out/'scores.csv',scores)
    if scales:write_csv(args.out/'implied_scale_by_bin.csv',scales)
    metadata=dict(mode='production ROOT null inputs' if args.current_root else 'synthetic only',
        full_basis=True,regularization_parameter_scan=False,sterile_fit=False,
        run_status='FAIL' if failure else 'PASS',error=str(failure) if failure else None,
        rtol=args.rtol,runtime=runtime,
        subset_size=args.nprof,cases=args.cases,tag_suffix=args.tag_suffix,
        assumptions='Fixed response, positive definite C, unrestricted coefficients, identical solve/score bins. No claim of coverage or preferred K.',
        note='ROOT covariance comparison is diagnostic: no assumption that a robust/IQR estimator equals a centered second moment.')
    (args.out/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    write_tex_report(args.out,scores,assertions,scales,metadata)
    if scores and not args.no_plots:plots(args.out,scores)
    print(f'Saved {len(assertions)} checks and diagnostics and test_report.tex to {args.out.resolve()}')
    if failure:raise failure


if __name__=='__main__':main()

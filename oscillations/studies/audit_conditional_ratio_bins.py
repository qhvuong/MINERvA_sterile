#!/usr/bin/env python3
"""Locate the v4.3.1 null ratio tension after conditioning on core data."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("arrays", type=Path)
    p.add_argument("--config", type=Path, required=True, help="Ratio-file HIST_CONFIG JSON")
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()

    cfg = json.loads(args.config.read_text())
    names = ("fhc_ratio", "rhc_ratio")
    blocks = [np.arange(cfg[n]["start"], cfg[n]["end"]+1) for n in names]
    ratio = np.concatenate(blocks)
    prefix = "predNueel_p8_v4p3p1__"
    with np.load(args.arrays) as a:
        A, C, data, mc = (a[prefix+k] for k in ("A","C_nonflux","data","mc"))
        core = a["ratios_omitted__solve_bins_zero_based"]
        np.testing.assert_array_equal(np.sort(ratio), a["ratios_only__solve_bins_zero_based"])
    np.testing.assert_array_equal(np.sort(np.r_[core, ratio]), np.arange(len(data)))
    if len(A) != 1000:
        p.error("expected 1000 flux universes")
    r = data-mc
    B = (A-A.mean(axis=0)).T/np.sqrt(len(A))
    cc, rr, rc = np.ix_(core, core), np.ix_(ratio, ratio), np.ix_(ratio, core)
    T = np.linalg.solve(C[cc], C[rc].T).T
    N = C[rr] - T @ C[rc].T
    # Compact factor of the same prior; avoid a 1000-dimensional solve.
    U, s, _ = np.linalg.svd(B, full_matrices=False)
    B = U*s
    H = B[ratio] - T @ B[core]
    W = np.linalg.solve(np.linalg.cholesky(C[cc]), B[core])
    L = np.linalg.cholesky(np.eye(B.shape[1])+W.T@W)
    F = np.linalg.solve(L, H.T).T
    S = C+B@B.T
    pred = mc[ratio] + S[rc] @ np.linalg.solve(S[cc], r[core])
    V = N+F@F.T
    np.testing.assert_allclose(V, S[rr]-S[rc]@np.linalg.solve(S[cc],S[rc].T), atol=1e-9)
    e = data[ratio]-pred
    P = np.linalg.solve(V, np.eye(len(ratio)))
    q = e@P@e
    sig = np.sqrt(np.diag(V))
    loo = (P@e)/np.sqrt(np.diag(P))
    # Diagonalize the remaining flux covariance relative to the conditional non-flux noise.
    LN = np.linalg.cholesky(N)
    U, s, _ = np.linalg.svd(np.linalg.solve(LN,F), full_matrices=False)
    z = U.T@np.linalg.solve(LN,e)
    mode_q = z*z/(1+s*s)
    np.testing.assert_allclose(q, mode_q.sum())
    mode = (LN@U[:,0])*s[0]
    if mode[np.argmax(abs(mode))] < 0:
        mode = -mode
    amplitude = mode@np.linalg.solve(N,e)/(1+s[0]**2)
    print(f"Conditional ratio chi2 = {q:.4f} ({len(ratio)} bins)")
    print(f"Leading mode chi2 = {mode_q[0]:.4f}/{q:.4f}; fitted pull = {amplitude:+.3f}")
    for k in range(min(5,len(s))):
        print(f"  mode {k+1}: eigenvalue={s[k]**2:.4f}, chi2={mode_q[k]:.4f}, "
              f"chi2 reduction={mode_q[k]*s[k]**2:.4f}")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2,2,figsize=(11,6.5),sharex="col",layout="constrained",
                            gridspec_kw={"height_ratios":[2.2,1]})
    rows = []
    start = 0
    for col,(name,block) in enumerate(zip(names, blocks)):
        j = np.arange(start,start+len(block))
        start += len(block)
        x = np.arange(1,len(j)+1)
        marginal = e[j]@np.linalg.solve(V[np.ix_(j,j)],e[j])
        other = np.setdiff1d(np.arange(len(e)), j)
        other_score = e[other] @ np.linalg.solve(V[np.ix_(other,other)], e[other])
        print(f"{name}: marginal chi2 {marginal:.4f}; extra given other beam {q-other_score:.4f}")
        top, bot = axes[:,col]
        top.fill_between(x,pred[j]-sig[j],pred[j]+sig[j],color="C0",alpha=.18,label="Conditional 1σ (per bin)")
        top.plot(x,pred[j],color="C0",label="Core-conditioned prediction")
        top.plot(x,pred[j]+amplitude*mode[j],ls="--",color="C1",label="With fitted leading mode only")
        top.plot(x,data[ratio[j]],"ko",ms=4,label="Ratio data")
        top.set_title(f"{name[:3].upper()} ratios • marginal χ² = {marginal:.2f}")
        bot.axhline(0,color=".5",lw=1)
        bot.axhspan(-2,2,color=".92")
        bot.plot(x,loo[j],"o",color="C3")
        bot.set(xlabel="Local ratio bin",xticks=x)
        for k,i in enumerate(j,1):
            rows.append(dict(sample=name,sample_bin=k,global_bin=int(ratio[i])+1,
                             data=data[ratio[i]],cv=mc[ratio[i]],conditional_prediction=pred[i],
                             conditional_sigma=sig[i],loo_pull=loo[i],leading_flux_1sigma=mode[i],
                             fitted_leading_shift=amplitude*mode[i]))
    axes[0,0].set_ylabel("CC ratio")
    axes[1,0].set_ylabel("Leave-one-bin-out pull")
    axes[0,0].legend(fontsize=8)
    fig.suptitle("MINERvA LE v4.3.1 • ratios predicted from the core at the null")
    args.out.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.out/"v4p3p1_conditional_ratio_bins.png",dpi=180)
    with (args.out/"conditional_ratio_bins.csv").open("w",newline="") as f:
        writer = csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for row in sorted(rows,key=lambda row:abs(row["loo_pull"]),reverse=True)[:5]:
        print(f"  {row['sample']} bin {row['sample_bin']}: leave-one-out pull {row['loo_pull']:+.3f}")


if __name__ == "__main__":
    main()

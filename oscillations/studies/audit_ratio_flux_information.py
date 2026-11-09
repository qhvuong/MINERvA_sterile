#!/usr/bin/env python3
"""Null-only v4.3.1 flux information and normalized L-curves from saved arrays."""
import argparse
from pathlib import Path

import numpy as np

from compare_normalized_flux_profile import TAGS, mean_center, write_csv


def whiten(cov, response, residual):
    L = np.linalg.cholesky(cov)
    return np.linalg.solve(L, response), np.linalg.solve(L, residual)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("arrays", type=Path)
    p.add_argument("--out", type=Path, help="Optional directory for curves and information CSVs")
    args = p.parse_args()
    rows, curves = [], []

    def scan(label, W, y):
        U, s, _ = np.linalg.svd(W, full_matrices=False)
        s2 = s*s
        z = U.T @ y
        total_at_one = None
        for lam in np.logspace(-4, 4, 161):
            pull2 = np.sum((s*z/(s2+lam))**2)
            left = y - U @ ((s2/(s2+lam))*z)
            residual = left @ left
            total = residual + lam*pull2
            curves.append(dict(case=label, lam=lam, pull2=pull2,
                               residual_chi2=residual, total=total))
            if lam == 1:
                total_at_one = total
        gain = s2/(1+s2)
        rows.append(dict(case=label, bins=len(y), leading_eigenvalue=s2[0],
                         strongest_variance_reduction=gain[0],
                         modes_halved=int(np.sum(s2 >= 1)), effective_flux_directions=gain.sum(),
                         leading_eigenvalues=" ".join(f"{v:.6g}" for v in s2[:6])))
        return total_at_one

    with np.load(args.arrays) as archive:
        def get(case):
            A, C, d, m = [archive[TAGS[case]+"__"+key]
                          for key in ("A", "C_nonflux", "data", "mc")]
            if len(A) != 1000:
                p.error("expected 1000 P8 universes")
            _, s, Vt = np.linalg.svd(mean_center(A)/np.sqrt(len(A)), full_matrices=False)
            return C, (s[:, None]*Vt).T, d-m

        C, B, r = get("direct_all")
        scan("direct_all", *whiten(C, B, r))
        C, B, r = get("ratio_all")
        core = archive["ratios_omitted__solve_bins_zero_based"]
        ratio = archive["ratios_only__solve_bins_zero_based"]
    np.testing.assert_array_equal(np.sort(np.r_[core, ratio]), np.arange(len(r)))
    cc, rr = np.ix_(core, core), np.ix_(ratio, ratio)
    rc = np.ix_(ratio, core)
    W, y = whiten(C[cc], B[core], r[core])
    qcore = scan("core_only", W, y)
    scan("ratios_only", *whiten(C[rr], B[ratio], r[ratio]))
    qall = scan("ratio_all", *whiten(C, B, r))

    # Condition the non-flux noise on the core, then propagate the core flux
    # posterior to that conditional ratio measurement.
    T = np.linalg.solve(C[cc], C[rc].T).T
    Ccond = C[rr] - T @ C[rc].T
    Bcond = B[ratio] - T @ B[core]
    rcond = r[ratio] - T @ r[core]
    L = np.linalg.cholesky(np.eye(B.shape[1]) + W.T @ W)
    post_mean = np.linalg.solve(L.T, np.linalg.solve(L, W.T @ y))
    post_response = np.linalg.solve(L, Bcond.T).T
    V = Ccond + post_response @ post_response.T
    err = rcond - Bcond @ post_mean
    qcond = err @ np.linalg.solve(V, err)
    S = C + B @ B.T
    err_ind = r[ratio] - S[rc] @ np.linalg.solve(S[cc], r[core])
    V_ind = S[rr] - S[rc] @ np.linalg.solve(S[cc], S[rc].T)
    np.testing.assert_allclose(err, err_ind, atol=1e-9)
    np.testing.assert_allclose(V, V_ind, atol=1e-9)
    np.testing.assert_allclose(qall, qcore + qcond, atol=1e-9)
    # Vary the strength of the *remaining core-posterior* covariance only.
    # This is an information diagnostic, not another scan of the original prior.
    scan("ratios_given_core", *whiten(Ccond, post_response, err))
    for row in rows:
        print(f"{row['case']:19s}  leading={row['leading_eigenvalue']:.4g}"
              f"  strongest reduction={row['strongest_variance_reduction']:.1%}"
              f"  halved={row['modes_halved']}  effective={row['effective_flux_directions']:.3f}"
              f"  top eigenvalues: {row['leading_eigenvalues']}")
    print(f"Null conditional ratio score = {qcond:.4f} ({len(ratio)} bins); "
          f"core {qcore:.4f} + conditional {qcond:.4f} = joint {qall:.4f}")
    if args.out:
        args.out.mkdir(parents=True, exist_ok=True)
        write_csv(args.out/"normalized_lcurves.csv", curves)
        write_csv(args.out/"flux_information.csv", rows)


if __name__ == "__main__":
    main()

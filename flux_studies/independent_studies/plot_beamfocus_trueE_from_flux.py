#!/usr/bin/env python3

import argparse
import math

import ROOT
import numpy as np
from matplotlib import pyplot as plt

ROOT.TH1.AddDirectory(False)


def parse_args():
    p = argparse.ArgumentParser(
        description="Study ME BeamFocus true-E uncertainty directly from produced flux ROOT files."
    )
    p.add_argument("--numu", required=True)
    p.add_argument("--nue", required=True)
    p.add_argument("--hist", default="flux_E_cvweighted")
    p.add_argument("--band", default="Flux_BeamFocus")
    p.add_argument("--emin", type=float, default=0.0)
    p.add_argument("--emax", type=float, default=20.0)
    p.add_argument("--bin-width", type=float, default=0.5)
    p.add_argument("--output", default="beamfocus_trueE_from_flux.png")
    p.add_argument("--npz", default=None)
    return p.parse_args()


def open_hist(path, hname):
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        raise RuntimeError("Could not open {}".format(path))

    h = f.Get(hname)
    if not h:
        keys = [k.GetName() for k in f.GetListOfKeys()]
        raise RuntimeError(
            "Could not find '{}' in {}\nAvailable keys include:\n{}".format(
                hname, path, "\n".join(keys[:50])
            )
        )

    h.SetDirectory(0)
    f.Close()
    return h


def get_band(hist, band_name):
    names = [str(x) for x in hist.GetVertErrorBandNames()]
    if band_name not in names:
        raise RuntimeError(
            "Histogram '{}' has no vertical error band '{}'. Available: {}".format(
                hist.GetName(), band_name, ", ".join(names)
            )
        )
    return hist.GetVertErrorBand(band_name)


def make_edges(emin, emax, width):
    n = int(round((emax - emin) / width))
    edges = emin + width * np.arange(n + 1, dtype=float)

    if not np.isclose(edges[-1], emax):
        raise ValueError("Requested range is not divisible by bin width")

    return edges


def integrate_density_hist(hist, out_edges, tol=1e-9):
    out = np.zeros(len(out_edges) - 1, dtype=float)

    for ibin in range(1, hist.GetNbinsX() + 1):
        lo = hist.GetXaxis().GetBinLowEdge(ibin)
        hi = hist.GetXaxis().GetBinUpEdge(ibin)
        value = hist.GetBinContent(ibin)

        if hi <= out_edges[0] + tol or lo >= out_edges[-1] - tol:
            continue

        for j in range(len(out)):
            overlap = min(hi, out_edges[j + 1]) - max(lo, out_edges[j])
            if overlap > tol:
                out[j] += value * overlap

    return out


def extract_rebinned(hist, band_name, out_edges):
    cv = integrate_density_hist(hist, out_edges)

    band = get_band(hist, band_name)
    nuniv = band.GetNHists()

    universes = np.zeros((len(out_edges) - 1, nuniv), dtype=float)

    for u in range(nuniv):
        universes[:, u] = integrate_density_hist(
            band.GetHist(u), out_edges
        )

    return cv, universes


def fractional_shifts(cv, universes):
    delta = np.full_like(universes, np.nan, dtype=float)
    good = np.isfinite(cv) & (cv != 0.0)
    delta[good, :] = universes[good, :] / cv[good, np.newaxis] - 1.0
    return delta


def rms_about_cv(delta):
    return np.sqrt(np.nanmean(delta * delta, axis=1))


def covariance_and_correlation(delta_a, delta_b):
    if delta_a.shape != delta_b.shape:
        raise ValueError("Universe arrays do not match")

    nbins = delta_a.shape[0]
    cov = np.full(nbins, np.nan)
    rho = np.full(nbins, np.nan)

    for i in range(nbins):
        a = delta_a[i]
        b = delta_b[i]
        mask = np.isfinite(a) & np.isfinite(b)

        if not np.any(mask):
            continue

        aa = a[mask]
        bb = b[mask]

        cov[i] = np.mean(aa * bb)

        sa = math.sqrt(np.mean(aa * aa))
        sb = math.sqrt(np.mean(bb * bb))

        if sa > 0.0 and sb > 0.0:
            rho[i] = cov[i] / (sa * sb)

    return cov, rho


def ratio_fractional_shifts(cv_numu, u_numu, cv_nue, u_nue):
    nbins, nuniv = u_numu.shape
    out = np.full((nbins, nuniv), np.nan, dtype=float)

    for i in range(nbins):
        if cv_numu[i] == 0.0 or cv_nue[i] == 0.0:
            continue

        cv_ratio = cv_numu[i] / cv_nue[i]

        for u in range(nuniv):
            if u_nue[i, u] == 0.0:
                continue

            ratio_u = u_numu[i, u] / u_nue[i, u]
            out[i, u] = ratio_u / cv_ratio - 1.0

    return out


def main():
    args = parse_args()

    edges = make_edges(args.emin, args.emax, args.bin_width)
    centers = 0.5 * (edges[:-1] + edges[1:])

    h_numu = open_hist(args.numu, args.hist)
    h_nue = open_hist(args.nue, args.hist)

    cv_mu, u_mu = extract_rebinned(h_numu, args.band, edges)
    cv_e, u_e = extract_rebinned(h_nue, args.band, edges)

    if u_mu.shape[1] != u_e.shape[1]:
        raise RuntimeError(
            "numu and nue have different numbers of universes: {} vs {}".format(
                u_mu.shape[1], u_e.shape[1]
            )
        )

    d_mu = fractional_shifts(cv_mu, u_mu)
    d_e = fractional_shifts(cv_e, u_e)

    sigma_mu = rms_about_cv(d_mu)
    sigma_e = rms_about_cv(d_e)

    cov_mu_e, rho_mu_e = covariance_and_correlation(d_mu, d_e)

    d_ratio = ratio_fractional_shifts(cv_mu, u_mu, cv_e, u_e)
    sigma_ratio = rms_about_cv(d_ratio)

    print("numu file :", args.numu)
    print("nue file  :", args.nue)
    print("hist      :", args.hist)
    print("band      :", args.band)
    print("universes :", u_mu.shape[1])
    print()
    print("{:>12s} {:>10s} {:>10s} {:>10s} {:>12s}".format(
        "Etrue [GeV]", "sig_mu", "sig_e", "rho", "sig_mu/e"
    ))

    for i in range(len(centers)):
        print(
            "{:5.1f}-{:<5.1f} {:9.3f}% {:9.3f}% {:10.4f} {:11.3f}%".format(
                edges[i],
                edges[i + 1],
                100.0 * sigma_mu[i],
                100.0 * sigma_e[i],
                rho_mu_e[i],
                100.0 * sigma_ratio[i],
            )
        )

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(centers, 100.0 * sigma_mu, marker="o", label=r"$\nu_\mu$")
    ax.plot(centers, 100.0 * sigma_e, marker="o", label=r"$\nu_e$")
    ax.plot(centers, 100.0 * sigma_ratio, marker="o", label=r"$\nu_\mu/\nu_e$")

    ax.set_xlabel(r"True $E_\nu$ [GeV]")
    ax.set_ylabel("BeamFocus fractional uncertainty [%]")
    ax.set_xlim(args.emin, args.emax)
    ax.grid(True, alpha=0.25)
    ax.legend()

    fig.tight_layout()
    fig.savefig(args.output, dpi=160)
    plt.close(fig)

    print()
    print("Saved plot:", args.output)

    if args.npz:
        np.savez(
            args.npz,
            edges=edges,
            centers=centers,
            cv_mu=cv_mu,
            cv_e=cv_e,
            universes_mu=u_mu,
            universes_e=u_e,
            delta_mu=d_mu,
            delta_e=d_e,
            sigma_mu=sigma_mu,
            sigma_e=sigma_e,
            covariance_mu_e=cov_mu_e,
            rho_mu_e=rho_mu_e,
            delta_ratio=d_ratio,
            sigma_ratio=sigma_ratio,
        )
        print("Saved NPZ :", args.npz)


if __name__ == "__main__":
    main()
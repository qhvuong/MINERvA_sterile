#!/usr/bin/env python3

import argparse
import math

import ROOT
import numpy as np
from matplotlib import pyplot as plt

ROOT.TH1.AddDirectory(False)


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Study ME BeamFocus true-E uncertainty directly from produced flux ROOT files "
            "using the native histogram binning."
        )
    )
    p.add_argument("--numu", required=True, help="Produced ME numu flux ROOT file")
    p.add_argument("--nue", required=True, help="Produced ME nue flux ROOT file")
    p.add_argument("--hist", default="flux_E_cvweighted")
    p.add_argument("--band", default="Flux_BeamFocus")
    p.add_argument("--output", default="beamfocus_trueE_from_flux_native.png")
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


def native_edges(hist):
    ax = hist.GetXaxis()
    return np.array(
        [ax.GetBinLowEdge(i) for i in range(1, hist.GetNbinsX() + 1)]
        + [ax.GetBinUpEdge(hist.GetNbinsX())],
        dtype=float,
    )


def assert_same_binning(h1, h2, tol=1e-12):
    e1 = native_edges(h1)
    e2 = native_edges(h2)

    if len(e1) != len(e2) or not np.allclose(e1, e2, rtol=0.0, atol=tol):
        raise RuntimeError(
            "numu and nue histograms do not have identical native binning.\n"
            "numu edges: {}\n"
            "nue edges : {}".format(e1, e2)
        )

    return e1


def extract_native(hist, band_name):
    nbins = hist.GetNbinsX()
    cv = np.array(
        [hist.GetBinContent(i) for i in range(1, nbins + 1)],
        dtype=float,
    )

    band = get_band(hist, band_name)
    nuniv = band.GetNHists()

    universes = np.zeros((nbins, nuniv), dtype=float)

    for u in range(nuniv):
        hu = band.GetHist(u)
        universes[:, u] = [
            hu.GetBinContent(i) for i in range(1, nbins + 1)
        ]

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
        raise ValueError(
            "Universe arrays do not match: {} vs {}".format(
                delta_a.shape, delta_b.shape
            )
        )

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
    if u_numu.shape != u_nue.shape:
        raise ValueError("numu and nue universe arrays have different shapes")

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

    h_numu = open_hist(args.numu, args.hist)
    h_nue = open_hist(args.nue, args.hist)

    edges = assert_same_binning(h_numu, h_nue)
    centers = 0.5 * (edges[:-1] + edges[1:])

    cv_mu, u_mu = extract_native(h_numu, args.band)
    cv_e, u_e = extract_native(h_nue, args.band)

    if u_mu.shape[1] != u_e.shape[1]:
        raise RuntimeError(
            "numu and nue have different numbers of {} universes: {} vs {}".format(
                args.band, u_mu.shape[1], u_e.shape[1]
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
    print("nbins     :", len(centers))
    print("universes :", u_mu.shape[1])
    print()
    print("{:>18s} {:>10s} {:>10s} {:>10s} {:>12s}".format(
        "Etrue bin [GeV]", "sig_mu", "sig_e", "rho", "sig_mu/e"
    ))

    for i in range(len(centers)):
        print(
            "{:7.3f}-{:7.3f} {:9.3f}% {:9.3f}% {:10.4f} {:11.3f}%".format(
                edges[i],
                edges[i + 1],
                100.0 * sigma_mu[i],
                100.0 * sigma_e[i],
                rho_mu_e[i],
                100.0 * sigma_ratio[i],
            )
        )

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.plot(
        centers,
        100.0 * sigma_mu,
        marker="o",
        label=r"$\nu_\mu$",
    )
    ax.plot(
        centers,
        100.0 * sigma_e,
        marker="o",
        label=r"$\nu_e$",
    )
    ax.plot(
        centers,
        100.0 * sigma_ratio,
        marker="o",
        label=r"$\nu_\mu/\nu_e$",
    )

    ax.set_xlabel(r"True $E_\nu$ [GeV]")
    ax.set_ylabel("BeamFocus fractional uncertainty [%]")
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

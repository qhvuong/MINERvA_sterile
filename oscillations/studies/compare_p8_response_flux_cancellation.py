#!/usr/bin/env python3

import argparse
import ROOT
import numpy as np
import matplotlib.pyplot as plt

ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)


# ============================================================
# STITCHED BINNING
# ============================================================
#
# Global ROOT bins:
#
#   1-4   : nu+e
#   5-16  : FHC mu/e
#   17-28 : RHC mu/e
#   29-40 : FHC numu
#   41-52 : RHC numubar
#
# For cancellation studies:
#
#   ratio bin 5+i  <-> electron component bin 5+i
#                    muon component bin 29+i
#
#   ratio bin 17+i <-> electron component bin 17+i
#                    muon component bin 41+i
#
# For cancellation studies:
#
#   Direct configuration:
#       bins  5-16  -> FHC nue
#       bins 17-28  -> RHC nuebar
#       bins 29-40  -> FHC numu
#       bins 41-52  -> RHC numubar
#
#   Ratio configuration:
#       bins  5-16  -> FHC mu/e
#       bins 17-28  -> RHC mu/e
#
# Electron and muon shifts are therefore taken from the direct file,
# while the stored mu/e shift is taken from the ratio file.


N_CC_BINS = 12


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare flux cancellation for p8 new vs old response."
    )

    parser.add_argument(
        "--new-ratio-file",
        default="rootfiles/NuE_stitched_hists_prodNueel_p8.root",
    )

    parser.add_argument(
        "--new-direct-file",
        default="rootfiles/NuE_stitched_hists_prodNueel_noRatio_p8.root",
    )

    parser.add_argument(
        "--old-ratio-file",
        default="rootfiles/NuE_stitched_hists_prodNueel_p8_oldResponse.root",
    )

    parser.add_argument(
        "--old-direct-file",
        default="rootfiles/NuE_stitched_hists_prodNueel_noRatio_p8_oldResponse.root",
    )

    parser.add_argument(
        "--error-band",
        default="Flux",
        help="Flux error band [default: Flux]",
    )

    parser.add_argument(
        "--out-prefix",
        default="p8_response_flux_cancellation",
    )

    return parser.parse_args()


def open_root(path):
    f = ROOT.TFile.Open(path, "READ")

    if not f or f.IsZombie():
        raise RuntimeError("Could not open {}".format(path))

    return f


def get_hist(f, name):
    h = f.Get(name)

    if not h:
        raise RuntimeError(
            "Could not find histogram '{}' in {}".format(
                name,
                f.GetName(),
            )
        )

    return h


def get_cv(hist):
    return np.array(
        [
            hist.GetBinContent(i)
            for i in range(1, hist.GetNbinsX() + 1)
        ],
        dtype=float,
    )


def get_universes(hist, band_name):
    if not hist.HasVertErrorBand(band_name):
        raise RuntimeError(
            "Histogram '{}' has no vertical error band '{}'. "
            "Available bands: {}".format(
                hist.GetName(),
                band_name,
                list(hist.GetVertErrorBandNames()),
            )
        )

    band = hist.GetVertErrorBand(band_name)

    nuniv = band.GetNHists()
    nbins = hist.GetNbinsX()

    values = np.zeros((nuniv, nbins), dtype=float)

    for u in range(nuniv):
        hu = band.GetHist(u)

        for i in range(1, nbins + 1):
            values[u, i - 1] = hu.GetBinContent(i)

    return values


def fractional_shifts(hist, band_name):
    cv = get_cv(hist)
    universes = get_universes(hist, band_name)

    delta = np.full_like(universes, np.nan, dtype=float)

    good = cv != 0.0

    delta[:, good] = (
        universes[:, good] - cv[good]
    ) / cv[good]

    return cv, delta


def std(x):
    return np.std(x, ddof=1)


def correlation(x, y):
    if std(x) == 0.0 or std(y) == 0.0:
        return 0.0

    return np.corrcoef(x, y)[0, 1]


def analyze_response_case(
    ratio_file,
    direct_file,
    band_name,
    beam,
):
    h_ratio = get_hist(ratio_file, "sample_mc")
    h_direct = get_hist(direct_file, "sample_mc")

    cv_ratio, d_ratio = fractional_shifts(
        h_ratio,
        band_name,
    )

    cv_direct, d_direct = fractional_shifts(
        h_direct,
        band_name,
    )

    if d_ratio.shape[0] != d_direct.shape[0]:
        raise RuntimeError(
            "Universe-count mismatch: ratio={}, direct={}".format(
                d_ratio.shape[0],
                d_direct.shape[0],
            )
        )

    if beam == "fhc":
        idx_e0 = 4
        idx_ratio0 = 4
        idx_mu0 = 28

    elif beam == "rhc":
        idx_e0 = 16
        idx_ratio0 = 16
        idx_mu0 = 40

    else:
        raise ValueError("beam must be 'fhc' or 'rhc'")

    sigma_mu = []
    sigma_e = []
    sigma_ratio = []
    sigma_diff = []
    rho = []

    for ibin in range(12):

        de = d_direct[:, idx_e0 + ibin]
        dm = d_direct[:, idx_mu0 + ibin]
        dr = d_ratio[:, idx_ratio0 + ibin]

        good = (
            np.isfinite(de)
            & np.isfinite(dm)
            & np.isfinite(dr)
        )

        de = de[good]
        dm = dm[good]
        dr = dr[good]

        se = np.std(de, ddof=1)
        sm = np.std(dm, ddof=1)
        sr = np.std(dr, ddof=1)

        r = np.corrcoef(dm, de)[0, 1]

        # Linearized ratio shift:
        # delta(mu/e) ~= delta_mu - delta_e
        sd = np.std(dm - de, ddof=1)

        sigma_e.append(se)
        sigma_mu.append(sm)
        sigma_ratio.append(sr)
        sigma_diff.append(sd)
        rho.append(r)

    return {
        "sigma_mu": np.asarray(sigma_mu),
        "sigma_e": np.asarray(sigma_e),
        "sigma_ratio": np.asarray(sigma_ratio),
        "sigma_diff": np.asarray(sigma_diff),
        "rho": np.asarray(rho),
    }

def print_table(label, result):
    print()
    print("=" * 100)
    print(label)
    print("=" * 100)

    print(
        "{:>4s} {:>12s} {:>12s} {:>12s} {:>14s} {:>14s}".format(
            "bin",
            "sigma_mu",
            "sigma_e",
            "rho",
            "stored mu/e",
            "mu-e approx",
        )
    )

    for i in range(len(result["sigma_ratio"])):
        print(
            "{:4d} {:11.4f}% {:11.4f}% {:12.5f} {:13.4f}% {:13.4f}%".format(
                i + 1,
                100.0 * result["sigma_mu"][i],
                100.0 * result["sigma_e"][i],
                result["rho"][i],
                100.0 * result["sigma_ratio"][i],
                100.0 * result["sigma_diff"][i],
            )
        )


def plot_uncertainty_comparison(
    old,
    new,
    beam_label,
    outname,
):
    x = np.arange(1, len(new["sigma_ratio"]) + 1)

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.plot(
        x,
        100.0 * old["sigma_mu"],
        marker="o",
        linestyle="--",
        label=r"Old response: $\sigma_\mu$",
    )

    ax.plot(
        x,
        100.0 * old["sigma_e"],
        marker="s",
        linestyle="--",
        label=r"Old response: $\sigma_e$",
    )

    ax.plot(
        x,
        100.0 * old["sigma_ratio"],
        marker="^",
        linestyle="--",
        label=r"Old response: $\sigma_{\mu/e}$",
    )

    ax.plot(
        x,
        100.0 * new["sigma_mu"],
        marker="o",
        label=r"New response: $\sigma_\mu$",
    )

    ax.plot(
        x,
        100.0 * new["sigma_e"],
        marker="s",
        label=r"New response: $\sigma_e$",
    )

    ax.plot(
        x,
        100.0 * new["sigma_ratio"],
        marker="^",
        label=r"New response: $\sigma_{\mu/e}$",
    )

    ax.set_xlabel("CC reconstructed-energy bin")
    ax.set_ylabel("Fractional flux uncertainty [%]")
    ax.set_title("{} flux cancellation".format(beam_label))

    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=9)

    fig.tight_layout()
    fig.savefig(outname)
    plt.close(fig)


def plot_ratio_only(
    old,
    new,
    beam_label,
    outname,
):
    x = np.arange(1, len(new["sigma_ratio"]) + 1)

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(
        x,
        100.0 * old["sigma_ratio"],
        marker="o",
        label="Old response",
    )

    ax.plot(
        x,
        100.0 * new["sigma_ratio"],
        marker="o",
        label="New response",
    )

    ax.set_xlabel("CC reconstructed-energy bin")
    ax.set_ylabel(r"$\sigma_{\mu/e}$ [%]")
    ax.set_title("{} residual flux uncertainty".format(beam_label))

    ax.grid(True, alpha=0.25)
    ax.legend()

    fig.tight_layout()
    fig.savefig(outname)
    plt.close(fig)


def plot_rho(
    old,
    new,
    beam_label,
    outname,
):
    x = np.arange(1, len(new["rho"]) + 1)

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(
        x,
        old["rho"],
        marker="o",
        label="Old response",
    )

    ax.plot(
        x,
        new["rho"],
        marker="o",
        label="New response",
    )

    ax.axhline(0.0, linewidth=1)

    ax.set_xlabel("CC reconstructed-energy bin")
    ax.set_ylabel(r"$\rho_{\mu e}$")
    ax.set_title("{} flux correlation".format(beam_label))

    ax.set_ylim(-1.05, 1.05)

    ax.grid(True, alpha=0.25)
    ax.legend()

    fig.tight_layout()
    fig.savefig(outname)
    plt.close(fig)


def plot_change(
    old,
    new,
    beam_label,
    outname,
):
    x = np.arange(1, len(new["sigma_ratio"]) + 1)

    diff = 100.0 * (
        new["sigma_ratio"] - old["sigma_ratio"]
    )

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.axhline(0.0, linewidth=1)

    ax.plot(
        x,
        diff,
        marker="o",
    )

    ax.set_xlabel("CC reconstructed-energy bin")
    ax.set_ylabel(
        r"New $-$ old $\sigma_{\mu/e}$ [percentage points]"
    )

    ax.set_title(
        "{} change in residual flux uncertainty".format(
            beam_label
        )
    )

    ax.grid(True, alpha=0.25)

    fig.tight_layout()
    fig.savefig(outname)
    plt.close(fig)


def main():
    args = parse_args()

    f_new_ratio = open_root(args.new_ratio_file)
    f_new_direct = open_root(args.new_direct_file)

    f_old_ratio = open_root(args.old_ratio_file)
    f_old_direct = open_root(args.old_direct_file)

    new_fhc = analyze_response_case(
        f_new_ratio,
        f_new_direct,
        args.error_band,
        "fhc",
    )

    old_fhc = analyze_response_case(
        f_old_ratio,
        f_old_direct,
        args.error_band,
        "fhc",
    )

    new_rhc = analyze_response_case(
        f_new_ratio,
        f_new_direct,
        args.error_band,
        "rhc",
    )

    old_rhc = analyze_response_case(
        f_old_ratio,
        f_old_direct,
        args.error_band,
        "rhc",
    )

    print_table(
        "FHC — OLD RESPONSE",
        old_fhc,
    )

    print_table(
        "FHC — NEW RESPONSE",
        new_fhc,
    )

    print_table(
        "RHC — OLD RESPONSE",
        old_rhc,
    )

    print_table(
        "RHC — NEW RESPONSE",
        new_rhc,
    )

    plot_uncertainty_comparison(
        old_fhc,
        new_fhc,
        "FHC",
        args.out_prefix + "_fhc_all.png",
    )

    plot_uncertainty_comparison(
        old_rhc,
        new_rhc,
        "RHC",
        args.out_prefix + "_rhc_all.png",
    )

    plot_ratio_only(
        old_fhc,
        new_fhc,
        "FHC",
        args.out_prefix + "_fhc_ratio.png",
    )

    plot_ratio_only(
        old_rhc,
        new_rhc,
        "RHC",
        args.out_prefix + "_rhc_ratio.png",
    )

    plot_rho(
        old_fhc,
        new_fhc,
        "FHC",
        args.out_prefix + "_fhc_rho.png",
    )

    plot_rho(
        old_rhc,
        new_rhc,
        "RHC",
        args.out_prefix + "_rhc_rho.png",
    )

    plot_change(
        old_fhc,
        new_fhc,
        "FHC",
        args.out_prefix + "_fhc_change.png",
    )

    plot_change(
        old_rhc,
        new_rhc,
        "RHC",
        args.out_prefix + "_rhc_change.png",
    )

    f_new_ratio.Close()
    f_new_direct.Close()
    f_old_ratio.Close()
    f_old_direct.Close()

    print("\nDone.")


if __name__ == "__main__":
    main()
#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt


INPUT = (
    "/exp/minerva/app/users/qvuong/MAT_AL9/CC-NuE-XSec/"
    "flux_studies/beamfocus_fullstats_FHC.npz"
)

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

SHIFT_MIN = 0.0
SHIFT_MAX = 4.0
SHIFT_STEP = 0.05

# Regions used to estimate flat baselines.
MU_BASELINE_REGIONS = [
    (2.5, 6.5),
    (14.0, 20.0),
]

E_BASELINE_REGIONS = [
    (2.5, 6.0),
    (8.0, 10.0),
]

# Region over which to compare shifted focusing structures.
COMPARE_REGION = (5.0, 13.0)


def region_mask(centers, regions):
    mask = np.zeros(len(centers), dtype=bool)

    for low, high in regions:
        mask |= (
            (centers >= low)
            & (centers < high)
        )

    return mask


def main():

    data = np.load(INPUT)

    edges = data["le_bin_edges"]
    sigma = data["le_sigma_centered"]

    centers = 0.5 * (
        edges[:-1]
        + edges[1:]
    )

    # species:
    # 0 = numu
    # 1 = nue

    sigma_mu = 100.0 * sigma[0]
    sigma_e = 100.0 * sigma[1]

    # --------------------------------------------------------
    # Restrict to useful range
    # --------------------------------------------------------

    useful = (
        np.isfinite(sigma_mu)
        & np.isfinite(sigma_e)
        & (centers < 20.0)
    )

    # --------------------------------------------------------
    # Baselines
    # --------------------------------------------------------

    mu_base_mask = (
        region_mask(
            centers,
            MU_BASELINE_REGIONS,
        )
        & useful
    )

    e_base_mask = (
        region_mask(
            centers,
            E_BASELINE_REGIONS,
        )
        & useful
    )

    mu_baseline = np.mean(
        sigma_mu[mu_base_mask]
    )

    e_baseline = np.mean(
        sigma_e[e_base_mask]
    )

    print()
    print("Muon baseline     = {:.4f}%".format(
        mu_baseline
    ))

    print("Electron baseline = {:.4f}%".format(
        e_baseline
    ))

    # --------------------------------------------------------
    # Excess over flat baseline
    # --------------------------------------------------------

    excess_mu = sigma_mu - mu_baseline
    excess_e = sigma_e - e_baseline

    # For the focusing-feature comparison, negative deviations
    # below baseline are not interpreted as focusing excess.
    excess_mu_positive = np.maximum(
        excess_mu,
        0.0,
    )

    excess_e_positive = np.maximum(
        excess_e,
        0.0,
    )

    # --------------------------------------------------------
    # Scan energy shift
    #
    # Positive shift means:
    #
    #   electron feature at E
    #     compared to
    #   muon feature at E + shift
    #
    # i.e. move the muon peak toward lower energy.
    # --------------------------------------------------------

    shifts = np.arange(
        SHIFT_MIN,
        SHIFT_MAX + 0.5 * SHIFT_STEP,
        SHIFT_STEP,
    )

    scores = []

    compare_mask = (
        (centers >= COMPARE_REGION[0])
        & (centers < COMPARE_REGION[1])
        & useful
    )

    x_compare = centers[compare_mask]
    e_compare = excess_e_positive[
        compare_mask
    ]

    for shift in shifts:

        shifted_mu = np.interp(
            x_compare + shift,
            centers,
            excess_mu_positive,
            left=np.nan,
            right=np.nan,
        )

        valid = (
            np.isfinite(shifted_mu)
            & np.isfinite(e_compare)
        )

        if np.count_nonzero(valid) < 2:
            scores.append(np.nan)
            continue

        diff = (
            shifted_mu[valid]
            - e_compare[valid]
        )

        rmse = np.sqrt(
            np.mean(
                diff * diff
            )
        )

        scores.append(rmse)

    scores = np.array(scores)

    best_index = np.nanargmin(scores)
    best_shift = shifts[best_index]
    best_rmse = scores[best_index]

    print()
    print(
        "Best muon -> electron energy shift = "
        "{:.2f} GeV".format(
            best_shift
        )
    )

    print(
        "RMSE of excess-shape comparison    = "
        "{:.4f}%".format(
            best_rmse
        )
    )

    # --------------------------------------------------------
    # Build shifted muon focusing excess
    # --------------------------------------------------------

    shifted_mu_excess = np.interp(
        centers + best_shift,
        centers,
        excess_mu_positive,
        left=0.0,
        right=0.0,
    )

    predicted_e_sigma = (
        e_baseline
        + shifted_mu_excess
    )

    # --------------------------------------------------------
    # Print useful comparison
    # --------------------------------------------------------

    print()
    print(
        "{:>12s} "
        "{:>10s} "
        "{:>10s} "
        "{:>10s} "
        "{:>10s}".format(
            "E_true",
            "mu",
            "e",
            "mu excess",
            "e excess",
        )
    )

    for i in range(len(centers)):

        if centers[i] < 4.0:
            continue

        if centers[i] >= 15.0:
            break

        print(
            "{:5.1f}-{:5.1f} "
            "{:9.3f}% "
            "{:9.3f}% "
            "{:9.3f}% "
            "{:9.3f}%".format(
                edges[i],
                edges[i + 1],
                sigma_mu[i],
                sigma_e[i],
                excess_mu_positive[i],
                excess_e_positive[i],
            )
        )

    # ========================================================
    # Plot 1: original uncertainty curves
    # ========================================================

    fig, ax = plt.subplots(
        figsize=(9, 6)
    )

    ax.plot(
        centers,
        sigma_mu,
        marker="o",
        label=r"$\nu_\mu$",
    )

    ax.plot(
        centers,
        sigma_e,
        marker="o",
        label=r"$\nu_e$",
    )

    ax.axhline(
        mu_baseline,
        linestyle="--",
        linewidth=1,
        label=(
            r"$\nu_\mu$ baseline"
        ),
    )

    ax.axhline(
        e_baseline,
        linestyle="--",
        linewidth=1,
        label=(
            r"$\nu_e$ baseline"
        ),
    )

    ax.set_xlim(
        2.0,
        20.0,
    )

    ax.set_xlabel(
        r"True $E_\nu$ [GeV]"
    )

    ax.set_ylabel(
        "BeamFocus fractional uncertainty [%]"
    )

    ax.set_title(
        "FHC BeamFocus uncertainty: muon vs electron"
    )

    ax.grid(
        True,
        alpha=0.25,
    )

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        "beamfocus_FHC_mu_e_baselines.png",
        dpi=160,
    )

    plt.close(fig)

    # ========================================================
    # Plot 2: excess focusing structure
    # ========================================================

    fig, ax = plt.subplots(
        figsize=(9, 6)
    )

    ax.plot(
        centers,
        excess_e_positive,
        marker="o",
        label=r"$\nu_e$ excess",
    )

    ax.plot(
        centers,
        shifted_mu_excess,
        marker="o",
        label=(
            r"shifted $\nu_\mu$ excess "
            + "({:.2f} GeV)".format(
                best_shift
            )
        ),
    )

    ax.set_xlim(
        4.0,
        14.0,
    )

    ax.set_xlabel(
        r"True $E_\nu$ [GeV]"
    )

    ax.set_ylabel(
        "Excess BeamFocus uncertainty above baseline [%]"
    )

    ax.set_title(
        "FHC focusing-feature alignment"
    )

    ax.grid(
        True,
        alpha=0.25,
    )

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        "beamfocus_FHC_mu_to_e_peak_alignment.png",
        dpi=160,
    )

    plt.close(fig)

    # ========================================================
    # Plot 3: actual electron vs simple prediction
    # ========================================================

    fig, ax = plt.subplots(
        figsize=(9, 6)
    )

    ax.plot(
        centers,
        sigma_e,
        marker="o",
        label=r"measured $\nu_e$",
    )

    ax.plot(
        centers,
        predicted_e_sigma,
        marker="o",
        label=(
            r"flat baseline + shifted $\nu_\mu$ peak"
        ),
    )

    ax.set_xlim(
        2.0,
        14.0,
    )

    ax.set_xlabel(
        r"True $E_\nu$ [GeV]"
    )

    ax.set_ylabel(
        "BeamFocus fractional uncertainty [%]"
    )

    ax.set_title(
        "Simple electron BeamFocus model"
    )

    ax.grid(
        True,
        alpha=0.25,
    )

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        "beamfocus_FHC_electron_shift_model.png",
        dpi=160,
    )

    plt.close(fig)

    # ========================================================
    # Plot 4: scan score
    # ========================================================

    fig, ax = plt.subplots(
        figsize=(8, 5)
    )

    ax.plot(
        shifts,
        scores,
        marker="o",
    )

    ax.axvline(
        best_shift,
        linestyle="--",
    )

    ax.set_xlabel(
        r"Muon-to-electron shift [GeV]"
    )

    ax.set_ylabel(
        "RMSE of focusing excess [%]"
    )

    ax.set_title(
        "Muon-to-electron focusing-peak shift scan"
    )

    ax.grid(
        True,
        alpha=0.25,
    )

    fig.tight_layout()

    fig.savefig(
        "beamfocus_FHC_mu_to_e_shift_scan.png",
        dpi=160,
    )

    plt.close(fig)

    print()
    print("Saved:")
    print("  beamfocus_FHC_mu_e_baselines.png")
    print("  beamfocus_FHC_mu_to_e_peak_alignment.png")
    print("  beamfocus_FHC_electron_shift_model.png")
    print("  beamfocus_FHC_mu_to_e_shift_scan.png")


if __name__ == "__main__":
    main()
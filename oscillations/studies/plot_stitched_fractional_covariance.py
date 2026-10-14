#!/usr/bin/env python3

import os
import argparse
import json

import ROOT
import PlotUtils
import numpy as np
import matplotlib.pyplot as plt


ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)


SAMPLE_LABELS = {
    "fhc_elastic": r"$\nu+e$",
    "fhc_ratio": r"FHC $\mu/e$",
    "rhc_ratio": r"RHC $\mu/e$",
    "fhc_numu": r"FHC $\nu_\mu$",
    "rhc_numu": r"RHC $\bar{\nu}_\mu$",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot fractional covariance and correlation matrices for grouped "
            "systematic uncertainties in a stitched MINERvA histogram."
        )
    )

    parser.add_argument(
        "input_file",
        help="Stitched ROOT file, e.g. NuE_stitched_hists_predNueel_p8.root",
    )

    parser.add_argument(
        "--hist",
        default="sample_mc",
        help="Stitched MnvH1D name [default: sample_mc]",
    )

    parser.add_argument(
        "--hist-config",
        default=None,
        help="Optional HIST_CONFIG_*.json used to draw sample boundaries and labels.",
    )

    parser.add_argument(
        "--outdir",
        default="plots/stitched_covariance",
    )

    parser.add_argument(
        "--print-bands",
        action="store_true",
        help="Print all error bands stored on the stitched histogram.",
    )

    return parser.parse_args()


def tmatrix_to_numpy(matrix):
    nrows = matrix.GetNrows()
    ncols = matrix.GetNcols()

    out = np.zeros((nrows, ncols), dtype=float)

    for i in range(nrows):
        for j in range(ncols):
            out[i, j] = matrix[i][j]

    return out


def root_string_vector_to_list(vec):
    return [str(x) for x in vec]


def get_available_error_bands(hist):
    vert = root_string_vector_to_list(hist.GetVertErrorBandNames())
    lat = root_string_vector_to_list(hist.GetLatErrorBandNames())

    return sorted(set(vert + lat))


def get_cv(hist):
    return np.array(
        [
            hist.GetBinContent(i)
            for i in range(1, hist.GetNbinsX() + 1)
        ],
        dtype=float,
    )


def get_sys_covariance(hist, band_name):
    """
    Return the absolute covariance matrix for one MnvH1D systematic error band.

    GetSysErrorMatrix includes ROOT underflow/overflow bins, so remove them
    to match the physical stitched analysis bins.
    """

    matrix = hist.GetSysErrorMatrix(band_name)
    cov = tmatrix_to_numpy(matrix)

    cov = cov[1:-1, 1:-1]

    return cov


def sum_covariances(hist, band_names, available):
    nbins = hist.GetNbinsX()
    total = np.zeros((nbins, nbins), dtype=float)

    used = []
    missing = []

    for name in band_names:
        if name not in available:
            missing.append(name)
            continue

        cov = get_sys_covariance(hist, name)

        if cov.shape != total.shape:
            raise RuntimeError(
                "Covariance for '{}' has shape {}, expected {}".format(
                    name,
                    cov.shape,
                    total.shape,
                )
            )

        total += cov
        used.append(name)

    return total, used, missing


def covariance_to_fractional(cov, cv):
    """
    Convert absolute covariance to fractional covariance:

        F_ij = V_ij / (N_i^CV N_j^CV)

    For bins with zero CV prediction, the corresponding row/column is set to zero.
    """

    denom = np.outer(cv, cv)

    frac_cov = np.divide(
        cov,
        denom,
        out=np.zeros_like(cov),
        where=(denom != 0.0),
    )

    return frac_cov


def covariance_to_correlation(cov):
    """
    Convert covariance to correlation:

        rho_ij = V_ij / sqrt(V_ii V_jj)
    """

    sigma = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    denom = np.outer(sigma, sigma)

    corr = np.divide(
        cov,
        denom,
        out=np.zeros_like(cov),
        where=(denom > 0.0),
    )

    return corr


def check_grouping(groups, available):
    available_set = set(available)

    assigned_to = {}

    for group_name, bands in groups.items():
        for band in bands:
            assigned_to.setdefault(band, []).append(group_name)

    duplicates = {
        band: names
        for band, names in assigned_to.items()
        if len(names) > 1
    }

    grouped_set = set(assigned_to)

    unassigned = sorted(available_set - grouped_set)
    missing = sorted(grouped_set - available_set)

    print("")
    print("===== systematic grouping check =====")
    print("available bands =", len(available_set))
    print("assigned bands  =", len(available_set & grouped_set))

    if duplicates:
        print("BANDS ASSIGNED TO MULTIPLE GROUPS:")
        for band, names in duplicates.items():
            print("  {}: {}".format(band, names))

    if unassigned:
        print("UNASSIGNED BANDS:")
        for band in unassigned:
            print("  ", band)

    if missing:
        print("GROUP DEFINITIONS NOT PRESENT IN FILE:")
        for band in missing:
            print("  ", band)

    if not duplicates and not unassigned and not missing:
        print("All available systematic bands are assigned exactly once.")

    print("===== end systematic grouping check =====")

    if duplicates or unassigned or missing:
        raise RuntimeError(
            "Systematic grouping is incomplete, contains overlaps, "
            "or requests missing bands."
        )


def get_group_definitions():
    flux = [
        "Flux",
    ]

    cross_section = [
        # Standard GENIE interaction uncertainties
        "GENIE_AGKYxF1pi",
        "GENIE_AhtBY",
        "GENIE_BhtBY",
        "GENIE_CCQEPauliSupViaKF",
        "GENIE_CV1uBY",
        "GENIE_CV2uBY",
        "GENIE_EP_MvRES",
        "GENIE_EtaNCEL",
        "GENIE_MaNCEL",
        "GENIE_MaRES",
        "GENIE_MaZExpCCQE",
        "GENIE_NormCCRES",
        "GENIE_NormDISCC",
        "GENIE_NormNCRES",
        "GENIE_RDecBR1gamma",
        "GENIE_Rvn1pi",
        "GENIE_Rvn2pi",
        "GENIE_Rvp1pi",
        "GENIE_Rvp2pi",
        "GENIE_Theta_Delta2Npi",
        "GENIE_VecFFCCQEshape",

        # GENIE final-state interaction uncertainties
        "GENIE_FrAbs_N",
        "GENIE_FrAbs_pi",
        "GENIE_FrCEx_N",
        "GENIE_FrCEx_pi",
        "GENIE_FrElas_N",
        "GENIE_FrElas_pi",
        "GENIE_FrInel_N",
        "GENIE_FrPiProd_N",
        "GENIE_FrPiProd_pi",
        "GENIE_MFP_N",
        "GENIE_MFP_pi",

        # MINERvA-specific interaction-model variations
        "RPA_HighQ2",
        "RPA_LowQ2",
        "Low_Recoil_2p2h_Tune",
        "fsi_weight",
        "SuSA_Valencia_Weight",
        "MK_model",
    ]

    detector = [
        # Beam direction
        "BeamAngleX",
        "BeamAngleY",

        # Electron reconstruction
        "eltheta",
        "elE_ECAL",
        "elE_HCAL",
        "electron_scale",
        "Leakage_Uncertainty",

        # Hadronic detector response
        "response_p",
        "response_meson",
        "response_em",
        "response_other",
        "response_xtalk",

        # Fiducial detector mass
        "Target_Mass_CH",

        # Muon reconstruction
        "Muon_Energy_MINERvA",
        "Muon_Energy_MINOS",
        "Muon_Energy_Resolution",
        "MuonAngleXResolution",
        "MuonAngleYResolution",
        "MINOS_Reconstruction_Efficiency",
    ]

    return {
        "flux": flux,
        "cross_section": cross_section,
        "detector": detector,
    }


def load_sample_layout(hist_config):
    if hist_config is None:
        return None

    with open(hist_config, "r") as f:
        cfg = json.load(f)

    samples = []

    for name, info in cfg.items():
        if "start" not in info or "end" not in info:
            continue

        samples.append(
            {
                "name": name,
                "start": int(info["start"]),
                "end": int(info["end"]),
            }
        )

    samples.sort(key=lambda x: x["start"])

    return samples


def add_sample_boundaries(ax, samples):
    if not samples:
        return

    for sample in samples[:-1]:
        boundary = sample["end"] + 1

        ax.axvline(
            boundary - 0.5,
            linewidth=1.2,
            color="black",
        )

        ax.axhline(
            boundary - 0.5,
            linewidth=1.2,
            color="black",
        )

    centers = [
        0.5 * (sample["start"] + sample["end"])
        for sample in samples
    ]

    labels = [
        SAMPLE_LABELS.get(sample["name"], sample["name"])
        for sample in samples
    ]

    ax.set_xticks(centers)
    ax.set_yticks(centers)

    ax.set_xticklabels(
        labels,
        rotation=45,
        ha="right",
        fontsize=8,
    )

    ax.set_yticklabels(
        labels,
        fontsize=8,
    )


def plot_matrix(
    matrix,
    title,
    outpath,
    samples=None,
    matrix_type="fractional_covariance",
):
    fig, ax = plt.subplots(figsize=(9, 8))

    if matrix_type == "correlation":
        im = ax.imshow(
            matrix,
            origin="upper",
            interpolation="nearest",
            vmin=-1.0,
            vmax=1.0,
            aspect="equal",
            cmap="RdBu_r",
        )

        cbar_label = "Correlation coefficient"

    elif matrix_type == "fractional_covariance":
        max_abs = np.max(np.abs(matrix))

        if max_abs > 0:
            im = ax.imshow(
                matrix,
                origin="upper",
                interpolation="nearest",
                vmin=-max_abs,
                vmax=max_abs,
                aspect="equal",
                cmap="RdBu_r",
            )
        else:
            im = ax.imshow(
                matrix,
                origin="upper",
                interpolation="nearest",
                aspect="equal",
                cmap="RdBu_r",
            )

        cbar_label = "Fractional covariance"

    else:
        raise ValueError(
            "Unknown matrix_type '{}'".format(matrix_type)
        )

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(cbar_label)

    ax.set_title(title)
    ax.set_xlabel("Stitched analysis bin")
    ax.set_ylabel("Stitched analysis bin")

    add_sample_boundaries(ax, samples)

    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)

    print("Wrote", outpath)


def main():
    args = parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    root_file = ROOT.TFile.Open(args.input_file)

    if not root_file or root_file.IsZombie():
        raise RuntimeError(
            "Could not open {}".format(args.input_file)
        )

    hist = root_file.Get(args.hist)

    if not hist:
        raise RuntimeError(
            "Could not find histogram '{}' in {}".format(
                args.hist,
                args.input_file,
            )
        )

    available = get_available_error_bands(hist)
    cv = get_cv(hist)

    print("")
    print("Input file :", args.input_file)
    print("Histogram  :", args.hist)
    print("Bins       :", hist.GetNbinsX())
    print("Error bands:", len(available))
    print("Zero-CV bins:", np.where(cv == 0.0)[0].tolist())

    if args.print_bands:
        print("")
        print("===== available error bands =====")

        for name in available:
            print(name)

        print("===== end error bands =====")

    samples = load_sample_layout(args.hist_config)
    groups = get_group_definitions()

    check_grouping(groups, available)

    print("")
    print("===== grouped covariance construction =====")

    for group_name, requested_bands in groups.items():
        cov, used, missing = sum_covariances(
            hist,
            requested_bands,
            available,
        )

        frac_cov = covariance_to_fractional(
            cov,
            cv,
        )

        corr = covariance_to_correlation(
            cov,
        )

        print("")
        print(group_name)

        print("  used bands:")

        for name in used:
            print("   ", name)

        if missing:
            print("  requested but absent:")

            for name in missing:
                print("   ", name)

        print(
            "  trace(abs cov) =",
            float(np.trace(cov)),
        )

        print(
            "  max |fractional cov| =",
            float(np.max(np.abs(frac_cov))),
        )

        print(
            "  max |correlation| =",
            float(np.max(np.abs(corr))),
        )

        # Save numerical matrices.
        np.save(
            os.path.join(
                args.outdir,
                "{}_covariance.npy".format(group_name),
            ),
            cov,
        )

        np.save(
            os.path.join(
                args.outdir,
                "{}_fractional_covariance.npy".format(group_name),
            ),
            frac_cov,
        )

        np.save(
            os.path.join(
                args.outdir,
                "{}_correlation.npy".format(group_name),
            ),
            corr,
        )

        # Plot fractional covariance.
        plot_matrix(
            frac_cov,
            "{} fractional covariance".format(
                group_name.replace("_", " ").title()
            ),
            os.path.join(
                args.outdir,
                "{}_fractional_covariance.png".format(group_name),
            ),
            samples=samples,
            matrix_type="fractional_covariance",
        )

        # Plot correlation.
        plot_matrix(
            corr,
            "{} correlation".format(
                group_name.replace("_", " ").title()
            ),
            os.path.join(
                args.outdir,
                "{}_correlation.png".format(group_name),
            ),
            samples=samples,
            matrix_type="correlation",
        )

    print("")
    print("===== end grouped covariance construction =====")

    root_file.Close()


if __name__ == "__main__":
    main()
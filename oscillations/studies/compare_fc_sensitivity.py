import os
import numpy as np

INPUT_DIR = "/exp/minerva/data/users/qvuong/surfaces/npys/p8"

MODES = {
    "ratio": "prodNueel_p8_profiledFlux_includeAll_Nprof30",
    "direct": "prodNueel_noRatio_p8_profiledFlux_includeAll_Nprof35",
}

# Must match makeSurface.py exactly.
UMU4 = 0.41 * np.logspace(-5, 0, 100)
UMU4[0] = 0.0

UE4 = 0.15 * np.logspace(-4, 0, 100)
UE4[0] = 0.0


def load_mode(mode):
    asimov_file = os.path.join(
        INPUT_DIR,
        "asimov_chi2s_{}.npy".format(mode),
    )

    fc_file = os.path.join(
        INPUT_DIR,
        "asimov_deltachi2s_{}.npy".format(mode),
    )

    dm_file = os.path.join(
        INPUT_DIR,
        "delta_m_values_{}.npy".format(mode),
    )

    asimov = np.load(asimov_file)
    toys = np.load(fc_file)
    dm2_values = np.load(dm_file)

    asimov_dchi2 = asimov - np.nanmin(asimov)
    critical95 = np.percentile(toys, 95.0)

    print("\nMode:", mode)
    print("  shape =", asimov.shape)
    print("  Asimov grid minimum =", np.nanmin(asimov))
    print("  FC 95% critical =", critical95)
    print("  toys =", toys.size)

    return {
        "asimov": asimov,
        "dchi2": asimov_dchi2,
        "critical95": critical95,
        "dm2": dm2_values,
    }


# ============================================================
# Part 1: old axis-crossing comparison
# ============================================================

def first_crossing(axis, values, critical):
    """
    Return the first axis value where dchi2 >= critical.

    Assumes movement from small mixing toward large mixing.
    Returns np.nan if threshold is never crossed.
    """

    above = np.where(values >= critical)[0]

    if len(above) == 0:
        return np.nan

    return axis[above[0]]


def safe_ratio(a, b):
    """
    Return a/b only when both are finite and denominator is nonzero.
    """

    if not np.isfinite(a):
        return np.nan

    if not np.isfinite(b):
        return np.nan

    if b == 0.0:
        return np.nan

    return a / b


def compare_fixed_ue4(results, dm2_indices, ue4_indices):
    print("\n")
    print("============================================")
    print("95% sensitivity: minimum Umu4 at fixed Ue4")
    print("============================================")

    for idm in dm2_indices:
        dm2 = results["ratio"]["dm2"][idm]

        print("\ndm2 = {:.5g}".format(dm2))

        for ie in ue4_indices:
            ue = UE4[ie]

            vals_ratio = results["ratio"]["dchi2"][idm, :, ie]
            vals_direct = results["direct"]["dchi2"][idm, :, ie]

            lim_ratio = first_crossing(
                UMU4,
                vals_ratio,
                results["ratio"]["critical95"],
            )

            lim_direct = first_crossing(
                UMU4,
                vals_direct,
                results["direct"]["critical95"],
            )

            comparison = safe_ratio(lim_ratio, lim_direct)

            print(
                "  Ue4={:9.5g}   "
                "ratio={:9.5g}   "
                "direct={:9.5g}   "
                "ratio/direct={:7.3f}".format(
                    ue,
                    lim_ratio,
                    lim_direct,
                    comparison,
                )
            )


def compare_fixed_umu4(results, dm2_indices, umu4_indices):
    print("\n")
    print("============================================")
    print("95% sensitivity: minimum Ue4 at fixed Umu4")
    print("============================================")

    for idm in dm2_indices:
        dm2 = results["ratio"]["dm2"][idm]

        print("\ndm2 = {:.5g}".format(dm2))

        for imu in umu4_indices:
            umu = UMU4[imu]

            vals_ratio = results["ratio"]["dchi2"][idm, imu, :]
            vals_direct = results["direct"]["dchi2"][idm, imu, :]

            lim_ratio = first_crossing(
                UE4,
                vals_ratio,
                results["ratio"]["critical95"],
            )

            lim_direct = first_crossing(
                UE4,
                vals_direct,
                results["direct"]["critical95"],
            )

            comparison = safe_ratio(lim_ratio, lim_direct)

            print(
                "  Umu4={:9.5g}   "
                "ratio={:9.5g}   "
                "direct={:9.5g}   "
                "ratio/direct={:7.3f}".format(
                    umu,
                    lim_ratio,
                    lim_direct,
                    comparison,
                )
            )


# ============================================================
# Part 2: 2D excluded-area / overlap comparison
# ============================================================

def build_log_cell_weights():
    """
    Construct area weights in:

        dlog10(Umu4) dlog10(Ue4)

    Zero-mixing row and column are excluded because log10(0)
    is undefined.

    Returned shape:
        (99, 99)

    corresponding to:
        UMU4[1:] x UE4[1:]
    """

    log_umu = np.log10(UMU4[1:])
    log_ue = np.log10(UE4[1:])

    umu_edges = np.empty(len(log_umu) + 1)
    ue_edges = np.empty(len(log_ue) + 1)

    umu_edges[1:-1] = (
        0.5 * (log_umu[:-1] + log_umu[1:])
    )

    ue_edges[1:-1] = (
        0.5 * (log_ue[:-1] + log_ue[1:])
    )

    umu_edges[0] = (
        log_umu[0]
        - 0.5 * (log_umu[1] - log_umu[0])
    )

    umu_edges[-1] = (
        log_umu[-1]
        + 0.5 * (log_umu[-1] - log_umu[-2])
    )

    ue_edges[0] = (
        log_ue[0]
        - 0.5 * (log_ue[1] - log_ue[0])
    )

    ue_edges[-1] = (
        log_ue[-1]
        + 0.5 * (log_ue[-1] - log_ue[-2])
    )

    dlog_umu = np.diff(umu_edges)
    dlog_ue = np.diff(ue_edges)

    return np.outer(dlog_umu, dlog_ue)


def compare_excluded_regions(results):
    weights = build_log_cell_weights()

    ratio = results["ratio"]
    direct = results["direct"]

    dm2_values = ratio["dm2"]

    if not np.allclose(dm2_values, direct["dm2"]):
        raise RuntimeError(
            "Ratio and direct configurations use different dm2 grids."
        )

    rows = []

    print("\n")
    print(
        "============================================================================"
    )
    print(
        "95% FC expected exclusion comparison in log10(Umu4)-log10(Ue4) area"
    )
    print(
        "============================================================================"
    )

    print(
        "{:>10s}  {:>9s}  {:>9s}  {:>9s}  {:>10s}  {:>10s}".format(
            "dm2",
            "ratio",
            "direct",
            "overlap",
            "ratioOnly",
            "directOnly",
        )
    )

    for idm, dm2 in enumerate(dm2_values):

        # Exclude zero-mixing row/column for log-space area.
        dr = ratio["dchi2"][idm, 1:, 1:]
        dd = direct["dchi2"][idm, 1:, 1:]

        finite = np.isfinite(dr) & np.isfinite(dd)

        ratio_excl = (
            dr >= ratio["critical95"]
        ) & finite

        direct_excl = (
            dd >= direct["critical95"]
        ) & finite

        overlap = ratio_excl & direct_excl
        ratio_only = ratio_excl & (~direct_excl)
        direct_only = direct_excl & (~ratio_excl)

        valid_weights = np.where(
            finite,
            weights,
            0.0,
        )

        total_area = np.sum(valid_weights)

        if total_area == 0.0:
            f_ratio = np.nan
            f_direct = np.nan
            f_overlap = np.nan
            f_ratio_only = np.nan
            f_direct_only = np.nan

        else:
            f_ratio = (
                np.sum(valid_weights * ratio_excl)
                / total_area
            )

            f_direct = (
                np.sum(valid_weights * direct_excl)
                / total_area
            )

            f_overlap = (
                np.sum(valid_weights * overlap)
                / total_area
            )

            f_ratio_only = (
                np.sum(valid_weights * ratio_only)
                / total_area
            )

            f_direct_only = (
                np.sum(valid_weights * direct_only)
                / total_area
            )

        rows.append(
            [
                dm2,
                f_ratio,
                f_direct,
                f_overlap,
                f_ratio_only,
                f_direct_only,
            ]
        )

        print(
            "{:10.4g}  {:8.2f}%  {:8.2f}%  {:8.2f}%  {:9.2f}%  {:9.2f}%".format(
                dm2,
                100.0 * f_ratio,
                100.0 * f_direct,
                100.0 * f_overlap,
                100.0 * f_ratio_only,
                100.0 * f_direct_only,
            )
        )

    return np.asarray(rows)


# ============================================================
# Part 3: global summary
# ============================================================

def summarize_sensitive_region(rows):
    f_ratio = rows[:, 1]
    f_direct = rows[:, 2]
    overlap = rows[:, 3]
    ratio_only = rows[:, 4]
    direct_only = rows[:, 5]

    sensitive = (
        np.isfinite(f_ratio)
        & np.isfinite(f_direct)
        & (
            (f_ratio > 0.0)
            | (f_direct > 0.0)
        )
    )

    if not np.any(sensitive):
        print("\nNo dm2 slices have 95% expected exclusion.")
        return

    diff = (
        f_ratio[sensitive]
        - f_direct[sensitive]
    )

    print("\n")
    print("============================================")
    print("Summary over sensitive dm2 slices")
    print("============================================")

    print(
        "Number of sensitive dm2 slices =",
        np.sum(sensitive),
    )

    print(
        "\nMean excluded area:"
    )

    print(
        "  ratio  = {:.2f}%".format(
            100.0 * np.mean(f_ratio[sensitive])
        )
    )

    print(
        "  direct = {:.2f}%".format(
            100.0 * np.mean(f_direct[sensitive])
        )
    )

    print(
        "\nMean overlap area = {:.2f}%".format(
            100.0 * np.mean(overlap[sensitive])
        )
    )

    print(
        "Mean ratio-only area = {:.2f}%".format(
            100.0 * np.mean(ratio_only[sensitive])
        )
    )

    print(
        "Mean direct-only area = {:.2f}%".format(
            100.0 * np.mean(direct_only[sensitive])
        )
    )

    print(
        "\nMean ratio - direct excluded area = "
        "{:+.2f} percentage points".format(
            100.0 * np.mean(diff)
        )
    )

    print(
        "Maximum |ratio - direct| = "
        "{:.2f} percentage points".format(
            100.0 * np.max(np.abs(diff))
        )
    )


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":

    results = {}

    for label, mode in MODES.items():
        results[label] = load_mode(mode)

    # --------------------------------------------------------
    # Old axis-crossing diagnostics
    # --------------------------------------------------------

    dm2_indices = [
        20,
        32,
        44,
        51,
        61,
        71,
    ]

    ue4_indices = [
        53,
        70,
        78,
        84,
        92,
    ]

    umu4_indices = [
        50,
        58,
        65,
        72,
        87,
    ]

    compare_fixed_ue4(
        results,
        dm2_indices,
        ue4_indices,
    )

    compare_fixed_umu4(
        results,
        dm2_indices,
        umu4_indices,
    )

    # --------------------------------------------------------
    # New full-2D excluded-area comparison
    # --------------------------------------------------------

    rows = compare_excluded_regions(results)

    summarize_sensitive_region(rows)

    outfile = "fc_expected_area_comparison.csv"

    header = (
        "dm2,"
        "ratio_excluded_fraction,"
        "direct_excluded_fraction,"
        "overlap_fraction,"
        "ratio_only_fraction,"
        "direct_only_fraction"
    )

    np.savetxt(
        outfile,
        rows,
        delimiter=",",
        header=header,
        comments="",
    )

    print("\nSaved:", outfile)
#!/usr/bin/env python3

import os
import sys
import argparse
import csv
import shutil
import json
import numpy as np
import logging

def parse_scan_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--hist-config-tag", default=None)
    parser.add_argument("--exclude", default=None)
    parser.add_argument("--mask", default="none")
    parser.add_argument("--out", default=None)

    parser.add_argument("--n-toys", type=int, default=1)
    parser.add_argument("--seed", type=int, default=12345)

    parser.add_argument(
        "--toy-mode",
        default="gauss_poisson",
        choices=["asimov", "poisson", "gauss", "gauss_poisson"],
        help=(
            "asimov: data = MC CV; "
            "poisson: Poisson(MC CV); "
            "gauss: MC CV + Gaussian covariance fluctuation; "
            "gauss_poisson: Gaussian covariance fluctuation then Poisson throw"
        ),
    )

    parser.add_argument(
        "--cov-mode",
        default="full",
        choices=["full", "sansFlux"],
        help="Covariance used for Gaussian toy fluctuation.",
    )

    parser.add_argument(
        "--lambdas",
        default="0.01,0.03,0.1,0.2,0.3,0.5,0.7,0.85,1,1.15,1.3,1.5,2,3,5,10,30,100",
        help="Comma-separated lambda values",
    )

    args, remaining = parser.parse_known_args()

    # Hide lambda-scan-only args from AnalysisConfig/project parsers.
    sys.argv = [sys.argv[0]] + remaining

    return args


args = parse_scan_args()

import ROOT
import PlotUtils

from config.AnalysisConfig import AnalysisConfig
from tools.Helper import TMatrix_to_Numpy
from tools.StitchedHistogram import StitchedHistogram
from tools.Fitters import Statistics, OscillationFitter

ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

def throw_multivariate_psd(mean, cov, size=1, tol=1e-10):
    """
    Draw from N(mean, cov) for a symmetric positive-semidefinite covariance.
    Allows exact zero eigenvalues, unlike Cholesky.

    Returns:
      size == 1 : shape (nbins,)
      size > 1  : shape (size, nbins)
    """
    mean = np.asarray(mean, dtype=float)
    cov = np.asarray(cov, dtype=float)

    # Force exact symmetry.
    cov = 0.5 * (cov + cov.T)

    eigvals, eigvecs = np.linalg.eigh(cov)

    min_eig = np.min(eigvals)
    if min_eig < -tol:
        print("WARNING: covariance has negative eigenvalue:", min_eig)

    # Clip tiny negative numerical noise.
    eigvals = np.clip(eigvals, 0.0, None)

    z = np.random.normal(size=(size, len(mean)))
    sqrt_cov = eigvecs @ np.diag(np.sqrt(eigvals))

    throws = mean + z @ sqrt_cov.T

    if size == 1:
        return throws[0]

    return throws


def ThrowSystematics(
    histogram,
    throwFlux=False,
    useDataSubMCCov=True,
    n_samples=50,
    doDiagnostics=False,
):
    pred_vals = np.array(histogram.GetMCHistogram())[1:-1]

    # Full covariance currently includes stat + all systematics.
    V_full = histogram.GetCovarianceMatrix(False)
    V_sansFlux = histogram.GetCovarianceMatrix(sansFlux=True)
    V_flux = V_full - V_sansFlux

    # Get the stat covariance from the same residual histogram convention
    # used in StitchedHistogram.SetCovarianceMatrices().
    h_test = histogram.GetDataHistogram().Clone()
    h_test.Add(histogram.GetMCHistogram(), -1)

    V_stat = TMatrix_to_Numpy(h_test.GetStatErrorMatrix())[1:-1, 1:-1]

    # Remove stat from the Gaussian throw covariance.
    # Poisson will handle statistical fluctuations later.
    V_full_systOnly = V_full - V_stat
    V_sansFlux_systOnly = V_sansFlux - V_stat

    # Force symmetry.
    V_flux = 0.5 * (V_flux + V_flux.T)
    V_full_systOnly = 0.5 * (V_full_systOnly + V_full_systOnly.T)
    V_sansFlux_systOnly = 0.5 * (V_sansFlux_systOnly + V_sansFlux_systOnly.T)

    if doDiagnostics:
        print_cov_diagnostics("V_full = stat + all syst", V_full)
        print_cov_diagnostics("V_sansFlux = stat + nonflux syst", V_sansFlux)
        print_cov_diagnostics("V_stat", V_stat)
        print_cov_diagnostics("V_flux", V_flux)
        print_cov_diagnostics("V_full_systOnly", V_full_systOnly)
        print_cov_diagnostics("V_sansFlux_systOnly", V_sansFlux_systOnly)

    if throwFlux:
        # Same strategy as Ryan's script:
        # one independent flux throw per toy.
        flux_throws = throw_multivariate_psd(
            pred_vals,
            V_flux,
            size=n_samples,
        )

        flux_throws = np.asarray(flux_throws)
        if flux_throws.ndim == 1:
            flux_throws = flux_throws.reshape(1, -1)

        # Then one independent non-flux systematic throw per toy,
        # centered on that toy's flux-shifted mean.
        sys_throws = []
        for flux_mean in flux_throws:
            toy_mean = throw_multivariate_psd(
                flux_mean,
                V_sansFlux_systOnly,
                size=1,
            )
            sys_throws.append(toy_mean)

        sys_throws = np.asarray(sys_throws)

    else:
        # No explicit flux profiling:
        # throw all systematics together, but still no stat in the Gaussian covariance.
        sys_throws = throw_multivariate_psd(
            pred_vals,
            V_full_systOnly,
            size=n_samples,
        )

    sys_throws = np.asarray(sys_throws)
    if sys_throws.ndim == 1:
        sys_throws = sys_throws.reshape(1, -1)

    print("ThrowSystematics shape:", sys_throws.shape)
    return sys_throws

def ThrowPoissons(lambdas, histogram, throwFlux=False):
    lambdas = np.asarray(lambdas)

    # Force shape = (n_toys, n_bins)
    if lambdas.ndim == 1:
        lambdas = lambdas.reshape(1, -1)

    print("ThrowPoissons input shape:", lambdas.shape)

    throws = []
    for lam in lambdas:
        while lam[lam < 0].any():
            logging.error("Negative value in sys throws, rerunning this throw...")
            lam = ThrowSystematics(
                histogram,
                throwFlux=throwFlux,
                n_samples=1,
                doDiagnostics=False,
            )[0]

        throw = np.random.poisson(lam)
        throws.append(throw)

    throws = np.asarray(throws)

    if throws.ndim == 1:
        throws = throws.reshape(1, -1)

    print("ThrowPoissons output shape:", throws.shape)

    return throws

def write_rows_csv(rows, out_csv):
    if len(rows) == 0:
        return

    fieldnames = list(rows[0].keys())

    with open(out_csv, "w") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_mask(mask_name):
    if mask_name in [None, "", "none", "None"]:
        return None

    if mask_name == "ratio":
        return {
            "fhc_ratio": list(range(1, 13)),
            "rhc_ratio": list(range(1, 13)),
        }

    if mask_name == "bin1":
        return {
            "fhc_ratio": [1],
            "rhc_ratio": [1],
        }

    if mask_name == "bins9_10":
        return {
            "fhc_ratio": [9, 10],
            "rhc_ratio": [9, 10],
        }

    if mask_name == "bins1_9_10":
        return {
            "fhc_ratio": [1, 9, 10],
            "rhc_ratio": [1, 9, 10],
        }

    if mask_name == "fhc_bins1_9_10":
        return {
            "fhc_ratio": [1, 9, 10],
        }

    if mask_name == "rhc_bins1_9_10":
        return {
            "rhc_ratio": [1, 9, 10],
        }

    raise ValueError("Unknown mask_name: {}".format(mask_name))


def psd_gaussian_shift(cov, rng):
    """
    Draw one correlated Gaussian fluctuation from covariance using PSD eigen sqrt.
    Handles tiny negative eigenvalues from numerical noise.
    """
    cov = 0.5 * (cov + cov.T)

    vals, vecs = np.linalg.eigh(cov)

    min_eval = np.min(vals)
    if min_eval < -1e-6 * max(np.max(vals), 1.0):
        print("WARNING: covariance has significantly negative eigenvalue:", min_eval)

    vals = np.clip(vals, 0.0, None)

    z = rng.normal(size=len(vals))
    shift = vecs @ (np.sqrt(vals) * z)

    return shift


def generate_pseudo_vector(base_histogram, rng, toy_mode, cov_mode):
    """
    Generate one pseudo-data vector in stitched global-bin order.

    toy_mode:
        asimov        -> MC CV
        poisson       -> Poisson(MC CV)
        gauss         -> MC CV + Gaussian cov fluctuation
        gauss_poisson -> Gaussian cov fluctuation, then Poisson
    """
    base = np.array(base_histogram.GetMCHistogram())[1:-1].astype(float)

    mean = base.copy()

    if toy_mode in ["gauss", "gauss_poisson"]:
        if cov_mode == "full":
            cov = base_histogram.GetCovarianceMatrix(sansFlux=False)
        elif cov_mode == "sansFlux":
            cov = base_histogram.GetCovarianceMatrix(sansFlux=True)
        else:
            raise ValueError("Unknown cov_mode: {}".format(cov_mode))

        shift = psd_gaussian_shift(cov, rng)
        mean = mean + shift

    mean = np.clip(mean, 0.0, None)

    if toy_mode in ["poisson", "gauss_poisson"]:
        pseudo = rng.poisson(mean)
    elif toy_mode in ["asimov", "gauss"]:
        pseudo = mean
    else:
        raise ValueError("Unknown toy_mode: {}".format(toy_mode))

    return pseudo.astype(float)


def overwrite_data_hist_from_stitched_vector(histogram, stitched_vector):
    """
    Replace the stitched data histogram with a pseudo-data vector.

    This follows the same strategy as fitAsimovs_quinn.py:
    start from the existing stitched data histogram, build weights so that
    DivideSingle changes the bin contents to the toy values, then call
    SetDataHistogram().
    """
    stitched_data = histogram.GetDataHistogram()

    toy = np.asarray(stitched_vector, dtype=float)
    if toy.ndim != 1:
        raise RuntimeError("stitched_vector must be 1D")

    if stitched_data.GetNbinsX() != len(toy):
        raise RuntimeError(
            "Toy vector length {} does not match stitched data bins {}".format(
                len(toy),
                stitched_data.GetNbinsX()
            )
        )

    weights = stitched_data.Clone().GetCVHistoWithStatError()

    for i in range(1, weights.GetNbinsX() + 1):
        if toy[i - 1] != 0:
            weight = stitched_data.GetBinContent(i) / toy[i - 1]
        else:
            weight = stitched_data.GetBinContent(i)

        weights.SetBinContent(i, weight)
        weights.SetBinError(i, 0)

    data_histogram = stitched_data.Clone()
    data_histogram.DivideSingle(data_histogram, weights)

    histogram.SetDataHistogram(data_histogram)


def get_residual_and_penalty(stat, useOsc=False):
    chi2, penalty = stat.Chi2DataMC(
        marginalize=True,
        useOsc=useOsc,
    )

    flux_fitter = stat.GetFluxFitter(useOsc=useOsc)
    a = flux_fitter.GetFluxSolution()

    residual = chi2 - penalty
    norm_a = np.linalg.norm(a)
    max_abs_a = np.max(np.abs(a))

    return chi2, residual, penalty, norm_a, max_abs_a


def run_one_lambda(sample_histogram, lam, exclude_samples, mask_spec):
    """
    Run null profiled chi2 and oscillation best fit for one lambda value.
    This uses whatever is currently stored in sample_histogram.data_hists.
    """

    stat_null = Statistics(
        sample_histogram,
        exclude=exclude_samples,
        lam=lam,
        mask_spec=mask_spec,
    )

    chi2_null, resid_null, pen_null, norm_null, max_null = get_residual_and_penalty(
        stat_null,
        useOsc=False,
    )

    fitter = OscillationFitter(
        sample_histogram,
        exclude=exclude_samples,
        lam=lam,
        marginalize_flux=True,
        mask_spec=mask_spec,
    )

    chi2_bf_fit, res = fitter.DoFit()

    sample_histogram.OscillateHistogram(
        res["m"],
        res["ue4"],
        res["umu4"],
        res["utau4"],
    )

    stat_bf = Statistics(
        sample_histogram,
        exclude=exclude_samples,
        lam=lam,
        mask_spec=mask_spec,
    )

    chi2_bf, resid_bf, pen_bf, norm_bf, max_bf = get_residual_and_penalty(
        stat_bf,
        useOsc=True,
    )

    row = {
        "lambda": lam,
        "log_lambda": np.log10(lam),

        "chi2_null": chi2_null,
        "resid_null": resid_null,
        "penalty_null": pen_null,
        "norm_a_null": norm_null,
        "max_abs_a_null": max_null,

        "chi2_bf": chi2_bf,
        "resid_bf": resid_bf,
        "penalty_bf": pen_bf,
        "norm_a_bf": norm_bf,
        "max_abs_a_bf": max_bf,

        "delta_chi2": chi2_null - chi2_bf,

        "dm2": res["m"],
        "ue4": res["ue4"],
        "umu4": res["umu4"],
        "utau4": res["utau4"],
    }

    return row


def main():
    global args

    # Match fitAsimovs_quinn.py seeding convention closely enough.
    np.random.seed(args.seed)

    ccnueroot = os.environ.get("CCNUEROOT")
    if ccnueroot is None:
        raise RuntimeError("CCNUEROOT is not set")

    plot_tag = args.hist_config_tag
    if plot_tag is None:
        plot_tag = getattr(AnalysisConfig, "hist_config_tag", "default")

    if plot_tag in [None, "", "none"]:
        plot_tag = "default"

    exclude_samples = args.exclude
    if exclude_samples is None:
        exclude_samples = getattr(AnalysisConfig, "exclude", "")

    if exclude_samples is None:
        exclude_samples = ""

    if isinstance(exclude_samples, str):
        if exclude_samples.strip().lower() in ["none", ""]:
            exclude_samples = ""

    mask_spec = parse_mask(args.mask)

    filename = "rootfiles/NuE_stitched_hists_{}.root".format(plot_tag)
    file_path = "{}/oscillations/{}".format(ccnueroot, filename)

    hist_config = "HIST_CONFIG_{}.json".format(plot_tag)
    if not os.path.exists(hist_config):
        raise RuntimeError("Missing requested hist config file: {}".format(hist_config))

    shutil.copyfile(hist_config, "HIST_CONFIG.json")

    lambda_values = [float(x) for x in args.lambdas.split(",")]

    # For FC-consistent toy throwing:
    # profileFlux=True means throw flux separately, then sans-flux syst, then Poisson.
    # profileFlux=False means throw full syst-only covariance, then Poisson.
    profileFlux = True

    if args.out is None:
        mask_label = args.mask if args.mask not in [None, ""] else "none"
        exclude_label = exclude_samples if exclude_samples not in [None, ""] else "none"
        out_csv = "lambda_scan_asimov_{}_exclude_{}_mask_{}_fcThrow_ntoys{}.csv".format(
            plot_tag,
            exclude_label.replace(",", "_"),
            mask_label,
            args.n_toys,
        )
    else:
        out_csv = args.out

    out_dir = os.path.dirname(out_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    print("\n===== FC-style Asimov/toy lambda scan setup =====")
    print("plot_tag       =", plot_tag)
    print("file           =", file_path)
    print("hist_config    =", hist_config)
    print("exclude        =", exclude_samples)
    print("mask           =", args.mask)
    print("mask_spec      =", mask_spec)
    print("profileFlux    =", profileFlux)
    print("n_toys         =", args.n_toys)
    print("seed           =", args.seed)
    print("lambdas        =", lambda_values)
    print("output csv     =", out_csv)

    # Load once to generate all toys.
    toy_source_hist = StitchedHistogram("toy_source")
    toy_source_hist.Load(file_path)

    print("\n===== Loaded stitched covariance check =====")
    print("full covariance shape      =", toy_source_hist.GetCovarianceMatrix(False).shape)
    print("sans-flux covariance shape =", toy_source_hist.GetCovarianceMatrix(True).shape)

    # Generate toys using the same method as fitAsimovs_quinn.py:
    # MC CV -> Gaussian systematic throw -> Poisson throw.
    # If profileFlux=True:
    #   flux covariance throw first,
    #   then sans-flux systematic covariance throw around flux-shifted mean,
    #   then Poisson.
    if args.toy_mode == "asimov":
        base = np.array(toy_source_hist.GetMCHistogram())[1:-1].astype(float)
        experiments = np.asarray([base.copy() for _ in range(args.n_toys)])
    elif args.toy_mode == "poisson":
        base = np.array(toy_source_hist.GetMCHistogram())[1:-1].astype(float)
        experiments = np.random.poisson(base, size=(args.n_toys, len(base)))
    elif args.toy_mode in ["gauss", "gauss_poisson"]:
        sys_throws = ThrowSystematics(
            toy_source_hist,
            throwFlux=profileFlux,
            n_samples=args.n_toys,
            doDiagnostics=False,
        )

        if args.toy_mode == "gauss":
            experiments = sys_throws
        else:
            experiments = ThrowPoissons(
                sys_throws,
                toy_source_hist,
                throwFlux=profileFlux,
            )
    else:
        raise ValueError("Unknown toy_mode: {}".format(args.toy_mode))

    experiments = np.asarray(experiments)
    if experiments.ndim == 1:
        experiments = experiments.reshape(1, -1)

    print("\n===== Generated experiments =====")
    print("experiments shape =", experiments.shape)
    print("sum min/max       =", np.min(np.sum(experiments, axis=1)), np.max(np.sum(experiments, axis=1)))
    print("bin min/max       =", np.min(experiments), np.max(experiments))

    rows = []

    for toy_id, pseudo_vector in enumerate(experiments):
        print("\n\n################################################")
        print("Scanning toy =", toy_id)
        print("pseudo sum =", np.sum(pseudo_vector))
        print("pseudo min =", np.min(pseudo_vector))
        print("pseudo max =", np.max(pseudo_vector))
        print("################################################")

        for lam in lambda_values:
            print("\n\n==============================================")
            print("Toy =", toy_id, "lambda =", lam)
            print("==============================================")

            # Reload fresh histogram for each lambda to avoid stale oscillated state.
            sample_histogram = StitchedHistogram("sample")
            sample_histogram.Load(file_path)

            # Apply exactly the same pseudo-data vector for every lambda in this toy.
            overwrite_data_hist_from_stitched_vector(
                sample_histogram,
                pseudo_vector,
            )
            check_data = np.array(sample_histogram.GetDataHistogram())[1:-1]
            print("check data sum after SetDataHistogram =", np.sum(check_data))
            print("check max |data - toy| =", np.max(np.abs(check_data - pseudo_vector)))


            row = run_one_lambda(
                sample_histogram=sample_histogram,
                lam=lam,
                exclude_samples=exclude_samples,
                mask_spec=mask_spec,
            )

            row["toy_id"] = toy_id
            row["seed"] = args.seed
            row["toy_mode"] = args.toy_mode
            row["throw_style"] = "fitAsimovs_quinn"
            row["hist_config_tag"] = plot_tag
            row["exclude"] = exclude_samples if exclude_samples != "" else "none"
            row["mask"] = args.mask

            rows.append(row)

            print("\n===== toy lambda result =====")
            for k in [
                "toy_id",
                "lambda",
                "log_lambda",
                "chi2_null",
                "chi2_bf",
                "delta_chi2",
                "resid_null",
                "resid_bf",
                "penalty_null",
                "penalty_bf",
                "norm_a_null",
                "norm_a_bf",
                "max_abs_a_null",
                "max_abs_a_bf",
                "dm2",
                "ue4",
                "umu4",
                "utau4",
            ]:
                print("{:20s} {}".format(k, row[k]))

            write_rows_csv(rows, out_csv)
            print("Wrote partial results to", out_csv)

    write_rows_csv(rows, out_csv)

    print("\n===== final Asimov/toy lambda scan table =====")
    for row in rows:
        print(row)

    print("\nSaved:", out_csv)


if __name__ == "__main__":
    main()
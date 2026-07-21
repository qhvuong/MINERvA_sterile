import os
import time
import logging, sys
import argparse
import math
try:
    import psutil
except ImportError:
    psutil = None

from array import array

import numpy as np
from scipy import optimize

ccnueroot = os.environ.get('CCNUEROOT')

import ROOT
import ctypes
import PlotUtils
#insert path for modules of this package.
from tools.PlotLibrary import HistHolder
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS 

from tools.StitchedHistogram import *
from tools.Helper import *

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

MNVPLOTTER = PlotUtils.MnvPlotter()
MNVPLOTTER.error_summary_group_map.clear()
for k,v in CONSOLIDATED_ERROR_GROUPS.items():
    vec = ROOT.vector("std::string")()
    for vs in v :
        vec.push_back(vs)
    MNVPLOTTER.error_summary_group_map[k]= vec

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

_PRINTED_FLUXFITTER_DIAGNOSTIC = False
_PRINTED_FLUX_SVD = set()
_PRINTED_FLUX_WEIGHTED_SVD = set()
_PRINTED_RELATIVE_A_DIAGNOSTIC = False

def GetMaskedBinIndices(mask_spec, hist_config="HIST_CONFIG.json", verbose=False):
    """
    Convert a mask specification into zero-based global bin indices.

    mask_spec example:
        {
            "fhc_ratio": [1, 9, 10],
            "rhc_ratio": [1, 9, 10],
        }

    local bin numbers are 1-based within each sample.
    returned global bin indices are 0-based for numpy arrays.
    """
    if mask_spec is None:
        return []

    import json

    with open(hist_config, "r") as f:
        cfg = json.load(f)

    masked = []

    for sample, local_bins in mask_spec.items():
        if sample not in cfg:
            print("WARNING: sample {} not found in HIST_CONFIG.json".format(sample))
            continue

        start = cfg[sample]["start"]  # zero-based global start

        for local_bin in local_bins:
            global_idx0 = start + local_bin - 1
            masked.append(global_idx0)

    masked = sorted(set(masked))
    if verbose:
        print("Masked global bin indices, zero-based:", masked)
        print("Masked ROOT/global bin numbers, one-based:", [x + 1 for x in masked])

    return masked

def GetProfileOnlyBinIndices(profile_only, hist_keys, hist_config="HIST_CONFIG.json", verbose=False):
    """
    Return zero-based global bin indices to use in the flux solve only.

    profile_only="ratio" keeps only fhc_ratio/rhc_ratio bins.
    """
    if profile_only in [None, "", "none", "None"]:
        return None

    import json

    with open(hist_config, "r") as f:
        cfg = json.load(f)

    if profile_only == "ratio":
        keep_samples = ["fhc_ratio", "rhc_ratio"]
    else:
        raise ValueError("Unknown profile_only option: {}".format(profile_only))

    keep = []

    for sample in keep_samples:
        if sample not in cfg:
            print("WARNING: profile_only sample {} not found in HIST_CONFIG.json".format(sample))
            continue

        if sample not in hist_keys:
            print("WARNING: profile_only sample {} not found in histogram keys".format(sample))
            continue

        start = int(cfg[sample]["start"])
        end = int(cfg[sample]["end"])

        keep.extend(list(range(start, end + 1)))

    keep = sorted(set(keep))

    if len(keep) == 0:
        raise RuntimeError("profile_only={} selected zero bins".format(profile_only))

    if verbose:
        print("profile_only =", profile_only)
        print("profile-only global bin indices, zero-based:", keep)
        print("profile-only ROOT/global bin numbers, one-based:", [x + 1 for x in keep])

    return keep

def GetFluxSolveBinIndices(hist_keys, exclude, profile_only=None, mask_bins=None):
    """
    Return the zero-based global bin indices used in the flux solve.
    """
    profile_only_inds = GetProfileOnlyBinIndices(
        profile_only,
        hist_keys,
        hist_config="HIST_CONFIG.json",
        verbose=False,
    )

    if profile_only_inds is None:
        solve_inds = GetSliceIndices("HIST_CONFIG.json", exclude, hist_keys)
    else:
        solve_inds = profile_only_inds

    if mask_bins is not None and len(mask_bins) > 0:
        solve_inds = [i for i in solve_inds if i not in mask_bins]

    solve_inds = sorted(set(solve_inds))

    if len(solve_inds) == 0:
        raise RuntimeError(
            "No bins selected for flux solve. exclude={}, profile_only={}, mask_bins={}".format(
                exclude,
                profile_only,
                mask_bins,
            )
        )

    return solve_inds

def GetSampleBinRanges(hist_keys, hist_config="HIST_CONFIG.json"):
    """
    Return zero-based global bin ranges for each stitched sample.

    The returned indices correspond to the internal arrays after removing
    ROOT underflow/overflow, i.e. np.array(hist)[1:-1].
    """
    import json

    with open(hist_config, "r") as f:
        cfg = json.load(f)

    ranges = {}

    for key in hist_keys:
        if key not in cfg:
            print("WARNING: sample {} not found in {}".format(key, hist_config))
            continue

        start = int(cfg[key]["start"])
        end = int(cfg[key]["end"])

        ranges[key] = list(range(start, end + 1))

    return ranges


class OscillationFitter():
    def __init__(
        self,
        histogram,
        exclude="ratio",
        lam=1,
        marginalize_flux=True,
        mask_spec=None,
        profile_only=None,
        profile_n_universes=None,
    ):
        self.hist = histogram
        self.tol = 1e-12
        self.exclude = exclude
        self.lam = lam
        self.marginalize_flux = marginalize_flux
        self.mask_spec = mask_spec
        self.profile_only = profile_only
        self.profile_n_universes = profile_n_universes

        self.statistic = Statistics(
            self.hist,
            exclude=self.exclude,
            lam=self.lam,
            mask_spec=self.mask_spec,
            profile_only=self.profile_only,
            profile_n_universes=self.profile_n_universes,
        )


    def DoFit(self):
        x0 = np.array([0.0, 0.0, 0.0, 0.0], dtype=float)

        # Current convention: dm2 = x[0] * 100
        bounds = np.array([
            [0.0, 1.0],   # dm2 / 100
            [0.0, 0.394],  # Ue4^2
            [0.0, 0.489],  # Umu4^2
            [0.0, 0.718],  # Utau4^2
        ], dtype=float)

        cons = optimize.LinearConstraint([[0, 1, 1, 1]], -np.inf, 1)

        candidates = []

        print("\n===== robust OscillationFitter.DoFit: multistart =====")

        # Near-null check
        null = optimize.minimize(
            fun=self.CalChi2,
            x0=x0,
            method="SLSQP",
            bounds=bounds,
            constraints=cons,
            tol=1e-10,
            options={
                "maxiter": 300,
                "ftol": 1e-10,
                "disp": False,
            },
        )

        print("near-null chi2 = {:.6f}, x = {}".format(float(null.fun), null.x))
        candidates.append(("near_null", null))

        # One global differential-evolution search
        print("\n===== differential evolution global search =====")
        res_de = optimize.differential_evolution(
            func=self.CalChi2,
            bounds=bounds,
            constraints=cons,
            x0=x0,
            maxiter=200,
            popsize=15,
            tol=1e-7,
            atol=1e-7,
            polish=False,
            disp=True,
            updating="immediate",
            workers=1,
        )

        print("DE success =", res_de.success)
        print("DE message =", res_de.message)
        print("DE chi2    =", res_de.fun)
        print("DE x       =", res_de.x)
        print("DE dm2     =", res_de.x[0] * 100.0)
        candidates.append(("DE", res_de))

        # Polish the DE result
        res_de_polish = optimize.minimize(
            fun=self.CalChi2,
            x0=res_de.x,
            method="SLSQP",
            bounds=bounds,
            constraints=cons,
            tol=1e-10,
            options={
                "maxiter": 300,
                "ftol": 1e-10,
                "disp": False,
            },
        )

        print(
            "DE polish -> chi2={:.6f}, dm2={:.4f}, Ue4={:.5f}, Umu4={:.5f}, Utau4={:.5f}, success={}".format(
                float(res_de_polish.fun),
                res_de_polish.x[0] * 100.0,
                res_de_polish.x[1],
                res_de_polish.x[2],
                res_de_polish.x[3],
                res_de_polish.success,
            )
        )

        candidates.append(("DE_polish", res_de_polish))

        # Multistart local fits across the full dm2 range.
        # These are generic seeds, not tuned to known minima.
        dm2_seeds = np.linspace(0.5, 99.5, 10)

        # A few generic mixing seeds.
        # Keep them inside physical bounds and away from exactly zero
        # so SLSQP has nonzero gradients in oscillation parameters.
        mixing_seeds = [
            (0.02, 0.02, 0.00),
            (0.05, 0.02, 0.00),
            (0.10, 0.02, 0.00),
            (0.15, 0.02, 0.00),
        ]

        print("\n===== multistart SLSQP over generic dm2 seeds =====")

        for dm2_seed in dm2_seeds:
            for ue_seed, umu_seed, utau_seed in mixing_seeds:
                x_seed = np.array([
                    dm2_seed / 100.0,
                    ue_seed,
                    umu_seed,
                    utau_seed,
                ], dtype=float)

                # Clip to bounds
                x_seed[0] = min(max(x_seed[0], bounds[0, 0]), bounds[0, 1])
                x_seed[1] = min(max(x_seed[1], bounds[1, 0]), bounds[1, 1])
                x_seed[2] = min(max(x_seed[2], bounds[2, 0]), bounds[2, 1])
                x_seed[3] = min(max(x_seed[3], bounds[3, 0]), bounds[3, 1])

                r = optimize.minimize(
                    fun=self.CalChi2,
                    x0=x_seed,
                    method="SLSQP",
                    bounds=bounds,
                    constraints=cons,
                    tol=1e-10,
                    options={
                        "maxiter": 300,
                        "ftol": 1e-10,
                        "disp": False,
                    },
                )

                print(
                    "seed dm2={:7.2f}, ue={:6.3f}, umu={:6.3f} -> "
                    "chi2={:12.6f}, dm2={:10.4f}, Ue4={:9.5f}, Umu4={:9.5f}, Utau4={:9.5f}, success={}".format(
                        dm2_seed,
                        ue_seed,
                        umu_seed,
                        float(r.fun),
                        r.x[0] * 100.0,
                        r.x[1],
                        r.x[2],
                        r.x[3],
                        r.success,
                    )
                )

                candidates.append(
                    (
                        "seed_dm2_{:.2f}_ue_{:.3f}_umu_{:.3f}".format(
                            dm2_seed, ue_seed, umu_seed
                        ),
                        r,
                    )
                )

        best_name, best = min(candidates, key=lambda item: float(item[1].fun))
        chi2 = float(best.fun)

        print("\n===== selected best fit candidate =====")
        print("best candidate =", best_name)
        print("best chi2      =", chi2)
        print("best x         =", best.x)
        print("best dm2       =", best.x[0] * 100.0)
        print("best Ue4       =", best.x[1])
        print("best Umu4      =", best.x[2])
        print("best Utau4     =", best.x[3])

        return (
            chi2,
            {
                "m": best.x[0] * 100.0,
                "ue4": best.x[1],
                "umu4": best.x[2],
                "utau4": best.x[3],
            },
        )

    def CalChi2(self,x):
        ms = x[0]*100
        U_e4 = x[1]
        U_mu4 = x[2]
        U_tau4 = x[3]
        self.hist.OscillateHistogram(ms,U_e4,U_mu4,U_tau4,False,False)

        chi2, penalty = self.statistic.Chi2DataMC(
            marginalize=self.marginalize_flux,
            useOsc=True
        )

        return chi2


class FluxFitter():
    def __init__(
        self,
        histogram,
        exclude="",
        lam=1,
        usePseudo=False,
        useOsc=False,
        mask_bins=None,
        profile_only=None,
        profile_n_universes=None,
    ):
        self.hist = histogram
        self.exclude = exclude
        self.lam = lam
        self.usePseudo = usePseudo
        self.useOsc = useOsc
        self.mask_bins = mask_bins if mask_bins is not None else []
        self.profile_only = profile_only
        self.profile_n_universes = profile_n_universes

        if usePseudo:
            self.dataHist = histogram.GetPseudoHistogram()
        else:
            self.dataHist = histogram.GetDataHistogram()

        if useOsc:
            self.mcHist = histogram.GetOscillatedHistogram()
        else:
            self.mcHist = histogram.GetMCHistogram()

        self.SolveFluxSolution()

    def SetHistogram(self, hist):
        self.hist = hist

        if self.usePseudo:
            self.dataHist = hist.GetPseudoHistogram()
        else:
            self.dataHist = hist.GetDataHistogram()

        if self.useOsc:
            self.mcHist = hist.GetOscillatedHistogram()
        else:
            self.mcHist = hist.GetMCHistogram()

        self.SolveFluxSolution()

    def GetProfileAMatrixAndUniverses(self):
        """
        Return the flux-universe basis used for profiling.

        Relative-A diagnostic mode:
            R[k,i] = A_abs[k,i] / MC_CV[i]
            A_used[k,i] = R[k,i] * MC_current[i]

        At null, or when USE_RELATIVE_FLUX_A is disabled,
        use the original absolute A matrix.
        """
        universes = self.hist.GetFluxUniverses()

        A_abs = np.asarray(
            self.hist.GetAMatrix(),
            dtype=float,
        )

        use_relative_A = (
            os.environ.get("USE_RELATIVE_FLUX_A", "0") == "1"
        )

        unosc_cv = None
        current_mc = None
        relative_A = None

        if use_relative_A and self.useOsc:
            unosc_cv = np.asarray(
                self.hist.GetMCHistogram(),
                dtype=float,
            )[1:-1]

            current_mc = np.asarray(
                self.mcHist,
                dtype=float,
            )[1:-1]

            if A_abs.shape[1] != len(unosc_cv):
                raise RuntimeError(
                    "A columns {} do not match unoscillated CV bins {}".format(
                        A_abs.shape[1],
                        len(unosc_cv),
                    )
                )

            if len(current_mc) != len(unosc_cv):
                raise RuntimeError(
                    "Oscillated MC bins {} do not match unoscillated CV bins {}".format(
                        len(current_mc),
                        len(unosc_cv),
                    )
                )

            relative_A = np.divide(
                A_abs,
                unosc_cv[None, :],
                out=np.zeros_like(A_abs),
                where=np.abs(unosc_cv[None, :]) > 1e-12,
            )

            A = relative_A * current_mc[None, :]

        else:
            A = A_abs

        # Apply Nprof exactly once.
        if self.profile_n_universes is not None:
            nprof = int(self.profile_n_universes)

            if nprof <= 0:
                raise RuntimeError(
                    "profile_n_universes must be positive, got {}".format(
                        nprof
                    )
                )

            if nprof > A.shape[0]:
                raise RuntimeError(
                    "profile_n_universes={} requested, but A only has "
                    "{} universes".format(
                        nprof,
                        A.shape[0],
                    )
                )

            universes = universes[:nprof]
            A = A[:nprof, :]

        # Print only at the first genuinely oscillated point.
        global _PRINTED_RELATIVE_A_DIAGNOSTIC

        if (
            use_relative_A
            and self.useOsc
            and not _PRINTED_RELATIVE_A_DIAGNOSTIC
        ):
            osc_to_cv = np.divide(
                current_mc,
                unosc_cv,
                out=np.ones_like(current_mc),
                where=np.abs(unosc_cv) > 1e-12,
            )

            nontrivial = (
                np.max(np.abs(osc_to_cv - 1.0)) > 1e-6
            )

            if nontrivial:
                A_abs_used = A_abs

                if self.profile_n_universes is not None:
                    A_abs_used = A_abs_used[
                        :int(self.profile_n_universes),
                        :,
                    ]

                print("")
                print("===== relative flux-A nontrivial diagnostic =====")
                print("returned A shape        =", A.shape)
                print("absolute A used shape   =", A_abs_used.shape)
                print("min osc/CV ratio        =", np.min(osc_to_cv))
                print("max osc/CV ratio        =", np.max(osc_to_cv))
                print("mean osc/CV ratio       =", np.mean(osc_to_cv))
                print(
                    "max |A-A_abs_used|      =",
                    np.max(np.abs(A - A_abs_used)),
                )
                print(
                    "mean |A-A_abs_used|     =",
                    np.mean(np.abs(A - A_abs_used)),
                )
                print(
                    "fraction changed >1e-8  =",
                    np.mean(np.abs(A - A_abs_used) > 1e-8),
                )
                print(
                    "relative-A finite       =",
                    np.all(np.isfinite(relative_A)),
                )
                print(
                    "returned A finite       =",
                    np.all(np.isfinite(A)),
                )
                print("===== end relative flux-A diagnostic =====")
                print("")

                _PRINTED_RELATIVE_A_DIAGNOSTIC = True

        return universes, A

    def SolveFluxSolution(self):
        data = np.array(self.dataHist)[1:-1]
        mc = np.array(self.mcHist)[1:-1]

        # universes = self.hist.GetFluxUniverses()
        # A = self.hist.GetAMatrix()
        universes, A = self.GetProfileAMatrixAndUniverses()

        profile_only_inds = GetProfileOnlyBinIndices(
            self.profile_only,
            self.hist.keys,
            hist_config="HIST_CONFIG.json",
            verbose=False,
        )

        # if profile_only_inds is None:
        #     sliceInds = GetSliceIndices("HIST_CONFIG.json", self.exclude, self.hist.keys)
        # else:
        #     # profile_only overrides exclude for the flux solve.
        #     sliceInds = profile_only_inds

        # if len(self.mask_bins) > 0:
        #     sliceInds = [i for i in sliceInds if i not in self.mask_bins]

        sliceInds = GetFluxSolveBinIndices(
            self.hist.keys,
            self.exclude,
            profile_only=self.profile_only,
            mask_bins=self.mask_bins,
        )

        if len(sliceInds) == 0:
            raise RuntimeError(
                "No bins selected for flux solve after profile_only/exclude/mask. "
                "profile_only={}, exclude={}, mask_bins={}".format(
                    self.profile_only,
                    self.exclude,
                    self.mask_bins,
                )
            )

        data = slicer(data, sliceInds)
        mc   = slicer(mc, sliceInds)
        A    = slicer(A, sliceInds, axis=1)

        # -------------------------------------------------
        # SVD diagnostic for the actual profiling matrix
        # -------------------------------------------------
        global _PRINTED_FLUX_SVD

        svd_key = (
            self.profile_n_universes,
            str(self.profile_only),
            str(self.exclude),
            tuple(self.mask_bins),
            tuple(sliceInds),
        )

        if (
            os.environ.get("DEBUG_FLUX_SVD", "0") == "1"
            and svd_key not in _PRINTED_FLUX_SVD
        ):
            _PRINTED_FLUX_SVD.add(svd_key)

            singular_values = np.linalg.svd(
                A,
                full_matrices=False,
                compute_uv=False,
            )

            tolerance = (
                np.finfo(float).eps
                * max(A.shape)
                * singular_values[0]
                if len(singular_values) > 0
                else 0.0
            )

            numerical_rank = int(np.sum(singular_values > tolerance))

            variance = singular_values**2
            total_variance = np.sum(variance)

            if total_variance > 0:
                variance_fraction = variance / total_variance
                cumulative_fraction = np.cumsum(variance_fraction)
            else:
                variance_fraction = np.zeros_like(variance)
                cumulative_fraction = np.zeros_like(variance)

            print("")
            print("===== flux profiling A-matrix SVD =====")
            print("profile_n_universes =", self.profile_n_universes)
            print("profile_only        =", self.profile_only)
            print("exclude             =", self.exclude)
            print("mask_bins           =", self.mask_bins)
            print("A shape             =", A.shape)
            print("numerical rank      =", numerical_rank)
            print("SVD tolerance       =", tolerance)

            if len(singular_values) > 0:
                print("largest singular value  =", singular_values[0])
                print("smallest singular value =", singular_values[-1])

                if singular_values[-1] > 0:
                    print(
                        "condition number A       =",
                        singular_values[0] / singular_values[-1],
                    )
                else:
                    print("condition number A       = inf")

            for target in [0.90, 0.95, 0.99, 0.999]:
                if total_variance > 0:
                    n_required = (
                        int(np.searchsorted(cumulative_fraction, target)) + 1
                    )
                else:
                    n_required = 0

                print(
                    "{:5.1f}% cumulative S^2 -> {} modes".format(
                        100.0 * target,
                        n_required,
                    )
                )

            print("")
            print(" index      singular value       frac S^2       cumulative")
            for i, value in enumerate(singular_values):
                print(
                    "{:5d}  {:18.10e}  {:12.6e}  {:12.6e}".format(
                        i,
                        value,
                        variance_fraction[i],
                        cumulative_fraction[i],
                    )
                )

            print("===== end flux profiling A-matrix SVD =====")
            print("")

        global _PRINTED_FLUXFITTER_DIAGNOSTIC
        if not _PRINTED_FLUXFITTER_DIAGNOSTIC:
            print("\n===== FluxFitter universe-basis diagnostic =====")
            print("profile_n_universes =", self.profile_n_universes)
            print("A shape used in flux solve =", A.shape)
            print("number of flux nuisance parameters =", len(universes))
            print("number of bins used in flux solve =", len(sliceInds))
            _PRINTED_FLUXFITTER_DIAGNOSTIC = True

        cov_sans = self.hist.GetCovarianceMatrix(sansFlux=True)
        cov_sliced = slicer(cov_sans, sliceInds)
        V = np.linalg.pinv(cov_sliced)

        # -------------------------------------------------
        # Covariance-weighted SVD diagnostic
        # Uses the same metric as A @ V @ A.T in the fit.
        # -------------------------------------------------
        global _PRINTED_FLUX_WEIGHTED_SVD

        weighted_svd_key = (
            self.profile_n_universes,
            str(self.profile_only),
            str(self.exclude),
            tuple(self.mask_bins),
            tuple(sliceInds),
        )

        if (
            os.environ.get("DEBUG_FLUX_WEIGHTED_SVD", "0") == "1"
            and weighted_svd_key not in _PRINTED_FLUX_WEIGHTED_SVD
        ):
            _PRINTED_FLUX_WEIGHTED_SVD.add(weighted_svd_key)

            # V should be symmetric, but explicitly symmetrize against
            # small numerical asymmetries from the pseudoinverse.
            V_sym = 0.5 * (V + V.T)

            eigvals, eigvecs = np.linalg.eigh(V_sym)

            # Protect against tiny negative eigenvalues from numerical precision.
            eigvals = np.clip(eigvals, 0.0, None)

            sqrtV = (
                eigvecs
                @ np.diag(np.sqrt(eigvals))
                @ eigvecs.T
            )

            A_weighted = A @ sqrtV

            S_weighted = np.linalg.svd(
                A_weighted,
                full_matrices=False,
                compute_uv=False,
            )

            weighted_variance = S_weighted**2
            total_weighted_variance = np.sum(weighted_variance)

            if total_weighted_variance > 0:
                weighted_fraction = (
                    weighted_variance / total_weighted_variance
                )
                weighted_cumulative = np.cumsum(weighted_fraction)
            else:
                weighted_fraction = np.zeros_like(S_weighted)
                weighted_cumulative = np.zeros_like(S_weighted)

            tolerance_weighted = (
                np.finfo(float).eps
                * max(A_weighted.shape)
                * S_weighted[0]
                if len(S_weighted) > 0
                else 0.0
            )

            weighted_rank = int(
                np.sum(S_weighted > tolerance_weighted)
            )

            print("")
            print("===== covariance-weighted flux A-matrix SVD =====")
            print("profile_n_universes =", self.profile_n_universes)
            print("profile_only        =", self.profile_only)
            print("exclude             =", self.exclude)
            print("mask_bins           =", self.mask_bins)
            print("A shape             =", A.shape)
            print("A_weighted shape    =", A_weighted.shape)
            print("weighted rank       =", weighted_rank)
            print("weighted tolerance  =", tolerance_weighted)

            if len(S_weighted) > 0:
                print(
                    "largest weighted singular value  =",
                    S_weighted[0],
                )
                print(
                    "smallest weighted singular value =",
                    S_weighted[-1],
                )

                if S_weighted[-1] > 0:
                    print(
                        "condition number A_weighted     =",
                        S_weighted[0] / S_weighted[-1],
                    )
                else:
                    print("condition number A_weighted     = inf")

            for target in [0.90, 0.95, 0.99, 0.999]:
                if total_weighted_variance > 0:
                    n_required = (
                        int(
                            np.searchsorted(
                                weighted_cumulative,
                                target,
                            )
                        )
                        + 1
                    )
                else:
                    n_required = 0

                print(
                    "{:5.1f}% cumulative weighted S^2 -> {} modes".format(
                        100.0 * target,
                        n_required,
                    )
                )

            print("")
            print(
                " index      weighted singular value   "
                "frac S^2       cumulative"
            )

            for i, value in enumerate(S_weighted):
                print(
                    "{:5d}  {:22.10e}  {:12.6e}  {:12.6e}".format(
                        i,
                        value,
                        weighted_fraction[i],
                        weighted_cumulative[i],
                    )
                )

            print(
                "===== end covariance-weighted flux A-matrix SVD ====="
            )
            print("")

        C = data - mc
        I = np.identity(len(universes))

        L = 2 * A @ V @ C
        Q = A @ V @ A.T + I * self.lam

        self.fluxSolution = np.linalg.pinv(Q) @ L / 2

        if os.environ.get("DEBUG_FLUX_ANALYTIC_SOLUTION", "0") == "1":
            a = self.fluxSolution
            rhs = A @ V @ C

            def obj(a_test):
                r = C - a_test @ A
                return float(r.T @ V @ r + self.lam * (a_test @ a_test))

            def grad(a_test):
                # derivative of (C-aA)^T V (C-aA) + lambda a^T a
                return 2.0 * ((A @ V @ A.T + self.lam * np.identity(len(a_test))) @ a_test - rhs)

            chi2_zero = obj(np.zeros_like(a))
            chi2_sol = obj(a)
            g = grad(a)

            print("")
            print("===== flux analytic solution check =====")
            print("useOsc =", self.useOsc)
            print("usePseudo =", self.usePseudo)
            print("exclude =", self.exclude)
            print("profile_only =", self.profile_only)
            print("profile_n_universes =", self.profile_n_universes)
            print("lambda =", self.lam)
            print("A shape =", A.shape)
            print("V shape =", V.shape)
            print("Q shape =", Q.shape)
            print("")
            print("objective at a=0        =", chi2_zero)
            print("objective at analytic a =", chi2_sol)
            print("improvement             =", chi2_zero - chi2_sol)
            print("penalty lambda aTa      =", float(self.lam * (a @ a)))
            print("|a|                     =", float(np.linalg.norm(a)))
            print("max |a_i|               =", float(np.max(np.abs(a))))
            print("")
            print("|gradient|              =", float(np.linalg.norm(g)))
            print("max |gradient_i|        =", float(np.max(np.abs(g))))
            print("|Q a - rhs|             =", float(np.linalg.norm(Q @ a - rhs)))
            print("max |Q a - rhs|         =", float(np.max(np.abs(Q @ a - rhs))))
            print("cond(Q)                 =", float(np.linalg.cond(Q)))
            print("===== end flux analytic solution check =====")
            print("")


    def SetFluxSolution(self,solution):
        self.fluxSolution = solution

    def GetFluxSolution(self,):
        return(self.fluxSolution)

    def MarginalizeFlux(self):
        universes, A = self.GetProfileAMatrixAndUniverses()
        solution = self.fluxSolution

        penalty = solution @ solution * self.lam
        new_cv = np.array(self.mcHist)[1:-1] + solution @ A

        return(new_cv, penalty)

    def ReweightToFluxSolution(self,histogram):
        mc = np.array(histogram)[1:-1]
        band = histogram.GetVertErrorBand("Flux")
        nhists = self.mcHist.GetVertErrorBand("Flux").GetNHists()
        universes = np.array([np.array(band.GetHist(l))[1:-1] for l in range(nhists)])
        cv_table = np.array([mc for l in range(len(universes))])
        A = universes - cv_table
        if self.profile_n_universes is not None:
            nprof = int(self.profile_n_universes)

            if nprof <= 0:
                raise RuntimeError("profile_n_universes must be positive, got {}".format(nprof))

            if nprof > A.shape[0]:
                raise RuntimeError(
                    "profile_n_universes={} requested, but histogram Flux band only has {} universes".format(
                        nprof,
                        A.shape[0],
                    )
                )

            A = A[:nprof, :]

        weights = histogram.GetCVHistoWithStatError()
        new_cv = mc + self.fluxSolution @ A

        for j in range(1,weights.GetNbinsX()+1):
            weight = weights.GetBinContent(j) / new_cv[j-1] if new_cv[j-1] != 0 else weights.GetBinContent(j)
            weights.SetBinContent(j,weight)
            weights.SetBinError(j,0)

        histogram.DivideSingle(histogram,weights)

    def ReweightToFluxSolutionStoredA(self, histogram):
        mc = np.array(histogram)[1:-1]
        universes, A = self.GetProfileAMatrixAndUniverses()
        new_cv = mc + self.fluxSolution @ A

        weights = histogram.GetCVHistoWithStatError()
        for j in range(1, weights.GetNbinsX() + 1):
            weight = weights.GetBinContent(j) / new_cv[j-1] if new_cv[j-1] != 0 else weights.GetBinContent(j)
            weights.SetBinContent(j, weight)
            weights.SetBinError(j, 0)

        histogram.DivideSingle(histogram, weights)


class Statistics():
    def __init__(
        self,
        histogram,
        exclude="ratio",
        lam=1,
        mask_spec=None,
        profile_only=None,
        profile_n_universes=None,
    ):
        self.hist = histogram
        self.exclude = exclude
        self.lam = lam
        self.mask_spec = mask_spec
        self.profile_only = profile_only
        self.profile_n_universes = profile_n_universes
        self.mask_bins = GetMaskedBinIndices(mask_spec, verbose=False) if mask_spec is not None else []
        self.nulFluxFitter = None
        self.oscFluxFitter = None
        self._printed_chi2_by_sample = False

    def ApplyMaskToVectorAndMatrix(self, data, mc, invCov):
        if len(self.mask_bins) == 0:
            return data, mc, invCov

        n = len(data)
        keep = np.array([i for i in range(n) if i not in self.mask_bins], dtype=int)

        data_m = data[keep]
        mc_m = mc[keep]
        invCov_m = invCov[np.ix_(keep, keep)]

        return data_m, mc_m, invCov_m

    def ApplyMaskToVectorAndCovariance(self, data, mc, cov):
        if len(self.mask_bins) == 0:
            invCov = np.linalg.pinv(cov)
            return data, mc, invCov

        n = len(data)
        keep = np.array([i for i in range(n) if i not in self.mask_bins], dtype=int)

        data_m = data[keep]
        mc_m = mc[keep]
        cov_m = cov[np.ix_(keep, keep)]
        invCov_m = np.linalg.pinv(cov_m)

        return data_m, mc_m, invCov_m

    def GetFluxFitter(self,useOsc=False):
        if useOsc:
            return(self.oscFluxFitter)
        else:
            return(self.nulFluxFitter)

    def SetHistogram(self,histogram):
        self.hist = histogram


    def PrintChi2BinDiagnostics(self, label, data, mc, invCov, diff, max_bins=20):
        q = diff * (invCov @ diff)

        order = np.argsort(np.abs(q))[::-1]

        print("\n===== chi2 bin diagnostics: {} =====".format(label))
        print("residual chi2 from sum q_i =", np.sum(q))

        import json
        sample_ranges = []
        try:
            with open("HIST_CONFIG.json", "r") as f:
                cfg = json.load(f)

            for sample, info in cfg.items():
                sample_ranges.append((sample, info["start"], info["end"]))
        except Exception:
            sample_ranges = []

        def sample_name(idx0):
            for sample, start, end in sample_ranges:
                if start <= idx0 <= end:
                    return sample
            return "unknown"

        print("{:>5s} {:>25s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s}".format(
            "bin", "sample", "data", "mc", "diff", "q_i", "frac"
        ))

        total = np.sum(q)
        for idx0 in order[:max_bins]:
            frac = q[idx0] / total if total != 0 else 0.0
            print("{:5d} {:>25s} {:12.6g} {:12.6g} {:12.6g} {:12.6g} {:12.6f}".format(
                idx0 + 1,
                sample_name(idx0),
                data[idx0],
                mc[idx0],
                diff[idx0],
                q[idx0],
                frac
            ))

    def DebugChi2BySample(self, data_vec, mc_vec, cov, label=""):
        """
        Diagnostic: compute chi2 using each sample's sub-block covariance.

        This is not expected to add exactly to the full chi2 because the full
        covariance has off-diagonal correlations between samples.
        """
        ranges = GetSampleBinRanges(self.hist.keys)

        print("")
        print("===== Chi2-by-sample diagnostic {} =====".format(label))
        print("NOTE: sample chi2 values use sub-block covariances.")
        print("      They are diagnostic and do not necessarily sum to full chi2.")
        print("")

        total_diag = 0.0

        for key, inds in ranges.items():
            d = slicer(data_vec, inds)
            m = slicer(mc_vec, inds)
            c = slicer(cov, inds)

            try:
                v = np.linalg.pinv(c)
                r = d - m
                chi2 = float(r.T @ v @ r)
            except Exception as e:
                print("{:<24s} failed: {}".format(key, e))
                continue

            total_diag += chi2

            print(
                "{:<24s} nbins={:<3d} chi2_subblock={:12.6f}".format(
                    key,
                    len(inds),
                    chi2,
                )
            )

        print("")
        print("sum of sample sub-block chi2 = {:.6f}".format(total_diag))
        print("===== end chi2-by-sample diagnostic =====")
        print("")

    def DebugHybridFluxCovarianceChi2(self, data, mc, cov_full, cov_sans, penalty):
        """
        Diagnostic for excluded samples.

        Current profiled chi2 uses cov_sans for all bins.
        This diagnostic adds flux covariance back for bins that were NOT used
        in the flux solve, then recomputes chi2.

        Two variants are printed:
          block-only: add flux cov only inside excluded-excluded block
          rows/cols: add flux cov for any element touching an excluded bin
        """
        n = len(data)

        solve_inds = GetFluxSolveBinIndices(
            self.hist.keys,
            self.exclude,
            profile_only=self.profile_only,
            mask_bins=self.mask_bins,
        )

        solve_set = set(solve_inds)
        excluded_inds = [i for i in range(n) if i not in solve_set]

        cov_flux = cov_full - cov_sans

        cov_hybrid_block = np.array(cov_sans, copy=True)
        cov_hybrid_rowscols = np.array(cov_sans, copy=True)

        if len(excluded_inds) > 0:
            ix = np.ix_(excluded_inds, excluded_inds)
            cov_hybrid_block[ix] += cov_flux[ix]

        excluded_set = set(excluded_inds)
        for i in range(n):
            for j in range(n):
                if i in excluded_set or j in excluded_set:
                    cov_hybrid_rowscols[i, j] += cov_flux[i, j]

        def calc(label, cov):
            data_m, mc_m, inv_m = self.ApplyMaskToVectorAndCovariance(data, mc, cov)
            diff_m = data_m - mc_m
            resid = float(diff_m.T @ inv_m @ diff_m)
            total = resid + float(penalty)

            print("{:<35s} residual={:12.6f}  penalty={:12.6f}  total={:12.6f}".format(
                label,
                resid,
                float(penalty),
                total,
            ))

        print("")
        print("===== Hybrid flux-covariance chi2 diagnostic =====")
        print("exclude =", self.exclude)
        print("profile_only =", self.profile_only)
        print("profile_n_universes =", self.profile_n_universes)
        print("lambda =", self.lam)
        print("number of total bins =", n)
        print("number of flux-solve bins =", len(solve_inds))
        print("number of excluded/non-profiled bins =", len(excluded_inds))
        print("excluded/non-profiled global bins, one-based =", [i + 1 for i in excluded_inds])
        print("")

        calc("current: sansFlux all bins", cov_sans)
        calc("hybrid: add excluded block", cov_hybrid_block)
        calc("hybrid: add excluded rows/cols", cov_hybrid_rowscols)
        calc("profiled MC with full covariance", cov_full)

        print("===== end hybrid diagnostic =====")
        print("")

    def BuildHybridExcludedFluxCovariance(self, cov_full, cov_sans, mode="block"):
        """
        Build a hybrid covariance for test fits.

        Current nominal profiled chi2 uses cov_sans for all bins.

        Hybrid mode adds flux covariance back for bins that were NOT used
        in the flux solve.

        mode="block":
            add flux covariance only inside excluded-excluded block

        mode="rowscols":
            add flux covariance for any element touching an excluded bin
        """
        n = cov_sans.shape[0]

        solve_inds = GetFluxSolveBinIndices(
            self.hist.keys,
            self.exclude,
            profile_only=self.profile_only,
            mask_bins=self.mask_bins,
        )

        solve_set = set(solve_inds)
        excluded_inds = [i for i in range(n) if i not in solve_set]
        excluded_set = set(excluded_inds)

        cov_flux = cov_full - cov_sans
        cov_hybrid = np.array(cov_sans, copy=True)

        if mode == "block":
            if len(excluded_inds) > 0:
                ix = np.ix_(excluded_inds, excluded_inds)
                cov_hybrid[ix] += cov_flux[ix]

        elif mode == "rowscols":
            for i in range(n):
                for j in range(n):
                    if i in excluded_set or j in excluded_set:
                        cov_hybrid[i, j] += cov_flux[i, j]

        else:
            raise RuntimeError("Unknown hybrid covariance mode: {}".format(mode))

        return cov_hybrid

    def Chi2DataMC(self, marginalize=True, usePseudo=False, useOsc=False, printBinDiagnostics=False): 
        if useOsc:
            mcHist = self.hist.GetOscillatedHistogram()
        else:
            mcHist = self.hist.GetMCHistogram()

        if usePseudo:
            dataHist = self.hist.GetPseudoHistogram()
        else:
            dataHist = self.hist.GetDataHistogram()

        if dataHist.GetNbinsX() != mcHist.GetNbinsX():
            logging.error("breaking error in Chi2DataMC")
            logging.error(
                "The number of bins from Data ({}) and MC ({}) histograms differ. Returning -1.".format(
                    dataHist.GetNbinsX(), mcHist.GetNbinsX()
                )
            )
            return(-1)

        mc = np.array(mcHist)[1:-1]
        data = np.array(dataHist)[1:-1]

        cov_full = self.hist.GetCovarianceMatrix(sansFlux=False)
        cov_sans = self.hist.GetCovarianceMatrix(sansFlux=True)

        cov = cov_full
        penalty = 0

        if marginalize:
            fluxFitter = FluxFitter(
                self.hist,
                self.exclude,
                self.lam,
                usePseudo,
                useOsc,
                mask_bins=self.mask_bins,
                profile_only=self.profile_only,
                profile_n_universes=self.profile_n_universes,
            )
            mc, penalty = fluxFitter.MarginalizeFlux()
            # cov = cov_sans
            if os.environ.get("USE_HYBRID_EXCLUDED_FLUX_COV", "0") == "1":
                hybrid_mode = os.environ.get("HYBRID_EXCLUDED_FLUX_COV_MODE", "block")

                cov = self.BuildHybridExcludedFluxCovariance(
                    cov_full,
                    cov_sans,
                    mode=hybrid_mode,
                )

                if os.environ.get("PRINT_HYBRID_TEST_MODE", "1") == "1" and not useOsc:
                    print("")
                    print("===== USING HYBRID EXCLUDED FLUX COVARIANCE FOR FIT =====")
                    print("mode =", hybrid_mode)
                    print("exclude =", self.exclude)
                    print("profile_only =", self.profile_only)
                    print("This is a diagnostic/test fit, not nominal.")
                    print("===== END HYBRID TEST MODE MESSAGE =====")
                    print("")
            else:
                cov = cov_sans

            if useOsc:
                self.oscFluxFitter = fluxFitter
            else:
                self.nulFluxFitter = fluxFitter

        data_chi2, mc_chi2, invCov_chi2 = self.ApplyMaskToVectorAndCovariance(
            data,
            mc,
            cov
        )

        diff = data_chi2 - mc_chi2
        chi2 = diff.T @ invCov_chi2 @ diff + penalty

        if (
            os.environ.get("DEBUG_CHI2_HYBRID", "0") == "1"
            and marginalize
            and not useOsc
        ):
            self.DebugHybridFluxCovarianceChi2(
                data,
                mc,
                cov_full,
                cov_sans,
                penalty,
            )

        if os.environ.get("DEBUG_CHI2_BY_SAMPLE", "0") == "1" and not self._printed_chi2_by_sample:
            self._printed_chi2_by_sample = True
            if useOsc:
                label = "oscillated"
            else:
                label = "null"

            if marginalize:
                label += "_profiled"
            else:
                label += "_raw"

            resid_chi2 = float(diff.T @ invCov_chi2 @ diff)

            print("")
            print("===== Global chi2 diagnostic: {} =====".format(label))
            print("exclude =", self.exclude)
            print("profile_only =", self.profile_only)
            print("profile_n_universes =", self.profile_n_universes)
            print("lambda =", self.lam)
            print("residual chi2 =", resid_chi2)
            print("flux penalty  =", float(penalty))
            print("total chi2    =", float(chi2))
            print("===== end global chi2 diagnostic =====")
            print("")

            self.DebugChi2BySample(
                data,
                mc,
                cov,
                label="{} exclude={} lambda={} profile_only={}".format(
                    label,
                    self.exclude,
                    self.lam,
                    self.profile_only,
                ),
            )

        if printBinDiagnostics:
            if useOsc:
                label = "oscillated"
            else:
                label = "null"

            if marginalize:
                label += "_profiled"
            else:
                label += "_raw"

            self.PrintChi2BinDiagnostics(
                label,
                data_chi2,
                mc_chi2,
                invCov_chi2,
                diff,
                max_bins=20
            )




        if abs(chi2) > 1e30:
            logging.error("chi2 has invalid value: {}".format(chi2))
            print("chi2 has invalid value: {}".format(chi2))
            return(-1)

        return(chi2, penalty)
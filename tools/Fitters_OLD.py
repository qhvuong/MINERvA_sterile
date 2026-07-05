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

class OscillationFitter():
    def __init__(self, histogram, exclude="ratio", lam=1, marginalize_flux=True):
        self.hist = histogram
        self.tol = 1e-12
        self.exclude = exclude
        self.lam = lam
        self.marginalize_flux = marginalize_flux
        self.statistic = Statistics(histogram, exclude, lam)


    def DoFit(self):
        x0 = np.array([0.0, 0.0, 0.0, 0.0], dtype=float)

        # Current convention: dm2 = x[0] * 100
        bounds = np.array([
            [0.0, 1.0],   # dm2 / 100
            [0.0, 0.15],  # Ue4^2
            [0.0, 0.41],  # Umu4^2
            [0.0, 0.66],  # Utau4^2
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
    # def DoFit(self):
    #     x0 = [0.0,0.0,0.0,0.0]
    #     bounds = np.array([[0.0,1.0],[0.0,0.15],[0.0,0.41],[0,0.66]], dtype = float)
    #     cons = optimize.LinearConstraint([[0,1,1,1]],-np.inf,1)

    #     null = optimize.minimize(fun=self.CalChi2,x0=x0,tol=self.tol,options={"maxiter":40},method="SLSQP",bounds=bounds,constraints=cons)
    #     null_chi2 = float(null.fun)
    #     print("fit near null: {}".format(null.fun))

    #     res = optimize.differential_evolution(func=self.CalChi2,tol=self.tol,bounds=bounds,polish=False,x0=x0,maxiter=100,disp=True,constraints=cons)
    #     new_x0 = res.x
    #     print("best fit: {}".format(res.fun))

    #     res = optimize.minimize(fun=self.CalChi2,x0=new_x0,tol=self.tol,options={"maxiter":20},method="SLSQP",bounds=bounds,constraints=cons)
    #     chi2 = float(res.fun)
    #     print("polish fit: {}".format(res.fun))

    #     if null_chi2 < chi2:
    #         chi2 = null_chi2
    #         res = null

    #     return(chi2,{"m":res.x[0]*100,"ue4":res.x[1],"umu4":res.x[2],"utau4":res.x[3]})

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

# class OscillationFitter():
#     def __init__(self, histogram, exclude="ratio", lam=1, marginalize_flux=True):
#         self.hist = histogram
#         self.tol = 1e-12
#         self.exclude = exclude
#         self.lam = lam
#         self.marginalize_flux = marginalize_flux
#         self.statistic = Statistics(histogram, exclude, lam)

#     def _RunSLSQP(self, x_seed, bounds, cons, label="", maxiter=300, ftol=1e-10, verbose=True):
#         x_seed = np.array(x_seed, dtype=float)

#         # Clip to bounds.
#         for i in range(len(x_seed)):
#             x_seed[i] = min(max(x_seed[i], bounds[i, 0]), bounds[i, 1])

#         res = optimize.minimize(
#             fun=self.CalChi2,
#             x0=x_seed,
#             method="SLSQP",
#             bounds=bounds,
#             constraints=cons,
#             tol=ftol,
#             options={
#                 "maxiter": maxiter,
#                 "ftol": ftol,
#                 "disp": False,
#             },
#         )

#         if verbose:
#             print(
#                 "{:35s} chi2={:12.6f}, dm2={:10.4f}, "
#                 "Ue4={:9.5f}, Umu4={:9.5f}, Utau4={:9.5f}, success={}".format(
#                     label,
#                     float(res.fun),
#                     res.x[0] * 100.0,
#                     res.x[1],
#                     res.x[2],
#                     res.x[3],
#                     res.success,
#                 )
#             )

#         return res

#     def DoFit(self):
#         x0 = np.array([0.0, 0.0, 0.0, 0.0], dtype=float)

#         # Current convention: dm2 = x[0] * 100
#         bounds = np.array([
#             [0.0, 1.0],   # dm2 / 100, so 0--100 eV^2
#             [0.0, 0.15],  # Ue4^2
#             [0.0, 0.41],  # Umu4^2
#             [0.0, 0.66],  # Utau4^2
#         ], dtype=float)

#         cons = optimize.LinearConstraint([[0, 1, 1, 1]], -np.inf, 1)

#         candidates = []

#         print("\n===== robust OscillationFitter.DoFit: dm2-focused multistart =====")

#         # -------------------------------------------------
#         # 0. Near-null check
#         # -------------------------------------------------
#         null = self._RunSLSQP(
#             x_seed=x0,
#             bounds=bounds,
#             cons=cons,
#             label="near_null",
#             maxiter=300,
#             ftol=1e-10,
#         )
#         candidates.append(("near_null", null))

#         # -------------------------------------------------
#         # 1. One global DE search
#         # -------------------------------------------------
#         print("\n===== differential evolution global search =====")

#         res_de = optimize.differential_evolution(
#             func=self.CalChi2,
#             bounds=bounds,
#             constraints=cons,
#             x0=x0,
#             maxiter=200,
#             popsize=15,
#             tol=1e-7,
#             atol=1e-7,
#             polish=False,
#             disp=True,
#             updating="immediate",
#             workers=1,
#         )

#         print("DE success =", res_de.success)
#         print("DE message =", res_de.message)
#         print("DE chi2    =", res_de.fun)
#         print("DE x       =", res_de.x)
#         print("DE dm2     =", res_de.x[0] * 100.0)

#         candidates.append(("DE", res_de))

#         # -------------------------------------------------
#         # 2. Polish DE result
#         # -------------------------------------------------
#         res_de_polish = self._RunSLSQP(
#             x_seed=res_de.x,
#             bounds=bounds,
#             cons=cons,
#             label="DE_polish",
#             maxiter=300,
#             ftol=1e-10,
#         )
#         candidates.append(("DE_polish", res_de_polish))

#         # -------------------------------------------------
#         # 3. Generic dm2-focused multistart
#         # -------------------------------------------------
#         # Dense enough in the low/mid dm2 region where phase structure matters,
#         # but still covers the full 0--100 eV^2 allowed range.
#         dm2_seeds = np.array([
#             0.5, 1.0, 2.0, 3.0, 5.0,
#             8.0, 10.0, 13.0, 15.0,
#             18.0, 20.0, 25.0, 30.0,
#             40.0, 50.0, 70.0, 100.0,
#         ], dtype=float)

#         # Only a couple of generic nonzero mixing seeds.
#         # The local minima are expected to be driven mostly by dm2.
#         mixing_seeds = [
#             (0.05, 0.02, 0.00),
#             (0.15, 0.02, 0.00),
#         ]

#         print("\n===== multistart SLSQP over dm2-focused generic seeds =====")

#         for dm2_seed in dm2_seeds:
#             for ue_seed, umu_seed, utau_seed in mixing_seeds:
#                 x_seed = np.array([
#                     dm2_seed / 100.0,
#                     ue_seed,
#                     umu_seed,
#                     utau_seed,
#                 ], dtype=float)

#                 label = "seed_dm2_{:.1f}_ue_{:.2f}_umu_{:.2f}".format(
#                     dm2_seed, ue_seed, umu_seed
#                 )

#                 res_seed = self._RunSLSQP(
#                     x_seed=x_seed,
#                     bounds=bounds,
#                     cons=cons,
#                     label=label,
#                     maxiter=300,
#                     ftol=1e-10,
#                 )

#                 candidates.append((label, res_seed))

#         # -------------------------------------------------
#         # 4. Choose best candidate
#         # -------------------------------------------------
#         best_name, best = min(candidates, key=lambda item: float(item[1].fun))
#         chi2 = float(best.fun)

#         print("\n===== selected best fit candidate =====")
#         print("best candidate =", best_name)
#         print("best chi2      =", chi2)
#         print("best x         =", best.x)
#         print("best dm2       =", best.x[0] * 100.0)
#         print("best Ue4       =", best.x[1])
#         print("best Umu4      =", best.x[2])
#         print("best Utau4     =", best.x[3])

#         return (
#             chi2,
#             {
#                 "m": best.x[0] * 100.0,
#                 "ue4": best.x[1],
#                 "umu4": best.x[2],
#                 "utau4": best.x[3],
#             },
#         )

#     def CalChi2(self, x):
#         ms = x[0] * 100
#         U_e4 = x[1]
#         U_mu4 = x[2]
#         U_tau4 = x[3]

#         self.hist.OscillateHistogram(ms, U_e4, U_mu4, U_tau4, False, False)

#         chi2, penalty = self.statistic.Chi2DataMC(
#             marginalize=self.marginalize_flux,
#             useOsc=True,
#         )

#         return chi2



class FluxFitter():
    def __init__(self,histogram,exclude,lam,usePseudo=False,useOsc=False):
        self.hist = histogram
        self.exclude = exclude
        self.lam = lam

        if usePseudo:
            self.dataHist = histogram.GetPseudoHistogram()
        else:
            self.dataHist = histogram.GetDataHistogram()

        if useOsc:
            self.mcHist = histogram.GetOscillatedHistogram()
        else:
            self.mcHist = histogram.GetMCHistogram()

        self.SolveFluxSolution()

    def SetHistogram(self,hist):
        self.hist = histogram
        if usePseudo:
            self.dataHist= histogram.GetPseudoHistogram()
        else:
            self.dataHist = histogram.GetDataHistogram()

        if useOsc:
            self.mcHist = histogram.GetOscillatedHistogram()
        else:
            self.mcHist = histogram.GetMCHistogram()

        self.SolveFluxSolution()



    def SolveFluxSolution(self):
        data = np.array(self.dataHist)[1:-1]
        mc = np.array(self.mcHist)[1:-1]
        universes = self.hist.GetFluxUniverses()
        A = self.hist.GetAMatrix()

        sliceInds = GetSliceIndices("HIST_CONFIG.json", self.exclude, self.hist.keys)

        data = slicer(data, sliceInds)
        mc   = slicer(mc, sliceInds)
        A    = slicer(A, sliceInds, axis=1)

        cov_sans = self.hist.GetCovarianceMatrix(sansFlux=True)
        cov_sliced = slicer(cov_sans, sliceInds)
        V = np.linalg.pinv(cov_sliced)

        C = data - mc
        I = np.identity(len(universes))

        L = 2 * A @ V @ C
        Q = A @ V @ A.T + I * self.lam

        self.fluxSolution = np.linalg.pinv(Q) @ L / 2

        # # After:
        # # self.fluxSolution = np.linalg.pinv(Q) @ L / 2

        # a = self.fluxSolution

        # # Objective at a = 0
        # chi2_zero = C.T @ V @ C

        # # Objective at fitted a
        # resid_fit = C - a @ A
        # chi2_fit_resid = resid_fit.T @ V @ resid_fit
        # penalty_fit = self.lam * (a @ a)
        # chi2_fit_total = chi2_fit_resid + penalty_fit

        # # Gradient check:
        # # d/da [ (C-aA)^T V (C-aA) + λ a^T a ]
        # # = -2 A V (C-aA) + 2 λ a
        # grad = -2 * A @ V @ resid_fit + 2 * self.lam * a

        # print("\n===== Flux profiling solve diagnostic =====")
        # print("selected bins       =", len(sliceInds))
        # print("A shape             =", A.shape)
        # print("V shape             =", V.shape)
        # print("lambda              =", self.lam)
        # print("Q shape             =", Q.shape)
        # print("Q condition number  =", np.linalg.cond(Q))
        # print("solution norm       =", np.sqrt(a @ a))
        # print("max |solution|      =", np.max(np.abs(a)))
        # print("chi2 at zero shift  =", chi2_zero)
        # print("chi2 prof residual  =", chi2_fit_resid)
        # print("flux penalty        =", penalty_fit)
        # print("chi2 prof total     =", chi2_fit_total)
        # print("improvement selected=", chi2_zero - chi2_fit_total)
        # print("gradient norm       =", np.linalg.norm(grad))
        # print("max |gradient|      =", np.max(np.abs(grad)))


    # def SolveFluxSolution(self):
    #     data = np.array(self.dataHist)[1:-1]
    #     mc = np.array(self.mcHist)[1:-1]
    #     universes = self.hist.GetFluxUniverses()
    #     invCov = self.hist.GetInverseCovarianceMatrix(sansFlux=True)
    #     A = self.hist.GetAMatrix()

    #     sliceInds = GetSliceIndices("HIST_CONFIG.json",self.exclude,self.hist.keys)

    #     invCov_full = self.hist.GetInverseCovarianceMatrix(sansFlux=True)
    #     V_current = slicer(invCov_full, sliceInds)

    #     cov_sans = self.hist.GetCovarianceMatrix(sansFlux=True)
    #     cov_sliced = slicer(cov_sans, sliceInds)
    #     V_alt = np.linalg.pinv(cov_sliced)

    #     print("\n===== selected covariance inverse test =====")
    #     print("V_current shape =", V_current.shape)
    #     print("V_alt shape     =", V_alt.shape)
    #     print("max abs diff    =", np.max(np.abs(V_current - V_alt)))
    #     print("rel norm diff   =", np.linalg.norm(V_current - V_alt) / np.linalg.norm(V_alt))

    #     global _PRINTED_FLUXFITTER_DIAGNOSTIC

    #     if not _PRINTED_FLUXFITTER_DIAGNOSTIC:
    #         print("\n===== FluxFitter profiling diagnostic, printed once =====")
    #         print("  exclude passed to FluxFitter =", self.exclude)
    #         print("  total stitched bins          =", len(data))
    #         print("  bins used for flux solve     =", len(sliceInds))
    #         print("  selected global bin indices  =", sliceInds)

    #         if hasattr(self.hist, "stitchKeys"):
    #             print("  stitched samples:")
    #             for k in self.hist.stitchKeys:
    #                 print("   ", k)

    #         _PRINTED_FLUXFITTER_DIAGNOSTIC = True

    #     # print("\n===== FluxFitter profiling diagnostic =====")
    #     # print("  exclude passed to FluxFitter =", self.exclude)
    #     # print("  total stitched bins          =", self.hist.GetMCHistogram().GetNbinsX())
    #     # print("  bins used for flux solve     =", len(sliceInds))

    #     # try:
    #     #     import json
    #     #     with open("HIST_CONFIG.json", "r") as f:
    #     #         hist_config = json.load(f)

    #     #     print("  samples available in HIST_CONFIG:")
    #     #     for key in self.hist.keys:
    #     #         if key not in hist_config:
    #     #             print("    {:25s}  not in HIST_CONFIG".format(key))
    #     #             continue

    #     #         # Handle a few likely config formats.
    #     #         cfg = hist_config[key]
    #     #         if isinstance(cfg, dict):
    #     #             start = cfg.get("start", cfg.get("bin_start", cfg.get("first", None)))
    #     #             end   = cfg.get("end",   cfg.get("bin_end",   cfg.get("last", None)))
    #     #         else:
    #     #             start, end = None, None

    #     #         print("    {:25s}  config = {}".format(key, cfg))

    #     #     print("  selected global bin indices for profiling:")
    #     #     print("   ", sliceInds)

    #     # except Exception as e:
    #     #     print("  Could not print HIST_CONFIG details:", e)

    #     data = slicer(data,sliceInds)
    #     mc   = slicer(mc,sliceInds)
    #     A    = slicer(A,sliceInds,axis=1)
    #     V    = slicer(invCov,sliceInds)

    #     # print("\n===== FluxFitter sliced arrays =====")
    #     # print("  sliced data bins =", len(data))
    #     # print("  sliced mc bins   =", len(mc))
    #     # print("  sliced A shape   =", A.shape)
    #     # print("  sliced V shape   =", V.shape)

    #     C = data - mc
    #     I = np.identity(len(universes))

    #     L = 2 * A @ V @ C
    #     Q = A @ V @ A.T + I * self.lam
    #     self.fluxSolution = np.linalg.inv(Q) @ L/2

    def SetFluxSolution(self,solution):
        self.fluxSolution = solution

    def GetFluxSolution(self,):
        return(self.fluxSolution)

    def MarginalizeFlux(self):
        A    = self.hist.GetAMatrix()
        solution = self.fluxSolution
        
        penalty = solution @ solution * self.lam
        new_cv = np.array(self.mcHist)[1:-1] + solution @ A
        return(new_cv,penalty)

    def ReweightToFluxSolution(self,histogram):
        mc = np.array(histogram)[1:-1]
        band = histogram.GetVertErrorBand("Flux")
        nhists = self.mcHist.GetVertErrorBand("Flux").GetNHists()
        universes = np.array([np.array(band.GetHist(l))[1:-1] for l in range(nhists)])
        cv_table = np.array([mc for l in range(len(universes))])
        A = universes - cv_table

        weights = histogram.GetCVHistoWithStatError()
        new_cv = mc + self.fluxSolution @ A

        for j in range(1,weights.GetNbinsX()+1):
            weight = weights.GetBinContent(j) / new_cv[j-1] if new_cv[j-1] != 0 else weights.GetBinContent(j)
            weights.SetBinContent(j,weight)
            weights.SetBinError(j,0)

        histogram.DivideSingle(histogram,weights)

class Statistics():
    def __init__(self,histogram,exclude="ratio",lam=1):
        self.hist = histogram
        self.exclude = exclude
        self.lam = lam
        self.nulFluxFitter = None
        self.oscFluxFitter = None

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

        inv_full = self.hist.GetInverseCovarianceMatrix(sansFlux=False)
        inv_sans = self.hist.GetInverseCovarianceMatrix(sansFlux=True)

        diff0 = data - mc
        chi2_full_no_profile = diff0.T @ inv_full @ diff0
        chi2_sans_no_profile = diff0.T @ inv_sans @ diff0

        # print("\n===== Chi2 baseline comparison =====")
        # print("useOsc =", useOsc)
        # print("marginalize =", marginalize)
        # print("chi2 full cov, no profile  =", chi2_full_no_profile)
        # print("chi2 sansFlux, no profile  =", chi2_sans_no_profile)

        invCov = inv_full
        penalty = 0

        if marginalize:
            fluxFitter = FluxFitter(self.hist, self.exclude, self.lam, usePseudo, useOsc)
            mc, penalty = fluxFitter.MarginalizeFlux()
            invCov = inv_sans

            diff_prof = data - mc
            chi2_prof_residual = diff_prof.T @ invCov @ diff_prof
            chi2_prof_total = chi2_prof_residual + penalty

            # print("\n===== Chi2 profiling comparison =====")
            # print("chi2 full cov, no profile    =", chi2_full_no_profile)
            # print("chi2 sansFlux, no profile    =", chi2_sans_no_profile)
            # print("chi2 profiled residual only  =", chi2_prof_residual)
            # print("flux penalty                 =", penalty)
            # print("chi2 profiled total          =", chi2_prof_total)
            # print("improvement vs sansFlux      =", chi2_sans_no_profile - chi2_prof_total)

            if useOsc:
                self.oscFluxFitter = fluxFitter
            else:
                self.nulFluxFitter = fluxFitter

        diff = data - mc
        chi2 = diff.T @ invCov @ diff + penalty


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
                data,
                mc,
                invCov,
                diff,
                max_bins=20
            )




        if abs(chi2) > 1e30:
            logging.error("chi2 has invalid value: {}".format(chi2))
            print("chi2 has invalid value: {}".format(chi2))
            return(-1)

        return(chi2, penalty)
    # def Chi2DataMC(self,marginalize=True,usePseudo=False,useOsc=False): 
    #     ##### Get self.hists to calculate chi2 between #####
    #     if useOsc:
    #         mcHist = self.hist.GetOscillatedHistogram()
    #     else:
    #         mcHist = self.hist.GetMCHistogram()

    #     if usePseudo:
    #         dataHist = self.hist.GetPseudoHistogram()
    #     else:
    #         dataHist = self.hist.GetDataHistogram()

    #     #We get the number of bins and make sure it's compatible with the NxN matrix given
    #     if dataHist.GetNbinsX() != mcHist.GetNbinsX():
    #         logging.error("breaking error in Chi2DataMC")
    #         logging.error("The number of bins from Data ({}) and MC ({}) histograms differ. Returning -1.".format(dataHist.GetNbinsX(),mcHist.GetNbinsX()))
    #         return(-1)

    #     mc = np.array(mcHist)[1:-1] # store MC bin contents excluding over/underflow bins
    #     data = np.array(dataHist)[1:-1] # store data bin contents excluding over/underflow bins 
       
    #     inv_full = self.hist.GetInverseCovarianceMatrix(sansFlux=False)
    #     inv_sans = self.hist.GetInverseCovarianceMatrix(sansFlux=True)

    #     diff0 = data - mc
    #     chi2_full_no_profile = diff0.T @ inv_full @ diff0
    #     chi2_sans_no_profile = diff0.T @ inv_sans @ diff0

    #     print("\n===== Chi2 baseline comparison =====")
    #     print("useOsc =", useOsc)
    #     print("marginalize =", marginalize)
    #     print("chi2 full cov, no profile  =", chi2_full_no_profile)
    #     print("chi2 sansFlux, no profile  =", chi2_sans_no_profile)

    #     invCov = self.hist.GetInverseCovarianceMatrix(sansFlux=False)
    #     penalty = 0

    #     # print("\n===== Chi2DataMC diagnostic =====")
    #     # print("  marginalize =", marginalize)
    #     # print("  useOsc      =", useOsc)
    #     # print("  usePseudo   =", usePseudo)
    #     # print("  exclude     =", self.exclude)
    #     # print("  full data bins used in final chi2 =", len(data))
    #     # print("  full mc bins used in final chi2   =", len(mc))
    #     # print("  invCov shape before profiling     =", invCov.shape)

    #     # Do we want to marginalize over the flux systematic before calculating chi2
    #     if marginalize:
    #         fluxFitter = FluxFitter(self.hist,self.exclude,self.lam,usePseudo,useOsc)
    #         oldMc = mc
    #         mc,penalty = fluxFitter.MarginalizeFlux()
    #         invCov = self.hist.GetInverseCovarianceMatrix(sansFlux=True)

    #         # print("\n===== Chi2DataMC after flux profiling =====")
    #         # print("  final chi2 still uses bins =", len(data))
    #         # print("  profiled mc bins           =", len(mc))
    #         # print("  sans-flux invCov shape     =", invCov.shape)
    #         # print("  flux penalty               =", penalty)
    #         # print("  total MC shift             =", np.sum(mc - oldMc))
    #         # print("  max abs MC shift           =", np.max(np.abs(mc - oldMc)))

    #         diff_prof = data - mc
    #         chi2_prof_residual = diff_prof.T @ invCov @ diff_prof
    #         chi2_prof_total = chi2_prof_residual + penalty

    #         print("\n===== Chi2 profiling comparison =====")
    #         print("chi2 full cov, no profile    =", chi2_full_no_profile)
    #         print("chi2 sansFlux, no profile    =", chi2_sans_no_profile)
    #         print("chi2 profiled residual only  =", chi2_prof_residual)
    #         print("flux penalty                 =", penalty)
    #         print("chi2 profiled total          =", chi2_prof_total)
    #         print("improvement vs sansFlux      =", chi2_sans_no_profile - chi2_prof_total)

    #         if useOsc:
    #             self.oscFluxFitter = fluxFitter
    #         else:
    #             self.nulFluxFitter = fluxFitter

    #     # ----- Calculate chi2 value ----= #
    #     diff = data - mc
    #     chi2 = diff.T @ invCov @ diff + penalty # @ is numpy efficient matrix multiplication

    #     # print("\n===== Final chi2 diagnostic =====")
    #     # print("  diff bins used in final chi2 =", len(diff))
    #     # print("  chi2 without/with penalty    =", chi2 - penalty, "+", penalty)
    #     # print("  final chi2                  =", chi2)

    #     if abs(chi2) > 1e30:
    #         logging.error("chi2 has invalid value: {}".format(chi2))
    #         print("chi2 has invalid value: {}".format(chi2))
    #         return(-1)

    #     return(chi2,penalty)

import os
import sys
import shutil
import logging

import ROOT
import PlotUtils
import numpy as np
from scipy import optimize

from config.AnalysisConfig import AnalysisConfig
from tools.StitchedHistogram import StitchedHistogram
from tools.Fitters import Statistics

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

ccnueroot = os.environ.get("CCNUEROOT")


class DiagnosticOscillationFitter:
    def __init__(
        self,
        histogram,
        exclude="ratio",
        lam=1,
        marginalize_flux=True,
        dm2_max=100.0,
        de_maxiter=300,
        de_popsize=15,
    ):
        self.hist = histogram
        self.exclude = exclude
        self.lam = lam
        self.marginalize_flux = marginalize_flux
        self.dm2_max = dm2_max
        self.de_maxiter = de_maxiter
        self.de_popsize = de_popsize

        self.statistic = Statistics(histogram, exclude=exclude, lam=lam)

    def CalChi2(self, x):
        dm2 = x[0] * self.dm2_max
        ue4 = x[1]
        umu4 = x[2]
        utau4 = x[3]

        self.hist.OscillateHistogram(dm2, ue4, umu4, utau4, False, False)

        chi2, penalty = self.statistic.Chi2DataMC(
            marginalize=self.marginalize_flux,
            useOsc=True,
        )

        return float(chi2)

    def DoFit(self):
        x0 = [0.0, 0.0, 0.0, 0.0]

        # x[0] is scaled to dm2 = x[0] * dm2_max
        bounds = np.array(
            [
                [0.0, 1.0],   # scaled dm2
                [0.0, 0.15],  # |Ue4|^2
                [0.0, 0.41],  # |Umu4|^2
                [0.0, 0.66],  # |Utau4|^2
            ],
            dtype=float,
        )

        # |Ue4|^2 + |Umu4|^2 + |Utau4|^2 <= 1
        cons = optimize.LinearConstraint([[0, 1, 1, 1]], -np.inf, 1)

        print("\n===== near-null SLSQP fit =====")
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

        print("near-null success =", null.success)
        print("near-null message =", null.message)
        print("near-null niter   =", getattr(null, "nit", None))
        print("near-null chi2    =", null.fun)
        print("near-null x       =", null.x)

        print("\n===== differential evolution global search =====")
        res_de = optimize.differential_evolution(
            func=self.CalChi2,
            bounds=bounds,
            constraints=cons,
            x0=x0,
            maxiter=self.de_maxiter,
            popsize=self.de_popsize,
            tol=1e-7,
            atol=1e-7,
            polish=False,
            disp=True,
            updating="immediate",
            workers=1,
        )

        print("DE success =", res_de.success)
        print("DE message =", res_de.message)
        print("DE niter   =", getattr(res_de, "nit", None))
        print("DE chi2    =", res_de.fun)
        print("DE x       =", res_de.x)

        print("\n===== polish SLSQP fit =====")
        res = optimize.minimize(
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

        print("polish success =", res.success)
        print("polish message =", res.message)
        print("polish niter   =", getattr(res, "nit", None))
        print("polish chi2    =", res.fun)
        print("polish x       =", res.x)

        best = res

        if float(null.fun) < float(best.fun):
            print("\nNear-null fit is better than global+polish.")
            best = null

        out = {
            "chi2": float(best.fun),
            "dm2": float(best.x[0] * self.dm2_max),
            "ue4": float(best.x[1]),
            "umu4": float(best.x[2]),
            "utau4": float(best.x[3]),
            "raw_x": best.x,
            "success": bool(best.success),
            "message": str(best.message),
        }

        return out


def load_sample():
    plot_tag = getattr(AnalysisConfig, "hist_config_tag", "default")

    if plot_tag in [None, "", "none"]:
        plot_tag = "default"

    filename = "rootfiles/NuE_stitched_hists_{}.root".format(plot_tag)
    file_path = "{}/oscillations/{}".format(ccnueroot, filename)

    hist_config = "HIST_CONFIG_{}.json".format(plot_tag)

    if not os.path.exists(hist_config):
        raise RuntimeError("Missing requested hist config file: {}".format(hist_config))

    shutil.copyfile(hist_config, "HIST_CONFIG.json")

    print("plot_tag    =", plot_tag)
    print("file        =", file_path)
    print("hist_config =", hist_config)
    print("Copied {} -> HIST_CONFIG.json".format(hist_config))

    sample_histogram = StitchedHistogram("sample")
    sample_histogram.Load(file_path)

    print("\n===== loaded sample =====")
    print("nbins data =", sample_histogram.GetDataHistogram().GetNbinsX())
    print("nbins mc   =", sample_histogram.GetMCHistogram().GetNbinsX())
    print("data int   =", sample_histogram.GetDataHistogram().Integral())
    print("mc int     =", sample_histogram.GetMCHistogram().Integral())

    return sample_histogram, plot_tag


def run_one_case(base_hist, case_name, exclude, dm2_max, lam, profile_flux):
    print("\n")
    print("=" * 100)
    print("CASE:", case_name)
    print("exclude       =", repr(exclude))
    print("dm2_max       =", dm2_max)
    print("lambda        =", lam)
    print("profile_flux  =", profile_flux)
    print("=" * 100)

    # Use a deep copy so each fit starts from a clean object.
    hist = base_hist.Copy()

    stat = Statistics(hist, exclude=exclude, lam=lam)
    chi2_null, penalty_null = stat.Chi2DataMC(marginalize=profile_flux)

    print("\n===== null result =====")
    print("null chi2 total =", chi2_null)
    print("null penalty    =", penalty_null)
    print("null residual   =", chi2_null - penalty_null)

    fitter = DiagnosticOscillationFitter(
        hist,
        exclude=exclude,
        lam=lam,
        marginalize_flux=profile_flux,
        dm2_max=dm2_max,
        de_maxiter=300,
        de_popsize=15,
    )

    fit = fitter.DoFit()

    delta_chi2 = chi2_null - fit["chi2"]

    print("\n===== case summary =====")
    print("case              =", case_name)
    print("exclude           =", repr(exclude))
    print("dm2_max           =", dm2_max)
    print("null chi2 total   =", chi2_null)
    print("best chi2 total   =", fit["chi2"])
    print("delta chi2        =", delta_chi2)
    print("best dm2          =", fit["dm2"])
    print("best Ue4          =", fit["ue4"])
    print("best Umu4         =", fit["umu4"])
    print("best Utau4        =", fit["utau4"])

    return {
        "case": case_name,
        "exclude": exclude,
        "dm2_max": dm2_max,
        "null_chi2": float(chi2_null),
        "null_penalty": float(penalty_null),
        "best_chi2": fit["chi2"],
        "delta_chi2": float(delta_chi2),
        "dm2": fit["dm2"],
        "ue4": fit["ue4"],
        "umu4": fit["umu4"],
        "utau4": fit["utau4"],
    }


if __name__ == "__main__":
    base_hist, plot_tag = load_sample()

    lambda_value = getattr(AnalysisConfig, "lambdaValue", 1)
    profile_flux = getattr(AnalysisConfig, "profileFlux", True)

    # Force this diagnostic to use flux profiling unless you intentionally disable it.
    profile_flux = True


    # cases = [
    #     {
    #         "case_name": "noRatio_profileAll_dm2max10",
    #         "exclude": "",
    #         "dm2_max": 10.0,
    #     },
    #     {
    #         "case_name": "noRatio_profileAll_dm2max20",
    #         "exclude": "",
    #         "dm2_max": 20.0,
    #     },
    #     {
    #         "case_name": "noRatio_profileAll_dm2max50",
    #         "exclude": "",
    #         "dm2_max": 50.0,
    #     },
    #     {
    #         "case_name": "noRatio_profileAll_dm2max100",
    #         "exclude": "",
    #         "dm2_max": 100.0,
    #     },
    # ]
    cases = [
        # Clean/nominal interpretation:
        # ratio samples are included in final chi2 but excluded from flux solve.
        {
            "case_name": "exclude_ratio_dm2max10",
            "exclude": "ratio",
            "dm2_max": 10.0,
        },
        {
            "case_name": "exclude_ratio_dm2max20",
            "exclude": "ratio",
            "dm2_max": 20.0,
        },
        {
            "case_name": "exclude_ratio_dm2max100",
            "exclude": "ratio",
            "dm2_max": 100.0,
        },

        # Diagnostic:
        # allow ratio samples to pull the flux profiling solve.
        {
            "case_name": "include_ratio_dm2max10",
            "exclude": "",
            "dm2_max": 10.0,
        },
        {
            "case_name": "include_ratio_dm2max20",
            "exclude": "",
            "dm2_max": 20.0,
        },
        {
            "case_name": "include_ratio_dm2max100",
            "exclude": "",
            "dm2_max": 100.0,
        },
    ]

    results = []

    for c in cases:
        result = run_one_case(
            base_hist=base_hist,
            case_name=c["case_name"],
            exclude=c["exclude"],
            dm2_max=c["dm2_max"],
            lam=lambda_value,
            profile_flux=profile_flux,
        )
        results.append(result)

    print("\n")
    print("=" * 100)
    print("FINAL DIAGNOSTIC SUMMARY")
    print("=" * 100)

    header = (
        "{:28s} {:12s} {:8s} {:12s} {:12s} {:12s} {:12s} {:10s} {:10s} {:10s}"
        .format(
            "case",
            "exclude",
            "dm2max",
            "null",
            "best",
            "delta",
            "dm2",
            "Ue4",
            "Umu4",
            "Utau4",
        )
    )
    print(header)
    print("-" * len(header))

    for r in results:
        print(
            "{:28s} {:12s} {:8.1f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:10.4g} {:10.4g} {:10.4g}"
            .format(
                r["case"],
                repr(r["exclude"]),
                r["dm2_max"],
                r["null_chi2"],
                r["best_chi2"],
                r["delta_chi2"],
                r["dm2"],
                r["ue4"],
                r["umu4"],
                r["utau4"],
            )
        )

    outname = "fit_diagnostic_summary_{}.csv".format(plot_tag)
    # outname = "fit_diagnostic_noRatio_profileAll_{}.csv".format(plot_tag)

    with open(outname, "w") as f:
        f.write(
            "case,exclude,dm2_max,null_chi2,null_penalty,best_chi2,delta_chi2,dm2,ue4,umu4,utau4\n"
        )
        for r in results:
            f.write(
                "{case},{exclude},{dm2_max},{null_chi2},{null_penalty},{best_chi2},{delta_chi2},{dm2},{ue4},{umu4},{utau4}\n"
                .format(**r)
            )

    print("\nWrote summary CSV:", outname)
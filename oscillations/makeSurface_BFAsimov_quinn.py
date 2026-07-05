import os
import sys
import ROOT
import PlotUtils
import numpy as np
import shutil
import argparse
import logging

np.set_printoptions(precision=5)
np.set_printoptions(linewidth=1520)
np.set_printoptions(threshold=sys.maxsize)

ccnueroot = os.environ.get("CCNUEROOT")

# ----------------------------------------------------------------------
# Parse BF-Asimov-specific args before AnalysisConfig sees sys.argv.
# This avoids modifying AnalysisConfig.
# ----------------------------------------------------------------------
_bf_parser = argparse.ArgumentParser(add_help=False)
_bf_parser.add_argument("--bf-dm2", type=float, default=None)
_bf_parser.add_argument("--bf-ue4", type=float, default=None)
_bf_parser.add_argument("--bf-umu4", type=float, default=None)
_bf_parser.add_argument("--bf-utau4", type=float, default=None)

_bf_args, _remaining_argv = _bf_parser.parse_known_args()
sys.argv = [sys.argv[0]] + _remaining_argv

from tools.PlotLibrary import HistHolder
from tools.Fitters import *
from tools.OscillationPlotTools import *
from config.AnalysisConfig import AnalysisConfig

ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)


def sanitize_tag(x):
    x = str(x)
    return (
        x.replace(".", "p")
         .replace("-", "m")
         .replace("+", "p")
         .replace("/", "_")
         .replace(",", "_")
         .replace(" ", "")
    )


def GetDefaultBestFitPoint(hist_config_tag, exclude):
    """
    Known continuous best-fit points from nominal no-throw fits.
    Used only if --bf-dm2/--bf-ue4/--bf-umu4/--bf-utau4 are not passed.
    """

    if hist_config_tag == "prodNueel_noRatio":
        return {
            "dm2": 18.967113,
            "ue4": 0.140490,
            "umu4": 0.013035,
            "utau4": 0.0,
            "chi2": 35.139265,
        }

    if hist_config_tag == "prodNueel":
        return {
            "dm2": 13.904505,
            "ue4": 0.055414,
            "umu4": 0.026505,
            "utau4": 0.0,
            "chi2": 38.769660,
        }

    raise RuntimeError(
        "No default BF point known for hist_config_tag = {}. "
        "Pass --bf-dm2 --bf-ue4 --bf-umu4 --bf-utau4 explicitly.".format(hist_config_tag)
    )


def ResolveBestFitPoint(hist_config_tag, exclude):
    default = GetDefaultBestFitPoint(hist_config_tag, exclude)

    bf_dm2 = _bf_args.bf_dm2 if _bf_args.bf_dm2 is not None else default["dm2"]
    bf_ue4 = _bf_args.bf_ue4 if _bf_args.bf_ue4 is not None else default["ue4"]
    bf_umu4 = _bf_args.bf_umu4 if _bf_args.bf_umu4 is not None else default["umu4"]
    bf_utau4 = _bf_args.bf_utau4 if _bf_args.bf_utau4 is not None else default["utau4"]

    return bf_dm2, bf_ue4, bf_umu4, bf_utau4


def BuildBFAsimovPseudoData(histogram, bf_dm2, bf_ue4, bf_umu4, bf_utau4):
    """
    Make usePseudo=True correspond to MC(best-fit), not MC(null).

    After this:
      chi2_asimov(theta) = chi2(MC(best-fit), MC(theta))
    """

    print("\n===== Building BF-Asimov pseudo-data =====")
    print("BF dm2   =", bf_dm2)
    print("BF ue4   =", bf_ue4)
    print("BF umu4  =", bf_umu4)
    print("BF utau4 =", bf_utau4)

    nominal_mc = histogram.GetMCHistogram()
    nominal_sum = nominal_mc.Integral()

    histogram.OscillateHistogram(bf_dm2, bf_ue4, bf_umu4, bf_utau4)
    bf_pseudo = histogram.GetOscillatedHistogram()
    bf_pseudo.SetName("bf_asimov_pseudodata")
    bf_pseudo.SetDirectory(0)

    histogram.SetPseudoHistogram(bf_pseudo)

    print("nominal MC sum        =", nominal_sum)
    print("BF-Asimov pseudo sum  =", bf_pseudo.Integral())
    print("BF/nominal sum ratio  =", bf_pseudo.Integral() / nominal_sum if nominal_sum != 0 else 0.0)

    # Sanity check: pseudo-data should equal the oscillated histogram
    # at the same BF hypothesis, so chi2 at BF should be ~0 up to profiling/numerics.
    stat = Statistics(histogram, exclude=AnalysisConfig.exclude, lam=AnalysisConfig.lambdaValue)
    chi2_bf, penalty_bf = stat.Chi2DataMC(
        marginalize=AnalysisConfig.profileFlux,
        useOsc=True,
        usePseudo=True
    )

    print("BF-Asimov self-check chi2 =", chi2_bf)
    print("BF-Asimov self-check penalty =", penalty_bf)


def EvalOnePoint(histogram, deltam, U_e4, U_mu4, U_tau4=0.0,
                 usePseudo=False, profileFlux=False, exclude="", lam=1):

    statistic = Statistics(histogram, exclude=exclude, lam=lam)

    chi2_null, penalty_null = statistic.Chi2DataMC(
        marginalize=profileFlux,
        useOsc=False,
        usePseudo=usePseudo
    )

    histogram.OscillateHistogram(deltam, U_e4, U_mu4, U_tau4)

    chi2_test, penalty_test = statistic.Chi2DataMC(
        marginalize=profileFlux,
        useOsc=True,
        usePseudo=usePseudo
    )

    print("\n===== Single-point oscillation test =====")
    print("profileFlux =", profileFlux)
    print("usePseudo   =", usePseudo)
    print("dm2         =", deltam)
    print("Ue4^2       =", U_e4)
    print("Umu4^2      =", U_mu4)
    print("Utau4^2     =", U_tau4)
    print("")
    print("chi2_null   = {:.6f}  penalty = {:.6f}".format(chi2_null, penalty_null))
    print("chi2_test   = {:.6f}  penalty = {:.6f}".format(chi2_test, penalty_test))
    print("test - null = {:.6f}".format(chi2_test - chi2_null))
    print("null - test = {:.6f}".format(chi2_null - chi2_test))

    return {
        "chi2_null": chi2_null,
        "chi2_test": chi2_test,
        "penalty_null": penalty_null,
        "penalty_test": penalty_test,
        "dchi2_test_minus_null": chi2_test - chi2_null,
        "dchi2_null_minus_test": chi2_null - chi2_test,
    }


def MakeSurface(histogram, outdir, deltam=1, U_tau4=0, makePlot=False,
                exclude="", lam=1, profileFlux=True, tag="default"):
    U_mu4s = 0.41 * np.logspace(-5, 0, 100)
    U_mu4s[0] = 0
    U_e4s = 0.15 * np.logspace(-4, 0, 100)
    U_e4s[0] = 0

    arrShape = (np.shape(U_mu4s)[0], np.shape(U_e4s)[0])

    asimov_surface    = np.zeros(arrShape, dtype="f")
    data_surface      = np.zeros(arrShape, dtype="f")
    data_penalties    = np.zeros(arrShape, dtype="f")
    asimov_penalties  = np.zeros(arrShape, dtype="f")
    count = 0

    statistic = Statistics(histogram, exclude=exclude, lam=lam)

    for i in range(U_mu4s.shape[0]):
        count += 1
        for j in range(U_e4s.shape[0]):
            U_mu4 = U_mu4s[i]
            U_e4  = U_e4s[j]

            histogram.OscillateHistogram(deltam, U_e4, U_mu4, U_tau4)

            # Real data surface:
            #   chi2_data(theta) = chi2(data, MC(theta))
            chi2_data, data_penalty = statistic.Chi2DataMC(
                marginalize=profileFlux,
                useOsc=True,
                usePseudo=False
            )

            # BF-Asimov surface:
            #   chi2_asimov(theta) = chi2(MC(best-fit), MC(theta))
            chi2_asimov, asimov_penalty = statistic.Chi2DataMC(
                marginalize=profileFlux,
                useOsc=True,
                usePseudo=True
            )

            data_surface[i, j] = chi2_data
            data_penalties[i, j] = data_penalty
            asimov_surface[i, j] = chi2_asimov
            asimov_penalties[i, j] = asimov_penalty

        logging.info(
            "{:.2f}% done with chi2s. Current data, BF-Asimov chi2s = {:.4f}, {:.4f}".format(
                100 * count / (U_mu4s.shape[0]),
                data_surface[i, -1],
                asimov_surface[i, -1]
            )
        )

    mode = "profiledFlux" if profileFlux else "noFluxProfile"

    np.save("{}/chi2_surface_data_{}_{}_m_{}.npy".format(outdir, tag, mode, deltam), data_surface)
    np.save("{}/chi2_surface_pseudodata_{}_{}_m_{}.npy".format(outdir, tag, mode, deltam), asimov_surface)
    np.save("{}/chi2_penalty_data_{}_{}_m_{}.npy".format(outdir, tag, mode, deltam), data_penalties)
    np.save("{}/chi2_penalty_pseudodata_{}_{}_m_{}.npy".format(outdir, tag, mode, deltam), asimov_penalties)


if __name__ == "__main__":
    hist_config_tag = getattr(AnalysisConfig, "hist_config_tag", "default")
    if hist_config_tag in [None, "", "none"]:
        hist_config_tag = "default"

    filename = "NuE_stitched_hists_{}.root".format(hist_config_tag)
    file_path = "{}/oscillations/rootfiles/{}".format(ccnueroot, filename)

    hist_config = "HIST_CONFIG_{}.json".format(hist_config_tag)
    if not os.path.exists(hist_config):
        raise RuntimeError("Missing requested hist config file: {}".format(hist_config))

    shutil.copyfile(hist_config, "HIST_CONFIG.json")

    bf_dm2, bf_ue4, bf_umu4, bf_utau4 = ResolveBestFitPoint(
        hist_config_tag,
        AnalysisConfig.exclude
    )

    bf_tag = "BFAsimov_dm2_{}_ue4_{}_umu4_{}_utau4_{}".format(
        sanitize_tag(bf_dm2),
        sanitize_tag(bf_ue4),
        sanitize_tag(bf_umu4),
        sanitize_tag(bf_utau4),
    )

    output_tag = "{}_{}".format(hist_config_tag, bf_tag)

    print("\n===== makeSurface BF-Asimov setup =====")
    print("hist_config_tag =", hist_config_tag)
    print("file            =", file_path)
    print("hist_config     =", hist_config)
    print("Copied {} -> HIST_CONFIG.json".format(hist_config))
    print("profileFlux     =", AnalysisConfig.profileFlux)
    print("exclude         =", AnalysisConfig.exclude)
    print("lambda          =", AnalysisConfig.lambdaValue)
    print("output tag      =", output_tag)
    print("BF dm2          =", bf_dm2)
    print("BF ue4          =", bf_ue4)
    print("BF umu4         =", bf_umu4)
    print("BF utau4        =", bf_utau4)

    sample_histogram = StitchedHistogram("sample")
    sample_histogram.Load(file_path)

    cat_to_exclude = AnalysisConfig.exclude_systematic
    if len(cat_to_exclude) > 0:
        sample_histogram.RemoveSystematics(cat_to_exclude)

    BuildBFAsimovPseudoData(
        sample_histogram,
        bf_dm2,
        bf_ue4,
        bf_umu4,
        bf_utau4
    )

    if AnalysisConfig.single_point:
        if AnalysisConfig.single_dm2 is None:
            raise ValueError("Need --dm2 when using --single-point")

        result = EvalOnePoint(
            sample_histogram,
            deltam=AnalysisConfig.single_dm2,
            U_e4=AnalysisConfig.single_ue4,
            U_mu4=AnalysisConfig.single_umu4,
            U_tau4=AnalysisConfig.single_utau4,
            usePseudo=AnalysisConfig.pseudodata,
            profileFlux=AnalysisConfig.profileFlux,
            lam=AnalysisConfig.lambdaValue,
            exclude=AnalysisConfig.exclude,
        )

        outname = "{}/single_point_{}_dm2_{:.5g}_ue4_{:.5g}_umu4_{:.5g}_utau4_{:.5g}_{}.txt".format(
            AnalysisConfig.output_dir,
            output_tag,
            AnalysisConfig.single_dm2,
            AnalysisConfig.single_ue4,
            AnalysisConfig.single_umu4,
            AnalysisConfig.single_utau4,
            "profiledFlux" if AnalysisConfig.profileFlux else "noFluxProfile",
        )

        with open(outname, "w") as f:
            for k, v in result.items():
                f.write("{} {}\n".format(k, v))

        print("Wrote:", outname)
        sys.exit(0)

    delta_ms = np.logspace(-1, 2, 100)
    delta_m = delta_ms[AnalysisConfig.delta_m]

    print("grid delta_m index =", AnalysisConfig.delta_m)
    print("grid delta_m value =", delta_m)

    if not AnalysisConfig.grid:
        m_toloop = np.logspace(-1, 2, 100)
        for m in m_toloop:
            print("running over delta_m^2 = {}".format(m))
            MakeSurface(
                sample_histogram,
                AnalysisConfig.output_dir,
                m,
                AnalysisConfig.U_tau4,
                lam=AnalysisConfig.lambdaValue,
                exclude=AnalysisConfig.exclude,
                profileFlux=AnalysisConfig.profileFlux,
                tag=output_tag,
            )
    else:
        MakeSurface(
            sample_histogram,
            AnalysisConfig.output_dir,
            delta_m,
            AnalysisConfig.U_tau4,
            lam=AnalysisConfig.lambdaValue,
            exclude=AnalysisConfig.exclude,
            profileFlux=AnalysisConfig.profileFlux,
            tag=output_tag,
        )
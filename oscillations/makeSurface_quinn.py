import os
import logging, sys
import ROOT
import PlotUtils
import numpy as np
from scipy import optimize, integrate
np.set_printoptions(precision=5)
np.set_printoptions(linewidth=1520)
np.set_printoptions(threshold=sys.maxsize)
ccnueroot = os.environ.get('CCNUEROOT')
import shutil
import argparse
import math
from array import array

def parse_surface_args():
    parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument("--hist-config-tag", default=None)
    parser.add_argument("--exclude", default=None)
    parser.add_argument("--lam", default=None, type=float)

    parser.add_argument("--profile-flux", action="store_true")
    parser.add_argument("--no-profile-flux", action="store_true")

    parser.add_argument("--profile-only", default=None)
    parser.add_argument("--profile-n-universes", default=None, type=int)

    args, remaining = parser.parse_known_args()
    sys.argv = [sys.argv[0]] + remaining
    return args


surface_args = parse_surface_args()


def make_safe_tag(x):
    if x is None:
        return "none"
    x = str(x).strip()
    if x == "":
        return "none"
    return x.replace(",", "-").replace("/", "-")


def build_analysis_tag(hist_config_tag, profileFlux, exclude, profile_only, profile_n_universes):
    profile_tag = "profiledFlux" if profileFlux else "noFluxProfile"

    if exclude in [None, "", "none"]:
        exclude_tag = "includeAll"
    else:
        exclude_tag = "exclude{}".format(make_safe_tag(exclude))

    if profile_n_universes is None:
        nprof_tag = "NprofAll"
    else:
        nprof_tag = "Nprof{}".format(profile_n_universes)

    profile_only_tag = ""
    if profile_only not in [None, "", "none"]:
        profile_only_tag = "_profileOnly{}".format(make_safe_tag(profile_only))

    return "{}_{}_{}_{}{}".format(
        hist_config_tag,
        profile_tag,
        exclude_tag,
        nprof_tag,
        profile_only_tag,
    )


from tools.PlotLibrary import HistHolder
from tools.Fitters import *
from tools.OscillationPlotTools import *
from config.AnalysisConfig import AnalysisConfig

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

def EvalOnePoint(
    histogram,
    deltam,
    U_e4,
    U_mu4,
    U_tau4=0.0,
    usePseudo=False,
    profileFlux=False,
    exclude="",
    lam=1,
    mask_spec=None,
    profile_only=None,
    profile_n_universes=None,
):

    # statistic = Statistics(histogram, exclude=exclude, lam=lam)
    statistic = Statistics(
        histogram,
        exclude=exclude,
        lam=lam,
        mask_spec=mask_spec,
        profile_only=profile_only,
        profile_n_universes=profile_n_universes,
    )

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
    print("profile_only =", profile_only)
    print("profile_n_universes =", profile_n_universes)
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

def MakeSurface(
    histogram,
    outdir,
    deltam=1,
    U_tau4=0,
    makePlot=False,
    exclude="",
    lam=1,
    profileFlux=True,
    tag="default",
    mask_spec=None,
    profile_only=None,
    profile_n_universes=None,
):
    U_mu4s = 0.41*np.logspace(-5,0,100)
    U_mu4s[0] = 0
    U_e4s = 0.15*np.logspace(-4,0,100)
    U_e4s[0] = 0

    arrShape = (np.shape(U_mu4s)[0],np.shape(U_e4s)[0])

    asimov_surface    = np.zeros(arrShape,dtype='f')
    data_surface      = np.zeros(arrShape,dtype='f')
    data_penalties    = np.zeros(arrShape,dtype='f')
    asimov_penalties  = np.zeros(arrShape,dtype='f')
    count = 0

    # statistic = Statistics(histogram, exclude=exclude, lam=lam)
    statistic = Statistics(
        histogram,
        exclude=exclude,
        lam=lam,
        mask_spec=mask_spec,
        profile_only=profile_only,
        profile_n_universes=profile_n_universes,
    )

    for i in range(U_mu4s.shape[0]):
        count += 1
        for j in range(U_e4s.shape[0]):
            U_mu4 = U_mu4s[i]
            U_e4  = U_e4s[j]

            histogram.OscillateHistogram(deltam, U_e4, U_mu4, U_tau4)

            chi2_data, data_penalty = statistic.Chi2DataMC(
                marginalize=profileFlux,
                useOsc=True
            )

            chi2_asimov, asimov_penalty = statistic.Chi2DataMC(
                marginalize=profileFlux,
                useOsc=True,
                usePseudo=True
            )

            if makePlot:
                from tools.PlotHistogram import PlottingContainer

                res = {"m": deltam, "ue4": U_e4, "umu4": U_mu4, "utau4": U_tau4}
                plotter = PlottingContainer("fitted_histogram", histogram)
                plotter.SetExclude(exclude)
                plotter.SetLambda(lam)
                plotter.SetMaskSpec(mask_spec)
                plotter.SetProfileOnly(profile_only)
                plotter.SetProfileNUniverses(profile_n_universes)

                plotter.PlotOscillationEffects(
                    res, "asimov",
                    useMarg=False,
                    plotSamples=False,
                    usePseudo=True
                )
                plotter.PlotOscillationEffects(
                    res, "asimov_marg",
                    useMarg=True,
                    plotSamples=False,
                    usePseudo=True
                )

            data_surface[i,j] = chi2_data
            data_penalties[i,j] = data_penalty
            asimov_surface[i,j] = chi2_asimov
            asimov_penalties[i,j] = asimov_penalty

        logging.info(
            "{:.2f}% done with chi2s. Current data, asimov chi2s = {:.4f}, {:.4f}".format(
                100*count/(U_mu4s.shape[0]),
                data_surface[i,-1],
                asimov_surface[i,-1]
            )
        )

    os.makedirs(outdir, exist_ok=True)

    np.save(
        "{}/chi2_surface_data_{}_m_{}.npy".format(outdir, tag, deltam),
        data_surface,
    )
    np.save(
        "{}/chi2_surface_pseudodata_{}_m_{}.npy".format(outdir, tag, deltam),
        asimov_surface,
    )
    np.save(
        "{}/chi2_penalty_data_{}_m_{}.npy".format(outdir, tag, deltam),
        data_penalties,
    )
    np.save(
        "{}/chi2_penalty_pseudodata_{}_m_{}.npy".format(outdir, tag, deltam),
        asimov_penalties,
    )

if __name__ == "__main__":
    hist_config_tag = surface_args.hist_config_tag
    if hist_config_tag is None:
        hist_config_tag = getattr(AnalysisConfig, "hist_config_tag", "default")

    if hist_config_tag in [None, "", "none"]:
        hist_config_tag = "default"

    filename = "NuE_stitched_hists_{}.root".format(hist_config_tag)
    file_path = "{}/oscillations/rootfiles/{}".format(ccnueroot, filename)

    hist_config = "HIST_CONFIG_{}.json".format(hist_config_tag)
    if not os.path.exists(hist_config):
        raise RuntimeError("Missing requested hist config file: {}".format(hist_config))

    shutil.copyfile(hist_config, "HIST_CONFIG.json")

    if surface_args.profile_flux:
        profileFlux = True
    elif surface_args.no_profile_flux:
        profileFlux = False
    else:
        profileFlux = getattr(AnalysisConfig, "profileFlux", True)

    lam = surface_args.lam
    if lam is None:
        lam = getattr(AnalysisConfig, "lambdaValue", 1)

    exclude = surface_args.exclude
    if exclude is None:
        exclude = getattr(AnalysisConfig, "exclude", "")

    if exclude is None:
        exclude = ""
    elif isinstance(exclude, str):
        if exclude.strip().lower() in ["none", ""]:
            exclude = ""

    profile_only = surface_args.profile_only
    if profile_only is not None and profile_only.strip().lower() in ["", "none"]:
        profile_only = None

    profile_n_universes = surface_args.profile_n_universes

    mask_spec = None

    surface_tag = build_analysis_tag(
        hist_config_tag,
        profileFlux,
        exclude,
        profile_only,
        profile_n_universes,
    )

    print("\n===== makeSurface setup =====")
    print("hist_config_tag =", hist_config_tag)
    print("file            =", file_path)
    print("hist_config     =", hist_config)
    print("Copied {} -> HIST_CONFIG.json".format(hist_config))
    print("profileFlux     =", profileFlux)
    print("exclude         =", exclude)
    print("lambda          =", lam)
    print("profile_only    =", profile_only)
    print("profile_n_universes =", profile_n_universes)
    print("surface_tag     =", surface_tag)

    sample_histogram = StitchedHistogram("sample")
    sample_histogram.Load(file_path)

    cat_to_exclude = AnalysisConfig.exclude_systematic
    if len(cat_to_exclude) > 0:
        sample_histogram.RemoveSystematics(cat_to_exclude)

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
            profileFlux=profileFlux,
            lam=lam,
            exclude=exclude,
            mask_spec=mask_spec,
            profile_only=profile_only,
            profile_n_universes=profile_n_universes,
        )

        outname = "{}/single_point_{}_dm2_{:.5g}_ue4_{:.5g}_umu4_{:.5g}_utau4_{:.5g}_{}.txt".format(
            AnalysisConfig.output_dir,
            surface_tag,
            AnalysisConfig.single_dm2,
            AnalysisConfig.single_ue4,
            AnalysisConfig.single_umu4,
            AnalysisConfig.single_utau4,
            "profiledFlux" if profileFlux else "noFluxProfile",
        )

        with open(outname, "w") as f:
            for k, v in result.items():
                f.write("{} {}\n".format(k, v))

        print("Wrote:", outname)
        sys.exit(0)

    # Existing surface behavior
    delta_ms = np.logspace(-1, 2, 100)
    delta_m = delta_ms[AnalysisConfig.delta_m]

    if not AnalysisConfig.grid:
        m_toloop = np.logspace(-1, 2, 100)
        for m in m_toloop:
            print("running over delta_m^2 = {}".format(m))
            MakeSurface(
                sample_histogram,
                AnalysisConfig.output_dir,
                m,
                AnalysisConfig.U_tau4,
                lam=lam,
                profileFlux=profileFlux,
                tag=surface_tag,
                exclude=exclude,
                mask_spec=mask_spec,
                profile_only=profile_only,
                profile_n_universes=profile_n_universes,
            )
    else:
        MakeSurface(
            sample_histogram,
            AnalysisConfig.output_dir,
            delta_m,
            AnalysisConfig.U_tau4,
            lam=lam,
            profileFlux=profileFlux,
            tag=surface_tag,
            exclude=exclude,
            mask_spec=mask_spec,
            profile_only=profile_only,
            profile_n_universes=profile_n_universes,
        )
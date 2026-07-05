import datetime as dt
import os, time, sys, math
from config.GRIDConfig import gridargs, anaargs
from tools import Utilities
import numpy as np

baseDir = os.path.dirname(os.path.abspath(__file__)) + "/../"
MacroName = baseDir.split("/")[-4]
MAT = os.environ["PLOTUTILSROOT"].split("/")[-2]
CONFIG = os.environ["PYTHONPATH"].split(":")[0].split("/")[-1]


def get_arg_value(args, key, default=None):
    if key in args:
        i = args.index(key)
        if i + 1 < len(args):
            return args[i + 1]
    return default


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


def ensure_arg(args, key, value):
    """
    If key is absent, append key/value to args.
    If key exists, leave user-provided value untouched.
    """
    if key not in args:
        args.extend([key, str(value)])
    return args


def get_default_bf(hist_config_tag):
    if hist_config_tag == "prodNueel_noRatio":
        return 18.967113, 0.140490, 0.013035, 0.0

    if hist_config_tag == "prodNueel":
        return 13.904505, 0.055414, 0.026505, 0.0

    raise RuntimeError(
        "No default BF point known for hist_config_tag = {}. "
        "Pass --bf-dm2 --bf-ue4 --bf-umu4 --bf-utau4 explicitly.".format(hist_config_tag)
    )


def createTarball(outDir):
    print("I'm inside createTarball()")
    found = os.path.isfile(outDir)
    if not found:
        cmd = "tar -czf %s -C %s %s" % (
            outDir,
            baseDir + "../",
            "{} {}".format(MacroName, MAT)
        )
        print(cmd)
        os.system(cmd)

    print("I'm done creating the tarballs")


def unpackTarball(mywrapper):
    mywrapper.write("cd $CONDOR_DIR_INPUT\n")
    mywrapper.write("source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh\n")
    mywrapper.write("spack load root@6.28.12\n")
    mywrapper.write("spack load cmake\n")
    mywrapper.write("spack load gcc\n")
    mywrapper.write("spack load fife-utils@3.7.4\n")
    mywrapper.write("spack load py-numpy\n")
    mywrapper.write("tar -xvzf {}\n".format(outdir_tarball.split("/")[-1]))
    mywrapper.write("export MINERVA_PREFIX=`pwd`/{}\n".format(MAT))
    mywrapper.write("pushd {}/bin\n".format(MAT))
    mywrapper.write("source setup.sh\n")
    mywrapper.write("popd\n")
    mywrapper.write("pushd {}\n".format(MacroName))
    mywrapper.write("source setup_ccnue.sh {}\n".format(CONFIG))
    mywrapper.write("export LD_LIBRARY_PATH=${ROOTSYS}/lib/root:${LD_LIBRARY_PATH}\n")
    mywrapper.write("source oscillations/py3env/bin/activate\n")
    mywrapper.write("popd\n")


def submitJob(tupleName, tag):
    wrapper_name = "grid_wrappers/%s/%s_wrapper.sh" % (processingID, tupleName)

    my_wrapper = open(wrapper_name, "w")
    my_wrapper.write("#!/bin/sh\n")
    unpackTarball(my_wrapper)

    my_wrapper.write("cd $CCNUEROOT/oscillations\n")
    my_wrapper.write("export USER=$(whoami)\n")
    my_wrapper.write("source py3env/bin/activate\n")
    my_wrapper.write("export PROCESS=${PROCESS:-0}\n")
    my_wrapper.write("echo HIST_CONFIG_TAG={}\n".format(hist_config_tag))
    my_wrapper.write("echo BF_TAG={}\n".format(bf_tag))
    my_wrapper.write("echo PROCESS=${PROCESS}\n")

    command = (
        "start_time=$(date +%%s)\n"
        "echo \"START_TIME=$(date)\" >> $CONDOR_DIR_LOGS/%s-%s-${PROCESS}.log\n"
        "py3env/bin/python3 makeSurface_BFAsimov_quinn.py "
        "--grid "
        "--delta_m ${PROCESS} "
        "--output $CONDOR_DIR_HISTS "
        "2>> $CONDOR_DIR_LOGS/%s-%s-${PROCESS}.err "
        "1>> $CONDOR_DIR_LOGS/%s-%s-${PROCESS}.log "
        "%s\n"
        "status=$?\n"
        "end_time=$(date +%%s)\n"
        "runtime=$((end_time - start_time))\n"
        "echo \"END_TIME=$(date)\" >> $CONDOR_DIR_LOGS/%s-%s-${PROCESS}.log\n"
        "echo \"RUNTIME_SECONDS=${runtime}\" >> $CONDOR_DIR_LOGS/%s-%s-${PROCESS}.log\n"
        "exit $status\n"
    ) % (
        tupleName, tag,
        tupleName, tag,
        tupleName, tag,
        argstring,
        tupleName, tag,
        tupleName, tag,
    )

    my_wrapper.write("echo 'Running command:'\n")
    my_wrapper.write("echo \"%s\"\n" % command.replace('"', '\\"').strip())
    my_wrapper.write(command)
    my_wrapper.close()

    print(command)
    os.system("chmod 777 %s" % wrapper_name)

    cmd = (
        "jobsub_submit "
        "--group=minerva "
        "-l '+SingularityImage=\\\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9:latest\\\"' "
        "--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC "
        "--append_condor_requirements='CpuFamily != 6' "
        "--role=Analysis "
        "--memory %dMB "
        "-f %s "
        "-d HISTS %s "
        "-d LOGS %s "
        "-N %d "
        "--expected-lifetime=%dh "
        "file://%s/%s"
    ) % (
        memory,
        outdir_tarball,
        outdir_hists,
        outdir_logs,
        njobs,
        expected_lifetime,
        os.environ["PWD"],
        wrapper_name
    )

    print(cmd)
    os.system(cmd)


if __name__ == "__main__":
    PNFS_switch = gridargs.PNFS_switch

    # Make a mutable copy because we may append default BF args.
    anaargs_local = list(anaargs)

    hist_config_tag = get_arg_value(anaargs_local, "--hist-config-tag", "default")
    if hist_config_tag in [None, "", "none"]:
        hist_config_tag = "default"

    default_bf_dm2, default_bf_ue4, default_bf_umu4, default_bf_utau4 = get_default_bf(hist_config_tag)

    bf_dm2 = get_arg_value(anaargs_local, "--bf-dm2", default_bf_dm2)
    bf_ue4 = get_arg_value(anaargs_local, "--bf-ue4", default_bf_ue4)
    bf_umu4 = get_arg_value(anaargs_local, "--bf-umu4", default_bf_umu4)
    bf_utau4 = get_arg_value(anaargs_local, "--bf-utau4", default_bf_utau4)

    # Ensure the worker receives explicit BF point.
    anaargs_local = ensure_arg(anaargs_local, "--bf-dm2", bf_dm2)
    anaargs_local = ensure_arg(anaargs_local, "--bf-ue4", bf_ue4)
    anaargs_local = ensure_arg(anaargs_local, "--bf-umu4", bf_umu4)
    anaargs_local = ensure_arg(anaargs_local, "--bf-utau4", bf_utau4)

    argstring = " ".join(anaargs_local)

    print("hist_config_tag =", hist_config_tag)
    print("argstring =", argstring)

    if "--profile-flux" in argstring:
        profile_tag = "profiledFlux"
    elif "--no-profile-flux" in argstring:
        profile_tag = "noFluxProfile"
    else:
        profile_tag = "defaultFluxMode"

    exclude = get_arg_value(anaargs_local, "--exclude", "none")
    exclude_tag = "exclude{}".format(exclude) if exclude not in ["none", "", None] else "excludeNone"

    bf_tag = "BFAsimov_dm2_{}_ue4_{}_umu4_{}_utau4_{}".format(
        sanitize_tag(bf_dm2),
        sanitize_tag(bf_ue4),
        sanitize_tag(bf_umu4),
        sanitize_tag(bf_utau4),
    )

    # Start with 1 for smoke test. Change to 100 for full run.
    # njobs = 1
    njobs = 100

    expected_lifetime = 3
    memory = gridargs.memory if gridargs.memory is not None else 1000

    processingID = "%s_%s_%s_%s_%s_%s-%s" % (
        "Oscillation_Surface",
        hist_config_tag,
        profile_tag,
        exclude_tag,
        bf_tag,
        dt.date.today(),
        dt.datetime.today().strftime("%H%M%S")
    )

    outdir_hists = "/pnfs/minerva/scratch/users/%s/%s_texts" % (
        os.environ["USER"],
        processingID
    )
    os.system("mkdir -p %s" % outdir_hists)

    outdir_logs = "/pnfs/minerva/scratch/users/%s/%s_logs" % (
        os.environ["USER"],
        processingID
    )
    os.system("mkdir -p %s" % outdir_logs)

    os.system("mkdir -p grid_wrappers/%s" % processingID)

    outdir_tarball = gridargs.tarball if gridargs.tarball else "/pnfs/minerva/resilient/tarballs/%s-%s.tar.gz" % (
        os.environ["USER"],
        processingID
    )

    createTarball(outdir_tarball)

    print("njobs =", njobs)
    print("expected_lifetime =", expected_lifetime)
    print("memory =", memory)
    print("bf_tag =", bf_tag)
    print("outdir_hists =", outdir_hists)
    print("outdir_logs =", outdir_logs)

    cmdString = "FitSpaceBFAsimov"
    tag = gridargs.ntuple_tag
    submitJob(cmdString, tag)
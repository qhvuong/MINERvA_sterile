import datetime as dt
import os, time, sys, math
from config.GRIDConfig import gridargs, anaargs
from tools import Utilities

baseDir = os.path.dirname(os.path.abspath(__file__)) + "/../"
MacroName = baseDir.split("/")[-4]
MAT = os.environ["PLOTUTILSROOT"].split("/")[-2]
CONFIG = os.environ["PYTHONPATH"].split(":")[0].split("/")[-1]
print(CONFIG)

def make_safe_tag(x):
  if x is None:
    return "none"
  x = str(x).strip()
  if x == "":
    return "none"
  return x.replace(",", "-").replace("/", "-")

def get_flag_value(args, key):
  return key in args

def get_arg_value(args, key, default="default"):
  """
  Extract value after a command-line option from anaargs.
  Example: get_arg_value(anaargs, "--hist-config-tag")
  """
  if key in args:
    i = args.index(key)
    if i + 1 < len(args):
      return args[i + 1]
  return default

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

def addBashLine(wrapper, command):
  wrapper.write("echo '---------------'\n")
  wrapper.write("echo '-------%s'\n" % command)
  wrapper.write("%s\n" % command)
  wrapper.write("echo '---------------'\n")

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
  my_wrapper.write("echo ANAARGS='{}'\n".format(argstring.replace("'", "'\"'\"'")))
  my_wrapper.write("echo PROFILE_N_UNIVERSES={}\n".format(nprof_arg))
  my_wrapper.write("echo EXCLUDE_TAG={}\n".format(exclude_tag))
  my_wrapper.write("echo ANALYSIS_TAG={}\n".format(analysis_tag))
  my_wrapper.write("echo PROCESS=${PROCESS}\n")
  my_wrapper.write("export N_TOYS=20\n")
  my_wrapper.write("echo N_TOYS=${N_TOYS}\n")
  my_wrapper.write("echo PROCESS=${PROCESS}\n")


  command = (
    "start_time=$(date +%%s)\n"
    "echo \"START_TIME=$(date)\" >> $CONDOR_DIR_LOGS/%s-%s-${PROCESS}.log\n"
    "py3env/bin/python3 fitAsimovs_quinn.py "
    "--grid "
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

  # my_wrapper.write("exit $?\n")
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

if __name__ == '__main__':
  PNFS_switch = gridargs.PNFS_switch

  argstring = " ".join(anaargs)

  hist_config_tag = get_arg_value(anaargs, "--hist-config-tag", "default")

  if hist_config_tag in [None, "", "none"]:
    hist_config_tag = "default"

  print("hist_config_tag =", hist_config_tag)

  exclude_arg = get_arg_value(anaargs, "--exclude", "none")
  if exclude_arg in [None, "", "none"]:
    exclude_tag = "includeAll"
  else:
    exclude_tag = "exclude{}".format(make_safe_tag(exclude_arg))

  profile_only_arg = get_arg_value(anaargs, "--profile-only", "none")
  if profile_only_arg in [None, "", "none"]:
    profile_only_tag = ""
  else:
    profile_only_tag = "_profileOnly{}".format(make_safe_tag(profile_only_arg))

  nprof_arg = get_arg_value(anaargs, "--profile-n-universes", "All")
  if nprof_arg in [None, "", "none", "All"]:
    nprof_tag = "NprofAll"
  else:
    nprof_tag = "Nprof{}".format(make_safe_tag(nprof_arg))

  if "--profile-flux" in argstring:
    profile_tag = "profiledFlux"
  elif "--no-profile-flux" in argstring:
    profile_tag = "noFluxProfile"
  else:
    profile_tag = "defaultFluxMode"

  analysis_tag = "{}_{}_{}{}".format(
    profile_tag,
    exclude_tag,
    nprof_tag,
    profile_only_tag
  )

  # Start with 1 for smoke test. Change to 400 for production.
  njobs = 500
  # njobs = 2

  expected_lifetime = 15

  processingID = '%s_%s_%s_%s-%s' % (
    "Asimovs",
    hist_config_tag,
    analysis_tag,
    dt.date.today(),
    dt.datetime.today().strftime("%H%M%S")
  )

  os.system("mkdir -p grid_wrappers/%s" % processingID)

  outdir_tarball = gridargs.tarball if gridargs.tarball else "/pnfs/minerva/resilient/tarballs/%s-%s.tar.gz" % (
    os.environ["USER"],
    processingID
  )

  createTarball(outdir_tarball)

  if gridargs.memory is None:
    memory = 1000
  else:
    memory = gridargs.memory

  if gridargs.count is None:
    count = 1000
  else:
    count = gridargs.count

  outdir_hists = "/pnfs/minerva/scratch/users/%s/%s_%s_%s_Asimov_dchi2s_texts" % (
    os.environ["USER"],
    hist_config_tag,
    analysis_tag,
    processingID
  )
  os.system("mkdir -p %s" % outdir_hists)

  outdir_logs = "/pnfs/minerva/scratch/users/%s/%s_%s_%s_Asimov_dchi2s_logs" % (
    os.environ["USER"],
    hist_config_tag,
    analysis_tag,
    processingID
  )
  os.system("mkdir -p %s" % outdir_logs)

  cmdString = "Asimovs"
  tag = gridargs.ntuple_tag
  submitJob(cmdString, tag)
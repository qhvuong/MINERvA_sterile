import datetime as dt
import os, time, sys, math
from config.GRIDConfig import gridargs, anaargs
from tools import Utilities

baseDir = os.path.dirname(os.path.abspath(__file__)) + "/../"
MacroName = baseDir.split("/")[-4]
MAT = os.environ["PLOTUTILSROOT"].split("/")[-2]
CONFIG = os.environ["PYTHONPATH"].split(":")[0].split("/")[-1]
print(CONFIG)


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

def remove_arg_with_value(args, keys):
  cleaned = []
  skip_next = False

  for i, a in enumerate(args):
    if skip_next:
      skip_next = False
      continue

    if a in keys:
      skip_next = True
      continue

    cleaned.append(a)

  return cleaned

def sanitize_tag(x):
  if x is None:
    return "none"
  x = str(x)
  return (
    x.replace(".", "p")
     .replace("-", "m")
     .replace("+", "p")
     .replace("/", "_")
     .replace(",", "_")
     .replace(" ", "")
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
  my_wrapper.write("echo OSC_TAG={}\n".format(osc_tag))
  my_wrapper.write("echo PROCESS=${PROCESS}\n")
  my_wrapper.write("export N_TOYS=20\n")
  my_wrapper.write("echo N_TOYS=${N_TOYS}\n")
  my_wrapper.write("echo PROCESS=${PROCESS}\n")
  my_wrapper.write("export TRUE_DM2={}\n".format(dm2))
  my_wrapper.write("export TRUE_UE4={}\n".format(ue4))
  my_wrapper.write("export TRUE_UMU4={}\n".format(umu4))
  my_wrapper.write("export TRUE_UTAU4={}\n".format(utau4))

  my_wrapper.write("echo TRUE_DM2=${TRUE_DM2}\n")
  my_wrapper.write("echo TRUE_UE4=${TRUE_UE4}\n")
  my_wrapper.write("echo TRUE_UMU4=${TRUE_UMU4}\n")
  my_wrapper.write("echo TRUE_UTAU4=${TRUE_UTAU4}\n")
  my_wrapper.write("echo ANAARGS='{}'\n".format(argstring))

  command = (
    "start_time=$(date +%%s)\n"
    "echo \"START_TIME=$(date)\" >> $CONDOR_DIR_LOGS/%s-%s-${PROCESS}.log\n"
    "py3env/bin/python3 fitOscToys_quinn.py "
    "--grid "
    "--output $CONDOR_DIR_HISTS "
    "%s "
    "2>> $CONDOR_DIR_LOGS/%s-%s-${PROCESS}.err "
    "1>> $CONDOR_DIR_LOGS/%s-%s-${PROCESS}.log\n"
    "status=$?\n"
    "end_time=$(date +%%s)\n"
    "runtime=$((end_time - start_time))\n"
    "echo \"END_TIME=$(date)\" >> $CONDOR_DIR_LOGS/%s-%s-${PROCESS}.log\n"
    "echo \"RUNTIME_SECONDS=${runtime}\" >> $CONDOR_DIR_LOGS/%s-%s-${PROCESS}.log\n"
    "exit $status\n"
  ) % (
    tupleName, tag,
    argstring,
    tupleName, tag,
    tupleName, tag,
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


if __name__ == '__main__':
  PNFS_switch = gridargs.PNFS_switch

  truth_keys = ["--dm2", "--ue4", "--umu4", "--utau4", "--truth-label"]
  fit_args_for_script = remove_arg_with_value(anaargs, truth_keys)
  argstring = " ".join(fit_args_for_script)

  hist_config_tag = get_arg_value(anaargs, "--hist-config-tag", "default")
  if hist_config_tag in [None, "", "none"]:
    hist_config_tag = "default"

  exclude = get_arg_value(anaargs, "--exclude", "none")
  if exclude in [None, "", "none"]:
    exclude_tag = "includeAll"
  else:
    exclude_tag = "exclude{}".format(sanitize_tag(exclude))

  nprof = get_arg_value(anaargs, "--profile-n-universes", "All")
  nprof_tag = "Nprof{}".format(sanitize_tag(nprof))

  profile_only = get_arg_value(anaargs, "--profile-only", "none")
  profile_only_tag = ""
  if profile_only not in [None, "", "none"]:
    profile_only_tag = "_profileOnly{}".format(sanitize_tag(profile_only))

  if "--profile-flux" in anaargs:
    profile_tag = "profiledFlux"
  elif "--no-profile-flux" in anaargs:
    profile_tag = "noFluxProfile"
  else:
    profile_tag = "defaultFluxMode"

  dm2 = get_arg_value(anaargs, "--dm2", "0")
  ue4 = get_arg_value(anaargs, "--ue4", "0")
  umu4 = get_arg_value(anaargs, "--umu4", "0")
  utau4 = get_arg_value(anaargs, "--utau4", "0")

  truth_label = get_arg_value(anaargs, "--truth-label", None)

  if truth_label not in [None, "", "none", "None"]:
    osc_tag = sanitize_tag(truth_label)
  elif float(dm2) == 0 and float(ue4) == 0 and float(umu4) == 0 and float(utau4) == 0:
    osc_tag = "null"
  else:
    osc_tag = "dm2_%s_ue4_%s_umu4_%s_utau4_%s" % (
      sanitize_tag(dm2),
      sanitize_tag(ue4),
      sanitize_tag(umu4),
      sanitize_tag(utau4),
    )

  analysis_tag = "{}_{}_{}{}".format(
    profile_tag,
    exclude_tag,
    nprof_tag,
    profile_only_tag,
  )

  # Start with 2 for smoke test. Change to 500 for production.
  # njobs = 2
  njobs = 500
  expected_lifetime = 15




  processingID = '%s_%s_%s_%s_%s-%s' % (
    "OscToys",
    hist_config_tag,
    analysis_tag,
    osc_tag,
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

  cmdString = "OscToys"
  tag = gridargs.ntuple_tag
  submitJob(cmdString, tag)
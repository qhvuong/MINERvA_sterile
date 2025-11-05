import datetime as dt
import os, time, sys, math
from config.GRIDConfig import gridargs,anaargs
from tools import Utilities

baseDir = os.path.dirname(os.path.abspath(__file__))+"/../"
MacroName = baseDir.split("/")[-4]
MAT = os.environ["PLOTUTILSROOT"].split("/")[-2]
CONFIG=os.environ["PYTHONPATH"].split(":")[0].split("/")[-1]

def createTarball(outDir):
  found = os.path.isfile(outDir)
  if(not found):
    cmd = "tar -czf %s -C %s %s"%(outDir, baseDir+"../", "{} {}".format(MacroName,MAT))# remove /Ana/CCNuE from baseDir bc want to tar the release.
    print(cmd)
    os.system(cmd)

def unpackTarball( mywrapper):
  mywrapper.write("cd $CONDOR_DIR_INPUT\n")

  mywrapper.write("source /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh\n")
  mywrapper.write("spack load root@6.28.12 arch=linux-almalinux9-x86_64_v3\n")
  mywrapper.write("spack load cmake@3.27.9%gcc@11.4.1 arch=linux-almalinux9-x86_64_v2\n")
  mywrapper.write("spack load gcc\n")
  mywrapper.write("spack load python@3.9.15\n")
  mywrapper.write("spack load ifdhc-config@2.6.20%gcc@11.4.1 arch=linux-almalinux9-x86_64_v2\n")
  # mywrapper.write("spack load ifdhc@2.8.0\n")
  mywrapper.write("spack load py-numpy@1.24.3%gcc@12.2.0\n")
  # mywrapper.write("htgettoken -i minerva --a htvaultprod.fnal.gov\n")
  mywrapper.write("htgettoken -i minerva --vaultserver htvaultprod.fnal.gov\n")
  # mywrapper.write("export BEARER_TOKEN_FILE=/run/user/`id -u`/bt_u`id -u`\n")
  # mywrapper.write("export BEARER_TOKEN_FILE=/tmp/bt_token_minerva_Analysis_`id -u`\n")
  mywrapper.write("tar -xvzf {}\n".format(outdir_tarball.split("/")[-1]))
  mywrapper.write("export MINERVA_PREFIX=`pwd`/{}\n".format(MAT))
  mywrapper.write("pushd {}/bin\n".format(MAT))
  # mywrapper.write("source setup.sh\n")
  mywrapper.write("popd\n")
  mywrapper.write("pushd {}\n".format(MacroName))
  mywrapper.write("source setup.sh\n")
  mywrapper.write("source setup_ccnue.sh {}\n".format(CONFIG))
  mywrapper.write("export LD_LIBRARY_PATH=${ROOTSYS}/lib/root:${LD_LIBRARY_PATH}\n")
  
  mywrapper.write("popd\n")

def addBashLine( wrapper , command ):

  wrapper.write("echo '---------------'\n")
  wrapper.write("echo '-------%s'\n" % command)
  wrapper.write("%s\n" % command)
  wrapper.write("echo '---------------'\n")


def submitJob( tupleName):

  wrapper_name = "grid_wrappers/%s/%s_wrapper.sh" % ( processingID , tupleName ) 
  
  my_wrapper = open(wrapper_name,"w")
  my_wrapper.write("#!/bin/sh\n")
  unpackTarball(my_wrapper)

  my_wrapper.write( "cd $CCNUEROOT/selection\n")
  my_wrapper.write( "export USER=$(whoami)\n")

  my_wrapper.write( "python eventSelection.py -p %s --grid --%s-only --ntuple_tag %s --count %d %d  --output $CONDOR_DIR_HISTS %s &> $CONDOR_DIR_LOGS/%s-${PROCESS}.log\n" % (playlist, dataSwitch, gridargs.ntuple_tag, start, count, argstring,tupleName) )
  # Add a check for expected output files
  my_wrapper.write("status=$?\n")
  my_wrapper.write("if [ $status -ne 0 ]; then\n")
  my_wrapper.write("  echo '[ERROR] eventSelection.py crashed with status $status' >&2\n")
  my_wrapper.write("  exit $status\n")
  my_wrapper.write("fi\n")

  # Check if output dir contains anything
  my_wrapper.write("n_files=$(ls $CONDOR_DIR_HISTS | wc -l)\n")
  my_wrapper.write("if [ $n_files -eq 0 ]; then\n")
  my_wrapper.write("  echo '[ERROR] No output files created in $CONDOR_DIR_HISTS' >&2\n")
  my_wrapper.write("  exit 1\n")
  my_wrapper.write("fi\n")

  my_wrapper.write("exit 0\n")

  # my_wrapper.write("exit $?\n")
  my_wrapper.close()
  
  os.system( "chmod 777 %s" % wrapper_name )
  
  # cmd = "jobsub_submit --group=minerva -l '+SingularityImage=\\\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9:latest\\\"' --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --role=Analysis --memory %dMB -f %s -d HISTS %s -d LOGS %s -N %d --expected-lifetime=%dh  file://%s/%s" % ( memory , outdir_tarball , outdir_hists , outdir_logs , njobs, 36, os.environ["PWD"] , wrapper_name )
  # --singularity-image /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9:latest
  cmd = "jobsub_submit -c has_avx2==True --group=minerva --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --role=Analysis --memory %dMB -f %s -d HISTS %s -d LOGS %s -N %d --expected-lifetime=%dh  file://%s/%s" % ( memory , outdir_tarball , outdir_hists , outdir_logs , njobs, 36, os.environ["PWD"] , wrapper_name )
  os.system(cmd)

  #cmdname = "jobsub_commands/submit_wrapper.sh"
  #jobsubcmd = open(cmdname, "w")
  #jobsubcmd.write("#!/bin/sh\n")

  #jobsubcmd.write(cmd+"\n")
  #jobsubcmd.close()
  #os.system( "chmod 777 %s" % cmdname)

if __name__ == '__main__':
  PNFS_switch = gridargs.PNFS_switch
  processingID = '%s_%s-%s' % ("CCNUE_selection", dt.date.today() , dt.datetime.today().strftime("%H%M%S") )
  outdir_hists = "/pnfs/minerva/%s/users/%s/%s_hists" % (PNFS_switch,os.environ["USER"],processingID)
  os.system( "mkdir -p %s" % outdir_hists )
  outdir_logs = "/pnfs/minerva/%s/users/%s/%s_logs" % (PNFS_switch,os.environ["USER"],processingID)
  os.system( "mkdir -p %s" % outdir_logs )
  os.system( "mkdir -p grid_wrappers/%s" % processingID )
  outdir_tarball=gridargs.tarball if gridargs.tarball else "/pnfs/minerva/resilient/tarballs/%s-%s.tar.gz" % (os.environ["USER"],processingID)
  createTarball(outdir_tarball)

  for playlist in gridargs.playlists:
    for dataSwitch in ["mc","data"]:
      if (gridargs.data_only and dataSwitch == "mc" ) or (gridargs.mc_only and dataSwitch == "data"):
        continue

      if gridargs.memory is None:
        memory = 1000 if dataSwitch == "data" else 8000
      else:
        memory = gridargs.memory
    
      if gridargs.count is None:
        count = 1000 if dataSwitch == "data" else 1
      else:
        count = gridargs.count

      argstring=" ".join(anaargs)

      start = gridargs.start
      if gridargs.njobs is None:
        end = Utilities.countLines(playlist,dataSwitch,gridargs.ntuple_tag) if gridargs.end is None else gridargs.end
        njobs = math.ceil(1.0*(end-start)/count)
      else:
        njobs = gridargs.njobs
        end = njobs*count + start

      if njobs > 1000:
        print("You should not submit more than a 1k jobs for 1 playlist.")
        print("Please merge ntuples or run more subruns per job.")
        print("there are %d ntuple files and you process %d of them per job" % ( end-start ,count))
        sys.exit(1)

      cmdString = "CCNuE-%s-%s" % (playlist,dataSwitch)
      submitJob(cmdString)

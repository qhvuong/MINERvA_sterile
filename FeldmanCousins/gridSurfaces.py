import datetime as dt
import os, time, sys, math
from config.GRIDConfig import gridargs,anaargs
from tools import Utilities
import numpy as np

baseDir = os.path.dirname(os.path.abspath(__file__))+"/../"
MacroName = baseDir.split("/")[-4]
MAT = os.environ["PLOTUTILSROOT"].split("/")[-2]
CONFIG=os.environ["PYTHONPATH"].split(":")[0].split("/")[-1]
print(CONFIG)

def createTarball(outDir):
  print("I'm inside createTarball()")
  found = os.path.isfile(outDir)
  if(not found):
    cmd = "tar -czf %s -C %s %s"%(outDir, baseDir+"../", "{} {}".format(MacroName,MAT))# remove /Ana/CCNuE from baseDir bc want to tar the release.
    print(cmd)
    os.system(cmd)

  print("I'm done creating the tarballs")

def unpackTarball( mywrapper):
  # Add lines to wrapper that wil unpack tarball; add additional setup steps here if necessary  
  mywrapper.write("cd $CONDOR_DIR_INPUT\n")
  mywrapper.write("source /cvmfs/larsoft.opensciencegrid.org/products/setup\n")
  mywrapper.write("setup root v6_22_06a -q e19:p383b:prof\n")
  mywrapper.write("tar -xvzf {}\n".format(outdir_tarball.split("/")[-1]))
  mywrapper.write("export MINERVA_PREFIX=`pwd`/{}\n".format(MAT))
  # Set up test release
  #mywrapper.write("pushd Tools/ProductionScriptsLite/cmt/\n")#assuming only one test release
  #mywrapper.write("cmt config\n")
  #mywrapper.write(". setup.sh\n")
  #mywrapper.write("popd\n")
  mywrapper.write("pushd {}/bin\n".format(MAT))
  mywrapper.write("source setup.sh\n")
  #mywrapper.write("cmt make\n")
  #mywrapper.write(". setup.sh\n")
  mywrapper.write("popd\n")
  #mywrapper.write("export LD_LIBRARY_PATH=/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/GSL/1.10/x86_64-slc7-gcc49-opt/lib:$LD_LIBRARY_PATH\n")
  mywrapper.write("pushd {}\n".format(MacroName))
  mywrapper.write("source setup_ccnue.sh {}\n".format(CONFIG))
  mywrapper.write("source FeldmanCousins/py3env/bin/activate\n")
  #mywrapper.write("make\n")
  
  mywrapper.write("popd\n")


# def copyLocalFilesToPNFS( tupleName , dest_dir):

#   if not os.path.isfile('%s/makeHists.py' % dest_dir):
#     os.system( "cp makeHists.py %s" % dest_dir )

# def clearLocalFilesFromPNFS():

#   os.system( "rm -rf /pnfs/minerva/%s/users/finer/temp/" % PNFS_switch)

def addBashLine( wrapper , command ):

  wrapper.write("echo '---------------'\n")
  wrapper.write("echo '-------%s'\n" % command)
  wrapper.write("%s\n" % command)
  wrapper.write("echo '---------------'\n")


def submitJob(tupleName,dm2,ue4):

  # Create wrapper
  wrapper_name = "grid_wrappers/%s/%s_wrapper.sh" % ( processingID , tupleName ) 
  
  my_wrapper = open(wrapper_name,"w")
  my_wrapper.write("#!/bin/sh\n")
  unpackTarball(my_wrapper)
  #addBashLine(my_wrapper,'env')
  #addBashLine(my_wrapper,'pwd')
  #addBashLine(my_wrapper,'ls')

  # Write dummy output files to get around dumb copy-back bug
  #makeDummyFile(my_wrapper,'$CONDOR_DIR_LOGS')
  #makeDummyFile(my_wrapper,'$CONDOR_DIR_HISTS')


  # This is the bash line that will be executed on the grid
  my_wrapper.write( "cd $CCNUEROOT/FeldmanCousins\n")
  my_wrapper.write( "export USER=$(whoami)\n")
  #my_wrapper.write( "export XRD_LOGLEVEL=\"Debug\"\n")
  my_wrapper.write( "source py3env/bin/activate\n")
  if gridargs.pseudodata:
    my_wrapper.write( "py3env/bin/python3 makeSurface.py --pseudodata --grid --delta_m %s --U_e4 %s --output $CONDOR_DIR_HISTS 2>> $CONDOR_DIR_LOGS/%s-%s.err 1>> $CONDOR_DIR_LOGS/%s-%s.log\n" % (dm2,ue4,argstring,tupleName,argstring,tupleName) )
  else:
    my_wrapper.write( "py3env/bin/python3 makeSurface.py --grid --delta_m %s --U_e4 %s --output $CONDOR_DIR_HISTS 2>> $CONDOR_DIR_LOGS/%s-%s.err 1>> $CONDOR_DIR_LOGS/%s-%s.log\n" % (dm2,ue4,argstring,tupleName,argstring,tupleName) )

  my_wrapper.write("exit $?\n")
  #my_wrapper.write( "python eventSelection.py -p %s --grid --%s-only --ntuple_tag %s --count %d %d  --output $CONDOR_DIR_HISTS %s \n" % (playlist, dataSwitch, gridargs.ntuple_tag, start, count, argstring) )
 
  #addBashLine(my_wrapper,'date')
  my_wrapper.close()
  
  os.system( "chmod 777 %s" % wrapper_name )
  
  #cmd = "jobsub_submit --group=minerva -l '+SingularityImage=\\\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest\\\"' --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --role=Analysis --memory %dMB -f %s/testRelease.tar.gz -d HISTS %s -d LOGS %s -r %s -N %d --expected-lifetime=%dh --cmtconfig=%s -i /cvmfs/minerva.opensciencegrid.org/minerva/software_releases/%s/ file://%s/%s" % ( memory , outdir_tarball , outdir_hists , outdir_logs , os.environ["MINERVA_RELEASE"], njobs, 12, os.environ["CMTCONFIG"],os.environ["MINERVA_RELEASE"], os.environ["PWD"] , wrapper_name )
  cmd = "jobsub_submit --group=minerva -l '+SingularityImage=\\\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest\\\"' --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --role=Analysis --memory %dMB -f %s -d HISTS %s -d LOGS %s -N %d --expected-lifetime=%dh  file://%s/%s" % ( memory , outdir_tarball , outdir_hists , outdir_logs , njobs, 36, os.environ["PWD"] , wrapper_name )
  # Copy local files to PNFS, they aren't there already
  #copyLocalFilesToPNFS(tupleName,outdir_logs) 
 
  print(cmd)

  if os.path.isfile(cmdname):
    jobsubcmd = open(cmdname, 'a')
  else:
    jobsubcmd = open(cmdname, "w")
    jobsubcmd.write("#!/bin/sh\n")

  jobsubcmd.write(cmd+"\n")
  jobsubcmd.close()
  os.system( "chmod 777 %s" % cmdname)
  #os.system(cmd)
 
  #sleepTime = 2 
  #print("Sleeping for %i seconds\n" % sleepTime)
  #time.sleep(sleepTime)

  ## # Clear files from PNFS
  ## clearLocalFilesFromPNFS()

if __name__ == '__main__':
  PNFS_switch = gridargs.PNFS_switch
  # Automatically generate unique output directory
  processingID = '%s_%s-%s' % ("Oscillation_Surface", dt.date.today() , dt.datetime.today().strftime("%H%M%S") )
  if gridargs.pseudodata:
      processingID = '%s_%s-%s-pseudodata' % ("Oscillation_Surface", dt.date.today() , dt.datetime.today().strftime("%H%M%S") )
  outdir_hists = "/pnfs/minerva/scratch/users/%s/%s_Oscillation_Surface_texts" % (os.environ["USER"],processingID)
  os.system( "mkdir -p %s" % outdir_hists )
  outdir_logs = "/pnfs/minerva/scratch/users/%s/%s_Oscillation_Surface_logs" % (os.environ["USER"],processingID)
  os.system( "mkdir -p %s" % outdir_logs )
  os.system( "mkdir -p grid_wrappers/%s" % processingID )
  outdir_tarball=gridargs.tarball if gridargs.tarball else "/pnfs/minerva/resilient/tarballs/rhowell-%s.tar.gz" % (processingID)
  createTarball(outdir_tarball)

  if gridargs.memory is None:
    memory = 1000
  else:
    memory = gridargs.memory

  if gridargs.count is None:
    count = 1000
  else:
    count = gridargs.count

  argstring=" ".join(anaargs)

  njobs = 1
  cmdname = "grid_surface.sh"
  if gridargs.pseudodata:
      cmdname = "grid_surface_pseudodata.sh"

  if os.path.exists(cmdname):
    os.system( "rm {}".format(cmdname))

  m_toloop = np.logspace(0,2,60)
  U_e4s = 0.15*np.logspace(-2.2,0,60)
  U_e4s[0] = 0

  for m in m_toloop:
    for ue4 in U_e4s:
      cmdString = "FitSpace_{:.3f}_{:.4f}".format(m,ue4)
      submitJob(cmdString,str(m),str(ue4))

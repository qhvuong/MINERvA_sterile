# source /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh

# spack load root@6.28.12 arch=linux-almalinux9-x86_64_v3

# #spack load cmake@3.27.9%gcc@11.4.1 arch=linux-almalinux9-x86_64_v2
# spack load cmake@3.27.9%gcc@11.4.1 arch=linux-almalinux9-x86_64_v3
# spack load gcc
# spack load python@3.9.15
# #spack load ifdhc-config@2.6.20%gcc@11.4.1 arch=linux-almalinux9-x86_64_v2
# spack load ifdhc-config@2.6.20%gcc@11.4.1 arch=linux-almalinux9-x86_64_v3
# spack load py-numpy@1.24.3%gcc@12.2.0
# export JOBSUB_GROUP=minerva

# # htgettoken -a htvaultprod.fnal.gov -i minerva
# htgettoken -i minerva --vaultserver htvaultprod.fnal.gov
# export BEARER_TOKEN_FILE=/run/user/`id -u`/bt_u`id -u`
# #export BEARER_TOKEN_FILE=/tmp/bt_token_minerva_Analysis_`id -u`

# source /exp/minerva/app/users/$USER/MAT_AL9/opt/bin/setup.sh
# # source /exp/minerva/app/users/$USER/MAT_AL9/CC-NuE-XSec/setup.sh
# export LD_LIBRARY_PATH=${ROOTSYS}/lib/root:${LD_LIBRARY_PATH}
source /exp/minerva/app/users/$USER/MAT_AL9/CC-NuE-XSec/setup_ccnumu.sh

cd $CCNUEROOT

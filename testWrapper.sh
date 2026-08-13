#!/bin/bash
set -euo pipefail

export XRD_NETWORKSTACK=IPv4
export JOBSUB_GROUP=minerva
umask 002

config=${1:?need config}
playlist=${2:?need playlist}
sample=${3:?need sample type: mc or data}
count_n=${4:?need count number per job}
selection_tag=${5:?need selection tag}
out_tag=${6:?need output tag}

# shift 6
# extra_args=("$@")
# count_start=$(( PROCESS * count_n ))

start_override=""
if [[ $# -ge 7 && "$7" != --* ]]; then
  start_override="$7"
  shift 7
else
  shift 6
fi
extra_args=("$@")

if [[ -n "$start_override" ]]; then
  count_start="$start_override"
else
  count_start=$(( PROCESS * count_n ))
fi

proc_id=$(( count_start / count_n ))
export PROCESS="$proc_id"

hist_out="$CONDOR_DIR_HISTS"
log_out="$CONDOR_DIR_LOGS"
# logfile="${log_out}/selection-${playlist}-${sample}-${PROCESS}.log"
logfile="${log_out}/selection-${playlist}-${sample}-${proc_id}.log"

echo "=== JOB START ==="
echo "Start time: $(date)"
echo "config=$config"
echo "playlist=$playlist"
echo "sample=$sample"
echo "count_start=$count_start"
echo "count_n=$count_n"
echo "selection_tag=$selection_tag"
echo "out_tag=$out_tag"
echo "logfile=$logfile"
echo

cd "$CONDOR_DIR_INPUT"

source /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh
spack load root@6.28.12 arch=linux-almalinux9-x86_64_v3
spack load cmake@3.27.9%gcc@11.4.1 arch=linux-almalinux9-x86_64_v3
spack load gcc
spack load python@3.9.15
spack load ifdhc-config@2.6.20%gcc@11.4.1 arch=linux-almalinux9-x86_64_v3
spack load py-numpy@1.24.3%gcc@12.2.0

cd "$(dirname "$INPUT_TAR_FILE")"

export MINERVA_PREFIX="$(pwd)/opt"
source opt/bin/setup.sh
export LD_LIBRARY_PATH="${ROOTSYS}/lib/root:${LD_LIBRARY_PATH}"
source CC-NuE-XSec/setup_ccnue.sh "$config"

cd "$CCNUEROOT/selection"
export USER=$(whoami)

case "$sample" in
  mc)   sample_flag="--mc_only" ;;
  data) sample_flag="--data_only" ;;
  *) echo "[ERROR] sample must be mc or data, got '$sample'" >&2; exit 2 ;;
esac

( /usr/bin/time -v python eventSelection.py \
    --playlist "$playlist" \
    --grid \
    $sample_flag \
    --ntuple_tag MAD \
    --count "$count_start" "$count_n" \
    --output "$hist_out" \
    --selection_tag "$selection_tag" \
    "${extra_args[@]}" \
) &> "$logfile"

status=$?
if [ $status -ne 0 ]; then
  echo "[ERROR] eventSelection.py crashed with status $status" >&2
  exit $status
fi

echo "=== JOB END ==="
echo "End time: $(date)"
#!/usr/bin/env bash
set -euo pipefail

tag=${1:?You must provide a selection_tag}
count=${2:-30}
exclude=${3:-None}
sideband=${4:-None}

## Running format: ./submitFHC.sh selection_tag count exclude sideband(s)

logdir="/exp/minerva/data/users/qvuong/runningNotes"
mkdir -p "${logdir}"

logfile="${logdir}/${tag}_$(date +%Y-%m-%d_%H%M%S).txt"
echo "Logging to $logfile"
exec > >(tee -a "$logfile") 2>&1

timestamp=$(date +%Y-%m-%d_%H%M%S)
tarball="/pnfs/minerva/resilient/tarballs/${USER}-CCNUE_selection_${tag}_${timestamp}.tar.gz"

echo "Using shared tarball: ${tarball}"
echo "Selection tag: ${tag}"
echo "MC count/job : ${count}"
echo "Exclude      : ${exclude}"
echo "Sideband     : ${sideband}"

EXCLUDE_ARGS=()
if [[ -n "${exclude}" && "${exclude}" != "None" ]]; then
  EXCLUDE_ARGS=(--exclude "${exclude}")
fi

USE_SIDEBAND_ARGS=(--use-sideband)
if [[ -n "${sideband}" && "${sideband}" != "None" ]]; then
  USE_SIDEBAND_ARGS=(--use-sideband "${@:4}")
fi

COMMON_ARGS=(
  --ntuple_tag MAD
  "${EXCLUDE_ARGS[@]}"
  "${USE_SIDEBAND_ARGS[@]}"
  --truth
  --cal_POT
  --selection_tag "${tag}"
  --tarball "${tarball}"
)

for name in 1 7 9 13C 13C_2p2h; do
# for name in 7; do
  cmd=(
    python selection/gridSelection.py
    --playlist le${name}_p6
    "${COMMON_ARGS[@]}"
    --count "${count}"
    --mc_only
  )

  echo "--------------------------------------------------"
  echo "Running command:"
  printf ' %q' "${cmd[@]}"
  echo
  echo "--------------------------------------------------"

  "${cmd[@]}"
done

# for name in 1 7 9 13A 13B 13C 13D 13E; do
# # for name in 7; do
#   cmd=(
#     python selection/gridSelection.py
#     --playlist le${name}_p6
#     "${COMMON_ARGS[@]}"
#     --count 200
#     --data_only
#   )

#   echo "--------------------------------------------------"
#   echo "Running command:"
#   printf ' %q' "${cmd[@]}"
#   echo
#   echo "--------------------------------------------------"

#   "${cmd[@]}"
# done
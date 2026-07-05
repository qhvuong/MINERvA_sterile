#!/usr/bin/env bash
set -euo pipefail

tag=${1:?You must provide a selection_tag}
count=${2:-30}
exclude_universes=${3:-None}
sideband=${4:-None}

## Running format:
## ./submitRHC.sh selection_tag count exclude_universes sideband(s)
##
## Example:
## ./submitRHC.sh test 30 all
## ./submitRHC.sh test 30 None my_sideband

logdir="/exp/minerva/data/users/qvuong/runningNotes"
mkdir -p "${logdir}"

logfile="${logdir}/${tag}_$(date +%Y-%m-%d_%H%M%S).txt"
echo "Logging to $logfile"
exec > >(tee -a "$logfile") 2>&1

timestamp=$(date +%Y-%m-%d_%H%M%S)
tarball="/pnfs/minerva/resilient/tarballs/${USER}-CCNUE_selection_${tag}_${timestamp}.tar.gz"

echo "Using shared tarball: ${tarball}"
echo "Selection tag      : ${tag}"
echo "MC count/job       : ${count}"
echo "Exclude universes  : ${exclude_universes}"
echo "Sideband           : ${sideband}"

EXCLUDE_UNIVERSE_ARGS=()
if [[ -n "${exclude_universes}" && "${exclude_universes}" != "None" ]]; then
  EXCLUDE_UNIVERSE_ARGS=(--exclude-universes "${exclude_universes}")
fi

USE_SIDEBAND_ARGS=(--use-sideband)
if [[ -n "${sideband}" && "${sideband}" != "None" ]]; then
  USE_SIDEBAND_ARGS=(--use-sideband "${@:4}")
fi

COMMON_ARGS=(
  --ntuple_tag MAD
  "${USE_SIDEBAND_ARGS[@]}"
  --truth
  --cal_POT
  --selection_tag "${tag}"
  "${EXCLUDE_UNIVERSE_ARGS[@]}"
  --tarball "${tarball}"
)

for name in 5; do
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

# for name in 5; do
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
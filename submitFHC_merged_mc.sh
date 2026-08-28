#!/usr/bin/env bash
set -euo pipefail

tag=${1:?You must provide a selection_tag}
exclude_universes=${2:-None}
sideband=${3:-None}

## Running format:
## ./submitFHC.sh selection_tag exclude_universes sideband(s)
##
## Example:
## ./submitFHC.sh test all
## ./submitFHC.sh test None my_sideband

logdir="/exp/minerva/data/users/qvuong/runningNotes"
mkdir -p "${logdir}"

logfile="${logdir}/${tag}_$(date +%Y-%m-%d_%H%M%S).txt"
echo "Logging to $logfile"
exec > >(tee -a "$logfile") 2>&1

timestamp=$(date +%Y-%m-%d_%H%M%S)
tarball="/pnfs/minerva/resilient/tarballs/${USER}-CCNUE_selection_${tag}_${timestamp}.tar.gz"

echo "Using shared tarball: ${tarball}"
echo "Selection tag      : ${tag}"
echo "Exclude universes  : ${exclude_universes}"
echo "Sideband           : ${sideband}"

EXCLUDE_UNIVERSE_ARGS=()
if [[ -n "${exclude_universes}" && "${exclude_universes}" != "None" ]]; then
  EXCLUDE_UNIVERSE_ARGS=(--exclude-universes "${exclude_universes}")
fi

USE_SIDEBAND_ARGS=()
if [[ -n "${sideband}" && "${sideband}" != "None" ]]; then
  USE_SIDEBAND_ARGS=(--use-sideband "${@:3}")
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

for name in 1 7 9 13C 13C_2p2h; do
# for name in 1 7 9 13C; do
# for name in 7 9 13C; do
# for name in 1; do
  cmd=(
    python selection/gridSelection.py
    --playlist le${name}_p8
    "${COMMON_ARGS[@]}"
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
#     --playlist le${name}_p8
#     "${COMMON_ARGS[@]}"
#     --count 20
#     --data_only
#   )

#   echo "--------------------------------------------------"
#   echo "Running command:"
#   printf ' %q' "${cmd[@]}"
#   echo
#   echo "--------------------------------------------------"

#   "${cmd[@]}"
# done
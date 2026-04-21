#!/usr/bin/env bash
set -euo pipefail

tag=${1:?You must provide a selection_tag}
count=${2:-50}
sideband=${3:-None}

logfile="runningNotes/${tag}_$(date +%Y-%m-%d_%H%M%S).txt"
echo "Logging to $logfile"
exec > >(tee -a "$logfile") 2>&1

USE_SIDEBAND_ARGS=(--use-sideband)
if [[ -n "${sideband}" && "${sideband}" != "None" ]]; then
  USE_SIDEBAND_ARGS=(--use-sideband "${@:3}")
fi

for name in 1 7 9 13C 13C_2p2h; do
  cmd=(
    python selection/gridSelection.py
    --playlist le${name}_p6
    --ntuple_tag MAD
    "${USE_SIDEBAND_ARGS[@]}"
    --truth
    --cal_POT
    --selection_tag "${tag}"
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
#   cmd=(
#     python selection/gridSelection.py
#     --playlist le${name}_p6
#     --ntuple_tag MAD
#     "${USE_SIDEBAND_ARGS[@]}"
#     --truth
#     --cal_POT
#     --selection_tag "${tag}"
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
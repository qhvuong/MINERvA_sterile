#!/usr/bin/env bash
set -euo pipefail

tag=${1:?You must provide a selection_tag}
sideband=${2:-None}

logfile="runningNotes/${tag}_$(date +%Y-%m-%d_%H%M%S).txt"
echo "Logging to $logfile"
exec > >(tee -a "$logfile") 2>&1

# Build the sideband part exactly how you want:
# - Always include --use-sideband
# - If sideband is None/empty, add nothing after it
# - Otherwise append the provided sideband tokens
USE_SIDEBAND_ARGS=(--use-sideband)
if [[ -n "${sideband}" && "${sideband}" != "None" ]]; then
  # If you pass multiple sidebands, put them after tag like:
  #   ./submitAll.sh test dEdX Eavail Etheta
  # so we take all args from $2 onward
  USE_SIDEBAND_ARGS=(--use-sideband "${@:2}")
fi

for name in 1 7 9 13C; do
# for name in 7 9 13C; do
  cmd=(
    python selection/gridSelection.py
    --playlist le${name}_p6
    --ntuple_tag MAD
    "${USE_SIDEBAND_ARGS[@]}"
    --truth
    --cal_POT
    --selection_tag "${tag}"
    --count 12
  )

  echo "--------------------------------------------------"
  echo "Running command:"
  printf ' %q' "${cmd[@]}"
  echo
  echo "--------------------------------------------------"

  "${cmd[@]}"
done

for name in 13A 13B 13D 13E; do
  cmd=(
    python selection/gridSelection.py
    --playlist le${name}_p6
    --ntuple_tag MAD
    "${USE_SIDEBAND_ARGS[@]}"
    --truth
    --cal_POT
    --selection_tag "${tag}"
    --count 12
    --data_only
  )

  echo "--------------------------------------------------"
  echo "Running command:"
  printf ' %q' "${cmd[@]}"
  echo
  echo "--------------------------------------------------"

  "${cmd[@]}"
done

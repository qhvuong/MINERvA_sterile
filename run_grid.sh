#!/bin/bash

MODE=""
SELECTION_TAG=""

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --mc)
      MODE="mc"
      shift
      ;;
    --data)
      MODE="data"
      shift
      ;;
    --selection_tag)
      SELECTION_TAG="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

if [[ -z "$MODE" ]]; then
  echo "ERROR: must specify --mc or --data"
  exit 1
fi

if [[ -z "$SELECTION_TAG" ]]; then
  echo "ERROR: must specify --selection_tag"
  exit 1
fi

# Add "no" prefix safely
if [[ "$SELECTION_TAG" == no* ]]; then
  FINAL_SELECTION_TAG="$SELECTION_TAG"
else
  FINAL_SELECTION_TAG="no${SELECTION_TAG}"
#   FINAL_SELECTION_TAG="${SELECTION_TAG}"
fi

# Playlist lists
if [[ "$MODE" == "mc" ]]; then
  PLAYLISTS=(1 7 9)
  MODE_FLAG="--mc_only"
else
  PLAYLISTS=(1 7 9 13A 13B 13C 13D 13E)
  MODE_FLAG="--data_only"
fi

# Print commands (do NOT run)
for name in "${PLAYLISTS[@]}"; do
  CMD="python selection/gridSelection.py \
    --playlist le${name}_p6 \
    --ntuple_tag MAD \
    --use-sideband \
    --truth \
    --cal_POT \
    --selection_tag ${FINAL_SELECTION_TAG} \
    ${MODE_FLAG} \
    --count 12"

  echo "Running command:"
  echo "${CMD}"
  echo "=================================================="

  eval ${CMD}
done

echo "All gridSelection jobs for ${SELECTION_TAG} in ${MODE} finished."
echo "=================================================="

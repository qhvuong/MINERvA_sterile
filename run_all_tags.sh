#!/bin/bash

SELECTION_TAGS=(
  HasNoBackExitingTracks
  HasTracks
  Vertex_Z
  Vertex_Apothem
  EMLikeTrackScore
  DSCalVisE
  ODCalVisE
  DeadTime
  Afterpulsing
  NonMIPClusFrac
  TransverseGapScore
  HasNoVertexMismatch
  StartPointVertexMultiplicity
  VertexTrackMultiplicity
  LeptonEnergy
  Q2
  InverseEtheta
  MeanFrontdEdX
  Eavail
)

# Pause 30 mins (1800 seconds)
PAUSE_SECS=$((30 * 60))

for tag in "${SELECTION_TAGS[@]}"; do
  
  echo "=============================================="
  echo "Running selection tag: $tag"
  echo "=============================================="

  for mode in mc data; do
    ./run_grid.sh --$mode --selection_tag "$tag"
  done

  echo "Finished tag = $tag"
  echo "Sleeping for ${PAUSE_SECS} seconds (30 minutes)..."
  sleep "$PAUSE_SECS"

done

echo "=============================================="
echo "All SELECTION_TAGS completed!"
echo "=============================================="

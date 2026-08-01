#!/bin/bash

script=make_flat_beamfocus_weights.py
filedir=../file_option
outdir=flat_beamfocus_outputs/FHC

mkdir -p "$outdir"

for playlist in 1A 1B 1C 1D 1E 1F 1G 1L 1M 1N 1O 1P; do
  filelist="${filedir}/playlist_mcme${playlist}p8.txt"
  csv="${outdir}/diagnostic_ME_FHC_${playlist}.csv"
  npz="${outdir}/accum_ME_FHC_${playlist}.npz"
  log="${outdir}/make_flat_ME_FHC_${playlist}.log"

  echo "Running playlist ${playlist}"

  python "$script" \
    --mode FHC \
    --filelist "$filelist" \
    --n_universes 1000 \
    --accumulator_output "$npz" \
    --output "$csv" \
    > "$log" 2>&1

  status=$?

  if [ "$status" -ne 0 ]; then
    echo "ERROR: playlist ${playlist} failed. See ${log}"
  else
    echo "Finished playlist ${playlist}"
  fi
done
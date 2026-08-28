#!/bin/bash

set -euo pipefail

SCRIPT="make_flat_beamfocus_weights.py"
NUNIV=1000

BASE_OUT="/exp/minerva/data/users/qvuong/flux_studies/flat_beamfocus_outputs_FV"
FHC_OUT="${BASE_OUT}/FHC"
RHC_OUT="${BASE_OUT}/RHC"

mkdir -p "${FHC_OUT}" "${RHC_OUT}"

# ============================================================
# FHC playlists
# ============================================================

FHC_PLAYLISTS=(
    # 1A
    # 1B
    # 1C
    # 1D
    # 1E
    # 1F
    # 1G
    # 1L
    # 1M
    # 1N
    # 1O
    # 1P
)

for pl in "${FHC_PLAYLISTS[@]}"; do
    echo
    echo "============================================"
    echo "Processing ME FHC ${pl}"
    echo "============================================"

    FILELIST="../file_option/playlist_mcme${pl}p8.txt"
    CSV="${FHC_OUT}/diagnostic_ME_FHC_${pl}.csv"
    NPZ="${FHC_OUT}/accum_ME_FHC_${pl}.npz"
    LOG="${FHC_OUT}/make_flat_ME_FHC_${pl}.log"

    if python -u "${SCRIPT}" \
        --mode FHC \
        --filelist "${FILELIST}" \
        --n_universes "${NUNIV}" \
        --accumulator_output "${NPZ}" \
        --output "${CSV}" \
        2>&1 | tee "${LOG}"
    then
        echo "Finished playlist ${pl}"
    else
        echo "ERROR: playlist ${pl} failed. See ${LOG}"
    fi
done


# ============================================================
# RHC playlists
# ============================================================

RHC_PLAYLISTS=(
    # 5A
    # 6A
    # 6B
    # 6C
    # 6D
    # 6E
    # 6F
    # 6G
    6H
    6I
    6J
)

for pl in "${RHC_PLAYLISTS[@]}"; do
    echo
    echo "============================================"
    echo "Processing ME RHC ${pl}"
    echo "============================================"

    FILELIST="../file_option/playlist_mcme${pl}p8.txt"
    CSV="${RHC_OUT}/diagnostic_ME_RHC_${pl}.csv"
    NPZ="${RHC_OUT}/accum_ME_RHC_${pl}.npz"
    LOG="${RHC_OUT}/make_flat_ME_RHC_${pl}.log"

    if python -u "${SCRIPT}" \
        --mode RHC \
        --filelist "${FILELIST}" \
        --n_universes "${NUNIV}" \
        --accumulator_output "${NPZ}" \
        --output "${CSV}" \
        2>&1 | tee "${LOG}"
    then
        echo "Finished playlist ${pl}"
    else
        echo "ERROR: playlist ${pl} failed. See ${LOG}"
    fi
done

echo
echo "============================================"
echo "All BeamFocus playlist jobs completed"
echo "============================================"
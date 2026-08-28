#!/bin/bash

set -euo pipefail

# PLAYLISTS=("LE1" "LE5" "LE13C")
PLAYLISTS=("LE13C")

SPECIES=(
    # "mu -1"
    # "mu 1"
    # "e -1"
    "e 1"
)

for PLAYLIST in "${PLAYLISTS[@]}"; do
    for SPEC in "${SPECIES[@]}"; do

        read -r FLAVOR HELICITY <<< "${SPEC}"

        echo
        echo "============================================================"
        echo "Running ${PLAYLIST} ${FLAVOR} ${HELICITY}"
        echo "============================================================"
        echo

        ./run_compute_flux.sh "${PLAYLIST}" "${FLAVOR}" "${HELICITY}"

        echo
        echo "============================================================"
        echo "Refreshing token after ${PLAYLIST} ${FLAVOR} ${HELICITY}"
        echo "============================================================"
        echo

        . ../../getToken.sh

    done
done

echo
echo "============================================================"
echo "All flux productions completed."
echo "============================================================"
#!/bin/bash

set -euo pipefail

# Get a fresh token before starting anything
. ../../getToken.sh

# ------------------------------------------------------------
# Keep the token alive while this script is running
# ------------------------------------------------------------
refresh_token_loop() {
    while true; do
        sleep 1800   # 30 minutes

        echo
        echo "============================================================"
        echo "Refreshing token in background"
        echo "============================================================"
        echo

        . ../../getToken.sh || {
            echo "[WARNING] Token refresh failed at $(date)"
        }
    done
}

refresh_token_loop &
TOKEN_REFRESH_PID=$!

# Kill the background refresher when this script exits/stops/fails
trap 'kill ${TOKEN_REFRESH_PID} 2>/dev/null || true' EXIT INT TERM


PLAYLISTS=(
    # "ME1A"
    # "ME1B"
    # "ME1C"
    # "ME1D"
    # "ME1E"
    # "ME1F"
    # "ME1G"
    # "ME1L"
    # "ME1M"
    "ME1N"
    # "ME1O"
    # "ME1P"

    # "ME5A"

    # "ME6A"
    # "ME6B"
    # "ME6C"
    # "ME6D"
    # "ME6E"
    # "ME6F"
    # "ME6G"
    # "ME6H"
    # "ME6I"
    # "ME6J"
)

SPECIES=(
    # "mu -1"
    "mu 1"
    # "e -1"
    # "e 1"
)

for PLAYLIST in "${PLAYLISTS[@]}"; do
    for SPEC in "${SPECIES[@]}"; do

        read -r FLAVOR HELICITY <<< "${SPEC}"

        echo
        echo "============================================================"
        echo "Running ${PLAYLIST} ${FLAVOR} ${HELICITY}"
        echo "============================================================"
        echo

        ./run_compute_flux.sh \
            "${PLAYLIST}" \
            "${FLAVOR}" \
            "${HELICITY}"

    done
done

echo
echo "============================================================"
echo "All ME flux productions completed."
echo "============================================================"
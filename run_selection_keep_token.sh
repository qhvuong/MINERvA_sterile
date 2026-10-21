#!/usr/bin/env bash
set -u

TOKEN_INTERVAL=1800

refresh_token() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Refreshing token..."
    . ../getToken.sh
}

keep_token_alive() {
    while true; do
        sleep "${TOKEN_INTERVAL}"
        refresh_token
    done
}

# Get a fresh token before starting.
refresh_token

# Start periodic refresh in the background.
keep_token_alive &
TOKEN_PID=$!

# Make sure the background loop is killed when the job ends.
cleanup() {
    kill "${TOKEN_PID}" 2>/dev/null || true
    wait "${TOKEN_PID}" 2>/dev/null || true
}
trap cleanup EXIT INT TERM

python -u selection/eventSelection.py \
    --playlist le9_p8 \
    --ntuple_tag MAD \
    --use-sideband \
    --truth \
    --cal_POT \
    --selection_tag test \
    --mc_only \
    --exclude-universes none \
    --count 3 1
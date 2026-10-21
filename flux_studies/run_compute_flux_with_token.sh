#!/bin/bash

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 PLAYLIST FLAVOR HELICITY"
    echo
    echo "Examples:"
    echo "  $0 LE1 e -1"
    echo "  $0 LE5 e 1"
    echo "  $0 LE13C mu -1"
    exit 1
fi

PLAYLIST="$1"
FLAVOR="$2"
HELICITY="$3"

TOKEN_REFRESH_SECONDS=3600

cleanup() {
    if [[ -n "${TOKEN_PID:-}" ]]; then
        kill "${TOKEN_PID}" 2>/dev/null || true
        wait "${TOKEN_PID}" 2>/dev/null || true
    fi
}

trap cleanup EXIT INT TERM

echo
echo "============================================================"
echo "Refreshing token before starting"
echo "============================================================"
echo

. ../../getToken.sh

(
    while true; do

        sleep "${TOKEN_REFRESH_SECONDS}"

        echo
        echo "============================================================"
        echo "Refreshing token"
        echo "============================================================"
        echo

        . ../../getToken.sh

    done
) &

TOKEN_PID=$!

echo
echo "============================================================"
echo "Running flux production"
echo "Playlist : ${PLAYLIST}"
echo "Flavor   : ${FLAVOR}"
echo "Helicity : ${HELICITY}"
echo "============================================================"
echo

./run_compute_flux.sh \
    "${PLAYLIST}" \
    "${FLAVOR}" \
    "${HELICITY}"

echo
echo "============================================================"
echo "Flux production completed"
echo "============================================================"
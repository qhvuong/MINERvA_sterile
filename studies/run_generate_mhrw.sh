#!/bin/bash

PLAYLIST=$1

if [ -z "$PLAYLIST" ]; then
    echo "Usage: $0 <playlist>"
    exit 1
fi

# Initial token
. ../getToken.sh

# Keep refreshing while the job runs
(
    while true; do
        sleep 1800
        echo "[TOKEN] refreshing at $(date)"
        . ../getToken.sh
    done
) &

TOKEN_PID=$!

python studies/generate_mhrw.py \
    --playlist "$PLAYLIST" \
    --ntuple_tag MAD

STATUS=$?

kill $TOKEN_PID 2>/dev/null
wait $TOKEN_PID 2>/dev/null

exit $STATUS
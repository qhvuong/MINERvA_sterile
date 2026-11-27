#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -lt 2 ]; then
  echo 'Usage: bash selection/run_fsi_probability_balance.sh OUTPUT_PREFIX [normal selection arguments]' >&2
  exit 2
fi
balance_prefix=$1
shift
for balance_item in "${balance_prefix}_reco" "${balance_prefix}_truth" \
                    "${balance_prefix}_reco.log" "${balance_prefix}_truth.log"; do
  if [ -e "$balance_item" ]; then
    echo "Output exists; choose a new prefix: $balance_item" >&2
    exit 2
  fi
done
mkdir -p "$(dirname "$balance_prefix")"
python -u selection/study_fsi_probability_balance.py "$@" \
  --balance-tree reco --balance-output "${balance_prefix}_reco" \
  2>&1 | tee "${balance_prefix}_reco.log"
python -u selection/study_fsi_probability_balance.py "$@" \
  --balance-tree truth --balance-output "${balance_prefix}_truth" \
  2>&1 | tee "${balance_prefix}_truth.log"

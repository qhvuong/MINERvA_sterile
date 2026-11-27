#!/usr/bin/env bash
# Run from the analysis repository root with the intended beam config active.
set -euo pipefail
if [ "$#" -lt 2 ]; then
  echo 'Usage: bash selection/run_dis_weight_pair.sh OUTPUT_PREFIX [normal selection options]' >&2
  exit 2
fi
pair_prefix=$1
shift
for pair_path in "${pair_prefix}_430" "${pair_prefix}_431" "${pair_prefix}_comparison" \
                 "${pair_prefix}_430.log" "${pair_prefix}_431.log"; do
  if [ -e "$pair_path" ]; then
    echo "Output already exists; use a new prefix: $pair_path" >&2
    exit 2
  fi
done
mkdir -p "$(dirname "$pair_prefix")"
python -u selection/dump_dis_weight_components.py "$@" \
  --pair-tune v4p3p0 --pair-output "${pair_prefix}_430" \
  2>&1 | tee "${pair_prefix}_430.log"
python -u selection/dump_dis_weight_components.py "$@" \
  --pair-tune v4p3p1 --pair-output "${pair_prefix}_431" \
  2>&1 | tee "${pair_prefix}_431.log"
python selection/compare_dis_weight_components.py \
  --v430 "${pair_prefix}_430" --v431 "${pair_prefix}_431" \
  --output "${pair_prefix}_comparison"

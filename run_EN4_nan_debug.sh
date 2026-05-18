#!/usr/bin/env bash
set -euo pipefail

start=${1:-1440}
end=${2:-1469}
tag=${3:-debug_EN4_single}
playlist=${4:-le5_p6}

logdir="/exp/minerva/data/users/qvuong/runningNotes/EN4_nan_debug_${playlist}_${start}_${end}_$(date +%Y-%m-%d_%H%M%S)"
mkdir -p "${logdir}"

summary="${logdir}/SUMMARY.txt"
: > "${summary}"

echo "Debug log dir: ${logdir}"
echo "Playlist: ${playlist}" | tee -a "${summary}"
echo "Range: ${start}-${end}" | tee -a "${summary}"
echo "Tag base: ${tag}" | tee -a "${summary}"
echo "" | tee -a "${summary}"

for i in $(seq "${start}" "${end}"); do
  logfile="${logdir}/file_${i}.log"

  echo "=================================================="
  echo "Running file index ${i}"
  echo "Log: ${logfile}"
  echo "=================================================="

  set +e
  python selection/eventSelection.py \
    --playlist "${playlist}" \
    --ntuple_tag MAD \
    --use-sideband \
    --truth \
    --cal_POT \
    --selection_tag "${tag}_${i}" \
    --mc_only \
    --count "${i}" 1 \
    > "${logfile}" 2>&1
  status=$?
  set -e

  bad_fill_count=$(grep -c '^\[EN4_BAD_FILL\]' "${logfile}" || true)
  bad_bin_count=$(grep -c '^\[EN4_BAD_BIN\]' "${logfile}" || true)
  branch_error_count=$(grep -c "Can't find branch" "${logfile}" || true)

  if [[ "${status}" -ne 0 || "${bad_fill_count}" -ne 0 || "${bad_bin_count}" -ne 0 || "${branch_error_count}" -ne 0 ]]; then
    echo "[BAD] index=${i} status=${status} bad_fill=${bad_fill_count} bad_bin=${bad_bin_count} branch_errors=${branch_error_count} log=${logfile}" | tee -a "${summary}"

    {
      echo "---- first EN4 bad lines for index ${i} ----"
      grep -E '^\[EN4_BAD_FILL\]|^\[EN4_BAD_BIN\]|Can'\''t find branch|Traceback|Error|Exception' "${logfile}" | head -80
      echo ""
    } >> "${summary}"
  else
    echo "[OK] index=${i} status=${status}" | tee -a "${summary}"
  fi
done

echo ""
echo "================ FINAL SUMMARY ================"
cat "${summary}"
echo "==============================================="
echo "Wrote summary: ${summary}"
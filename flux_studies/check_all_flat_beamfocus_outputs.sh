#!/bin/bash

set -u

OUTDIR="/exp/minerva/data/users/qvuong/flux_studies/producedFluxes_p8_LE"

FHC_CSV="flat_beamfocus_outputs/FHC/flat_beamfocus_ME_FHC.csv"
RHC_CSV="flat_beamfocus_outputs/RHC/flat_beamfocus_ME_RHC.csv"

CHECKER="check_flat_beamfocus_output.py"
HIST="flux_E_cvweighted"

failures=0

run_check() {
  local root_file="$1"
  local csv_file="$2"
  local column="$3"
  local offset="$4"

  echo
  echo "============================================================"
  echo "Checking: ${root_file}"
  echo "CSV:      ${csv_file}"
  echo "Column:   ${column}"
  echo "Offset:   ${offset}"
  echo "============================================================"

  python "${CHECKER}" \
    --root-file "${root_file}" \
    --csv "${csv_file}" \
    --column "${column}" \
    --offset "${offset}" \
    --hist-name "${HIST}" \
    --print-universes 3

  status=$?

  if [ ${status} -ne 0 ]; then
    echo "FAILED: ${root_file}"
    failures=$((failures + 1))
  else
    echo "PASSED: ${root_file}"
  fi
}


# ------------------------------------------------------------
# FHC LE1: right sign is numu, so test nue, nuebar, and numubar.
# ------------------------------------------------------------

run_check \
  "${OUTDIR}/LE1_nue.root" \
  "${FHC_CSV}" \
  "nu_e" \
  509

run_check \
  "${OUTDIR}/LE1_nuebar.root" \
  "${FHC_CSV}" \
  "nubar_e" \
  509

run_check \
  "${OUTDIR}/LE1_numubar.root" \
  "${FHC_CSV}" \
  "nubar_mu" \
  509


# ------------------------------------------------------------
# FHC LE13: right sign is numu.
# ------------------------------------------------------------

run_check \
  "${OUTDIR}/LE13_nue.root" \
  "${FHC_CSV}" \
  "nu_e" \
  509

run_check \
  "${OUTDIR}/LE13_nuebar.root" \
  "${FHC_CSV}" \
  "nubar_e" \
  509

run_check \
  "${OUTDIR}/LE13_numubar.root" \
  "${FHC_CSV}" \
  "nubar_mu" \
  509


# ------------------------------------------------------------
# RHC LE5: right sign is numubar, so test numu, nue, and nuebar.
# ------------------------------------------------------------

run_check \
  "${OUTDIR}/LE5_numu.root" \
  "${RHC_CSV}" \
  "nu_mu" \
  811

run_check \
  "${OUTDIR}/LE5_nue.root" \
  "${RHC_CSV}" \
  "nu_e" \
  811

run_check \
  "${OUTDIR}/LE5_nuebar.root" \
  "${RHC_CSV}" \
  "nubar_e" \
  811


echo
echo "============================================================"

if [ ${failures} -eq 0 ]; then
  echo "ALL NON-RIGHT-SIGN BEAMFOCUS CHECKS PASSED"
else
  echo "${failures} CHECK(S) FAILED"
fi

echo "============================================================"

exit ${failures}
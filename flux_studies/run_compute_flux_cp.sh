#!/bin/bash

set -euo pipefail

COMPUTE_FLUX="compute_flux_cp.py"
NUNIV=1000

FHC_OFFSET=509
RHC_OFFSET=811

FHC_FLAT="/exp/minerva/data/users/qvuong/flux_studies/flat_beamfocus_outputs_FV/FHC/flat_beamfocus_FHC.csv"
RHC_FLAT="/exp/minerva/data/users/qvuong/flux_studies/flat_beamfocus_outputs_FV/RHC/flat_beamfocus_RHC.csv"

FILEDIR="../file_option"

usage() {
    echo "Usage:"
    echo "  $0 PLAYLIST FLAVOR HELICITY"
    echo
    echo "PLAYLIST: LE1, LE5, LE13C, ME1*, ME5*, ME6*"
    echo "FLAVOR  : mu, e"
    echo "HELICITY: -1 for neutrino, +1 for antineutrino"
    exit 1
}

if [ "$#" -ne 3 ]; then
    usage
fi

PLAYLIST="$1"
FLAVOR="$2"
HELICITY="$3"

# ============================================================
# Playlist configuration
# ============================================================

case "${PLAYLIST}" in

    # -------------------------
    # Low-energy FHC
    # -------------------------
    LE1|LE13C)
        ERA="LE"
        MODE="FHC"
        MODE_FLAG="--FHC"
        OFFSET="${FHC_OFFSET}"
        FLAT_CSV="${FHC_FLAT}"
        OUTDIR="/exp/minerva/data/users/qvuong/flux_studies/producedFluxes_p8_LE_FV"
        ;;

    # -------------------------
    # Low-energy RHC
    # -------------------------
    LE5)
        ERA="LE"
        MODE="RHC"
        MODE_FLAG="--RHC"
        OFFSET="${RHC_OFFSET}"
        FLAT_CSV="${RHC_FLAT}"
        OUTDIR="/exp/minerva/data/users/qvuong/flux_studies/producedFluxes_p8_LE_FV"
        ;;

    # -------------------------
    # Medium-energy FHC
    # -------------------------
    ME1*)
        ERA="ME"
        MODE="FHC"
        MODE_FLAG="--FHC"
        OFFSET="${FHC_OFFSET}"
        FLAT_CSV=""
        OUTDIR="/exp/minerva/data/users/qvuong/flux_studies/producedFluxes_p8_ME_offsetPPFX"
        ;;

    # -------------------------
    # Medium-energy RHC
    # -------------------------
    ME5*|ME6*)
        ERA="ME"
        MODE="RHC"
        MODE_FLAG="--RHC"
        OFFSET="${RHC_OFFSET}"
        FLAT_CSV=""
        OUTDIR="/exp/minerva/data/users/qvuong/flux_studies/producedFluxes_p8_ME_offsetPPFX"
        ;;

    *)
        echo "ERROR: unsupported playlist: ${PLAYLIST}"
        exit 1
        ;;
esac

mkdir -p "${OUTDIR}"


# ============================================================
# Filelist
# ============================================================

FILELIST=$(find "${FILEDIR}" \
    -maxdepth 1 \
    -type f \
    -iname "playlist_mc${PLAYLIST}p8.txt" \
    -print -quit)

if [ -z "${FILELIST}" ]; then
    echo "ERROR: could not find filelist for ${PLAYLIST}"
    exit 1
fi



# ============================================================
# Flavor/helicity -> PDG
# ============================================================

case "${FLAVOR}" in
    mu)
        ABS_PDG=14
        ;;
    e)
        ABS_PDG=12
        ;;
    *)
        echo "ERROR: FLAVOR must be mu or e"
        exit 1
        ;;
esac

if [ "${HELICITY}" -eq -1 ]; then
    PDG="${ABS_PDG}"
elif [ "${HELICITY}" -eq 1 ]; then
    PDG="-$ABS_PDG"
else
    echo "ERROR: HELICITY must be -1 or +1"
    exit 1
fi


# ============================================================
# Determine right-sign muon species
#
# FHC RS = nu_mu
# RHC RS = nubar_mu
# ============================================================

if [ "${MODE}" = "FHC" ] && [ "${PDG}" -eq 14 ]; then
    RIGHT_SIGN=true
elif [ "${MODE}" = "RHC" ] && [ "${PDG}" -eq -14 ]; then
    RIGHT_SIGN=true
else
    RIGHT_SIGN=false
fi



# ============================================================
# Determine BeamFocus source
#
# LE right-sign muon:
#   native tuple BeamFocus
#
# LE non-right-sign:
#   flat ME-derived BeamFocus
#
# ME, all species:
#   native tuple BeamFocus
# ============================================================

USE_FLAT_BEAMFOCUS=false

if [ "${ERA}" = "LE" ] && [ "${RIGHT_SIGN}" = false ]; then
    USE_FLAT_BEAMFOCUS=true
fi




# ============================================================
# Species name
# ============================================================

case "${PDG}" in
    14)
        SPECIES="numu"
        ;;
    -14)
        SPECIES="numubar"
        ;;
    12)
        SPECIES="nue"
        ;;
    -12)
        SPECIES="nuebar"
        ;;
    *)
        echo "ERROR: unsupported PDG ${PDG}"
        exit 1
        ;;
esac


OUTPUT="${OUTDIR}/${PLAYLIST}_${SPECIES}.root"

# ============================================================
# Print configuration
# ============================================================

echo
echo "============================================"
echo "compute_flux configuration"
echo "============================================"
echo "Playlist          : ${PLAYLIST}"
echo "Mode              : ${MODE}"
echo "Flavor            : ${FLAVOR}"
echo "Helicity          : ${HELICITY}"
echo "PDG               : ${PDG}"
echo "Beam offset       : ${OFFSET}"
echo "Right-sign muon   : ${RIGHT_SIGN}"
echo "Filelist          : ${FILELIST}"
echo "Output            : ${OUTPUT}"

echo "Era               : ${ERA}"

if [ "${USE_FLAT_BEAMFOCUS}" = true ]; then
    echo "BeamFocus         : flat ME-derived weights"
    echo "Flat BF CSV       : ${FLAT_CSV}"
else
    echo "BeamFocus         : native tuple weights"
fi

echo "============================================"
echo


# ============================================================
# Build compute_flux command
# ============================================================

CMD=(
    python -u "${COMPUTE_FLUX}"
    "${MODE_FLAG}"
    --nu_flavor "${FLAVOR}"
    --nu_helicity "${HELICITY}"
    --filelist "${FILELIST}"
    --use_meta_tree
    --calc_errors
    --n_universes "${NUNIV}"
    --beam_universe_offset "${OFFSET}"
)

# Flat BeamFocus is used only for LE non-right-sign species.
# All ME species use their native tuple BeamFocus weights.
if [ "${USE_FLAT_BEAMFOCUS}" = true ]; then
    CMD+=(
        --flat_beamfocus_csv "${FLAT_CSV}"
    )
fi

CMD+=(
    -o "${OUTPUT}"
)


# ============================================================
# Run
# ============================================================

echo "Running:"
printf ' %q' "${CMD[@]}"
echo
echo

LOGDIR="${OUTDIR}/logs"
mkdir -p "${LOGDIR}"

LOG="${LOGDIR}/${PLAYLIST}_${SPECIES}.log"

echo "Log: ${LOG}"
echo

"${CMD[@]}" 2>&1 | tee "${LOG}"
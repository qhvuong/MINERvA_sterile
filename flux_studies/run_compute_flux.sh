#!/bin/bash

set -euo pipefail

COMPUTE_FLUX="compute_flux.py"
NUNIV=1000

FHC_OFFSET=509
RHC_OFFSET=811

FHC_FLAT="/exp/minerva/data/users/qvuong/flux_studies/flat_beamfocus_outputs_FV/FHC/flat_beamfocus_FHC.csv"
RHC_FLAT="/exp/minerva/data/users/qvuong/flux_studies/flat_beamfocus_outputs_FV/RHC/flat_beamfocus_RHC.csv"

BEAMFOCUS_RESPONSE="/exp/minerva/app/users/qvuong/FluxAna/LE_BeamFocus_universes_nue.root"

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
        OUTDIR="/exp/minerva/data/users/qvuong/flux_studies/producedFluxes_p8_LE_newElectronResponse"
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
        OUTDIR="/exp/minerva/data/users/qvuong/flux_studies/producedFluxes_p8_LE_newElectronResponse"
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
# Determine BeamFocus source
#
# LE:
#   right-sign muon:
#       FHC numu / RHC numubar
#       -> native tuple BeamFocus
#
#   right-sign electron:
#       FHC nue / RHC nuebar
#       -> regenerated energy-dependent BeamFocus response
#
#   wrong-sign species:
#       -> flat BeamFocus weights
#
# ME:
#   all species
#       -> native tuple BeamFocus
# ============================================================

RIGHT_SIGN_MUON=false
RIGHT_SIGN_ELECTRON=false

if [ "${MODE}" = "FHC" ]; then

    if [ "${PDG}" -eq 14 ]; then
        RIGHT_SIGN_MUON=true
    elif [ "${PDG}" -eq 12 ]; then
        RIGHT_SIGN_ELECTRON=true
    fi

elif [ "${MODE}" = "RHC" ]; then

    if [ "${PDG}" -eq -14 ]; then
        RIGHT_SIGN_MUON=true
    elif [ "${PDG}" -eq -12 ]; then
        RIGHT_SIGN_ELECTRON=true
    fi

fi


BEAMFOCUS_SOURCE="tuple"

if [ "${ERA}" = "LE" ]; then

    if [ "${RIGHT_SIGN_MUON}" = true ]; then
        BEAMFOCUS_SOURCE="tuple"

    elif [ "${RIGHT_SIGN_ELECTRON}" = true ]; then
        BEAMFOCUS_SOURCE="response"

    else
        BEAMFOCUS_SOURCE="flat"
    fi

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
echo "Right-sign muon   : ${RIGHT_SIGN_MUON}"
echo "Right-sign e      : ${RIGHT_SIGN_ELECTRON}"
echo "Filelist          : ${FILELIST}"
echo "Output            : ${OUTPUT}"

echo "Era               : ${ERA}"

case "${BEAMFOCUS_SOURCE}" in

    tuple)
        echo "BeamFocus         : native tuple weights"
        ;;

    flat)
        echo "BeamFocus         : flat ME-derived weights"
        echo "Flat BF CSV       : ${FLAT_CSV}"
        ;;

    response)
        echo "BeamFocus         : regenerated energy-dependent response"
        echo "BF response file  : ${BEAMFOCUS_RESPONSE}"
        ;;

esac

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
    --output_file "${OUTPUT}"
)

# LE BeamFocus:
#   right-sign muon     -> native tuple
#   right-sign electron -> regenerated energy-dependent response
#   wrong-sign species  -> flat weights
#
# ME:
#   all species -> native tuple
case "${BEAMFOCUS_SOURCE}" in

    tuple)
        # No external BeamFocus input needed.
        ;;

    flat)
        CMD+=(
            --flat_beamfocus_csv "${FLAT_CSV}"
        )
        ;;

    response)
        CMD+=(
            --beamfocus_response_file "${BEAMFOCUS_RESPONSE}"
        )
        ;;

esac

if [ "${BEAMFOCUS_SOURCE}" = "response" ]; then
    if [ ! -f "${BEAMFOCUS_RESPONSE}" ]; then
        echo "ERROR: BeamFocus response file not found:"
        echo "  ${BEAMFOCUS_RESPONSE}"
        exit 1
    fi
fi

if [ "${BEAMFOCUS_SOURCE}" = "flat" ]; then
    if [ ! -f "${FLAT_CSV}" ]; then
        echo "ERROR: Flat BeamFocus CSV not found:"
        echo "  ${FLAT_CSV}"
        exit 1
    fi
fi

case "${BEAMFOCUS_SOURCE}" in
    tuple|flat|response)
        ;;
    *)
        echo "ERROR: unknown BeamFocus source: ${BEAMFOCUS_SOURCE}"
        exit 1
        ;;
esac

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
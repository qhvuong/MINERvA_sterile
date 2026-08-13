#!/bin/bash
set -euo pipefail

# -----------------------------
# User settings
# -----------------------------
selection_tag=${1:?You must provide a selection_tag}
config=${2:?You must provide a config}
count=${3:-200}
sideband=${4:-None}
tarball_tag=${5:-}
start_override=${6:-}
rerun_playlist=${7:-}
rerun_sample=${8:-}

USE_SIDEBAND_ARGS=()
case "$sideband" in
  None|none)
    ;;
  dEdX)
    USE_SIDEBAND_ARGS=(--use-sideband dEdX)
    ;;
  all)
    USE_SIDEBAND_ARGS=(--use-sideband dEdX Etheta Eavail)
    ;;
  *)
    echo "[ERROR] sideband must be None, dEdX, or all" >&2
    exit 1
    ;;
esac

ntuple_tag="MAD"
wrapper="$PWD/testWrapper.sh"


# optional per-playlist/sample overrides
declare -A COUNT_MAP=(
  ["le1_p6|data"]=513       # *6
  ["le5_p6|data"]=528       # *12
  ["le7_p6|data"]=500       # *1
  ["le9_p6|data"]=600       # *1
  ["le13A_p6|data"]=800     # *1
  ["le13B_p6|data"]=511     # *4
  ["le13C_p6|data"]=512     # *14
  ["le13D_p6|data"]=400     # *1
  ["le13E_p6|data"]=578     # *5


  ["le1_p6|mc"]=326         # *30
  ["le5_p6|mc"]=328         # *30
  ["le7_p6|mc"]=316         # *3
  ["le9_p6|mc"]=330         # *3
  ["le13C_p6|mc"]=327       # *35
)


mkdir -p /exp/minerva/data/users/qvuong/tarballs

if [[ -n "$tarball_tag" ]]; then
  out_tag="$tarball_tag"
  tarball_local="/exp/minerva/data/users/qvuong/tarballs/GridTarball_${out_tag}.tar.gz"

  [[ -f "$tarball_local" ]] || {
    echo "[ERROR] requested existing tarball not found: $tarball_local" >&2
    exit 1
  }

  echo "[INFO] Reusing existing tarball: $tarball_local"
else
  out_tag=$(date +%Y-%m-%d_%H%M%S)
  tarball_local="/exp/minerva/data/users/qvuong/tarballs/GridTarball_${out_tag}.tar.gz"

  echo "[INFO] Creating new tarball: $tarball_local"
  tar -czf "$tarball_local" \
    -C /exp/minerva/app/users/qvuong/MAT_AL9 \
    CC-NuE-XSec opt
fi

tarball="dropbox://${tarball_local}"
echo "[INFO] out_tag       = $out_tag"
echo "[INFO] tarball_local = $tarball_local"
echo "[INFO] tarball       = $tarball"


log_suffix="_logs"
if [[ -n "$start_override" ]]; then
  log_suffix="_logs_rerun"
fi

logdir="/exp/minerva/data/users/qvuong/runningNotes"
mkdir -p "$logdir"
logfile="${logdir}/${config}_${selection_tag}_${out_tag}${log_suffix}.txt"
echo "Logging to $logfile"
exec > >(tee -a "$logfile") 2>&1


DRY_RUN=0
# out_tag="DRYRUN"
# tarball="dropbox:///DUMMY/TARBALL.tar.gz"


# optional global start/end
start=0
end_override=""   # leave empty to use full count from Utilities.countLines

# optional forced njobs
force_njobs=""    # leave empty to auto-compute

# optional forced memory
force_memory=""   # leave empty for default: data=1500, mc=2000


# extra args passed through to eventSelection.py
extra_args=(
  "${USE_SIDEBAND_ARGS[@]}"
  --truth
  --cal_POT
)


# playlists with both mc and data
both_names=(
  1
  7
  9
  13C
)

# playlists with data only
data_only_names=(
  13A
  13B
  13D
  13E
)

REPO_DIR="/exp/minerva/app/users/qvuong/MAT_AL9/CC-NuE-XSec"
POT_JSON="${REPO_DIR}/configs/POT.json"
# -----------------------------
# helper: count files
# -----------------------------
get_playlist_file() {
  local playlist="$1"
  local sample="$2"

  python - <<PY
import json
import os
import sys

pot_json = "${POT_JSON}"
repo_dir = "${REPO_DIR}"
playlist = "${playlist}"
sample = "${sample}"
ntuple_tag = "${ntuple_tag}"

with open(pot_json) as f:
    data = json.load(f)

try:
    relpath = data[playlist][sample][ntuple_tag]["playlist_location"]
except KeyError as e:
    print(f"ERROR: missing key in POT.json for playlist={playlist}, sample={sample}, ntuple_tag={ntuple_tag}: {e}", file=sys.stderr)
    sys.exit(1)

print(os.path.join(repo_dir, relpath))
PY
}

count_lines() {
  local playlist="$1"
  local sample="$2"
  local playlist_file

  playlist_file=$(get_playlist_file "$playlist" "$sample")

  [[ -f "$playlist_file" ]] || {
    echo "[ERROR] playlist file not found: $playlist_file" >&2
    exit 1
  }

  wc -l < "$playlist_file"
}

# -----------------------------
# helper: submit one sample
# -----------------------------
submit_one() {
  local playlist="$1"
  local sample="$2"

  local memory end n_to_process njobs cmd_string
  local key="${playlist}|${sample}"
  # local count_per_job="$count"
  # local count_per_job="${COUNT_MAP[$key]:-$count}"
  local count_per_job

  if [[ "$config" == "CCNuMu" || "$config" == "CCAntiNuMu" ]]; then
    count_per_job="$count"
  else
    count_per_job="${COUNT_MAP[$key]:-$count}"
  fi

  if [[ -n "$force_memory" ]]; then
    memory="$force_memory"
  else
    if [[ "$sample" == "data" ]]; then
      memory=1500
    else
      memory=2000
    fi
  fi

  if [[ -n "$end_override" ]]; then
    end="$end_override"
  else
    end=$(count_lines "$playlist" "$sample")
  fi

  [[ "$end" =~ ^[0-9]+$ ]] || {
    echo "[ERROR] bad count returned for playlist=$playlist sample=$sample : $end" >&2
    exit 1
  }

  n_to_process=$(( end - start ))

  if (( n_to_process <= 0 )); then
    echo "[WARN] skipping playlist=$playlist sample=$sample start=$start end=$end"
    return 0
  fi

  if [[ -n "$force_njobs" ]]; then
    njobs="$force_njobs"
    end=$(( start + njobs * count_per_job ))
  else
    njobs=$(( (n_to_process + count_per_job - 1) / count_per_job ))
  fi

  if [[ -n "$rerun_playlist" && "$playlist" != "$rerun_playlist" ]]; then
    return 0
  fi
  if [[ -n "$rerun_sample" && "$sample" != "$rerun_sample" ]]; then
    return 0
  fi

  if [[ -n "$start_override" ]]; then
    njobs=1
  fi

  if (( njobs > 1000 )); then
    echo "You should not submit more than 1000 jobs for one playlist/sample."
    echo "Please merge ntuples or run more files per job."
    echo "There are $(( end - start )) ntuple files and you process $count_per_job per job."
    exit 1
  fi

  cmd_string="${config}-${playlist}-${sample}"

  echo "=================================================="
  echo "Submitting: $cmd_string"
  echo "  playlist  = $playlist"
  echo "  sample    = $sample"
  echo "  total     = $(( end - start ))"
  echo "  count/job = $count_per_job"
  echo "  njobs     = $njobs"
  echo "  memory    = ${memory}MB"
  echo "  out_tag   = $out_tag"
  echo "=================================================="

  # cmd=(
  #   jobsub_submit
  #   -G minerva
  #   --memory="${memory}MB"
  #   --disk=5GB
  #   --expected-lifetime=8h
  #   -N "$njobs"
  #   -d HISTS "/pnfs/minerva/scratch/users/qvuong/${config}/${out_tag}_${selection_tag}_hists"
  #   -d LOGS "/pnfs/minerva/scratch/users/qvuong/${config}/${out_tag}_${selection_tag}_logs"
  #   --tar-file-name "$tarball"
  #   --use-cvmfs-dropbox
  #   file://"$wrapper"
  #   "$config" "$playlist" "$sample" "$count_per_job" "$selection_tag" "$out_tag"
  #   "${extra_args[@]}"
  # )

  cmd=(
    jobsub_submit
    -G minerva
    --memory="${memory}MB"
    --disk=5GB
    --expected-lifetime=8h
    -N "$njobs"
    -d HISTS "/pnfs/minerva/scratch/users/qvuong/${config}/${out_tag}_${selection_tag}_hists"
    # -d LOGS "/pnfs/minerva/scratch/users/qvuong/${config}/${out_tag}_${selection_tag}_logs_rerun"
    -d LOGS "/pnfs/minerva/scratch/users/qvuong/${config}/${out_tag}_${selection_tag}${log_suffix}"
    --tar-file-name "$tarball"
    --use-cvmfs-dropbox
    file://"$wrapper"
    "$config" "$playlist" "$sample" "$count_per_job" "$selection_tag" "$out_tag"
  )

  if [[ -n "$start_override" ]]; then
    cmd+=("$start_override")
  fi

  cmd+=("${extra_args[@]}")


  echo "--------------------------------------------------"
  echo "Running command:"
  printf ' %q' "${cmd[@]}"
  echo
  echo "--------------------------------------------------"

  if [[ "${DRY_RUN}" -eq 0 ]]; then
    "${cmd[@]}"
  fi
}

# -----------------------------
# main loop
# -----------------------------

# playlists with both mc and data
for name in "${both_names[@]}"; do
  playlist="le${name}_p6"
  submit_one "$playlist" mc
  submit_one "$playlist" data
done

# playlists with data only
for name in "${data_only_names[@]}"; do
  playlist="le${name}_p6"
  submit_one "$playlist" data
done
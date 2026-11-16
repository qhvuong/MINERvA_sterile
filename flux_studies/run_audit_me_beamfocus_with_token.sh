#!/usr/bin/env bash
# Pass all arguments through to independent_studies/audit_me_beamfocus_bins.py.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
TOKEN_SCRIPT="${TOKEN_SCRIPT:-${SCRIPT_DIR}/../../getToken.sh}"
TOKEN_REFRESH_SECONDS="${TOKEN_REFRESH_SECONDS:-3600}"
PYTHON_BIN="${PYTHON_BIN:-python}"

if [[ $# -eq 0 || "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    "${PYTHON_BIN}" independent_studies/audit_me_beamfocus_bins.py --help
    exit 0
fi
if [[ ! "${TOKEN_REFRESH_SECONDS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: TOKEN_REFRESH_SECONDS must be a positive integer." >&2
    exit 1
fi
if [[ ! -r "${TOKEN_SCRIPT}" ]]; then
    echo "ERROR: Cannot read token script: ${TOKEN_SCRIPT}" >&2
    exit 1
fi

echo "Refreshing token before audit; renewal interval: ${TOKEN_REFRESH_SECONDS} seconds"
# Source once here so the audit inherits any token-file environment settings.
. "${TOKEN_SCRIPT}"

refresh_token_loop() {
    local sleep_pid=""
    trap 'if [[ -n "${sleep_pid}" ]]; then kill "${sleep_pid}" 2>/dev/null || true; wait "${sleep_pid}" 2>/dev/null || true; fi; exit 0' TERM INT
    while true; do
        sleep "${TOKEN_REFRESH_SECONDS}" &
        sleep_pid=$!
        wait "${sleep_pid}"
        sleep_pid=""
        echo "Refreshing token during audit: $(date -Is)"
        # Renew the same token file the running audit uses. Environment-only
        # changes here cannot change the environment of an existing process.
        if ! . "${TOKEN_SCRIPT}"; then
            echo "WARNING: token renewal failed; retrying at next interval." >&2
        fi
    done
}

TOKEN_PID=""
AUDIT_PID=""
cleanup() {
    if [[ -n "${AUDIT_PID}" ]]; then
        kill "${AUDIT_PID}" 2>/dev/null || true
        wait "${AUDIT_PID}" 2>/dev/null || true
    fi
    if [[ -n "${TOKEN_PID}" ]]; then
        kill "${TOKEN_PID}" 2>/dev/null || true
        wait "${TOKEN_PID}" 2>/dev/null || true
    fi
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

refresh_token_loop &
TOKEN_PID=$!
"${PYTHON_BIN}" independent_studies/audit_me_beamfocus_bins.py "$@" &
AUDIT_PID=$!
status=0
wait "${AUDIT_PID}" || status=$?
AUDIT_PID=""
exit "${status}"

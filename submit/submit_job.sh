#!/usr/bin/env bash
# Submit any worker through the same optional uger: YAML interface.
# bash submit/submit_job.sh --config run.yaml --profile kallisto_reference \
#   [--dry-run] -- [qsub options] submit/submit_kallisto_reference.sh [builder args]
set -euo pipefail
SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${CSM_RESOURCE_CONFIG:-${CSM_CONFIG:-${SUBMIT_DIR}/../config.yaml}}"
PROFILE=""
DRY_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config) CONFIG="$2"; shift 2 ;;
        --profile) PROFILE="$2"; shift 2 ;;
        --dry-run) DRY_ARGS=(--dry-run); shift ;;
        --) shift; break ;;
        *) echo 'Usage: submit_job.sh --config YAML --profile PROFILE [--dry-run] -- [qsub options] WORKER [args]' >&2; exit 2 ;;
    esac
done
: "${PROFILE:?--profile is required}"
# An activated Python/PyYAML environment or CSM_UGER_PYTHON is sufficient.
PYTHON_BIN="${CSM_UGER_PYTHON:-python3}"
TIMEOUT_ARGS=()
if command -v timeout >/dev/null 2>&1; then TIMEOUT_ARGS=(timeout 30s); fi
exec ${TIMEOUT_ARGS[@]+"${TIMEOUT_ARGS[@]}"} "${PYTHON_BIN}" "${SUBMIT_DIR}/../scripts/uger_resources.py" \
    --config "${CONFIG}" --submit-profile "${PROFILE}" ${DRY_ARGS[@]+"${DRY_ARGS[@]}"} -- "$@"

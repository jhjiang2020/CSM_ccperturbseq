#!/bin/bash
# ==============================================================================
# Master Submission Script — ccPerturbSeq Pipeline
# ==============================================================================
# Reads the config to determine the number of libraries, then submits all
# pipeline steps with correct array ranges and dependency chains.
#
# Usage:
#   bash submit/submit_pipeline.sh [--config CONFIG] [--step STEP] [--dry-run]
#
# Examples:
#   bash submit/submit_pipeline.sh                  # Submit all steps
#   bash submit/submit_pipeline.sh --step 1         # Step 1a/1b + optional Step 1c
#   bash submit/submit_pipeline.sh --step 2,3       # Submit Steps 2 and 3
#   bash submit/submit_pipeline.sh --dry-run        # Show what would be submitted
#   bash submit/submit_pipeline.sh --config /path/to/config.yaml
# ==============================================================================

set -euo pipefail

# --- Defaults (derived from script location) ---
SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$(dirname "${SUBMIT_DIR}")/config.yaml"
STEPS="all"
DRY_RUN=false

# --- Parse arguments ---
while [[ $# -gt 0 ]]; do
    case $1 in
        --config)   CONFIG="$2"; shift 2 ;;
        --step)     STEPS="$2"; shift 2 ;;
        --dry-run)  DRY_RUN=true; shift ;;
        -h|--help)
            echo "Usage: bash submit_pipeline.sh [--config CONFIG] [--step STEP] [--dry-run]"
            echo "  --config CONFIG   Path to config.yaml (default: ${CONFIG})"
            echo "  --step STEP       Comma-separated steps to submit: 1,2,3,4,5 or 'all'"
            echo "  --dry-run         Print qsub commands without executing"
            exit 0
            ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [[ "${DRY_RUN}" != true ]] && ! command -v qsub >/dev/null 2>&1; then
    echo "ERROR: qsub is unavailable. On Broad run:" >&2
    echo "       source /broad/software/scripts/useuse; reuse -q UGER || true; command -v qsub" >&2
    exit 1
fi

# --- Bootstrap only the lightweight YAML reader on the login node ---
if [[ ! -r "${CONFIG}" ]]; then
    echo "ERROR: Config is not readable: ${CONFIG}" >&2
    exit 1
fi
_CONDA_PREFIX=$(grep '^ *conda_prefix:' "${CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' || true)
_CONDA_ENV=$(grep '^ *conda_env:' "${CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' || true)
CONDA_ROOT="${_CONDA_PREFIX%/}"
PYTHON_BIN="$(command -v python3 2>/dev/null || true)"
if [[ -z "${PYTHON_BIN}" ]] || ! "${PYTHON_BIN}" -c 'import yaml' >/dev/null 2>&1; then
    if [[ -z "${CONDA_ROOT}" || -z "${_CONDA_ENV}" || ! -x "${CONDA_ROOT}/condabin/conda" ]]; then
        echo "ERROR: Python/PyYAML is unavailable and the configured Conda environment cannot be activated" >&2
        echo "       Check environment.conda_prefix and environment.conda_env in ${CONFIG}" >&2
        exit 1
    fi
    eval "$("${CONDA_ROOT}/condabin/conda" shell.bash hook)"
    conda activate "${_CONDA_ENV}"
    PYTHON_BIN="$(command -v python3 2>/dev/null || true)"
fi
if [[ -z "${PYTHON_BIN}" ]] || ! "${PYTHON_BIN}" -c 'import yaml' >/dev/null 2>&1; then
    echo "ERROR: Python with PyYAML is unavailable after activating ${_CONDA_ENV}" >&2
    exit 1
fi

# Resource resolution/submission reads only small YAML and worker paths. All
# analysis runs in compute jobs. Bound login-node bookkeeping to 30 seconds.
DRY_ARGS=()
if [[ "${DRY_RUN}" == true ]]; then DRY_ARGS+=(--dry-run); fi
TIMEOUT_ARGS=()
if command -v timeout >/dev/null 2>&1; then TIMEOUT_ARGS=(timeout 30s); fi
exec ${TIMEOUT_ARGS[@]+"${TIMEOUT_ARGS[@]}"} "${PYTHON_BIN}" "${SUBMIT_DIR}/../scripts/uger_resources.py" \
    --config "${CONFIG}" --pipeline --steps "${STEPS}" ${DRY_ARGS[@]+"${DRY_ARGS[@]}"}

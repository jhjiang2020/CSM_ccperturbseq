#!/usr/bin/env bash
# Shared submission adapter. Reads only YAML/text on the login node; no data.
# The calling launcher supplies CONFIG/QC_CONFIG or CSM_RESOURCE_CONFIG.
# Child audit stages also source this helper; their tools use explicit NSLOTS.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
_CSM_RESOURCE_SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../scripts" && pwd)/uger_resources.py"
csm_qsub() {
    local profile="$1" resource_config python_bin
    local -a timeout_args=() record_args=()
    shift
    resource_config="${CSM_RESOURCE_CONFIG:-${CONFIG:-${QC_CONFIG:-${AUDIT_CONFIG:-${PROJECT_DIR}/config.yaml}}}}"
    python_bin="${CSM_UGER_PYTHON:-${PYTHON_BIN:-python3}}"
    if command -v timeout >/dev/null 2>&1; then timeout_args=(timeout 30s); fi
    if [[ -n "${LOG_DIR:-}" ]]; then record_args=(--record "${LOG_DIR}/uger_submissions.jsonl"); fi
    ${timeout_args[@]+"${timeout_args[@]}"} "${python_bin}" "${_CSM_RESOURCE_SCRIPT}" \
        --config "${resource_config}" --submit-profile "${profile}" ${record_args[@]+"${record_args[@]}"} -- "$@"
}

#!/bin/bash
# Login-safe scheduler bookkeeping: bounded YAML/resource parsing, no scientific I/O.
set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="${CSM_CONFIG:-${PROJECT_DIR}/config.yaml}"
WORK_DIR=""
REFERENCE=""
KALLISTO=""
MAX_CONCURRENT=2
DRY_RUN=false
PAIRS=50000000
FROM_MODE=init
HOLD_JOB=""
while (( $# )); do
    case "$1" in
        --config) CONFIG="$2"; shift 2;;
        --work-dir) WORK_DIR="$2"; shift 2;;
        --reference) REFERENCE="$2"; shift 2;;
        --kallisto) KALLISTO="$2"; shift 2;;
        --max-concurrent) MAX_CONCURRENT="$2"; shift 2;;
        --pairs) PAIRS="$2"; shift 2;;
        --from-mode) FROM_MODE="$2"; shift 2;;
        --hold-job) HOLD_JOB="$2"; shift 2;;
        --dry-run) DRY_RUN=true; shift;;
        *) echo "Usage: $0 --config FILE --work-dir DIR --reference DIR --kallisto FILE [--max-concurrent 2] [--pairs 50000000] [--dry-run]" >&2; exit 2;;
    esac
done
[[ -n "${WORK_DIR}" && -n "${REFERENCE}" && -n "${KALLISTO}" ]]
[[ "${PAIRS}" =~ ^[1-9][0-9]*$ && "${MAX_CONCURRENT}" =~ ^[1-9][0-9]*$ ]]
[[ "${FROM_MODE}" =~ ^(init|prepare|align|assign|run|report|publish)$ ]]
[[ -z "${HOLD_JOB}" || "${HOLD_JOB}" =~ ^[0-9]+$ ]]
for value in "${CONFIG}" "${WORK_DIR}" "${REFERENCE}" "${KALLISTO}" "${PROJECT_DIR}"; do
    [[ "${value}" == /* && "${value}" != *','* && "${value}" != *$'\n'* ]]
done
LOG_DIR="${WORK_DIR}/logs"
mkdir -p "${LOG_DIR}"
PYTHON_BIN="${CSM_UGER_PYTHON:-python3}"
read -r slots memory hours runtime < <(timeout 30s "${PYTHON_BIN}" "${PROJECT_DIR}/scripts/uger_resources.py" \
    --config "${CONFIG}" --step kallisto_benchmark --format tsv)
common="CSM_PROJECT_DIR=${PROJECT_DIR},CSM_CONFIG=${CONFIG},CSM_BENCH_DIR=${WORK_DIR},CSM_BENCH_REFERENCE=${REFERENCE},CSM_BENCH_KALLISTO=${KALLISTO},CSM_BENCH_PAIRS=${PAIRS},CSM_BENCH_MEMORY_PER_SLOT_GB=${memory},CSM_BENCH_RUNTIME_HOURS=${hours}"
previous="${HOLD_JOB}"
started=false
for mode in init prepare align assign run report publish; do
    if [[ "${mode}" == "${FROM_MODE}" ]]; then started=true; fi
    ${started} || continue
    options=(-N "csm_v2_bench_${mode}" -j y -o "${LOG_DIR}" -v "${common},CSM_BENCH_MODE=${mode}")
    worker="${PROJECT_DIR}/submit/run_gex_benchmark_v2.sh"
    if [[ "${mode}" == publish ]]; then
        worker="${PROJECT_DIR}/submit/publish_gex_benchmark_v2_report.sh"
        options+=(-pe smp 1 -binding linear:1 -l h_vmem=8G,h_rt=1:00:00)
    fi
    if [[ -n "${previous}" ]]; then options+=(-hold_jid "${previous}"); fi
    case "${mode}" in
        prepare) options+=(-t 1-4 -tc "${MAX_CONCURRENT}");;
        align|assign) options+=(-t 1-8 -tc "${MAX_CONCURRENT}");;
        run) options+=(-t 1-48 -tc "${MAX_CONCURRENT}");;
    esac
    cmd=(timeout 30s "${PYTHON_BIN}" "${PROJECT_DIR}/scripts/uger_resources.py" --config "${CONFIG}"
         --submit-profile kallisto_benchmark --record "${LOG_DIR}/uger_submissions.jsonl")
    if ${DRY_RUN}; then
        "${cmd[@]}" --dry-run -- "${options[@]}" "${worker}"
        previous="csm_v2_bench_${mode}"
    else
        submitted=$("${cmd[@]}" -- "${options[@]}" "${worker}")
        [[ "${submitted}" =~ ^[0-9]+(\.[0-9,:-]+)?$ ]]
        previous="${submitted%%.*}"
        echo "${mode}: ${previous}"
    fi
done

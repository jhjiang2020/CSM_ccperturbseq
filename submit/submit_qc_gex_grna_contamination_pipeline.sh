#!/bin/bash
# Submit the GEX gRNA contamination audit array and dependent aggregation job.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/uger_common.sh"

SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SUBMIT_DIR}")"
QC_CONFIG="${CSM_GEX_CONTAM_CONFIG:-}"
MODE="all"
MAX_CONCURRENT=2
ALIGNMENT_MAX_CONCURRENT=""
RAW_WORKERS=8
DRY_RUN=false
RUN_AGGREGATE=true
RUN_ALIGNMENT=true

usage() {
    echo "Usage: $0 [--config FILE] [--mode benchmark|bam|guide|all|filter|alignment|aggregate] [--workers N] [--max-concurrent N] [--alignment-max-concurrent N] [--no-alignment] [--no-aggregate] [--dry-run]"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            QC_CONFIG="$2"
            shift 2
            ;;
        --mode)
            MODE="$2"
            shift 2
            ;;
        --max-concurrent)
            MAX_CONCURRENT="$2"
            shift 2
            ;;
        --workers)
            RAW_WORKERS="$2"
            shift 2
            ;;
        --alignment-max-concurrent)
            ALIGNMENT_MAX_CONCURRENT="$2"
            shift 2
            ;;
        --no-alignment)
            RUN_ALIGNMENT=false
            shift
            ;;
        --no-aggregate)
            RUN_AGGREGATE=false
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: Unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done
: "${QC_CONFIG:?Pass --config /absolute/path/to/my_qc.yaml; see examples/qc/README.md}"


if [[ "${DRY_RUN}" != true ]] && ! command -v qsub >/dev/null 2>&1; then
    echo "ERROR: qsub is unavailable. On Broad run:" >&2
    echo "       source /broad/software/scripts/useuse; reuse -q UGER || true; command -v qsub" >&2
    exit 1
fi

case "${MODE}" in
    benchmark|bam|guide|all|filter|alignment|aggregate) ;;
    *)
        echo "ERROR: Invalid mode: ${MODE}" >&2
        exit 2
        ;;
esac

if [[ ! "${MAX_CONCURRENT}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --max-concurrent must be a positive integer" >&2
    exit 2
fi
if [[ ! "${RAW_WORKERS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --workers must be a positive integer" >&2
    exit 2
fi
ALIGNMENT_MAX_CONCURRENT="${ALIGNMENT_MAX_CONCURRENT:-${MAX_CONCURRENT}}"
if [[ ! "${ALIGNMENT_MAX_CONCURRENT}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --alignment-max-concurrent must be a positive integer" >&2
    exit 2
fi

QC_CONFIG="$(cd "$(dirname "${QC_CONFIG}")" && pwd)/$(basename "${QC_CONFIG}")"
_BASE_CONFIG_VALUE=$(grep '^ *base_config:' "${QC_CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' || true)
_BASE_CONFIG_VALUE="${_BASE_CONFIG_VALUE:-config.yaml}"
if [[ "${_BASE_CONFIG_VALUE}" = /* ]]; then
    BASE_CONFIG="${_BASE_CONFIG_VALUE}"
else
    BASE_CONFIG="$(dirname "${QC_CONFIG}")/${_BASE_CONFIG_VALUE}"
fi
if [[ ! -r "${BASE_CONFIG}" ]]; then
    echo "ERROR: Base config is not readable: ${BASE_CONFIG}" >&2
    exit 1
fi
_CONDA_PREFIX=$(grep '^ *conda_prefix:' "${BASE_CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' || true)
_CONDA_ENV=$(grep '^ *conda_env:' "${BASE_CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' || true)
CONDA_ROOT="${_CONDA_PREFIX%/}"
CONDA_ENV_NAME="${CSM_GEX_CONTAM_ENV:-${_CONDA_ENV}}"
if [[ -z "${CONDA_ROOT}" || -z "${CONDA_ENV_NAME}" ]]; then
    echo "ERROR: Configure environment.conda_prefix and environment.conda_env" >&2
    exit 1
fi
PYTHON_BIN="$(command -v python3 2>/dev/null || command -v python 2>/dev/null || true)"
# Like the established master launcher, only YAML-capable Python is needed to
# submit. If the login shell has none, bootstrap the configured environment;
# every compute job activates the environment independently below.
if [[ -z "${PYTHON_BIN}" ]] || ! "${PYTHON_BIN}" -c 'import yaml' >/dev/null 2>&1; then
    if [[ ! -x "${CONDA_ROOT}/condabin/conda" ]]; then
        echo "ERROR: Python/PyYAML is unavailable and Conda is not at ${CONDA_ROOT}/condabin/conda" >&2
        echo "       Check environment.conda_prefix in ${BASE_CONFIG}" >&2
        exit 1
    fi
    eval "$("${CONDA_ROOT}/condabin/conda" shell.bash hook)"
    conda activate "${CONDA_ENV_NAME}"
    PYTHON_BIN="$(command -v python3 2>/dev/null || command -v python 2>/dev/null || true)"
fi
if [[ -z "${PYTHON_BIN}" ]] || ! "${PYTHON_BIN}" -c 'import yaml' >/dev/null 2>&1; then
    echo "ERROR: Python with PyYAML is unavailable after activating ${CONDA_ENV_NAME}" >&2
    exit 1
fi

N_LIBS=$("${PYTHON_BIN}" -c '
import sys, yaml
cfg = yaml.safe_load(open(sys.argv[1])) or {}
print(len(cfg.get("libraries", [])))
' "${BASE_CONFIG}")
LOG_DIR=$("${PYTHON_BIN}" -c '
import os, sys, yaml
cfg = yaml.safe_load(open(sys.argv[1])) or {}
print(cfg.get("paths", {}).get("log_dir", os.path.join(os.path.dirname(sys.argv[1]), "logs")))
' "${BASE_CONFIG}")
LOG_DIR="${LOG_DIR%/}"
HPC_PROJECT_DIR=$("${PYTHON_BIN}" -c '
import os, sys, yaml
cfg = yaml.safe_load(open(sys.argv[1])) or {}
print(cfg.get("paths", {}).get("project_dir", os.path.dirname(sys.argv[1])))
' "${BASE_CONFIG}")
read -r ALIGNMENT_CONFIG_ENABLED N_ALIGNMENT_VARIANTS ALIGNMENT_THREADS ALIGNMENT_MEMORY ALIGNMENT_RUNTIME < <("${PYTHON_BIN}" -c '
import sys, yaml
cfg = yaml.safe_load(open(sys.argv[1])) or {}
qc = cfg.get("gex_grna_contamination_qc", {})
align = qc.get("alignment_benchmark", {})
read_length = int(qc.get("short_read_length", 54))
expected_short = f"short{read_length}"
variants = [str(value) for value in align.get("variants", ["full", expected_short])]
variants = [expected_short if value == "short54" else value for value in variants]
allowed = {"full", expected_short}
if not variants:
    raise ValueError("alignment_benchmark.variants must not be empty")
if len(set(variants)) != len(variants) or not set(variants).issubset(allowed):
    raise ValueError(
        f"alignment_benchmark.variants must be a unique subset of {sorted(allowed)}"
    )
print(
    str(bool(align.get("enabled", True))).lower(),
    len(variants),
    int(align.get("star_threads", 8)),
    str(align.get("star_memory", "16G")),
    str(align.get("star_runtime", "12:00:00")),
)
' "${QC_CONFIG}")

# Reject invalid optional resources before dispatching any part of the graph.
RESOURCE_CONFIG="${CSM_RESOURCE_CONFIG:-${QC_CONFIG}}"
"${PYTHON_BIN}" "${_CSM_RESOURCE_SCRIPT}" --config "${RESOURCE_CONFIG}" --validate >/dev/null
show_resource_command() {
    local profile="$1"; shift
    "${PYTHON_BIN}" "${_CSM_RESOURCE_SCRIPT}" --config "${RESOURCE_CONFIG}" \
        --submit-profile "${profile}" --dry-run -- "$@"
}

# Render paths as their shared UGER locations. This also makes a local --dry-run
# useful without leaking workstation-specific mount paths into qsub commands.
if [[ "${QC_CONFIG}" == "${PROJECT_DIR}/"* ]]; then
    QC_RELATIVE="${QC_CONFIG#${PROJECT_DIR}/}"
    JOB_QC_CONFIG="${HPC_PROJECT_DIR}/${QC_RELATIVE}"
else
    JOB_QC_CONFIG="${QC_CONFIG}"
fi
JOB_SUBMIT_DIR="${HPC_PROJECT_DIR}/submit"

if [[ "${N_LIBS}" -lt 1 ]]; then
    echo "ERROR: No libraries configured in ${BASE_CONFIG}" >&2
    exit 1
fi

SUBMIT_RAW=true
case "${MODE}" in
    alignment|aggregate) SUBMIT_RAW=false ;;
esac

SUBMIT_ALIGNMENT=false
if [[ "${RUN_ALIGNMENT}" == true && "${ALIGNMENT_CONFIG_ENABLED}" == true ]]; then
    case "${MODE}" in
        benchmark|all|alignment) SUBMIT_ALIGNMENT=true ;;
    esac
fi
if [[ "${MODE}" == "alignment" && "${SUBMIT_ALIGNMENT}" != true ]]; then
    echo "ERROR: --mode alignment requires alignment_benchmark.enabled: true" >&2
    exit 1
fi
if [[ "${MODE}" == "aggregate" ]]; then
    RUN_AGGREGATE=true
fi
N_ALIGNMENT_TASKS=$((N_LIBS * N_ALIGNMENT_VARIANTS))

ARRAY_CMD=(
    qsub -terse
    -t "1-${N_LIBS}"
    -tc "${MAX_CONCURRENT}"
    -pe smp "${RAW_WORKERS}"
    -binding "linear:${RAW_WORKERS}"
    -v "CSM_GEX_CONTAM_CONFIG=${JOB_QC_CONFIG},CSM_GEX_CONTAM_MODE=${MODE},CSM_PROJECT_DIR=${HPC_PROJECT_DIR}"
    -o "${LOG_DIR}/gex_grna_qc.\$JOB_ID.\$TASK_ID.log"
    "${JOB_SUBMIT_DIR}/submit_qc_gex_grna_contamination.sh"
)

echo "GEX gRNA contamination QC"
echo "  Config:         ${JOB_QC_CONFIG}"
echo "  Mode:           ${MODE}"
echo "  Libraries:      ${N_LIBS}"
echo "  Raw workers:    ${RAW_WORKERS} per library (fallback; uger overrides apply)"
echo "  Max concurrent: ${MAX_CONCURRENT} (fallback; see resolved requests)"
echo "  Logs:           ${LOG_DIR}"
if [[ "${SUBMIT_ALIGNMENT}" == true ]]; then
    echo "  Controlled STAR: ${N_ALIGNMENT_TASKS} tasks, max ${ALIGNMENT_MAX_CONCURRENT} concurrent"
else
    echo "  Controlled STAR: skipped (disabled or not applicable to mode ${MODE})"
fi

if [[ "${DRY_RUN}" == true ]]; then
    FINAL_DRY_HOLD=""
    if [[ "${SUBMIT_RAW}" == true ]]; then
        echo '  Array command:'
        show_resource_command gex_contamination "${ARRAY_CMD[@]:1}"
        FINAL_DRY_HOLD="<qc_array_job_id>"
    fi
    if [[ "${SUBMIT_ALIGNMENT}" == true ]]; then
        ALIGNMENT_CMD=(
            qsub -terse
            -t "1-${N_ALIGNMENT_TASKS}"
            -tc "${ALIGNMENT_MAX_CONCURRENT}"
            -pe smp "${ALIGNMENT_THREADS}"
            -binding "linear:${ALIGNMENT_THREADS}"
            -l "h_vmem=${ALIGNMENT_MEMORY}"
            -l "h_rt=${ALIGNMENT_RUNTIME}"
            -v "CSM_GEX_CONTAM_CONFIG=${JOB_QC_CONFIG},CSM_PROJECT_DIR=${HPC_PROJECT_DIR}"
            -o "${LOG_DIR}/gex_grna_star.\$JOB_ID.\$TASK_ID.log"
            "${JOB_SUBMIT_DIR}/submit_qc_gex_grna_alignment_array.sh"
        )
        if [[ -n "${FINAL_DRY_HOLD}" ]]; then
            ALIGNMENT_CMD=(qsub -terse -hold_jid "${FINAL_DRY_HOLD}" "${ALIGNMENT_CMD[@]:2}")
        fi
        echo '  Alignment command:'
        show_resource_command gex_contamination_alignment "${ALIGNMENT_CMD[@]:1}"
        FINAL_DRY_HOLD="<alignment_array_job_id>"
    fi
    if [[ "${RUN_AGGREGATE}" == true ]]; then
        echo '  Aggregate command:'
        AGG_CMD=(-terse
            -v "CSM_GEX_CONTAM_CONFIG=${JOB_QC_CONFIG},CSM_PROJECT_DIR=${HPC_PROJECT_DIR}"
            -o "${LOG_DIR}/gex_grna_qc_aggregate.\$JOB_ID.log")
        if [[ -n "${FINAL_DRY_HOLD}" ]]; then AGG_CMD+=(-hold_jid "${FINAL_DRY_HOLD}"); fi
        show_resource_command gex_contamination_aggregate "${AGG_CMD[@]}" \
            "${JOB_SUBMIT_DIR}/submit_qc_gex_grna_contamination_aggregate.sh"
    fi
    exit 0
fi

mkdir -p "${LOG_DIR}"
FINAL_HOLD_ID=""
if [[ "${SUBMIT_RAW}" == true ]]; then
    ARRAY_JOB_ID=$(csm_qsub gex_contamination "${ARRAY_CMD[@]:1}")
    echo "  Submitted array job: ${ARRAY_JOB_ID}"
    FINAL_HOLD_ID="${ARRAY_JOB_ID%%.*}"
fi

if [[ "${SUBMIT_ALIGNMENT}" == true ]]; then
    ALIGNMENT_CMD=(qsub -terse \
        -t "1-${N_ALIGNMENT_TASKS}" \
        -tc "${ALIGNMENT_MAX_CONCURRENT}" \
        -pe smp "${ALIGNMENT_THREADS}" \
        -binding "linear:${ALIGNMENT_THREADS}" \
        -l "h_vmem=${ALIGNMENT_MEMORY}" \
        -l "h_rt=${ALIGNMENT_RUNTIME}" \
        -v "CSM_GEX_CONTAM_CONFIG=${JOB_QC_CONFIG},CSM_PROJECT_DIR=${HPC_PROJECT_DIR}" \
        -o "${LOG_DIR}/gex_grna_star.\$JOB_ID.\$TASK_ID.log" \
        "${JOB_SUBMIT_DIR}/submit_qc_gex_grna_alignment_array.sh")
    if [[ -n "${FINAL_HOLD_ID}" ]]; then
        ALIGNMENT_CMD=(qsub -terse -hold_jid "${FINAL_HOLD_ID}" "${ALIGNMENT_CMD[@]:2}")
    fi
    ALIGNMENT_JOB_ID=$(csm_qsub gex_contamination_alignment "${ALIGNMENT_CMD[@]:1}")
    echo "  Submitted controlled STAR array: ${ALIGNMENT_JOB_ID}"
    FINAL_HOLD_ID="${ALIGNMENT_JOB_ID%%.*}"
fi

if [[ "${RUN_AGGREGATE}" == true ]]; then
    AGG_CMD=(qsub -terse \
        -v "CSM_GEX_CONTAM_CONFIG=${JOB_QC_CONFIG},CSM_PROJECT_DIR=${HPC_PROJECT_DIR}" \
        -o "${LOG_DIR}/gex_grna_qc_aggregate.\$JOB_ID.log" \
        "${JOB_SUBMIT_DIR}/submit_qc_gex_grna_contamination_aggregate.sh")
    if [[ -n "${FINAL_HOLD_ID}" ]]; then
        AGG_CMD=(qsub -terse -hold_jid "${FINAL_HOLD_ID}" "${AGG_CMD[@]:2}")
    fi
    AGG_JOB_ID=$(csm_qsub gex_contamination_aggregate "${AGG_CMD[@]:1}")
    echo "  Submitted aggregate job: ${AGG_JOB_ID}"
fi

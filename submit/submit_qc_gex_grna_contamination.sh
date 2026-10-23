#!/bin/bash
# ==============================================================================
# GEX gRNA cross-contamination QC (UGE array; one task per GEX sublibrary)
# ==============================================================================
# qsub -t 1-N \
#   -v CSM_GEX_CONTAM_CONFIG=/abs/path/qc_config_gex_grna_contamination.yaml \
#   submit/submit_qc_gex_grna_contamination.sh
#
# Optional: CSM_GEX_CONTAM_MODE=benchmark|bam|all|filter (default: all)
# ==============================================================================
#$ -q broad
#$ -N csm_gex_grna_qc
#$ -j y
#$ -l h_vmem=8G
#$ -l h_rt=48:00:00
#$ -pe smp 8
#$ -binding linear:8

set -euo pipefail
# Explicit tools use NSLOTS; nested numerical pools use one thread.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SUBMIT_DIR}")"
: "${CSM_GEX_CONTAM_CONFIG:?Set an absolute QC YAML path; see examples/qc/README.md}"
QC_CONFIG="${CSM_GEX_CONTAM_CONFIG}"
MODE="${CSM_GEX_CONTAM_MODE:-all}"

# Resolve the base config without Python: compute nodes can start with a
# minimal PATH, and Python is supplied by the configured Conda environment.
_BASE_CONFIG_VALUE=$(grep '^ *base_config:' "${QC_CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' || true)
_BASE_CONFIG_VALUE="${_BASE_CONFIG_VALUE:-config.yaml}"
if [[ "${_BASE_CONFIG_VALUE}" = /* ]]; then
    BASE_CONFIG="${_BASE_CONFIG_VALUE}"
else
    BASE_CONFIG="$(cd "$(dirname "${QC_CONFIG}")" && pwd)/${_BASE_CONFIG_VALUE}"
fi
if [[ ! -r "${BASE_CONFIG}" ]]; then
    echo "ERROR: Base config is not readable: ${BASE_CONFIG}" >&2
    exit 1
fi

# UGER runs a spooled copy of this file, so BASH_SOURCE is not a reliable
# project root inside a job. Match the established Step 2-5 wrapper contract:
# qsub environment > base config > directory containing the base config.
_CONFIG_PROJECT_DIR=$(grep '^ *project_dir:' "${BASE_CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' || true)
PROJECT_DIR="${CSM_PROJECT_DIR:-${_CONFIG_PROJECT_DIR:-$(dirname "${BASE_CONFIG}")}}"

_CONDA_PREFIX=$(grep '^ *conda_prefix:' "${BASE_CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' || true)
_CONDA_ENV=$(grep '^ *conda_env:' "${BASE_CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' || true)
CONDA_ROOT="${_CONDA_PREFIX%/}"
CONDA_ENV_NAME="${CSM_GEX_CONTAM_ENV:-${_CONDA_ENV}}"
if [[ -z "${CONDA_ROOT}" || -z "${CONDA_ENV_NAME}" ]]; then
    echo "ERROR: Configure environment.conda_prefix and environment.conda_env" >&2
    exit 1
fi
if [[ ! -x "${CONDA_ROOT}/condabin/conda" ]]; then
    echo "ERROR: Conda executable is not available at ${CONDA_ROOT}/condabin/conda" >&2
    echo "       Check environment.conda_prefix in ${BASE_CONFIG}" >&2
    exit 1
fi
eval "$("${CONDA_ROOT}/condabin/conda" shell.bash hook)"
conda activate "${CONDA_ENV_NAME}"
PYTHON_BIN="$(command -v python3 2>/dev/null || command -v python 2>/dev/null || true)"
if [[ -z "${PYTHON_BIN}" ]]; then
    echo "ERROR: Neither python3 nor python is available after activating ${CONDA_ENV_NAME}" >&2
    exit 1
fi
REQUIRED_IMPORTS='import edlib, pandas, yaml'
if [[ "${MODE}" == "all" || "${MODE}" == "bam" ]]; then
    REQUIRED_IMPORTS='import edlib, pandas, pysam, yaml'
fi
if ! "${PYTHON_BIN}" -c "${REQUIRED_IMPORTS}" >/dev/null 2>&1; then
    echo "ERROR: ${CONDA_ENV_NAME} is incomplete for QC mode ${MODE}" >&2
    echo "       Install requirements_extra.txt into the production splitpipe environment" >&2
    exit 1
fi

TASK_ID="${SGE_TASK_ID:-1}"
LIB_NAME=$("${PYTHON_BIN}" -c "
import yaml
c=yaml.safe_load(open('${BASE_CONFIG}'))
libs=c['libraries']
i=int('${TASK_ID}')-1
if i < 0 or i >= len(libs): raise IndexError(f'Task {i+1} outside 1-{len(libs)}')
print(libs[i]['sublibrary_name'])
")

echo "=========================================="
echo " GEX gRNA contamination QC"
echo " Host:       $(hostname)"
echo " Date:       $(date)"
echo " Job:        ${JOB_ID:-local}.${TASK_ID}"
echo " Library:    ${LIB_NAME}"
echo " Mode:       ${MODE}"
echo " Config:     ${QC_CONFIG}"
echo " Environment:${CONDA_ENV_NAME}"
echo " Python:     ${PYTHON_BIN}"
echo " Threads:    ${NSLOTS:-8}"
echo "=========================================="

"${PYTHON_BIN}" "${PROJECT_DIR}/scripts/qc_gex_grna_contamination.py" \
    --config "${QC_CONFIG}" \
    --library "${LIB_NAME}" \
    --mode "${MODE}" \
    --threads "${NSLOTS:-8}"

echo "GEX gRNA contamination QC complete: $(date)"

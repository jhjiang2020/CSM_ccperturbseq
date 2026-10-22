#!/bin/bash
# Aggregate completed per-library GEX gRNA contamination QC outputs.
#$ -q broad
#$ -N csm_gex_grna_qc_agg
#$ -j y
#$ -l h_vmem=4G
#$ -l h_rt=01:00:00
#$ -pe smp 1

set -euo pipefail
# Explicit tools use NSLOTS; nested numerical pools use one thread.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SUBMIT_DIR}")"
QC_CONFIG="${CSM_GEX_CONTAM_CONFIG:-${PROJECT_DIR}/qc_config_gex_grna_contamination.yaml}"

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
if ! "${PYTHON_BIN}" -c 'import matplotlib, pandas, yaml' >/dev/null 2>&1; then
    echo "ERROR: ${CONDA_ENV_NAME} is incomplete for QC aggregation" >&2
    echo "       Install requirements_extra.txt into the production splitpipe environment" >&2
    exit 1
fi

"${PYTHON_BIN}" "${PROJECT_DIR}/scripts/qc_gex_grna_contamination.py" \
    --config "${QC_CONFIG}" --mode aggregate

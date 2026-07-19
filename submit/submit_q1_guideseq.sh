#!/bin/bash
# ==============================================================================
# Q1: Guide-Seq QC — Decoupling & Pre-Transduction Recombination
# ==============================================================================
# Submit via:
#   qsub -v CSM_QC_CONFIG=/abs/path/to/qc_config_guideseq.yaml -o /path/to/logs/ submit/submit_q1_guideseq.sh
# ==============================================================================
# UGE Directives:
#$ -q broad
#$ -N csm_q1_guideseq
#$ -j y
#$ -l h_vmem=4G
#$ -l h_rt=12:00:00
#$ -pe smp 8
#$ -binding linear:8

set -euo pipefail

# --- Configuration (must be passed as absolute path via -v CSM_QC_CONFIG=...) ---
if [[ -z "${CSM_QC_CONFIG:-}" ]]; then
    echo "ERROR: CSM_QC_CONFIG not set. Submit with:"
    echo "  qsub -v CSM_QC_CONFIG=/abs/path/to/qc_config_guideseq.yaml submit/submit_q1_guideseq.sh"
    exit 1
fi
QC_CONFIG="$(cd "$(dirname "${CSM_QC_CONFIG}")" && pwd)/$(basename "${CSM_QC_CONFIG}")"
_SC_CONFIG_VALUE=$(grep '^ *sc_config:' "${QC_CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' || true)
if [[ -n "${_SC_CONFIG_VALUE}" && "${_SC_CONFIG_VALUE}" != /* ]]; then
    BASE_CONFIG="$(dirname "${QC_CONFIG}")/${_SC_CONFIG_VALUE}"
else
    BASE_CONFIG="${_SC_CONFIG_VALUE}"
fi
_CONFIG_PROJECT_DIR=""
if [[ -n "${BASE_CONFIG}" && -r "${BASE_CONFIG}" ]]; then
    _CONFIG_PROJECT_DIR=$(grep '^ *project_dir:' "${BASE_CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' || true)
fi
PROJECT_DIR="${CSM_PROJECT_DIR:-${_CONFIG_PROJECT_DIR:-$(dirname "${QC_CONFIG}")}}"
SCRIPTS_DIR="${PROJECT_DIR}/scripts"

# Initialize the configured QC environment; no user-specific path/name fallback.
_CONDA_PREFIX=$(grep '^ *qc_conda_prefix:' "${QC_CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' || true)
if [[ -z "${_CONDA_PREFIX}" && -n "${BASE_CONFIG}" && -r "${BASE_CONFIG}" ]]; then
    _CONDA_PREFIX=$(grep '^ *conda_prefix:' "${BASE_CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' || true)
fi
_QC_CONDA_ENV=$(grep '^ *qc_conda_env:' "${QC_CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' || true)
CONDA_ROOT="${_CONDA_PREFIX%/}"
CONDA_ENV_NAME="${CSM_QC_ENV:-${_QC_CONDA_ENV}}"
if [[ -z "${CONDA_ROOT}" || -z "${CONDA_ENV_NAME}" || ! -x "${CONDA_ROOT}/condabin/conda" ]]; then
    echo "ERROR: Configure qc_conda_prefix (or sc_config environment.conda_prefix) and qc_conda_env" >&2
    exit 1
fi
eval "$("${CONDA_ROOT}/condabin/conda" shell.bash hook)"
conda activate "${CONDA_ENV_NAME}"
PYTHON_BIN=$(command -v python3 2>/dev/null || command -v python 2>/dev/null || true)
if [[ -z "${PYTHON_BIN}" ]] || ! "${PYTHON_BIN}" -c 'import edlib, matplotlib, numpy, pandas, tqdm, yaml' >/dev/null 2>&1; then
    echo "ERROR: ${CONDA_ENV_NAME} is incomplete for Q1" >&2
    echo "       Create/update it with: conda env update -n ${CONDA_ENV_NAME} -f environment_qc.yml --prune" >&2
    exit 1
fi

# --- Log directory from base config ---
LOG_DIR=$("${PYTHON_BIN}" -c "
import yaml, os
with open('${QC_CONFIG}') as f:
    qc = yaml.safe_load(f)
base = qc.get('base_config', '')
if base and not os.path.isabs(base):
    base = os.path.join('${PROJECT_DIR}', base)
if base and os.path.exists(base):
    with open(base) as f:
        cfg = yaml.safe_load(f)
    log_dir = cfg.get('paths', {}).get('log_dir', '${PROJECT_DIR}/logs')
else:
    log_dir = '${PROJECT_DIR}/logs'
print(log_dir)
")
mkdir -p "${LOG_DIR}"

echo "=========================================="
echo "  Q1: Guide-Seq QC"
echo "  Host:       $(hostname)"
echo "  Date:       $(date)"
echo "  Job:        ${JOB_ID:-local}"
echo "  Cores:      ${NSLOTS:-8}"
echo "  QC Config:  ${QC_CONFIG}"
echo "  Project:    ${PROJECT_DIR}"
echo "  Logs:       ${LOG_DIR}"
echo "=========================================="

# --- Run ---
"${PYTHON_BIN}" "${SCRIPTS_DIR}/qc_q1_guideseq.py" --config "${QC_CONFIG}" --n-workers "${NSLOTS:-8}"

echo ""
echo "Q1 Guide-Seq QC complete: $(date)"

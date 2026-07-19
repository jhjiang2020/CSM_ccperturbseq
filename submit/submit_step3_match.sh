#!/bin/bash
# ==============================================================================
# Step 3: Match gRNA + iBAR (ARRAY JOB — one task per sublibrary)
# ==============================================================================
# Submit via submit_pipeline.sh (recommended), or manually:
#   qsub -hold_jid csm_step2_bc -t 1-N submit/submit_step3_match.sh
# ==============================================================================
# UGE Directives:
#$ -q broad
#$ -N csm_step3_match
#$ -j y
#$ -l h_vmem=4G
#$ -l h_rt=12:00:00
#$ -pe smp 8
#$ -binding linear:8
# NOTE: -t (array range) and -o (log dir) are set dynamically by submit_pipeline.sh

set -euo pipefail

# --- Configuration ---
if [[ -n "${CSM_CONFIG:-}" ]]; then
    CONFIG="${CSM_CONFIG}"
else
    SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    CONFIG="$(dirname "${SUBMIT_DIR}")/config.yaml"
fi

# --- Resolve project dir: env var > config (grep) > dirname(config) ---
if [[ -n "${CSM_PROJECT_DIR:-}" ]]; then
    PROJECT_DIR="${CSM_PROJECT_DIR}"
else
    PROJECT_DIR=$(grep '^ *project_dir:' "${CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/')
fi
SCRIPTS_DIR="${PROJECT_DIR}/scripts"

# --- Initialize conda (parsed from config via grep — no python3 needed yet) ---
_CONDA_PREFIX=$(grep '^ *conda_prefix:' "${CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/')
_CONDA_ENV=$(grep '^ *conda_env:' "${CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/')
[[ -n "${_CONDA_PREFIX}" && -n "${_CONDA_ENV}" ]] || { echo "ERROR: Configure environment.conda_prefix and environment.conda_env" >&2; exit 1; }
CONDA_PREFIX="${_CONDA_PREFIX}"
CONDA_ENV="${_CONDA_ENV}"
[[ -x "${CONDA_PREFIX}/condabin/conda" ]] || { echo "ERROR: Missing configured Conda: ${CONDA_PREFIX}/condabin/conda" >&2; exit 1; }
eval "$(${CONDA_PREFIX}/condabin/conda shell.bash hook)"
conda activate "${CONDA_ENV}"

LOG_DIR=$(python3 -c "import yaml; print(yaml.safe_load(open('${CONFIG}'))['paths'].get('log_dir', '${PROJECT_DIR}/logs'))")
mkdir -p "${LOG_DIR}"

# --- Map task ID to sublibrary ---
LIB_NAME=$(python3 -c "
import yaml
with open('${CONFIG}') as f:
    cfg = yaml.safe_load(f)
libs = cfg['libraries']
idx = ${SGE_TASK_ID} - 1
if idx >= len(libs):
    raise IndexError(f'Task ID ${SGE_TASK_ID} exceeds number of libraries ({len(libs)})')
print(libs[idx]['sublibrary_name'])
")

echo "=========================================="
echo "  Step 3: Match gRNA + iBAR (task ${SGE_TASK_ID})"
echo "  Host:       $(hostname)"
echo "  Date:       $(date)"
echo "  Job:        ${JOB_ID}.${SGE_TASK_ID}"
echo "  Sublibrary: ${LIB_NAME}"
echo "  Cores:      ${NSLOTS:-8}"
echo "=========================================="

# --- Run ---
python "${SCRIPTS_DIR}/03_match_grna_ibar.py" \
    --config "${CONFIG}" \
    --library "${LIB_NAME}" \
    --n_processes "${NSLOTS:-8}"

echo ""
echo "Step 3 complete for ${LIB_NAME}: $(date)"

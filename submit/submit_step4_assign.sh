#!/bin/bash
# ==============================================================================
# Step 4: Assign gRNAs to Cells (all libraries — run AFTER step 3 completes)
# ==============================================================================
# Submit via submit_pipeline.sh (recommended), or manually:
#   qsub -hold_jid csm_step3_match submit/submit_step4_assign.sh
# ==============================================================================
# UGE Directives:
#$ -q broad
#$ -N csm_step4_assign
#$ -j y
#$ -l h_vmem=32G
#$ -l h_rt=4:00:00
#$ -pe smp 4

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
CONDA_PREFIX="${_CONDA_PREFIX:-$HOME/my-conda}"
CONDA_ENV="${_CONDA_ENV:-splitpipe}"
eval "$(${CONDA_PREFIX}/condabin/conda shell.bash hook)"
conda activate "${CONDA_ENV}"

LOG_DIR=$(python3 -c "import yaml; print(yaml.safe_load(open('${CONFIG}'))['paths'].get('log_dir', '${PROJECT_DIR}/logs'))")
mkdir -p "${LOG_DIR}"

echo "=========================================="
echo "  Step 4: Assign gRNAs to Cells"
echo "  Host: $(hostname)"
echo "  Date: $(date)"
echo "  Job:  ${JOB_ID:-local}"
echo "=========================================="

# --- Run ---
python "${SCRIPTS_DIR}/04_assign_grna_to_cells.py" --config "${CONFIG}" --ncpus ${NSLOTS:-4}

echo ""
echo "Step 4 complete: $(date)"

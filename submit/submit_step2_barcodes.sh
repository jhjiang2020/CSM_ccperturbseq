#!/bin/bash
# ==============================================================================
# Step 2: Extract gRNA Barcodes (ARRAY JOB — one task per sublibrary)
# ==============================================================================
# Submit via submit_pipeline.sh (recommended), or manually:
#   qsub -t 1-N submit/submit_step2_barcodes.sh
# ==============================================================================
# UGE Directives:
#$ -q broad
#$ -N csm_step2_bc
#$ -j y
#$ -l h_vmem=4G
#$ -l h_rt=12:00:00
#$ -pe smp 1
# NOTE: -t (array range) and -o (log dir) are set dynamically by submit_pipeline.sh

set -euo pipefail

# --- Configuration ---
if [[ -n "${CSM_CONFIG:-}" ]]; then
    CONFIG="${CSM_CONFIG}"
else
    SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    CONFIG="$(dirname "${SUBMIT_DIR}")/config.yaml"
fi
SCRIPTS_DIR="$(dirname "${CONFIG}")/scripts"

# Initialize conda (must happen before any python3 calls)
eval "$($HOME/my-conda/condabin/conda shell.bash hook)"
conda activate splitpipe

LOG_DIR=$(python3 -c "import yaml; print(yaml.safe_load(open('${CONFIG}'))['paths'].get('log_dir', '$(dirname "${CONFIG}")/logs'))")
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
echo "  Step 2: Extract gRNA Barcodes (task ${SGE_TASK_ID})"
echo "  Host:       $(hostname)"
echo "  Date:       $(date)"
echo "  Job:        ${JOB_ID}.${SGE_TASK_ID}"
echo "  Sublibrary: ${LIB_NAME}"
echo "=========================================="

# --- Run ---
python "${SCRIPTS_DIR}/02_extract_grna_barcodes.py" \
    --config "${CONFIG}" \
    --library "${LIB_NAME}"

echo ""
echo "Step 2 complete for ${LIB_NAME}: $(date)"

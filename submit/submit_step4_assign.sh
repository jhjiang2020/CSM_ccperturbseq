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
SCRIPTS_DIR="$(dirname "${CONFIG}")/scripts"

# Initialize conda (must happen before any python3 calls)
eval "$($HOME/my-conda/condabin/conda shell.bash hook)"
conda activate splitpipe

LOG_DIR=$(python3 -c "import yaml; print(yaml.safe_load(open('${CONFIG}'))['paths'].get('log_dir', '$(dirname "${CONFIG}")/logs'))")
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

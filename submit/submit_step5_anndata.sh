#!/bin/bash
# ==============================================================================
# Step 5: Build Combined AnnData (run AFTER steps 1b and 4 complete)
# ==============================================================================
# Submit via submit_pipeline.sh (recommended), or manually:
#   qsub -hold_jid csm_step1b_combine,csm_step4_assign submit/submit_step5_anndata.sh
# ==============================================================================
# UGE Directives:
#$ -q broad
#$ -N csm_step5_anndata
#$ -j y
#$ -l h_vmem=32G
#$ -l h_rt=4:00:00
#$ -pe smp 1

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
echo "  Step 5: Build Combined AnnData"
echo "  Host: $(hostname)"
echo "  Date: $(date)"
echo "  Job:  ${JOB_ID:-local}"
echo "=========================================="

# --- Run ---
python "${SCRIPTS_DIR}/05_build_anndata.py" --config "${CONFIG}"

echo ""
echo "Step 5 complete: $(date)"

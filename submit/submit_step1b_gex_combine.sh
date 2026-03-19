#!/bin/bash
# ==============================================================================
# Step 1b: Combine GEX sublibraries (run AFTER step1a array completes)
# ==============================================================================
# Submit via submit_pipeline.sh (recommended), or manually:
#   qsub -hold_jid csm_step1a_gex submit/submit_step1b_gex_combine.sh
# ==============================================================================
# UGE Directives:
#$ -q broad
#$ -N csm_step1b_combine
#$ -j y
#$ -l h_vmem=16G
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
echo "  Step 1b: Combine GEX Sublibraries"
echo "  Host: $(hostname)"
echo "  Date: $(date)"
echo "  Job:  ${JOB_ID:-local}"
echo "=========================================="

# Read config
SPLITPIPE_DIR=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); print(cfg['paths']['splitpipe_dir'])")
OUTPUT_DIR=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); print(cfg['paths']['output_dir'])")
SPLITPIPE="${SPLITPIPE_DIR}/split-pipe"

# Build list of sublibrary output directories
SUBLIB_DIRS=$(python3 -c "
import yaml, os
with open('${CONFIG}') as f:
    cfg = yaml.safe_load(f)
out = cfg['paths']['output_dir']
dirs = [os.path.join(out, f\"gex_{lib['sublibrary_name']}\") for lib in cfg['libraries']]
print(' '.join(dirs))
")

COMB_OUT="${OUTPUT_DIR}/gex_combined"

echo "Sublibraries: ${SUBLIB_DIRS}"
echo "Output:       ${COMB_OUT}"
echo ""

# Verify all sublibrary outputs exist
for d in ${SUBLIB_DIRS}; do
    if [[ ! -d "$d" ]]; then
        echo "ERROR: Sublibrary output not found: $d"
        echo "Ensure Step 1a array job completed successfully."
        exit 1
    fi
done

# Execute the custom python wrapper that patches split-pipe's list-index bug at runtime
python3 "${SCRIPTS_DIR}/splitpipe_wrapper.py" \
    --mode comb \
    --sublibraries ${SUBLIB_DIRS} \
    --output_dir "${COMB_OUT}" \
    --save_anndata

echo ""
echo "Step 1b complete: $(date)"
echo "Combined output: ${COMB_OUT}"

#!/bin/bash
# ==============================================================================
# Step 1a: Run GEX split-pipe per sublibrary (ARRAY JOB)
# ==============================================================================
# Submit via submit_pipeline.sh (recommended), or manually:
#   qsub -t 1-N submit/submit_step1a_gex_array.sh
# ==============================================================================
# UGE Directives:
#$ -q broad
#$ -N csm_step1a_gex
#$ -j y
#$ -l h_vmem=8G
#$ -l h_rt=24:00:00
#$ -pe smp 16
# NOTE: -t (array range) and -o (log dir) are set dynamically by submit_pipeline.sh

set -euo pipefail

# --- Configuration ---
# CSM_CONFIG is passed via qsub -v from submit_pipeline.sh.
# Falls back to BASH_SOURCE derivation for local/interactive runs.
if [[ -n "${CSM_CONFIG:-}" ]]; then
    CONFIG="${CSM_CONFIG}"
else
    SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    CONFIG="$(dirname "${SUBMIT_DIR}")/config.yaml"
fi

# Initialize conda (must happen before any python3 calls)
eval "$($HOME/my-conda/condabin/conda shell.bash hook)"
conda activate splitpipe

LOG_DIR=$(python3 -c "import yaml; print(yaml.safe_load(open('${CONFIG}'))['paths'].get('log_dir', '$(dirname "${CONFIG}")/logs'))")
mkdir -p "${LOG_DIR}"

# --- Map array task ID to sublibrary ---
# Read library info from config
LIB_NAME=$(python3 -c "
import yaml
with open('${CONFIG}') as f:
    cfg = yaml.safe_load(f)
libs = cfg['libraries']
idx = ${SGE_TASK_ID} - 1
if idx >= len(libs):
    raise IndexError(f'Task ID {${SGE_TASK_ID}} exceeds number of libraries ({len(libs)})')
print(libs[idx]['sublibrary_name'])
")

GEX_PREFIX=$(python3 -c "
import yaml
with open('${CONFIG}') as f:
    cfg = yaml.safe_load(f)
print(cfg['libraries'][${SGE_TASK_ID}-1]['gex_prefix'])
")

SAMPLE_NAME=$(python3 -c "
import yaml
with open('${CONFIG}') as f:
    cfg = yaml.safe_load(f)
print(cfg['libraries'][${SGE_TASK_ID}-1]['sample_name'])
")

# Read other config values
FASTQ_DIR=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); print(cfg['paths']['fastq_dir'])")
SPLITPIPE_DIR=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); print(cfg['paths']['splitpipe_dir'])")
GENOME_DIR=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); print(cfg['paths']['genome_dir'])")
OUTPUT_DIR=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); print(cfg['paths']['output_dir'])")
KIT=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); print(cfg['splitpipe']['kit'])")
CHEMISTRY=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); print(cfg['splitpipe']['chemistry'])")
SAMPLE_WELLS=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); print(cfg['splitpipe']['sample_wells'])")
BYPASS_KIT=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); print(cfg['splitpipe'].get('bypass_kit_score', False))")

SPLITPIPE="${SPLITPIPE_DIR}/split-pipe"
OUT_SUBLIB="${OUTPUT_DIR}/gex_${LIB_NAME}"

# Build optional flags
EXTRA_FLAGS=""
if [[ "$BYPASS_KIT" == "True" ]]; then
    EXTRA_FLAGS="--kit_score_skip"
    echo "  (kit score check will be skipped)"
fi

echo "=========================================="
echo "  Step 1a: GEX split-pipe (task ${SGE_TASK_ID})"
echo "  Host:       $(hostname)"
echo "  Date:       $(date)"
echo "  Job:        ${JOB_ID}.${SGE_TASK_ID}"
echo "  Sublibrary: ${LIB_NAME}"
echo "  GEX prefix: ${GEX_PREFIX}"
echo "  Sample:     ${SAMPLE_NAME}"
echo "  Cores:      ${NSLOTS:-16}"
echo "=========================================="

# Find FASTQ files
FQ1=$(find "${FASTQ_DIR}" -name "*${GEX_PREFIX}*_R1_*fastq.gz" | head -1)
FQ2=$(find "${FASTQ_DIR}" -name "*${GEX_PREFIX}*_R2_*fastq.gz" | head -1)

if [[ -z "$FQ1" || -z "$FQ2" ]]; then
    echo "ERROR: Could not find FASTQ files for GEX prefix '${GEX_PREFIX}'"
    exit 1
fi

echo "R1: ${FQ1}"
echo "R2: ${FQ2}"
echo "Output: ${OUT_SUBLIB}"
echo ""

# Run split-pipe
"${SPLITPIPE}" \
    --mode all \
    --kit "${KIT}" \
    --chemistry "${CHEMISTRY}" \
    --genome_dir "${GENOME_DIR}" \
    --fq1 "${FQ1}" \
    --fq2 "${FQ2}" \
    --output_dir "${OUT_SUBLIB}" \
    --sample "${SAMPLE_NAME}" "${SAMPLE_WELLS}" \
    --nthreads "${NSLOTS:-16}" \
    --save_anndata \
    ${EXTRA_FLAGS}

echo ""
echo "Step 1a complete for ${LIB_NAME}: $(date)"

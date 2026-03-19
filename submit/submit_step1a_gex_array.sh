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

# --- Initialize conda (parsed from config via grep — no python3 needed yet) ---
_CONDA_PREFIX=$(grep '^ *conda_prefix:' "${CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/')
_CONDA_ENV=$(grep '^ *conda_env:' "${CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/')
CONDA_PREFIX="${_CONDA_PREFIX:-$HOME/my-conda}"
CONDA_ENV="${_CONDA_ENV:-splitpipe}"
eval "$(${CONDA_PREFIX}/condabin/conda shell.bash hook)"
conda activate "${CONDA_ENV}"

LOG_DIR=$(python3 -c "import yaml; print(yaml.safe_load(open('${CONFIG}'))['paths'].get('log_dir', '$(dirname \"${CONFIG}\")/logs'))")
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
FQ1_LIST=$(find "${FASTQ_DIR}" -name "*${GEX_PREFIX}*_R1_*fastq.gz" | sort | tr '\n' ' ')
FQ2_LIST=$(find "${FASTQ_DIR}" -name "*${GEX_PREFIX}*_R2_*fastq.gz" | sort | tr '\n' ' ')

if [[ -z "${FQ1_LIST// /}" || -z "${FQ2_LIST// /}" ]]; then
    echo "ERROR: Could not find FASTQ files for GEX prefix '${GEX_PREFIX}' in ${FASTQ_DIR}"
    exit 1
fi

echo "R1 files: ${FQ1_LIST}"
echo "R2 files: ${FQ2_LIST}"
echo "Output: ${OUT_SUBLIB}"
echo ""

# Run split-pipe
# Note: FQ1_LIST and FQ2_LIST are unquoted so bash splits them into multiple arguments
"${SPLITPIPE}" \
    --mode all \
    --kit "${KIT}" \
    --chemistry "${CHEMISTRY}" \
    --genome_dir "${GENOME_DIR}" \
    --fq1 ${FQ1_LIST} \
    --fq2 ${FQ2_LIST} \
    --output_dir "${OUT_SUBLIB}" \
    --sample "${SAMPLE_NAME}" "${SAMPLE_WELLS}" \
    --nthreads "${NSLOTS:-16}" \
    --save_anndata \
    ${EXTRA_FLAGS}

echo ""
echo "Step 1a complete for ${LIB_NAME}: $(date)"

#!/bin/bash
# ==============================================================================
# Test: v0.3 Stem Ambiguity Impact (ARRAY JOB — one task per sublibrary)
# ==============================================================================
# Submit via:
#   qsub -v CSM_CONFIG=/abs/path/to/config.yaml -t 1-4 submit/submit_test_v03_impact.sh
# Or for a single sample:
#   qsub -v CSM_CONFIG=/abs/path/to/config.yaml -t 1 submit/submit_test_v03_impact.sh
# ==============================================================================
# UGE Directives:
#$ -q broad
#$ -N csm_v03_impact
#$ -j y
#$ -l h_vmem=8G
#$ -l h_rt=6:00:00
#$ -pe smp 1
#$ -o /medpop/rgupta/Jiahao/logs/CSM_ccperturbseq/v03_impact.$JOB_ID.$TASK_ID.log

set -euo pipefail

# --- Configuration (must be passed as absolute path via -v CSM_CONFIG=...) ---
if [[ -z "${CSM_CONFIG:-}" ]]; then
    echo "ERROR: CSM_CONFIG not set. Submit with:"
    echo "  qsub -v CSM_CONFIG=/abs/path/to/config.yaml -t 1-4 submit/submit_test_v03_impact.sh"
    exit 1
fi
CONFIG="${CSM_CONFIG}"
SCRIPTS_DIR="$(dirname "${CONFIG}")/scripts"

# Initialize conda
eval "$($HOME/my-conda/condabin/conda shell.bash hook)"
conda activate splitpipe

# Map task ID to sample name
SAMPLES=(C1 C2 T1 T2)
IDX=$(( ${SGE_TASK_ID:-1} - 1 ))
if [[ $IDX -ge ${#SAMPLES[@]} ]]; then
    echo "ERROR: Task ID ${SGE_TASK_ID} exceeds number of samples (${#SAMPLES[@]})"
    exit 1
fi
SAMPLE="${SAMPLES[$IDX]}"

echo "=========================================="
echo "  v0.3 Stem Impact Test (task ${SGE_TASK_ID:-1})"
echo "  Host:     $(hostname)"
echo "  Date:     $(date)"
echo "  Job:      ${JOB_ID:-local}.${SGE_TASK_ID:-1}"
echo "  Sample:   ${SAMPLE}"
echo "  Config:   ${CONFIG}"
echo "=========================================="

# --- Run (process all reads) ---
python3 "${SCRIPTS_DIR}/test_v03_stem_impact.py" \
    --config "${CONFIG}" \
    --sample "${SAMPLE}"

echo ""
echo "v0.3 impact test complete for ${SAMPLE}: $(date)"

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
#$ -o /medpop/rgupta/Jiahao/logs/CSM_ccperturbseq/q1_guideseq.$JOB_ID.log

set -euo pipefail

# --- Configuration (must be passed as absolute path via -v CSM_QC_CONFIG=...) ---
if [[ -z "${CSM_QC_CONFIG:-}" ]]; then
    echo "ERROR: CSM_QC_CONFIG not set. Submit with:"
    echo "  qsub -v CSM_QC_CONFIG=/abs/path/to/qc_config_guideseq.yaml submit/submit_q1_guideseq.sh"
    exit 1
fi
QC_CONFIG="${CSM_QC_CONFIG}"
PROJECT_DIR="$(dirname "${QC_CONFIG}")"
SCRIPTS_DIR="${PROJECT_DIR}/scripts"

# Initialize conda
eval "$($HOME/my-conda/condabin/conda shell.bash hook)"
conda activate csm_analysis

# --- Log directory from base config ---
LOG_DIR=$(python3 -c "
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
python3 "${SCRIPTS_DIR}/qc_q1_guideseq.py" --config "${QC_CONFIG}" --n-workers "${NSLOTS:-8}"

echo ""
echo "Q1 Guide-Seq QC complete: $(date)"

#!/usr/bin/env bash
# Archived project-specific worker retained for provenance.
#$ -N thp11_rrna_final
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_vmem=8G
#$ -l h_rt=04:00:00

set -euo pipefail

PROJECT_DIR="${CSM_PROJECT_DIR:-/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq}"
CONFIG="${THP_RRNA_CONFIG:-${PROJECT_DIR}/config/THP_1_sublib11_rrna_audit.yaml}"
AUDIT_ROOT="${THP_RRNA_AUDIT_ROOT:-${PROJECT_DIR}/output/THP_1_sublib11/rrna_barcode_diagnostic}"

CONDA_PREFIX_CONFIG=$(awk '
    /^environment:/ { in_environment=1; next }
    in_environment && /^[^[:space:]]/ { exit }
    in_environment && /^[[:space:]]+conda_prefix:/ {
        line=$0
        sub(/^[^:]*:[[:space:]]*/, "", line)
        gsub(/"/, "", line)
        print line
        exit
    }
' "${CONFIG}")
source "${CONDA_PREFIX_CONFIG}/etc/profile.d/conda.sh"
conda activate splitpipe

python "${PROJECT_DIR}/scripts/external_rrna_audit.py" finalize \
    --project THP_1_sublib11 \
    --output-dir "${AUDIT_ROOT}"

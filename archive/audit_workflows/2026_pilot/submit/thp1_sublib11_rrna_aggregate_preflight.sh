#!/usr/bin/env bash
# Archived project-specific worker retained for provenance.
#$ -N thp11_rrna_preagg
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_vmem=4G
#$ -l h_rt=02:00:00

set -euo pipefail

PROJECT_DIR="${CSM_PROJECT_DIR:-/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq}"
CONFIG="${THP_RRNA_CONFIG:-${PROJECT_DIR}/config/THP_1_sublib11_rrna_audit.yaml}"
FASTQ_MANIFEST="${THP_RRNA_FASTQ_MANIFEST:-${PROJECT_DIR}/config/THP_1_sublib11_fastq_manifest.csv}"
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

python "${PROJECT_DIR}/scripts/external_rrna_audit.py" aggregate-preflight \
    --project THP_1_sublib11 \
    --fastq-manifest "${FASTQ_MANIFEST}" \
    --preflight-dir "${AUDIT_ROOT}/preflight" \
    --output-dir "${AUDIT_ROOT}" \
    --expected-pairs 16 \
    --expected-r1-length 54 \
    --expected-r2-length 61

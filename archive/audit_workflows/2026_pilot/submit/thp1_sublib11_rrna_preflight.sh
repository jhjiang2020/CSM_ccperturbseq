#!/usr/bin/env bash
# Archived project-specific worker retained for provenance.
#$ -N thp11_rrna_pre
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -binding linear:4
#$ -l h_vmem=4G
#$ -l h_rt=48:00:00
#$ -t 1-16
#$ -tc 4

set -euo pipefail

PROJECT_DIR="${CSM_PROJECT_DIR:-/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq}"
CONFIG="${THP_RRNA_CONFIG:-${PROJECT_DIR}/config/THP_1_sublib11_rrna_audit.yaml}"
FASTQ_MANIFEST="${THP_RRNA_FASTQ_MANIFEST:-${PROJECT_DIR}/config/THP_1_sublib11_fastq_manifest.csv}"
AUDIT_ROOT="${THP_RRNA_AUDIT_ROOT:-${PROJECT_DIR}/output/THP_1_sublib11/rrna_barcode_diagnostic}"
TASK_ID="${SGE_TASK_ID:?Run as a 16-task UGER array}"
THREADS="${NSLOTS:?UGER did not provide NSLOTS}"

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
CONDA_ENV=$(awk '
    /^environment:/ { in_environment=1; next }
    in_environment && /^[^[:space:]]/ { exit }
    in_environment && /^[[:space:]]+conda_env:/ {
        line=$0
        sub(/^[^:]*:[[:space:]]*/, "", line)
        gsub(/"/, "", line)
        print line
        exit
    }
' "${CONFIG}")
source "${CONDA_PREFIX_CONFIG}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

echo "THP_1_sublib11 FASTQ preflight task ${TASK_ID}; slots=${THREADS}"
grep Cpus_allowed_list /proc/self/status

python "${PROJECT_DIR}/scripts/external_rrna_audit.py" preflight-fastq \
    --fastq-manifest "${FASTQ_MANIFEST}" \
    --task-index "${TASK_ID}" \
    --output-dir "${AUDIT_ROOT}/preflight" \
    --threads "${THREADS}" \
    --expected-r1-length 54 \
    --expected-r2-length 61

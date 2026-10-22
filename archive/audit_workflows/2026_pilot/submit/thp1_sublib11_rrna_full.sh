#!/usr/bin/env bash
# Archived project-specific worker retained for provenance.
#$ -N thp11_rrna_full
#$ -cwd
#$ -j y
#$ -pe smp 16
#$ -binding linear:16
#$ -l h_vmem=8G
#$ -l h_rt=168:00:00

set -euo pipefail

PROJECT_DIR="${CSM_PROJECT_DIR:-/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq}"
CONFIG="${THP_RRNA_CONFIG:-${PROJECT_DIR}/config/THP_1_sublib11_rrna_audit.yaml}"
FASTQ_MANIFEST="${THP_RRNA_FASTQ_MANIFEST:-${PROJECT_DIR}/config/THP_1_sublib11_fastq_manifest.csv}"
AUDIT_ROOT="${THP_RRNA_AUDIT_ROOT:-${PROJECT_DIR}/output/THP_1_sublib11/rrna_barcode_diagnostic}"
REFERENCE_DIR="${THP_RRNA_REFERENCE_DIR:-${PROJECT_DIR}/output/pilot_4_sublib_audit/rrna_barcode_diagnostic/reference}"
SCRATCH_ROOT="${THP_RRNA_SCRATCH_ROOT:-/broad/hptmp/jj/THP_1_sublib11_rrna}"
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
source "${CONDA_PREFIX_CONFIG}/etc/profile.d/conda.sh"
conda activate splitpipe
export PYTHONHASHSEED=0
python -c 'import numpy, pandas, pysam, yaml'
command -v STAR >/dev/null
mkdir -p "${SCRATCH_ROOT}"

echo "THP_1_sublib11 full rRNA audit; slots=${THREADS}"
grep Cpus_allowed_list /proc/self/status

python "${PROJECT_DIR}/scripts/rrna_barcode_diagnostic.py" full-audit \
    --config "${CONFIG}" \
    --library THP_1_sublib11 \
    --fastq-manifest "${FASTQ_MANIFEST}" \
    --expected-records-manifest "${AUDIT_ROOT}/fastq_summary.json" \
    --barcode-contract-mode runtime-kit \
    --sample-wells A1-H12 \
    --reference-dir "${REFERENCE_DIR}" \
    --output-dir "${AUDIT_ROOT}" \
    --direct-output \
    --scratch-dir "${SCRATCH_ROOT}" \
    --threads "${THREADS}" \
    --minimum-r1-length 45 \
    --r2-rrna-mode disabled-short-read \
    --overlap-mode disabled \
    --force

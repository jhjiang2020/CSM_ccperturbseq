#!/usr/bin/env bash
# Archived single-use smoke-report workflow.
# Publish the aggregate-only interim report from the retained C1 100k smoke run.
#$ -N csm_rrna_sr
#$ -cwd
#$ -j y
#$ -o /medpop/rgupta/Jiahao/logs/CSM_ccperturbseq/
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_vmem=8G
#$ -l h_rt=01:00:00

set -euo pipefail

PROJECT_DIR="${CSM_PROJECT_DIR:-/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq}"
CONFIG="${CSM_PILOT_AUDIT_CONFIG:-${PROJECT_DIR}/config.yaml}"
SMOKE_ROOT="${CSM_RRNA_SMOKE_ROOT:-/broad/hptmp/jj/rrna_polyt_examples_20260729}"
SOURCE_DIR="${SMOKE_ROOT}/output/C1"
OUTPUT_DIR="${PROJECT_DIR}/output/pilot_4_sublib_audit/rrna_barcode_diagnostic/smoke_C1_100k"
REPORT_PATH="${PROJECT_DIR}/docs/PILOT_100K_RRNA_BARCODE_SMOKE_REPORT.md"

set +u
source /broad/software/scripts/useuse
set -u
reuse -q UGER || true

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

python "${PROJECT_DIR}/scripts/rrna_barcode_diagnostic.py" smoke-report \
    --smoke-dir "${SOURCE_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --report-path "${REPORT_PATH}" \
    --force

test -s "${OUTPUT_DIR}/smoke_validation.json"
test -s "${REPORT_PATH}"
echo "Interim C1 100k rRNA/barcode report published."

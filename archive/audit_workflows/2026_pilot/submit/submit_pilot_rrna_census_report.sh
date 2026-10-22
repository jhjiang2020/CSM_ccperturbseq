#!/usr/bin/env bash
# Archived worker retained for audit provenance.
# Publish the standalone four-sublibrary rRNA/barcode census report from the
# reconciled aggregate tables.
#$ -N csm_rrna_cr
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
RRNA_ROOT="${CSM_RRNA_AUDIT_ROOT:-${PROJECT_DIR}/output/pilot_4_sublib_audit/rrna_barcode_diagnostic}"
AGGREGATE_DIR="${RRNA_ROOT}/aggregate"
REPORT_PATH="${PROJECT_DIR}/docs/PILOT_4_SUBLIB_RRNA_BARCODE_REPORT.md"

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
[[ -n "${CONDA_PREFIX_CONFIG}" ]]
source "${CONDA_PREFIX_CONFIG}/etc/profile.d/conda.sh"
conda activate splitpipe
python -c 'import numpy, pandas'

python "${PROJECT_DIR}/scripts/rrna_barcode_diagnostic.py" census-report \
    --aggregate-dir "${AGGREGATE_DIR}" \
    --library-dir "${RRNA_ROOT}" \
    --report-path "${REPORT_PATH}" \
    --force

test -s "${AGGREGATE_DIR}/census_validation.json"
test -s "${AGGREGATE_DIR}/census_report_manifest.json"
test -s "${REPORT_PATH}"
echo "Four-sublibrary rRNA/barcode census report published."

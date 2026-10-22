#!/usr/bin/env bash
# Archived worker retained for audit provenance.
# Reconcile the four completed tasks, publish compact assets, and update Markdown.
#$ -N csm_rrna_final
#$ -cwd
#$ -j y
#$ -o /medpop/rgupta/Jiahao/logs/CSM_ccperturbseq/
#$ -pe smp 2
#$ -binding linear:2
#$ -l h_vmem=12G
#$ -l h_rt=08:00:00

set -euo pipefail

PROJECT_DIR="${CSM_PROJECT_DIR:-/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq}"
CONFIG="${CSM_PILOT_AUDIT_CONFIG:-${PROJECT_DIR}/config.yaml}"
RRNA_ROOT="${CSM_RRNA_AUDIT_ROOT:-${PROJECT_DIR}/output/pilot_4_sublib_audit/rrna_barcode_diagnostic}"
THREADS="${NSLOTS:?UGER did not provide NSLOTS}"

set +u
source /broad/software/scripts/useuse
set -u
reuse -q UGER || true
command -v qsub >/dev/null

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
python -c 'import matplotlib, numpy, pandas, yaml'

echo "Finalizer slots/binding workers: ${THREADS}"
python "${PROJECT_DIR}/scripts/rrna_barcode_diagnostic.py" aggregate \
    --output-dir "${RRNA_ROOT}" \
    --asset-dir "${PROJECT_DIR}/docs/assets/pilot_4_sublib_audit" \
    --report-path "${PROJECT_DIR}/docs/PILOT_4_SUBLIB_TECHNICAL_AUDIT.md" \
    --library-metrics "${PROJECT_DIR}/docs/assets/pilot_4_sublib_audit/library_metrics.csv" \
    --downsampling-dir "${PROJECT_DIR}/docs/assets/pilot_4_sublib_audit/downsampling"

python "${PROJECT_DIR}/scripts/rrna_barcode_diagnostic.py" validate \
    --output-dir "${RRNA_ROOT}" \
    --asset-dir "${PROJECT_DIR}/docs/assets/pilot_4_sublib_audit" \
    --report-path "${PROJECT_DIR}/docs/PILOT_4_SUBLIB_TECHNICAL_AUDIT.md"

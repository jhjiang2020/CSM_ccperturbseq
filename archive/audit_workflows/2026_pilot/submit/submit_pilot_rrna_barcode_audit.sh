#!/usr/bin/env bash
# Archived entry point; use submit/submit_gex_rrna_audit.sh.
# Submit the reference -> four-library array -> aggregate/validate dependency chain.
# Run on a Broad login node with:
#   bash submit/submit_pilot_rrna_barcode_audit.sh

set -euo pipefail

PROJECT_DIR="${CSM_PROJECT_DIR:-/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq}"
LOG_DIR="/medpop/rgupta/Jiahao/logs/CSM_ccperturbseq"
mkdir -p "${LOG_DIR}"

set +u
source /broad/software/scripts/useuse
set -u
reuse -q UGER || true
command -v qsub >/dev/null

REFERENCE_JOB=$(qsub -terse "${PROJECT_DIR}/submit/submit_pilot_rrna_reference.sh")
REFERENCE_ID="${REFERENCE_JOB%%.*}"
ARRAY_JOB=$(qsub -terse -hold_jid "${REFERENCE_ID}" "${PROJECT_DIR}/submit/submit_pilot_rrna_full_array.sh")
ARRAY_ID="${ARRAY_JOB%%.*}"
FINAL_JOB=$(qsub -terse -hold_jid "${ARRAY_ID}" "${PROJECT_DIR}/submit/submit_pilot_rrna_finalize.sh")
FINAL_ID="${FINAL_JOB%%.*}"

echo "Reference job: ${REFERENCE_JOB}"
echo "Full audit:    ${ARRAY_JOB}"
echo "Finalizer:     ${FINAL_JOB}"
echo
echo "Binding request recorded for full array:"
qstat -j "${ARRAY_ID}" | awk '/parallel environment|binding|job_number/ {print}'

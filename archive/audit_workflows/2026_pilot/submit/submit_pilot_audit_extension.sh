#!/usr/bin/env bash
# Archived single-use entry point; use submit/submit_gex_grna_read_recovery.sh.
# Submit preparation, eight rarefaction tasks, four raw scans, and held finalizer.

set -euo pipefail

PROJECT_DIR="${CSM_PROJECT_DIR:-/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq}"

set +u
source /broad/software/scripts/useuse
set -u
reuse -q UGER || true
command -v qsub >/dev/null

submit_job() {
    local output
    output=$("$@")
    echo "${output}" >&2
    awk '{ for (i=1; i<=NF; i++) if ($i ~ /^[0-9]+([.][0-9:-]+)?$/) { split($i,a,"."); print a[1]; exit } }' <<<"${output}"
}

PREPARE_ID=$(submit_job qsub "${PROJECT_DIR}/submit/submit_pilot_screen_audit_prepare.sh")
[[ "${PREPARE_ID}" =~ ^[0-9]+$ ]]

RAREFACTION_ID=$(submit_job qsub \
    -hold_jid "${PREPARE_ID}" \
    -v "CSM_PILOT_AUDIT_SKIP_FASTQ_COUNT=true" \
    "${PROJECT_DIR}/submit/submit_pilot_screen_audit.sh")
[[ "${RAREFACTION_ID}" =~ ^[0-9]+$ ]]

RAW_ID=$(submit_job qsub \
    -hold_jid "${PREPARE_ID}" \
    "${PROJECT_DIR}/submit/submit_pilot_raw_guide_diagnostic.sh")
[[ "${RAW_ID}" =~ ^[0-9]+$ ]]

FINAL_ID=$(submit_job qsub \
    -hold_jid "${RAREFACTION_ID},${RAW_ID}" \
    "${PROJECT_DIR}/submit/submit_pilot_screen_audit_finalize.sh")
[[ "${FINAL_ID}" =~ ^[0-9]+$ ]]

echo "prepare=${PREPARE_ID}"
echo "rarefaction=${RAREFACTION_ID}"
echo "raw_diagnostic=${RAW_ID}"
echo "finalizer=${FINAL_ID}"

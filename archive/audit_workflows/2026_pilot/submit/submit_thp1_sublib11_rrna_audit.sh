#!/usr/bin/env bash
# Archived project-specific entry point; use the manifest profile of submit/submit_gex_rrna_audit.sh.
set -euo pipefail

PROJECT_DIR="${CSM_PROJECT_DIR:-/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq}"
AUDIT_ROOT="${THP_RRNA_AUDIT_ROOT:-${PROJECT_DIR}/output/THP_1_sublib11/rrna_barcode_diagnostic}"
LOG_DIR="${AUDIT_ROOT}/logs"
FASTQ_MANIFEST="${THP_RRNA_FASTQ_MANIFEST:-${PROJECT_DIR}/config/THP_1_sublib11_fastq_manifest.csv}"
REFERENCE_DIR="${THP_RRNA_REFERENCE_DIR:-${PROJECT_DIR}/output/pilot_4_sublib_audit/rrna_barcode_diagnostic/reference}"

set +u
source /broad/software/scripts/useuse
set -u
reuse -q UGER || true
command -v qsub >/dev/null

[[ -f "${FASTQ_MANIFEST}" ]]
[[ -f "${REFERENCE_DIR}/reference_manifest.json" ]]
[[ ! -f "${AUDIT_ROOT}/run_manifest.json" ]]
mkdir -p "${AUDIT_ROOT}/preflight" "${LOG_DIR}"

PREFLIGHT_SPEC=$(qsub -terse -o "${LOG_DIR}" \
    "${PROJECT_DIR}/submit/thp1_sublib11_rrna_preflight.sh")
PREFLIGHT_JOB="${PREFLIGHT_SPEC%%.*}"

AGGREGATE_SPEC=$(qsub -terse -hold_jid "${PREFLIGHT_JOB}" -o "${LOG_DIR}" \
    "${PROJECT_DIR}/submit/thp1_sublib11_rrna_aggregate_preflight.sh")
AGGREGATE_JOB="${AGGREGATE_SPEC%%.*}"

AUDIT_SPEC=$(qsub -terse -hold_jid "${AGGREGATE_JOB}" -o "${LOG_DIR}" \
    "${PROJECT_DIR}/submit/thp1_sublib11_rrna_full.sh")
AUDIT_JOB="${AUDIT_SPEC%%.*}"

FINALIZE_SPEC=$(qsub -terse -hold_jid "${AUDIT_JOB}" -o "${LOG_DIR}" \
    "${PROJECT_DIR}/submit/thp1_sublib11_rrna_finalize.sh")
FINALIZE_JOB="${FINALIZE_SPEC%%.*}"

python - "${AUDIT_ROOT}/submission_manifest.json" \
    "${PREFLIGHT_SPEC}" "${AGGREGATE_SPEC}" "${AUDIT_SPEC}" "${FINALIZE_SPEC}" <<'PY'
import json
import sys
import time
from pathlib import Path

destination, preflight, aggregate, audit, finalize = sys.argv[1:]
payload = {
    "schema_version": 1,
    "project": "THP_1_sublib11",
    "submitted_unix": time.time(),
    "preflight_array_job": preflight,
    "aggregate_preflight_job": aggregate,
    "full_audit_job": audit,
    "finalizer_job": finalize,
    "dependency_order": [
        "preflight_array_job",
        "aggregate_preflight_job",
        "full_audit_job",
        "finalizer_job",
    ],
}
Path(destination).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY

echo "preflight=${PREFLIGHT_SPEC}"
echo "aggregate=${AGGREGATE_SPEC}"
echo "audit=${AUDIT_SPEC}"
echo "finalize=${FINALIZE_SPEC}"
echo "submission_manifest=${AUDIT_ROOT}/submission_manifest.json"

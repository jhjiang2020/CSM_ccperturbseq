#!/usr/bin/env bash
# Archived single-use example extraction workflow.
# Extract anonymized sequence examples from a kept 100k-pair smoke run.
#$ -N csm_rrna_ex
#$ -cwd
#$ -j y
#$ -o /medpop/rgupta/Jiahao/logs/CSM_ccperturbseq/
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_vmem=8G
#$ -l h_rt=02:00:00

set -euo pipefail

PROJECT_DIR="${CSM_PROJECT_DIR:-/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq}"
CONFIG="${CSM_PILOT_AUDIT_CONFIG:-${PROJECT_DIR}/config.yaml}"
SMOKE_ROOT="${CSM_RRNA_SMOKE_ROOT:?Set CSM_RRNA_SMOKE_ROOT}"
REFERENCE_DIR="${PROJECT_DIR}/output/pilot_4_sublib_audit/rrna_barcode_diagnostic/reference"

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

python "${PROJECT_DIR}/scripts/rrna_probe_examples.py" \
    --config "${CONFIG}" \
    --library C1 \
    --smoke-root "${SMOKE_ROOT}" \
    --reference-dir "${REFERENCE_DIR}" \
    --examples 4

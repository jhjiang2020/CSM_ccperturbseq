#!/usr/bin/env bash
# Archived single-use smoke workflow.
# End-to-end 100k-pair C1 smoke test using only UGER scratch outputs.
#$ -N csm_rrna_smoke
#$ -cwd
#$ -j y
#$ -o /medpop/rgupta/Jiahao/logs/CSM_ccperturbseq/
#$ -pe smp 4
#$ -binding linear:4
#$ -l h_vmem=8G
#$ -l h_rt=04:00:00

set -euo pipefail

PROJECT_DIR="${CSM_PROJECT_DIR:-/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq}"
CONFIG="${CSM_PILOT_AUDIT_CONFIG:-${PROJECT_DIR}/config.yaml}"
REFERENCE_DIR="${PROJECT_DIR}/output/pilot_4_sublib_audit/rrna_barcode_diagnostic/reference"
SMOKE_ROOT="${CSM_RRNA_SMOKE_ROOT:-/broad/hptmp/jj/rrna_barcode_smoke}"
THREADS="${NSLOTS:?UGER did not provide NSLOTS}"
EXTRA_ARGS=()
if [[ "${KEEP_SMOKE_SCRATCH:-0}" == "1" ]]; then
    EXTRA_ARGS+=(--keep-scratch)
fi

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
export PYTHONHASHSEED=0

python "${PROJECT_DIR}/scripts/rrna_barcode_diagnostic.py" full-audit \
    --config "${CONFIG}" \
    --library C1 \
    --reference-dir "${REFERENCE_DIR}" \
    --output-dir "${SMOKE_ROOT}/output" \
    --scratch-dir "${SMOKE_ROOT}/scratch" \
    --threads "${THREADS}" \
    --max-records 100000 \
    --progress-every 25000 \
    --valid-r2-sample-modulus 20 \
    --quality-sample-modulus 100 \
    --reconciliation-sample-modulus 50 \
    --overlap-sample-modulus 100 \
    --force \
    "${EXTRA_ARGS[@]}"

test -s "${SMOKE_ROOT}/output/C1/population_summary.csv"
test -s "${SMOKE_ROOT}/output/C1/star_multimap_reconciliation.csv"
test -s "${SMOKE_ROOT}/output/C1/paired_overlap_summary.csv"
echo "End-to-end rRNA smoke test passed."

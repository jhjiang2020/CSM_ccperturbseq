#!/usr/bin/env bash
# Archived single-use representative-read finalizer.
# Aggregate selected examples after the four scan tasks finish.
# Submit with:
#   qsub -hold_jid ARRAY_JOB_ID submit/submit_bc1_representative_reads_finalize.sh
#$ -N csm_bc1_examples_final
#$ -cwd
#$ -j y
#$ -o /medpop/rgupta/Jiahao/logs/CSM_ccperturbseq/
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_vmem=4G
#$ -l h_rt=01:00:00

set -euo pipefail

PROJECT_DIR="${CSM_PROJECT_DIR:-/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq}"
CONFIG="${CSM_PILOT_AUDIT_CONFIG:-${PROJECT_DIR}/config.yaml}"
WORK_DIR="${CSM_BC1_EXAMPLE_WORK_DIR:-${PROJECT_DIR}/output/pilot_4_sublib_audit/bc1_examples/work}"
OUTPUT_DIR="${CSM_BC1_EXAMPLE_OUTPUT_DIR:-${PROJECT_DIR}/output/pilot_4_sublib_audit/bc1_examples}"

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

python "${PROJECT_DIR}/scripts/extract_representative_bc1_failures.py" aggregate \
    --config "${CONFIG}" \
    --work-dir "${WORK_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --top-sequences 25 \
    --top-examples 8

echo "Representative BC1-failure examples aggregated."

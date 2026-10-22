#!/usr/bin/env bash
# Archived single-use barcode-rejection workflow.
# Lightweight paired-arm barcode rejection replay (one sublibrary per task).
#$ -N csm_barcode_reject
#$ -cwd
#$ -j y
#$ -o /medpop/rgupta/Jiahao/logs/CSM_ccperturbseq/
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_vmem=8G
#$ -l h_rt=08:00:00
#$ -t 1-4
#$ -tc 4

set -euo pipefail

PROJECT_DIR="${CSM_PROJECT_DIR:-/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq}"
CONFIG="${CSM_PILOT_AUDIT_CONFIG:-${PROJECT_DIR}/config.yaml}"
WORK_DIR="${CSM_PILOT_AUDIT_WORK_DIR:-${PROJECT_DIR}/output/pilot_4_sublib_audit/work}"
TASK_ID="${SGE_TASK_ID:?Run as a UGER array task (1-4)}"
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
CONDA_ENV=$(awk '
    /^environment:/ { in_environment=1; next }
    in_environment && /^[^[:space:]]/ { exit }
    in_environment && /^[[:space:]]+conda_env:/ {
        line=$0
        sub(/^[^:]*:[[:space:]]*/, "", line)
        gsub(/"/, "", line)
        print line
        exit
    }
' "${CONFIG}")
[[ -n "${CONDA_PREFIX_CONFIG}" && -n "${CONDA_ENV}" ]]
source "${CONDA_PREFIX_CONFIG}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"
python -c 'import matplotlib, numpy, pandas, yaml'

LIB_INDEX=$(( TASK_ID - 1 ))
LIBRARY=$(python -c "import yaml; c=yaml.safe_load(open('${CONFIG}')); print(c['libraries'][${LIB_INDEX}]['sublibrary_name'])")

echo "Barcode rejection task ${TASK_ID}: ${LIBRARY}"
echo "Slots/binding workers: ${THREADS}"

python "${PROJECT_DIR}/scripts/barcode_rejection_qc.py" audit \
    --config "${CONFIG}" \
    --work-dir "${WORK_DIR}" \
    --library "${LIBRARY}"

echo "Barcode rejection task ${TASK_ID} complete."

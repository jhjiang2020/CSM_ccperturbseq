#!/usr/bin/env bash
# Archived worker retained for audit provenance.
# Run one heavy arm × sublibrary audit task (tasks 1-4 GEX, 5-8 gRNA).
#$ -N csm_pilot_audit
#$ -cwd
#$ -j y
#$ -o /medpop/rgupta/Jiahao/logs/CSM_ccperturbseq/
#$ -pe smp 4
#$ -binding linear:4
#$ -l h_vmem=8G
#$ -l h_rt=48:00:00
#$ -t 1-8
#$ -tc 4

set -euo pipefail

PROJECT_DIR="${CSM_PROJECT_DIR:-/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq}"
CONFIG="${CSM_PILOT_AUDIT_CONFIG:-${PROJECT_DIR}/config.yaml}"
REPORT_DIR="${CSM_PILOT_AUDIT_REPORT_DIR:-${PROJECT_DIR}/docs/assets/pilot_4_sublib_audit}"
WORK_DIR="${CSM_PILOT_AUDIT_WORK_DIR:-${PROJECT_DIR}/output/pilot_4_sublib_audit/work}"
TASK_ID="${SGE_TASK_ID:?Run as a UGER array task (1-8)}"
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

LIB_INDEX=$(( (TASK_ID - 1) % 4 ))
ARM_INDEX=$(( (TASK_ID - 1) / 4 ))
LIBRARY=$(python -c "import yaml; c=yaml.safe_load(open('${CONFIG}')); print(c['libraries'][${LIB_INDEX}]['sublibrary_name'])")
if [[ "${ARM_INDEX}" -eq 0 ]]; then
    ARM="GEX"
else
    ARM="gRNA"
fi

echo "Pilot audit task ${TASK_ID}: ${LIBRARY} ${ARM}"
echo "Slots/binding workers: ${THREADS}"

if [[ "${CSM_PILOT_AUDIT_SKIP_FASTQ_COUNT:-false}" != "true" ]]; then
    python "${PROJECT_DIR}/scripts/pilot_screen_audit.py" count-fastq \
        --config "${CONFIG}" \
        --report-dir "${REPORT_DIR}" \
        --library "${LIBRARY}" \
        --arm "${ARM}" \
        --threads "${THREADS}"
fi

python "${PROJECT_DIR}/scripts/pilot_screen_audit.py" downsample \
    --config "${CONFIG}" \
    --work-dir "${WORK_DIR}" \
    --report-dir "${REPORT_DIR}" \
    --library "${LIBRARY}" \
    --arm "${ARM}"

echo "Pilot audit task ${TASK_ID} complete."

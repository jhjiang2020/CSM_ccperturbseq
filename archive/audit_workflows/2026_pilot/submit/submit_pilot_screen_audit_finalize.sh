#!/usr/bin/env bash
# Archived worker retained for audit provenance.
# Finalize the compact report after all FASTQ/downsampling array tasks complete.
# Submit with: qsub -hold_jid JOB_ID submit/submit_pilot_screen_audit_finalize.sh
#$ -N csm_pilot_audit_final
#$ -cwd
#$ -j y
#$ -o /medpop/rgupta/Jiahao/logs/CSM_ccperturbseq/
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_vmem=16G
#$ -l h_rt=04:00:00

set -euo pipefail

PROJECT_DIR="${CSM_PROJECT_DIR:-/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq}"
CONFIG="${CSM_PILOT_AUDIT_CONFIG:-${PROJECT_DIR}/config.yaml}"
REPORT_DIR="${CSM_PILOT_AUDIT_REPORT_DIR:-${PROJECT_DIR}/docs/assets/pilot_4_sublib_audit}"
WORK_DIR="${CSM_PILOT_AUDIT_WORK_DIR:-${PROJECT_DIR}/output/pilot_4_sublib_audit/work}"
REPORT_PATH="${CSM_PILOT_AUDIT_REPORT_PATH:-${PROJECT_DIR}/docs/PILOT_4_SUBLIB_TECHNICAL_AUDIT.md}"

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
python -c 'import h5py, matplotlib, numpy, pandas, yaml'

python "${PROJECT_DIR}/scripts/pilot_screen_audit.py" summarize \
    --config "${CONFIG}" \
    --report-dir "${REPORT_DIR}" \
    --work-dir "${WORK_DIR}" \
    --report-path "${REPORT_PATH}" \
    --require-matrix-validation \
    --require-downsampling \
    --require-raw-diagnostic

echo "Pilot audit finalization complete: ${REPORT_PATH}"

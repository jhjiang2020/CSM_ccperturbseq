#!/bin/bash
#$ -N csm_asym_v2
#$ -q broad
#$ -l os=RedHat8,h_vmem=8G,h_rt=1:00:00
#$ -pe smp 2
#$ -binding linear:2
#$ -j y

set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the shared checkout}"
: "${RUN_DIR:?Set RUN_DIR to a fresh durable validation directory}"
: "${FIXTURE_DIR:?Set FIXTURE_DIR to a fresh fixture directory, or an existing one for validate}"
: "${MODE:=validate}"
: "${CONDA_ENV:=splitpipe}"
: "${CONDA_EXE:?Set the absolute conda executable}"
mkdir -p "${RUN_DIR}"
trap 'status=$?; if (( status != 0 )); then printf "%s exit=%s\n" "$(date -u +%FT%TZ)" "$status" > "${RUN_DIR}/_FAILED.${JOB_ID}"; fi' EXIT
eval "$("${CONDA_EXE}" shell.bash hook)"
conda activate "${CONDA_ENV}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
export MPLCONFIGDIR="${TMPDIR:-/tmp}/csm-asym-mpl-${JOB_ID}"
export XDG_CACHE_HOME="${TMPDIR:-/tmp}/csm-asym-cache-${JOB_ID}"
cd "${PROJECT_ROOT}"
echo "job=${JOB_ID} host=$(hostname) python=$(command -v python) env=${CONDA_PREFIX} mode=${MODE}"
grep Cpus_allowed_list /proc/self/status
if [[ "${MODE}" == build ]]; then
    python -B scripts/validate_asymmetric_pilot.py build --config "${PROJECT_ROOT}/config.yaml" --destination "${FIXTURE_DIR}"
elif [[ "${MODE}" == validate ]]; then
    python -B scripts/validate_asymmetric_pilot.py validate --fixture "${FIXTURE_DIR}" --output "${RUN_DIR}/validation" --workers "${NSLOTS:-1}"
else
    echo "Unknown MODE=${MODE}" >&2
    exit 2
fi
date -u +%FT%TZ > "${RUN_DIR}/_SUCCESS.${JOB_ID}"

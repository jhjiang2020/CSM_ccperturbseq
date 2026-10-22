#!/bin/bash
# Extract compact C2 evidence for the gRNA-to-GEX misalignment case study.
#$ -q broad
#$ -N csm_c2_gex_case
#$ -j y
#$ -l h_vmem=4G
#$ -l h_rt=04:00:00
#$ -pe smp 4
#$ -binding linear:4

set -euo pipefail
# Explicit tools use NSLOTS; nested numerical pools use one thread.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${CSM_PROJECT_DIR:-$(dirname "${SUBMIT_DIR}")}"
QC_CONFIG="${CSM_C2_CASE_CONFIG:-${PROJECT_DIR}/output/v1_4_performance_audit/v14_qc_completion_20260720_02/configs/qc_C2.yaml}"
CASE_OUTPUT="${CSM_C2_CASE_OUTPUT:-${PROJECT_DIR}/output/case_studies/c2_gex_grna_misalignment/data}"
SCRATCH_ROOT="${CSM_C2_CASE_SCRATCH:-/broad/hptmp/jj/CSM_ccperturbseq/c2_gex_grna_case_study_${JOB_ID:-manual}}"

if [[ ! -r "${QC_CONFIG}" ]]; then
    echo "ERROR: Completed C2 QC config is not readable: ${QC_CONFIG}" >&2
    exit 1
fi

BASE_CONFIG_VALUE=$(grep '^ *base_config:' "${QC_CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' || true)
if [[ -z "${BASE_CONFIG_VALUE}" ]]; then
    echo "ERROR: base_config is not set in ${QC_CONFIG}" >&2
    exit 1
fi
if [[ "${BASE_CONFIG_VALUE}" = /* ]]; then
    BASE_CONFIG="${BASE_CONFIG_VALUE}"
else
    BASE_CONFIG="$(cd "$(dirname "${QC_CONFIG}")" && pwd)/${BASE_CONFIG_VALUE}"
fi
if [[ ! -r "${BASE_CONFIG}" ]]; then
    echo "ERROR: Base config is not readable: ${BASE_CONFIG}" >&2
    exit 1
fi

CONDA_PREFIX=$(grep '^ *conda_prefix:' "${BASE_CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' | sed 's:/*$::' || true)
CONDA_ENV=$(grep '^ *conda_env:' "${BASE_CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' || true)

if [[ -z "${CONDA_PREFIX}" || -z "${CONDA_ENV}" ]]; then
    echo "ERROR: environment.conda_prefix and environment.conda_env are required" >&2
    exit 1
fi
if [[ ! -x "${CONDA_PREFIX}/condabin/conda" ]]; then
    echo "ERROR: Conda executable is unavailable under ${CONDA_PREFIX}" >&2
    exit 1
fi

eval "$("${CONDA_PREFIX}/condabin/conda" shell.bash hook)"
conda activate "${CONDA_ENV}"

set +u
source /broad/software/scripts/useuse
set -u
reuse -q Samtools || true
command -v samtools
command -v pigz
PYTHON_BIN=$(command -v python3 2>/dev/null || command -v python)
"${PYTHON_BIN}" -c 'import pandas, pysam, splitpipe, yaml'

mkdir -p "${SCRATCH_ROOT}" "${CASE_OUTPUT}"

echo "=========================================="
echo " C2 gRNA-to-GEX case-study extraction"
echo " Host:       $(hostname)"
echo " Date:       $(date)"
echo " Config:     ${QC_CONFIG}"
echo " Output:     ${CASE_OUTPUT}"
echo " Scratch:    ${SCRATCH_ROOT}"
echo " Slots:      ${NSLOTS:-4}"
grep -E '^(Cpus_allowed_list|Mems_allowed_list):' /proc/self/status || true
echo "=========================================="

"${PYTHON_BIN}" "${PROJECT_DIR}/scripts/c2_gex_grna_case_study.py" extract \
    --config "${QC_CONFIG}" \
    --output-dir "${CASE_OUTPUT}" \
    --scratch-dir "${SCRATCH_ROOT}" \
    --threads "${NSLOTS:-4}" \
    --samtools "$(command -v samtools)" \
    --pigz "$(command -v pigz)"

echo "Case-study extraction complete: $(date)"

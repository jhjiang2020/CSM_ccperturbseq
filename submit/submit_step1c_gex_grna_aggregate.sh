#!/bin/bash
# ==============================================================================
# Step 1c aggregate: validate every library and build the GEX-only identity arm
# ==============================================================================
#$ -q broad
#$ -N csm_step1c_aggregate
#$ -j y
#$ -l h_vmem=32G
#$ -l h_rt=8:00:00
#$ -pe smp 4
#$ -binding linear:4

set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

if [[ -n "${CSM_CONFIG:-}" ]]; then
    CONFIG="${CSM_CONFIG}"
else
    SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    CONFIG="$(dirname "${SUBMIT_DIR}")/config.yaml"
fi

if [[ -n "${CSM_PROJECT_DIR:-}" ]]; then
    PROJECT_DIR="${CSM_PROJECT_DIR}"
else
    PROJECT_DIR=$(grep '^ *project_dir:' "${CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/')
fi

_CONDA_PREFIX=$(grep '^ *conda_prefix:' "${CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/')
_CONDA_ENV=$(grep '^ *conda_env:' "${CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/')
[[ -n "${_CONDA_PREFIX}" && -n "${_CONDA_ENV}" ]] || { echo "ERROR: Configure environment.conda_prefix and environment.conda_env" >&2; exit 1; }
CONDA_PREFIX="${_CONDA_PREFIX}"
CONDA_ENV="${_CONDA_ENV}"
[[ -x "${CONDA_PREFIX}/condabin/conda" ]] || { echo "ERROR: Missing configured Conda: ${CONDA_PREFIX}/condabin/conda" >&2; exit 1; }
eval "$("${CONDA_PREFIX}/condabin/conda" shell.bash hook)"
conda activate "${CONDA_ENV}"

echo "=========================================="
echo "  Step 1c aggregate: GEX guide identity"
echo "  Host: $(hostname)"
echo "  Date: $(date)"
echo "  Job:  ${JOB_ID:-local}"
echo "=========================================="

python3 "${PROJECT_DIR}/scripts/validate_step1c_manifests.py" --config "${CONFIG}"
python3 "${PROJECT_DIR}/scripts/04_assign_grna_to_cells.py" \
    --config "${CONFIG}" \
    --source-mode gex \
    --ncpus "${NSLOTS:-4}"

echo "Step 1c aggregate complete: $(date)"

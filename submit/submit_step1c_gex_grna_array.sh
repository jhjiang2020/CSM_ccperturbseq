#!/bin/bash
# ==============================================================================
# Step 1c: run ordinary Step 3 on compact Step 1a evidence (one task/library)
# ==============================================================================
#$ -q broad
#$ -N csm_step1c_grna
#$ -j y
#$ -l h_vmem=4G
#$ -l h_rt=12:00:00
#$ -pe smp 1
#$ -binding linear:1
# Array range, dependency on Step 1a, and log path are supplied by qsub.

set -euo pipefail

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

LIB_NAME=$(python3 -c "
import yaml
with open('${CONFIG}') as handle:
    libraries = yaml.safe_load(handle)['libraries']
index = ${SGE_TASK_ID} - 1
if index < 0 or index >= len(libraries):
    raise IndexError('SGE_TASK_ID is outside configured libraries')
print(libraries[index]['sublibrary_name'])
")

OUTPUT_DIR=$(python3 -c "
import yaml
with open('${CONFIG}') as handle:
    print(yaml.safe_load(handle)['paths']['output_dir'])
")
MANIFEST_ROOT="${CSM_STEP1C_MANIFEST_ROOT:-${OUTPUT_DIR}/gex_grna_filter}"
OUTPUT_NAMESPACE="${CSM_STEP1C_OUTPUT_NAMESPACE:-gex_grna_identity}"
INPUT_MANIFEST="${MANIFEST_ROOT}/${LIB_NAME}/rejected_pairs.complete.json"

echo "=========================================="
echo "  Step 1c: GEX-derived gRNA identity"
echo "  Host:       $(hostname)"
echo "  Date:       $(date)"
echo "  Job:        ${JOB_ID}.${SGE_TASK_ID}"
echo "  Sublibrary: ${LIB_NAME}"
echo "  Manifest:   ${INPUT_MANIFEST}"
echo "  Namespace:  ${OUTPUT_NAMESPACE}"
echo "=========================================="

python "${PROJECT_DIR}/scripts/01c_process_gex_grna.py" \
    --config "${CONFIG}" \
    --library "${LIB_NAME}" \
    --input-manifest "${INPUT_MANIFEST}" \
    --output-namespace "${OUTPUT_NAMESPACE}" \
    --n_processes "${NSLOTS:-1}"

echo "Step 1c complete for ${LIB_NAME}: $(date)"

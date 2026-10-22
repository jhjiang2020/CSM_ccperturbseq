#!/bin/bash
# ==============================================================================
# Step 4: Assign gRNAs to Cells (all libraries — run AFTER step 3 completes)
# ==============================================================================
# Submit via submit_pipeline.sh (recommended), or manually:
#   qsub -hold_jid csm_step3_match submit/submit_step4_assign.sh
# ==============================================================================
# UGE Directives:
#$ -q broad
#$ -N csm_step4_assign
#$ -j y
#$ -l h_vmem=32G
#$ -l h_rt=4:00:00
#$ -pe smp 4
#$ -binding linear:4

set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

# --- Configuration ---
if [[ -n "${CSM_CONFIG:-}" ]]; then
    CONFIG="${CSM_CONFIG}"
else
    SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    CONFIG="$(dirname "${SUBMIT_DIR}")/config.yaml"
fi

# --- Resolve project dir: env var > config (grep) > dirname(config) ---
if [[ -n "${CSM_PROJECT_DIR:-}" ]]; then
    PROJECT_DIR="${CSM_PROJECT_DIR}"
else
    PROJECT_DIR=$(grep '^ *project_dir:' "${CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/')
fi
SCRIPTS_DIR="${PROJECT_DIR}/scripts"

# --- Initialize conda (parsed from config via grep — no python3 needed yet) ---
_CONDA_PREFIX=$(grep '^ *conda_prefix:' "${CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/')
_CONDA_ENV=$(grep '^ *conda_env:' "${CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/')
[[ -n "${_CONDA_PREFIX}" && -n "${_CONDA_ENV}" ]] || { echo "ERROR: Configure environment.conda_prefix and environment.conda_env" >&2; exit 1; }
CONDA_PREFIX="${_CONDA_PREFIX}"
CONDA_ENV="${_CONDA_ENV}"
[[ -x "${CONDA_PREFIX}/condabin/conda" ]] || { echo "ERROR: Missing configured Conda: ${CONDA_PREFIX}/condabin/conda" >&2; exit 1; }
eval "$(${CONDA_PREFIX}/condabin/conda shell.bash hook)"
conda activate "${CONDA_ENV}"

LOG_DIR=$(python3 -c "import yaml; print(yaml.safe_load(open('${CONFIG}'))['paths'].get('log_dir', '${PROJECT_DIR}/logs'))")
mkdir -p "${LOG_DIR}"

echo "=========================================="
echo "  Step 4: Assign gRNAs to Cells"
echo "  Host: $(hostname)"
echo "  Date: $(date)"
echo "  Job:  ${JOB_ID:-local}"
echo "=========================================="

# --- Run ---
# When the current master submission already completed the Step 1c aggregate,
# avoid recomputing the GEX-only matrix. A standalone Step 4 submission keeps
# auto mode and reconstructs every configured product safely.
if [[ "${CSM_STEP1C_GEX_READY:-false}" == "true" ]]; then
    python "${SCRIPTS_DIR}/04_assign_grna_to_cells.py" \
        --config "${CONFIG}" --ncpus "${NSLOTS:-4}" --source-mode dedicated
    MERGE_ENABLED=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); print(cfg.get('gex_decontamination', {}).get('combine_with_grna_arm', {}).get('enabled', False))")
    if [[ "${MERGE_ENABLED}" == "True" ]]; then
        python "${SCRIPTS_DIR}/04_assign_grna_to_cells.py" \
            --config "${CONFIG}" --ncpus "${NSLOTS:-4}" --source-mode combined
    fi
else
    python "${SCRIPTS_DIR}/04_assign_grna_to_cells.py" \
        --config "${CONFIG}" --ncpus "${NSLOTS:-4}" --source-mode auto
fi

echo ""
echo "Step 4 complete: $(date)"

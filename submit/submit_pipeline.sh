#!/bin/bash
# ==============================================================================
# Master Submission Script — ccPerturbSeq Pipeline
# ==============================================================================
# Reads the config to determine the number of libraries, then submits all
# pipeline steps with correct array ranges and dependency chains.
#
# Usage:
#   bash submit/submit_pipeline.sh [--config CONFIG] [--step STEP] [--dry-run]
#
# Examples:
#   bash submit/submit_pipeline.sh                  # Submit all steps
#   bash submit/submit_pipeline.sh --step 1         # Step 1a/1b + optional Step 1c
#   bash submit/submit_pipeline.sh --step 2,3       # Submit Steps 2 and 3
#   bash submit/submit_pipeline.sh --dry-run        # Show what would be submitted
#   bash submit/submit_pipeline.sh --config /path/to/config.yaml
# ==============================================================================

set -euo pipefail

# --- Defaults (derived from script location) ---
SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$(dirname "${SUBMIT_DIR}")/config.yaml"
STEPS="all"
DRY_RUN=false

# --- Parse arguments ---
while [[ $# -gt 0 ]]; do
    case $1 in
        --config)   CONFIG="$2"; shift 2 ;;
        --step)     STEPS="$2"; shift 2 ;;
        --dry-run)  DRY_RUN=true; shift ;;
        -h|--help)
            echo "Usage: bash submit_pipeline.sh [--config CONFIG] [--step STEP] [--dry-run]"
            echo "  --config CONFIG   Path to config.yaml (default: ${CONFIG})"
            echo "  --step STEP       Comma-separated steps to submit: 1,2,3,4,5 or 'all'"
            echo "  --dry-run         Print qsub commands without executing"
            exit 0
            ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [[ "${DRY_RUN}" != true ]] && ! command -v qsub >/dev/null 2>&1; then
    echo "ERROR: qsub is unavailable. On Broad run:" >&2
    echo "       source /broad/software/scripts/useuse; reuse -q UGER || true; command -v qsub" >&2
    exit 1
fi

# --- Bootstrap only the lightweight YAML reader on the login node ---
if [[ ! -r "${CONFIG}" ]]; then
    echo "ERROR: Config is not readable: ${CONFIG}" >&2
    exit 1
fi
_CONDA_PREFIX=$(grep '^ *conda_prefix:' "${CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' || true)
_CONDA_ENV=$(grep '^ *conda_env:' "${CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' || true)
CONDA_ROOT="${_CONDA_PREFIX%/}"
PYTHON_BIN="$(command -v python3 2>/dev/null || true)"
if [[ -z "${PYTHON_BIN}" ]] || ! "${PYTHON_BIN}" -c 'import yaml' >/dev/null 2>&1; then
    if [[ -z "${CONDA_ROOT}" || -z "${_CONDA_ENV}" || ! -x "${CONDA_ROOT}/condabin/conda" ]]; then
        echo "ERROR: Python/PyYAML is unavailable and the configured Conda environment cannot be activated" >&2
        echo "       Check environment.conda_prefix and environment.conda_env in ${CONFIG}" >&2
        exit 1
    fi
    eval "$("${CONDA_ROOT}/condabin/conda" shell.bash hook)"
    conda activate "${_CONDA_ENV}"
    PYTHON_BIN="$(command -v python3 2>/dev/null || true)"
fi
if [[ -z "${PYTHON_BIN}" ]] || ! "${PYTHON_BIN}" -c 'import yaml' >/dev/null 2>&1; then
    echo "ERROR: Python with PyYAML is unavailable after activating ${_CONDA_ENV}" >&2
    exit 1
fi

# --- Resolve project dir from config (grep — no pyyaml needed) ---
PROJECT_DIR=$(grep '^ *project_dir:' "${CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/')
JOB_SUBMIT_DIR="${PROJECT_DIR%/}/submit"

N_LIBS=$(python3 -c "
import yaml
with open('${CONFIG}') as f:
    cfg = yaml.safe_load(f)
print(len(cfg['libraries']))
")

LIB_NAMES=$(python3 -c "
import yaml
with open('${CONFIG}') as f:
    cfg = yaml.safe_load(f)
for lib in cfg['libraries']:
    print(f\"  {lib['sublibrary_name']}\")
")

LOG_DIR=$(python3 -c "
import yaml
with open('${CONFIG}') as f:
    cfg = yaml.safe_load(f)
print(cfg['paths'].get('log_dir', '${PROJECT_DIR}/logs'))
")
LOG_DIR="${LOG_DIR%/}"
STEP1C_ENABLED=$(python3 -c "
import yaml
with open('${CONFIG}') as f:
    cfg = yaml.safe_load(f)
print(cfg.get('gex_decontamination', {}).get('step1c', {}).get('enabled', False))
")
GEX_SOURCE_MERGE_ENABLED=$(python3 -c "
import yaml
with open('${CONFIG}') as f:
    cfg = yaml.safe_load(f)
print(cfg.get('gex_decontamination', {}).get('combine_with_grna_arm', {}).get('enabled', False))
")
RESOURCE_RESOLVER="${PROJECT_DIR}/scripts/uger_resources.py"
if [[ ! -r "${RESOURCE_RESOLVER}" ]]; then
    RESOURCE_RESOLVER="${SUBMIT_DIR}/../scripts/uger_resources.py"
fi
IFS=$'\t' read -r STEP1A_SLOTS STEP1A_MEMORY_PER_SLOT_GB STEP1A_RUNTIME_HOURS STEP1A_H_RT < <(
    python3 "${RESOURCE_RESOLVER}" --config "${CONFIG}" --format tsv
)
if [[ "${DRY_RUN}" != true ]]; then
    mkdir -p "${LOG_DIR}"
fi

echo "=========================================="
echo "  ccPerturbSeq Pipeline Submission"
echo "=========================================="
echo "  Config:      ${CONFIG}"
echo "  Project dir: ${PROJECT_DIR}"
echo "  Libraries:   ${N_LIBS}"
echo "${LIB_NAMES}"
echo "  Log dir:    ${LOG_DIR}"
echo "  Steps:      ${STEPS}"
echo "  Dry run:    ${DRY_RUN}"
echo "  Step 1c:    ${STEP1C_ENABLED}"
echo "  GEX merge:  ${GEX_SOURCE_MERGE_ENABLED}"
echo "  Step 1a:    ${STEP1A_SLOTS} slots, ${STEP1A_MEMORY_PER_SLOT_GB}G/slot, ${STEP1A_H_RT}"
echo "=========================================="
echo ""

# --- Helper to submit or print ---
submit() {
    local cmd="$1"
    if [[ "$DRY_RUN" == true ]]; then
        echo "  [DRY RUN] ${cmd}"
    else
        echo "  $ ${cmd}"
        eval "${cmd}"
    fi
}

should_run() {
    local step="$1"
    [[ "$STEPS" == "all" ]] || [[ ",$STEPS," == *",$step,"* ]]
}

# --- Submit Steps ---

# Step 1a: GEX split-pipe (array job)
if should_run 1; then
    echo "--- Step 1a: GEX split-pipe (array ×${N_LIBS}) ---"
    submit "qsub -pe smp ${STEP1A_SLOTS} -binding linear:${STEP1A_SLOTS} -l h_vmem=${STEP1A_MEMORY_PER_SLOT_GB}G,h_rt=${STEP1A_H_RT} -v CSM_CONFIG=${CONFIG},CSM_PROJECT_DIR=${PROJECT_DIR} -t 1-${N_LIBS} -o ${LOG_DIR}/step1a_gex.\\\$JOB_ID.\\\$TASK_ID.log ${JOB_SUBMIT_DIR}/submit_step1a_gex_array.sh"
    echo ""

    # Step 1b: Combine (after 1a)
    echo "--- Step 1b: Combine GEX sublibraries ---"
    submit "qsub -v CSM_CONFIG=${CONFIG},CSM_PROJECT_DIR=${PROJECT_DIR} -hold_jid csm_step1a_gex -o ${LOG_DIR}/step1b_combine.\\\$JOB_ID.log ${JOB_SUBMIT_DIR}/submit_step1b_gex_combine.sh"
    echo ""

    if [[ "${STEP1C_ENABLED}" == "True" ]]; then
        echo "--- Step 1c: Process GEX-derived gRNA reads (array ×${N_LIBS}) ---"
        submit "qsub -v CSM_CONFIG=${CONFIG},CSM_PROJECT_DIR=${PROJECT_DIR} -hold_jid csm_step1a_gex -t 1-${N_LIBS} -o ${LOG_DIR}/step1c_grna.\\\$JOB_ID.\\\$TASK_ID.log ${JOB_SUBMIT_DIR}/submit_step1c_gex_grna_array.sh"
        submit "qsub -v CSM_CONFIG=${CONFIG},CSM_PROJECT_DIR=${PROJECT_DIR} -hold_jid csm_step1c_grna -o ${LOG_DIR}/step1c_aggregate.\\\$JOB_ID.log ${JOB_SUBMIT_DIR}/submit_step1c_gex_grna_aggregate.sh"
        echo ""
    fi
fi

# Step 2: Barcode extraction (array job)
if should_run 2; then
    echo "--- Step 2: Extract gRNA barcodes (array ×${N_LIBS}) ---"
    submit "qsub -v CSM_CONFIG=${CONFIG},CSM_PROJECT_DIR=${PROJECT_DIR} -t 1-${N_LIBS} -o ${LOG_DIR}/step2_bc.\\\$JOB_ID.\\\$TASK_ID.log ${JOB_SUBMIT_DIR}/submit_step2_barcodes.sh"
    echo ""
fi

# Step 3: gRNA matching (array job, after Step 2)
if should_run 3; then
    echo "--- Step 3: Match gRNA + iBAR (array ×${N_LIBS}) ---"
    if should_run 2; then
        submit "qsub -v CSM_CONFIG=${CONFIG},CSM_PROJECT_DIR=${PROJECT_DIR} -hold_jid csm_step2_bc -t 1-${N_LIBS} -o ${LOG_DIR}/step3_match.\\\$JOB_ID.\\\$TASK_ID.log ${JOB_SUBMIT_DIR}/submit_step3_match.sh"
    else
        submit "qsub -v CSM_CONFIG=${CONFIG},CSM_PROJECT_DIR=${PROJECT_DIR} -t 1-${N_LIBS} -o ${LOG_DIR}/step3_match.\\\$JOB_ID.\\\$TASK_ID.log ${JOB_SUBMIT_DIR}/submit_step3_match.sh"
    fi
    echo ""
fi

# Step 4: Assign gRNAs to cells (single job, after Step 3)
if should_run 4; then
    echo "--- Step 4: Assign gRNAs to cells ---"
    step4_holds=""
    if should_run 3; then
        step4_holds="csm_step3_match"
    fi
    if [[ "${STEP1C_ENABLED}" == "True" ]] && should_run 1; then
        if [[ -n "${step4_holds}" ]]; then
            step4_holds="${step4_holds},csm_step1c_aggregate"
        else
            step4_holds="csm_step1c_aggregate"
        fi
    fi
    step4_extra_env=""
    if [[ "${STEP1C_ENABLED}" == "True" ]] && should_run 1; then
        step4_extra_env=",CSM_STEP1C_GEX_READY=true"
    fi
    if [[ -n "${step4_holds}" ]]; then
        submit "qsub -v CSM_CONFIG=${CONFIG},CSM_PROJECT_DIR=${PROJECT_DIR}${step4_extra_env} -hold_jid ${step4_holds} -o ${LOG_DIR}/step4_assign.\\\$JOB_ID.log ${JOB_SUBMIT_DIR}/submit_step4_assign.sh"
    else
        submit "qsub -v CSM_CONFIG=${CONFIG},CSM_PROJECT_DIR=${PROJECT_DIR}${step4_extra_env} -o ${LOG_DIR}/step4_assign.\\\$JOB_ID.log ${JOB_SUBMIT_DIR}/submit_step4_assign.sh"
    fi
    echo ""
fi

# Step 5: Build AnnData (after Steps 1b and 4)
if should_run 5; then
    echo "--- Step 5: Build combined AnnData ---"
    local_holds=""
    if should_run 1; then local_holds="csm_step1b_combine"; fi
    if should_run 4; then
        if [[ -n "$local_holds" ]]; then
            local_holds="${local_holds},csm_step4_assign"
        else
            local_holds="csm_step4_assign"
        fi
    fi
    if [[ -n "$local_holds" ]]; then
        submit "qsub -v CSM_CONFIG=${CONFIG},CSM_PROJECT_DIR=${PROJECT_DIR} -hold_jid ${local_holds} -o ${LOG_DIR}/step5_anndata.\\\$JOB_ID.log ${JOB_SUBMIT_DIR}/submit_step5_anndata.sh"
    else
        submit "qsub -v CSM_CONFIG=${CONFIG},CSM_PROJECT_DIR=${PROJECT_DIR} -o ${LOG_DIR}/step5_anndata.\\\$JOB_ID.log ${JOB_SUBMIT_DIR}/submit_step5_anndata.sh"
    fi
    echo ""
fi

echo "=========================================="
echo "  Submission complete!"
echo "  Log dir:  ${LOG_DIR}"
echo "  Monitor:  qstat -u \$USER"
echo "=========================================="

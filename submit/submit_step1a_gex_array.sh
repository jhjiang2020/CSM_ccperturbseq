#!/bin/bash
# ==============================================================================
# Step 1a: Run GEX split-pipe per sublibrary (ARRAY JOB)
# ==============================================================================
# Submit via submit_pipeline.sh (recommended), or manually:
#   qsub -t 1-N submit/submit_step1a_gex_array.sh
# ==============================================================================
# UGE Directives:
#$ -q broad
#$ -N csm_step1a_gex
#$ -j y
#$ -l h_vmem=8G
#$ -l h_rt=48:00:00
#$ -pe smp 16
#$ -binding linear:16
# NOTE: -t (array range) and -o (log dir) are set dynamically by submit_pipeline.sh

set -euo pipefail

# --- Configuration ---
# CSM_CONFIG is passed via qsub -v from submit_pipeline.sh.
# Falls back to BASH_SOURCE derivation for local/interactive runs.
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
SCRIPTS_DIR="${PROJECT_DIR}/scripts"

# --- Initialize conda (parsed from config via grep — no python3 needed yet) ---
_CONDA_PREFIX=$(grep '^ *conda_prefix:' "${CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/')
_CONDA_ENV=$(grep '^ *conda_env:' "${CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/')
if [[ -z "${_CONDA_PREFIX}" || -z "${_CONDA_ENV}" ]]; then
    echo "ERROR: environment.conda_prefix and environment.conda_env are required in ${CONFIG}" >&2
    exit 1
fi
CONDA_PREFIX="${_CONDA_PREFIX}"
CONDA_ENV="${_CONDA_ENV}"
if [[ ! -x "${CONDA_PREFIX}/condabin/conda" ]]; then
    echo "ERROR: configured conda executable is missing: ${CONDA_PREFIX}/condabin/conda" >&2
    exit 1
fi
eval "$(${CONDA_PREFIX}/condabin/conda shell.bash hook)"
conda activate "${CONDA_ENV}"

LOG_DIR=$(python3 -c "import yaml; print(yaml.safe_load(open('${CONFIG}'))['paths'].get('log_dir', '$(dirname \"${CONFIG}\")/logs'))")
mkdir -p "${LOG_DIR}"

# --- Map array task ID to sublibrary ---
# Read library info from config
LIB_NAME=$(python3 -c "
import yaml
with open('${CONFIG}') as f:
    cfg = yaml.safe_load(f)
libs = cfg['libraries']
idx = ${SGE_TASK_ID} - 1
if idx >= len(libs):
    raise IndexError(f'Task ID {${SGE_TASK_ID}} exceeds number of libraries ({len(libs)})')
print(libs[idx]['sublibrary_name'])
")

GEX_PREFIX=$(python3 -c "
import yaml
with open('${CONFIG}') as f:
    cfg = yaml.safe_load(f)
print(cfg['libraries'][${SGE_TASK_ID}-1]['gex_prefix'])
")

SAMPLE_NAME=$(python3 -c "
import yaml
with open('${CONFIG}') as f:
    cfg = yaml.safe_load(f)
print(cfg['libraries'][${SGE_TASK_ID}-1]['sample_name'])
")

# Read other config values
FASTQ_DIR=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); print(cfg['paths']['fastq_dir'])")
SPLITPIPE_DIR=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); print(cfg['paths']['splitpipe_dir'])")
GENOME_DIR=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); print(cfg['paths']['genome_dir'])")
OUTPUT_DIR=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); print(cfg['paths']['output_dir'])")
KIT=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); print(cfg['splitpipe']['kit'])")
CHEMISTRY=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); print(cfg['splitpipe']['chemistry'])")
SAMPLE_WELLS=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); print(cfg['splitpipe']['sample_wells'])")
BYPASS_KIT=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); print(cfg['splitpipe'].get('bypass_kit_score', False))")
GEX_DECONTAM_ENABLED=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); print(cfg.get('gex_decontamination', {}).get('enabled', False))")
EXISTING_ALIGNMENT_ENABLED=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); print(cfg.get('gex_decontamination', {}).get('existing_alignment', {}).get('enabled', False))")

SPLITPIPE="${SPLITPIPE_DIR}/split-pipe"
SPLITPIPE_GEX_WRAPPER="${SCRIPTS_DIR}/splitpipe_gex_wrapper.py"
SPLITPIPE_GEX_PATCH="${PROJECT_DIR}/patches/parsebio_splitpipe_1.6.4_gex_decontamination.patch"
OUT_SUBLIB="${OUTPUT_DIR}/gex_${LIB_NAME}"

if [[ "${EXISTING_ALIGNMENT_ENABLED}" == "True" ]]; then
    SOURCE_SUBLIB=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); section=cfg['gex_decontamination']['existing_alignment']; print(section.get('source_output_dir_template', '{output_dir}/gex_{library}').format(output_dir=cfg['paths']['output_dir'], library='${LIB_NAME}'))")
    OUT_SUBLIB=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); section=cfg['gex_decontamination']['existing_alignment']; print(section.get('sandbox_output_dir_template', '{output_dir}/gex_{library}_v1_4').format(output_dir=cfg['paths']['output_dir'], library='${LIB_NAME}'))")
    BAM_BACKEND=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${CONFIG}')); print(cfg['gex_decontamination']['existing_alignment'].get('bam_backend', 'auto'))")
fi

# The runtime overlay is entirely environment-gated. With the option disabled,
# Step 1a executes pristine Parse v1.6.4 directly. With it enabled, the wrapper
# validates and transforms preprocess.py/molinfo.py only in memory.
if [[ "${GEX_DECONTAM_ENABLED}" == "True" ]]; then
    if [[ ! -f "${SPLITPIPE_GEX_WRAPPER}" || ! -f "${SPLITPIPE_GEX_PATCH}" ]]; then
        echo "ERROR: Missing v1.4 split-pipe runtime overlay wrapper or patch."
        echo "  Wrapper: ${SPLITPIPE_GEX_WRAPPER}"
        echo "  Patch:   ${SPLITPIPE_GEX_PATCH}"
        exit 1
    fi
    export CSM_GEX_DECONTAM_CONFIG="${CONFIG}"
    export CSM_GEX_DECONTAM_LIBRARY="${LIB_NAME}"
    export CSM_GEX_DECONTAM_OUTPUT_DIR="${OUTPUT_DIR}/gex_grna_filter/${LIB_NAME}"
    export CSM_GEX_DECONTAM_HOOK="${SCRIPTS_DIR}/gex_grna_preprocess_hook.py"
    export CSM_GEX_DECONTAM_RUN_ID="${JOB_ID:-local}.${SGE_TASK_ID:-1}.$(date +%s)"
    export CSM_PRESERVE_HARMONIZED_BAM=1
    echo "  GEX gRNA decontamination: enabled (runtime overlay, packed inline, before STAR)"
else
    unset CSM_GEX_DECONTAM_CONFIG CSM_GEX_DECONTAM_LIBRARY
    unset CSM_GEX_DECONTAM_OUTPUT_DIR CSM_GEX_DECONTAM_HOOK CSM_GEX_DECONTAM_RUN_ID
    unset CSM_PRESERVE_HARMONIZED_BAM
    echo "  GEX gRNA decontamination: disabled"
fi

run_splitpipe() {
    if [[ "${GEX_DECONTAM_ENABLED}" == "True" ]]; then
        python3 "${SPLITPIPE_GEX_WRAPPER}" \
            --wrapper_splitpipe_dir "${SPLITPIPE_DIR}" \
            --wrapper_patch "${SPLITPIPE_GEX_PATCH}" \
            "$@"
    else
        "${SPLITPIPE}" "$@"
    fi
}

# Build optional flags
EXTRA_FLAGS=""
if [[ "$BYPASS_KIT" == "True" ]]; then
    EXTRA_FLAGS="--kit_score_skip"
    echo "  (kit score check will be skipped)"
fi

echo "=========================================="
echo "  Step 1a: GEX split-pipe (task ${SGE_TASK_ID})"
echo "  Host:       $(hostname)"
echo "  Date:       $(date)"
echo "  Job:        ${JOB_ID}.${SGE_TASK_ID}"
echo "  Sublibrary: ${LIB_NAME}"
echo "  GEX prefix: ${GEX_PREFIX}"
echo "  Sample:     ${SAMPLE_NAME}"
echo "  Cores:      ${NSLOTS:-16}"
echo "=========================================="

# Find FASTQ files
mapfile -t FQ1_FILES < <(find "${FASTQ_DIR}" -name "*${GEX_PREFIX}*_R1_*fastq.gz" -print | sort -V)
mapfile -t FQ2_FILES < <(find "${FASTQ_DIR}" -name "*${GEX_PREFIX}*_R2_*fastq.gz" -print | sort -V)

if [[ ${#FQ1_FILES[@]} -eq 0 || ${#FQ1_FILES[@]} -ne ${#FQ2_FILES[@]} ]]; then
    echo "ERROR: Could not find FASTQ files for GEX prefix '${GEX_PREFIX}' in ${FASTQ_DIR}"
    exit 1
fi

printf 'R1 file: %s\n' "${FQ1_FILES[@]}"
printf 'R2 file: %s\n' "${FQ2_FILES[@]}"
echo "Output: ${OUT_SUBLIB}"
echo ""

if [[ "${GEX_DECONTAM_ENABLED}" == "True" ]]; then
    INPUT_IDENTITY="${CSM_GEX_DECONTAM_OUTPUT_DIR}/step1a_input_identity.json"
    identity_args=(
        --library "${LIB_NAME}"
        --output "${INPUT_IDENTITY}"
    )
    for path in "${FQ1_FILES[@]}"; do identity_args+=(--r1 "${path}"); done
    for path in "${FQ2_FILES[@]}"; do identity_args+=(--r2 "${path}"); done
    python3 "${SCRIPTS_DIR}/write_fastq_input_identity.py" "${identity_args[@]}"
    export CSM_GEX_DECONTAM_INPUT_MANIFEST="${INPUT_IDENTITY}"
fi

if [[ "${EXISTING_ALIGNMENT_ENABLED}" == "True" ]]; then
    if [[ "${GEX_DECONTAM_ENABLED}" != "True" ]]; then
        echo "ERROR: existing_alignment.enabled requires gex_decontamination.enabled"
        exit 1
    fi
    if [[ "${SOURCE_SUBLIB}" == "${OUT_SUBLIB}" ]]; then
        echo "ERROR: source_output_dir_template and sandbox_output_dir_template resolve to the same path"
        exit 1
    fi
    SOURCE_BAM=""
    for candidate in \
        "${SOURCE_SUBLIB}/process/barcode_headAligned_anno.bam" \
        "${SOURCE_SUBLIB}/process/barcode_headAligned_sorted.bam"; do
        if [[ -f "${candidate}" ]]; then
            SOURCE_BAM="${candidate}"
            break
        fi
    done
    if [[ -z "${SOURCE_BAM}" ]]; then
        echo "ERROR: No compatible name-sorted/annotated BAM found under ${SOURCE_SUBLIB}/process"
        exit 1
    fi
    SOURCE_PRE_ALIGN_STATS="${SOURCE_SUBLIB}/process/pre_align_stats.csv"
    if [[ ! -f "${SOURCE_PRE_ALIGN_STATS}" ]]; then
        echo "ERROR: Historical alignment statistics are missing: ${SOURCE_PRE_ALIGN_STATS}"
        exit 1
    fi
    if [[ -e "${OUT_SUBLIB}/process/barcode_headAligned_sorted.bam" ]]; then
        echo "ERROR: Existing-alignment sandbox already contains a harmonized BAM: ${OUT_SUBLIB}"
        exit 1
    fi

    echo "Legacy source (read-only): ${SOURCE_BAM}"
    echo "Harmonized sandbox:        ${OUT_SUBLIB}"

    # Re-run only preprocessing in the new sandbox to obtain the same positive
    # QNAME evidence and filtered barcode_head FASTQ. STAR is not rerun.
    run_splitpipe \
        --mode pre --one_step \
        --kit "${KIT}" \
        --chemistry "${CHEMISTRY}" \
        --genome_dir "${GENOME_DIR}" \
        --fq1 "${FQ1_FILES[@]}" \
        --fq2 "${FQ2_FILES[@]}" \
        --output_dir "${OUT_SUBLIB}" \
        --sample "${SAMPLE_NAME}" "${SAMPLE_WELLS}" \
        --nthreads "${NSLOTS:-16}" \
        ${EXTRA_FLAGS}

    python3 "${SCRIPTS_DIR}/gex_grna_bam_filter.py" \
        --source-bam "${SOURCE_BAM}" \
        --source-pre-align-stats "${SOURCE_PRE_ALIGN_STATS}" \
        --output-dir "${OUT_SUBLIB}" \
        --filter-manifest "${OUTPUT_DIR}/gex_grna_filter/${LIB_NAME}/rejected_pairs.complete.json" \
        --library "${LIB_NAME}" \
        --backend "${BAM_BACKEND}"

    # The sandbox now satisfies split-pipe's ordinary mol input contract.
    run_splitpipe \
        --mode mol \
        --output_dir "${OUT_SUBLIB}" \
        --nthreads "${NSLOTS:-16}" \
        --save_anndata \
        ${EXTRA_FLAGS}
else
    run_splitpipe \
        --mode all \
        --kit "${KIT}" \
        --chemistry "${CHEMISTRY}" \
        --genome_dir "${GENOME_DIR}" \
        --fq1 "${FQ1_FILES[@]}" \
        --fq2 "${FQ2_FILES[@]}" \
        --output_dir "${OUT_SUBLIB}" \
        --sample "${SAMPLE_NAME}" "${SAMPLE_WELLS}" \
        --nthreads "${NSLOTS:-16}" \
        --save_anndata \
        ${EXTRA_FLAGS}
    if [[ "${GEX_DECONTAM_ENABLED}" == "True" ]]; then
        python3 "${SCRIPTS_DIR}/gex_grna_bam_filter.py" \
            --register-fresh \
            --source-bam "${OUT_SUBLIB}/process/barcode_headAligned_sorted.bam" \
            --output-dir "${OUT_SUBLIB}" \
            --filter-manifest "${OUTPUT_DIR}/gex_grna_filter/${LIB_NAME}/rejected_pairs.complete.json" \
            --library "${LIB_NAME}"
    fi
fi

echo ""
echo "Step 1a complete for ${LIB_NAME}: $(date)"

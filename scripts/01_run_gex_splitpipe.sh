#!/bin/bash
# ==============================================================================
# Step 1: Run standard split-pipe GEX processing
# ==============================================================================
# Usage: bash scripts/01_run_gex_splitpipe.sh --config config.yaml
#
# This script:
#   1. Runs split-pipe in --mode all for each GEX sublibrary
#   2. Combines results with split-pipe --mode comb
# ==============================================================================

set -euo pipefail

# --- Parse arguments ---
CONFIG_FILE=""
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --config) CONFIG_FILE="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [[ -z "$CONFIG_FILE" ]]; then
    echo "Usage: bash scripts/01_run_gex_splitpipe.sh --config config.yaml"
    exit 1
fi

# v1.4 safety: this historical sequential driver called split-pipe directly
# and would bypass the preprocessing hook, Step 1c, manifest validation, and
# existing-BAM sandbox. Delegate to the config-driven UGER orchestrator.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_SUBMIT="$(dirname "${SCRIPT_DIR}")/submit/submit_pipeline.sh"
DELEGATE_ARGS=(--config "${CONFIG_FILE}" --step 1)
if [[ "${DRY_RUN}" == true ]]; then
    DELEGATE_ARGS+=(--dry-run)
fi
echo "NOTE: scripts/01_run_gex_splitpipe.sh is deprecated in v1.4; delegating to submit_pipeline.sh."
exec bash "${PIPELINE_SUBMIT}" "${DELEGATE_ARGS[@]}"

# --- Read config with Python (handles YAML properly) ---
read_config() {
    python3 -c "
import yaml, sys
with open('$CONFIG_FILE') as f:
    cfg = yaml.safe_load(f)
# Navigate dotted key path
keys = '$1'.split('.')
val = cfg
for k in keys:
    val = val[k]
print(val)
"
}

read_config_list() {
    python3 -c "
import yaml, sys
with open('$CONFIG_FILE') as f:
    cfg = yaml.safe_load(f)
keys = '$1'.split('.')
val = cfg
for k in keys:
    if isinstance(val, list):
        val = val[int(k)]
    else:
        val = val[k]
print(val)
"
}

# --- Extract config values ---
PROJECT_DIR=$(read_config "paths.project_dir")
FASTQ_DIR=$(read_config "paths.fastq_dir")
SPLITPIPE_DIR=$(read_config "paths.splitpipe_dir")
GENOME_DIR=$(read_config "paths.genome_dir")
OUTPUT_DIR=$(read_config "paths.output_dir")
KIT=$(read_config "splitpipe.kit")
CHEMISTRY=$(read_config "splitpipe.chemistry")
SAMPLE_WELLS=$(read_config "splitpipe.sample_wells")
THREADS=$(read_config "splitpipe.threads")

SPLITPIPE="${SPLITPIPE_DIR}/split-pipe"

echo "=============================================="
echo "  Step 1: GEX Processing with split-pipe"
echo "=============================================="
echo "Project:    $PROJECT_DIR"
echo "FASTQ dir:  $FASTQ_DIR"
echo "Genome:     $GENOME_DIR"
echo "Kit:        $KIT"
echo "Chemistry:  $CHEMISTRY"
echo "Wells:      $SAMPLE_WELLS"
echo "Threads:    $THREADS"
echo ""

# Optional flags
BYPASS_KIT=$(read_config "splitpipe.bypass_kit_score" 2>/dev/null || echo "false")
EXTRA_FLAGS=""
if [[ "$BYPASS_KIT" == "True" || "$BYPASS_KIT" == "true" ]]; then
    EXTRA_FLAGS="--kit_score_skip"
    echo "Kit score check: SKIPPED (bypass_kit_score=true)"
fi
echo ""

# --- Get library info ---
N_LIBS=$(python3 -c "
import yaml
with open('$CONFIG_FILE') as f:
    cfg = yaml.safe_load(f)
print(len(cfg['libraries']))
")

echo "Found $N_LIBS GEX sublibraries to process."
echo ""

# --- Process each GEX sublibrary ---
SUBLIBRARY_DIRS=()

for ((i=0; i<N_LIBS; i++)); do
    LIB_NAME=$(python3 -c "
import yaml
with open('$CONFIG_FILE') as f:
    cfg = yaml.safe_load(f)
print(cfg['libraries'][$i]['sublibrary_name'])
")
    GEX_PREFIX=$(python3 -c "
import yaml
with open('$CONFIG_FILE') as f:
    cfg = yaml.safe_load(f)
print(cfg['libraries'][$i]['gex_prefix'])
")
    SAMPLE_NAME=$(python3 -c "
import yaml
with open('$CONFIG_FILE') as f:
    cfg = yaml.safe_load(f)
print(cfg['libraries'][$i]['sample_name'])
")

    # Find FASTQ files matching the prefix
    FQ1=$(find "$FASTQ_DIR" -name "*${GEX_PREFIX}*_R1_*fastq.gz" | head -1)
    FQ2=$(find "$FASTQ_DIR" -name "*${GEX_PREFIX}*_R2_*fastq.gz" | head -1)

    if [[ -z "$FQ1" || -z "$FQ2" ]]; then
        echo "ERROR: Could not find FASTQ files for GEX prefix '$GEX_PREFIX' in $FASTQ_DIR"
        echo "  Looking for: *${GEX_PREFIX}*_R1_*fastq.gz"
        exit 1
    fi

    OUT_SUBLIB="${OUTPUT_DIR}/gex_${LIB_NAME}"
    SUBLIBRARY_DIRS+=("$OUT_SUBLIB")

    echo "----------------------------------------------"
    echo "Processing sublibrary: $LIB_NAME"
    echo "  GEX prefix:   $GEX_PREFIX"
    echo "  Sample:        $SAMPLE_NAME"
    echo "  R1:            $FQ1"
    echo "  R2:            $FQ2"
    echo "  Output:        $OUT_SUBLIB"
    echo ""

    if [[ "$DRY_RUN" == true ]]; then
        echo "  [DRY RUN] Would run split-pipe --mode all"
        continue
    fi

    # Run split-pipe for this sublibrary
    "$SPLITPIPE" \
        --mode all \
        --kit "$KIT" \
        --chemistry "$CHEMISTRY" \
        --genome_dir "$GENOME_DIR" \
        --fq1 "$FQ1" \
        --fq2 "$FQ2" \
        --output_dir "$OUT_SUBLIB" \
        --sample "$SAMPLE_NAME" "$SAMPLE_WELLS" \
        --nthreads "$THREADS" \
        ${EXTRA_FLAGS}

    echo "  ✓ Completed sublibrary $LIB_NAME"
    echo ""
done

# --- Combine all sublibraries ---
echo "=============================================="
echo "  Combining GEX sublibraries"
echo "=============================================="
COMB_OUT="${OUTPUT_DIR}/gex_combined"

echo "  Sublibraries: ${SUBLIBRARY_DIRS[*]}"
echo "  Output:       $COMB_OUT"
echo ""

if [[ "$DRY_RUN" == true ]]; then
    echo "  [DRY RUN] Would run split-pipe --mode comb"
else
    "$SPLITPIPE" \
        --mode comb \
        --sublibraries "${SUBLIBRARY_DIRS[@]}" \
        --output_dir "$COMB_OUT"

    echo "  ✓ Combined output saved to $COMB_OUT"
fi

echo ""
echo "=============================================="
echo "  Step 1 complete!"
echo "=============================================="

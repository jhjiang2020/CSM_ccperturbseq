#!/bin/bash
# Controlled STAR benchmark: library x {full, short54} ground-truth array.
# Submit only after the raw benchmark array has written alignment_input FASTQs.
#$ -q broad
#$ -N csm_gex_grna_star
#$ -j y
#$ -l h_vmem=16G
#$ -l h_rt=12:00:00
#$ -pe smp 8
#$ -binding linear:8

set -euo pipefail

SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SUBMIT_DIR}")"
QC_CONFIG="${CSM_GEX_CONTAM_CONFIG:-${PROJECT_DIR}/qc_config_gex_grna_contamination.yaml}"

_BASE_CONFIG_VALUE=$(grep '^ *base_config:' "${QC_CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' || true)
_BASE_CONFIG_VALUE="${_BASE_CONFIG_VALUE:-config.yaml}"
if [[ "${_BASE_CONFIG_VALUE}" = /* ]]; then
    BASE_CONFIG="${_BASE_CONFIG_VALUE}"
else
    BASE_CONFIG="$(cd "$(dirname "${QC_CONFIG}")" && pwd)/${_BASE_CONFIG_VALUE}"
fi
if [[ ! -r "${BASE_CONFIG}" ]]; then
    echo "ERROR: Base config is not readable: ${BASE_CONFIG}" >&2
    exit 1
fi
_CONFIG_PROJECT_DIR=$(grep '^ *project_dir:' "${BASE_CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' || true)
PROJECT_DIR="${CSM_PROJECT_DIR:-${_CONFIG_PROJECT_DIR:-$(dirname "${BASE_CONFIG}")}}"
_CONDA_PREFIX=$(grep '^ *conda_prefix:' "${BASE_CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' || true)
_CONDA_ENV=$(grep '^ *conda_env:' "${BASE_CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/' || true)
CONDA_ROOT="${_CONDA_PREFIX%/}"
CONDA_ENV_NAME="${CSM_GEX_CONTAM_ENV:-${_CONDA_ENV}}"
if [[ -z "${CONDA_ROOT}" || -z "${CONDA_ENV_NAME}" ]]; then
    echo "ERROR: Configure environment.conda_prefix and environment.conda_env" >&2
    exit 1
fi
if [[ ! -x "${CONDA_ROOT}/condabin/conda" ]]; then
    echo "ERROR: Conda executable is not available at ${CONDA_ROOT}/condabin/conda" >&2
    echo "       Check environment.conda_prefix in ${BASE_CONFIG}" >&2
    exit 1
fi
eval "$("${CONDA_ROOT}/condabin/conda" shell.bash hook)"
conda activate "${CONDA_ENV_NAME}"
PYTHON_BIN="$(command -v python3 2>/dev/null || command -v python 2>/dev/null || true)"
if [[ -z "${PYTHON_BIN}" ]]; then
    echo "ERROR: Neither python3 nor python is available after activating ${CONDA_ENV_NAME}" >&2
    exit 1
fi
if ! "${PYTHON_BIN}" -c 'import edlib, pandas, pysam, yaml' >/dev/null 2>&1; then
    echo "ERROR: ${CONDA_ENV_NAME} is incomplete for controlled STAR QC" >&2
    echo "       Install requirements_extra.txt into the production splitpipe environment" >&2
    exit 1
fi

TASK_ID="${SGE_TASK_ID:-1}"
IFS=$'\t' read -r LIBRARY VARIANT INPUT_FASTQ OUTPUT_DIR GENOME_DIR ADAPTER < <(
    "${PYTHON_BIN}" - "${PROJECT_DIR}" "${QC_CONFIG}" "${TASK_ID}" <<'PY'
import os
import sys
import yaml
import hashlib

project, config_path, task = sys.argv[1], sys.argv[2], int(sys.argv[3]) - 1
sys.path.insert(0, os.path.join(project, "scripts"))
from qc_config_loader import load_qc_config
from qc_gex_grna_contamination import (
    benchmark_config_sha256,
    benchmark_summary_status,
)

cfg = load_qc_config(config_path)
qc = cfg["gex_grna_contamination_qc"]
align = qc["alignment_benchmark"]
libraries = [entry["sublibrary_name"] for entry in cfg["libraries"]]
read_length = int(qc.get("short_read_length", 54))
expected_short = f"short{read_length}"
variants = [str(value) for value in align.get("variants", ["full", expected_short])]
variants = [expected_short if value == "short54" else value for value in variants]
allowed = {"full", expected_short}
if not variants:
    raise ValueError("alignment_benchmark.variants must not be empty")
if len(set(variants)) != len(variants) or not set(variants).issubset(allowed):
    raise ValueError(
        f"alignment_benchmark.variants must be a unique subset of {sorted(allowed)}"
    )
total = len(libraries) * len(variants)
if task < 0 or task >= total:
    raise IndexError(f"Task {task + 1} outside 1-{total}")
library = libraries[task // len(variants)]
variant = variants[task % len(variants)]
suffix = "full" if variant == "full" else str(read_length)
root = qc["output_dir"]
input_fastq = os.path.join(root, library, "alignment_input", variant, f"{library}_ground_truth_R1_{suffix}.fastq.gz")
output_dir = os.path.join(root, library, align.get("output_subdir", "alignment_star"), variant)
adapter = align.get("star_clip3p_adapter", "GAGGTGGTTGGA")
summary_path = os.path.join(root, library, "benchmark_summary.yaml")
if not os.path.exists(summary_path):
    raise FileNotFoundError(f"Benchmark summary not found: {summary_path}")
with open(summary_path) as handle:
    summary = yaml.safe_load(handle) or {}
truth_input = summary.get("alignment_ground_truth_input", {})
module_path = os.path.join(project, "scripts", "qc_gex_grna_contamination.py")
with open(module_path, "rb") as handle:
    current_module_sha256 = hashlib.sha256(handle.read()).hexdigest()
benchmark_current, benchmark_reason = benchmark_summary_status(
    summary_path, cfg, library, current_module_sha256
)
if not benchmark_current:
    raise RuntimeError(
        f"Refusing controlled STAR: benchmark is stale ({benchmark_reason}); "
        "rerun the raw benchmark."
    )
if truth_input.get("qc_module_sha256", "") != current_module_sha256:
    raise RuntimeError(
        "Refusing controlled STAR: QC code changed after truth generation; "
        "rerun the raw benchmark."
    )
if truth_input.get("benchmark_config_sha256", "") != benchmark_config_sha256(
    cfg, library
):
    raise RuntimeError(
        "Refusing controlled STAR: QC configuration changed after truth "
        "generation; rerun the raw benchmark."
    )
if truth_input.get("safety_cap_reached", False):
    raise RuntimeError(
        "Refusing controlled STAR: alignment truth safety cap was reached. "
        "Lower alignment_benchmark.sample_fraction and rerun the raw benchmark."
    )
if int(truth_input.get("reads_written_per_variant", 0)) <= 0:
    raise RuntimeError(
        "Refusing controlled STAR: no ground-truth structural contaminant reads "
        f"were written for {library}."
    )
print("\t".join([library, variant, input_fastq, output_dir, cfg["paths"]["genome_dir"], adapter]))
PY
)

if [[ ! -s "${INPUT_FASTQ}" ]]; then
    echo "ERROR: Missing or empty ground-truth STAR input: ${INPUT_FASTQ}" >&2
    exit 1
fi
if STAR_BIN=$(command -v STAR 2>/dev/null); then
    :
elif [[ -x "${CONDA_ROOT}/bin/STAR" ]]; then
    # Parse installations commonly provide STAR in the configured Conda root
    # rather than the separate QC environment.  This remains relocatable
    # because CONDA_ROOT comes from environment.conda_prefix.
    STAR_BIN="${CONDA_ROOT}/bin/STAR"
else
    echo "ERROR: STAR is unavailable in ${CONDA_ENV_NAME} and ${CONDA_ROOT}/bin" >&2
    exit 1
fi
mkdir -p "${OUTPUT_DIR}"

echo "=========================================="
echo " Controlled GEX-contaminant STAR benchmark"
echo " Host:      $(hostname)"
echo " Date:      $(date)"
echo " Library:   ${LIBRARY}"
echo " Variant:   ${VARIANT}"
echo " Input:     ${INPUT_FASTQ}"
echo " Output:    ${OUTPUT_DIR}"
echo " STAR:      ${STAR_BIN}"
echo " Threads:   ${NSLOTS:-8}"
echo "=========================================="

"${STAR_BIN}" \
    --genomeDir "${GENOME_DIR}" \
    --runThreadN "${NSLOTS:-8}" \
    --readFilesIn "${INPUT_FASTQ}" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${OUTPUT_DIR}/${LIBRARY}_" \
    --outSAMtype BAM Unsorted \
    --outFilterMultimapNmax 3 \
    --clip3pAdapterSeq "${ADAPTER}"

"${PYTHON_BIN}" "${PROJECT_DIR}/scripts/qc_gex_grna_contamination.py" \
    --config "${QC_CONFIG}" \
    --library "${LIBRARY}" \
    --mode alignment \
    --alignment-variant "${VARIANT}" \
    --threads "${NSLOTS:-8}"

echo "Controlled alignment audit complete: $(date)"

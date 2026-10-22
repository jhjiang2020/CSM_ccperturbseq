#!/usr/bin/env bash
# Archived worker retained for audit provenance.
# Build the pinned small rRNA reference used by the full-raw pilot audit.
#$ -N csm_rrna_ref
#$ -cwd
#$ -j y
#$ -o /medpop/rgupta/Jiahao/logs/CSM_ccperturbseq/
#$ -pe smp 4
#$ -binding linear:4
#$ -l h_vmem=8G
#$ -l h_rt=04:00:00

set -euo pipefail

PROJECT_DIR="${CSM_PROJECT_DIR:-/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq}"
CONFIG="${CSM_PILOT_AUDIT_CONFIG:-${PROJECT_DIR}/config.yaml}"
RRNA_ROOT="${CSM_RRNA_AUDIT_ROOT:-${PROJECT_DIR}/output/pilot_4_sublib_audit/rrna_barcode_diagnostic}"
REFERENCE_DIR="${RRNA_ROOT}/reference"
SOURCE_DIR="${REFERENCE_DIR}/sources"
GENCODE_GTF="${SOURCE_DIR}/gencode.v44.primary_assembly.annotation.gtf.gz"
THREADS="${NSLOTS:?UGER did not provide NSLOTS}"

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
[[ -n "${CONDA_PREFIX_CONFIG}" ]]
source "${CONDA_PREFIX_CONFIG}/etc/profile.d/conda.sh"
conda activate splitpipe
python -c 'import numpy, pandas, pysam, yaml'
command -v STAR >/dev/null

mkdir -p "${SOURCE_DIR}"
if [[ ! -s "${GENCODE_GTF}" ]]; then
    DOWNLOAD="${GENCODE_GTF}.download"
    curl --fail --location --retry 5 --retry-delay 10 \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz" \
        --output "${DOWNLOAD}"
    gzip -t "${DOWNLOAD}"
    mv "${DOWNLOAD}" "${GENCODE_GTF}"
fi

if [[ -s "${REFERENCE_DIR}/star/Genome" && -s "${REFERENCE_DIR}/reference_manifest.json" ]]; then
    echo "Pinned rRNA reference already exists: ${REFERENCE_DIR}"
    exit 0
fi

python "${PROJECT_DIR}/scripts/rrna_barcode_diagnostic.py" prepare-reference \
    --config "${CONFIG}" \
    --gencode-gtf "${GENCODE_GTF}" \
    --reference-dir "${REFERENCE_DIR}" \
    --star "$(command -v STAR)" \
    --threads "${THREADS}" \
    --force

#!/usr/bin/env bash
# Archived worker retained for audit provenance.
# Full-raw GEX rRNA/barcode census: one sublibrary per bound 16-slot task.
#$ -N csm_rrna_full
#$ -cwd
#$ -j y
#$ -o /medpop/rgupta/Jiahao/logs/CSM_ccperturbseq/
#$ -pe smp 16
#$ -binding linear:16
#$ -l h_vmem=8G
#$ -l h_rt=168:00:00
#$ -t 1-4
#$ -tc 4

set -euo pipefail

PROJECT_DIR="${CSM_PROJECT_DIR:-/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq}"
CONFIG="${CSM_PILOT_AUDIT_CONFIG:-${PROJECT_DIR}/config.yaml}"
RRNA_ROOT="${CSM_RRNA_AUDIT_ROOT:-${PROJECT_DIR}/output/pilot_4_sublib_audit/rrna_barcode_diagnostic}"
REFERENCE_DIR="${RRNA_ROOT}/reference"
SCRATCH_ROOT="${CSM_RRNA_SCRATCH_ROOT:-/broad/hptmp/jj}"
TASK_ID="${SGE_TASK_ID:?Run as a UGER array task (1-4)}"
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
export PYTHONHASHSEED=0
python -c 'import numpy, pandas, pysam, yaml'
command -v STAR >/dev/null

LIB_INDEX=$(( TASK_ID - 1 ))
LIBRARY=$(python -c "import yaml; c=yaml.safe_load(open('${CONFIG}')); print(c['libraries'][${LIB_INDEX}]['sublibrary_name'])")

echo "Full-raw rRNA audit task ${TASK_ID}: ${LIBRARY}"
echo "Slots/binding workers: ${THREADS}"
grep Cpus_allowed_list /proc/self/status

python "${PROJECT_DIR}/scripts/rrna_barcode_diagnostic.py" full-audit \
    --config "${CONFIG}" \
    --library "${LIBRARY}" \
    --reference-dir "${REFERENCE_DIR}" \
    --output-dir "${RRNA_ROOT}" \
    --scratch-dir "${SCRATCH_ROOT}" \
    --threads "${THREADS}" \
    --valid-r2-sample-modulus 200 \
    --force

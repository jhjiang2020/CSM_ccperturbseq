#!/bin/bash
# Build a reference on a compute node; forward the builder's FASTA/GTF/output args.
# qsub -v CSM_CONFIG=/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq/config.yaml submit/submit_kallisto_reference.sh --fasta ... --gtf ... --genome hg38 --output ... --kb ... --kallisto ...
#$ -q broad
#$ -N csm_kallisto_ref
#$ -j y
#$ -pe smp 8
#$ -binding linear:8
#$ -l h_vmem=8G
#$ -l h_rt=12:00:00
#$ -l os=RedHat8
set -euo pipefail
# Explicit tools use NSLOTS; nested numerical pools use one thread.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [[ -n "${CSM_CONFIG:-}" ]]; then
    CONFIG="${CSM_CONFIG}"
elif [[ -n "${CSM_PROJECT_DIR:-}" ]]; then
    CONFIG="${CSM_PROJECT_DIR}/config.yaml"
else
    SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    CONFIG="$(dirname "${SUBMIT_DIR}")/config.yaml"
fi
PROJECT_DIR="${CSM_PROJECT_DIR:-$(sed -n 's/^ *project_dir: *"\{0,1\}\([^"#]*\)"\{0,1\}.*/\1/p' "${CONFIG}" | head -1)}"
_CONDA_PREFIX=$(sed -n 's/^ *conda_prefix: *"\{0,1\}\([^"#]*\)"\{0,1\}.*/\1/p' "${CONFIG}" | head -1)
_CONDA_ENV=$(sed -n 's/^ *conda_env: *"\{0,1\}\([^"#]*\)"\{0,1\}.*/\1/p' "${CONFIG}" | head -1)
[[ -x "${_CONDA_PREFIX}/condabin/conda" && -n "${_CONDA_ENV}" ]]
eval "$("${_CONDA_PREFIX}/condabin/conda" shell.bash hook)"
conda activate "${_CONDA_ENV}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 "${PROJECT_DIR}/scripts/build_kallisto_reference.py" --threads "${NSLOTS:-8}" "$@"

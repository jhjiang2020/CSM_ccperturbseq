#!/bin/bash
# Scheduler worker for the matched four-library backend benchmark.
#$ -q broad
#$ -N csm_kallisto_bench
#$ -j y
#$ -pe smp 16
#$ -binding linear:16
#$ -l h_vmem=8G
#$ -l h_rt=12:00:00
#$ -l os=RedHat8
set -euo pipefail
: "${CSM_PROJECT_DIR:?}"
: "${CSM_BENCH_DIR:?}"
: "${CSM_BENCH_MODE:?reference, prepare, run, or report}"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/csm-kallisto-mpl-${JOB_ID:-local}"
cd "${CSM_PROJECT_DIR}"
BENCH_PYTHON="${CSM_BENCH_PYTHON:-/medpop/rgupta/Jiahao/conda/envs/splitpipe/bin/python}"
CONFIG="${CSM_CONFIG:-${CSM_PROJECT_DIR}/config.yaml}"
_CONDA_PREFIX=$(grep '^ *conda_prefix:' "${CONFIG}" | head -1 | sed 's/.*: *"\{0,1\}\([^"]*\)"\{0,1\}/\1/')
eval "$("${_CONDA_PREFIX}/condabin/conda" shell.bash hook)"
conda activate "$(dirname "$(dirname "${BENCH_PYTHON}")")"
mkdir -p "${CSM_BENCH_DIR}/logs" "${CSM_BENCH_DIR}/manifests"
if [[ "${CSM_BENCH_MODE}" == reference ]]; then
    rm -f "${CSM_BENCH_DIR}/REFERENCE_FAILED" "${CSM_BENCH_DIR}/REFERENCE_SUCCESS"
    trap 'touch "${CSM_BENCH_DIR}/REFERENCE_FAILED"' ERR
    tar -xzf "${CSM_BENCH_DIR}/tools/kallisto_linux-v0.52.0.tar.gz" -C "${CSM_BENCH_DIR}/tools"
    "${BENCH_PYTHON}" -m venv --system-site-packages "${CSM_BENCH_DIR}/tools/kb_env"
    "${CSM_BENCH_DIR}/tools/kb_env/bin/python" -m pip install --disable-pip-version-check 'kb-python==0.29.5'
    time "${BENCH_PYTHON}" scripts/build_kallisto_reference.py \
        --fasta /medpop/rgupta/Jiahao/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa \
        --gtf /medpop/rgupta/Jiahao/reference_genome/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz \
        --genome hg38 --output "${CSM_BENCH_DIR}/reference" \
        --kb "${CSM_BENCH_DIR}/tools/kb_env/bin/kb" \
        --kallisto "${CSM_BENCH_DIR}/tools/kallisto/kallisto" --threads "${NSLOTS:-16}"
    touch "${CSM_BENCH_DIR}/REFERENCE_SUCCESS"
    exit 0
fi
libraries=(C1 C2 T1 T2)
if [[ "${CSM_BENCH_MODE}" == prepare ]]; then
    library="${libraries[$((SGE_TASK_ID - 1))]}"
    rm -f "${CSM_BENCH_DIR}/manifests/${library}.prepare_FAILED"
    trap 'touch "${CSM_BENCH_DIR}/manifests/${library}.prepare_FAILED"' ERR
    "${BENCH_PYTHON}" scripts/benchmark_gex_backends.py prepare \
        --config "${CONFIG}" --work-dir "${CSM_BENCH_DIR}" --library "${library}" \
        --pairs 10000000 --threads "${NSLOTS:-16}"
elif [[ "${CSM_BENCH_MODE}" == run ]]; then
    arms=(star kallisto_unique kallisto_fractional star_directional kallisto_unique_directional kallisto_fractional_directional star_none kallisto_unique_none kallisto_fractional_none)
    library="${libraries[$(((SGE_TASK_ID - 1) % 4))]}"
    arm="${arms[$(((SGE_TASK_ID - 1) / 4))]}"
    rm -f "${CSM_BENCH_DIR}/manifests/${library}.${arm}_FAILED"
    trap 'touch "${CSM_BENCH_DIR}/manifests/${library}.${arm}_FAILED"' ERR
    if [[ "${arm}" != star* ]]; then test -f "${CSM_BENCH_DIR}/REFERENCE_SUCCESS"; fi
    "${BENCH_PYTHON}" scripts/benchmark_gex_backends.py run \
        --config "${CONFIG}" --work-dir "${CSM_BENCH_DIR}" --library "${library}" \
        --arm "${arm}" --pairs 10000000 --threads "${NSLOTS:-16}"
elif [[ "${CSM_BENCH_MODE}" == report ]]; then
    rm -f "${CSM_BENCH_DIR}/BENCHMARK_FAILED"
    trap 'touch "${CSM_BENCH_DIR}/BENCHMARK_FAILED"' ERR
    "${BENCH_PYTHON}" scripts/benchmark_gex_backends.py report \
        --config "${CONFIG}" --work-dir "${CSM_BENCH_DIR}" --pairs 10000000
else
    echo "Unknown benchmark mode: ${CSM_BENCH_MODE}" >&2
    exit 2
fi

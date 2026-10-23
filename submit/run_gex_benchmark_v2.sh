#!/bin/bash
# v2 cached Step 1a benchmark worker; use submit_gex_benchmark_v2.sh.
#$ -q broad
#$ -N csm_gex_v2_bench
#$ -j y
#$ -pe smp 16
#$ -binding linear:16
#$ -l h_vmem=8G
#$ -l h_rt=12:00:00
#$ -l os=RedHat8
set -euo pipefail
: "${CSM_PROJECT_DIR:?}"
: "${CSM_CONFIG:?}"
: "${CSM_BENCH_DIR:?}"
: "${CSM_BENCH_MODE:?init, prepare, align, assign, run, report}"
: "${CSM_BENCH_REFERENCE:?existing validated mature-plus-nascent reference}"
: "${CSM_BENCH_KALLISTO:?absolute kallisto executable}"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/csm-v2-bench-mpl-${JOB_ID:-local}-${SGE_TASK_ID:-1}"
cd "${CSM_PROJECT_DIR}"
BENCH_SOURCE_ROOT="${CSM_BENCH_SOURCE_ROOT:-${CSM_PROJECT_DIR}}"
export PYTHONPATH="${BENCH_SOURCE_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
BENCH_PYTHON="${CSM_BENCH_PYTHON:-/medpop/rgupta/Jiahao/conda/envs/splitpipe/bin/python}"
CONDA_BASE=$(sed -n 's/^ *conda_prefix: *"\{0,1\}\([^"#]*\)"\{0,1\}.*/\1/p' "${CSM_CONFIG}" | head -1)
eval "$("${CONDA_BASE%/}/condabin/conda" shell.bash hook)"
conda activate "$(dirname "$(dirname "${BENCH_PYTHON}")")"
test "$(command -v python)" = "${BENCH_PYTHON}"
mkdir -p "${CSM_BENCH_DIR}/logs" "${CSM_BENCH_DIR}/manifests"
sentinel="${CSM_BENCH_DIR}/manifests/${CSM_BENCH_MODE}.${JOB_ID:-local}.${SGE_TASK_ID:-1}"
trap 'touch "${sentinel}.FAILED"' ERR
if [[ "${CSM_BENCH_REQUIRE_PUBLICATION_REPAIR:-false}" == true ]]; then
    "${BENCH_PYTHON}" -c 'import json,sys; assert json.load(open(sys.argv[1])).get("complete") is True' \
        "${CSM_BENCH_DIR}/PUBLICATION_REPAIR_SUCCESS.json"
fi
if [[ "${CSM_BENCH_MODE}" == validate ]]; then
    export CSM_TEST_KALLISTO="${CSM_BENCH_KALLISTO}"
    export CSM_TEST_STAR="$(command -v STAR)"
    for pattern in test_benchmark_gex_v2.py test_kallisto_backend.py test_star_fractional.py test_uger_resources.py; do
        "${BENCH_PYTHON}" -m unittest discover -s tests -p "${pattern}" -v
    done
    date -u +%Y-%m-%dT%H:%M:%SZ > "${sentinel}.SUCCESS"
    exit 0
fi
args=("${CSM_BENCH_MODE}" --config "${CSM_CONFIG}" --work-dir "${CSM_BENCH_DIR}"
      --reference "${CSM_BENCH_REFERENCE}" --kallisto "${CSM_BENCH_KALLISTO}"
      --pairs "${CSM_BENCH_PAIRS:-50000000}" --threads "${NSLOTS:-16}")
if [[ "${CSM_BENCH_REQUIRE_BASELINE:-false}" == true ]]; then args+=(--require-count-baseline); fi
if [[ -n "${CSM_BENCH_ARCHIVE_ARMS:-}" ]]; then
    IFS=: read -r -a archive_arms <<< "${CSM_BENCH_ARCHIVE_ARMS}"
    for archive_arm in "${archive_arms[@]}"; do args+=(--archive-arm "${archive_arm}"); done
fi
libraries=(C1 C2 T1 T2)
backends=(star kallisto)
arms=(star_unique_parse star_unique_directional star_unique_none
      star_fractional_parse star_fractional_directional star_fractional_none
      kallisto_unique_parse kallisto_unique_directional kallisto_unique_none
      kallisto_fractional_parse kallisto_fractional_directional kallisto_fractional_none)
case "${CSM_BENCH_MODE}" in
    prepare|align|assign|run|retry-star)
        : "${SGE_TASK_ID:?array task required}"
        (( SGE_TASK_ID >= 1 ))
        case "${CSM_BENCH_MODE}" in
            prepare|retry-star) (( SGE_TASK_ID <= 4 ));;
            align|assign) (( SGE_TASK_ID <= 8 ));;
            run) (( SGE_TASK_ID <= 48 ));;
        esac
        args+=(--library "${libraries[$(((SGE_TASK_ID - 1) % 4))]}")
        if [[ "${CSM_BENCH_MODE}" == align || "${CSM_BENCH_MODE}" == assign ]]; then
            args+=(--backend "${backends[$(((SGE_TASK_ID - 1) / 4))]}")
        elif [[ "${CSM_BENCH_MODE}" == run ]]; then
            args+=(--arm "${arms[$(((SGE_TASK_ID - 1) / 4))]}")
        elif [[ "${CSM_BENCH_MODE}" == retry-star ]]; then
            args+=(--repair-note "${CSM_BENCH_REPAIR_NOTE:?}")
        fi
        ;;
    report) args+=(--verify-cache);;
    init) ;;
    accept-revision|archive-kallisto-counting) args+=(--repair-note "${CSM_BENCH_REPAIR_NOTE:?}");;
    *) echo "Unknown mode: ${CSM_BENCH_MODE}" >&2; exit 2;;
esac
"${BENCH_PYTHON}" "${BENCH_SOURCE_ROOT}/scripts/benchmark_gex_v2.py" "${args[@]}"
date -u +%Y-%m-%dT%H:%M:%SZ > "${sentinel}.SUCCESS"

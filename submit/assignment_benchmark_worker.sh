#!/bin/bash
set -euo pipefail
: "${PROJECT_DIR:?}" "${CORE_DIR:?}" "${RUN_DIR:?}" "${BENCH_STAGE:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
task="${SGE_TASK_ID:-0}"
[[ "$task" =~ ^[0-9]+$ ]] || task=0
trap 's=$?; if [[ $s == 0 ]]; then touch "$RUN_DIR/${task}._SUCCESS"; else echo "$s" > "$RUN_DIR/${task}._FAILED"; fi' EXIT
source "$HOME/my-conda/etc/profile.d/conda.sh"
conda activate splitpipe
cd "$PROJECT_DIR"
cfg="$CORE_DIR/recovered.yaml"
libs=(C1 C2 T1 T2)
lib="${libs[$(( task > 0 ? task - 1 : 0 ))]:-}"
require_stage() {
  local stage="$1" count="$2" prior i
  prior="$(awk -F '\t' -v s="$stage" '$2==s {p=$4} END {print p}' "$CORE_DIR/submissions.tsv")"
  [[ -n "$prior" ]] || { echo "Missing dependency submission: $stage" >&2; exit 1; }
  if [[ "$count" == 0 ]]; then [[ -f "$prior/0._SUCCESS" ]];
  else for ((i=1;i<=count;i++)); do [[ -f "$prior/$i._SUCCESS" ]] || return 1; done; fi
}
case "$BENCH_STAGE" in
  prepare) python -u scripts/grna_assignment_benchmark.py prepare --core "$CORE_DIR" ;;
  capture) python -u scripts/grna_assignment_benchmark.py capture --core "$CORE_DIR" --library "$lib" ;;
  barcodes)
    [[ -f "$CORE_DIR/prepared/_SUCCESS" ]]
    python -u scripts/02_extract_grna_barcodes.py --config "$cfg" --library "$lib"
    ;;
  pilot) require_stage barcodes 4; python -u scripts/03_match_grna_ibar.py --config "$cfg" --preflight-only ;;
  match|standard_match)
    if [[ "$BENCH_STAGE" == match ]]; then require_stage pilot 0; else require_stage barcodes 4; fi
    [[ "$BENCH_STAGE" != standard_match ]] || cfg="$CORE_DIR/standard.yaml"
    python -u scripts/03_match_grna_ibar.py --config "$cfg" --library "$lib" --n_processes "${NSLOTS:-1}"
    ;;
  step1c)
    require_stage match 4
    [[ -f "$CORE_DIR/recovered/gex_grna_filter/$lib/_SUCCESS" ]]
    python -u scripts/01c_process_gex_grna.py --config "$cfg" --library "$lib" \
      --input-manifest "$CORE_DIR/recovered/gex_grna_filter/$lib/rejected_pairs.complete.json" --n_processes "${NSLOTS:-1}"
    ;;
  count|count_recovered)
    require_stage match 4; require_stage standard_match 4; require_stage step1c 4
    require_stage tests 0
    python scripts/validate_step1c_manifests.py --config "$cfg"
    if [[ "$BENCH_STAGE" == count ]]; then
      python -u scripts/04_assign_grna_to_cells.py --config "$CORE_DIR/standard.yaml" --source-mode dedicated --ncpus "${NSLOTS:-1}"
    else
      python -u scripts/grna_assignment_benchmark.py validate-standard --core "$CORE_DIR"
    fi
    python -u scripts/04_assign_grna_to_cells.py --config "$cfg" --source-mode dedicated --ncpus "${NSLOTS:-1}"
    python -u scripts/04_assign_grna_to_cells.py --config "$cfg" --source-mode gex --ncpus "${NSLOTS:-1}"
    python -u scripts/04_assign_grna_to_cells.py --config "$cfg" --source-mode combined --ncpus "${NSLOTS:-1}"
    touch "$CORE_DIR/_COUNTS_SUCCESS"
    ;;
  environment)
    conda activate csm_analysis
    bash "$PROJECT_DIR/submit/setup_assignment_benchmark_env.sh"
    ;;
  export|methods|concordance|tests)
    if [[ "$BENCH_STAGE" == tests ]]; then
      python -m unittest discover -s tests -p test_step4_match_stream.py -v
      python -m unittest discover -s tests -p test_guide_design.py -v
      python -m unittest discover -s tests -p test_grna_recovery.py -v
    fi
    source "$ASSIGN_PREFIX/bin/activate"
    [[ -f "$ASSIGN_PREFIX/SEALED" ]]
    if [[ "$BENCH_STAGE" == tests ]]; then
      python -m unittest discover -s tests -p test_assignment_benchmark.py -v
      python -m unittest discover -s analysis/tests -p test_assignment_report.py -v
      /medpop/rgupta/Jiahao/conda/envs/csm_sceptre_20261116/bin/Rscript --vanilla \
        -e 'invisible(parse(file="analysis/scripts/a8b_assignment_kd.R")); cat("KD R syntax passed\n")'
    else
      [[ "$BENCH_STAGE" != methods ]] || require_stage tests 0
      python -u scripts/benchmark_guide_methods.py "$BENCH_STAGE" --core "$CORE_DIR" --task "$task"
    fi
    ;;
  *) echo "Unknown stage $BENCH_STAGE" >&2; exit 2 ;;
esac

#!/bin/bash
# Login-safe: scheduler bookkeeping only. Scientific work enters UGER.
set -euo pipefail
stage="${1:?stage}"; core="${2:?core directory}"; dependency="${3:-}"; task_range="${4:-}"
self="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJECT_DIR="$(cd "$self/.." && pwd)"
source "$PROJECT_DIR/analysis/scripts/submit_lib.sh"
case "$core" in /*) ;; *) core="$PROJECT_DIR/$core" ;; esac
export CORE_DIR="$core" BENCH_STAGE="$stage"
export ASSIGN_PREFIX="${ASSIGN_PREFIX:-/medpop/rgupta/Jiahao/venvs/csm_assign_20261116_v2}"
export RUN_DIR="$core/jobs/${stage}_$(date -u +%Y%m%dT%H%M%S)"
[[ ! -e "$RUN_DIR" ]] || csm_die 'Submission directory exists'
mkdir -p "$RUN_DIR"
resources=(-binding linear:1 -l h_vmem=24G -l h_rt=24:00:00)
array=()
case "$stage" in
 prepare|environment|export|concordance|tests) ;;
 capture) resources=(-pe smp 8 -binding linear:8 -l h_vmem=8G -l h_rt=48:00:00); array=(-t 1-4 -tc 2) ;;
 barcodes) array=(-t 1-4 -tc 4) ;;
 pilot) resources=(-pe smp 4 -binding linear:4 -l h_vmem=8G -l h_rt=12:00:00) ;;
 match|standard_match|step1c) resources=(-pe smp 8 -binding linear:8 -l h_vmem=4G -l h_rt=24:00:00); array=(-t 1-4 -tc 2) ;;
 count|count_recovered) resources=(-pe smp 4 -binding linear:4 -l h_vmem=16G -l h_rt=24:00:00) ;;
 methods) array=(-t 1-11 -tc 3) ;;
 *) csm_die "Unknown stage: $stage" ;;
esac
if [[ -n "$task_range" ]]; then
  [[ "$stage" == methods && "$task_range" =~ ^([1-9]|1[01])(-([1-9]|1[01]))?$ ]] || csm_die 'Invalid methods task range'
  array=(-t "$task_range" -tc 3)
fi
csm_require_uger
holds=(); [[ -z "$dependency" ]] || holds=(-hold_jid "$dependency")
job="$(csm_qsub -N "ga_$stage" -q broad -hard -l os=RedHat8 "${resources[@]}" \
 ${array[@]+"${array[@]}"} ${holds[@]+"${holds[@]}"} -j y -o "$RUN_DIR" \
 -v PROJECT_DIR,CORE_DIR,RUN_DIR,BENCH_STAGE,ASSIGN_PREFIX "$self/assignment_benchmark_worker.sh")"
[[ "$job" =~ ^[0-9]+ ]] || csm_die "Bad job ID: $job"
printf '%s\n' "$job" > "$RUN_DIR/job_id.txt"
printf '%s\t%s\t%s\t%s\n' "$(date -u +%FT%TZ)" "$stage" "$job" "$RUN_DIR" >> "$core/submissions.tsv"
printf '%s\n' "$job"

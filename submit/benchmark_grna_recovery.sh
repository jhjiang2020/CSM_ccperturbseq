#!/bin/bash
#$ -N v21_benchmark
#$ -q broad
#$ -hard
#$ -l os=RedHat8,h_vmem=8G,h_rt=8:00:00
#$ -pe smp 4
#$ -binding linear:4
#$ -t 1-4
#$ -j y
set -euo pipefail
: "${PROJECT_ROOT:?}"
eval "$("${HOME}/my-conda/condabin/conda" shell.bash hook)"
conda activate splitpipe
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
cd "$PROJECT_ROOT"
base="$PROJECT_ROOT/output/grna_spacer_scan_20261023/fiveprime_20261104/v2_1"
libs=(C1 C2 T1 T2)
lib="${libs[$SGE_TASK_ID-1]}"
trap 'status=$?; if ((status)); then echo "$status" > "${base}/jobs/_FAILED.${JOB_ID}.${SGE_TASK_ID}"; fi' EXIT
[[ -f "$base/pilot/_SUCCESS.$lib" ]]
python -B scripts/benchmark_grna_recovery.py --library "$lib" --output "$base" --workers "${NSLOTS:-1}"

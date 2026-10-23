#!/bin/bash
# Publish only a fully validated v2 report; never copy experimental matrices.
#$ -q broad
#$ -N csm_v2_bench_publish
#$ -j y
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_vmem=8G
#$ -l h_rt=1:00:00
#$ -l os=RedHat8
set -euo pipefail
: "${CSM_PROJECT_DIR:?}"
: "${CSM_BENCH_DIR:?}"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
sentinel="${CSM_BENCH_DIR}/manifests/publish.${JOB_ID:-local}"
trap 'touch "${sentinel}.FAILED"' ERR
"${CSM_BENCH_PYTHON:-/medpop/rgupta/Jiahao/conda/envs/splitpipe/bin/python}" - "${CSM_PROJECT_DIR}" "${CSM_BENCH_DIR}" <<'PY'
import json
from pathlib import Path
import shutil
import sys

project, run = map(Path, sys.argv[1:])
state = json.loads((run / 'BENCHMARK_SUCCESS.json').read_text())
if (state.get('complete') is not True or state.get('pairs_per_library') != 50_000_000
        or state.get('libraries') != ['C1','C2','T1','T2'] or len(state.get('arms', [])) != 12
        or state.get('configuration_pairs_per_library') != 66
        or state.get('full_cache_checksums_verified') is not True):
    raise ValueError('Refusing to publish an incomplete, unchecked or non-50M v2 experiment')
source = run / 'comparison'
names = ['stage_resources.csv', 'runtime.csv', 'library_summary.csv', 'agreement.csv',
         'umi_gains.csv', 'stage_summary.csv', 'failed_attempts.csv', 'cache_storage_baseline.csv',
         'count_staging_equality.csv', 'umi_gains.png', 'agreement.png', 'runtime.png']
names += [name for name in ('resource_peaks.csv', 'input_read_flow.csv', 'fixed_cohort_zero_counts.csv')
          if (source / name).is_file()]
accounting = source / 'scheduler_accounting_summary.json'
if accounting.is_file() and json.loads(accounting.read_text()).get('missing_accounting') == []:
    names += ['scheduler_accounting.csv', 'scheduler_coverage.csv', 'scheduler_accounting_summary.json']
for name in ['REPORT.md', *names]:
    if not (source / name).is_file():
        raise FileNotFoundError(source / name)
docs = project / 'docs'
assets_relative = Path('assets') / run.name
assets = docs / assets_relative
assets.mkdir(parents=True, exist_ok=True)
for name in names:
    temporary = assets / (name + '.tmp')
    shutil.copy2(source / name, temporary)
    temporary.replace(assets / name)
text = (source / 'REPORT.md').read_text()
for name in names:
    text = text.replace(f']({name})', f']({(assets_relative / name).as_posix()})')
text += f'\nExperiment: `{run.relative_to(project)}`. [Design and accounting definitions](GEX_BENCHMARK_V2.md).\n'
destination = docs / 'GEX_BENCHMARK_V2_REPORT.md'
temporary = destination.with_suffix('.md.tmp')
temporary.write_text(text)
temporary.replace(destination)
design = docs / 'GEX_BENCHMARK_V2.md'
text = design.read_text()
start, end = '<!-- benchmark-v2-run-status:start -->', '<!-- benchmark-v2-run-status:end -->'
if start in text and end in text:
    before, rest = text.split(start, 1)
    _, after = rest.split(end, 1)
    text = before + start + '\n\nThe 50M-pair experiment completed and passed full cache validation. See the\n[measured performance, agreement and UMI-gain report](GEX_BENCHMARK_V2_REPORT.md).\nRun provenance and all scheduler job IDs remain in `' + str(run.relative_to(project)) + '/logs/uger_submissions.jsonl`.\n\n' + end + after
    temporary = design.with_suffix('.md.tmp')
    temporary.write_text(text)
    temporary.replace(design)
PY
date -u +%Y-%m-%dT%H:%M:%SZ > "${sentinel}.SUCCESS"

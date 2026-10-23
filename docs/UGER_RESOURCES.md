# UGER resources and portability — v2.0

All resource settings are optional. A configuration with no `uger:` block keeps
the preceding workflow's CPU, memory and time defaults. The same resource
resolver serves production Steps 1–5, reference builds, and maintained QC/audit
submission drivers. STAR and kallisto use the same Step 1a reservation.

Use `submit/submit_pipeline.sh` for the main workflow. Use
`submit/submit_job.sh` for an individual worker. Direct `qsub worker.sh` remains
supported with the worker's historical `#$` defaults, but **does not read YAML
resource settings**: scheduler requests must be resolved before a job starts.
Editing YAML after submission does not change a queued or running job.

## YAML interface

Merge this into your existing configuration's single `uger:` block; do not add
a second block with the same key:

```yaml
uger:
  queue: broad
  parallel_environment: smp
  binding: linear
  # os: RedHat8
  # project: your_uger_account
  # extra_resources:
  #   scratch: "100G"  # Only if this complex exists at your site.
  defaults:
    max_concurrent: 2
  step1a:
    slots: 16
    memory_per_slot_gb: 8
    runtime_hours: 48
  step1b:
    slots: 4
    memory_per_slot_gb: 24
    runtime_hours: 8
  step3:
    slots: 12
    memory_per_slot_gb: 4
    runtime_hours: 24
    max_concurrent: 3
  step5:
    slots: 1
    memory_per_slot_gb: 64
    runtime_hours: 8
```

This example increases Step 1b's total reservation to 96 GB and Step 5's to
64 GB. It limits arrays to two tasks at once, except Step 3, which allows three.
It is an example, not a recommendation that every run needs larger requests.

| Key | Meaning and accepted values |
| --- | --- |
| `slots` | CPU slots, positive integer up to 256. The resolver requests matching `-pe` and CPU binding. |
| `memory_per_slot_gb` | Positive number up to 4096; rendered as `h_vmem=<value>G`. Fractions such as `2.5` are accepted. |
| `runtime_hours` | Positive hours up to 744, including fractions such as `0.5`; rendered as `h_rt=HH:MM:SS`. Minimum one second. |
| `max_concurrent` | Positive integer or `null`; emits `-tc` only for an array. `null` removes the resolver's concurrency cap. |
| `queue` | UGER queue or queue/host selection, such as `broad` or `queue@host`; a comma-separated queue list is accepted. |
| `parallel_environment` | Name of the site's shared-memory PE; default `smp`. |
| `binding` | `linear` (default, binds the requested number of cores) or `none` (explicitly disables binding). |
| `os` | Optional OS complex, e.g. `RedHat8`. |
| `project` | Optional scheduler accounting project, passed with `-P`; unrelated to `paths.project_dir`. |
| `extra_resources` | Mapping of other site-defined `-l` complexes to scalar values. Use the named fields for `h_vmem`, `h_rt` and `os`. |

Scheduler fields (`queue`, `parallel_environment`, `binding`, `os`, `project`,
`extra_resources`) may appear directly under `uger:` or inside a stage.
All fields may appear under `uger.defaults` or a stage. Stage values take
precedence over `defaults`, which take precedence over common scheduler values
and the built-in stage defaults. A stage's `extra_resources` mapping replaces
the common mapping. Null queue/OS/project values omit that command-line flag;
standalone worker directives or site defaults may still supply it.

Unknown profiles/fields, malformed numbers, nonpositive requests and unsafe
scheduler tokens fail before the first main-workflow submission. Resource
settings never change biological thresholds or UMI correction.

## Production defaults and where CPUs are used

| Profile | Slots | GB/slot | Total GB/task | Hours | Parallel work |
| --- | ---: | ---: | ---: | ---: | --- |
| `step1a` | 16 | 8 | 128 | 48 | Parse preprocessing, STAR/kallisto alignment and molecule counting in successive phases. |
| `step1b` | 4 | 16 | 64 | 4 | Parse combine uses allocated threads; native kallisto combine is largely sparse I/O/aggregation. |
| `step1c` | 1 | 4 | 4 | 12 | Guide matching consumes `NSLOTS`; one array task per library. |
| `step1c_aggregate` | 4 | 32 | 128 | 8 | Optional guide evidence aggregation/assignment consumes `NSLOTS`. |
| `step2` | 1 | 4 | 4 | 12 | Current dedicated-guide barcode FASTQ loop is serial. More slots alone will not speed it up. |
| `step3` | 8 | 4 | 32 | 12 | Guide/iBAR matching workers consume `NSLOTS`. |
| `step4` | 4 | 32 | 128 | 4 | Guide assignment/rarefaction consumes `NSLOTS`. |
| `step5` | 1 | 32 | 32 | 4 | AnnData integration is principally a serial matrix/I/O stage. Increase RAM before requesting extra cores. |

Array stages are 1a, 1c, 2 and 3. Other production profiles submit a single job.
`max_concurrent` limits tasks **within each array**, not the total number of
jobs across every stage or every user. Peak capacity per array is
`max_concurrent × slots × memory_per_slot_gb`, plus overlapping jobs.

The worker's allocation (`NSLOTS`) is authoritative for its parallel tools.
`gex_alignment.kallisto.count_threads: 0` uses the allocated Step 1a cores;
a positive value caps native counting, and `1` forces serial counting. It does
not reserve additional cores. Parse preprocessing currently caps its own workers
at eight. Kallisto alignment and native molecule workers run in sequence.
Numerical thread pools are capped in multiprocessing workers to avoid nested
oversubscription.

`splitpipe.threads`, `splitpipe.memory_gb` and
`barcode_extraction.n_processes` are retained legacy config fields; they do not
reserve UGER resources. Step 2's current loop does not use that process count.
On another scheduler, `uger:` does not allocate anything: set `NSLOTS` from that
scheduler's actual CPU reservation and use the [worker recipe](RUNNING.md#another-scheduler-or-an-existing-linux-allocation).

## Reference and standalone QC workers

The `kallisto_benchmark` profile now controls the cached twelve-choice v2 GEX
experiment. Launch it with `submit/submit_gex_benchmark_v2.sh`; see the
[50M-pair benchmark design](GEX_BENCHMARK_V2.md) and the completed
[measured CPU, memory and disk report](GEX_BENCHMARK_V2_REPORT.md). It requests 128 GiB total by
default (16 × 8 GiB), with two concurrent array tasks unless overridden.

| Profile | Slots | GB/slot | Hours | Worker |
| --- | ---: | ---: | ---: | --- |
| `kallisto_reference` | 8 | 8 | 12 | `submit/submit_kallisto_reference.sh` |
| `kallisto_benchmark` | 16 | 8 | 12 | `submit/run_gex_benchmark_v2.sh` |
| `shared_kallisto_environment` | 4 | 8 | 3 | `submit/submit_shared_kallisto_environment.sh` |
| `q1` | 8 | 4 | 12 | `submit/submit_q1_guideseq.sh` |
| `gex_contamination` | 8 | 8 | 48 | Raw GEX carry-over QC array |
| `gex_contamination_alignment` | 8 | 16 | 12 | Controlled STAR QC array |
| `gex_contamination_aggregate` | 1 | 4 | 1 | GEX carry-over QC aggregation |
| `c2_case_study` | 4 | 4 | 4 | `submit/submit_c2_gex_grna_case_study.sh` |

Reference, benchmark and environment-build profiles retain their RedHat8 default.
Q2–Q4 notebook/analysis code does not have production `qsub` workers in this
repository. Allocate its resources through the analysis environment/scheduler;
these profiles do not silently control notebook kernels or the `analysis/`
submodule.

To submit one worker with YAML resources, use absolute paths in worker arguments
and exported scheduler variables:

```bash
CONFIG="$(realpath my_project.v2.yaml)"
PROJECT_DIR="$(pwd -P)"
bash submit/submit_job.sh --config "${CONFIG}" --profile kallisto_reference \
  --dry-run -- -v "CSM_CONFIG=${CONFIG},CSM_PROJECT_DIR=${PROJECT_DIR}" \
  -o "${PROJECT_DIR}/logs/kallisto-ref.\$JOB_ID.log" \
  "${PROJECT_DIR}/submit/submit_kallisto_reference.sh" \
  --fasta /path/to/genome.fa --gtf /path/to/annotation.gtf \
  --genome hg38 --output "${PROJECT_DIR}/references/kallisto_GRCh38" \
  --kb /path/to/tools/bin/kb --kallisto /path/to/tools/bin/kallisto
```

Create the log directory and remove `--dry-run` to submit. For Q1, select
`--profile q1` and pass `-v CSM_QC_CONFIG=/absolute/qc.yaml` plus its worker.
Select the correct profile explicitly: the generic launcher does not guess it
from the script filename. It accepts `--record /absolute/file.jsonl` through the
Python resolver CLI; otherwise the generic shell launcher records under the
resource config's `paths.log_dir`, or a sibling `logs/` directory.

## Standalone QC driver

`submit/submit_qc_gex_grna_contamination_pipeline.sh` uses the same resolver
with profiles `gex_contamination`, `gex_contamination_alignment` and
`gex_contamination_aggregate`. Pass `--config /absolute/path/to/my_qc.yaml`;
see the [QC examples](../examples/qc/README.md). Default raw/alignment array
concurrency is two. Set `CSM_UGER_PYTHON` when Python with PyYAML is not on PATH.

Older experiment-specific resource profiles remain readable for saved YAMLs;
their one-use drivers are local archive material, not release entry points.

## Inspect requests and tune a rerun

```bash
python scripts/uger_resources.py --config "${CONFIG}" --validate
python scripts/uger_resources.py --config "${CONFIG}" --step step5
bash submit/submit_pipeline.sh --config "${CONFIG}" --dry-run
qstat -j JOB_ID
qacct -j JOB_ID
```

The resolver only reads small config/text files. Submission bookkeeping uses a
30-second timeout when GNU `timeout` is available, and each `qsub` call has a
25-second timeout. It imports no numerical libraries or assay data.

Successful submissions first save an atomic per-job record under
`uger_submissions.records/<job-id>.json`, then append under a file lock to
`uger_submissions.jsonl` in the configured log folder (multi-stage QC drivers use
their own `LOG_DIR`). If shared-filesystem locking or the JSONL append fails,
the resolver returns the submitted job ID with a warning and preserves the
per-job record. Check that record before retrying a submission. Each entry records UTC
time, exact job ID, command arguments and resolved resource values including
total GB. Main-workflow dependencies use the exact returned job IDs so two runs
do not share holds just because their job names match. A scheduler hold still
means ordering, not success; every downstream stage must validate its inputs.
A failure midway through submission leaves the already-recorded jobs visible;
inspect them before resubmitting.

For a memory failure, increase that profile's `memory_per_slot_gb`; total memory
changes when `slots` changes too. For a wall-time failure, increase
`runtime_hours`. For CPU-bound matching/alignment/counting, increase `slots`
only if the stage uses parallel workers. For storage contention, reduce
`max_concurrent` before adding more cores. Use measured resident memory and
phase timings; `qacct maxvmem` is virtual memory, not RSS. Larger requests can
increase queue wait, and more workers can increase counting RAM/scratch I/O.
Keep biological config, read sampling and output policy fixed when comparing
resource changes.

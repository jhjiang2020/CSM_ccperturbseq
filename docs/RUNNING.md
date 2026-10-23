# Execute, monitor and resume the v2.0 workflow

Use a complete configuration prepared with the [quick start](QUICKSTART.md).
The workflow is Steps 0–5; Step 0 and genome-reference building are separate
preparations, not jobs submitted by the Steps 1–5 master launcher.

## Broad UGER

On the login node, initialize scheduler commands and use an absolute config:

```bash
cd CSM_ccperturbseq
PROJECT_DIR="$(pwd -P)"
CONFIG="$(realpath my_project.v2.yaml)"
source /broad/software/scripts/useuse
reuse -q UGER
bash submit/submit_pipeline.sh --config "${CONFIG}" --dry-run
bash submit/submit_pipeline.sh --config "${CONFIG}"
```

Review library count, script paths, CPU binding, memory **per slot**, wall time,
config path and log destinations in the dry run. The launcher bootstraps only
lightweight YAML/resource parsing on the login node. Its workers initialize
the configured production Conda environment inside each compute job. Use simple absolute paths. The launcher passes command arguments directly;
config/project paths cannot contain commas/newlines (`qsub -v`) or apostrophes
(legacy worker config readers).

Select a subset with the same selector in the dry run and actual submission:

```bash
bash submit/submit_pipeline.sh --config "${CONFIG}" --step 1 --dry-run
bash submit/submit_pipeline.sh --config "${CONFIG}" --step 1

bash submit/submit_pipeline.sh --config "${CONFIG}" --step 2,3,4 --dry-run
bash submit/submit_pipeline.sh --config "${CONFIG}" --step 2,3,4

bash submit/submit_pipeline.sh --config "${CONFIG}" --step 4,5 --dry-run
bash submit/submit_pipeline.sh --config "${CONFIG}" --step 4,5
```

Omitted prerequisites must already be complete under the configured output
root. `--step 1` means the Step 1a array, Step 1b, and optional Step 1c array and
aggregate. There is no `--step 1a` selector. Steps 2–4 remain the dedicated-guide
path regardless of GEX backend.

The actual dependency graph is in the [detailed workflow](DETAILED_WORKFLOW.md#uger-job-dependency-chain).
Step 2 depends on Step 1a when both are submitted; Step 1b and later guide work
can proceed along separate branches after their prerequisites. Optional Step 1c
feeds Step 4 only when source merging is enabled.

Monitor with short, fresh queries:

```bash
qstat -u "$USER"
qacct -j JOB_ID
tail -n 60 "${PROJECT_DIR}/logs/step1a.JOB_ID.TASK_ID.log"
```

Use the configured `paths.log_dir` if it differs. Do not run Python/H5AD/BAM
probes on the login node. Inspect those in an allocation. Dependency holds mean
that jobs finished, not that all prerequisite jobs succeeded; check exit status,
logs and completion validation. Submitted job IDs and resolved resources are recorded in
`paths.log_dir/uger_submissions.jsonl`; dependencies use those exact IDs, so
independent runs can share readable job names.

## Another scheduler or an existing Linux allocation

The master launcher is Broad-specific. Its worker scripts perform the stage
work without invoking `qsub`; their `#$` directives are comments when run by
Bash. This lets another scheduler supply resources and dependencies without
reimplementing the biological workflow.

The following is a **sequential recipe inside an already allocated Linux job**.
It does not allocate resources or submit jobs by itself. Use enough RAM/time for
the largest stage, or schedule the individual workers separately using the
resource/dependency tables. Bash 4+ and GNU utilities are expected.

```bash
set -euo pipefail
cd CSM_ccperturbseq
export CSM_PROJECT_DIR="$(pwd -P)"
export CSM_CONFIG="$(realpath my_project.v2.yaml)"
# Activate your configured production environment for the small setup probes.
conda activate splitpipe
# On Slurm, use --cpus-per-task and set NSLOTS from that allocation.
export NSLOTS="${SLURM_CPUS_PER_TASK:-16}"
export JOB_ID="${SLURM_JOB_ID:-manual-v2}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1

N_LIBRARIES="$(python -c 'import sys,yaml; print(len(yaml.safe_load(open(sys.argv[1]))["libraries"]))' "${CSM_CONFIG}")"
STEP1C="$(python -c 'import sys,yaml; c=yaml.safe_load(open(sys.argv[1])); print("yes" if c.get("gex_decontamination",{}).get("step1c",{}).get("enabled",False) else "no")' "${CSM_CONFIG}")"

for task in $(seq 1 "${N_LIBRARIES}"); do
  SGE_TASK_ID="${task}" bash submit/submit_step1a_gex_array.sh
done
bash submit/submit_step1b_gex_combine.sh

if [[ "${STEP1C}" == yes ]]; then
  for task in $(seq 1 "${N_LIBRARIES}"); do
    SGE_TASK_ID="${task}" bash submit/submit_step1c_gex_grna_array.sh
  done
  bash submit/submit_step1c_gex_grna_aggregate.sh
fi

for task in $(seq 1 "${N_LIBRARIES}"); do
  SGE_TASK_ID="${task}" bash submit/submit_step2_barcodes.sh
done
for task in $(seq 1 "${N_LIBRARIES}"); do
  SGE_TASK_ID="${task}" bash submit/submit_step3_match.sh
done
bash submit/submit_step4_assign.sh
bash submit/submit_step5_anndata.sh
```

Replace the initial activation if you use an absolute/shared prefix. Every
worker subsequently activates `environment.conda_env` from your YAML. Keep
`JOB_ID` defined: some array workers use it without a fallback.
`SGE_TASK_ID` is a **1-based index into `libraries`**, not a sample name. For a
scheduler array, map its task ID to that index explicitly and set `NSLOTS` to
the number of CPUs allocated to that task. Do not run the loops on a login node.

`scripts/01_run_gex_splitpipe.sh` is a deprecated compatibility entry point that
executes the UGER master launcher. It is not the portable sequential runner.
Invoking `run_kallisto_gex.py` alone also does not replace Step 1a: it expects
completed Parse preprocessing and its metadata/evidence contracts.

## Resources

Defaults in the current workers/master are:

| Stage | CPUs | UGER memory per slot | Total at default | Wall time | Shape |
| --- | ---: | ---: | ---: | --- | --- |
| 1a: preprocessing + GEX backend | 16 | 8 GB | 128 GB | 48 h | One task per library |
| 1b: GEX combine | 4 | 16 GB | 64 GB | 4 h | One job |
| 1c: optional guide match | 1 | 4 GB | 4 GB | 12 h | One task per library |
| 1c: optional aggregate | 4 | 32 GB | 128 GB | 8 h | One job |
| 2: guide barcode extraction | 1 | 4 GB | 4 GB | 12 h | One task per library |
| 3: guide/iBAR match | 8 | 4 GB | 32 GB | 12 h | One task per library |
| 4: guide assignment | 4 | 32 GB | 128 GB | 4 h | One job |
| 5: AnnData integration | 1 | 32 GB | 32 GB | 4 h | One job |

Every profile in this table is tunable through optional `uger.<profile>`
settings. See [the complete resource interface](UGER_RESOURCES.md) for defaults,
array concurrency, queue/OS/account settings, QC profiles and worked examples.
Use the master or generic `submit_job.sh` launcher to apply YAML resources.
Direct `qsub` workers use fallback directives. Kallisto alignment and molecule
counting use allocated cores in successive phases; `count_threads` can cap
native counting. Resource demand depends on read depth and unique molecules.

## Outputs and success checks

Paths below are relative to `paths.output_dir`:

| Product | Interpretation |
| --- | --- |
| `gex_<library>/process/` | Parse preprocessing identity, per-library count evidence and backend QC. |
| `gex_<library>/<sample>/DGE_filtered/` | Called-cell expression product; kallisto also writes `DGE_unfiltered`. |
| `gex_combined/<sample>/DGE_filtered/` | Step 1b product read by Step 5. |
| `grna_barcodes/` | Step 2 dedicated-guide read/barcode TSVs. |
| `grna_matches/` | Step 3 matches, preflight diagnostics and optional failed reads. |
| `grna_assignment/` | Dedicated or selected combined guide assignment/matrix products. |
| `gex_grna_filter/<library>/` | Enabled filter's completion manifest, counts and corrected-barcode evidence. |
| `gex_grna_identity/` | Default optional Step 1c namespace; configurable. |
| `combined_adata.h5ad` | Step 5 expression plus aligned guide counts/annotations. |

Kallisto's `process/kallisto.complete.json` and combined
`gex_combined/kallisto.complete.json` seal settings, source and output identities.
`process/kallisto_qc.json` records molecule/mapping counts, cell calls, phase
timings and `counting_execution` workers/shards. Its assignment CSV retains
molecule read support and optional fractional gene weights. Genomic BAMs and
exonic fractions are not fabricated for this backend. H5AD provenance records
the feature level and counting policy.

Unique STAR runs have `process/gex_umi_policy.json` and a combined policy marker.
Fractional STAR uses `process/star.complete.json`, `process/star_qc.json` and
`gex_combined/star.complete.json`. Its complete pre-molinfo BAM is retained and
bound to provenance; it uses the shared sparse combine.
[Counting choices and BAM requirements](GEX_COUNTING.md).
Legacy unmarked STAR outputs are accepted only with unique genes and Parse-default correction.
The filter has a separate `rejected_pairs.complete.json`; Step 1c also has its
own validated completion artifacts. Keep these with the output files.

Run the common filter/backend handoff validator in a compute allocation after
Step 1a, before manual combination:

```bash
python scripts/validate_gex_filter_manifests.py --config "${CSM_CONFIG}"
```

Step 1b already runs it, and Step 5 validates the combined expression contract.
For unique STAR, `python scripts/gex_alignment.py --config "${CSM_CONFIG}"
--validate-star-policies` additionally checks per-library policy/output identity.
Do not create marker files manually to bypass a failed validation.

## Resume and change parameters

| Change | Required work |
| --- | --- |
| STAR ↔ kallisto, GEX UMI flavor, gene/transcript policy or expression reference | New output namespace; rerun Step 1a/1b and Step 5. Preserve/recreate compatible guide products in the new namespace. |
| Guide reference/construct interpretation or GEX removal policy | Regenerate affected Step 1 filtering/evidence and guide Steps 2–4 as appropriate, then Step 5. Reference hashes intentionally reject stale evidence. |
| Step 3 matching settings | Rerun 3, 4 and 5; include Step 1c if its matching/acceptance interpretation changed. |
| Step 4 guide threshold or gRNA UMI settings | Rerun 4 and 5 from complete compatible read-match evidence. |
| Enable source merging | Require complete Step 1c plus dedicated Step 3; rebuild 4 and 5 with both source options enabled. |
| Failed Step 1a task | Inspect the error and its files; preserve the failed attempt for diagnosis before retrying in a clean destination when state is uncertain. |

Avoid editing backend/UMI settings on an old output root and merely rerunning
Step 5. The manifests are designed to reject that relabeling. See
[migration details](MIGRATING_TO_V2.md) for historical STAR BAM reuse and count
comparisons. For fresh comparisons, hold FASTQ sampling, preprocessing,
reference annotation and cell cohorts explicit; timing should separate shared
preprocessing, alignment, counting and publishing.

# CSM ccPerturbSeq Pipeline

This pipeline processes Parse Biosciences GEX FASTQs and dedicated
CROPseq-multi guide libraries into analysis-ready AnnData. It supports
single- and dual-guide/iBAR designs, any number of sublibraries, and optional
v1.4 removal of guide products that escaped GEX/gRNA size selection.

For architecture, see [pipeline_overview.md](pipeline_overview.md). Detailed
logic is in [docs/DETAILED_WORKFLOW.md](docs/DETAILED_WORKFLOW.md); v1.4 is
covered by the [operating guide](docs/GEX_DECONTAMINATION_V1_4.md) and
[QC guide](docs/GEX_GRNA_CONTAMINATION_QC.md).

v1.5 adds early iBAR reference diagnostics and Step 3 scaling controls. See the
[v1.5 release notes](RELEASE_NOTES.md#v15--ibar-preflight-and-step-3-scaling-2026-10-14)
and [preflight/scaling guide](docs/SCALABILITY_AND_PREFLIGHT.md).

## Quick start

### 1. Install the production environment

This repository does not provide the Parse split-pipe source. Obtain split-pipe
1.6.4 directly from Parse Biosciences and extract it beside this repository's
scripts under the expected directory name:

```text
CSM_ccperturbseq/
├── ParseBiosciences-Pipeline.1.6.4/
│   ├── split-pipe
│   ├── splitpipe/
│   ├── setup.py
│   ├── install_dependencies_conda.sh
│   └── install_dependencies.sh
├── scripts/
└── submit/
```

An example Broad project tree, including a locally obtained copy of the Parse
distribution, is available at
`/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq`. Access to that example does not
replace obtaining and accepting the Parse software through Parse Biosciences
for a new installation.

Parse split-pipe 1.6.4 requires Python 3.12.8 or newer. From the project root:

```bash
conda create -n splitpipe python=3.12.8
conda activate splitpipe
cd ParseBiosciences-Pipeline.1.6.4
bash install_dependencies_conda.sh --install --yes
python -m pip install .
cd ..
python -m pip install -r requirements_extra.txt
```

`requirements_extra.txt` installs the required C-backed `edlib` matcher for
v1.4. All production pipeline, decontamination, controlled-alignment, and QC
jobs run in `splitpipe`; STAR and `pysam` must be available there. To update and
verify an existing production environment:

```bash
conda activate splitpipe
python -m pip install -r requirements_extra.txt
python -c 'import edlib, pandas, pysam, splitpipe, yaml'
STAR --version
```

`environment_qc.yml` is optional and creates the separate `csm_analysis`
environment for downstream analysis modules and notebooks. Production wrappers
never activate it. On Broad, command-line `samtools` is also available after
`use Samtools`; this optional module does not replace the production Python
`pysam` dependency.

Configure the actual Conda installation rather than relying on a personal
fallback. Production QC reads the same settings:

```yaml
environment:
  conda_prefix: "/path/reported/by/conda-info-base"
  conda_env: "splitpipe"
```

### Optional installation-free shared environment on Broad

Members of `rguptalab` may reuse the published v1.4 production environment.
Point both fields to its absolute prefix; do not use a personal `/home/unix/`
path or the parent `envs/` directory:

```yaml
environment:
  conda_prefix: "/medpop/rgupta/Jiahao/conda/envs/splitpipe-share"
  conda_env: "/medpop/rgupta/Jiahao/conda/envs/splitpipe-share"
```

The Step 1-5 and production QC wrappers activate that prefix directly. Before
the first submission, run this read-only check:

```bash
set -euo pipefail
SHARED_ENV=/medpop/rgupta/Jiahao/conda/envs/splitpipe-share
CONDA="${SHARED_ENV}/condabin/conda"

printf 'shared environment:\n'
test -r "${SHARED_ENV}/csm_shared_environment.complete.json"
printf '  completion manifest: OK\n'
"${CONDA}" run --prefix "${SHARED_ENV}" \
  python -c 'import edlib, pandas, pysam, splitpipe, yaml' >/dev/null
printf '  production imports: OK\n'
"${CONDA}" run --prefix "${SHARED_ENV}" python -m pip check >/dev/null
printf '  dependencies: OK\n'
test -z "$(find "${SHARED_ENV}" -xdev -not -user jjiahao -print -quit)"
printf '  ownership: OK\n'
test -z "$(find "${SHARED_ENV}" -xdev -not -type l -perm /022 -print -quit)"
printf '  permissions: OK\n'

# Expected result:
# shared environment:
#   completion manifest: OK
#   production imports: OK
#   dependencies: OK
#   ownership: OK
#   permissions: OK
```

The prefix is owner-write and `rguptalab`-read/execute; its parent directories
use sticky group-sharing mode `1770`. Consumers must configure writable project,
output, log, and scratch paths and must not change the shared prefix. Maintainers
republish the completion manifest after updates. Never store credentials or
private configuration in the environment.

The scripts under `submit/` are written specifically for the Broad Institute's
UGER environment: they use `qsub`, the `broad` queue, UGER job arrays and
dependencies, and Broad linear core binding. On another scheduler, retain the
pipeline scripts but translate the submission layer and resource directives.

### 2. Configure inputs

Copy the closest public configuration and edit the copy:

```bash
cp config_dual.yaml my_project.yaml       # dual-guide/iBAR
# or: cp config_single.yaml my_project.yaml
```

Set `paths.project_dir`, `fastq_dir`, `splitpipe_dir`, `genome_dir`,
`output_dir`, `log_dir`, and `grna_reference`. Define every sublibrary under
`libraries`; array dimensions are derived from this list. Prefix discovery is
token-bounded, and explicit `gex_fastq_pairs` can be used for unusual names or
multiple lanes.

For Step 1a, the defaults are deliberately conservative:

```yaml
uger:
  step1a:
    slots: 16
    memory_per_slot_gb: 8
    runtime_hours: 48
```

Large or multi-lane GEX inputs should retain the 48-hour limit. More than 16
slots showed no benefit in the release benchmark, but 16 is not a hardcoded
pipeline ceiling. `submit_pipeline.sh` validates `uger.step1a.slots`, requests
the same number with `-pe smp` and `-binding linear`, and Step 1a passes the
allocated `NSLOTS` to split-pipe as `--nthreads`. Users can raise or lower
`slots` in YAML and inspect the resolved request with `--dry-run`.

The reason to keep 16 as the default is measured diminishing return: Parse
1.6.4 internally caps its preprocessing/barcode-correction worker pool at
eight (`pre_num_multiproc=8`), and the full-C2 run averaged about 5.2 cores.
STAR, sorting, and later split-pipe work may use more of the allocation, so a
different dataset can justify a measured override. More slots also reserve
more total memory because `memory_per_slot_gb` is per slot. The legacy
`splitpipe.threads` config field does not control the current UGER Step 1a
wrapper; `uger.step1a.slots` is authoritative.

Step 3 now checks a bounded read sample for possible iBAR strand, length, and
reference misspecification before starting its matching workers. Existing
configs enable this automatically; warnings and sample counts appear in the
job log and `grna_matches/<library>_preflight.yaml`. A standalone raw-R1 check
can run before any production processing with `--preflight-only
--preflight-fastq`. Optional `uger.step3` overrides use the same resource keys
as Step 1a; omitted values retain 8 slots, 4 GB per slot, and 12 hours. See
[preflight usage and the pipeline scaling review](docs/SCALABILITY_AND_PREFLIGHT.md)
for early checks, resource examples, and the limits of adding more cores.

### 3. Initialize UGER and dry-run

On a Broad login node, set the project and config paths once, initialize UGER,
and inspect the complete submission before running it:

```bash
PROJECT_DIR=/medpop/rgupta/Jiahao/path/to/CSM_ccperturbseq
cd "${PROJECT_DIR}"
CONFIG="$(realpath my_project.yaml)"

source /broad/software/scripts/useuse
reuse -q UGER || true
command -v qsub

bash submit/submit_pipeline.sh --config "${CONFIG}" --dry-run
```

Heavy commands, including split-pipe, must run through `qsub`. The dry run
prints the exact resource, binding, dependency, config, and log arguments.
The config passed through `qsub` must be an absolute HPC path: use `realpath`
as above or assign the full `/medpop/rgupta/Jiahao/...` path. A bare relative
`my_project.yaml` will not resolve reliably on a compute node. Run the dry-run
again whenever the config or requested step set changes.

### 4. Submit and monitor

To submit Steps 1–5 as one dependency graph:

```bash
bash submit/submit_pipeline.sh --config "${CONFIG}" --dry-run
bash submit/submit_pipeline.sh --config "${CONFIG}"
```

To submit only one step or a comma-separated set, use the same selector for the
dry run and production submission. Prerequisite outputs for omitted steps must
already be complete.

```bash
# Step 1a/1b and optional Step 1c only
bash submit/submit_pipeline.sh --config "${CONFIG}" --step 1 --dry-run
bash submit/submit_pipeline.sh --config "${CONFIG}" --step 1

# Dedicated-guide extraction, matching, and assignment only
bash submit/submit_pipeline.sh --config "${CONFIG}" --step 2,3,4 --dry-run
bash submit/submit_pipeline.sh --config "${CONFIG}" --step 2,3,4

# Rebuild assignment and final AnnData from completed upstream products
bash submit/submit_pipeline.sh --config "${CONFIG}" --step 4,5 --dry-run
bash submit/submit_pipeline.sh --config "${CONFIG}" --step 4,5
```

Monitor and inspect accounting separately:

```bash
qstat -u "$USER"
qacct -j JOB_ID
tail -f /path/to/logs/step1a_gex.JOB_ID.TASK_ID.log
```

Step 1 submits per-library Step 1a followed by Step 1b. Optional Step 1c is
submitted only when enabled. Steps 2–4 process the dedicated guide arm; Step 5
joins guide identities to GEX AnnData.

### 5. Find the main outputs

```text
output/
├── gex_<library>/                         # per-library split-pipe output
├── gex_combined/<sample>/                 # Step 1b combined GEX products
├── grna_barcodes/                         # dedicated guide Step 2
├── grna_matches/                          # dedicated guide Step 3
├── grna_assignment/                       # dedicated/combined Step 4 matrices
├── gex_grna_filter/<library>/             # v1.4 manifest and Step 1a evidence
├── gex_grna_identity/                     # optional Step 1c identities
└── combined_adata.h5ad                    # Step 5 analysis-ready AnnData
```

The Step 1a schema-4 completion manifest records the synchronized rejected
R1/R2 hashes and counts, ordered read-ID digest, exact FASTQ identities,
reference/settings hashes, and a compact valid-positive barcode-evidence TSV.
For historical aligned-BAM reuse, the separate schema-2 harmonized manifest
also proves that source alignment statistics match the source BAM, records the
unique/multimap removals, and installs the adjusted counts required by Parse
`--mode mol`. Invalid-barcode and unaligned filter positives may be absent from
an aligned-only source BAM; validation requires exact accounting of the
positive QNAMEs that are present and removed rather than requiring all raw
positive IDs to appear in that BAM.

Fresh v1.4 setup and the complete decontamination configuration now live in the
[v1.4 GEX decontamination guide](docs/GEX_DECONTAMINATION_V1_4.md#fresh-v14-fastq-start).
For short-read experiments, use the consolidated
[100-cycle scenario guide](docs/ILLUMINA_100_CYCLE_SEQUENCING.md#choose-the-scenario)
instead of copying a long-read configuration.

## v1.4 migration after a successful v1.3-or-earlier Step 1

Do not edit or overwrite the old run. The v1.4 runtime overlay, Step 1c, and
manifest checks require updated code; changing only the old YAML is not
sufficient.

### Update the v1.3 code and environment

The preferred route is to preserve the old config, update a clean project
checkout to the v1.4 release, and copy the config under a new name:

```bash
PROJECT_DIR=/medpop/rgupta/Jiahao/path/to/CSM_ccperturbseq
cd "${PROJECT_DIR}"
git status --short
git pull --ff-only

cp /medpop/rgupta/Jiahao/path/to/old_project.yaml \
   /medpop/rgupta/Jiahao/path/to/old_project.v1_4.yaml
conda activate splitpipe
python -m pip install -r requirements_extra.txt
python -c 'import edlib, pandas, pysam, splitpipe, yaml'
```

If `git pull` is not possible, replace the v1.3 Step 1 runtime as one unit; do
not copy only the filter module. The minimum v1.4 set is:

```text
requirements_extra.txt
patches/parsebio_splitpipe_1.6.4_gex_decontamination.patch
submit/submit_pipeline.sh
submit/submit_step1a_gex_array.sh
submit/submit_step1b_gex_combine.sh
submit/submit_step1c_gex_grna_array.sh
submit/submit_step1c_gex_grna_aggregate.sh
scripts/01_run_gex_splitpipe.sh
scripts/01c_process_gex_grna.py
scripts/03_match_grna_ibar.py
scripts/04_assign_grna_to_cells.py
scripts/gex_grna_bam_filter.py
scripts/gex_grna_filter.py
scripts/gex_grna_preprocess_hook.py
scripts/splitpipe_gex_wrapper.py
scripts/splitpipe_wrapper.py
scripts/uger_resources.py
scripts/validate_gex_filter_manifests.py
scripts/validate_step1c_manifests.py
scripts/write_fastq_input_identity.py
```

Also replace `submit/submit_step4_assign.sh`,
`submit/submit_step5_anndata.sh`, and `scripts/05_build_anndata.py` before using
the optional GEX/dedicated-guide merge. Replace
`submit/submit_step2_barcodes.sh`, `submit/submit_step3_match.sh`, and
`scripts/02_extract_grna_barcodes.py` as well if the dedicated arm will be
rerun. A full release checkout is safer than maintaining this manual file set.

### Migrate the v1.3 configuration

In the copied YAML, verify `environment.conda_prefix` and
`environment.conda_env`. Use the named `splitpipe` environment for a
conventional Conda installation or the absolute shared-prefix configuration in
the quick start above. Keep the existing `libraries`, FASTQ paths, Parse
kit/chemistry, guide structure, and guide reference, but update:

- `paths.project_dir` and `paths.splitpipe_dir` to the v1.4 checkout and its
  pristine `ParseBiosciences-Pipeline.1.6.4` tree;
- `paths.output_dir` and `paths.log_dir` to new v1.4 locations;
- `uger.step1a` with the new bound resource controls;
- the complete `gex_decontamination` section from `config_dual.yaml` or
  `config_single.yaml`, using the canonical `splitpipe_pre_router` integration;
- `existing_alignment.source_output_dir_template` to the old read-only Step 1
  output, and its sandbox template to the new output root.

The essential historical-BAM settings are:

```yaml
paths:
  project_dir: "/medpop/rgupta/Jiahao/path/to/CSM_ccperturbseq"
  splitpipe_dir: "/medpop/rgupta/Jiahao/path/to/CSM_ccperturbseq/ParseBiosciences-Pipeline.1.6.4"
  output_dir: "/medpop/rgupta/Jiahao/path/to/output_v1_4"
  log_dir: "/medpop/rgupta/Jiahao/path/to/logs_v1_4"

uger:
  step1a:
    slots: 16
    memory_per_slot_gb: 8
    runtime_hours: 48

gex_decontamination:
  enabled: true
  integration: "splitpipe_pre_router"
  tso_matching:
    backend: "edlib"
  step1c:
    enabled: false
  combine_with_grna_arm:
    enabled: false
  existing_alignment:
    enabled: true
    source_output_dir_template: "/medpop/rgupta/Jiahao/path/to/old_output/gex_{library}"
    sandbox_output_dir_template: "{output_dir}/gex_{library}_v1_4"
```

Copy the remaining matching, I/O, and strict-threshold fields from the matching
public v1.4 config rather than relying on implicit values. The old raw FASTQs
must still be readable: the sandbox route scans them once to reproduce the
pre-STAR decisions, but reuses the old aligned BAM and does not rerun STAR.

### Start decontamination of the existing Step 1 run

On the Broad login node, dry-run the exact Step 1 migration before submitting:

```bash
PROJECT_DIR=/medpop/rgupta/Jiahao/path/to/CSM_ccperturbseq
CONFIG=/medpop/rgupta/Jiahao/path/to/old_project.v1_4.yaml
cd "${PROJECT_DIR}"

source /broad/software/scripts/useuse
reuse -q UGER || true
command -v qsub

bash submit/submit_pipeline.sh --config "${CONFIG}" --step 1 --dry-run
bash submit/submit_pipeline.sh --config "${CONFIG}" --step 1
```

This sandbox route avoids rerunning human STAR and never modifies the source
BAM. It cannot retroactively apply v1.4’s improved full/fuzzy-TSO trimming to
surviving reads; use a fresh Step 1/STAR run when that change matters.

### Optional Step 1c and guide-arm integration

Step 1c is a conservative identity-recovery option, not part of required GEX
decontamination. It reads only the compact corrected-barcode evidence written
by Step 1a, so it does not rescan the full GEX FASTQs. Enable
`gex_decontamination.step1c.enabled: true` before the Step 1 submission if the
GEX R1 is long enough to resolve the required guide and, for a dual design,
iBAR. Keep it disabled for 54–60-nt dual-guide/iBAR R1, where complete construct
identity is not observed.

First leave `combine_with_grna_arm.enabled: false`, run Step 1, and inspect the
GEX-only products under `gex_grna_identity/`. To add accepted GEX guide evidence
to the dedicated-guide arm, set both Step 1c and
`combine_with_grna_arm.enabled` to true and keep
`combine_before_umi_collapse: true`.

Because the migration deliberately uses a new `paths.output_dir`, that new
namespace must contain complete dedicated Step 3 products before merging. The
safest route is to rerun the much smaller dedicated Steps 2 and 3 there:

```bash
bash submit/submit_pipeline.sh --config "${CONFIG}" --step 2,3 --dry-run
bash submit/submit_pipeline.sh --config "${CONFIG}" --step 2,3
```

After Step 1c and dedicated Step 3 are complete, rebuild Steps 4 and 5:

```bash
bash submit/submit_pipeline.sh --config "${CONFIG}" --step 4,5 --dry-run
bash submit/submit_pipeline.sh --config "${CONFIG}" --step 4,5
```

This standalone Step 4 reads the existing dedicated Step 3 and Step 1c
products, preserves dedicated-only and GEX-only matrices, combines source reads
before read-depth filtering and directional UMI collapse, and counts an exact
cross-source UMI once. A verified copy of old dedicated Step 3 products may be
used instead, but do not make the new output root writable over the historical
production run. For a wholly new run, enabling both options and submitting all
steps creates the required dependency chain automatically.

## GEX carry-over QC quick start

The QC is independent of production filtering. Eight bound workers are the
recommended full-scan setting:

```bash
bash submit/submit_qc_gex_grna_contamination_pipeline.sh \
  --config "$(realpath qc_config_gex_grna_contamination.yaml)" \
  --mode all --workers 8 --dry-run
```

The raw scan, controlled full/short54 STAR arms, historical BAM audit, and
optional guide-assignment arm are separately completeness-tracked. See the
[QC guide](docs/GEX_GRNA_CONTAMINATION_QC.md) before enabling decontamination
for a new design.

## More documentation

- [Architecture and data flow](pipeline_overview.md)
- [Detailed workflow and resources](docs/DETAILED_WORKFLOW.md)
- [v1.4 operating guide](docs/GEX_DECONTAMINATION_V1_4.md)
- [GEX carry-over QC](docs/GEX_GRNA_CONTAMINATION_QC.md)
- [100-cycle recommendations](docs/ILLUMINA_100_CYCLE_SEQUENCING.md)
- [Release notes](RELEASE_NOTES.md)

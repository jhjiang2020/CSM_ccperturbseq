# CSM ccPerturbSeq Pipeline

This pipeline processes Parse Biosciences GEX FASTQs and dedicated
CROPseq-multi guide libraries into analysis-ready AnnData. It supports
single- and dual-guide/iBAR designs, any number of sublibraries, and optional
v1.4 removal of guide products that escaped GEX/gRNA size selection.

For architecture, see [pipeline_overview.md](pipeline_overview.md). Detailed
logic is in [docs/DETAILED_WORKFLOW.md](docs/DETAILED_WORKFLOW.md); v1.4 is
covered by the [operating guide](docs/GEX_DECONTAMINATION_V1_4.md) and
[QC guide](docs/GEX_GRNA_CONTAMINATION_QC.md).

## Quick start

### 1. Install the production environment

Parse split-pipe 1.6.4 requires Python 3.12.8 or newer. From the repository
root:

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
slots is not normally useful: Parse preprocessing uses at most eight workers
in this setup. Every multi-slot wrapper requests matching SMP and linear core
binding, and workers derive their count from `NSLOTS`.

### 3. Initialize UGER and dry-run

On a Broad login node:

```bash
source /broad/software/scripts/useuse
reuse -q UGER || true
command -v qsub

bash submit/submit_pipeline.sh \
  --config "$(realpath my_project.yaml)" --dry-run
```

Heavy commands, including split-pipe, must run through `qsub`. The dry run
prints the exact resource, binding, dependency, config, and log arguments.

### 4. Submit and monitor

```bash
# All production steps
bash submit/submit_pipeline.sh --config "$(realpath my_project.yaml)"

# Selected steps
bash submit/submit_pipeline.sh \
  --config "$(realpath my_project.yaml)" --step 1

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

## v1.4 quick start for a fresh run

Enable the pre-STAR filter while leaving identity recovery opt-in:

```yaml
gex_decontamination:
  enabled: true
  integration: "splitpipe_pre_router"
  tso_matching:
    backend: "edlib"
  io:
    write_rejected_raw_pairs: true
  step1c:
    enabled: false
    strict_matching:
      enabled: true
      max_guide_edit_distance: 1
      max_ibar_edit_distance: 1
      require_resolved_ibar_in_dual_mode: true
      reject_ibar_decoupled: true
  combine_with_grna_arm:
    enabled: false
```

Step 1a removes only reference-supported TSO–guide reads. It stores full
paired positives by default because the observed C2 burden was only about 30
per million reads. Filter-negative full/fuzzy-TSO reads are trimmed at the
detected full TSO endpoint; reads without one retain Parse’s exact-17-mer
fallback.

Optional Step 1c does **not** rescan the GEX FASTQs or recalibrate barcodes. It
runs ordinary Step 3 only on Step 1a’s small valid-positive evidence table.
The default job therefore requests one core; if a caller allocates more slots,
Step 1c starts at most one worker per 100,000-row chunk and never starts idle
workers.
The unfiltered Step 3 result is retained, while conservative configurable
guide/iBAR thresholds determine which rows may augment the guide matrix. The
release defaults accept edit distance at most one and, in dual-guide mode,
require a resolved, guide-compatible corrected iBAR; decoupled guide/iBAR calls
are rejected. This intentionally favors specificity over recovery because the
optional arm can only add guide UMIs.

## v1.4 migration after a successful v1.3-or-earlier Step 1

Do not edit or overwrite the old run. Copy its configuration and select a new
output root:

```bash
cp old_project.yaml old_project.v1_4.yaml
conda activate splitpipe
python -m pip install -r requirements_extra.txt
```

Add `uger.step1a` and `gex_decontamination`, then configure the historical BAM
route:

```yaml
paths:
  output_dir: "/new/project/output_v1_4"

uger:
  step1a:
    slots: 16
    memory_per_slot_gb: 8
    runtime_hours: 48

gex_decontamination:
  enabled: true
  existing_alignment:
    enabled: true
    source_output_dir_template: "/old/read_only/output/gex_{library}"
    sandbox_output_dir_template: "{output_dir}/gex_{library}_v1_4"
```

Verify and submit only Step 1:

```bash
bash submit/submit_pipeline.sh \
  --config "$(realpath old_project.v1_4.yaml)" --step 1 --dry-run
bash submit/submit_pipeline.sh \
  --config "$(realpath old_project.v1_4.yaml)" --step 1
```

This sandbox route avoids rerunning human STAR and never modifies the source
BAM. It cannot retroactively apply v1.4’s improved full/fuzzy-TSO trimming to
surviving reads; use a fresh Step 1/STAR run when that change matters.

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

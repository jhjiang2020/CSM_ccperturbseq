# v2.0 quick start

This guide gets paired GEX and dedicated-guide FASTQs to
`combined_adata.h5ad`. The defaults are **STAR, Parse GEX UMI correction and
 gene-level output**. Kallisto is an optional Step 1a backend. Both retain
Parse split-pipe 1.6.4.

## 1. Choose your installation route

| Your access | Route |
| --- | --- |
| Broad with access to the shared rgupta software/reference folders | Clone your own project checkout; optionally reuse the read-only production environment and compatible references. See [shared setup](INSTALLATION.md#optional-shared-rgupta-setup). |
| Broad without that access | Install Parse 1.6.4 and dependencies in your own storage and prepare your own references. Use the same UGER launcher with your own absolute paths. |
| Another institution or Linux server | Install your own dependencies/references. Use your scheduler to allocate resources, then invoke the existing worker scripts with the [portable recipe](RUNNING.md#another-scheduler-or-an-existing-linux-allocation). |

No workflow input is fetched automatically from rgupta. Raw reads, guide design,
references and output storage must be available to your jobs. On a cluster,
installation, environment solving, reference building and data inspection belong
on compute nodes. Only cloning, editing, scheduler submission and lightweight
job/log checks belong on a login host.

## 2. Obtain the code and dependencies

```bash
git clone https://github.com/jhjiang2020/CSM_ccperturbseq.git
cd CSM_ccperturbseq
PROJECT_DIR="$(pwd -P)"
git rev-parse HEAD
```

Use a checkout containing `scripts/run_kallisto_gex.py` and
`scripts/kallisto_count_parallel.py` for the v2.0 backend. Do not assume a
`v2.0` Git tag has already been published. The downstream `analysis/` submodule
is optional for pipeline execution.

[Install Parse 1.6.4 and `requirements_extra.txt`](INSTALLATION.md#own-production-environment)
in the production `splitpipe` environment, or configure the shared prefix.
Kallisto users can use the [published shared companion](INSTALLATION.md#shared-kallisto-companion-published-2026-10-22),
or create the [separate tools environment and reference](INSTALLATION.md#optional-kallisto-tools-and-reference).
Keep the Parse distribution pristine.

## 3. Make your run configuration

```bash
cp config_dual.yaml my_project.v2.yaml
# For single-guide designs, start with config_single.yaml instead.
CONFIG="$(realpath my_project.v2.yaml)"
```

Edit the copy. **Every example project path must be replaced** if it is not
readable and appropriate for your experiment. Do not edit the shared lab config.

| Required edit | Meaning |
| --- | --- |
| `environment.conda_prefix` | Your Conda base containing `condabin/conda`, or the documented shared standalone prefix. |
| `environment.conda_env` | Production environment name or absolute prefix. |
| `paths.project_dir` | Absolute path to this checkout. |
| `paths.fastq_dir` | Absolute directory containing your paired, demultiplexed FASTQs. |
| `paths.splitpipe_dir` | Pristine Parse 1.6.4 distribution, including its launcher and source modules. |
| `paths.genome_dir` | Compatible Parse reference; required even when kallisto is selected. |
| `paths.grna_reference` | Your formatted guide/iBAR reference CSV. |
| `paths.output_dir`, `paths.log_dir` | New writable locations for this run. |
| `libraries` | One entry for each paired GEX/guide sublibrary; replace the pilot examples. |
| `splitpipe.kit`, `chemistry`, `sample_wells` | Match the actual experiment and sample wells. |
| `grna_structure`, `guide_assignment` | Match single/dual-guide mode, anchors, guide lengths and iBAR presence. |

Use absolute values for production paths. YAML does not expand shell variables
such as `${PROJECT_DIR}`. See [configuration details](CONFIGURATION.md) for
matching options and naming constraints.

For two sublibraries from one biological sample, the library section might be:

```yaml
libraries:
  - sublibrary_name: S1
    gex_prefix: MY_GEX_S1
    grna_prefix: MY_GUIDE_S1
    sample_name: SampleA
  - sublibrary_name: S2
    gex_prefix: MY_GEX_S2
    grna_prefix: MY_GUIDE_S2
    sample_name: SampleA
```

The GEX worker finds `*<gex_prefix>*_R1_*fastq.gz` and corresponding R2 files.
Use an unambiguous prefix and matching lane names; inspect the resolved list in
the Step 1a log. Do not assume explicit-pair settings supported by QC utilities
also override this production shell discovery. All library ordering must stay
fixed from Step 1b through Step 5.

## 4. Prepare the guide reference and check its orientation

If the reference CSV already exists, use it directly. Otherwise configure
`design.design_file` and run Step 0 **in a compute allocation**:

```bash
conda activate splitpipe
python scripts/00_format_reference.py --config "${CONFIG}"
```

Steps 1–5 read `paths.grna_reference`, not the original design CSV. Step 0
normally reverse-complements design-strand iBARs into the sequencing-read
orientation. Do not reverse-complement an already formatted reference again.
A single-guide/no-iBAR run still needs an `ibar_seq` column, with empty values.
See [reference formats](CONFIGURATION.md#step-0-guide-reference-and-identity-model).

Run the bounded guide preflight in that allocation:

```bash
python scripts/03_match_grna_ibar.py --config "${CONFIG}" \
  --preflight-only --preflight-fastq
```

Review orientation/length warnings before a large run. The check diagnoses the
reference; it does not repair it. [Details](SCALABILITY_AND_PREFLIGHT.md).

## 5. Select the GEX backend and UMI policy

For the default workflow:

```yaml
gex_alignment:
  backend: star
  umi_correction: parse
```

For native kallisto, replace that section with:

```yaml
gex_alignment:
  backend: kallisto
  umi_correction: parse
  kallisto:
    executable: kallisto
    reference_dir: references/kallisto_GRCh38
    gene_assignment: unique
    output_level: gene
    count_threads: 0
```

Set `executable` to the absolute kallisto binary in your tools environment when
it is not on the production environment's PATH. `reference_dir` here is relative
to the YAML file; it must contain the manifest produced by the repository's
reference builder. A plain kallisto index alone is insufficient.

`parse`, `directional` and `none` are available on **both** backends. Fractional
gene assignment and transcript output are kallisto options. Use
`output_level: both` to keep genes primary while also exporting transcript/target
counts. `output_level: transcript` changes the primary expression feature level.
[Counting semantics and output limitations](KALLISTO_BACKEND.md).

## 6. Decide whether to remove gRNA carry-over from GEX

The public configs enable the production filter for the standard CSM construct.
Confirm the guide reference and TSO sequence before retaining that setting.
Filtering takes place during Parse preprocessing and works with either backend.

```yaml
gex_decontamination:
  enabled: true
  integration: splitpipe_pre_router
  # Retain the rest of this section from the public config.
```

Keep `step1c.enabled` and `combine_with_grna_arm.enabled` false for a first run
unless the GEX R1 reaches the guide/iBAR sequence needed for your identity model.
Removal is useful even when reads are too short for that optional recovery.
For short-GEX experiments, choose the [sequencing scenario](ILLUMINA_100_CYCLE_SEQUENCING.md#choose-the-scenario)
before enabling guide-source merging. The [decontamination guide](GEX_DECONTAMINATION_V1_4.md)
explains the full configuration, evidence and historical-STAR BAM route.

## 7. Dry-run and submit

On a Broad login node, after finishing setup:

```bash
cd "${PROJECT_DIR}"
source /broad/software/scripts/useuse
reuse -q UGER
CONFIG="$(realpath my_project.v2.yaml)"
bash submit/submit_pipeline.sh --config "${CONFIG}" --dry-run
bash submit/submit_pipeline.sh --config "${CONFIG}"
qstat -u "$USER"
```

The launcher submits Steps 1–5 with dependencies; Step 0 is separate. Arrays
follow the number of configured libraries. Default Step 1a allocation is 16
slots, 8 GB **per slot** and 48 hours. Kallisto counting uses those allocated
cores unless capped by `count_threads`; `splitpipe.threads` is a legacy field,
not the current Step 1a scheduler control.

For another scheduler, use the [portable worker recipe](RUNNING.md#another-scheduler-or-an-existing-linux-allocation).
Do not use the deprecated `scripts/01_run_gex_splitpipe.sh` as a local runner:
it delegates to the UGER launcher.

## 8. Check the result

The final file is `paths.output_dir/combined_adata.h5ad`. Check scheduler exit
status and logs as well as file presence. For kallisto, Step 1a/1b completion
manifests validate reference, settings and output hashes; Step 5 validates the
combined product before using it. For STAR, UMI policy markers protect new runs.

Inspect called-cell yield, UMIs/genes per cell, guide assignment and the chosen
expression feature level. Read H5AD or matrices on a compute node. Use
[the runbook](RUNNING.md#outputs-and-success-checks) for the output map,
validation commands and restart guidance. Keep gene counts primary for QC or
downstream methods that expect gene identifiers; do not silently round
fractional counts to satisfy an integer-only method.

# v2.0 quick start

This guide gets paired GEX and dedicated-guide FASTQs to
`combined_adata.h5ad`. The defaults are **STAR, Parse GEX UMI correction and
 gene-level output**. Kallisto is an optional Step 1a backend. Both retain
Parse split-pipe 1.6.4.

## 1. Choose your installation route

| Your access | Route |
| --- | --- |
| Broad with access to the shared rgupta software/reference folders | Clone your own checkout and copy [`config_rgupta.yaml`](../config_rgupta.yaml); the shared environment, kallisto executable and both GRCh38 references are preconfigured. |
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

**rgupta lab users on UGER can reuse this complete shared setup.** The
[`config_rgupta.yaml`](../config_rgupta.yaml) preset contains all these settings:

```yaml
environment:
  conda_prefix: "/medpop/rgupta/Jiahao/conda/envs/splitpipe-share"
  conda_env: "/medpop/rgupta/Jiahao/conda/envs/splitpipe-share"
paths:
  splitpipe_dir: "/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq/ParseBiosciences-Pipeline.1.6.4"
  # split-pipe-generated metadata for both backends; STAR also uses its index.
  genome_dir: "/medpop/rgupta/Jiahao/reference_genome/parse_refdata-gex-GRCh38-2024-A"
gex_alignment:
  backend: star                 # change to kallisto to use the shared stack below
  gene_assignment: unique        # or fractional
  umi_correction: parse          # or directional / none
  output_level: gene
  count_threads: 0
  kallisto:
    executable: "/medpop/rgupta/Jiahao/conda/envs/csm_kallisto-share/bin/kallisto"
    reference_dir: "/medpop/rgupta/Jiahao/reference_genome/kallisto_refdata-gex-GRCh38-2024-A_nac_20261023"
    strand: forward
```

No personal Conda/kallisto installation, kallisto module or reference build is
needed. Both references use GRCh38 2024-A / GENCODE v44; the kallisto reference
was validated with the shared 0.52.0 executable. Keep these dependency paths
and put your run's inputs, output and logs in your own storage. See
[shared setup and access requirements](INSTALLATION.md#optional-shared-rgupta-setup)
and [reference validation and permissions](INSTALLATION.md#shared-grch38-kallisto-reference).

`paths.genome_dir` provides Parse's canonical gene IDs, names, feature order and
genome/species labels for the integrated GEX outputs and QC. Kallisto's
target-to-gene assignments are checked against that metadata; pseudoalignment
uses the separate kallisto index above. The Parse directory is generated
directly with split-pipe outside this pipeline and can be shared as shown.
See [reference preparation and future metadata/index separation](INSTALLATION.md#parse-reference-for-both-backends).

Without shared access, [install Parse 1.6.4 and `requirements_extra.txt`](INSTALLATION.md#own-production-environment)
in your production `splitpipe` environment and prepare compatible references.
Kallisto users can create the [separate tools environment and reference](INSTALLATION.md#optional-kallisto-tools-and-reference).
Tools and references may be mixed: a compatible personal kallisto 0.52.0 build
can use the shared reference, or the shared executable can build/use your own
matching reference. Check index compatibility as described in the installation guide.
Keep the Parse distribution pristine.

## 3. Make your run configuration

Choose a complete rgupta preset for your single-target/one-guide design. Both
use the same shared software and references, and support either GEX backend:

| Guide design | iBAR | Preset | Step 4 output |
| --- | --- | --- | --- |
| Single-target / one guide | Present | [`config_rgupta_single.yaml`](../config_rgupta_single.yaml) | Individual guide calls, with iBAR matching. |
| Single-target / one guide | Absent | [`config_rgupta_single_no_ibar.yaml`](../config_rgupta_single_no_ibar.yaml) | Individual guide calls, with the no-iBAR parser. |

`config_rgupta.yaml` has the same settings as the single-target/iBAR preset.
The shared presets cover single-target inputs. For short-read experiments that
do not reach iBAR, also review the [guide-only recipe](ILLUMINA_100_CYCLE_SEQUENCING.md#scenario-2-gex-and-grna-100-cycle-no-ibar).

```bash
# rgupta lab / shared UGER setup:
cp config_rgupta_single.yaml my_project.v2.yaml  # or config_rgupta_single_no_ibar.yaml
# Own installation: copy config_dual.yaml or config_single.yaml instead.
CONFIG="$(realpath my_project.v2.yaml)"
```

Edit the copy. In the lab preset, replace **every `EDIT_ME` value** with your own
absolute cluster path or FASTQ prefix. Set `analysis_name`, libraries and the
experiment settings below; keep the prefilled shared software/reference paths.
The chosen preset supplies the guide parsing, iBAR and assignment mode; confirm
these against the actual construct. For the other templates, replace every
example path that is not appropriate for your experiment. Keep edits in your
own run config.

| Required edit | Meaning |
| --- | --- |
| `analysis_name` | Name for your analysis run. |
| `environment.conda_prefix` | Keep the lab preset; otherwise set your Conda base containing `condabin/conda`. |
| `environment.conda_env` | Keep the lab preset; otherwise set your production environment name or absolute prefix. |
| `paths.project_dir` | Absolute path to this checkout. |
| `paths.fastq_dir` | Absolute directory containing your paired, demultiplexed FASTQs. |
| `paths.splitpipe_dir` | Keep the lab preset; otherwise use your pristine Parse 1.6.4 distribution. |
| `paths.genome_dir` | Keep the lab preset; otherwise prepare a reference directly with split-pipe. Supplies canonical gene metadata for both backends and the index used by STAR. |
| `paths.grna_reference` | Your formatted guide/iBAR reference CSV. |
| `paths.output_dir`, `paths.log_dir` | New writable locations for this run. |
| `design.design_file` | Your design CSV, only if you run Step 0; otherwise unused. |
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

In a **compute allocation**, activate the production environment. For the
shared lab preset:

```bash
source /medpop/rgupta/Jiahao/conda/envs/splitpipe-share/etc/profile.d/conda.sh
conda activate /medpop/rgupta/Jiahao/conda/envs/splitpipe-share
```

For a personal installation, activate the environment in your config instead.
If the reference CSV already exists, use it directly and skip Step 0 (the
`design.design_file` placeholder is unused). Otherwise set `design.design_file`
and run Step 0 in that allocation:

```bash
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

## 5. Select the GEX backend, assignment and UMI policy

For the default workflow:

```yaml
gex_alignment:
  backend: star
  gene_assignment: unique
  umi_correction: parse
```

**With `config_rgupta.yaml`, change only `gex_alignment.backend` to `kallisto`**;
its executable and reference are already set. Use a new output/log root when
switching backend or counting policies. For a personal kallisto setup, use:

```yaml
gex_alignment:
  backend: kallisto
  gene_assignment: unique
  umi_correction: parse
  output_level: gene
  count_threads: 0
  kallisto:
    executable: kallisto
    reference_dir: references/kallisto_GRCh38
```

Set `executable` to the absolute kallisto binary in your tools environment when
it is not on the production environment's PATH. `reference_dir` here is relative
to the YAML file; it must contain the manifest produced by the repository's
reference builder. A plain kallisto index alone is insufficient.

`parse`, `directional` and `none` are available on **both** backends. Fractional
gene assignment is available on both aligners; transcript output requires kallisto. Use
`output_level: both` to keep genes primary while also exporting transcript/target
counts. `output_level: transcript` changes the primary expression feature level.
[Counting choices, Parse defaults and Tony’s workflow](GEX_COUNTING.md).

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

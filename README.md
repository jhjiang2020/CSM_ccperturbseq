# CSM ccPerturbSeq Pipeline — v2.0

Process Parse Biosciences GEX and CROPseq-multi guide FASTQs into an
analysis-ready `combined_adata.h5ad`. The workflow supports single-guide and
dual-guide/iBAR designs, multiple lanes and sublibraries, optional removal of
gRNA carry-over from GEX, and optional recovery of those guide identities.

**Start with the [quick start](docs/QUICKSTART.md).** It covers a new installation,
the optional shared rgupta setup, and running on Broad UGER or another Linux
compute environment. Access to the rgupta folder is not required.

## What changes in v2.0

- **STAR stays the default.** The optional kallisto backend changes transcriptomic
  alignment and counting in Step 1a. Parse split-pipe **1.6.4 remains required**
  for preprocessing, barcode/chemistry definitions, sample handling and cell
  calling; its installation stays unchanged.
- **GEX UMI correction is configurable for both backends:**
  `gex_alignment.umi_correction: parse` (default), `directional`, or `none`.
  The directional condition is `n_high >= 2*n_low - 1`. This setting does not
  change dedicated or GEX-derived **gRNA** UMI correction in Step 4.
- **Both aligners support unique/fractional gene assignment** through
  `gex_alignment.gene_assignment`. STAR fractional counting retains secondary
  alignments and uses the same molecule/output core as kallisto.
- **Native kallisto counting writes compatible sparse products** for Step 1b
  and downstream integration. Gene counts are the default; equal-fractional
  gene assignment and transcript/target output are optional. No BAM projection
  is used to feed pseudoalignments into Parse gene counting.
- **Parallel molecule counting** keeps each corrected barcode intact and merges
  results in barcode order. `gex_alignment.count_threads: 0` uses the allocated cores;
  `1` selects serial counting. Fractional values are retained.

The [parallel counting benchmark](docs/KALLISTO_PARALLEL_BENCHMARK_REPORT.md)
compares the optimized counter with saved baseline counts and runtime. Counting
uses allocated cores by default; `gex_alignment.count_threads` can cap it.

See [counting choices and Parse/Tony presets](docs/GEX_COUNTING.md), the [release notes](RELEASE_NOTES.md), [backend guide](docs/KALLISTO_BACKEND.md)
and [migration guide](docs/MIGRATING_TO_V2.md). Use a separate output directory
when changing backend, UMI policy, reference or feature level.

## Quick start

```bash
git clone https://github.com/jhjiang2020/CSM_ccperturbseq.git
cd CSM_ccperturbseq
cp config_dual.yaml my_project.v2.yaml    # or config_single.yaml
```

Use the checkout containing the v2.0 implementation; record its commit. A Git
release tag is not assumed by these instructions. For your own installation,
obtain Parse 1.6.4 separately and follow [installation and references](docs/INSTALLATION.md).
Edit the copied YAML's environment, paths, libraries, kit/chemistry and guide
reference. The rgupta shared preset below already fills in the dependency paths.

**rgupta lab users on UGER:** start with
[`config_rgupta.yaml`](config_rgupta.yaml):
`cp config_rgupta.yaml my_project.v2.yaml`. It already configures the shared
production environment, kallisto executable and **both precomputed GRCh38
references** together:

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

No personal Conda/kallisto installation, module loading or reference build is
needed with this preset. Replace every `EDIT_ME` value with your own absolute
cluster path or FASTQ prefix, and confirm the libraries, kit/chemistry, wells
and guide construct. Keep the shared dependency paths; use your own writable
checkout, output and logs. Change only `gex_alignment.backend` to select
kallisto; keep `paths.genome_dir` with either backend. Complete shared examples cover
[single-target with iBAR](config_rgupta_single.yaml) and
[single-target without iBAR](config_rgupta_single_no_ibar.yaml).
`config_rgupta.yaml` uses the single-target/iBAR settings.
See the [preset table](docs/QUICKSTART.md#3-make-your-run-configuration)
for read-layout requirements and the [shared setup](docs/INSTALLATION.md#optional-shared-rgupta-setup)
for validation, permissions and [personal installation options](docs/INSTALLATION.md#optional-kallisto-tools-and-reference).

`paths.genome_dir` supplies Parse's canonical gene table—IDs, names, feature
order and genome/species labels—for the integrated GEX outputs and QC.
Kallisto's target-to-gene assignments are checked against that metadata;
pseudoalignment uses its own kallisto index. Generate the Parse directory
directly with split-pipe outside this pipeline, or reuse the shared directory.
See [reference preparation and future metadata/index separation](docs/INSTALLATION.md#parse-reference-for-both-backends).

The default GEX choice is:

```yaml
gex_alignment:
  backend: star
  gene_assignment: unique
  umi_correction: parse
```

For a personal kallisto setup, configure its executable and a validated reference
built from the same FASTA/GTF release as the Parse reference. The lab preset
already supplies these fields; retain its absolute paths:

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

`executable: kallisto` works only if the executable is on the **production job's**
PATH. An absolute executable path is preferable for a separate tools environment.
Other production path fields should be absolute; see the
[configuration reference](docs/CONFIGURATION.md).

On Broad UGER, initialize the scheduler and submit from the project checkout:

```bash
source /broad/software/scripts/useuse
reuse -q UGER
CONFIG="$(realpath my_project.v2.yaml)"
bash submit/submit_pipeline.sh --config "${CONFIG}" --dry-run
bash submit/submit_pipeline.sh --config "${CONFIG}"
qstat -u "$USER"
```

Run installation, reference construction and data processing in compute
allocations. The master launcher's YAML parsing/dry run is scheduler bookkeeping.
For another scheduler, follow the [portable worker recipe](docs/RUNNING.md#another-scheduler-or-an-existing-linux-allocation).

## Documentation

| Need | Guide |
| --- | --- |
| First successful run | [Quick start](docs/QUICKSTART.md) |
| Own installation or optional shared Broad setup | [Installation and reference preparation](docs/INSTALLATION.md) |
| YAML options, defaults and which stage consumes them | [Configuration reference](docs/CONFIGURATION.md) |
| Submit, monitor, resume and check outputs | [Execution runbook](docs/RUNNING.md) |
| CPU, memory, time, queues and array limits | [UGER resources](docs/UGER_RESOURCES.md) |
| Overall data flow | [Architecture](pipeline_overview.md) |
| Stage inputs, outputs and counting logic | [Detailed workflow](docs/DETAILED_WORKFLOW.md) |
| Kallisto, fractional genes, transcript output, GEX UMI policies | [Backend guide](docs/KALLISTO_BACKEND.md) |
| GEX gRNA removal, Step 1c and source merging | [Decontamination guide](docs/GEX_DECONTAMINATION_V1_4.md) |
| Read length, iBAR reach and short-read configurations | [Sequencing scenarios](docs/ILLUMINA_100_CYCLE_SEQUENCING.md) |
| Reference preflight and resource scaling | [Preflight/scaling guide](docs/SCALABILITY_AND_PREFLIGHT.md) |
| Existing run → v2.0 | [Migration guide](docs/MIGRATING_TO_V2.md) |
| QC and historical validation reports | [Documentation index](docs/README.md) |

The final product is `paths.output_dir/combined_adata.h5ad`. Per-library GEX,
combined GEX, guide counts and optional contaminant evidence remain available
under that output root; [the runbook](docs/RUNNING.md#outputs-and-success-checks)
explains their contracts.

This repository owns FASTQ-to-AnnData processing and count/QC validation.
`analysis/` is the separate `CSM_analysis` Git submodule for downstream analyses;
it is not required to run Steps 0–5. Historical reports record their original
methods and datasets and are not current installation recipes.

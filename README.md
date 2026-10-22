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
- **Native kallisto counting writes compatible sparse products** for Step 1b
  and downstream integration. Gene counts are the default; equal-fractional
  gene assignment and transcript/target output are optional. No BAM projection
  is used to feed pseudoalignments into Parse gene counting.
- **Parallel molecule counting** keeps each corrected barcode intact and merges
  results in barcode order. `kallisto.count_threads: 0` uses the allocated cores;
  `1` selects serial counting. Fractional values are retained.

The [parallel counting benchmark](docs/KALLISTO_PARALLEL_BENCHMARK_REPORT.md)
compares the optimized counter with saved baseline counts and runtime. Counting
uses allocated cores by default; `gex_alignment.kallisto.count_threads` can cap it.

See the [release notes](RELEASE_NOTES.md), [backend guide](docs/KALLISTO_BACKEND.md)
and [migration guide](docs/MIGRATING_TO_V2.md). Use a separate output directory
when changing backend, UMI policy, reference or feature level.

## Quick start

```bash
git clone https://github.com/jhjiang2020/CSM_ccperturbseq.git
cd CSM_ccperturbseq
cp config_dual.yaml my_project.v2.yaml    # or config_single.yaml
```

Use the checkout containing the v2.0 implementation; record its commit. A Git
release tag is not assumed by these instructions. Obtain Parse 1.6.4 separately,
then follow [installation and references](docs/INSTALLATION.md). Edit the copied
YAML's environment, all paths, libraries, kit/chemistry and guide reference.
The public configs contain example project paths and must be customized.

For authorized lab members, an **installation-free shared environment** is
available. In the copied YAML, use:

```yaml
environment:
  conda_prefix: "/medpop/rgupta/Jiahao/conda/envs/splitpipe-share"
  conda_env: "/medpop/rgupta/Jiahao/conda/envs/splitpipe-share"
```

The [shared setup](docs/INSTALLATION.md#optional-shared-rgupta-setup) also lists
the validated dated kallisto companion and optional lab references. Configure
your own writable project, output and logs; the shared environment is a
dependency. Without shared-folder access, follow the own-installation route.

The default GEX choice is:

```yaml
gex_alignment:
  backend: star
  umi_correction: parse
```

To use kallisto, add its executable and the validated reference built from the
same FASTA/GTF release as the Parse reference:

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

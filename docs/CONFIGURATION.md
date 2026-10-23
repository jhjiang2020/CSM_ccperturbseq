# v2.0 configuration reference

Start with `config_dual.yaml` or `config_single.yaml`, copy it under a run-specific
name, and edit the copy. The public files are complete examples with lab-specific
paths; they are not ready-to-run configs for another project. Use
[the quick start](QUICKSTART.md) for the shortest setup route.

## Paths, environment and libraries

| Section/key | Used by | What to set |
| --- | --- | --- |
| `environment.conda_prefix` | Workers | Conda base containing `condabin/conda`; shared standalone prefix is a documented exception. |
| `environment.conda_env` | Workers | Production `splitpipe` environment name or absolute prefix. |
| `paths.project_dir` | Workers | Absolute v2.0 checkout path; determines the scripts actually executed. |
| `paths.fastq_dir` | Steps 1a/2 and raw-read QC | Your paired FASTQ directory. |
| `paths.splitpipe_dir` | Steps 1a/1b, barcode decoding and QC | Pristine Parse 1.6.4 distribution. |
| `paths.genome_dir` | Step 1 and controlled QC | Complete Parse reference, including when kallisto is selected. |
| `paths.grna_reference` | Guide matching, GEX filtering, assignment | Formatted guide/iBAR CSV. |
| `paths.output_dir` | All processing stages | New writable output root for this parameter set. |
| `paths.log_dir` | Launcher/workers | Writable scheduler log directory. |
| `libraries[].sublibrary_name` | All stages | Unique stable identifier for a paired GEX/dedicated-guide sublibrary. |
| `libraries[].gex_prefix` | Step 1a | Unambiguous substring/glob in conventional paired GEX filenames. |
| `libraries[].grna_prefix` | Step 2 | Dedicated-guide FASTQ identifier. |
| `libraries[].sample_name` | Steps 1/5 | Biological sample label; equal labels pool sublibraries in Step 1b. |

Use absolute production paths, particularly in copied configs passed to compute
nodes. Shell variables in YAML are not expanded. Kallisto's `reference_dir` and
slash-containing `executable` are resolved relative to the YAML when relative;
this behavior is not a promise that every legacy path key is config-relative.

The Step 1a shell worker searches
`*<gex_prefix>*_R1_*fastq.gz` and matching R2 names. It accepts multiple lanes,
but the broad glob can overmatch: choose specific prefixes, keep paired lane
names/counts consistent, and inspect the printed file list. Explicit FASTQ-pair
options in QC resolvers do not override that shell worker. Keep `libraries`
order unchanged after counting: Step 1b's sample-local suffixes (`__s1`, etc.)
and Step 5's sublibrary mapping depend on that order.

## Step 0: guide reference and identity model

Steps 1–5 read the formatted CSV at `paths.grna_reference`. Only Step 0 reads
`design.design_file`. `design.rev_comp_ibar: true` converts CSM design-strand
iBARs to read orientation during formatting; changing it does not modify an
already generated reference.

The minimum matching columns are `grna_name`, `grna_seq`, `ibar_seq` and
`target_gene`. A no-iBAR reference retains an empty `ibar_seq` column. Dual-guide
inference also needs correct position/pair annotations such as `position`,
`stem_variant`, `pair_id`, `target_combo` and `target_version`; generate these
with Step 0 from the actual CSM design rather than inventing pair IDs from
individual guide names. The formatter's required design fields are documented
in [the detailed workflow](DETAILED_WORKFLOW.md#step-0-format-reference) and the
public YAML comments.

| Identity model | `grna_structure.single_guide_mode` | `grna_structure.ibar_length` | `guide_assignment.dual_guide_counting_mode` |
| --- | --- | --- | --- |
| Standard dual CSM with iBAR | `false` | Actual iBAR length (example: `7`) | `true` |
| Single guide, no iBAR | `true` | `0` | `false` |
| Guide-only analysis of reads too short for iBAR | Choose the guide-only parsing configuration in the sequencing guide | `0` | `false` |

The standard examples use 20/21-nt guides. Confirm `grna_lengths`, TSO,
stem/tracr anchors and iBAR length against the construct. Do not turn off iBARs
merely to suppress orientation warnings in a dual-guide analysis.
[Read-length scenarios](ILLUMINA_100_CYCLE_SEQUENCING.md) and
[bounded reference preflight](SCALABILITY_AND_PREFLIGHT.md) explain those choices.

## Step 1a: GEX backend, assignment and UMI correction

```yaml
gex_alignment:
  backend: star
  gene_assignment: unique
  umi_correction: parse
```

Omitting the entire section preserves STAR and Parse-default correction.

| Key | Default | Effect |
| --- | --- | --- |
| `gex_alignment.backend` | `star` | `star` or `kallisto`; only the transcriptomic alignment/counting path changes. |
| `gex_alignment.gene_assignment` | `unique` | `unique` or `fractional`, on either aligner. |
| `gex_alignment.output_level` | `gene` | STAR: gene only; kallisto: gene, both, transcript. |
| `gex_alignment.transcript_assignment` | `fractional` | Kallisto target assignment: unique or fractional. |
| `gex_alignment.count_threads` | `0` | Shared counter uses allocated cores; positive integer caps workers. |
| `gex_alignment.star.max_multimaps` | `3` | Maximum STAR loci per read; all reported hits are retained for fractional counting. |
| `gex_alignment.umi_correction` | `parse` | Applies to **both** GEX backends; independent of Step 4 gRNA correction. |
| `splitpipe.kit`, `chemistry` | Values in copied template | Must match the Parse experiment, not the aligner choice. |
| `splitpipe.sample_wells` | Example: `A1-B12` | Sample well range used for each configured sublibrary/sample. |
| `splitpipe.bypass_kit_score` | Template: `true` | Adds `--kit_score_skip` when enabled. |
| `splitpipe.threads` | Legacy example: `16` | Compatibility field; current workers use allocated `NSLOTS`. |
| `splitpipe.memory_gb` | Legacy example: `64` | Does not set the launcher's UGER reservation; use the resource controls below. |

| GEX UMI flavor | Molecule semantics |
| --- | --- |
| `parse` | Original Parse 1.6.4 strict-abundance Hamming-one, non-transitive correction. Equal-abundance neighbors do not merge. |
| `directional` | Hamming-one graph edges satisfy `n_high >= 2*n_low - 1` using original read support; transitive traversal follows bustools-compatible root ownership. Singleton ties can merge. |
| `none` | No correction between distinct UMI strings; exact duplicates within a counting group still form one molecule. |

Corrected barcode-index/RT identities remain separate while counting and are
aggregated to physical wells afterward. Parse's cell-calling algorithm is
retained; its selected cells can differ when a backend or assignment policy
changes the counts. A new output root is required when comparing policies.

### Kallisto-specific keys

```yaml
gex_alignment:
  backend: kallisto
  gene_assignment: unique
  umi_correction: parse
  output_level: gene
  count_threads: 0
  transcript_assignment: fractional
  kallisto:
    executable: kallisto
    reference_dir: references/kallisto_GRCh38
    strand: forward
    cell_umi_threshold: null
```

| Key under `gex_alignment.kallisto` | Default | Meaning |
| --- | --- | --- |
| `executable` | `kallisto` | Executable on production PATH or an explicit path into a separate tools environment. |
| `reference_dir` | Required | Reference directory with the repository builder's validated manifest and index. |
| `strand` | `forward` | `forward`, `reverse` or `unstranded`, according to the assay. |
| `cell_umi_threshold` | `null` | Initialized Parse automatic cell calling; a positive number explicitly fixes the cutoff. |

For rgupta users, the [shared GRCh38 setup](INSTALLATION.md#shared-grch38-kallisto-reference)
provides exact production-environment, executable and reference paths. The
[personal installation options](INSTALLATION.md#optional-kallisto-tools-and-reference)
also cover mixing a shared reference with your own tools, or shared tools with
your own reference.

Legacy nested policy keys remain accepted; conflicting common/nested values
are rejected. See [counting choices](GEX_COUNTING.md) for STAR fractional BAM
requirements and the Parse/Tony settings comparison.

Fractional gene assignment is optional and is **not** an EM or TPM estimate.
It changes the definition of the gene count, so compare it as a separate policy.
Transcript output is UMI-derived target support; a mature-plus-nascent reference
also contains nascent targets. It is not a bulk kallisto abundance table or a
separate RNA-velocity result. Some gene-resolved molecules have unresolved
transcript intersections and cannot appear in the target matrix. Cell calls
use whole gene-molecule totals, with float roundoff excluded from cutoffs.

The backend writes raw sparse counts and preserves fractional values through
Step 1b/5. `output_level: both` is the convenient choice when colleagues need
gene-based QC plus transcript exploration. Kallisto does not create the genomic
BAM required by coordinate-specific STAR QC. See [the backend guide](KALLISTO_BACKEND.md)
for collision resolution, family membership, QC and completion contracts.

## Step 1a/1c: GEX gRNA decontamination

The complete block in the public templates is the starting point. These options
apply before **either** transcriptomic backend:

| Key | Template/default behavior |
| --- | --- |
| `gex_decontamination.enabled` | Public templates: `true`; enables reference-supported guide removal during Parse preprocessing. |
| `integration` | `splitpipe_pre_router` production route. |
| `tso_max_edit_dist`, `max_tso_start` | Template: `3`, `2`; full TSO detection near the beginning of GEX R1. |
| `tso_matching.backend` | Required production matcher: `edlib`. |
| `filter_matching` | Packed full-guide Hamming matching, maximum one mismatch, no guide indels. |
| `io.write_rejected_raw_pairs` | Template: `true`; retain synchronized raw evidence. |
| `step1c.enabled` | `false`; optional guide identity recovery from compact corrected-barcode evidence. |
| `step1c.strict_matching` | Conservative acceptance after ordinary matching; actual guide/iBAR reach is required. |
| `combine_with_grna_arm.enabled` | `false`; optional dedicated/GEX guide-source union in Step 4. |
| `combine_before_umi_collapse` | `true`; union reads before depth filtering/collapse, so exact overlap is counted once. |
| `existing_alignment.enabled` | `false`; historical compatible STAR BAM reuse only, incompatible with kallisto. |

TSO alone never justifies removal. Removal and optional identity recovery are
separate decisions; leave Step 1c/source merging off for short dual-guide GEX
reads that cannot reveal the iBAR. Full policies and manifest fields are in the
[decontamination guide](GEX_DECONTAMINATION_V1_4.md).

## Steps 2–4: dedicated-guide processing

| Section/key | Scope and example default |
| --- | --- |
| `barcode_extraction.bc_max_edit_dist` | Step 2 Parse barcode correction; example: `1`. |
| `barcode_extraction.min_r1_length`, `min_r2_length` | Step 2 extraction limits, not replacements for chemistry-specific GEX/R2 requirements. |
| `barcode_extraction.max_reads` | `0` means all reads; positive values are diagnostic subsampling. |
| `barcode_extraction.n_processes` | Legacy field; current per-library Step 2 FASTQ loop remains serial. |
| `grna_matching.max_lev_dist`, `max_lev_dist_ibar` | Step 3 fuzzy guide/iBAR distances; examples: `3`, `2`. |
| `grna_matching.save_failed_reads` | `false`; optional reason-annotated failed-read TSV. |
| `grna_matching.preflight` | Enabled by default with a bounded 10,000-read check. |
| `guide_assignment.min_reads_per_umi` | `1`; gRNA read-depth filter before directional collapse. |
| `guide_assignment.max_umi_edits` | `1`; gRNA Hamming-one correction, or `0` for exact-only. |
| `guide_assignment.umi_dedup_ibar` | Include corrected/unresolved iBAR identity in appropriate groups; disabled in no-iBAR mode. |
| `guide_assignment.min_umi_threshold` | Example: `5`; pair-sum threshold for complete dual pairs, single-side threshold for half pairs. |
| `multiplet_strategy`, `min_dominance_fraction` | Call handling; examples: `all`, `0.0`. |
| `compute_saturation_curve` | `true`; guide sequencing-depth QC. |

`gex_alignment.umi_correction` does not alter any row in this table. GEX-derived
**guide** evidence from Step 1c also uses the Step 4 policy. The GEX directional
implementation and the existing gRNA traversal share an abundance condition,
but are not interchangeable algorithms. GEX read-support cutoffs come from
initialized Parse molecule parameters, not `guide_assignment.min_reads_per_umi`.

## Resources and output

All UGER resource options are optional and use one interface across Steps
1a/1b/1c/2/3/4/5, reference builds and maintained QC drivers. Old configs keep
the previous requests. A focused override is sufficient:

```yaml
uger:
  queue: broad
  defaults:
    max_concurrent: 2
  step1a: {slots: 16, memory_per_slot_gb: 8, runtime_hours: 48}
  step1b: {slots: 4, memory_per_slot_gb: 24, runtime_hours: 8}
  step5: {slots: 1, memory_per_slot_gb: 64, runtime_hours: 8}
```

Memory is **per slot**. The Step 1a example reserves 128 GB per array task;
Step 1b reserves 96 GB. See [UGER resources](UGER_RESOURCES.md) for every
profile/default, scheduler-wide fields, precedence, QC compatibility and tuning
examples. The master's dry run shows resolved resources; direct `qsub` worker
submissions use their fallback directives. `NSLOTS` controls parallel workers;
`gex_alignment.count_threads` caps counting within that allocation.

`output` contains plotting/intermediate preferences, but not every historical
key is a global switch: Step 5 currently writes `combined_adata.h5ad` directly.
Setting `output.output_format: loom` does not convert that deliverable. Preserve
intermediates and manifests needed for provenance/restarts. See
[outputs and success checks](RUNNING.md#outputs-and-success-checks).

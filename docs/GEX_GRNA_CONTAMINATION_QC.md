# GEX gRNA cross-contamination QC — v2.0 usage

Resource requests for every maintained stage are optional YAML settings; see
[UGER resources](UGER_RESOURCES.md) for CPU/RAM/time defaults, array limits,
site configuration and examples of increasing resources.

This QC module, introduced in v1.4, measures gRNA-construct carry-over in GEX libraries,
benchmarks removal rules for a 54-nt R1, and quantifies whether contaminant
reads can become false gene-expression alignments or molecules. It is an audit
workflow, not the production Step 1 decontamination path: it reads the original
FASTQs and existing split-pipe BAMs without modifying them.

## v2.0 scope and current entry points

Use the [quick start](QUICKSTART.md) and [installation guide](INSTALLATION.md)
for your own software/reference paths; a shared rgupta environment is optional.
This audit uses production `splitpipe`, not the Q1–Q4 notebook environment.
The maintained launcher is
`submit/submit_qc_gex_grna_contamination_pipeline.sh`; production removal is
configured separately with `gex_decontamination` in the Steps 1–5 config.

Raw structural scans and corrected-barcode evidence remain relevant with either
GEX backend. Controlled full/short54 alignment arms and historical annotated-BAM
audits are **STAR experiments**. A native kallisto production run does not
provide a genomic BAM; do not point those audit arms at its BUS/target products
or interpret their old STAR alignment percentages as kallisto estimates.
The measured C1/C2 values below retain their original input/policy scope.

`gex_alignment.umi_correction` controls production GEX transcriptomic counting.
It does not retroactively change a historical audit's molecule proxy/collapse
definition or the Step 4 gRNA correction policy. Report those definitions when
comparing results. [Current execution/dependency guide](RUNNING.md).

## Rationale

### The contamination risk

The assay begins with one reverse-transcription reaction containing both
poly(A)-directed probes and gRNA-directed probes. It therefore creates long
gene-transcript cDNA and much shorter gRNA-derived products in the same
reaction. The GEX and dedicated gRNA libraries are separated afterward by size
selection. A small fraction of the short gRNA product can escape that
separation and remain in the GEX fraction. The pilot 300-cycle GEX libraries
contain this low-frequency population with the expected gRNA construct:

```text
full TSO (~30 nt) - spacer (~20/21 nt) - CSM stem - iBAR - partial CSM tracr
```

In both the GEX and dedicated gRNA libraries, R1 contains the biological insert
(cDNA or gRNA) and R2 contains the Parse UMI and cell-barcoding information. For
Parse chemistry v3, the required 58-nt R2 layout is:

```text
10-nt UMI + 8-nt BC3 + 12-nt linker + 8-nt BC2 + 12-nt linker + 8-nt BC1
```

R2 is decoded and corrected but is not aligned. Cell identity is the corrected
combination of all three 8-nt barcodes, not one index sequence.

The measured raw-read frequency is low: the completed C1 and C2 pilot scans
found 36.6 and 32.1 structural-contaminant reads per million, respectively.
That frequency does not by itself bound biological impact. If STAR aligns a
guide-derived sequence to the reference—especially after a short 100-cycle
read is reduced to a guide-dominated insert—the resulting false count can be
large relative to a target gene expressed at only tens of counts per million
or less.

This carry-over is therefore especially concerning for the proposed 100-cycle
sequencing configuration. The kit supplies 138 total cycles; after reserving
10 cycles for each of i7 and i5, 118 cycles remain for R1 and R2. The preferred
configuration is `R1=60, R2=58`; the `R1=54, R2=64` validation configuration is
also accepted, although split-pipe does not decode R2 beyond base 58. See the
[100-cycle sequencing note](ILLUMINA_100_CYCLE_SEQUENCING.md).

The gRNAs in this CRISPRi screen target promoter regions. After TSO trimming, a
short contaminant R1 may be dominated by its spacer and can therefore align
uniquely within or near the genomic sequence it was designed to target. If that
alignment falls on the correct strand inside an annotated gene span, it can be
assigned to a gene and counted as expression, creating a target-concordant
artifact in the same cells used to estimate the perturbation effect. A promoter
hit outside an annotated gene span is not necessarily assigned, so this risk
must be measured rather than assumed from guide identity alone.

### Molecular assumptions used by the audit

An intact contaminant carried into a GEX library is expected to start with the
full TSO at R1 offset 0. Offsets 1-2 are accepted only as low-frequency tolerance
for 5-prime damage or sequencing error; they are not an intentional GEX prefix.
The dedicated gRNA library is different because it receives a separate
circularization-capture amplification after size selection. That chemistry
introduces a random N7 segment between the Illumina sequencing-primer/index
side of the library molecule and the TSO, so the demultiplexed gRNA R1 begins
with N7 and the intact TSO normally starts at offset 7. The instrument's i7 and
i5 reads remain separate and do not appear in demultiplexed R1; no additional
N3 prefix is assumed. GEX carry-over escaped size selection before this
dedicated amplification and therefore lacks the circularization-capture N7.
The GEX classifier does not assume it.

Genuine short cDNA can also retain all or part of the TSO. TSO detection alone
is therefore not evidence of gRNA contamination. The headline truth definition
requires the downstream spacer-sized segment, paired CSM stem, 7-nt iBAR-sized
block, and paired partial tracr structure.

### Relevant split-pipe behavior

The local Parse pipeline is v1.6.4; pilot logs report STAR v2.7.11b.

- split-pipe rejects raw GEX R1 shorter than 45 nt and chemistry-v3 R2 shorter
  than 58 nt. Barcode correction permits a maximum edit distance of 2; this is
  independent of genomic-alignment mismatch handling.
- Before STAR, split-pipe searches R1 offsets 0-15 for the exact 17-mer
  `AACGCAGAGTGAATGGG`, bases 14-30 of the full 30-nt CSM Split_TSO, and removes
  sequence through the end of the match.
- For a raw 54-nt contaminant with an intact TSO at offset 0, exact trimming
  leaves approximately 24 nt: a 20-nt spacer plus four stem bases, or a 21-nt
  spacer plus three stem bases. Offsets 1-2 leave one or two fewer stem bases.
- split-pipe does not otherwise remove spacer, CSM stem, iBAR, or tracr sequence
  before alignment. STAR is also passed
  `--clip3pAdapterSeq GAGGTGGTTGGA`.
- The pilot command sets `--outFilterMultimapNmax 3`; other mapping filters use
  the local STAR defaults. STAR may emit reads mapping to up to three loci, but
  unique alignments have MAPQ 255 and split-pipe assigns genes only to
  non-secondary MAPQ-255 records. Multimappers may remain in BAM but do not
  become counted split-pipe molecules.

This makes exact TSO trimming important: an untrimmed 54-nt read needs much more
matched sequence to pass STAR's score/matched-length fractions than the
22-24-nt sequence passed to STAR after trimming. The audit therefore evaluates
the actual post-trim sequence state in both full-length and synthetic short-R1
arms.

Local implementation evidence is in:

- `ParseBiosciences-Pipeline.1.6.4/splitpipe/pb_const.py`
- `ParseBiosciences-Pipeline.1.6.4/splitpipe/preprocess.py`
- `ParseBiosciences-Pipeline.1.6.4/splitpipe/align.py`
- `ParseBiosciences-Pipeline.1.6.4/splitpipe/molinfo.py`
- `output/gex_C1/process/STAR_align_sdtout.txt`
- `output/gex_C1/process/barcode_headLog.out`

## Quick Start

Install the production additions into the same `splitpipe` environment used by
the main pipeline:

```bash
conda activate splitpipe
python -m pip install -r requirements_extra.txt
python -c 'import edlib, matplotlib, pandas, pysam, splitpipe, yaml'
STAR --version
```

`environment_qc.yml` creates the separate optional `csm_analysis` environment
for downstream analysis modules and notebooks. It is not a production pipeline
runtime. On Broad, `use Samtools` exposes the optional command-line `samtools`
module; production Python jobs still require `pysam` in `splitpipe`.

On a Broad login node, initialize the scheduler once per shell before using the
launcher:

```bash
source /broad/software/scripts/useuse
reuse -q UGER || true
command -v qsub
```

Copy `examples/qc/qc_config_gex_grna_contamination.example.yaml` to a run
configuration and set `base_config` to your pipeline YAML. Pass that file
explicitly through `--config` (launcher) or `CSM_GEX_CONTAM_CONFIG` (worker). Every production
launcher activates `environment.conda_env` from the base config, either the
named `splitpipe` environment or the absolute shared `splitpipe-share` prefix;
it never activates `csm_analysis`. The optional
`CSM_GEX_CONTAM_ENV` submission variable is an explicit testing override.
Conda's root is always read from `environment.conda_prefix`, with no
user-specific path fallback. Manual activation before submission is not
required after the one-time installation. The launcher also exports
`CSM_PROJECT_DIR` explicitly because UGER executes a spooled copy of each job
script.

The production FASTQ scan defaults to `auto`: one requested thread selects the
serial backend and more than one selects bounded multiprocessing. `--threads`
controls worker count, `--execution-backend serial|parallel|auto` permits an
explicit comparison, and `--batch-size` overrides
`performance.batch_size` (default 20,000 records). The parent keeps no more
than `2 × workers` batches, or `2 × workers × batch_size` records, pending.
UGER wrappers derive workers from `NSLOTS` and bind the same number of physical
cores as their `-pe smp` request.

The pipeline launcher uses eight raw-scan workers by default and passes the
bound `NSLOTS` value to the QC engine. This is the recommended starting point
for full, large GEX scans when the queue and memory budget permit. The completed four-worker
performance screen showed essentially linear scaling from one to four workers;
do a bounded 4-versus-8 comparison before requesting more than eight because
gzip decoding and parent aggregation eventually limit scaling.

The parallel backend uses the production packed exact/three-seed distance-2
guide index and a validated binary paired-FASTQ parser. Workers return aggregate
truth, filter, barcode, HLL/molecule, and diagnostic state plus bounded artifact
selections. The parent writes `read_evidence.tsv.gz`, synthetic paired FASTQs,
and controlled-alignment panels in original input order, so selecting the
parallel backend does not remove or reorder existing QC artifacts.

From the HPC project directory, first inspect the generated qsub commands:

```bash
cd /medpop/rgupta/Jiahao/my_project/CSM_ccperturbseq
bash submit/submit_qc_gex_grna_contamination_pipeline.sh \
  --config qc_config_gex_grna_contamination.yaml \
  --mode all \
  --workers 8 \
  --max-concurrent 2 \
  --dry-run
```

Then submit the complete audit:

```bash
bash submit/submit_qc_gex_grna_contamination_pipeline.sh \
  --config qc_config_gex_grna_contamination.yaml \
  --mode all \
  --workers 8 \
  --max-concurrent 2 \
  --alignment-max-concurrent 2
```

`--mode all` and `--mode benchmark` create a three-stage dependency graph:

1. one raw FASTQ benchmark/BAM-audit task per GEX sublibrary;
2. a held `library x {full, short54}` controlled STAR array; and
3. a held aggregation job after all alignment audits finish.

Other useful modes are:

| Mode or option | Use |
|---|---|
| `--mode bam` | Re-run only the historical annotated-BAM audit. |
| `--mode guide` | Compare strict Step 1c evidence with the dedicated guide arm; requires `guide_assignment_impact.enabled: true`. |
| `--mode alignment` | Re-run the controlled STAR impact stage from an existing complete truth panel. |
| `--mode aggregate` | Rebuild cross-library tables, plots, and the generated report. |
| `--mode filter --no-aggregate` | Materialize full-library paired clean/rejected 54-nt FASTQs after a filter has been accepted. |
| `--no-alignment --no-aggregate` | Fast code/data smoke test; it does not answer the alignment-impact question. |

For a smoke test, temporarily cap `raw_benchmark.max_reads` and
`bam_audit.max_primary_reads`. Such output is deliberately marked incomplete
and must not be reported as prevalence. For the final audit, use
`raw_benchmark.sample_fraction: 1.0`, `raw_benchmark.max_reads: 0`, and
`bam_audit.max_primary_reads: 0`; the controlled truth-panel safety cap must not
be reached.

> [!IMPORTANT]
> Begin every review with `audit_completeness.csv` and interpret a library only
> when `audit_complete` is true. Raw FASTQ, controlled STAR, and annotated-BAM
> fractions use different denominators and are not interchangeable estimates.

The fastest review path is:

1. verify completeness and provenance;
2. report raw-read truth tiers and pre-alignment UMI-family bounds;
3. require the hard-negative false-positive-rate gate for any 54-nt filter;
4. compare the controlled full and short54 STAR arms; and
5. use the historical BAM audit to quantify impact already visible in the
   pilot result.

## Detailed Step-by-Step QC Workflow

### Step 0 — Validate inputs, provenance, and audit completeness

**What it does.** The workflow records FASTQ/BAM metadata, the QC module digest,
the relevant configuration digest, truth-panel token, expected record counts,
and scan-completion state. Aggregation rejects stale or incomplete evidence
instead of silently combining it.

**Expected output.** The primary cross-library gate is:

```text
output/qc_gex_grna_contamination/audit_completeness.csv
```

Important columns include `benchmark_complete`, `benchmark_status`,
`truth_safety_cap_reached`, `bam_audit_complete`,
`alignment_full_complete`, `alignment_short54_complete`, and
`audit_complete`.

**How to interpret it.** Require `audit_complete=true`. `benchmark_status` must
be `current`; module/configuration digests and input file metadata must still
match; both alignment arms must have the same truth token and expected read
count; and
`truth_safety_cap_reached` must be false. A partial scan, changed input,
digest mismatch, missing arm, or capped truth panel invalidates downstream
prevalence and comparison results. UGER job holds enforce order, not successful
completion, which is why this file is the first review gate.

The historical BAM is optional and disabled by default. When
`bam_audit.enabled: true`, its status is still reported. It gates
`audit_complete` only when `bam_audit.required_for_audit: true`; use that stricter
setting only for a BAM that retained unmapped records or for a specifically
defined mapped-survivor question.

**v1.4 release evidence.** Complete raw scans are available for C1 and C2. The
corrected full C2 production scan used eight bound workers and generated the
all-contaminant controlled panel. Both full and short54 STAR arms completed on
all 10,559 corrected-barcode structural contaminants. The historical mapped-BAM
scan is supplementary: split-pipe normally discards STAR-unmapped reads, so
that BAM cannot estimate raw contamination and is not required to interpret
the controlled experiment.

### Step 1 — Establish full-length contaminant truth and raw-read prevalence

**What it does.** The original paired pilot FASTQs are scanned before any R1
truncation or split-pipe TSO trimming. TSO-only reads are retained as a
diagnostic/hard-negative population, never as contaminant truth. Four nested
confidence tiers are reported:

1. `definitive_reference_construct`: TSO, complete paired structure, known
   spacer, and compatible reference iBAR.
2. `reference_supported_structure`: TSO, complete structure, and known spacer.
3. `primary_structural_surrogate`: TSO, spacer-sized segment, complete paired
   stem, 7-nt iBAR-sized block, and paired partial tracr. Guide/iBAR identity
   may be unknown or recombined. This is the headline tier.
4. `inclusive_possible_contaminant`: the primary structure or TSO plus a known
   spacer when downstream structure is incomplete. This is the
   sensitivity-oriented upper tier.

**Expected outputs.** Per library:

- `benchmark_summary.yaml`
- `read_evidence.tsv.gz`
- `contaminant_guide_counts.csv`
- `contaminant_target_counts.csv`
- `contaminant_tso_offsets.csv`
- `per_cycle_base_composition.csv`

In `benchmark_summary.yaml`, the relevant sections are
`raw_read_prevalence`, `full_length_surrogate_truth.confidence_tiers`, and scan
provenance/completeness fields.

**How to interpret them.** `all_reads` is the denominator for the raw-library
prevalence. `tso_reads` includes genuine short cDNA and must not be called
contamination. `structure_contaminant_reads` is the primary numerator;
`contaminant_fraction_of_all_reads` and its Wilson 95% interval are the headline
read-level estimate. `contaminant_fraction_of_tso_reads` explains how unsafe a
TSO-only removal rule would be.

Report the primary tier together with the definitive and inclusive tiers as a
plausible evidence range. A large separation between tiers warrants review of
`read_evidence.tsv.gz`, guide/iBAR compatibility, and base composition. True GEX
carry-over should be concentrated at full-TSO offset 0 with, at most, a small
tail at damage/error-tolerant offsets 1-2. Require `raw_scan_complete=true`; a
fractional or record-limited development run is not the final prevalence audit.

The complete pilot scans found a low but reproducible signal:

| Library | Read pairs | Definitive (CPM) | Primary structural (CPM, 95% CI) | Inclusive (CPM) |
|---|---:|---:|---:|---:|
| C1 | 452,146,811 | 33.08 | 36.61 (36.05-37.17) | 39.25 |
| C2 | 363,646,678 | 29.55 | 32.13 (31.55-32.72) | 34.80 |

Thus the observed carry-over is approximately 30-40 reads per million, not a
large fraction of the GEX library. It can still matter biologically: an
incorrect STAR assignment can rival or exceed the true count of a lowly
expressed target gene, particularly with a 100-cycle read that ends within the
guide-derived insert.

### Step 2 — Estimate cell-associated reads and pre-alignment UMI families

**What it does.** R2 is decoded using the chemistry-v3 layout and corrected
against retained `barcode_data.csv` tables using split-pipe-compatible
calibration: perfect combinations from logical records 3,000,000-3,999,999 in
the configured ordered R2 stream (continuing across lane boundaries), a cumulative
92% cell threshold with minimum count 5, and a uniquely
supported correction at Hamming distance 0-2. Both corrected cell identity
(`CB`) and barcode-index identity (`pB`) are retained because different
round-1 barcodes can map to the same well.

Distinct pre-alignment families are estimated with HyperLogLog at precision 18
(nominal relative error approximately 0.20%) using two keys:

- `cell_pb_umi = (CB, pB, exact UMI)`, a lower molecular-diversity bound that
  can merge UMI collisions across inserts; and
- `cell_pb_umi_insert = (CB, pB, exact UMI, full raw R1)`, an upper bound that
  separates insert collisions but can split one PCR family at R1 sequencing
  errors.

**Expected outputs.** `benchmark_summary.yaml` contains
`prealignment_molecule_prevalence`, including `valid_barcode_read_counts`,
`estimated_distinct_keys`, fractions, and nominal sketch intervals.
`cell_read_contamination.csv` provides the corresponding cell-level read
summary.

**How to interpret them.** Use valid-barcode counts to determine how much raw
contamination can be associated with corrected cells. Report
`fractions.cell_pb_umi.contaminant_of_all` and
`fractions.cell_pb_umi_insert.contaminant_of_all` together as lower/upper
pre-alignment family estimates; the analogous `contaminant_of_tso` fields use
TSO-associated families as their denominator. The nominal interval represents
sketch error only, not uncertainty in structural truth or barcode correction.

These are not final split-pipe molecules or a literal count of cDNAs before PCR.
Gene identity is unavailable before alignment, exact-key collisions/errors can
merge or split families, and split-pipe later performs abundance-directed
Hamming-1 UMI collapse within assigned genes.

For the complete pilot scans, the cell-associated and family estimates were:

| Library | Structural reads with valid barcodes | Valid-barcode reads | Pre-alignment contaminant-family range (CPM of all families) |
|---|---:|---:|---:|
| C1 | 14,643 | 300,115,343 | 48.05-75.62 |
| C2 | 10,559 | 244,524,486 | 42.47-69.66 |

The lower endpoint includes the full raw R1 in the family key and the upper
endpoint does not; use the pair as bounds rather than choosing one as a final
molecule count.

### Step 3 — Benchmark 54-nt contaminant-removal rules in silico

**What it does.** R1 sequence and quality are truncated to 54 nt in memory; R2
is unchanged. Each candidate rule is compared with the full-length truth from
Step 1. The rule panel includes:

- TSO-only detection, as an unsafe diagnostic rather than a production option;
- TSO plus exact/fuzzy full spacer;
- TSO plus the expected `spacer + stem` prefix;
- TSO plus a complete spacer and at least one compatible stem base
  (`tso_template_boundary_d0/d1/d2`);
- exact 12-, 15-, and 18-nt guide-prefix rules;
- simulated split-pipe exact-17-mer trimming followed by guide matching; and
- wider-offset template detection as a diagnostic.

Primary rules expect the TSO at offset 0 and allow offsets 1-2 only for
low-frequency damage/error rescue.

**Expected outputs.** Per library:

- `short54_filter_metrics.csv`
- `short54_metrics_by_tso_offset.csv`
- `benchmark_summary.yaml` sections `short_read_hazard`, `selected_filter`, and
  `recommended_filter_from_benchmark`

**How to interpret them.** Start with rows where
`truth_definition == structure_truth`.

| Metric | Interpretation |
|---|---|
| `sensitivity` | Fraction of full-length structural contaminants removed at 54 nt. |
| `specificity` | Fraction of all structural negatives retained. |
| `precision` | Fraction of calls that are truth-positive; this depends on prevalence. |
| `false_positive_rate` | False calls among all structural negatives. |
| `hard_negative_fpr` | False calls among full-TSO, structure-negative GEX-like reads; the key cDNA-safety metric. |
| `hard_negative_fpr_95_high` | Upper Wilson bound used as the primary safety gate. |
| `f1` | Sensitivity/precision summary; it never replaces the hard-negative gate. |

By default a non-TSO-only rule is recommendation-eligible only when
`hard_negative_fpr_95_high <= 0.001`. Even zero observed false positives has a
nonzero upper bound when the negative panel is small. Eligible rules are ranked
by sensitivity, then precision, then global FPR.

`selected_filter` is copied from the QC configuration and generated the example
FASTQs; it is not automatically accepted. An empty
`recommended_filter_from_benchmark` means no rule passed the gate. Inspect
offset 0 first, then offsets 1-2 as rescue strata, and require acceptable
performance for both 20- and 21-nt spacers and both construct variants. TSO-only
or broad-offset matching should not be promoted merely because it is sensitive.

In both pilot libraries, the selected template-boundary diagnostic achieved
about 82% sensitivity and 95% precision against structural truth. Its
hard-negative FPR upper bound was 14.26 CPM for C1 and 8.80 CPM for C2. The
more conservative production-like `tso_guide_d1` rule had 82.84%/82.50%
sensitivity and 94.64%/94.33% precision for C1/C2, with hard-negative upper
bounds of 16.67 and 10.88 CPM. By contrast, a TSO-only rule had approximately
0.03% precision in C1 and 0.018% in C2 because it classified nearly every
TSO-positive hard negative as contamination. Never use TSO alone.

The QC ranker preferred a distance-2 guide rule for sensitivity, but production
v1.4 deliberately retains distance 1 as the conservative default. On full C2,
that production rule removed 10,267 pairs (28.23 CPM), including 9,511 pairs
with valid corrected barcodes, while preserving the validated removal and
survivor invariants.

### Step 4 — Generate synthetic 54-nt and optional filtered paired FASTQs

**What it does.** A deterministic hash sample distributed across input lanes is
written in three paired forms: `all`, `clean`, and `rejected`. R1 sequence and
quality are truncated together to 54 nt; R2 remains unchanged. Removing an R1
always removes its paired R2, and pair IDs and record counts are validated.
After a rule passes Step 3, `--mode filter` can apply the same operation to the
whole library.

**Expected outputs.** Per library:

```text
synthetic_54/   # sampled all, clean, and rejected paired FASTQs
filtered_54/    # full-library output, present only after --mode filter
```

Paths and counts are also recorded in `benchmark_summary.yaml`.

**How to interpret them.** These files are reproducible test material, not
independent evidence that a filter is safe. Verify synchronized R1/R2 record
counts, 54-nt R1 lengths, unchanged R2, and `all = clean + rejected`. Only
materialize full-library output after the chosen rule passes the hard-negative
gate and stratified sensitivity review. An optional split-pipe comparison of the
synthetic `all` and `clean` pairs provides an end-to-end check before production
adoption.

For full C2, the four-worker implementation reproduced the serial result and
all 17 detailed artifacts, including the sampled `all`, `clean`, and `rejected`
paired FASTQs. This establishes deterministic artifact generation across the
two backends for this input; it is not a substitute for the biological safety
checks in Step 3.

### Step 5 — Measure false alignment in controlled full and short54 STAR arms

**What it does.** Full-length structural truth is retained before truncation.
The release configuration includes every corrected-barcode structural
contaminant (`sample_fraction: 1.0`) and keeps PCR-family abundance and UMI
collapse auditable. It writes two single-end STAR inputs with the same truth
IDs:

1. full pilot R1 after split-pipe's exact 17-mer TSO preprocessing; and
2. raw R1 truncated to 54 nt, followed by the same preprocessing.

**Both arms align to the same `paths.genome_dir` and use identical STAR
parameters**, including `--outFilterMultimapNmax 3` and the 3-prime adapter.
The only arm-specific transformation is that `short54` truncates raw R1
sequence and quality to 54 nt before the same Parse-equivalent TSO
preprocessing. Every input read
is contaminant truth, so every gene assignment in these arms is a false GEX
assignment by construction. The module reproduces split-pipe's directional
Hamming-1 collapse for input `(CB,pB,guide)` families and false output
`(CB,pB,assigned gene)` molecules.

**Expected outputs.** Per library:

```text
alignment_input/
├── ground_truth_manifest.tsv.gz
├── full/<library>_ground_truth_R1_full.fastq.gz
└── short54/<library>_ground_truth_R1_54.fastq.gz

alignment_star/full/
├── <library>_Aligned.out.bam
├── alignment_impact_summary.yaml
├── false_gene_assignments.csv
├── cell_false_molecule_impact.csv
└── cb_pb_cluster_alignment_outcomes.csv

alignment_star/short54/
└── the same five output types
```

The cross-library comparison is
`all_libraries_full_vs_short54_alignment_delta.csv`.

**How to interpret them.** In each `alignment_impact_summary.yaml`, report
mapped, multimapping, unique MAPQ-255, false-gene-assigned, and guide-target
fractions. `false_gene_assignment_fraction_of_input_reads` is the conditional
probability that a known contaminant becomes a false gene-expression read;
`guide_target_assignment_fraction_of_input_reads` is the target-concordant
subset. Prefer the corrected-`(CB,pB)` cluster-bootstrap 95% intervals; naive
read-level Wilson intervals do not account for dependence within PCR families.

Use `splitpipe_equivalent_hamming1_molecule_impact` for the controlled molecule
result and `outcomes_by_predicted_tso_trim` to separate exact-TSO-trimmed from
untrimmed inputs. In the paired delta table, an absolute delta above zero means
54-nt R1 increased false assignment. Treat the change as supported only when
the paired cluster-bootstrap interval excludes zero; fold change is unstable
when the full-length rate is near zero.

The controlled rate is conditional on contaminant truth and is not the library
prevalence. Aggregate burden projections multiply this rate by raw prevalence,
but mixing with genuine GEX can merge or change gene/UMI clusters. The
historical BAM in Step 6 measures that competition directly for the current
pilot. If the truth safety cap is reached, raise the cap or intentionally define
a new sampled analysis; never report a silently split/capped panel as the
all-read C2 census. The corrected panel size must match the current
split-pipe-compatible barcode decoder rather than an earlier Hamming-only
approximation.

**C2 controlled-alignment result.** The corrected release panel contains
10,559 reads in each arm. UGER array job `60030670` aligned both arms with STAR
2.7.11b, eight linearly bound cores per task, the same genome, and the same
options. The only experimental difference was truncation to 54 nt before
Parse-equivalent TSO preprocessing.

| Outcome | Full length | Short54 |
|---|---:|---:|
| Input structural contaminants | 10,559 | 10,559 |
| Mapped reads | 0 | 8,113 (76.83%) |
| Unique MAPQ-255 reads | 0 | 6,965 (65.96%) |
| Multimapping primary reads | 0 | 1,148 (10.87%) |
| False gene-assigned reads | 0 | 4,008 (37.96%) |
| Intended-target gene assignments | 0 | 2,431 (23.02% of inputs; 60.65% of false assignments) |
| Split-pipe-equivalent false gene molecules | 0 | 3,099/8,238 (37.62%) |
| Intended-target false molecules | 0 | 1,894/8,238 (22.99%; 61.12% of false molecules) |

The corrected-`(CB,pB)` cluster-bootstrap 95% interval was 36.09-39.88% for
false assigned reads and 35.88-39.41% for false molecules. All 8,113 mapped
short54 reads came from the 9,551 reads predicted to receive TSO trimming. That
trimming left mostly 23-27 nt of spacer plus a few stem bases; 84.94% of this
trimmed stratum mapped. In contrast, the longer post-trim inputs in the full
arm were all reported by STAR as too short to map relative to their total read
length. This is consistent with the downstream stem-iBAR-tracr context
protecting the 150-cycle library from spurious alignment.

The production v1.4 decision would remove 3,833/4,008 (95.63%) of the short54
false assigned reads and 2,963/3,099 (95.61%) of the corresponding collapsed
molecules. The retained stratum contained 175 false assigned reads and 136
false molecules, which is why the rule remains conservative rather than being
described as complete removal.

Post-STAR annotation job `60030674` normalized the literal quote characters in
Parse's `gene_info.json`. It found 2,431 target-concordant reads and 1,894
target-concordant collapsed molecules. Every one was in the production removal
stratum; none of the 175 retained false assigned reads was target-concordant.
Major target-concordant assignments included KLF2, NFKB1, FDX1, LAMB2,
ADAMTS7, KLF4, CEBPB, JUN, and CCM2. This is the biologically important result:
short guide-derived reads can create false expression precisely at the intended
perturbation target, while the v1.4 filter captured every such event observed
in this C2 panel.

A read- and locus-level follow-up for the two strongest complementary examples
is available in the [C2 KLF2/NFKB1 misalignment case study](C2_GEX_GRNA_MISALIGNMENT_CASE_STUDY.md).

### Step 6 — Audit contaminant impact in existing split-pipe BAMs

**What it does.** The module scans the existing annotated BAM, restores
reverse-strand records to sequencing orientation, detects complete constructs
whether the TSO was retained or trimmed, and reconstructs split-pipe-style
gene/UMI groups using abundance-directed Hamming-1 collapse. It reports
alignment uniqueness, gene assignment, guide-target concordance, and molecules
containing any, all, or a majority of contaminant reads.

**Expected outputs.** Per library:

- `bam_audit_summary.yaml`
- `bam_contaminant_gene_assignments.csv`
- `bam_guide_target_concordance.csv`
- `bam_contaminant_mapq_counts.csv`
- `bam_cell_molecule_contamination.csv`

**How to interpret read-level results.** `primary_reads` is the BAM denominator,
not the raw FASTQ denominator. The annotated BAM normally excludes unmapped
reads, so `fraction_of_primary_records` and
`fraction_of_mapped_primary_reads` describe visible alignment impact, not raw
contamination prevalence. `unique_contaminants` uses MAPQ 255;
`gene_assigned_contaminants` requires split-pipe-compatible `GX`/`GN` tags; and
`assigned_to_guide_target` is the most concerning target-concordant subset.
Review the assignment and concordance tables because a promoter guide need not
be assigned to its intended target, while another unique genic hit may be.

**How to interpret molecule results.** In `counted_molecule_impact`:

- `reconstructed_all_molecules` is the audit denominator.
- `vendor_number_of_tscp` is split-pipe's reported transcript count when
  available. `reconstruction_difference` should be small for a complete scan.
  A tiny discrepancy can reflect the v1.6.4 final-chunk `[:-1]` EOF quirk; a
  large discrepancy indicates an incompatible BAM/tag/order assumption.
- `molecules_with_any_guide_read_upper_bound` includes mixed genuine/guide
  clusters and is the upper bound.
- `molecules_with_all_guide_reads_lower_bound` is the conservative lower bound.
- `molecules_with_majority_guide_reads` is a useful middle summary, not a new
  truth definition.
- `guide_target_molecules_with_any_guide_read` is the expected-target subset.

Report the any/all interval rather than selecting one endpoint when the bounds
are wide. The BAM post-trim classifier is a cross-check; the controlled arms are
the primary full-versus-short ground-truth comparison because they retain reads
that disappear before the historical annotated BAM.

The C2 historical-BAM audit job `60030663` was stopped after 115 million primary
records with zero detected contaminants. It ran under the superseded analysis
environment and is not release evidence. More importantly, it is not a
substitute for the controlled arms: the source BAM normally omits STAR-unmapped
reads, so even a complete zero count would mean only that no recognizable
contaminant survived into that mapped-only artifact. The job is retained as a
supplementary compatibility observation and does not gate v1.4.

### Step 7 — Recover conservative Step 1c guide/iBAR evidence

**What it does.** Production Step 1a emits one compact row for each
filter-positive read with a valid split-pipe barcode. Each row contains the
full-run corrected cell barcode (`CB`), barcode-index identity (`pB`), UMI, R1
sequence, and R1 quality. Step 1c runs ordinary Step 3 directly on this table;
it does not scan the GEX FASTQs again and does not recalibrate barcodes from
the small rejected subset.

Ordinary Step 3 matches remain in
`<library>_grna_matches_unfiltered.tsv.gz`. The matrix input is a conservative
subset: guide edit distance at most one, iBAR edit distance at most one, and a
resolved guide-compatible corrected iBAR in dual-guide mode; decoupled
guide/iBAR calls are excluded. This is intentionally stricter
than the dedicated-arm Step 3 limits because Step 1c is optional and can only
add guide UMIs.

**Expected outputs.** Step 1a and Step 1c provide:

```text
output/gex_grna_filter/<library>/valid_positive_barcode_evidence.tsv.gz
output/gex_grna_identity/grna_matches/<library>_grna_matches_unfiltered.tsv.gz
output/gex_grna_identity/grna_matches/<library>_grna_matches.tsv.gz
output/gex_grna_identity/manifests/<library>.step1c.complete.json
```

**How to interpret them.** Require the Step 1a schema-4 evidence count and hash
to match the Step 1c completion manifest. For C2 the release gate expects
9,511 evidence rows from 10,267 filter-positive pairs. A guide-only match
without the expected iBAR is useful diagnostic evidence but must not be used to
augment a dual-guide identity matrix.

**C2 release result.** Step 1c read exactly the 9,511 Step 1a evidence rows; it
did not reopen the 363,646,678-pair GEX input. Ordinary Step 3 extracted 9,070
constructs and matched 9,060 (8,190 exact, 866 fuzzy, and 4 guide-only). The
strict dual-guide release policy retained 9,014 rows. It excluded 3 guide
distance failures, 11 iBAR distance failures, 4 unresolved iBARs, and 28
guide/iBAR-decoupled calls. The one-core UGER task completed in about 5.4
seconds. These are read-evidence counts before Step 4 depth filtering and
directional UMI collapse, not added guide-matrix UMIs.

The guide-impact job also fits the ordinary full-GEX 3–4M barcode calibration
interval, then applies the independent QC decoder only to the 10,267 stored
Step 1a positive pairs. It requires every field of all 9,511 valid calls—read
ID, corrected cell barcode, pB, UMI, R1 sequence, and quality—to equal the Step
1a evidence in order. This is a small release check, not another full-library
scan.

### Step 8 — Quantify dedicated-versus-GEX guide-matrix impact

**What it does.** `--mode guide` passes dedicated and strict GEX matches through
the production Step 4 accumulator. Reads from both sources are combined before
`min_reads_per_umi` filtering and directional Hamming-1 UMI collapse. It keeps
dedicated-only, GEX-only, and combined counts, then applies the production
single- or dual-guide assignment threshold to identify cells whose assignments
are added, removed, changed, or unchanged.

The exact cross-source molecular key is `(sublibrary, corrected cell barcode,
guide sequence, corrected iBAR, UMI)`. A stricter supporting key also requires
the same pB. When GEX R1 is too short to resolve iBAR, guide-only overlap is
reported separately and never substituted for the construct-aware key.

**Expected outputs.** Per library:

```text
guide_assignment_impact/guide_assignment_impact_summary.yaml
guide_assignment_impact/combined_umi_counts_by_source.csv.gz
guide_assignment_impact/per_cell_guide_umi_delta.csv
guide_assignment_impact/per_guide_umi_delta.csv
guide_assignment_impact/cell_assignment_changes.csv
guide_assignment_impact/cross_source_molecule_overlap.csv
```

**How to interpret them.** Report exact and collapsed UMI increments, the
number of changed cells, and both overlap definitions. Exact overlap is
evidence—not proof—that the same RT molecule was amplified into both
size-selected libraries, because cell-barcode/UMI collisions remain possible.
pB agreement strengthens that evidence. A zero overlap does not invalidate
GEX carry-over; the GEX contaminant may not have entered or survived the
separate circularization-capture amplification used for the dedicated library.

**C2 production-merge result.** The 9,014 strict GEX reads contained 6,847
exact UMIs and produced 6,788 collapsed GEX UMI clusters. C2 dedicated reads
contained 1,000,411 exact UMIs and produced 789,746 collapsed clusters. The
combined C2 stream produced 790,439 collapsed clusters, only 693 more than the
dedicated-only stream. Exact `(CB,guide,iBAR,UMI)` cross-source overlap was
5,995 and all 5,995 also agreed on pB. At the configured pair-assignment
threshold, 7,619 C2 cells were unchanged, two changed assignment, and four
gained an assignment.

The isolated final-AnnData integration test retained all four dedicated
libraries while adding the C2 GEX source. Across that four-library matrix,
the dedicated-only total was reproduced exactly at 3,272,782 UMIs and the
combined total was 3,273,475. Assigned cells increased from 38,325 to 38,329
and assignment rows from 45,759 to 45,765. The earlier 6,041/6,098 values are
collapsed source-union cluster diagnostics, not the stricter exact/pB-aware
cross-source key counts.

These small net changes support keeping Step 1c opt-in and conservative. Exact
overlap is strong evidence of shared molecules across the two size-selected
products, but it is not proof because barcode/UMI collisions are possible.

### Step 9 — Aggregate libraries and make the decision

**What it does.** Once all required per-library stages pass completeness checks,
the held aggregation job combines evidence without mixing denominators and
creates a concise report and figures.

**Expected outputs.** At `output/qc_gex_grna_contamination/`:

- `all_libraries_contamination_summary.csv`
- `all_libraries_short54_filter_metrics.csv`
- `all_libraries_short54_metrics_by_tso_offset.csv`
- `all_libraries_bam_audit_summary.csv`
- `all_libraries_alignment_impact.csv`
- `all_libraries_full_vs_short54_alignment_delta.csv`
- `audit_completeness.csv`
- `GEX_GRNA_CONTAMINATION_AUDIT_REPORT.md`
- `plots/contamination_fraction.png`
- `plots/short_read_filter_tradeoff.png`
- `plots/controlled_alignment_false_assignment.png`

**How to interpret them.** The generated report is the quickest overview, but
the YAML/CSV files remain authoritative for denominators, provenance, intervals,
and strata. Review in this order:

1. Require `audit_complete=true` for every compared library.
2. Report primary, definitive, and inclusive raw-read prevalence with intervals.
3. Report the `cell_pb_umi` and `cell_pb_umi_insert` pre-alignment family bounds.
4. Accept a short filter only if the hard-negative FPR upper bound passes the
   predeclared gate and offset-0, guide-length, and variant sensitivity is
   adequate.
5. Compare full versus short54 false-read and Hamming-1 molecule rates using
   paired corrected-`(CB,pB)` cluster-bootstrap intervals.
6. Report historical-BAM any/all molecule bounds, affected genes, and
   guide-target concordance.
7. Keep `selected_filter` unchanged if no method passes. Revise the assay or
   collect more validation negatives instead of relaxing the gate post hoc.

For experiments using the same universal RT, size-selection, and
circularization-capture design, the C1/C2 evidence supports enabling
`gex_decontamination.enabled: true` when the complete guide reference is
available. The signal is reproducible, the production rule is intentionally
conservative, and prep correctness/performance passed the release gates. Run
this QC first or leave decontamination disabled for different chemistries,
incomplete guide references, or experiments without a shared RT product and
size-selection opportunity. The v1.4 release summary covers C1/C2 raw evidence
and the completed C2 controlled full/short54 experiment; it does not claim
T1/T2 or a four-library controlled-alignment result.

### Step 10 — Validate the module before reporting results

**What it does.** Unit and miniature integration tests exercise reference
constructs at TSO offsets 0-2, both stem variants, 20/21-nt spacers, the
21-mer-ending-G edge case, shared spacers, TSO/spacer errors, TSO-positive cDNA
negatives, BAM post-trim structure, paired FASTQ integrity, and benchmark/filter
execution.

**Expected output.** Run on the HPC environment:

```bash
conda activate splitpipe
python3 -m unittest discover -s tests -v
```

The command should finish with all tests passing and no paired-FASTQ integrity
errors.

The miniature production test compares serial and four-worker outputs by
decompressed content, including per-read tables, synthetic R1/R2 files,
barcode/molecule summaries, and full/short controlled-alignment panels. The
generality regression executes complete 1-, 4-, and 12-library aggregates,
multi-lane and explicit-pair resolution, overlapping prefixes, duplicate/shared
input rejection, and 20k/50k-guide packed-index checks.

**How to interpret it.** Passing tests establish software consistency for the
covered cases; they do not validate the biological truth definition or approve
a production filter. The pilot evidence, completeness gates, hard-negative FPR,
controlled STAR comparison, and historical-BAM audit remain the acceptance
criteria.

The corrected full C2 release scan processed all 363,646,678 pairs with eight
bound workers at 68,089 pairs/s (5,340.7 seconds of scan time; 5,368.1 seconds
scheduler wall time). The earlier bounded one-versus-four-worker benchmark
measured a 4.019x scan speedup and reproduced all detailed artifacts. Eight
workers is therefore the recommended full-scan default for this scale;
benchmark before requesting more because gzip decoding and parent-process
aggregation eventually limit scaling.

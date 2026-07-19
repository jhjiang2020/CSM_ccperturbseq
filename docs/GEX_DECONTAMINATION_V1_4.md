# v1.4 GEX gRNA decontamination

v1.4 removes reference-supported gRNA contamination during split-pipe
preprocessing, before the human STAR input FASTQ is finalized. Guide identity
is deliberately handled in a separate optional Step 1c.

The contamination originates upstream of sequencing. Universal reverse
transcription uses both poly(A) probes and gRNA probes, producing long
gene-transcript cDNA and short gRNA-derived products together. Size selection
separates the GEX and gRNA fractions, but a small amount of the short product
can remain in GEX. The completed C1/C2 pilot scans measured approximately
32–37 structural-contaminant reads per million. This low frequency can still
matter when a guide-derived read is assigned by STAR to a target gene whose
true expression is at the same scale or lower, particularly with short
100-cycle read recipes.

The completed C2 controlled experiment makes that short-read risk quantitative.
None of 10,559 full-length structural contaminants mapped, whereas after raw R1
was truncated to 54 nt and passed through the same TSO preprocessing, 8,113
mapped and 4,008 (37.96%) received a false gene assignment. Directional
Hamming-1 collapse yielded 3,099 false gene molecules. The production v1.4
decision removed 95.63% of those false assigned reads while retaining a
conservative rule for protecting genuine TSO-bearing cDNA.
After normalizing Parse gene-name quoting, 2,431 false reads and 1,894 false
molecules were assigned to their own intended guide target. All of these
target-concordant events were removed by the production decision; the small
retained false-assignment stratum was off-target in this C2 panel.

The pilot QC rationale, execution workflow, and interpretation of its raw,
controlled-alignment, and historical-BAM summaries are documented in the
consolidated [GEX gRNA contamination QC guide](GEX_GRNA_CONTAMINATION_QC.md).

## Runtime overlay; no Parse source installation

The Parse 1.6.4 distribution remains pristine. When decontamination is enabled,
Step 1a starts `scripts/splitpipe_gex_wrapper.py`, which validates the exact
Parse version and source checksums, applies the tracked unified diff in memory,
and then executes the original `split-pipe` launcher. The wrapper intercepts
only `splitpipe.preprocess` and `splitpipe.molinfo`; it never writes either
source file or the installed package.

This works because the Parse launcher imports its processing modules only after
argument parsing. UGER/Linux uses the `fork` multiprocessing method, so the
preloaded packed guide index and overlaid preprocessing functions are inherited
by workers. A non-`fork` start method is rejected when the hook is enabled.

The checksum/context checks are intentionally strict. If Parse is upgraded,
locally edited, or still contains the earlier physical v1.4 patch, Step 1a
stops before processing and reports the unexpected checksum. Restore a pristine
Parse 1.6.4 tree rather than patching it. With
`gex_decontamination.enabled: false`, Step 1a invokes the original executable
directly and no overlay is installed.

## Fresh-run data flow

```text
raw paired GEX FASTQ
  → split-pipe R2 barcode correction
  → full TSO expected at offset 0; offsets 1–2 tolerate rare 5' damage/error
  → complete guide prefix matched against the reference
       exact first; otherwise verified Hamming distance 1
  → filter positive: omit from STAR; retain corrected-barcode evidence for Step 1c
       and, by default, synchronized raw R1/R2 for provenance and QC
  → filter negative with full TSO: trim at detected full/fuzzy TSO endpoint
  → filter negative without full TSO: retain Parse exact-17-mer fallback trim
  → surviving reads: barcode_head.fastq → human STAR
  → ordinary split-pipe post → mol → dge → sub → ana → rep → clean
```

Step 1a does not assign a guide name, iBAR, pair, target, cell, or construct.
Its anonymous sequence-group IDs are diagnostics only. This separation allows
removal to use a narrow, sensitive policy while Step 1c keeps the ordinary
gRNA-arm identity thresholds.

The packed index stores each full guide as two-bit DNA plus two complementary
seed tables. Any query with at most one substitution preserves at least one
seed; every candidate is then verified across the complete guide. It does not
materialize all Hamming-1 neighbors and is practical for 20–50k unique guides.
The implemented v1.4 backend is `packed_inline`; indels and mismatch distance
greater than one are rejected by configuration rather than silently changing
semantics.

The production router avoids constructing rich per-read decision objects in
the normal path. It precomputes the TSO pattern and pigeonhole blocks, uses
`edlib` for one constrained C-backed fuzzy alignment, and retains the
dependency-free Myers route as its executable correctness oracle. The common
guide-negative path uses exact strings and complementary half-seeds before any
packed Hamming verification. Classification, survivor trimming, and fixed
worker counters are fused so ordinary reads do not parse unused read IDs or R2
artifact fields. The full `classify()` result remains available when
`io.write_read_decisions` is enabled.

The canonical integration setting is:

```yaml
gex_decontamination:
  integration: splitpipe_pre_router
  tso_matching:
    backend: edlib
```

Install the tracked add-ons once in the configured splitpipe environment:

```bash
conda activate splitpipe
python -m pip install -r requirements_extra.txt
python -c 'import edlib, yaml'
```

`edlib>=1.3` is required when decontamination is enabled. Step 1a fails during
preload rather than silently using the much slower Python matcher. The master
launcher activates `environment.conda_env` inside every UGER job, so no manual
activation is needed at submission time after this one-time installation.

`splitpipe_pre_candidate_router` remains a deprecated read-only alias so an old
audit configuration fails visibly with a warning rather than changing routing
semantics. New configurations and the historical-BAM test use only
`splitpipe_pre_router`.

### Survivor TSO trimming

The removal decision and STAR-survivor trimming use one ordered contract:

```text
full TSO + reference-supported guide  → remove before STAR
full TSO + no supported guide         → trim at the full-TSO detector endpoint
no full TSO + Parse exact 17-mer      → trim with Parse's original rule
neither                               → leave R1 unchanged
```

Thus a genuine short cDNA with a fuzzy full TSO is retained but no longer sends
that TSO to STAR merely because an error occurred within Parse's terminal
17-mer. Partial-TSO cDNA retains Parse's historical exact-suffix behavior. The
full-TSO call remains insufficient for removal by itself.

This behavior is recorded as manifest schema 4 with survivor trim policy
`full_tso_endpoint_then_parse_exact_17mer`. Earlier v1.4 preprocessing
manifests are deliberately incompatible and force preprocessing to rerun.

## Read-length behavior

For complete R1/R2 structures, the 138-cycle budget, and recommended GEX versus
gRNA recipes, see the
[100-cycle Illumina sequencing note](ILLUMINA_100_CYCLE_SEQUENCING.md).

The removal requirement is approximately:

```text
TSO start offset + full TSO + shortest configured guide
```

With the current 30-bp TSO, 20/21-bp guides, and accepted start offsets 0–2,
this is 50–53 bp. Offset 0 is the intact GEX expectation; offsets 1–2 are only
damage/error tolerance. Therefore a 54-bp R1 exposes a complete guide for
removal. A typical
full dual-guide identity requires about `30 + 20 + 19 + 7 = 76 bp` before
allowing for start offset, so 54-bp R1 generally cannot expose a complete iBAR.
Step 1c emits a warning and the shared Step 3 logic rejects rather than invents
dual identity. `single_guide_mode: true` plus `ibar_length: 0` activates the
existing v1.2 short-read identity path and works without an iBAR.

## Configuration

The default is removal enabled, identity recovery disabled:

```yaml
gex_decontamination:
  enabled: true
  tso_max_edit_dist: 3
  max_tso_start: 2
  filter_matching:
    backend: packed_inline
    guide_distance: hamming
    max_guide_mismatches: 1
    allow_guide_indels: false
    full_guide_required: true
    remove_reference_supported_ambiguous_hits: true
  io:
    write_rejected_raw_pairs: true
    keep_rejected_raw_pairs_after_step1c: true
    write_read_decisions: false
  step1c:
    enabled: false
    output_namespace: gex_grna_identity
    strict_matching:
      enabled: true
      max_guide_edit_distance: 1
      max_ibar_edit_distance: 1
      require_resolved_ibar_in_dual_mode: true
      reject_ibar_decoupled: true
  combine_with_grna_arm:
    enabled: false
    combine_before_umi_collapse: true
    retain_source_provenance: true
```

When decontamination is enabled, Step 1a always writes the compact
valid-positive evidence needed by Step 1c:
the full-run corrected cell barcode, pB, UMI, R1 sequence, and R1 quality. It
does so during the existing preprocessing pass, so Step 1c never rescans the
GEX FASTQs and never estimates barcode corrections from the small rejected
subset. In the completed C2 run this reduced Step 1c input from 363,646,678
FASTQ pairs to 9,511 evidence rows.

Because this table is normally small, the production Step 1c array requests
one core by default. If its UGER request is explicitly raised for an unusually
large evidence table, the runner starts at most one worker per 100,000-row
Step-3 chunk and records both allocated slots and effective workers.

`write_rejected_raw_pairs: true` additionally retains synchronized full-length
R1/R2 positives for provenance, controlled QC, and safe historical-BAM reruns.
Those paired FASTQs are not the Step 1c matching input and can be disabled when
those audit artifacts are not required.

For release validation, QC independently decodes only these stored positive
pairs after fitting the normal full-GEX calibration interval and requires an
ordered, field-for-field match to Step 1a evidence. Production Step 1c itself
still reads only the evidence TSV.

On C2, the production Step 1c task consumed 9,511 evidence rows and completed
in about 5.4 seconds on one UGER core. Ordinary Step 3 produced 9,060 matches;
the conservative dual-guide policy retained 9,014 and excluded 46
lower-confidence or incompatible rows. This is intentionally a high-specificity
augmentation path, not an attempt to recover every structural contaminant.

In the isolated release merge, those 9,014 reads collapsed to 6,788 GEX UMI
clusters. Exact construct-aware overlap accounted for 6,041; after joint
directional collapse only 693 clusters were incremental to the dedicated
matrix. Dedicated-only recomputation remained exactly 3,272,782 UMIs and the
combined matrix contained 3,273,475. This high overlap and small net change are
why source merging remains disabled by default despite being supported.

Set both `step1c.enabled` and `combine_with_grna_arm.enabled` to include
GEX-derived evidence in final calls. Ordinary Step 3 output is retained as an
unfiltered audit product; only the configurable strict subset enters Step 4.
The defaults accept guide and iBAR edit distance at most one and require a
resolved, guide-compatible corrected iBAR in dual-guide mode; decoupled calls
are excluded. This deliberately favors specificity over recovery because this
optional arm can only add guide UMIs. Dedicated and
GEX rows then enter one per-library UMI accumulator before
`min_reads_per_umi` and directional collapse. An exact UMI observed in both
sources counts once.

## Existing split-pipe runs

The safe route never overwrites the historical run or reruns STAR:

```yaml
gex_decontamination:
  existing_alignment:
    enabled: true
    source_output_dir_template: "{output_dir}/gex_{library}"
    sandbox_output_dir_template: "{output_dir}/gex_{library}_v1_4"
    bam_backend: auto
```

For each library Step 1a then:

1. runs split-pipe `pre --one_step` on raw FASTQs in the new sandbox;
2. obtains positive IDs from the completed Step 1a manifest;
3. removes every alignment record for those QNAMEs from the historical
   `barcode_headAligned_anno.bam` or `barcode_headAligned_sorted.bam`;
4. atomically installs `process/barcode_headAligned_sorted.bam` plus
   `csm_harmonized_bam.complete.json` in the sandbox;
5. resumes the sandbox from split-pipe `--mode mol`.

The environment-gated mol cleanup seam preserves only the harmonized sorted BAM
referenced by its completion manifest; ordinary chunk/temp cleanup still runs.

The legacy-BAM route can remove filter-positive QNAMEs but cannot retroactively
change the sequence that was aligned. Consequently the new fuzzy full-TSO
survivor trimming applies to fresh STAR runs. Obtaining that trimming benefit
for a historical BAM requires a new alignment; the sandboxed BAM-only route
otherwise remains intentionally alignment-free.

The source BAM must retain split-pipe's `__OH_` original-header token and be in
query-name order. Missing tokens, incompatible ordering, unavailable positive
IDs, in-place paths, partial outputs, and record-count invariant failures are
hard errors. Step 1b automatically combines the configured sandbox directories.

## Outputs to inspect

Per GEX sublibrary:

```text
output/gex_grna_filter/<library>/
  filter_summary.json
  rejected_pairs.complete.json
  valid_positive_barcode_evidence.tsv.gz
  rejected_R1.fastq.gz             # when retained
  rejected_R2.fastq.gz
  filter_decisions.jsonl.gz        # only when explicitly enabled
```

The corresponding fresh-run or legacy sandbox also retains:

```text
<gex run>/process/barcode_headAligned_sorted.bam
<gex run>/process/csm_harmonized_bam.complete.json
```

Only the harmonized root BAM is exempted from normal mol cleanup; worker/chunk
temporaries are still removed.

The schema-4 manifest reports reads examined, TSO candidates, exact/fuzzy TSO evidence,
filter positives, unique/ambiguous sequence support, valid-barcode positives
removed from STAR, invalid-barcode positives, Step 1c evidence rows, synchronized
raw-pair counts and hashes, and an ordered paired-read-ID digest. It
also separates full-TSO survivor trims, Parse exact-17-mer survivor trims, and
untrimmed survivors, and records reference/removal-policy checksums so
split-pipe cannot reuse a stale preprocessed FASTQ.

With Step 1c:

```text
output/gex_grna_identity/grna_matches/
output/gex_grna_identity/grna_assignment/
  gex_grna_assignments.csv.gz
  gex_grna_umi_counts.csv.gz
  gex_grna_umi_matrix.mtx
  step1c_summary.yaml
output/gex_grna_identity/gex_grna_decontamination_summary.yaml
```

Production v1.4 points Step 3 directly at
`valid_positive_barcode_evidence.tsv.gz`; it does not create or recalibrate a
second barcode table.

With source merging:

```text
output/grna_assignment/
  grna_umi_counts_by_source.csv.gz
  grna_umi_matrix_combined.mtx
  grna_cell_assignment_combined.csv
  grna_gex_merge_summary.yaml
```

The combined long table and summary retain dedicated-only, GEX-only, overlap,
and incremental-GEX UMI-cluster statistics. Step 5 records which source drove
final calls and keeps dedicated, GEX, and combined sparse matrices available
for audit.

## Validation

Local regression tests:

```bash
python3 -m unittest discover -s tests -v
bash -n submit/submit_step1a_gex_array.sh
bash -n submit/submit_step1c_gex_grna_array.sh
bash -n submit/submit_step1c_gex_grna_aggregate.sh
bash -n submit/submit_pipeline.sh
```

`tests/test_splitpipe_gex_wrapper.py` also checks that the in-memory diff
compiles, refuses modified Parse sources, and leaves the vendor file hashes
unchanged.

Every multi-slot UGER wrapper requests a linear binding equal to its SMP slot
count. The worker count comes from `NSLOTS`; a submission must not request
`-pe smp N` without the matching `-binding linear:N`.

Step 1a resources are configurable in the main YAML:

```yaml
uger:
  step1a:
    slots: 16
    memory_per_slot_gb: 8
    runtime_hours: 48
```

The master and historical-BAM launchers validate these integer values and emit
matching `-pe`, `-binding`, `h_vmem`, and `h_rt` arguments. Use `--dry-run` to
inspect the exact command. More than 16 slots is not recommended for this
workflow: full-C2 preprocessing used eight Parse workers and averaged about
5.2 cores, while the observed memory footprint already fit comfortably.

Before production, run the pilot-derived QC module and enable
`require_qc_policy` only after its filter policy has been reviewed and marked
approved.

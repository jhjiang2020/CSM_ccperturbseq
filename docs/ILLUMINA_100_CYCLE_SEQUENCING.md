# Illumina sequencing-cycle recommendations

This note compares the 300-cycle baseline with two shorter-run strategies for
CSM ccPerturb-seq. The central constraint is simple: GEX can be shortened while
retaining v1.4 contaminant removal, but a 100-cycle dedicated-gRNA read cannot
reach the iBAR while also preserving the complete Parse R2 barcode amplicon.

## Choose the scenario

| Scenario | GEX recipe | Dedicated-gRNA recipe | Identity recovered | Step 1c merge recommendation |
|---|---|---|---|---|
| Baseline | 300-cycle, pilot PE150 | 300-cycle, pilot PE150 | Guide + stem/position + iBAR | Optional; supported when GEX R1 also reaches iBAR |
| 1. Short GEX, long gRNA | 100-cycle, preferably R1 60/R2 58 | 300-cycle PE150, or R1 >=100/R2 >=58 | Full guide/iBAR identity from the dedicated arm | Do not merge short-GEX Step 1c evidence into the iBAR-aware matrix |
| 2. Both arms short | 100-cycle R1 60/R2 58 | 100-cycle R1 60/R2 58 | Guide-only identity; no iBAR, position, or pair recovery | Optional and compatible at the `(cell, guide, UMI)` level |

Scenario 1 preserves the original dual-guide/iBAR experiment and is the
recommended cost-saving design when construct or pair identity matters.
Scenario 2 is appropriate only when guide-level assignment is scientifically
sufficient. It is not an iBAR-aware dual-guide analysis.

Read recipes apply to an entire sequencing run. GEX and gRNA libraries can use
different recipes only when they are sequenced separately.

## Why 100-cycle sequencing cannot recover iBAR

For the 100-cycle service considered here, 138 total cycles are available, with
10 cycles assigned to i7 and 10 to i5:

```text
R1 + R2 <= 138 - 10(i7) - 10(i5) = 118 cycles
```

The index reads demultiplex sequencing libraries; they are not Parse cell
barcodes and do not appear in demultiplexed R1. Chemistry-v3 cell and molecule
identity occupies the first 58 bases of R2:

```text
cycles       1-10    11-18       19-30       31-38       39-50       51-58
R2 field     10N UMI 8-nt BC3    12-nt linker 8-nt BC2   12-nt linker 8-nt BC1
```

R2 should therefore remain at least 58 nt. That leaves at most 60 nt for R1.

The relevant R1 structures are:

```text
GEX cDNA:
  [cDNA insert ...]

GEX gRNA carry-over:
  [Split_TSO, 30 nt] [guide, 20/21 nt] [stem ...]

Dedicated gRNA after circularization capture:
  [N7] [Split_TSO, 30 nt] [guide, 20/21 nt]
  [stem, 19 nt] [iBAR, 7 nt] [partial tracr ...]
```

Universal reverse transcription produces long poly(A)-primed cDNA and short
gRNA-probe products together. Size selection separates the GEX and gRNA
fractions, but a small amount of short gRNA product can remain in GEX. Those GEX
contaminants do not contain the leading N7. N7 is introduced later, when the
dedicated gRNA fraction undergoes circularization-capture amplification.

The resulting R1 boundaries are:

| Use | Required R1 bases | Interpretation |
|---|---:|---|
| split-pipe GEX input | 45 | Shorter GEX R1 is rejected by split-pipe. |
| v1.4 GEX contaminant removal | 50–53 | Full TSO + full 20/21-nt guide at accepted TSO offsets 0–2. |
| Dedicated guide without iBAR | 57–58 | N7 + TSO + full guide. R1 60 leaves 2–3 stem bases. |
| GEX-derived dual identity through iBAR | Up to 79 | TSO offset + TSO + guide + stem + iBAR. |
| Dedicated dual identity through iBAR | About 84 | N7 + TSO + guide + stem + iBAR. |
| Dedicated identity plus full configured tracr anchor | About 100 | Recommended for robust dual-guide/iBAR parsing. |

No allocation within R1+R2 <=118 can provide both the complete 58-nt Parse R2
and a complete dedicated-library iBAR in R1. Changing mismatch thresholds or
setting `ibar_length: 0` changes the identity model; it cannot recover bases
that were not sequenced.

## Baseline: GEX 300-cycle and gRNA 300-cycle

The pilot baseline used long paired reads for both arms. Long GEX R1 maximizes
cDNA sequence available to STAR and can expose guide, stem, and iBAR on the rare
GEX contaminants. Long dedicated-gRNA R1 supports guide, stem/position, iBAR,
and downstream-anchor validation.

For the dual CSM design, retain:

```yaml
grna_structure:
  single_guide_mode: false
  grna_lengths: [20, 21]
  ibar_length: 7
  tso_max_edit_dist: 3
  fixed_seq_max_mismatches: 2

grna_matching:
  max_lev_dist: 3
  max_lev_dist_ibar: 2

guide_assignment:
  dual_guide_counting_mode: true
  umi_dedup_ibar: true

barcode_extraction:
  min_r1_length: 84
  min_r2_length: 58
```

### Step 1c in the baseline

Step 1c is optional. With sufficiently long GEX R1, it can recover a complete
guide/iBAR identity from the compact valid-positive evidence emitted by Step
1a. It does not rescan the full GEX FASTQs.

For the validated C2 PE150 data, strict Step 1c evidence was compatible with
the iBAR-aware dedicated arm. Joint source processing added only 693 collapsed
UMI clusters beyond the dedicated matrix, changed two existing cell
assignments, and added four assignments. This small effect supports leaving
merging disabled by default, but the route is scientifically compatible when
the complete construct is observed.

Start with `step1c.enabled: true` and
`combine_with_grna_arm.enabled: false`; inspect the GEX-only product before
enabling the merge.

## Scenario 1: GEX 100-cycle, gRNA 300-cycle, iBAR retained

### Read recipe and scientific output

Sequence GEX and the dedicated gRNA library separately:

- GEX: R1 60/R2 58/i7 10/i5 10 is preferred. R1 54/R2 64 is supported, but
  bases 59–64 of R2 do not add Parse barcode information.
- Dedicated gRNA: retain the 300-cycle PE150 baseline, or at minimum use R1
  around 100 and R2 >=58.

The short GEX read retains enough sequence for v1.4 to recognize and remove a
TSO-plus-reference-guide contaminant. The long dedicated-gRNA read retains the
full guide/iBAR identity needed for dual-guide pair assignment.

### Configuration

Keep the dual-guide/iBAR model and enable GEX removal:

```yaml
grna_structure:
  single_guide_mode: false
  fixed_seq_1: "AAGCAGTGGTATCAACGCAGAGTGAATGGG"
  fixed_seq_2: "GTTTGAGAGCTAAGCAGGA"
  fixed_seq_3: "GTTTCAGAGCTATGCTGGA"
  fixed_seq_4: "GACTGCTTAGCAAGTT"
  fixed_seq_5: "AACAGCATAGCAAGTT"
  grna_lengths: [20, 21]
  ibar_length: 7
  tso_max_edit_dist: 3
  fixed_seq_max_mismatches: 2

grna_matching:
  max_lev_dist: 3
  max_lev_dist_ibar: 2
  save_failed_reads: true

guide_assignment:
  dual_guide_counting_mode: true
  umi_dedup_ibar: true

barcode_extraction:
  min_r1_length: 84
  min_r2_length: 58

gex_decontamination:
  enabled: true
  integration: "splitpipe_pre_router"
  tso_max_edit_dist: 3
  max_tso_start: 2
  tso_matching:
    backend: "edlib"
  filter_matching:
    backend: "packed_inline"
    guide_distance: "hamming"
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
    output_namespace: "gex_grna_identity"
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
  existing_alignment:
    enabled: false
```

### Step 1c merge decision

Do not merge Step 1c evidence from 54–60-nt GEX R1 into this iBAR-aware
dedicated matrix. The GEX read exposes the guide needed for removal but not the
stem+iBAR needed to identify the construct or guide position. Enabling Step 1c
for diagnostics is harmless and still avoids a full FASTQ rescan, but ordinary
dual-guide Step 3 should reject unresolved iBAR evidence. Do not relax that
gate merely to increase guide counts.

The correct output for this scenario is therefore:

- decontaminated short-read GEX expression;
- guide+iBAR assignments from the long dedicated-gRNA library;
- no GEX-derived guide contribution to the final guide matrix.

## Scenario 2: GEX and gRNA 100-cycle, no iBAR

### Read recipe and scientific output

Use R1 60/R2 58/i7 10/i5 10 for both arms. A single pooled run can use this
shared recipe. The analysis reports guide identity only; it cannot recover
iBAR, construct, guide position, or dual-guide pairing.

This scenario is naturally appropriate for a true single-guide/no-iBAR
experiment. It can also be used as a deliberately guide-only analysis of a
physical dual/iBAR library, but only if loss of construct and pair identity is
acceptable. Guide sequences must be unique enough to stand on their own.

### Configuration

Use the short-read single-guide path consistently in Steps 1c, 3, and 4:

```yaml
grna_structure:
  single_guide_mode: true
  fixed_seq_1: "AAGCAGTGGTATCAACGCAGAGTGAATGGG"
  fixed_seq_2: "GTTTGAGAGCTAAGCAGGA"
  fixed_seq_3: ""
  fixed_seq_4: "GACTGCTTAGCAAGTT"
  fixed_seq_5: ""
  grna_lengths: [20, 21]
  ibar_length: 0
  tso_max_edit_dist: 3
  fixed_seq_max_mismatches: 2

grna_matching:
  max_lev_dist: 1
  max_lev_dist_ibar: 0
  save_failed_reads: true

guide_assignment:
  dual_guide_counting_mode: false
  umi_dedup_ibar: false
  max_umi_edits: 1
  min_reads_per_umi: 1
  min_umi_threshold: 5
  multiplet_strategy: "all"
  compute_saturation_curve: true

barcode_extraction:
  min_r1_length: 20
  min_r2_length: 58

gex_decontamination:
  enabled: true
  integration: "splitpipe_pre_router"
  tso_max_edit_dist: 3
  max_tso_start: 2
  tso_matching:
    backend: "edlib"
  filter_matching:
    backend: "packed_inline"
    guide_distance: "hamming"
    max_guide_mismatches: 1
    allow_guide_indels: false
    full_guide_required: true
    remove_reference_supported_ambiguous_hits: true
  io:
    write_rejected_raw_pairs: true
    keep_rejected_raw_pairs_after_step1c: true
    write_read_decisions: false
  step1c:
    enabled: true
    output_namespace: "gex_grna_identity"
    strict_matching:
      enabled: true
      max_guide_edit_distance: 1
      max_ibar_edit_distance: 0
      require_resolved_ibar_in_dual_mode: false
      reject_ibar_decoupled: false
  combine_with_grna_arm:
    enabled: false          # review GEX-only results before changing to true
    combine_before_umi_collapse: true
    retain_source_provenance: true
  existing_alignment:
    enabled: false
```

The guide reference must still contain an `ibar_seq` column, but it may be
empty. Generate it with Step 0 using `ibar_length: 0`. If multiple reference
records share a guide sequence, no-iBAR data cannot distinguish them: Step 0
reports and collapses duplicates, while an uncollapsed supplied reference can
produce ambiguous Step 3 calls.

`max_lev_dist: 1` is a conservative starting point, not a universal constant.
Validate it against the actual guide set. Large 20–50k-guide references should
not inherit a permissive distance without an overlap audit.

### Step 1c merge decision

Step 1c is scientifically compatible with this guide-only identity model. It
processes only Step 1a's compact valid-positive evidence and matches the guide
without requiring an iBAR. Step 4 automatically disables iBAR-aware grouping
when `ibar_length: 0`.

Recommended rollout:

1. Run with Step 1c enabled and `combine_with_grna_arm.enabled: false`.
2. Inspect the GEX-only guide/UMI matrix and strict-match summary.
3. If the added evidence is acceptable, set `combine_with_grna_arm.enabled:
   true` and rebuild Steps 4 and 5.

When merging is enabled, dedicated and GEX reads enter the same accumulator
before read-depth filtering and directional UMI collapse. An identical
`(cell, guide, UMI)` observed in both sources counts once; novel GEX-derived
UMIs can increase an existing guide call or move a cell across
`min_umi_threshold`. Dedicated-only, GEX-only, and combined products remain
available for comparison.

## Common UGER settings and submission

For large or multi-lane GEX inputs, retain the 48-hour Step 1a default:

```yaml
uger:
  step1a:
    slots: 16
    memory_per_slot_gb: 8
    runtime_hours: 48
```

The Broad UGER launcher requests matching SMP slots and linear binding, then
passes `NSLOTS` to split-pipe. Sixteen is an evidence-based default rather than
a hard ceiling. Parse 1.6.4 caps its preprocessing pool at eight workers, and
the full-C2 run averaged about 5.2 cores; larger allocations should be justified
by STAR or other phase accounting. Total reserved memory increases with slots
because memory is requested per slot.

Use an absolute config path and dry-run before submission:

```bash
PROJECT_DIR=/medpop/rgupta/Jiahao/path/to/CSM_ccperturbseq
cd "${PROJECT_DIR}"
CONFIG="$(realpath my_100_cycle_v1_4.yaml)"

source /broad/software/scripts/useuse
reuse -q UGER || true
command -v qsub

bash submit/submit_pipeline.sh --config "${CONFIG}" --dry-run
bash submit/submit_pipeline.sh --config "${CONFIG}"
```

These wrappers are specific to Broad UGER. Another scheduler requires an
adapted submission layer. Heavy and split-pipe commands must run through the
submitted jobs, not on the login node.

## Why decontamination matters more for short GEX R1

The controlled C2 experiment aligned full-length and synthetic 54-nt
contaminants to the same human genome with identical STAR settings. None of
10,559 full-length contaminants mapped. After 54-nt truncation and the same
Parse TSO preprocessing, 8,113/10,559 (76.83%) mapped and 4,008 (37.96%)
received a false gene assignment. Of those, 2,431 reads and 1,894 collapsed
molecules were assigned to the guide's intended target gene; the v1.4 filter
removed every target-concordant event in this panel.

After TSO trimming, short contaminant R1 consists mostly of guide sequence plus
a few stem bases, which can align spuriously. The longer stem+iBAR+tracr context
in PE150 contaminants suppressed mapping in C2. This supports enabling v1.4
GEX decontamination in both short-read scenarios even though the raw
contaminant frequency is low.

## Pre-run and post-run checks

Before submission:

1. Confirm whether the two arms are on the same run and therefore share one
   read recipe.
2. Preserve R2 >=58 for the complete Parse v3 barcode/UMI amplicon.
3. Decide explicitly whether the output must retain guide+iBAR/pair identity or
   guide identity alone.
4. For Scenario 2, confirm that guide sequences are unique enough without iBAR.
5. Test that each `gex_prefix` and `grna_prefix` resolves only the intended
   FASTQs.
6. Pass an absolute HPC config path and inspect the `--dry-run` output.

After processing, inspect:

- Step 1a `filter_summary.json`, manifest completeness, and contaminant rate;
- Step 3 extraction scenarios, ambiguous/no-match counts, and guide-distance
  distribution;
- Step 1c strict-match and GEX-only UMI summaries when enabled;
- dedicated-only versus combined UMI and cell-assignment deltas before adopting
  a Step 1c merge.

See the [v1.4 decontamination guide](GEX_DECONTAMINATION_V1_4.md) and
[contamination QC guide](GEX_GRNA_CONTAMINATION_QC.md) for output contracts and
the full C2 evidence.

# Special Note: 100-Cycle Illumina Sequencing

This note defines the supported sequencing recipes for the CSM ccPerturbSeq
workflow when the sequencing service permits **138 total cycles**, with
**10 cycles reserved for i7** and **10 cycles reserved for i5**.

> [!IMPORTANT]
> The usable paired-read budget is therefore
> `R1 + R2 <= 138 - 10(i7) - 10(i5) = 118 cycles`.
> The i7 and i5 reads demultiplex sequencing libraries; they are not the Parse
> cell barcodes. For both GEX and gRNA libraries, biological insert information
> is in R1 and the Parse UMI/cell-barcode amplicon is in R2.

## Read orientation confirmed by the protocol and pilot data

The protocol figures and the pilot 300-cycle FASTQs establish the following
orientation.

### R1: biological insert

**GEX library**

```text
R1 = cDNA insert

most fragmented molecules:
  cDNA sequence ......................................................>

short molecules retaining the 5' end:
  [partial or full Split_TSO] [cDNA] .................................>

gRNA carry-over constructed through the GEX library workflow:
  expected: [Split_TSO, 30 nt] [guide, 20/21 nt] [stem ...] ..........>
  tolerated: full-TSO alignment begins at offset 1-2 after rare 5' damage/error
```

Universal reverse transcription uses both poly(A)-directed and gRNA-directed
probes, so long cDNA and short gRNA-derived products are initially made
together. The GEX/gRNA size selection is not absolute: a low-frequency short
gRNA product can remain in the GEX fraction. Most GEX cDNA molecules were
fragmented upstream of sequencing and therefore do not show the TSO in R1;
short cDNA molecules can retain a partial or full TSO. Intact gRNA carry-over
starts with the full TSO at R1 offset 0 followed by a reference guide. The v1.4
filter accepts a detected full-TSO start at offset 1 or 2 only to tolerate
low-frequency 5' damage and/or sequencing error; those offsets are not an
expected library prefix. GEX-derived contaminants do **not** contain the
leading N7 random sequence added during the dedicated gRNA library's later
circularization-capture amplification.

After R2 barcode correction, split-pipe trims recognized TSO sequence from GEX
R1 and sends only R1 to human STAR. R2 is never aligned.

**Dedicated gRNA library**

```text
R1 = [N7 random prefix]
     [Split_TSO, 30 nt]
     [guide, 20/21 nt]
     [CSM stem 1/2, 19 nt]
     [iBAR, 7 nt]
     [CSM partial tracr 1/2, 16 nt ...]
```

After size selection, the dedicated gRNA fraction undergoes a separate
circularization-capture amplification. This chemistry places a random N7
between the Illumina sequencing-primer/index side of the molecule and the TSO.
Consequently, in an intact dedicated gRNA read the TSO begins exactly after N7,
at zero-based R1 offset 7. The 10-nt i7 and 10-nt i5 sequences consume
instrument cycles as separate index reads; they do not appear in demultiplexed
R1 and do not create an additional R1 offset. No separate N3 prefix is assumed.
Step 3 locates the expected N7-to-TSO boundary and parses guide, stem, iBAR, and
optional downstream tracr sequence relative to the TSO endpoint. The random
prefix is not part of the guide.

### R2: Parse UMI and cell identity

Both GEX and gRNA libraries use the same chemistry-v3 R2 layout:

```text
cycles       1-10    11-18       19-30       31-38       39-50       51-58
R2 field     10N UMI 8-nt BC3    12-nt linker 8-nt BC2   12-nt linker 8-nt BC1
```

The complete R2 amplicon is **58 nt**. Cell identity is the corrected
combination of BC1, BC2, and BC3; the 10-nt leading sequence is the molecular
UMI. Parse split-pipe 1.6.4 derives its R2 minimum from this 58-nt v3 amplicon,
and its independent hard minimum for raw GEX R1 is 45 nt.

## What fits within 118 paired-read cycles

The relevant R1 boundaries, including the structures above, are:

| Use case | R1 bases required | Consequence |
|---|---:|---|
| split-pipe GEX input gate | 45 | Shorter GEX R1 pairs are rejected. |
| v1.4 GEX contaminant removal | 50-53 | `0-2 accepted start + 30 TSO + 20/21 guide`; offset 0 is expected and offsets 1-2 are damage/error tolerance. |
| Dedicated single-guide/no-iBAR, full guide | 57-58 | `N7 prefix + 30 TSO + 20/21 guide`; 60-nt R1 leaves 2-3 stem bases. |
| GEX-derived dual identity through iBAR | up to 79 | `2 offset + 30 TSO + 21 guide + 19 stem + 7 iBAR`; impossible with R2 >=58. |
| Dedicated dual-guide identity through iBAR | about 84 | `N7 prefix + 30 TSO + 21 guide + 19 stem + 7 iBAR`; impossible with R2 >=58. |
| Dedicated dual identity plus full configured partial tracr | about 100 | Adds the 16-nt tracr anchor; use a longer-cycle kit. |

No allocation within `R1 + R2 <= 118` can provide both a complete 58-nt R2
and a complete dual-guide iBAR in R1. Lowering mismatch thresholds, minimum-read
settings, or iBAR length in the YAML cannot recover bases that were not
sequenced.

## Recommended sequencing recipes

Read recipes apply to an entire sequencing run. GEX and gRNA libraries pooled
on the same run cannot receive different R1/R2 lengths.

| Library and analysis goal | R1 | R2 | i7 | i5 | Recommendation |
|---|---:|---:|---:|---:|---|
| GEX, preferred 100-cycle configuration | 60 | 58 | 10 | 10 | Uses all 118 paired-read cycles, maximizes alignable cDNA, and gives a 7-nt margin above the longest v1.4 removal boundary. |
| GEX, planned 54/64 compatibility recipe | 54 | 64 | 10 | 10 | Supported and matches the synthetic 54-bp validation design. R2 bases 59-64 are not used by the v3 barcode decoder. |
| Dedicated single-guide, no iBAR | 60 | 58 | 10 | 10 | Supported by the v1.2 short-read path; normally captures the full guide and 2-3 stem bases after the N7 prefix. |
| Dedicated dual-guide with iBAR | — | — | 10 | 10 | **Not supported with this cycle budget.** Use a longer-cycle kit. R1 >=84 is the structural minimum; R1 around 100 is recommended to include the configured partial tracr anchor, with R2 >=58. |

For a single pooled 100-cycle run, use **R1=60 and R2=58**. This is suitable
for GEX and for a genuinely single-guide/no-iBAR gRNA arm. It does not make the
current dual-guide/iBAR arm identifiable. If the experiment uses the current
dual CSM design, sequence the dedicated gRNA library separately with a longer
recipe, for example **R1=100 and R2=58**, while the GEX library can use the
short recipe.

The 54/64 GEX recipe remains a valid conservative rollout target. If the core
allows the allocation to change, 60/58 is preferred because the extra six
cycles are more informative as cDNA R1 than beyond the completed 58-nt R2
barcode amplicon.

## Recommended pipeline configuration

### Common v3 barcode settings

Use the correct kit for the experiment and retain chemistry v3. The custom
dedicated-guide Step 2 barcode reader should explicitly enforce the complete
R2 length:

```yaml
splitpipe:
  kit: "WT"                 # change only if a different Parse kit was used
  chemistry: "v3"

barcode_extraction:
  bc_max_edit_dist: 1
  min_r1_length: 20         # preserves the v1.2 truncated-guide rescue path
  min_r2_length: 58         # complete v3 UMI + BC3 + BC2 + BC1 amplicon
  n_processes: 8
  max_reads: 0
```

`barcode_extraction` controls custom dedicated-gRNA Step 2. The production
Step 1c wrapper requests one core and starts additional workers only when an
explicitly larger allocation and multiple 100,000-row evidence chunks warrant
them. Step 1c does not decode GEX R2: it consumes the compact, full-run
corrected-barcode evidence emitted during Step 1a. These settings do not
override split-pipe's GEX gates: raw GEX R1 must be >=45 and chemistry-v3 R2
must be >=58.

### Short GEX with dual-guide contamination removal

For a fresh 54/64 or 60/58 GEX run using the current dual-guide reference:

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

gex_decontamination:
  enabled: true
  integration: "splitpipe_pre_router"
  tso_max_edit_dist: 3
  # Intact GEX contaminants start at 0; 1-2 tolerate rare 5' damage/error.
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
    enabled: false          # recommended for short dual/iBAR GEX R1
    output_namespace: "gex_grna_identity"
    strict_matching:
      enabled: true
      max_guide_edit_distance: 1
      max_ibar_edit_distance: 1
      require_resolved_ibar_in_dual_mode: true
      reject_ibar_decoupled: true
  combine_with_grna_arm:
    enabled: false          # do not merge identity that R1 cannot resolve
    combine_before_umi_collapse: true
    retain_source_provenance: true
  existing_alignment:
    enabled: false          # fresh FASTQ run; enable only for the sandboxed BAM route

uger:
  step1a:
    slots: 16
    memory_per_slot_gb: 8
    runtime_hours: 48
```

Step 1a can remove a short dual-guide contaminant from GEX because it requires
only full TSO plus full guide. Step 1c should remain disabled for 54-60-nt dual
GEX R1 because stem, iBAR, position, and pair identity are not fully observed.
If Step 1c is enabled for diagnostics, its short-iBAR warning and high identity
failure rate are expected; such evidence must not be merged into final calls.
Enabling Step 1c does not trigger a second full FASTQ scan: only Step 1a's
valid-positive evidence rows are examined.
The 48-hour Step 1a limit is recommended for large or multi-lane GEX inputs and
can be changed in YAML. More cores are not a substitute for matcher efficiency:
Parse preprocessing uses eight workers in the tested configuration.

### Short dedicated gRNA library: supported mode

The supported 60/58 dedicated-gRNA configuration is specifically
**single-guide mode without iBAR**:

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
  max_lev_dist: 1           # conservative starting point for a large short-read screen
  max_lev_dist_ibar: 0      # ignored when ibar_length is 0
  save_failed_reads: true   # recommended during rollout/QC

guide_assignment:
  dual_guide_counting_mode: false
  umi_dedup_ibar: false
  max_umi_edits: 1
  min_reads_per_umi: 1
  compute_saturation_curve: true
```

`max_lev_dist: 1` is a conservative starting point, not a universal biological
constant. Validate it against the actual guide reference and sequencing-error
profile. Small pilot references may support a larger validated value; large
20-50k-guide references should not inherit a permissive distance without an
overlap and performance audit.

### Dual-guide/iBAR gRNA library: longer-read configuration

For the current dual CSM library, keep the dual settings rather than pretending
the library has no iBAR:

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
  save_failed_reads: true

guide_assignment:
  dual_guide_counting_mode: true
  umi_dedup_ibar: true

barcode_extraction:
  min_r1_length: 84         # full dedicated-library guide + stem + iBAR
  min_r2_length: 58
```

Use R1 around 100 nt when possible so the ordinary Step 3 combined stem/tracr
checks see the full configured 16-nt partial tracr. Do not set `ibar_length: 0`
or enable single-guide mode merely to force 60-nt dual-library reads through
the pipeline; that would change the experiment's identity model and discard the
construct/pair information the dual workflow is designed to recover.

The C2 controlled STAR experiment directly tested the GEX consequence of a
54-nt R1. All 10,559 full-length contaminants failed to map, but 8,113/10,559
(76.83%) mapped after 54-nt truncation and the same Parse TSO preprocessing;
4,008 (37.96%) received a false gene assignment. TSO removal leaves mostly
23-27 nt of spacer plus a few stem bases, which explains why the short arm can
align while the full stem-iBAR-tracr context suppresses mapping. This result
supports enabling the conservative v1.4 GEX decontamination filter for the same
universal-RT/size-selection chemistry when using a short-read recipe.
Of those false assignments, 2,431 reads and 1,894 collapsed molecules were
assigned to the guide's intended target gene. The v1.4 decision removed every
target-concordant event observed in the C2 panel.

## Pre-run and post-run checks

Before submission:

1. Confirm the run recipe and actual FASTQ lengths, including `R2 >=58`.
2. Confirm whether GEX and gRNA libraries share one sequencing run and therefore
   one read recipe.
3. Confirm whether the guide design is genuinely single/no-iBAR or dual/iBAR.
4. Test that `gex_prefix` and `grna_prefix` select only the intended FASTQs.
5. Pass the configuration to UGER as an absolute HPC path, as described in the
   [quick start](../README.md#3-initialize-uger-and-dry-run).

After processing, inspect the split-pipe barcode-valid fraction, Step 1a
`filter_summary.json`, the full-TSO and guide-positive rates, Step 3 short-read
scenario/failure counts, and the expected iBAR-unobservable warning. See the
[v1.4 decontamination guide](GEX_DECONTAMINATION_V1_4.md) and
[consolidated contamination QC guide](GEX_GRNA_CONTAMINATION_QC.md) for the
corresponding outputs, acceptance checks, and field-by-field result
interpretation.

# Detailed Step-by-Step Workflow

This document describes each processing stage, its inputs and outputs, matching
and assignment decisions, and the UGER dependency graph. For setup and
submission commands, begin with the [quick start](../pipeline_overview.md).
Versioned behavior changes are tracked in the [release notes](../RELEASE_NOTES.md),
and operational details for GEX decontamination are in the
[v1.4 guide](GEX_DECONTAMINATION_V1_4.md). Read-recipe constraints and matching
configurations are summarized in the
[100-cycle Illumina sequencing note](ILLUMINA_100_CYCLE_SEQUENCING.md).

## UGER Job Dependency Chain

```mermaid
flowchart TD
    S0["Step 0 (manual)"] --> S3_arr

    subgraph GEX["GEX Branch"]
        S1a["Step 1a: split-pipe<br>(array ×N)"]
        S1b["Step 1b: combine"]
        S1c["Step 1c: GEX guide identity<br>(optional array + aggregate)"]
        S1a --> S1b
        S1a -.-> S1c
    end

    subgraph gRNA["gRNA Branch"]
        S2_arr["Step 2: Extract BCs<br>(array ×N)"]
        S3_arr["Step 3: Match gRNA+iBAR<br>(array ×N)"]
        S4["Step 4: Assign/merge to cells<br>(1 job, 4 CPUs)"]
        S2_arr --> S3_arr --> S4
    end

    S1b --> S5
    S1c -. optional source union .-> S4
    S4 --> S5["Step 5: Build AnnData"]

    style GEX fill:#fff3e0,stroke:#e65100
    style gRNA fill:#e8f5e9,stroke:#2e7d32
```

> [!TIP]
> GEX and gRNA branches run in **parallel** — they only merge at Step 5.
> You can test guide extraction and assignment with different parameters by re-runing Step 3, 4 and 5 with modified `config.yaml`.

---

## Step 0: Format Reference

**Script**: [00_format_reference.py](../scripts/00_format_reference.py)

**Run**: Manual (one-time)

Converts the dual-guide library design file into the gRNA/iBAR reference lookup used by downstream steps.

```mermaid
flowchart TD
    IN[("Library Design CSV<br>~700 dual-guide rows")]
    IN --> SPLIT["Split each row into<br>2 entries (pos1 + pos2)"]
    SPLIT --> RC{"rev_comp_ibar?"}
    RC -->|Yes| REVCOMP["Reverse-complement<br>iBAR sequences"]
    RC -->|No| KEEP["Keep iBAR as-is"]
    REVCOMP --> BUILD["Build reference table"]
    KEEP --> BUILD
    BUILD --> QC["QC: Check iBAR<br>uniqueness per pair_id"]
    QC --> OUT[("grna_ibar_reference.csv<br>~1400 rows × 14 cols")]

    style IN fill:#e3f2fd
    style OUT fill:#c8e6c9
```

| Input column | Output (Position 1) | Output (Position 2) |
|---|---|---|
| `spacer_1` | `grna_seq` | — |
| `spacer_2` | — | `grna_seq` |
| `iBAR_1` | `ibar_seq` (rev-comp) | — |
| `iBAR_2` | — | `ibar_seq` (rev-comp) |
| `target` | `target_combo` | `target_combo` |
| — | `stem_variant = A` | `stem_variant = B` |

#### Reference overlap check

When run with `--config`, Step 0 automatically checks whether the configured `max_lev_dist` and `max_lev_dist_ibar` could cause **ambiguous matching** in Step 3. The check computes all pairwise Levenshtein distances between gRNA sequences and between iBARs within each gRNA, classifying pairs into three zones:

| Zone | Condition | Impact |
|------|-----------|--------|
| **CRITICAL** | d(A,B) ≤ max_lev_dist | A read from ref A with ≤max_lev_dist errors matches ref B equally → **guaranteed ambiguity**. Read is **discarded** (not mis-assigned). |
| **WARNING** | max_lev_dist < d(A,B) ≤ 2×max_lev_dist | Error balls partially overlap — some error-prone reads may fall in the overlap region and be lost. |
| **SAFE** | d(A,B) > 2×max_lev_dist | No overlap possible. |

The check reports the minimum pairwise distance and a **safe limit** for each parameter (the largest distance that guarantees no overlap: `⌊(min_dist - 1) / 2⌋`).

To run the overlap check on an existing reference without reformatting:

```bash
# Using config distances
python scripts/00_format_reference.py --check-only references/grna_ibar_reference.csv --config config.yaml

# With custom distances
python scripts/00_format_reference.py --check-only references/grna_ibar_reference.csv --max-lev-dist 2 --max-lev-dist-ibar 2
```

> [!NOTE]
> Step 3's `NeighborhoodMatcher` always returns the **unique best match**. When multiple references tie at the same best distance, the read is marked **ambiguous** and discarded — there is no risk of wrong assignment, only data loss. Users can safely proceed with overlapping references if they accept some ambiguous-read loss (visible in Step 3 QC as "% ambiguous").


---

## Step 1: GEX Processing (split-pipe)

**Script**: [01_run_gex_splitpipe.sh](../scripts/01_run_gex_splitpipe.sh)

**UGER**: Step 1a (array ×N) → Step 1b (combine), plus optional Step 1c

The assay's universal reverse-transcription reaction contains both poly(A)
probes and gRNA probes. Long gene-transcript cDNA and short gRNA products are
therefore created together and separated later by size. Step 1a decontamination
targets the small amount of short gRNA product that escapes into the GEX size
fraction. Pilot C1/C2 raw scans measured approximately 32–37 structural
contaminants per million reads; this is low in absolute terms but can materially
bias a target gene expressed at a similar or lower rate if STAR assigns the
guide-derived sequence to that gene.

When enabled, Step 1a installs a checksum-validated runtime overlay over the
pristine split-pipe 1.6.4 preprocessing and molinfo modules; it does not edit
the Parse installation. TSO-bearing reads are screened with a packed guide
index. A complete guide match at Hamming distance 0 or 1 is removed before the
human STAR input is written. For filter-negative reads, a detected full/fuzzy
TSO is trimmed at its full endpoint; without a full-TSO call, Parse's original
exact terminal 17-mer rule remains the fallback. The full-TSO call alone never
causes removal.

Step 1a is deliberately removal-only, but it emits a compact table for every
valid-barcode filter positive using split-pipe's already-corrected cell
barcode, pB, UMI, R1 sequence, and quality. Optional Step 1c runs ordinary
Step 3 only on that table, applies stricter guide/iBAR acceptance thresholds,
and can pass the conservative subset to Step 4. It never rescans the full GEX
FASTQs or recalibrates barcodes from rejected reads. The array size N follows
`config.yaml`.

```mermaid
flowchart TD
    subgraph Step1a["Step 1a: Per-sublibrary (array ×N)"]
        FQ[("GEX FASTQ<br>R1 + R2")] --> PRE["split-pipe pre<br>R2 correction + TSO route"]
        PRE --> FILT{"full TSO + packed guide hit<br>Hamming 0/1?"}
        FILT -->|No| TRIM{"full/fuzzy TSO?"}
        FILT -->|Yes| EV[("corrected barcode evidence<br>+ optional raw R1/R2")]
        TRIM -->|Yes| FULLTRIM["trim at full-TSO endpoint"]
        TRIM -->|No| PARSETRIM["Parse exact-17-mer fallback"]
        FULLTRIM --> STAR["human STAR"]
        PARSETRIM --> STAR
        STAR --> SP["split-pipe<br>mol → downstream"]
        SP --> GEX[("gex_sublib_1/<br>gex_sublib_2/...")]
    end

    EV -. optional Step 1c .-> ID["ordinary Step 3 → strict subset → Step 4<br>GEX-only guide matrix"]

    subgraph Step1b["Step 1b: Combine"]
        GEX --> COMB["split-pipe<br>--mode comb"]
        COMB --> OUT[("gex_combined/<br>adata.h5ad")]
    end

    style FQ fill:#e3f2fd
    style OUT fill:#c8e6c9
```

At 54-bp R1, a 30-bp TSO plus a complete 20/21-bp guide is observable, so
removal remains supported. An intact GEX-derived contaminant starts at offset
0; the configured offsets 1–2 only tolerate rare 5' damage/sequencing error.
A standard dual-guide stem plus 7-bp iBAR is not observable at that length;
Step 1c warns and does not invent identity. Single-guide/no-iBAR mode uses the
existing v1.2 short-read logic.

For a historical compatible split-pipe BAM, Step 1a can instead preprocess raw
reads in a new sandbox, remove all records for filter-positive QNAMEs from the
name-sorted/annotated source BAM, and resume at `--mode mol`. This route avoids
STAR but cannot retroactively apply the new full/fuzzy-TSO survivor trimming;
that benefit requires a fresh alignment. See the
[v1.4 guide](GEX_DECONTAMINATION_V1_4.md#existing-split-pipe-runs).

---

## Step 2: Extract gRNA Barcodes

**Script**: [02_extract_grna_barcodes.py](../scripts/02_extract_grna_barcodes.py)

**UGER**: Array ×N, 32GB RAM, 8 processes

Extracts cell barcodes from gRNA library R2 using Parse Biosciences barcode correction logic, producing a TSV linking each read to its cell identity.

```mermaid
flowchart TD
    FQ1[("gRNA R1 FASTQ<br>(gRNA sequence)")] --> PAIR["Read paired<br>FASTQ entries"]
    FQ2[("gRNA R2 FASTQ<br>(cell barcodes)")] --> PAIR

    PAIR --> LENCK{"R1 ≥ 20bp?<br>R2 ≥ 50bp?"}
    LENCK -->|No| SKIP["Discard"]
    LENCK -->|Yes| EXTRACT["Extract 3 barcodes<br>+ UMI from R2"]

    EXTRACT --> AMBIG{"Ambiguous<br>N bases?"}
    AMBIG -->|Yes| SKIP
    AMBIG -->|No| CORRECT["Correct barcodes<br>(Levenshtein ≤ 1)"]

    CORRECT --> VALID{"All 3 BCs<br>corrected?"}
    VALID -->|No| SKIP
    VALID -->|Yes| BUILD["Build bc_wells<br>= w1_w2_w3"]

    BUILD --> OUT[("grna_barcodes/<br>sublib_1_grna_barcodes.tsv.gz<br>Cols: read_id, bc_wells,<br>bc_bci, umi, r1_seq, r1_qual")]

    style FQ1 fill:#e3f2fd
    style FQ2 fill:#e3f2fd
    style OUT fill:#c8e6c9
```

---

## Step 3: Match gRNA + iBAR

**Script**: [03_match_grna_ibar.py](../scripts/03_match_grna_ibar.py)

**UGER**: Array ×N, 64GB RAM, 8 processes (multiprocessing)

The most complex step — locates the gRNA and iBAR within each R1 read, then fuzzy-matches against the reference.

```mermaid
flowchart TD
    BC[("Barcode TSV<br>(from Step 2)")] --> CHUNK["Read in chunks<br>(100K rows)"]
    REF[("gRNA/iBAR<br>Reference")] --> MATCHERS["Build<br>NeighborhoodMatchers"]

    subgraph WorkerPool["Multiprocessing Pool (8 workers)"]
        CHUNK --> ROWS["Per-read processing"]
        MATCHERS --> ROWS

        ROWS --> TSO["1. Find Split_TSO<br>(edlib semi-global)"]
        TSO --> EXTRACT["2. Extract gRNA + iBAR<br>by position relative to TSO"]
        EXTRACT --> STEM{"3. Match stem<br>(variant A or B)"}
        STEM --> GRNA["4. Fuzzy-match gRNA<br>(NeighborhoodMatcher<br>Levenshtein ≤ 3)"]
        GRNA --> IBAR_G["5a. Global iBAR match"]
        GRNA --> IBAR_W["5b. Within-guide<br>iBAR match"]
        IBAR_G --> RESOLVE["6. Resolve corrected iBAR<br>+ decoupling flag"]
        IBAR_W --> RESOLVE
        ROWS -->|failed| FAIL["Record failure reason"]
    end

    RESOLVE --> OUT[("grna_matches/<br>sublib_1_grna_matches.tsv.gz<br>~17 columns per read")]
    FAIL -->|save_failed_reads| FOUT[("grna_matches/<br>sublib_1_grna_failed.tsv.gz<br>read + failure_reason")]

    style BC fill:#e3f2fd
    style REF fill:#e3f2fd
    style OUT fill:#c8e6c9
    style FOUT fill:#ffcdd2
```

### R1 Read Structure

The leading N7 is introduced only when the size-selected gRNA fraction is
further amplified by circularization capture. It lies between the Illumina
sequencing-primer/index side of the library molecule and the TSO, placing the
TSO at R1 offset 7. GEX carry-over escaped size selection before this separate
amplification and therefore begins at the TSO rather than N7.

```
┌──────────────────────────────────────────────────────────┐
│ N7 random prefix │ Split_TSO (30bp) │ gRNA (20-21bp) │ stem (19bp) │ iBAR (7bp) │ tracr │
│ var  │ AAGCAGT...GAATGGG │   <variable>  │ GTTTGAGAG.. │  <7bp>     │ GACT  │
│ off  │                   │               │     or      │            │  or   │
│      │                   │               │ GTTTCAGAG.. │            │ AACA  │
└──────────────────────────────────────────────────────────┘
```

### Dual iBAR Matching Logic
```mermaid
flowchart TD
    INPUT["Query iBAR"] --> GLOBAL["Global match<br>(all iBARs)"]
    INPUT --> WITHIN["Within-guide match<br>(iBARs for matched gRNA only)"]

    GLOBAL --> DEC{Both found?}
    WITHIN --> DEC

    DEC -->|Within-guide hit| USE_W["Use within-guide<br>as corrected iBAR"]
    DEC -->|Global only| USE_G["Use global<br>(flag decoupled)"]
    DEC -->|Both agree| USE_W
    DEC -->|Disagree| FLAG["Use within-guide<br>flag as decoupled"]
    DEC -->|Neither| EMPTY["ibar_seq_corrected = ''"]
```

> [!IMPORTANT]
> **iBAR correction requires unique-best matching by design.** The within-guide iBAR match (`_match_ibar_within_guide`) finds the single closest reference iBAR (by Levenshtein distance) among the iBARs associated with the resolved gRNA. If two or more reference iBARs tie at the same best distance, the match is rejected — no correction is applied. This conservative approach is intentional: since each iBAR is unique to a specific construct (pair_id), an ambiguous correction could silently assign reads to the wrong construct, introducing noise into downstream analysis. It is better to leave a read uncorrected (and flag it as decoupled) than to guess.

> [!NOTE]
> **Shared gRNA disambiguation also requires unique-best iBAR matching.** When a gRNA sequence is shared across positions (e.g., control guides), disambiguation tries iBAR matching per candidate. If exactly one candidate's iBAR set contains a unique-best match, that candidate is chosen. If multiple candidates match (both have an iBAR within the edit distance threshold), disambiguation falls back to stem variant identity. This prevents false assignment when the short iBAR length (7bp) and permissive edit distance threshold create overlapping Hamming balls across candidates.

---

## Step 4: Assign gRNAs to Cells

**Script**: [04_assign_grna_to_cells.py](../scripts/04_assign_grna_to_cells.py)

**UGER**: 1 job, 32GB RAM, 4 CPUs (parallelized per library)

Collapses read-level matches to cell-level gRNA assignments via **streaming UMI deduplication**.

```mermaid
flowchart TD
    subgraph Stream["Streaming UMI Dedup (per library, ×N parallel)"]
        MATCH[(".tsv.gz<br>match file")] --> CHUNKS["Read 500K-row<br>chunks"]
        CHUNKS --> ACCUM["Accumulate<br>{(cell, guide, ibar, target):<br>{umi: read_count}}"]
        ACCUM --> FILTER{"min_reads_per_umi<br>> 1?"}
        FILTER -->|Yes| PRUNE["Remove UMIs with<br>reads < threshold"]
        FILTER -->|No| COLLAPSE
        PRUNE --> COLLAPSE["Directional UMI<br>collapsing (Hamming ≤ 1)"]
        COLLAPSE --> DEDUP[("Per-library<br>deduped DataFrame")]
    end

    C1[("sublib_1 matches<br>273MB gz")] --> Stream
    C2[("sublib_2 matches<br>636MB gz")] --> Stream
    TN[("... matches<br>...MB gz")] --> Stream

    DEDUP --> CONCAT["Concatenate<br>deduped results"]
    CONCAT --> ASSIGN["Assign guides to cells<br>(UMI threshold ≥ 5)"]
    ASSIGN --> MATRIX["Build sparse<br>cells × guides matrix"]

    ASSIGN --> OUT1[("grna_cell_assignment.csv")]
    CONCAT --> OUT2[("grna_umi_counts.csv")]
    MATRIX --> OUT3[("grna_umi_matrix.mtx")]
    CONCAT --> OUT4[("umi_threshold_sensitivity.png")]
    ACCUM --> OUT5[("reads_per_umi_distribution.png")]
    ACCUM --> OUT6[("saturation_curve.png")]
    CONCAT --> OUT7[("median_reads_per_umi_per_cell.png")]

    style C1 fill:#e3f2fd
    style C2 fill:#e3f2fd
    style OUT1 fill:#c8e6c9
    style OUT2 fill:#c8e6c9
    style OUT3 fill:#c8e6c9
    style OUT4 fill:#c8e6c9
    style OUT5 fill:#fff9c4
    style OUT6 fill:#fff9c4
    style OUT7 fill:#fff9c4
```

### UMI Directional Collapsing Algorithm

Uses the UMI-tools directional method (Smith et al., 2017). A child UMI is only collapsed into a parent if its read count is consistent with being a PCR error offspring — i.e., the parent's count is at least 2× the child's count minus 1. This prevents over-collapsing of genuinely distinct molecules that happen to be close in sequence space.

```mermaid
flowchart TD
    UMIS["Sort UMIs by<br>read count (desc)"] --> ITER["For each parent UMI<br>(high count)"]
    ITER --> CHECK["For each child UMI<br>(lower count)"]
    CHECK --> HAM{"Hamming<br>dist ≤ 1?"}
    HAM -->|No| NEXT["Skip"]
    HAM -->|Yes| DIR{"count_high ≥<br>2×count_low − 1?"}
    DIR -->|Yes| MERGE["Collapse child<br>into parent"]
    DIR -->|No| NEXT["Keep as separate<br>molecule"]
    MERGE --> COUNT["Count remaining<br>unmerged clusters<br>= UMI count"]
```

### Cell Assignment Logic (Dual Guide Mode)

In dual-guide mode (`dual_guide_mode: true`), the pipeline determines pair completion rather than independent guides.

1.  **Pair Grouping**: Reads are assigned to a reference `pair_id` via the detected `pos1` or `pos2` guide and their associated `ibar_seq`.
2.  **iBAR Inference**: If a cell has reads for `pos1` but not `pos2`, the pipeline uses the detected `pos1_ibar` to infer the missing expected `pos2` guide based on the dual-guide library design. (Or vice versa for observed `pos2`).
3.  **Thresholding** (two-tier):
    -   **Complete pairs** (both positions detected): `pos1_umi + pos2_umi >= min_umi_threshold` (pair-level sum). Example: pos1=3, pos2=3 → sum=6 ≥ 5 → **passes**.
    -   **Half-detected pairs** (one position only): `single_position_umi >= min_umi_threshold` (guide-level). Example: pos1=3, pos2=undetected → 3 < 5 → **fails**.
4.  **Classification**:
    -   `combo_type` (**cell-level**, each cell counted once):
        -   **`singlet`**: Cell contains exactly one unique guide pair above threshold.
        -   **`multiplet`**: Cell contains >1 unique guide pair above threshold (could indicate true multiplet infection, or potential recombination).
    -   `detection_status` (**per pair-assignment**, one row per pair):
        -   **`complete`**: Both `pos1` and `pos2` gRNAs were explicitly detected.
        -   **`inferred`**: Only one position was detected; the partner is inferred via iBAR mapping from the known pair design.

---

## Step 5: Build Combined AnnData

**Script**: [05_build_anndata.py](../scripts/05_build_anndata.py)

**UGER**: 1 job, 16GB RAM

Merges GEX expression data with gRNA assignments into a single analysis-ready AnnData.

> [!NOTE]
> **Split-Pipe Combine Fixes (Wrapper Script)**
> The `split-pipe` tool has two bugs/limitations when combining sublibraries:
> 1. It initializes an empty well list with `[] * 3` instead of `[[], [], []]`, causing an `IndexError` when merging mismatching sublibrary setups.
> 2. Re-defining sample pools (e.g. pooling C1 and C2 into "Control") normally requires a complete Step 1a rerun; retain the configurable 48-hour UGER limit for large inputs.
>
> The pipeline bypasses both issues using the [splitpipe_wrapper.py](../scripts/splitpipe_wrapper.py) script. This wrapper dynamically patches `run_proc_def.json` files on the fly to sync sample names with `config.yaml`, and then monkey-patches `spsamp.SpSamp.__init__` in memory before launching `split-pipe`, saving roughly 24 hours of re-computation.

```mermaid
flowchart TD
    GEX[("gex_combined/<br>adata.h5ad")] --> LOAD_GEX["Load GEX AnnData<br>(cells × genes)"]
    ASSIGN[("grna_cell_<br>assignment.csv")] --> LOAD_ASSIGN["Load gRNA<br>assignments"]
    MTX[("grna_umi_<br>matrix.mtx")] --> LOAD_MTX["Load gRNA<br>UMI matrix"]

    LOAD_GEX --> MATCH["Match cell barcodes<br>(bc_wells overlap)"]
    LOAD_ASSIGN --> MATCH

    MATCH --> OBS["Add to adata.obs:<br>• combo_type / detection_status<br>• pair_id / target_combo<br>• pos1 & pos2 info (detected)<br>• inferred pos1 / pos2<br>• has_grna"]

    LOAD_MTX --> ALIGN["Align matrix to<br>adata.obs.index"]
    ALIGN --> OBSM["Add to adata.obsm:<br>• gRNA_counts (sparse, per-gRNA)"]

    OBS --> META["Add pipeline metadata<br>to adata.uns"]
    OBSM --> META

    META --> OUT[("combined_adata.h5ad<br>cells × genes<br>+ dual-guide metadata")]

    style GEX fill:#e3f2fd
    style ASSIGN fill:#e3f2fd
    style MTX fill:#e3f2fd
    style OUT fill:#c8e6c9
```

### AnnData Structure

| Layer | Content |
|---|---|
| `adata.X` | GEX count matrix (cells × genes) |
| `adata.obs` | Dual-guide metadata: `target_combo`, `pair_id`, `combo_type`, `detection_status`, `pos1_grna`, `pos2_grna`, `inferred_pos*_grna`, `has_grna` |
| `adata.obs` | Sublibrary metadata: `sublib_index` (Parse combinatorial ID), `sublibrary` (e.g., C1/T1/C2), `sample_specific` |
| `adata.obsm['gRNA_counts']` | Sparse UMI count matrix (cells × independent gRNAs) |
| `adata.uns['gRNA_names']` | Guide name list matching matrix columns |
| `adata.uns['guide_annotations']` | Design-level annotations per individual guide |
| `adata.uns['ccperturbseq_pipeline']` | Pipeline parameters snapshot |

---

## Configuration Reference

Key sections of [config.yaml](../config.yaml):

| Section | Key Parameters |
|---|---|
| `gex_decontamination` | pre-STAR removal, TSO/guide thresholds, raw reject retention, optional Step 1c/source merge, and existing-alignment sandbox |
| `grna_structure` | `fixed_seq_1` (TSO), `fixed_seq_2/3` (stems), `grna_lengths: [20, 21]`, `ibar_length: 7` |
| `grna_matching` | `max_lev_dist: 3`, `max_lev_dist_ibar: 2`, `save_failed_reads: false` |
| `guide_assignment` | `dual_guide_counting_mode: true`, `min_umi_threshold: 5`, `umi_dedup_ibar: true`, `max_umi_edits: 1`, `min_reads_per_umi: 1`, `compute_saturation_curve: true` |
| `barcode_extraction` | `bc_max_edit_dist: 1`, `n_processes: 8` |
| `libraries` | 4 sublibraries: C1, C2, T1, T2 (2 Control + 2 Treatment) |

## UGER Resource Allocations

| Step | Memory | CPUs | Runtime | Job Type |
|---|---|---|---|---|
| 1a (GEX) | 8GB/slot (128GB at default) | 16 | 48h | Array ×N |
| 1b (Combine) | 16GB/slot | 4 | 4h | Single |
| 1c (optional match) | 4GB/slot | 8 | 12h | Array ×N |
| 1c (optional aggregate) | 32GB/slot | 4 | 8h | Single |
| 2 (Barcodes) | 4GB | 1 | 12h | Array ×N |
| 3 (Matching) | 4GB/slot | 8 | 12h | Array ×N |
| 4 (Assignment) | 32GB/slot | 4 | 4h | Single |
| 5 (AnnData) | 32GB | 1 | 4h | Single |

Step 1a values come from `uger.step1a.slots`,
`uger.step1a.memory_per_slot_gb`, and `uger.step1a.runtime_hours`. The launcher
uses the slot count for both `-pe smp` and `-binding linear`, and split-pipe
derives its thread count from `NSLOTS`.

The 48-hour default is recommended for large or multi-lane GEX inputs. Verify
the resolved `-pe`, `-binding`, `h_vmem`, and `h_rt` values with
`submit/submit_pipeline.sh --dry-run`. More than 16 slots is not recommended:
the full-C2 evidence used eight Parse workers, averaged about 5.2 cores, and
already had ample memory headroom.

---

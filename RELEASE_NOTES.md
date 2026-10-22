# Release Notes

## v2.0 — Native kallisto and tunable GEX UMI correction (2026-10-22)

All maintained UGER stages now share optional YAML resources: CPU slots, memory
per slot, wall time, array concurrency, queue/OS/PE/account and site complexes.
Old configs retain their previous requests. Exact submitted IDs/resources are
recorded, and main-workflow holds use those IDs. See
[UGER resources](docs/UGER_RESOURCES.md).

The transcriptomic backend is selectable in Step 1a while Parse split-pipe
1.6.4 remains a dependency for preprocessing, barcode/sample handling and cell
calling. STAR remains the default. The Parse installation is unchanged.

- Added native kallisto BUS counting with Parse-compatible UMI families,
  sparse per-library products and a matching Step 1b combine path. Counts feed
  the existing guide/AnnData workflow without a pseudo-BAM projection.
- Added `gex_alignment.umi_correction: parse|directional|none` to **both STAR
  and kallisto**. Parse remains the default; directional uses the fixed-count
  `n_high >= 2*n_low - 1` rule and bustools-compatible traversal. Exact duplicate
  grouping remains when neighbor correction is disabled. Step 4 gRNA UMI
  correction remains independent.
- Gene output and unique gene assignment are default. Optional fractional
  gene assignment preserves noninteger values; `output_level: both` adds target
  matrices, and `transcript` selects targets as the primary expression level.
- Parallelized native counting by complete corrected-barcode shards, cached
  EC/gene annotation and ordered disk-backed merging. `count_threads: 0` uses
  allocated cores; a positive cap, including serial `1`, is supported. Workers
  read an immutable spool and merge into a separate database.
- Added reference/settings/output completion checks and STAR UMI policy markers
  to prevent mixed or stale quantification at combination/integration.
- GEX decontamination and optional Step 1c remain available with either backend;
  historical BAM reuse and coordinate-specific QC remain STAR-only.

The [quick start](docs/QUICKSTART.md), [installation](docs/INSTALLATION.md),
[configuration reference](docs/CONFIGURATION.md), [execution runbook](docs/RUNNING.md)
and [migration guide](docs/MIGRATING_TO_V2.md) cover colleagues with and without
shared rgupta access. Earlier operating guides now describe the backend boundary
and current worker entry points. Historical reports retain their measured
methods/results and point to the current guides.

The optional lab tools companion uses the stable public name
`/medpop/rgupta/Jiahao/conda/envs/csm_kallisto-share`, pointing to the validated
`csm_kallisto-share_20261022` build, with 60 checks
passed against the unchanged `splitpipe-share` environment. See
[shared installation](docs/INSTALLATION.md#shared-kallisto-companion-published-2026-10-22)
for the pinned versions and completion manifest. The targeted UGER regression
suite for resources, report templates and backend fixtures passed 144 tests and
49 subtests; shell/YAML examples and local documentation links were checked.

The [original four-library benchmark](docs/KALLISTO_PILOT_BENCHMARK_REPORT.md)
is the serial-counter baseline. The completed
[optimized benchmark](docs/KALLISTO_PARALLEL_BENCHMARK_REPORT.md) repeats all
36 backend/UMI-mode cases at 10 million pairs per library with 16 cores.
Kallisto counting takes 42–47 seconds per library, and preprocessing through
DGE output takes 6.59–7.16 minutes versus 8.37–8.70 minutes for STAR (medians
across four libraries; shared preprocessing is reused). All 96 tested count
products match their serial reference exactly. See
[the backend guide](docs/KALLISTO_BACKEND.md) for validation and output semantics.


This file records user-visible behavior changes. New users should begin with
the [README quick start](README.md); architecture is summarized in
[pipeline_overview.md](pipeline_overview.md), and implementation details are in
the [step-by-step workflow](docs/DETAILED_WORKFLOW.md).

## v1.5 — iBAR Preflight and Step 3 Scaling (2026-10-14)

This release adds early checks for possible iBAR reference misspecification
and improves Step 3 resource control and memory behavior for large libraries
and deep sequencing. Existing configs retain their matching settings and
resource defaults. See the [preflight and scaling guide](docs/SCALABILITY_AND_PREFLIGHT.md)
for usage, optional configuration, and the review of all pipeline stages.

**Early iBAR diagnostics**

- Step 3 automatically samples up to 10,000 reads per sublibrary before
  starting matching workers or opening production match files. It compares
  exact iBAR support for the stored and reverse-complemented reference both
  globally and conditional on an exact guide match.
- Flushed log warnings flag possible strand reversal, incompatible iBAR
  lengths, empty/invalid reference iBARs, low reference support, and low
  structure-extraction support. Shared guides, multiple iBARs per guide, and
  strand-neutral sequences are handled explicitly; sparse evidence is
  reported as inconclusive. No-iBAR mode skips the check.
- Each library writes `grna_matches/<library>_preflight.yaml` with sample
  counts, reference context, and diagnostic messages. Warnings are advisory:
  production continues using the configured reference and thresholds.
- Added `--preflight-only` for diagnostics without production matching and
  `--preflight-fastq` to check raw gRNA R1 before Step 2 or the full pipeline.
  Optional `grna_matching.preflight.enabled` and `sample_reads` settings
  control automatic checking; existing configs enable it by default.
- Diagnostics explain that Step 3 reads the existing reference CSV. Changing
  the Step 0 `design.rev_comp_ibar` flag alone does not repair that CSV, and
  no automatic reference correction is performed.

**Step 3 execution and resources**

- Bound submitted, unconsumed chunks to twice the worker count so completed
  results cannot accumulate with sequencing depth behind a slower gzip writer.
  Input parsing selects only the five required columns and preserves strings.
- Removed unused parent-process fuzzy matchers. Actual matchers still
  initialize once per worker. The distance-two fallback stops when two
  distinct references prove an ambiguity, preserving matching decisions.
- Progress is flushed after the first completed chunk and subsequently every
  ten chunks. The UGER worker uses unbuffered Python and pins numerical-library
  threads to one per matching process.
- The master launcher accepts optional `uger.step3.slots`,
  `memory_per_slot_gb`, and `runtime_hours`, with matching CPU binding and
  `NSLOTS`-based worker counts. Defaults remain 8 slots, 4 GB per slot
  (32 GB total), and 12 hours. Direct worker submissions retain their existing
  directives unless overridden at `qsub` submission.

**Compatibility and validation**

- Existing config files, library arrays, matching output columns, barcode
  correction, UMI grouping/collapse, and AnnData schemas are retained.
- Added a stage-by-stage scaling review identifying where additional cores,
  memory, or a future code refactor would help; larger Step 0/2/4/5 refactors
  are recommendations rather than changes in this release.
- Thirty synthetic regression tests passed in the production `splitpipe`
  environment on UGER, covering warning detection, standalone raw-R1 checks,
  bounded work submission, exhaustive fuzzy-match agreement, serial/parallel
  output agreement, resource overrides, and Step 1c compatibility.

## v1.4 — Pre-STAR GEX gRNA Decontamination (2026-07-17)

This release adds an optional, auditable GEX decontamination path while keeping
guide removal and guide identity assignment as separate decisions. See the
[v1.4 operating guide](docs/GEX_DECONTAMINATION_V1_4.md) for configuration,
outputs, and validation commands.

**Pilot contamination QC module**

- Documented the assay origin of GEX carry-over: universal reverse
  transcription uses both poly(A) and gRNA probes, then size selection
  separates long cDNA from short gRNA products. The dedicated gRNA fraction's
  later circularization-capture amplification adds the R1 N7 prefix; GEX
  carry-over lacks it. Completed C1/C2 scans measured 32.1–36.6 structural
  contaminant reads per million.
- Added a standalone, HPC-oriented QC workflow that scans the pilot 300-cycle
  paired FASTQs, estimates raw-read and corrected-cell/exact-UMI contamination,
  and benchmarks candidate filters after in-silico R1 truncation to 54 nt.
- Full-length structural truth uses TSO plus the downstream guide-cassette
  structure; TSO alone is never contaminant truth. Intact GEX-derived guide
  reads are expected at TSO offset 0, while offsets 1–2 only rescue rare 5' end
  damage/sequencing error. The dedicated gRNA N7 prefix is not assumed in GEX.
- Added matched full-versus-short54 controlled STAR arms. Truth is sampled by
  corrected `(CB,pB)` group to retain PCR/UMI families, and primary uncertainty
  intervals bootstrap those groups rather than treating reads as independent.
- Added an audit of historical split-pipe annotated BAMs, including unique and
  gene-assigned contaminant reads, guide-target concordance, and lower/upper
  molecule-impact bounds after split-pipe-equivalent Hamming-1 UMI collapse.
- Added code/config/input digests, truth-panel tokens, safety-cap checks, and
  `audit_completeness.csv` so partial or stale outputs cannot silently enter the
  aggregate report.
- Added per-library `benchmark_summary.yaml`, `bam_audit_summary.yaml`, filter
  metrics, controlled alignment summaries, cross-library CSVs/plots, and a
  generated audit report. See the consolidated
  [GEX gRNA contamination QC guide](docs/GEX_GRNA_CONTAMINATION_QC.md) for the
  rationale, execution steps, denominators, and field-by-field interpretation.
- The corrected eight-worker C2 QC scan processed all 363,646,678 pairs at
  68,089 pairs/s (5,340.7 seconds scan time), finding 11,683 structural
  contaminants and 10,559 with valid corrected barcodes. The bounded
  one-versus-four-worker benchmark reproduced all detailed artifacts and
  measured a 4.019x scan speedup. The launcher requests and binds the same
  number of cores; eight is the production default for full scans, while one
  core is used for the compact Step 1c arm.
- Controlled STAR job `60030670` aligned all 10,559 C2 valid-barcode
  structural contaminants in matched full and short54 arms. No full-length
  read mapped, whereas 8,113 short54 reads mapped, 6,965 mapped uniquely, and
  4,008 (37.96%) received a false gene assignment. Split-pipe-equivalent
  collapse produced 3,099 false gene molecules (37.62%). The v1.4 production
  decision removed 95.63% of false assigned reads and 95.61% of false
  molecules in this counterfactual.
- Normalized Parse gene-name serialization showed that 2,431 false reads and
  1,894 collapsed false molecules were assigned to their own intended guide
  target. Every target-concordant event was in the v1.4 removal stratum; the
  retained false assignments in this C2 panel were off-target.

**Runtime integration and removal**

- Added a checksum-validated, in-memory overlay for the split-pipe 1.6.4
  preprocessing and molinfo seams. The Parse source tree remains pristine;
  disabling `gex_decontamination.enabled` invokes the original executable.
- Added a packed, exact-first Hamming-1 guide filter designed for references
  containing 20–50k unique guides. It uses complementary seed tables plus
  full-guide verification rather than materializing every mismatch neighbor.
- A read is removed only when it contains a full/fuzzy 30-bp TSO plus a complete
  reference-supported guide. Offset 0 is expected for an intact GEX-derived
  contaminant; offsets 1–2 only tolerate rare 5' damage/sequencing error. Indels
  and more than one guide mismatch are not accepted by this removal backend.
- Filter-positive, valid-barcode reads are omitted before human STAR. Step 1a
  writes compact corrected-barcode evidence for Step 1c during the same pass;
  synchronized raw paired evidence is retained by default for provenance, QC,
  and the historical-BAM route.
- Promoted the precomputed edlib TSO matcher, exact/pigeonhole fast path,
  compact guide routing, and fused fixed-counter hook. Full-C2 production prep
  retained all removal invariants and completed in 3,215.796 seconds versus the
  immutable 2,923.947-second pristine baseline (9.98% overhead).
- Step 1a resources now default to 16 slots, 8 GB per slot, and 48 hours and are
  configurable under `uger.step1a`; `-pe smp` and linear binding are always
  emitted with the same slot count.

**Reproducible environments**

- `requirements_extra.txt` is the single add-on installation command for the
  splitpipe environment; `edlib>=1.3` is required when decontamination is
  enabled.
- All production pipeline, v1.4 decontamination, controlled-alignment, and QC
  jobs use the configured `splitpipe` environment. `csm_analysis` and
  `environment_qc.yml` are reserved for separate downstream analysis modules
  and notebooks. UGER wrappers activate the production environment inside each
  job, and both master launchers support dry-run verification.

**TSO survivor policy and auditability**

- Filter-negative reads with a full/fuzzy TSO are now trimmed at the detected
  full-TSO endpoint. Reads without a full-TSO call retain Parse's exact terminal
  17-mer fallback; reads matching neither rule remain unchanged.
- The full-TSO call by itself is never sufficient for removal, preserving short
  cDNA fragments that legitimately contain full or partial TSO sequence.
- Manifest schema 4 records the survivor policy
  `full_tso_endpoint_then_parse_exact_17mer`, full-TSO/guide statistics,
  survivor trim routes, exact input identities, synchronized artifact hashes,
  and reference/removal-policy checksums. Incompatible earlier preprocessing
  manifests are rejected and regenerated.

**Short reads, Step 1c, and source merging**

- A 54-bp R1 can support removal with the standard 30-bp TSO and 20/21-bp
  guides at expected offset 0 and damage/error-tolerant offsets 1–2. It is
  generally too short to recover the full
  dual-guide stem+iBAR structure; Step 1c warns and does not invent identity.
  Single-guide/no-iBAR identity continues to use the v1.2 short-read logic.
- Added optional Step 1c processing that applies ordinary Step 3 directly to
  Step 1a's compact full-run corrected-barcode evidence. It neither rescans the
  GEX FASTQs nor recalibrates barcodes from rejected reads. Ordinary matches are
  retained for audit; a stricter guide/iBAR subset feeds Step 4 in a separate
  output namespace.
- On full C2, Step 1c reduced 363,646,678 input pairs to 9,511 Step 1a evidence
  rows, obtained 9,060 ordinary Step 3 matches, and retained 9,014 under the
  conservative dual-guide policy. The one-core task completed in about 5.4
  seconds; these read-evidence counts precede UMI collapse.
- Added optional dedicated-plus-GEX evidence union before
  `min_reads_per_umi` filtering and directional UMI collapse. Source provenance,
  overlap, GEX-only, dedicated-only, and incremental-GEX statistics are kept.
- In C2-only QC, 6,788 collapsed GEX UMI clusters added only 693 clusters to
  the 789,746-cluster dedicated stream. Exactly 5,995
  `(CB,guide,iBAR,UMI)` keys overlapped across sources and all also matched pB.
  Four cells gained a guide-pair assignment and two changed. In the separate
  four-library AnnData integration test, the dedicated-only total was
  reproduced at 3,272,782 UMIs and the combined total was 3,273,475.
- Step 5 can retain GEX-only, dedicated-only, and combined sparse matrices,
  while preserving dedicated-only calls for comparison when combined calls
  drive the final result.

**Existing split-pipe runs**

- Added a sandboxed route for compatible historical runs: preprocess the raw
  FASTQs, remove every positive QNAME from the historical name-sorted/annotated
  BAM, install a harmonized BAM in a new output tree, and resume split-pipe at
  `--mode mol`. The historical output is never overwritten and STAR is not
  rerun.
- The BAM route cannot retroactively apply the new full/fuzzy-TSO survivor
  trimming because those sequences were already aligned. A fresh alignment is
  required when that trimming change is part of the test objective.
## v1.3 — Failed Read Recording & Sequencing Saturation QC (2026-03-31)

**Step 3 — Save Failed Reads:**
- Added `save_failed_reads` option under `grna_matching`. When enabled, reads that fail Step 3 matching are written to `grna_matches/<lib>_grna_failed.tsv.gz` with a `failure_reason` column. Default: `false`.

**Step 4 — UMI Read-Depth Filtering & Saturation QC:**
- Added `min_reads_per_umi` (default: `1`) under `guide_assignment`. Filters out UMIs supported by fewer reads than the threshold before directional collapsing — removes noise from singleton UMIs in deeply-sequenced libraries.
- Added `compute_saturation_curve` (default: `true`). Computes a CellRanger-style sequencing saturation downsampling curve via Poisson approximation.
- **New QC plots**:
  - `reads_per_umi_distribution.png` — histogram + reverse CDF of reads-per-UMI, with threshold annotation.
  - `saturation_curve.png` — UMI discovery & saturation % vs subsampled read depth, per-library + combined.
  - `median_reads_per_umi_per_cell.png` — violin plot of per-cell sequencing depth by sublibrary.
- **New QC stats** in `grna_qc_stats.yaml`: `sequencing_saturation`, `total_umis_raw`, `umis_filtered_low_reads`, `median_reads_per_umi`.

## v1.2 — Short Read Support (2026-03-20)

> [!IMPORTANT]
> This update applies only to **single-guide mode without iBAR** (`single_guide_mode: true`, `ibar_length: 0`). Dual-guide and iBAR-enabled pipelines are unaffected.

**Problem**: In libraries with shorter sequencing reads (e.g., R1 ≈ 50–55bp), R1 may be too short to capture the full `TSO + gRNA + stem` structure. Previously, these reads were silently discarded by the fixed-length `ibar_end > len(r1_seq)` boundary check, even though the gRNA sequence itself may be fully or partially recoverable.

**Fix — Four-scenario graceful degradation in `extract_grna_single_short()`:**

| Scenario | Condition | Behavior |
|---|---|---|
| 1. Full stem | R1 ≥ TSO + gRNA + stem | Normal extraction; highest confidence |
| 2. Partial stem | TSO + gRNA ≤ R1 < TSO + gRNA + stem | Partial stem boundary check with strict mismatch rules; shortest matching gRNA length preferred on ambiguity |
| 3. Truncated gRNA | TSO < R1 < TSO + min(gRNA) | Rescue all bases after TSO; matched against truncated reference |
| 4. Discard | No TSO or R1 ≤ TSO | Read discarded |

**Partial stem mismatch rules:**
- Partial stem ≤ 2bp: Exactly 0 mismatches allowed (overrides `fixed_seq_max_mismatches`)
- Partial stem > 2bp: 1 mismatch allowed (if user config ≥ 1)
- A one-time warning is printed if the applied thresholds differ from the user config

**Truncated gRNA matching:**
- The pipeline automatically samples the first 100,000 reads per library to detect all naturally-occurring gRNA lengths (no user input required)
- Detected lengths are logged alongside the expected lengths from `config.yaml`
- Pre-built `NeighborhoodMatcher` dictionaries are constructed for each detected truncated length during worker initialization, avoiding per-read computational overhead
- Truncated references may collapse uniqueness — reads matching multiple references at the same distance are flagged as `grna_ambiguous`

**New output columns** (single-guide short-read mode only):

| Column | Type | Description |
|---|---|---|
| `stem_truncated` | bool | `True` if stem was partially or fully absent |
| `grna_truncated` | bool | `True` if gRNA itself was truncated (Scenario 3) |
| `observed_grna_length` | int | Actual extracted gRNA length |

**New Step 3 stats:**

| Stat | Description |
|---|---|
| `short_read_scenario_1/2/3` | Read counts per scenario |
| `short_read_discarded` | Reads discarded (Scenario 4) |
| `truncated_matcher_used` | Reads routed to a truncated matcher |
| `grna_ambiguous_truncated` | Ambiguous matches caused by truncation |
| `grna_no_matcher_for_length` | Reads with no pre-built matcher for their length |

## v1.1 — Multi-FASTQ Support (2026-03-19)

**Changes:**
- **Graceful Multi-FASTQ Handling**: Added native support for processing multiple lane FASTQ files corresponding to a single sublibrary. Users can now provide a `gex_prefix` or `grna_prefix` that utilizes standard shell wildcard globbing (e.g., `library_S*_L00*` to match `L001`, `L002` lane files simultaneously). The pipeline locates matching files using the pattern `*{prefix}*_R1_*.fastq.gz`. To verify your pattern works, test it with: `find /path/to/fastq_dir -name "*<prefix>*_R1_*.fastq.gz"`. The pipeline aggregates all matched files seamlessly without overwriting outputs.
- **Step 1a (split-pipe) Update**: Replaced `head -1` discovery. Step 1a now groups the full list of matched GEX `_R1_` and `_R2_` FASTQ files and passes them simultaneously into the `split-pipe` execution command.
- **Step 2 (Barcode Extraction) Update**: Refactored the internal file loop. The script now establishes a single `.tsv.gz` writable stream per sublibrary and continuously appends the barcode extraction data across all matching gRNA FASTQ pairs, eliminating catastrophic overwrites.

## v1.0 — Single-Guide & No-iBAR Modes (2026-03-17)

**Changes:**
- **Single-Guide Mode**: Added robust support for single-guide libraries via `single_guide_mode: true` under `grna_structure` in `config.yaml`. The pipeline safely bypasses dual-specific matching steps (stem2, tracr2, inferred pairing) and employs a streamlined `assign_single_guides()` assignment algorithm.
- **Optional iBAR Support**: The pipeline now seamlessly supports processing libraries without iBARs. Setting `ibar_length: 0` natively skips all iBAR-based extraction, fuzzy matching, disambiguation, and UMI grouping logic.
- **Consolidated Metadata**: The final AnnData builder (`05_build_anndata.py`) intelligently aligns detection semantics between dual and single mode, populating compatible `combo_type`, `n_guides`, and unified single-guide metadata mappings.
- **Revised UGER Submission Wrapper**: Now support custom conda environment names, which can be set under the `environment` section in `config.yaml`.

## v0.4 — iBAR Unresolved Handling (2026-03-16)

**Changes:**
- **Retaining Unresolved iBARs**: For gRNAs where `ibar_seq_corrected` could not be resolved, the pipeline now keeps the raw iBAR sequence in a new mutually exclusive column `ibar_seq_unresolved`.
- **Accurate UMI Deduplication**: UMI deduplication now groups reads utilizing both `ibar_seq_corrected` and `ibar_seq_unresolved`. This ensures that reads with distinct, unresolved raw iBAR sequences are correctly maintained as separate molecules rather than being falsely collapsed together into an empty iBAR bin, thus retaining quantitative accuracy.
- **Output Format Update**: The `grna_umi_counts.csv` post-deduplication output file format is now saved directly as `grna_umi_counts.csv.gz` to improve storage efficiency across large single-cell datasets.

## v0.3 — Robust Matching (2026-03-10)

**Changes:**
- **Ambiguous Stem Handling**: When reads have equally good matches to both stem variants (`score_A == score_B`), they are now explicitly flagged as `AMBIGUOUS`. Previously, they arbitrarily tie-broke to `pos1`.
- **Strict Disambiguation**: For shared gRNAs, if the stem is `AMBIGUOUS` and `iBAR` cannot resolve it, the read is now marked unresolved (preventing silent mis-assignment).
- **Case-Insensitive Reference**: The pipeline now robustly handles user-provided reference CSVs (`grna_ibar_reference.csv`) containing lowercase sequences. All reference and query sequences are internally normalized to uppercase during extraction.

## v0.2 — Shared gRNA Disambiguation (2026-03-04)

> [!IMPORTANT]
> This release fixes a critical bug where gRNA sequences shared across positions (e.g., TLNRD1_TF vs CCM2_TLNRD1) were silently mis-assigned due to `NeighborhoodMatcher.exact` dict using last-write-wins.

**Problem**: `NeighborhoodMatcher.exact` stored `seq → single_index`, silently overwriting earlier entries for duplicate sequences. Exact matches for shared gRNAs returned an arbitrary position based on reference row order. This also caused false decoupling detection and inflated `cross_swap` rates in Q3.

**Fix — Three-layer disambiguation:**

| Layer | Signal | Method |
|---|---|---|
| 1. Exact combined lookup | `(grna_seq, ibar_seq)` | Already unique — resolves ~60-70% of reads |
| 2. iBAR-based (fuzzy path) | Within-guide iBAR match | For each candidate `grna_name`, try iBAR match against its iBAR set. 0 iBAR overlaps across positions → exactly one candidate matches |
| 3. Stem+tracr variant (fallback) | Combined stem+tracr Hamming | `extract_grna_ibar()` now scores all 4 stem+tracr combos jointly; filter by `position == variant` |

**Changes:**

- **`NeighborhoodMatcher`**: `exact` dict stores lists of indices; `match()` returns `list[str]` for position-ambiguous matches
- **`extract_grna_ibar()`**: Combined stem+tracr anchor scoring (`score = stem_mm + tracr_mm`); reports `anchor_combo` (A/B/cross) for QC stats
- **`_process_rows()`**: Uses `_resolve_shared_grna()` — tries iBAR disambiguation first, falls back to stem variant
- **`00_format_reference.py`**: Detects shared gRNA sequences across positions, checks iBAR uniqueness
- **`config.yaml`**: `fixed_seq_max_mismatches` documented as combined stem+tracr threshold (safe limit: ≤ `⌊(combined_hamming - 1) / 2⌋ = 2`)

**New Step 3 stats:**

| Stat | Description |
|---|---|
| `anchor_combo_A` / `_B` / `_cross` | Stem+tracr pairing distribution |
| `shared_grna_hit` | Reads matching a shared gRNA sequence |
| `shared_grna_ibar_resolved` | Resolved via iBAR disambiguation |
| `shared_grna_variant_resolved` | Resolved via stem variant fallback |
| `shared_grna_unresolved` | Could not be disambiguated (expected ~0) |

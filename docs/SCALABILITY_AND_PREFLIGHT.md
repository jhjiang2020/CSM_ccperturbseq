# iBAR preflight and pipeline scaling

Step 3 checks iBAR specification before starting matching workers or opening
the production match files. Old configs enable the check without adding keys.
The check uses exact set lookups and the existing structure extractor, rather
than the expensive global fuzzy matcher. It samples at most 10,000 reads per
sublibrary and writes `grna_matches/<library>_preflight.yaml` as well as a flushed
summary and any warnings in the job log. Counts are read-weighted.

## Check before a full production run

Run the standalone raw-R1 check in a compute allocation, before submitting the
full pipeline. It needs the formatted reference and raw gRNA FASTQs; Step 2
outputs and GEX processing are unnecessary. It samples the matching R1 files
in sorted order until the per-library cap is reached. FASTQ discovery follows
the same `fastq_dir`/`grna_prefix` glob as Step 2. It does not validate R2 cell
barcodes or sample every lane independently.

On UGER, after initializing the scheduler:

```bash
PROJECT=/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq
qrsh -q broad -l os=RedHat8,h_vmem=8G,h_rt=0:30:00 -binding linear:1 \
  /medpop/rgupta/Jiahao/conda/envs/splitpipe/bin/python -B -u \
  "$PROJECT/scripts/03_match_grna_ibar.py" --config "$PROJECT/config.yaml" \
  --preflight-only --preflight-fastq
```

Add `--library <name>` to check one sublibrary. To use an existing Step 2 TSV,
omit `--preflight-fastq`; `--barcode-input <path>` accepts an explicit TSV for
one library. `--output-namespace <path>` can separate diagnostic records from
production outputs. `--preflight-only` writes only diagnostic records and exits
without producing match files or matching statistics. Missing preflight input
is an error. Normal production also repeats the bounded TSV preflight at the
start of each Step 3 library, including Step 1c calls to the same matcher.

## What the warning means

The check compares observed iBARs with the stored reference and its reverse
complement, both globally and conditional on an exactly matched guide.
Shared guides and multiple iBARs per guide use sets. Same-guide evidence can
detect reversal even when the global pool contains both strands. Sequences
matching both strands, including palindromes, do not count as strand-specific
evidence.

A possible reverse-complement warning requires at least 100 informative
extractions or exact-guide reads, at least 20 reverse-complement-only hits,
at least a fivefold excess over stored-only hits, and a difference of at least
10 percentage points. These are conservative alert heuristics, not a
statistical test. Other warnings cover incompatible reference iBAR lengths,
empty/invalid iBAR sequences, fewer than 10% exact hits in either orientation,
and fewer than 10% successful structure extractions. Sparse or strand-neutral
evidence is reported as inconclusive. No-iBAR mode is skipped.

An orientation warning asks the user to check `ibar_seq` in the actual
`paths.grna_reference` CSV. Step 3 does not read `design.rev_comp_ibar`;
changing that Step 0 setting alone does not repair an existing CSV. The check
does not change sequences, matching thresholds, guide calls, or reference
files, and production continues after a warning. Use the standalone check if
you want to review warnings before committing production compute time.

First-read samples can be biased, and low support can also reflect read
quality, wrong fixed anchors, a different reference/library, or real
guide–iBAR uncoupling. The diagnostic does not establish a biological cause.

Optional settings, with defaults shown:

```yaml
grna_matching:
  # Keep existing matching keys here.
  preflight:
    enabled: true
    sample_reads: 10000
```

`sample_reads` accepts 1–100,000; samples below the minimum evidence thresholds
may remain inconclusive. `--preflight-only` explicitly runs the check even if
automatic preflight is disabled in the config.

## Resource controls with existing configs

The launcher now supports optional Step 3 resources alongside `uger.step1a`:

```yaml
uger:
  # Preserve existing step1a settings.
  step3:
    slots: 16
    memory_per_slot_gb: 4
    runtime_hours: 48
```

This example reserves **64 GB total**, requests matching `-pe smp 16` and
`-binding linear:16`, and runs 16 workers through `NSLOTS`. It is a tuning
example, not a new default. Omitting the block retains **8 workers, 32 GB total,
12 hours**. Validate the request with:

```bash
bash submit/submit_pipeline.sh --config config.yaml --step 3 --dry-run
```

The optional YAML resource override is applied by `submit_pipeline.sh`.
Direct `qsub submit/submit_step3_match.sh` retains the worker's directives;
use explicit `qsub -pe`, `-binding`, and `-l` overrides for a direct submission.
Existing library arrays, job dependencies, matching config keys, and output
columns remain in place. The default config file needs no edits.

## Scaling review

This is a review of the current code paths, not a new dataset benchmark.
Read depth, reference cardinality, unique molecules, and cell count stress
different resources. Additional requested slots help only phases that use them.

| Stage | Scaling constraint in the current code | Resource that can help | Smallest useful next change |
|---|---|---|---|
| 0: reference formatting/validation | Global guide and iBAR pairwise distance scans are quadratic in unique reference sequences; close-pair lists can also grow very large. These loops are serial. | More RAM may prevent exhaustion; extra slots alone do not parallelize the loops. | Offer a bounded/sample mode for descriptive global distance statistics or an indexed overlap check, while retaining exact within-guide validation. Keep this separate from production counts. |
| 1a: GEX split-pipe | One array task per sublibrary; the wrapper passes `NSLOTS` as `--nthreads`. Parse 1.6.4 preprocessing has an eight-worker cap; STAR/sorting may use more. | More RAM for alignment/sorting; extra cores when the active phase uses them. | Use existing `uger.step1a` controls and observed CPU/memory usage. The documented release benchmark already supports keeping 16 slots as the starting point. |
| 1b: GEX combine | Passes four allocated threads to split-pipe; combined sample data and temporary products increase with cell/feature counts. | RAM and filesystem throughput; extra cores only for threaded vendor phases. | Override resources at submission after examining utilization; do not assume linear CPU scaling. |
| 1c: GEX-derived guides | Reuses Step 3; effective workers are capped by evidence rows/chunk count. Compact evidence can be too small to occupy many workers. | More workers only for sufficiently large evidence; reference memory per worker. | Retain the existing worker cap and new Step 3 safeguards. |
| 2: cell-barcode extraction | One serial FASTQ-pair loop and gzip writer per array task; lane/file pairs are processed serially. `barcode_extraction.n_processes` does not parallelize this stage. | More concurrent sublibrary tasks help; additional slots in one task currently do not. RAM is mainly bounded sampling/buffers. | If needed, parallelize original lane/file pairs and merge TSVs/stats with a single header, preserving each original file's barcode-calibration sample. This needs a code change, not a CPU setting. |
| 3: guide/iBAR matching | Exact/one-edit lookup is indexed; remaining searches can scan the whole reference per read. Input decompression, TSV parsing, result serialization, and gzip writing run in the parent. Worker reference state is replicated. | More bound cores for CPU-heavy matching, with RAM increasing for worker/reference copies. More RAM alone does not accelerate a CPU scan. Serial gzip/I/O eventually limits CPU scaling. | First correct reference/anchor misspecification. Start with 8 workers; consider 16 when CPU utilization supports it. Use optional `uger.step3`, or lower existing CLI `--chunk_size` for memory/latency. Avoid requesting 64 workers by default. |
| 4: UMI deduplication/assignment | Input is chunked, but dictionaries retain all distinct group–UMI counts per library. Dedicated libraries run in parallel up to `min(ncpus, libraries)`; source-union paths and final assignment run serially. Directional UMI collapse compares UMIs within each group quadratically in the worst case. | RAM for unique molecules/groups; more cores up to library count for the dedicated path, with concurrent accumulators increasing peak RAM. | Raise `--ncpus` via allocated `NSLOTS` only for independent libraries and sufficient memory. A future disk-partitioned accumulator or per-library compact dedup stage would address memory beyond that. |
| 5: AnnData assembly | Loads sample AnnData objects and assignment tables, then concatenates; intermediate copies coexist. Guide matrices are sparse. | RAM sized to cells/nonzeros and intermediate copies; extra cores have limited benefit in this mostly serial wrapper. | Increase memory on the existing single-slot job when needed. Future work could release sample intermediates earlier or assemble backed data without changing output schema. |

Step 4 currently requests four slots at 32 GB/slot (**128 GB total**), while
Step 5 requests one slot at 32 GB. Increasing Step 4 from four to eight workers
at the same per-slot request would reserve 256 GB and may run twice as many
large UMI accumulators concurrently. Request memory from measured working set,
not from the slot count alone. In Step 5, an increased single-slot memory
request is more useful than reserving unused CPU slots.

## Changes included in this update

- Automatic bounded preflight, an optional raw-FASTQ preflight before production,
  flushed warnings, and structured per-library diagnostics.
- Optional Step 3 resource resolution, preserving current defaults and `NSLOTS`.
- At most twice the worker count of submitted, unconsumed chunks. This prevents
  completed results from accumulating with sequencing depth when output is
  slower than workers. Live memory still depends on chunk size, worker count,
  output size, and per-worker reference state.
- Input parsing selects only five required TSV columns, preserving string values.
- The parent no longer constructs an unused extra set of fuzzy matchers.
- The distance-two fallback stops after two distinct references prove a tie;
  exact matches, unique matches, and ambiguity decisions are preserved.
- First-chunk progress and explicit thread limits in the Step 3 worker wrapper.

The larger Step 0/2/4/5 refactors above remain recommendations. This update
does not change barcode correction, UMI grouping/collapse, downstream analysis,
or the matching output schema. Synthetic regression tests compare fuzzy
matches against exhaustive edit-distance search and compare serial/parallel
Step 3 outputs. No further collaborator-dataset benchmarks were run.

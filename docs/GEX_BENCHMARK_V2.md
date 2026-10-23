# v2 Step 1a benchmark: 50 million GEX read pairs per sublibrary

This is the current performance experiment, replacing the nine-configuration,
10-million-pair benchmarks. Superseded reports, plots and runners are retained
in Git history and the local UGER archive. Numerical performance and UMI-gain
claims below come from a completed `comparison/REPORT.md` and
`BENCHMARK_SUCCESS.json` from this runner.

<!-- benchmark-v2-run-status:start -->

The 50M-pair experiment completed and passed full cache validation. See the
[measured performance, agreement and UMI-gain report](GEX_BENCHMARK_V2_REPORT.md).
Run provenance and scheduler job IDs remain in the run's `logs/uger_submissions.jsonl` and atomic `logs/uger_submissions.records/` files.
The measured report includes complete scheduler coverage, memory/disk peaks and
a dedicated explanation of the STAR unique/fractional assignment caches.

<!-- benchmark-v2-run-status:end -->

## Dataset and controlled settings

- Use the four-library pilot, in order **C1, C2, T1, T2**: **50,000,000 paired
  GEX fragments per sublibrary**, 200 million pairs total. “50M reads” here means
  50M R1 records and their 50M R2 mates, before barcode extraction or filtering.
- Select equal-size deterministic prefixes from the resolved lane pairs, with
  the remainder allocated in lane order. Record names, lane quotas, source file
  sizes/mtimes and sampled FASTQ SHA256s. Pair-name mismatches, truncated records
  and lanes shorter than their quota fail. This is a reproducible matched sample,
  not a random subsample or an estimate of sampling uncertainty.
- Run the fixed **Parse split-pipe 1.6.4** barcode extraction, chemistry,
  corrected-barcode identities and sample masks once per library. Preserve the
  configured GEX guide-decontamination policy for all configurations. Step 1c,
  historical-BAM reuse and guide-source combination are disabled for this GEX
  experiment. Sampling time is measured separately from Step 1a.
- Keep the same GRCh38-2024-A FASTA/GTF annotation, gene output, forward strand,
  Parse cell-calling rules, minimum supporting-read threshold, and allocated
  cores. Reuse the validated full mature-plus-nascent kallisto reference; do not
  rebuild it for a different read depth. Reference construction is outside the
  benchmark. This controlled comparison does not reproduce Tony's unstranded
  reference/workflow.
- Use 16 bound CPU slots, 8 GiB per slot (**128 GiB total requested**) and
  12 hours per task initially; two concurrent tasks per array. Override the
  `uger.kallisto_benchmark` resource profile for another allocation. All phases
  use `NSLOTS`; nested numerical thread pools are pinned to one. Retain job IDs,
  host, affinity and actual allocation in the run record.

The twelve scientific configurations are the Cartesian product
`(star | kallisto) × (unique | fractional) × (parse | directional | none)`.
Their identifiers are explicit, e.g. `star_unique_parse` and
`kallisto_fractional_none`. Every configuration has its own output directory.

## Cache graph and measured steps

There are **two aligner executions per library**, one STAR and one kallisto:
eight alignment executions across the entire experiment. Twelve complete
Step 1a runs per library are unnecessary.

| Measured stage | Executions per library | Retained intermediate / reuse |
| --- | ---: | --- |
| Raw paired sampling | 1 | Matched 50M-pair compressed FASTQs; preparation cost only |
| Parse barcode extraction and configured guide filter | 1 | Immutable `barcode_head.fastq`, barcode data, run definition; used by both aligners |
| STAR alignment, including index load | 1 | Unsorted genomic BAM; retain all hits up to `max_multimaps: 3`, `--outSAMmultNmax -1` |
| STAR queryname sorting / postprocessing | 1 | Complete pre-molinfo sorted BAM, preserving secondary hits |
| STAR unique gene annotation | 1 | Primary-read BAM with Parse-selected `GX` gene and `RE` exon flag; used by all three unique UMI modes |
| STAR fractional gene annotation and evidence spool | 1 | Complete-read candidate gene sets in SQLite plus feature/barcode metadata; used by all three fractional UMI modes |
| Kallisto input formatting | 1 | cDNA/tag FASTQs, barcode lookup, BUS lengths and exclusion counters |
| Kallisto pseudoalignment, including index load | 1 | Persistent BUS, EC membership, target list and run information |
| Kallisto target/gene annotation and evidence spool | 1 | SQLite read evidence and cached target/gene/barcode metadata; used by all six kallisto policies |
| UMI correction, molecule aggregation and matrix assembly | 12 | Per-policy molecule assignments and counts; never rerun an aligner |
| Parse cell calling and DGE publication | 12 | Per-policy called/unfiltered and all-well gene products |

STAR unique annotation calls the installed Parse `get_gene` with its MAPQ gate,
strand/containment checks and overlapping-gene ranking. During subsequent native
Parse molecule runs, an in-memory `get_gene` adapter returns the cached `GX/RE`
values. Native molecule grouping, correction kernel, supporting-read threshold,
exon flags and DGE remain Parse's implementations. Secondary alignments are
excluded only from this unique cache. The complete sorted BAM is retained for
fractional assignment. Vendor files and the source alignment are unchanged.

The unique-gene cache is constructed by one streaming annotation process. In a
normal production STAR-unique run, Parse annotates reads inside its parallel
molecule workers. The reconstructed benchmark path therefore includes explicit
cache materialization and a different division of work; it is not a direct
production end-to-end runtime or a pure aligner speed comparison. The cached
policy timings measure the additional work after that cache exists.

STAR fractional annotation uses the production `spool_bam`: original-read
identity, primary/NH completeness checks, all retained loci and all overlapping
genes. Kallisto uses the production BUS reader and target annotation. Both feed
the production shared molecule counter, with exact-UMI/candidate resolution and
the selected neighbor-correction policy. Candidate geometry/EC-to-gene mapping
is cached; policy-dependent exact-UMI/collision resolution and correction remain
inside counting. Each policy stages a physical copy of the sealed read spool
into its job-local temporary directory, opens it read-only, and keeps count
databases and shards there too. The copy's time, bytes and temporary disk space
are part of the measured molecule phase. Direct shared-filesystem SQLite reads
were found to add a storage bottleneck; those initial successful measurements
are retained under `timing_baselines/shared_sqlite_v1` and excluded from the
final policy timings. Retiming requires exact matrices, ordered metadata,
scientific counters and molecule-assignment bytes against every completed
original molecule product. Gene weights stay
floating point. Native STAR accepts its usual UMI alphabet; BUS/fractional STAR
exclude non-ACGT UMIs with diagnostics, as in production.

The native STAR molecule phase includes its split/merge and all-well outputs;
the shared-counter phase includes its shard merge, matrix construction and
matrix-cache serialization. DGE publication includes native cell calling and
output I/O. Reference initialization and process startup occur within each
measured command. The adapter's cache creation cost is explicit. These are
measurements of implemented stages, not isolated algorithm timings.

## Runtime and resource accounting

Each stage uses **per-child `wait4` resource accounting**, avoiding the old
runner's cumulative `RUSAGE_CHILDREN` CPU/RSS measurements. Usage includes the
child and descendants it waited for; no external GNU time installation is needed. The resource CSV
records wall seconds, user/system CPU seconds, CPU-hours, reserved core-hours,
largest-process RSS, process-tree RSS sampled every two seconds, page faults,
filesystem block-I/O counters, output sizes, and sampled temporary/working-file
sizes. It also records requested slots/memory/time, hostname, affinity, command,
UTC start and scheduler job/task IDs.

Interpret the measurements precisely:

- `max_process_rss_kb` is the largest individual-process RSS high-water mark.
  Sampled `peak_tree_rss_gb` sums descendants and counts shared pages more than
  once; short peaks can be missed. Neither is scheduler `maxvmem`.
- `peak_temp_bytes` covers the stage's isolated `TMPDIR`.
  `peak_working_file_bytes` covers output directories, including linked inputs;
  it is apparent file size, not extra physical disk allocation. Hard links are
  deduplicated within a sample. Samtools temporaries alongside its output are
  included; these sampled peaks can miss short-lived files.
- Kernel filesystem I/O counts are block operations, not FASTQ/BAM logical bytes;
  filesystem caching and host/storage contention can affect them.
- `reconstructed_step1a_seconds` sums preprocessing + backend formatting,
  alignment/sort + the appropriate assignment cache + the selected UMI/matrix
  and DGE phases. It excludes raw sampling, reference build, Step 1b,
  clustering/report generation, checksums between stages and orchestration.
- `cached_policy_seconds` includes only the last two phases, the additional
  work for a setting after alignment/assignment caches exist.
- Actual experiment cost sums each executed shared stage **once**. Summing the
  twelve reconstructed paths would repeatedly bill the shared work. Summed task
  wall time differs from elapsed makespan under concurrent scheduling.

Collect `qacct -j <job-id>` after completion into the run's logs; verify every
array task's `failed` and `exit_status` and retain `ru_wallclock`, CPU and
`maxvmem` separately. One observation per stage/library yields observed times
without timing confidence intervals. Report library-specific measurements and
medians; the four different libraries are not technical timing replicates.

After completion, collect scheduler metadata and check coverage separately from
scientific-stage timing (both helpers read only compact bookkeeping files):

```bash
CSM_BENCH_DIR="${RUN}" bash submit/collect_gex_benchmark_accounting.sh
timeout 30s python scripts/benchmark_gex_scheduler_report.py "${RUN}" --require-complete
```

The scheduler export merges atomic per-job records and JSONL records by job ID.
It includes executed failures, successful cache-hit workers, setup, validation,
reporting and publication. Explicitly recorded tasks canceled before starting
are listed separately; unresolved missing accounting fails the final check.
Scheduler `wallclock × slots` measures reserved core-hours and includes overhead
outside the per-stage commands. `maxvmem` remains virtual memory.


## Agreement and gains

Use **every cell called by STAR + unique + Parse correction** as the primary,
fixed cohort within each sublibrary. Retrieve every other configuration's counts
from its all-well matrix; zero-fill absent cells and align gene IDs explicitly.
Never estimate gains by comparing medians from different called-cell sets.

The report produces:

1. **All 66 configuration pairs per library** (264 rows): cell-call intersection,
   union and Jaccard; fixed-cohort UMI Pearson/Spearman; per-cell log-expression
   Pearson including shared zeros; fixed-cohort pseudobulk log-expression
   agreement. A separate common-called-cell expression metric shows cohort
   sensitivity. Do not concatenate libraries without preserving library identity.
2. Called-cell yield, fixed-cohort UMI totals/genes and signed UMI changes versus
   the baseline for all 48 library/configuration products.
3. Factor-matched transitions: unique → fractional holding aligner and UMI mode
   fixed; Parse → directional, Parse → none and directional → none holding
   aligner/assignment fixed. Report absolute UMI change, percentage change in
   total fixed-cohort UMIs, median per-cell delta/ratio, and cells gaining,
   losing or unchanged. No correction still deduplicates exact UMIs.
4. Per-cell UMI tables and per-gene gains/losses against the baseline, preserving
   fractional weights. These expose changes concentrated in a subset of cells
   or genes. Agreement and gain heatmaps accompany the machine-readable tables.

Parse and directional correction are **not a guaranteed relaxation ordering**:
they use different abundance/tie/traversal rules. Retain negative deltas.
Likewise unique → fractional changes assignment and collision resolution and
need not increase every cell/gene. Additional counted UMIs can include ambiguous
assignments or retained sequencing errors; they are not established biological
recovery. Do not multiply 50M-pair gains into a full-depth prediction.

## Launch, resume and review

Run regression checks on a compute node before launching the full experiment.
The new tests cover the cache graph, stale/partial cache rejection, fixed-cohort
alignment, signed gains, all-pairs reporting and noncumulative resource records.
The production backend fixtures cover the twelve scientific policy choices.

On UGER, use an absolute project path and existing validated tools:

```bash
cd /medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq
source /broad/software/scripts/useuse
reuse -q UGER
BASE="${PWD}/output/kallisto_benchmark_10m_20261022"
RUN="${PWD}/output/gex_benchmark_v2_50m_20261023"
bash submit/submit_gex_benchmark_v2.sh \
  --config "${PWD}/config.yaml" --work-dir "${RUN}" \
  --reference "${BASE}/reference" \
  --kallisto "${BASE}/tools/kallisto/kallisto" --dry-run
# Repeat without --dry-run to submit.
```

The dependency chain is `init → prepare[1–4] → align[1–8] → assign[1–8] →
run[1–48] → report → publish`. Publication copies validated tables/figures and
the measured report to `docs/GEX_BENCHMARK_V2_REPORT.md` and updates the status
above. It refuses incomplete runs or reports without full cache checksums.
The `kallisto_benchmark` resource profile is retained for
compatibility; it now controls the complete v2 benchmark. The compute worker is `submit/run_gex_benchmark_v2.sh`.

Initialization seals settings, transitive GEX dependency/vendor hashes, annotation identity,
kallisto executable and reference manifest, and verifies the reference artifacts.
Local imports and executed Python helpers are followed recursively; unrelated
reference-building or gRNA workflow edits do not change this seal. Each successful
stage records dependency-manifest hashes and artifact
checksums/sizes/mtimes. Successful caches are reused without another timing
measurement; changed scientific settings require a new root. A reviewed
implementation repair can use `accept-revision --repair-note <reason>` with
the same sealed allocation and scientific settings. It preserves the previous
contract, snapshots the revised sources and records both code hashes in
`code_revisions.jsonl`; successful caches retain their original measurements.
The initial sources are preserved under `code_revisions/original`. To keep a
long-running experiment independent of edits in the shared checkout, workers
can set `CSM_BENCH_SOURCE_ROOT` to a saved revision directory. The worker runs
that snapshot and prepends it to `PYTHONPATH`; a saved raw YAML must reproduce
the sealed settings exactly. The 2026-10-23 continuation uses revision
`a936e58ba23d90ec` and `frozen_config.yaml`. Only known code-seal publication
failures are archived, after their workers terminate; successful molecule
products and measurements remain intact. Normal reuse checks
size/mtime; `--verify-cache` rehashes files. The report worker verifies checksums.
Every consumer validates completion records because scheduler dependencies do
not imply predecessor success. A partial `.started` stage refuses ordinary reuse.
The `retry-star` repair archives a confirmed failed alignment's output directory,
logs, resource metrics and marker under `failed_attempts/`, then restarts only
STAR from successful preprocessing. Other partial stages need an inspected,
recorded repair or a new experiment root. No automatic partial overwrite can turn
a skipped vendor stage into a misleading fresh timing.

Storage retiming uses a scoped `archive-kallisto-counting` operation so newly
completed local-storage policies remain untouched. Its completion record lists
the relocated library/configuration products. Retimed tasks explicitly require
that matching success record; scheduler ordering alone is insufficient. Failed
or superseded measurements remain separately accounted for in the report.

The launcher records job IDs in `logs/uger_submissions.jsonl` and atomic
`logs/uger_submissions.records/<job-id>.json` files. Accounting collection
merges both sources, deduplicating job IDs. Workers write
per-job/task `.SUCCESS` / `.FAILED` sentinels. The final success record is written
only after all twelve configurations in all four libraries validate and the
tables/plots/report are generated. A failed report removes any older success
record. Inspect `comparison/REPORT.md`, `stage_resources.csv`, `runtime.csv`,
`agreement.csv`, `umi_gains.csv` and `library_summary.csv` before publishing
results. Benchmark products are experimental DGE outputs, not independently
sealed production Step 1b inputs.

If submission stops after dispatching part of the dependency graph, continue
with `--from-mode <next-mode> --hold-job <last-job-id>`, using IDs from the
submission record. Completed data stages are reused; do not submit a duplicate
array over stages that are still running. The launcher strips array-task ranges
from scheduler IDs before constructing job holds.

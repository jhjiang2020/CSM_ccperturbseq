# GEX backends and shared molecule counting — v2.0

For the twelve gene-counting combinations, STAR BAM requirements and
Parse-default/Tony settings, see [Step 1a counting choices](GEX_COUNTING.md).

Start with the [quick start](QUICKSTART.md), [installation guide](INSTALLATION.md)
and [configuration reference](CONFIGURATION.md). Shared lab storage is optional.
Tune CPU/RAM/time through [UGER resources](UGER_RESOURCES.md).

STAR remains the default. The kallisto backend retains Parse 1.6.4 as a
dependency and leaves its installation unchanged. Parse owns preprocessing,
barcode correction, chemistry definitions, sample masks and the cell-calling
algorithm. Repository-owned code handles native transcript pseudoalignment,
gene resolution, UMI families, sparse output and Step 1b combination.

## Configuration

Use a new `paths.output_dir` when changing backends. Old configurations without
`gex_alignment` continue to select STAR.

```yaml
gex_alignment:
  backend: kallisto
  gene_assignment: unique
  umi_correction: parse
  output_level: gene
  count_threads: 0             # 0 uses allocated cores; 1 runs serially
  transcript_assignment: fractional
  kallisto:
    executable: kallisto  # Or an absolute path in your separate tools prefix.
    reference_dir: references/kallisto_GRCh38
    strand: forward                 # forward, reverse, or unstranded
    cell_umi_threshold: null         # null: initialized Parse cell calling
```

`output_level: both` retains gene counts as the primary matrix and writes an
additional transcript/target matrix. `transcript` selects that matrix as the
primary product consumed by Step 5. Gene totals are still computed internally
for the shared cell calls. An explicit positive `cell_umi_threshold` is useful
for controlled fixtures; production normally uses Parse's automatic policy.

Step 1 is submitted through the existing launcher:

```bash
CONFIG="$(realpath my_project.v2.yaml)"
bash submit/submit_pipeline.sh --config "${CONFIG}" --step 1
```

Step 1a runs Parse `pre --one_step`, reformats the resulting survivor FASTQ,
and invokes kallisto `bus`. It resumes a Parse `dge --one_step` entry point
with a repository-owned in-memory callback that obtains initialized reference,
barcode/sample information and cell-calling parameters. It writes the sparse
DGE products directly. Step 1b selects a fraction-safe sparse combine routine.
Neither stage projects pseudoalignments into a BAM.

The existing `gex_decontamination.existing_alignment` BAM-reuse route cannot be
combined with kallisto. Fresh preprocessing supports the usual v1.4 filtering
and Step 1c evidence. BAM/coordinate-specific QC remains specific to STAR.

## Reference and optional tools

Run environment creation and reference building on a compute node. The optional
`environment_kallisto.yml` provisions tools separately from production Python.
Set `executable` to that environment's kallisto binary; pipeline scripts still
run in `splitpipe`.

Build from the same genome FASTA and GTF release used by the Parse reference.
`--genome` must equal its genome label. The default `nac` reference includes
mature transcripts and nascent targets, retaining intronic sensitivity.

See the [reference-build recipe](INSTALLATION.md#optional-kallisto-tools-and-reference)
for a self-contained build or the optional shared tools prefix. The production
job must be able to execute the configured binary: activating a separate tools
environment interactively does not place it on a later production job's PATH.
Both backends still require `paths.genome_dir` to identify the Parse reference.


The builder refuses an existing destination. It uses `kb ref`, checks that the
runtime can read the index, and saves `reference.json`, `index.idx`, `t2g.txt`,
`features.tsv`, category lists, and its build log. The manifest records source
checksums, tool version, command and runtime-asset checksums. The builder maps
kb's versioned gene aliases back to the GTF's exact `gene_id` using its explicit
`gene_version` attribute, recording that policy in the manifest. Transcript and
nascent target IDs remain unchanged. Counting rejects unknown gene IDs or genome
labels rather than silently stripping versions.

## Molecule semantics

`gex_alignment.umi_correction` controls GEX correction for **both STAR and
kallisto**, independently of the gene/transcript output options:

| Value | Neighbor correction |
| --- | --- |
| `parse` (default) | Original Parse 1.6.4 strict-abundance, non-transitive collapse. |
| `directional` | Hamming-one edges satisfying `n_parent >= 2*n_child - 1`, followed transitively using original read counts. Singleton ties can merge. |
| `none` | No neighboring sequences merge; reads sharing an exact UMI within the same counting group still form one molecule. |

Directional roots claim reachable descendants in lexical UMI order, matching
`bustools umicorrect` for gene-resolved evidence. Correction remains scoped to
barcode-index identity and the resolved gene/candidate-gene set. The fractional
gene policy below is therefore still distinct from running the complete
bustools workflow. The gRNA code keeps its existing greedy traversal; sharing
the directional inequality does not make the two traversals identical.

The native unique STAR path uses its original kernel for `parse`. The other modes replace only its
in-memory `collapse_polyN` function through the existing runtime wrapper,
including forked molecule workers. No vendor source changes occur. Use a new
output directory when changing STAR correction policy, since native Parse can
otherwise reuse existing molecule files. Policy manifests prevent mixing
choices at Step 1b or labelling an old combined matrix with a new policy.

Native unique STAR policy records are `process/gex_umi_policy.json` and
`gex_combined/gex_umi_policy.json`. Kallisto records the choice in its existing
completion manifests and H5AD provenance. Step 5 records the common choice as
`ccperturbseq_pipeline.gex_umi_correction`. Historical STAR outputs without a
policy record remain usable with unique genes and the default Parse choice.
Fractional STAR instead uses `star.complete.json` manifests and the shared
sparse counting/combine implementation; see [STAR fractional assignment](GEX_COUNTING.md#star-fractional-assignment).

Preprocessed headers provide corrected barcode sequences, barcode indices,
physical-well identity, RT primer type and UMI. No second barcode correction or
random-primer replacement runs. Each barcode-index identity is counted
separately and only then aggregated into physical wells.

BUS read multiplicities are preserved. Disk-backed SQLite aggregation reunites
records across all lanes/chunks before deduplication. Within a candidate-gene
group, the default UMI families follow Parse 1.6.4's strict abundance/Hamming-distance-one
algorithm, including equal-abundance ties, neighbor order and non-transitive
merges. Each result is checked against the installed Parse collapse kernel.
The configured Parse minimum-read threshold is applied after the selected correction.

With `gene_assignment: unique`, a read's transcript compatibility set must map
to one gene; multiple isoforms of that gene remain acceptable. Multigene reads
are excluded and counted in QC. The final gene matrix contains integer UMI
counts. This reproduces UMI behavior for identical assigned evidence, while
pseudoalignment and STAR can assign different evidence to genes.

With `fractional`, compatible evidence for the same exact UMI is intersected.
Disjoint gene evidence retains independent UMI collisions; evidence compatible
with multiple existing anchors does not create an additional molecule. UMI
correction then operates within the remaining candidate-gene sets. Each
collapsed molecule contributes `1 / number_of_candidate_genes` to each gene.
This adopts the equal-assignment policy used by the Engreitz workflow, with
the custom Parse-compatible UMI policy; it does not claim bit-identical output
to `bustools count --mm`. Counts stay floating point through Step 1b and H5AD.

Transcript output retains compatible targets across the reads in each collapsed
UMI family under the selected gene-assignment policy. Set `gene_assignment:
fractional` to retain multigene evidence in those families. Transcript `unique`
retains singleton target assignments; `fractional` divides
one molecule equally among its compatible targets. This is an EC-derived UMI
matrix, not bulk read abundances, TPM, or an EM estimate. Contradictory target
support remains unresolved and is recorded in QC. Reference targets labelled
`nascent` remain explicitly labelled alongside annotated `transcript` targets.
Gene counts are computed at gene level; they are not obtained by summing
independently deduplicated transcript columns.

Native BUS cannot represent N-containing UMIs losslessly. Such reads are
excluded explicitly and counted in `invalid_umi_reads`. Empty cDNA reads are
likewise counted. No arbitrary base replacement occurs.

## Output contracts

Default primary output paths are unchanged:

```text
gex_<library>/process/{barcode_data.csv,run_proc_def.json}
gex_<library>/process/{tscp_assignment.csv,aw_uf_count_matrix.mtx,aw_uf_bc_counts.csv}
gex_<library>/<sample>/DGE_{filtered,unfiltered}/
gex_combined/<sample>/DGE_{filtered,unfiltered}/
```

Each DGE directory contains `count_matrix.mtx` (cells by features),
`cell_metadata.csv`, `features.csv`, and `anndata.h5ad`. Gene products also have
`all_genes.csv`; transcript products have `all_transcripts.csv`. Gene column
order follows Parse's reference. H5AD gene names follow the pipeline's existing
symbol/genome naming convention; transcript indices are reference target IDs.

`both` additionally writes `DGE_transcript_filtered` and
`DGE_transcript_unfiltered`. All levels use the same called-cell set and order.
Step 1b assigns sample-local `__s1`, `__s2`, ... suffixes in configuration order,
which Step 5 translates to the existing global library prefixes.

Kallisto H5AD `X` and `layers['counts']` contain unnormalized UMI values; no
clustering or expression transformation is performed. `tscp_count` and
`gene_count` describe the internal gene matrix; `matrix_count` and
`matrix_feature_count` describe the selected matrix. Fractional outputs are
identified in `uns['gex_quantification']`; consumers needing integer counts
must select the default unique-gene product. Gene-specific downstream analysis
requires a gene matrix.

Each gene-resolved molecule contributes one total UMI, including fractional
assignment. Cell calling uses whole-molecule totals, avoiding floating-point
roundoff at exact cutoffs; individual gene weights remain fractional.

`tscp_assignment.csv` retains Parse's molecule-table columns and adds
`assignment_weight`. Its `count` is supporting reads, not expression or the
fractional molecule contribution. Precise alignment-derived `exonic` values
are unavailable from native ECs and left missing. `kallisto_qc.json` and native
alignment logs replace STAR-specific alignment/report metrics for this branch.

Completion manifests bind configuration, reference, preprocessing artifacts,
decontamination evidence and output hashes. Step 1b and Step 5 reject stale or
incomplete products. Backend mixing is rejected before preprocessing. A failed
run never installs a successful completion manifest.

## Validation

Run from the project root in the production Python environment:

```bash
python -m unittest discover -s tests -p test_kallisto_backend.py -v
CSM_TEST_KALLISTO=/medpop/rgupta/Jiahao/conda/envs/csm_kallisto/bin/kallisto \
  python -m unittest discover -s tests -p test_kallisto_backend.py -v
```

The suite exercises Parse's original collapse function as an oracle, both
neighbor-search branches, UMI collisions, isoform ambiguity, primer-index
aggregation, fractional conservation, malformed BUS records, reference hashes,
matrix/H5AD round trips, combined barcode suffixes and stale-output rejection.
The optional executable test builds a tiny real index and pseudoaligns reads.
Set `CSM_TEST_BUSTOOLS` to the optional environment's `bin/bustools` to check directional families against
the real `bustools umicorrect` executable, including chains and competing roots.

Before production use, benchmark matched preprocessed inputs on the four
sublibraries, measuring molecule/gene yield, cell-call overlap, guide joins,
runtime, memory and temporary disk use. These small-fixture checks do not
establish full-depth resource requirements or biological equivalence to STAR.

The pilot benchmark uses 10 million raw paired reads per C1/C2/T1/T2 library,
then shares the exact preprocessed survivor FASTQ across all combinations.
`scripts/benchmark_gex_backends.py` and `submit/submit_kallisto_benchmark.sh`
cover nine combinations: STAR, kallisto unique genes, and kallisto fractional
genes, each with `parse`, `directional`, and `none`. Array tasks 1–12 cover the
Parse-default baseline, 13–24 directional correction, and 25–36 no correction;
each group contains four libraries per backend/assignment choice.

Each core runtime independently measures alignment through DGE, using 16 CPU
slots, with preprocessing added separately. Reference construction, sampling,
Step 1b validation and optional clustering/reporting are excluded from that
timing. Existing prepared reads can be used to add correction choices without
sampling or filtering again. The report compares matching correction modes
between backends and also compares correction modes within each backend. It
retains fractional counts and uses a fixed reference-called cohort as well as
each arm's own called cells.

The completed [pilot report](KALLISTO_PILOT_BENCHMARK_REPORT.md) includes all
36 library runs, UMI/cell yield, fixed-cohort and expression agreement, cell-call
overlap, component timings and figures. The validated mature-plus-nascent
reference built from the Parse GRCh38-2024-A FASTA/GTF is available at
`output/kallisto_benchmark_10m_20261022/reference`; its `reference.json` records
the source and runtime checksums. This is a full reference, independent of the
10-million-pair read sampling.

The shared STAR-fractional/kallisto molecule counter uses disk-backed, complete-barcode shards across
the allocated threads. `gex_alignment.count_threads` can cap the worker count
(0 uses all allocated threads; 1 runs serially). The legacy kallisto key is accepted. The counter caches EC/gene annotations
and merges shards in original barcode order, preserving molecule assignments
and fractional-weight accumulation exactly. Each Parse-corrected UMI group
retains the installed Parse kernel parity check. Threads are reused sequentially
for pseudoalignment and counting, rather than nested. Parallel counting requires
`fork`; use one counting thread on platforms without it. Temporary shards live
in the job temporary directory and are removed after success or worker failure.
`kallisto_qc.json` records workers, shards and reference/spool/correction/merge/
matrix timings under `counting_execution`. The benchmark reports
formatting, alignment, counting and publishing separately for kallisto. Its
mature-plus-nascent reference also retains intronic targets. Observed core wall
times describe this full implementation and reference, rather than an intrinsic
speed advantage of one aligner.

The optimized counter is checked against saved outputs by
`scripts/benchmark_kallisto_optimization.py`. Run it on a compute node after
all benchmark arms complete:

```bash
python scripts/benchmark_kallisto_optimization.py --config config.yaml \
  --baseline output/kallisto_benchmark_10m_20261022 \
  --work-dir output/kallisto_parallel_benchmark_10m_20261022
```

The comparison requires exact matrix entries (including fractional weights),
ordered cells/features, cell metadata, byte-identical molecule assignments,
scientific QC and H5AD raw counts. It also runs and checks the Step 1b handoff.
Runtime tables distinguish counting from the complete alignment-through-DGE
stage. The baseline report and outputs remain separate from the new results.
For the earlier benchmark's documented low-UMI metadata-repair rows, the
validator also reconstructs the original serial sums from saved assignments
using round-trip float parsing. It records any prior CSV serialization noise
separately and requires exact equality against that serial reference; baseline
files are never rewritten and no tolerance hides count differences.

Reference implementation:
[Engreitz GEX workflow](https://github.com/EngreitzLab/telohaec-genomewide-perturb-seq-analysis/blob/280a3747288b8cb30e5102a10c9777e391cf7756/workflows/00_fastq_to_sublibrary_h5ad.smk).
Directional traversal reference:
[bustools 0.45.1 implementation](https://github.com/BUStools/bustools/blob/v0.45.1/src/bustools_umicorrect.cpp).

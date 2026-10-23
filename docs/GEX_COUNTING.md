# Step 1a counting choices

Step 1a supports all twelve gene-counting combinations:
`(star | kallisto) × (unique | fractional) × (parse | directional | none)`.
Parse split-pipe 1.6.4 supplies preprocessing, corrected barcode identities,
reference annotation, sample masks and automatic cell calling. Its installation
is never modified. GEX UMI choices are independent of Step 4 guide correction.

```yaml
gex_alignment:
  backend: star                 # star or kallisto
  gene_assignment: unique       # unique or fractional
  umi_correction: parse         # parse, directional, none
  output_level: gene            # STAR: gene; kallisto: gene, both, transcript
  transcript_assignment: fractional  # kallisto target output only
  count_threads: 0              # shared counter: allocated cores; 1 = serial
  star:
    max_multimaps: 3            # STAR reports all alignments up to this limit
  kallisto:
    executable: kallisto
    reference_dir: references/kallisto_GRCh38
    strand: forward
```

Only configure the selected backend's tools/reference. Historical kallisto
configs can keep `gene_assignment`, `output_level`, `transcript_assignment` and
`count_threads` under `gex_alignment.kallisto`. New configs put these policies
directly under `gex_alignment`; contradictory values at both locations are an
error. Omitting the section retains STAR + unique + Parse correction. Use a
separate `paths.output_dir` for a different scientific policy.

## Common settings and the reference workflows

| Workflow | Backend | Gene assignment | Neighbor UMI correction | Other relevant choices |
| --- | --- | --- | --- | --- |
| Parse 1.6.4 default counting | STAR | `unique` | `parse` | MAPQ ≥255, stranded gene-span assignment; Parse ranks overlapping genes and picks one. |
| STAR with ambiguous molecules retained | STAR | `fractional` | Any of the three | All retained loci and all overlapping gene candidates; default maximum 3 loci per read. |
| Kallisto with Parse-style correction | kallisto | `unique` or `fractional` | `parse` | BUS compatibility sets, Parse preprocessing/cell calls. |
| Closest counting settings to Tony's published workflow | kallisto | `fractional` | `none` | Published config is **unstranded**, mature + nascent `nac`, `--sum total`, `--mm`. |
| Optional directional comparison | Either | Either | `directional` | Hamming-one graph with `n_parent >= 2*n_child - 1`; transitive traversal. |

Parse's default row describes its counting policy. This pipeline's optional
guide decontamination and kit-score settings are separate experiment choices.
In the unique STAR path, the native Parse molecule/DGE implementation remains
in use, with an in-memory UMI-kernel replacement for `directional` or `none`.
Consequently, `unique` means the legacy unique-alignment/selected-gene behavior
on STAR, whereas kallisto requires its read compatibility set to resolve to
one gene. Multiple isoforms of that gene are allowed.

Tony's source is pinned to
[`280a3747288b8cb30e5102a10c9777e391cf7756`](https://github.com/EngreitzLab/telohaec-genomewide-perturb-seq-analysis/tree/280a3747288b8cb30e5102a10c9777e391cf7756).
The [GEX rule](https://github.com/EngreitzLab/telohaec-genomewide-perturb-seq-analysis/blob/280a3747288b8cb30e5102a10c9777e391cf7756/workflows/00_fastq_to_sublibrary_h5ad.smk)
calls [its wrapper](https://github.com/EngreitzLab/telohaec-genomewide-perturb-seq-analysis/blob/280a3747288b8cb30e5102a10c9777e391cf7756/scripts/run_kb_count.py)
with `nac` and `--mm`; the wrapper adds `--sum total`. The
[configuration](https://github.com/EngreitzLab/telohaec-genomewide-perturb-seq-analysis/blob/280a3747288b8cb30e5102a10c9777e391cf7756/config/config.gw_combined.yaml)
sets `strand: unstranded`. The corresponding
[`kb-python` 0.29.5 counting code](https://github.com/pachterlab/kb_python/blob/v0.29.5/kb_python/count.py)
uses gene-level exact-UMI deduplication and does not call `bustools umicorrect`.
Thus `none` is the closest neighbor-correction setting; exact duplicates still
deduplicate. `directional` is an additional experiment, not a Tony default.

For that comparison, use:

```yaml
gex_alignment:
  backend: kallisto
  gene_assignment: fractional
  umi_correction: none
  output_level: gene
  kallisto:
    executable: kallisto
    reference_dir: references/kallisto_GRCh38
    strand: unstranded
```

These settings match the broad counting choices, not the entire workflow.
This pipeline retains Parse barcode preprocessing and cell calling, its own
candidate-set/collision resolution, and the configured reference. Tony's
published workflow uses its own barcode replacement, IGVF/GENCODE 43 reference,
nascent/mature handling and BarcodeRanks inflection cell calling. Reproducing
that workflow exactly requires its code and reference assets. Neither `--mm`
nor this pipeline's fractional mode estimates expression by EM or produces TPM.

## STAR fractional assignment

A genomic STAR BAM can support fractional gene assignment. The quantifier
must retain alternative alignments and annotate genes; STAR alignment alone
does not define fractional UMI counts. See the
[STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
for multimapping output settings.

Fresh fractional STAR processing runs Parse `pre --one_step`, `align
--one_step` and `post --one_step`, then the shared counting callback through
Parse's initialized DGE entry point. It bypasses native molinfo, whose BAM
splitter drops secondary alignments and whose MAPQ gate admits unique mappings.
The STAR command explicitly sets `--outFilterMultimapNmax` and
`--outSAMmultNmax -1`. Reads above the mapping limit remain excluded by STAR;
fractional counting cannot recover alignments STAR did not report.

The adapter processes one original read at a time from the queryname-sorted
`process/barcode_headAligned_sorted.bam`. It checks one primary record and an
`NH` value equal to the number of retained mapped records. It requires Parse's
original read ID (`pre_head_keep_original=True`, the default). Missing hits,
duplicate original IDs, paired-end BAM records, supplementary records or a
wrong sort order fail explicitly. Low STAR multimapper MAPQ is expected and is
not filtered with Parse's unique-only MAPQ 255 threshold.

For each alignment, candidate genes use Parse's reference bins, forward strand
and whole aligned-read containment in the gene span, including introns.
All overlapping genes survive; Parse's coding/name/exon preference is not used.
Candidate genes are unioned across the read's loci, with duplicate gene IDs
removed. Thus two loci in the same gene still support one gene. Intergenic loci
do not become output features. This is gene-set weighting, not `1/NH` weighting.

STAR and kallisto then share the following processing:

1. Spool complete read evidence into SQLite, keeping corrected barcode-index
   identities distinct until counting is finished.
2. Resolve compatible evidence within each exact UMI by candidate intersection;
   retain disjoint gene collisions as separate molecules. Ambiguous evidence
   supporting multiple anchored molecules does not add another molecule.
3. Apply the selected UMI correction within each remaining candidate-gene set.
   `parse` uses the vendor collapse kernel as an oracle; `none` retains exact
   deduplication; `directional` traverses the fixed-count graph.
4. Apply Parse's minimum supporting-read threshold after correction. Each
   surviving family contributes `1/k` to its `k` candidate genes. An unresolved
   A/B molecule contributes 0.5 to each; two alignments to A contribute 1 to A.
5. Aggregate physical wells, call cells with whole-molecule totals and write
   floating sparse counts. The shared Step 1b combine preserves those weights.

N-containing UMIs are excluded with an explicit counter, matching the native
BUS representation. Fractional STAR reports alignment/molecule diagnostics in
`process/star_qc.json`; per-family exonic fractions are left missing because the
shared evidence representation does not retain a resolved genomic origin.
The sorted genomic BAM remains available for coordinate-based QC.

The molecule table keeps Parse's historical `tscp_assignment.csv` filename.
Its `count` column is supporting reads, and `assignment_weight` is the gene's
fraction of a molecule. “tscp” in this filename does not imply isoform counts.
Ordinary Parse reports/clustering are not generated by the shared callback.

## Transcript output and provenance

Transcript `unique`/`fractional` output is supported by kallisto's target
compatibility sets. STAR's current Parse reference adapter annotates genes;
it does not reconstruct isoform compatibility from a genomic BAM. STAR with
`output_level: transcript` or `both` fails configuration validation. Transcript
quantification from STAR alignments is possible with a separate transcript-aware
quantifier, but is outside this adapter. See [target semantics](KALLISTO_BACKEND.md).

Fractional STAR has `process/star.complete.json` and
`gex_combined/star.complete.json`; kallisto retains its `kallisto.complete.json`
names. An early `process/gex_counting.settings.json` also binds partial shared-counter
runs before Parse can reuse intermediates. These contracts bind settings, preprocessing, annotation reference, outputs, and
for STAR the source BAM. Step 1b/5 reject changed inputs or policies. Native
unique STAR retains `gex_umi_policy.json`; legacy unmarked products are accepted
only for unique genes with Parse correction.

Historical BAM decontamination can feed fractional STAR only from a complete
pre-molinfo sorted BAM. The annotated BAM is not selected as a fallback because
Parse may already have removed secondary hits. The source remains read-only;
the existing sandbox route filters complete QNAME groups before counting.

Run the focused checks in a compute allocation:

```bash
CSM_TEST_STAR="$(command -v STAR)" \
  python -m unittest discover -s tests -p test_star_fractional.py -v
python -m unittest discover -s tests -p test_kallisto_backend.py -v
```

Fixtures cover twelve gene-policy configurations, overlapping/spliced genes,
multilocus ambiguity, same-gene duplicate hits, incomplete BAM rejection,
all UMI modes, read thresholds, parallel equality, fractional MTX/H5AD output,
combination, stale input rejection and cleanup after reference-validation failures.
The optional real STAR test builds a tiny repeated-locus genome, checks the
reported secondary alignments and fractional counts under all UMI modes, and
rejects a primary-only BAM. Full-depth STAR fractional resource and
biological comparisons require new measurements. The current
[50M-pair benchmark](GEX_BENCHMARK_V2.md) covers all twelve configurations with
shared alignment/assignment caches, stage resource accounting, all-pairs
agreement and fixed-cohort UMI gains. See the completed
[measured comparison and assignment-cache explanation](GEX_BENCHMARK_V2_REPORT.md).
The earlier 10M-pair kallisto reports are superseded by this experiment.

# Four-library guide assignment benchmark

Implemented and launched on 2026-11-16. This is an experimental benchmark;
production assignment defaults and expression products are unchanged. The
[analysis protocol](../analysis/docs/grna_assignment_benchmark.md) defines the
KD comparisons, missing-partner experiment, interpretation and limitations.

## Frozen inputs and processing

The cohort is the 37,086 GEX-good cells in
`output/pilot_4_sublib_audit/work/cohort_manifest.json`: C1 7,187; C2 6,994;
T1 11,284; T2 11,621. Cells with no guide counts remain in the denominator.
The expression source is `output/combined_adata.h5ad`, using the existing Parse
preprocessing, STAR and unique gene assignments. Neither alignment nor GEX
molecule/gene counting is resubmitted. No production Step 5 output is replaced.

`scripts/grna_assignment_benchmark.py prepare` snapshots a task-specific
reference and configurations. The reference has 500 constructs, 52 distinct
spacer sequences and 68 spacer/position/target identities. Six control-guide
labels each name two different sequences. Canonical benchmark guide names
therefore identify the **designed spacer, position and target**, retaining the
original labels as aliases. Counts retain iBAR-specific `guide_unit_id` and
construct candidates before aggregation. The same sequence in opposite
positions remains distinguishable; the reference permits subsequent pooling.

Explicit unpositioned `unit_*` evidence cannot identify a spacer/position guide.
It is retained in `matrices/unpositioned_units_<arm>.csv.gz` and summarized in
the recovery table, but is not distributed between canonical guides. Unexpected
identities with a known position remain fatal validation errors.

Both new standard and recovered runs use the current design-aware matcher and
the same canonical reference. Historical production pair calls are a separate
comparator, not a claim of exact historical read-processing reproduction.

The run performs:

1. Step 2 on full dedicated-arm FASTQs, enforcing the original Step 1a barcode
   contract and BC1 namespace.
2. Standard Step 3 and a separate recovered Step 3. Recovered matching selects
   one shared `relaxed_structure: auto` policy across all four libraries, with
   maximum 5-prime loss 3 nt, retained-length floor 17 nt, endpoint shift at
   most 3 nt and retained-base Q30. General matching edit budgets are unchanged.
3. Full GEX FASTQ Parse preprocessing in disposable scratch, with the current
   decontamination hook writing durable barcode-positive guide evidence. The
   captured barcode coordinates are checked against production coordinates.
   This is **not** the selected GEX archive used in the earlier pilot refresh.
   Preprocessing does not advance to STAR or expression counting.
4. Step 1c with that shared matching policy, its strict guide/iBAR filters and
   completion-manifest validation; then GEX-only aggregation.
5. Standard dedicated, recovered dedicated and recovered combined count
   products. The combined product merges read-level evidence **before**
   minimum-read filtering and directional UMI collapse. Collapsed matrices are
   never added. Source provenance remains in the production count ledger.
6. Guide-level matrices before assignment thresholds, all assignment methods,
   and cell-, guide-, target- and construct-level concordance summaries.

The decontamination filter still requires reference-supported full-guide
evidence. Relaxed Step 1c matching cannot recover GEX reads that never entered
that filter's evidence stream. The full FASTQ denominator does not imply an
exhaustive search for every shortened GEX guide.

## Methods and comparisons

The pinned package is crispat 0.9.11, commit
`a9372b9cf9b82b7409cc3be4d192056ae57a2673`. All eleven families are run:

| Family | Implemented benchmark settings |
| --- | --- |
| `ga_umi` | Guide UMI thresholds 1, 2, 3, 5, 8, 10, 15, 20 |
| `ga_max` | Native whole-cell maximum and separate per-position maximum |
| `ga_ratio` | Native and per-position; fractions >0.1, 0.2, 0.3, 0.5, 0.7, 0.9 |
| `ga_gauss` | Native VI with all counts and nonzero-only fitting |
| `ga_poisson_gauss` | Native nonzero log-count mixture |
| `ga_2beta`, `ga_3beta` | Native and per-position dominant-guide mixtures |
| `ga_poisson`, `ga_negative_binomial`, `ga_binomial` | Native batch-aware count mixtures |
| `ga_quantiles` | Retain top 50%, 70% or 90% of nonzero guide fractions per guide |

All non-UMI methods have an explicit minimum of one observed UMI. Count
mixtures pool libraries with a batch term; the other families run separately
by library. Zero-count rows stay in the shared output universe. Native package
iteration counts and initialization behavior are retained; no expression
outcome selects an assignment model. Fits are transductive on guide counts.
Fitted parameters, thresholds and package outputs remain under `fits/`.
Unavailable/failed fits are recorded, never converted to negative calls.

Maximum, ratio and beta methods primarily choose a dominant guide. Native
whole-cell use structurally misses a true second guide; per-position adaptation
still cannot represent arbitrary multi-construct cells. Both are useful
diagnostic comparators, not unrestricted multi-guide default candidates.
Per-position adaptation is explicitly labeled as a CSM wrapper, not an
unchanged crispat algorithm.

Additional comparators are:

- Current pair-sum threshold 5 on historical, new standard, recovered dedicated
  and recovered combined evidence. A 3+3 pair passes this rule; neither guide
  passes a per-guide threshold of 5. Only observed arms enter observed matrices.
- Combined pair-sum threshold 5 with the designed missing partner imputed.
  This explicit sensitivity matrix can contain positive labels with zero UMIs.
- Independent guide threshold 5 plus iBAR-supported rescue of an observed
  partner with at least 1 or 2 UMIs, when the opposite arm has at least 5 UMIs
  **in the same resolved, structurally concordant construct**. Zero-UMI arms
  are never rescued by these rules.

Concordance uses positive cell-guide Jaccard, directed disagreements, exact
cell guide-set equality, per-guide and per-target agreement, 0/1/2/>2 assigned
guide counts and complete/half-observed resolved constructs. Shared absences
are not the primary agreement statistic. Agreement and pair compatibility are
not experimental accuracy estimates, particularly when iBAR also informed a
method's calls.

## Implementation and execution

Processing and assignments belong to the parent repository:

- `scripts/grna_assignment_benchmark.py`: frozen cohort/configuration and full
  GEX evidence capture.
- `scripts/benchmark_guide_methods.py`: matrix export, crispat adapters, iBAR
  comparators and concordance.
- `submit/submit_assignment_benchmark.sh`: versioned UGER submissions.
- `tests/test_assignment_benchmark.py`: identity, zero preservation,
  pair-vs-guide thresholds, imputation, tie handling and all eleven APIs.

The task core is `output/grna_assignment_benchmark_20261116/`. Each submission
has a timestamped `jobs/` directory and an entry in `submissions.tsv`. Scientific
work runs on compute nodes; dependencies check completion artifacts as well as
scheduler ordering. The assignment venv uses the read-only `csm_analysis` base
with pinned crispat/Pyro/Dask additions and a recorded resolved package inventory:
`/medpop/rgupta/Jiahao/venvs/csm_assign_20261116_v2`.

Stages, in dependency order, are `prepare`; `capture` and `barcodes`; `pilot`
and `standard_match`; `match`; `step1c`; `count`; `export`; `methods`; and
`concordance`. `environment` and `tests` gate method execution. The launch
records contain the actual job IDs; the downstream submodule owns its own
submission graph and output core.

Do not run this benchmark's configurations against `output/` directly or
replace production reference names. Completed products refuse overwriting.
Retain failed-attempt logs and reuse completed inputs only after checking their
manifests. Large scratch Parse outputs are reproducible intermediates, not the
only copy of evidence; no cleanup was authorized as part of this launch.

Sources: [crispat source](https://github.com/velten-group/crispat/tree/a9372b9cf9b82b7409cc3be4d192056ae57a2673),
[Braunger and Velten, 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11399228/).

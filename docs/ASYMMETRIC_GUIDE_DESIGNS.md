# v2.0 optional asymmetric guide designs

The optional design-aware gRNA arm supports symmetric and asymmetric dual
constructs, including 7/7-, 7/1-, 7/0-, 0/7- and 0/0-nt iBAR layouts. Spacer
lengths are derived from the reference, including 19/20-nt spacers. Existing
configurations retain the legacy parser, matching schema and assignment path.
The GEX UMI correction policy is independent of these guide settings.

## Enable the feature

Start with a complete pipeline configuration and add:

```yaml
guide_design:
  enabled: true
  strict: false
  ibar_orientation: as_reference
  discovery_max_ibar_length: 32
  discovery_grna_lengths: [18, 19, 20, 21, 22, 23]
  scaffold_by_position:
    pos1: A
    pos2: B

grna_structure:
  # Retain your fixed_seq_1 through fixed_seq_5 anchor definitions.
  grna_lengths: [19, 20]
  ibar_lengths_by_position: auto

guide_assignment:
  dual_guide_counting_mode: true
  infer_missing_partner: auto
```

Merge these entries into existing blocks; do not replace the remaining
configuration. `ibar_lengths_by_position` also enables the feature and accepts
an explicit mapping such as `{pos1: [7], pos2: [0]}`. Length zero means an
intentionally absent iBAR. Distinct lengths within a position must occur in
the reference; they are not an unconstrained production search range. The
legacy scalar `ibar_length` remains valid for symmetric configurations.

Step 0 treats the design as authoritative and warns when explicit lengths
disagree. It records requested-setting conflicts and the effective reference
lengths in `<reference>.design.yaml`. `strict: true` turns those warnings into
errors, including requests for unsafe iBAR correction. Set an appropriate
`grna_matching.max_lev_dist_ibar` when using strict mode. Missing or invalid
anchors, invalid DNA and incomplete dual-pair annotations are errors.

Step 0 preserves the existing reference columns and adds observable guide-unit
IDs, iBAR lengths and design absence annotations. It also writes capabilities
for position identification, single-arm pair resolution, joint iBAR/spacer/full
signature resolution, and missing-partner inference coverage. Two distinct
construct IDs with identical observable sequences remain ambiguous.

## Run the bounded pilot before production

Run these commands inside a compute allocation with the production environment
active. Step 0 and the pilot open project data; they are not login-node probes.

```bash
python scripts/00_format_reference.py --config config_asymmetric.yaml
python scripts/03_match_grna_ibar.py --config config_asymmetric.yaml \
  --preflight-only --preflight-fastq
```

The pilot measures iBAR boundaries by independently anchoring the stem and
downstream tracr across a bounded discovery range. It can therefore flag lengths
that differ from the design instead of depending on the configured length to
extract them. Spacer discovery defaults to 18–23 nt plus the reference lengths;
it reports unexpected lengths and low exact-reference support. It also checks
exact spacers against scaffold positions and
compares iBAR orientations conditional on spacer and position.

The sample is bounded by `grna_matching.preflight.sample_reads` (default 10,000,
maximum 100,000), apportioned across input files. Alternating reads form
discovery and held-out validation subsets. Both subsets need sufficient,
consistent strand evidence before an orientation is resolved. Per-file evidence
is reported and conflicting file orientations prevent automatic reversal.
These are read-weighted diagnostics; PCR copies are not independent molecules.
Short, sparse or strand-neutral evidence remains inconclusive.

`ibar_orientation: as_reference` diagnoses orientation without changing the
reference. `ibar_orientation: auto` permits a validated per-position reversal
in a new effective contract. Source CSVs remain unchanged. No-iBAR positions
have orientation `not_applicable`; a 1-nt iBAR can be strand-inconclusive.
Wrong lengths or construct annotations are diagnosed, not silently rewritten.
Only configured scaffolds and the supplied reference can be evaluated; the
pilot cannot identify an arbitrary unknown library.

Each library writes:

- `grna_matches/<library>_preflight.yaml`: evidence, warnings and conclusions.
- `grna_matches/<library>_design.yaml`: effective sequences, capabilities and
  fingerprints used for matching.

To freeze a reviewed pilot result for all libraries and evidence sources, set:

```yaml
guide_design:
  enabled: true
  ibar_orientation: auto
  resolved_contract: output/grna_matches/C1_design.yaml
```

The contract path is relative to the configuration file unless absolute. Its
content and input fingerprint must match the reference/settings. A subsequent
pilot with confident contradictory orientation fails. This is particularly
useful before Step 1c, whose small GEX-derived evidence may not resolve strand
independently. Without a frozen contract, libraries resolve separately and
Step 4 refuses incompatible orientations rather than merging their counts.

## Extraction, matching and molecule counts

The production extractor scores candidate spacer, stem, iBAR and downstream
boundaries jointly. It records length-based position and scaffold-based position
separately. Length can identify position when allowed length sets are disjoint;
identical scaffold sequences do not provide independent position evidence.
Scaffold/length conflicts remain diagnostic guide observations with no assigned
pair. Ambiguous boundaries are retained in failed-read output when
`save_failed_reads: true`. A truncated read cannot establish iBAR absence.

The design-aware extractor requires a complete configured downstream anchor,
including at positions with no iBAR. Existing single-guide short-read parsing
remains available through the legacy configuration.

iBAR correction is conditional on spacer and position. Each codeword's allowed
edit radius is bounded by its distance to other codewords in that group. Zero-
and one-base iBARs use exact matching. Very large codebooks (over 1,024 iBARs
for one spacer/position) also use exact matching to bound validation cost.
An exact foreign iBAR remains a decoupling signal even if a fuzzy within-guide
correction could otherwise absorb it.

UMI correction remains available without iBARs. It operates within observable
guide-unit groups, separately from iBAR correction. Design-aware grouping
preserves iBAR identity even if the legacy `umi_dedup_ibar` setting is false;
Step 0 warns about that override. Shared evidence is never copied into multiple
candidate constructs, including when dedicated and GEX evidence are merged.
An identical UMI observed in the two positions does not establish physical
linkage between their transcripts.

Step 4 can resolve an ambiguous arm through compatibility with the other
detected arm. When several complete pairs compete for the same observation,
the observation remains unallocated. Missing partners are inferred only from
units uniquely identifying one designed pair, and can be disabled using
`infer_missing_partner: false`. Inferred guides never receive observed UMIs.

Read-level extensions include `guide_unit_id`, `position`, `ibar_length`,
`ibar_status`, `pair_candidates` (JSON), `pair_candidate_count`,
`resolution_status`, `position_from_length`, `position_from_scaffold`,
`scaffold_variant`, `tracr_variant`, `structure_status`,
`pair_resolution_basis`, `guide_design_fingerprint` and `legacy_model_eligible`.

`ibar_status` distinguishes `present_resolved`, `present_unresolved` and
`absent_by_design`. Read truncation is recorded as `unobserved_truncated` in
failed-read reasons when the minimum structure cannot fit. Cell assignments
also distinguish an unobserved arm and retain `unpositioned_umi` for detected
guides whose position cannot be resolved. Pair IDs are empty for unresolved calls;
candidate lists preserve their ambiguity. Assignment statuses add `ambiguous`,
`guide_only` and `discordant` alongside `complete` and `inferred`.

Step 4 validates the reference/settings fingerprints before consuming matches.
Use separate output roots when switching legacy/design-aware modes. Step 5
retains candidate and iBAR-state metadata, adds `has_resolved_pair`, and stores
per-library contract summaries as JSON in `adata.uns['guide_design']`.
`has_grna` remains true for detected guide evidence even without a resolved pair.

## QC availability

Q1 uses the same structural extractor with its sequencing adapters and analyzes
physically linked R1/R2 observations. Enable `guide_design` in the standalone
Q1 configuration or inherit it through `sc_config`. The maintained config
loader passes these settings into the `guideseq` block. Direct Python callers
pass `guide_design: {enabled: true}` within their `guideseq_cfg` mapping.

The design-aware Q1 path writes linked observations, per-library rates,
capabilities and frequency tables for uniquely resolved designed constructs.
`recombination_compatible_rate` is the observed incompatible fraction among
eligible linked units. Spacer information can make this measurable even when
the iBAR tuple is uninformative. Events producing another valid design signature
are invisible; sequence errors can also produce incompatible observations.
The rate is not a total biological recombination estimate. Design-compatible
observations include ambiguous construct IDs and are reported separately from
uniquely resolved pairs. Per-position decoupling rates use only known-iBAR
observations from spacer groups with distinguishable forbidden combinations;
their eligible denominators and unavailable reasons are explicit. The new Q1
path reports compatibility rather than the historical fine-grained event
decomposition, which remains available through legacy symmetric Q1 mode.
Right-primer anchoring on R2 verifies the boundary; its scaffold label is
`unobserved` for the tracr rather than claiming an independent tracr observation.

Q2 writes `q2_design_capabilities.yaml` and a design-aware UMI threshold sweep
when the legacy observation model is unsupported. The legacy likelihood and
Q3 inferred-source heuristic remain eligible only for symmetric, resolving
iBAR designs with distinct scaffolds and iBARs of at least two bases. The
capability is propagated through counting and assignment; it is not inferred
from an empty barcode column at QC time. Unavailable estimates are null with a
reason. Q3 can still analyze resolved efficacy cohorts when its inferred
dropout-versus-recombinant source classifier is unavailable.
Q4 excludes unresolved construct assignments from fitness comparisons.
Published experiment-specific pilot reports retain their original definitions;
they are not generic asymmetric-design QC reports.

## Validation

```bash
python -B -m unittest discover -s tests -p 'test_guide_design.py'
python -B -m unittest discover -s tests -p 'test_step3_matching.py'
python -B -m unittest discover -s tests -p 'test_grna_preflight.py'
python -B -m unittest discover -s tests -p 'test_step1c*.py'
```

Tests cover symmetric/asymmetric/no-iBAR extraction, shared signatures, pair
ambiguity, scaffold conflicts, truncation, conservative correction, strand and
length pilots, linked incompatibility QC, source deduplication, Step 0 output,
parallel Step 3 matching and Step 4 CLI/fingerprint validation.

The [UGER validation report](ASYMMETRIC_GUIDE_UGER_VALIDATION.md) records the
production `splitpipe` smoke test using a shipped pooled pilot-derived 7/4-nt
fixture, an unedited 7/7-nt control and the regression suite. It includes
reproduction instructions and explicitly separates artificial cell/UMI truth
from real Parse barcode extraction.

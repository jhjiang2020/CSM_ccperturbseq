# v2.1 plan: recover guide evidence without forcing identity

Status: implementation and four-library benchmarks completed, 4 November 2026. Automatic recovery adds 4,931,644 reads, 437 assigned cells and 694 missing partners, with a net gain of 127 singlets. The primary biological gate is inconclusive, so recovery remains optional and the recommended starting preset is audit-only. The consolidated [v2.1 guide](GRNA_RECOVERY_V2_1_EVIDENCE.md) records activation, parameters, results and limits.

Implementation scope: Steps 2–5 use the shared recovery policy; Step 1c reuses the dedicated library's frozen policy. Q1 amplicon QC keeps its separate primer-based parser. Partial TSO, scaffold indels and larger TSO edit limits remain deferred experimental extensions. Asymmetric/no-iBAR compatibility is tested with fixtures; the pilot supplies biological data only for its own layout. This document preserves the original design and release criteria below, including criteria that the pilot could not establish.

## Recommended direction

Make joint TSO–spacer–scaffold interpretation the main change. Add a separate, optional rule for 1–3 bases missing from the spacer's **own designed 5′ end**. Retain spacer-first identification, with optional iBAR corroboration. Select conservative settings from the reference and a bounded pilot, then freeze them before production.

Do not start by raising the general guide edit distance or TSO edit limit. The earlier endpoint audit and the new four-library probe recover reads while retaining the current TSO limit of three edits. A retained length of 18 bases is a useful conservative comparison; 17 bases is a candidate automatic setting when reference separation and held-out tests support it.

## 1. Extraction and identity rules

### Interpret boundaries together

1. Find candidate stem and downstream scaffold locations using the reference's layouts. Enumerate the upstream spacer lengths and TSO endpoints compatible with those locations.
2. Score the TSO at each candidate endpoint, including alternatives that are worse than the best isolated TSO alignment but remain within the allowed edit limit. Considering only equally optimal alignments would miss much of the earlier rescue.
3. Initially retain full-TSO alignment, the existing search region and the existing combined scaffold mismatch limit. Test scaffold indels, more than three TSO edits and partial TSO in separate experimental arms.
4. Preserve all credible guide identities and length interpretations. Resolve an identity only when the evidence supports it; record uncertainty rather than choosing the first layout.

A read can identify one guide while supporting both full-length and shortened layouts. Count it once as guide evidence, label its length uncertain, and exclude it from claims about abnormal-only activity. Do not automatically prefer a full-length interpretation simply because it is longer.

Start assignment-changing rescue with reads that failed extraction or matching. Preserve established calls and report conflicting alternative identities separately. Revisit already assigned reads in an explicit sensitivity analysis, not through silent reassignment.

### Make shortening reference-aware

For a guide of length `L`, test its exact suffixes after loss of 1, 2 or 3 bases. A 20-base design can therefore contribute a 17-base candidate; a 21-base design has an 18-base minimum under the same loss limit. A global `minimum_design_length − 3` cutoff would also admit four-base loss from the latter.

Build a suffix index over the **complete reference**, including intact sequences and competitors from other permitted loss lengths. Include a broader set of plausible short competitors for ambiguity checks without necessarily admitting them to assignment. Never remove rare or undetected reference guides to make a read appear unique.

Initially require exact retained-spacer sequence. Compute nearest-neighbor edit distances between the observable sequences, not just between full spacers. Any later fuzzy correction must have a unique candidate and a correction radius supported by the reference separation and held-out error tests. Where position is unresolved, compare candidates across positions.

The automatic search should test maximum loss 0–3 and retained-length floors of 18 and 17 bases. A shorter floor can be explicitly evaluated for other references, but is not a default justified by this pilot. Intact guides shorter than those rescue floors remain valid according to their own design.

### Support all existing construct layouts

Use the same design contract for single and dual guides, symmetric and asymmetric iBARs, and absent iBARs at either or both positions. Required fixtures include 7/7, 7/1, 7/0, 0/7 and 0/0, mixed reference spacer lengths, shared spacers and identical scaffolds.

Preserve the existing single-guide route for reads truncated at the 3′ end as a distinct capability. It must not be used as proof of 5′ loss or intentional iBAR absence. If a reference/read lacks enough downstream anchor information for the new rescue, report that rescue as unavailable while retaining the established parsing behavior.

- Identify the spacer before using iBAR. Default iBAR handling reports agreement or disagreement; an optional mode requires agreement.
- An exact foreign iBAR must remain visible. Do not overwrite a spacer identity with its expected iBAR owner or automatically learn a replacement linkage map.
- Distinguish absent-by-design iBAR from missing sequence caused by read truncation. Require sufficient downstream sequence to establish the layout.
- Keep guide, target and construct resolution separate. Shared spacers can identify a gene without identifying a construct. Ambiguous evidence must not be copied into multiple constructs.
- Keep iBAR correction, UMI correction and missing-partner inference separate. Identical UMIs in the two positions do not establish physical linkage.

## 2. Automatic settings from a bounded pilot

Extend the existing preflight rather than adding an unrelated configuration system.

1. **Compile reference limits.** Derive scaffold layouts, iBAR lengths, suffix collisions and allowable correction radii per guide/position. These are safety limits, not quantities to learn from observed abundance.
2. **Collect up to 100,000 reads per library.** A quick sample from the beginning of each input file can give an initial recommendation. For automatic activation, prefer a deterministic sample collected during the existing Step 2 pass across the files. Record file coverage, unique cell/UMI support and sampling bias. Do not describe early-read samples as random whole-library samples.
3. **Separate discovery and validation.** Split by cell/UMI family so PCR copies cannot appear in both. When UMIs are unavailable, raw-R1 preflight remains diagnostic; it must not claim molecule-level validation.
4. **Discover structure independently of the current TSO endpoint.** Measure endpoint alternatives, scaffold support, retained lengths, sequence quality, iBAR agreement and competing identities. The current preflight is conditioned on the original TSO finder and would otherwise miss the failure mode being fixed.
5. **Select the least permissive passing policy.** Use unique-molecule recovery, masked known-guide accuracy, ambiguity and negative-control rates. If two settings recover nearly the same usable evidence, choose the stricter one. Test 18-base retention before admitting 17-base candidates.
6. **Freeze and record.** Save effective per-guide limits, validation evidence, unavailable reasons and reference/settings fingerprints. Reuse the frozen contract across dedicated gRNA and optional GEX-derived guide evidence. An inconclusive pilot leaves the affected rescue class in audit mode.

Do not use target expression or knockdown to tune these settings. Biological validation follows the frozen choice. iBAR coupling is a user-selected assumption, not something silently enabled because a small pilot shows few discordant reads.

Implemented configuration interface:

```yaml
grna_recovery:
  mode: auto                    # off | audit | auto | manual
  tso:
    endpoint_policy: joint
    max_edit_distance: 3
    partial_tso: false
  spacer:
    max_5prime_loss: 3
    min_retained: auto           # evaluate 18, then 17 where supported
    retained_edit_distance: 0
  ibar_corroboration: report     # report | require
  pilot:
    sample_reads: 100000
    split_by: cell_umi
    seed: 21
  resolved_contract: null
```

Manual settings may narrow or explicitly test a policy; they must not force an ambiguous identity. Exact thresholds for admitting a newly recovered unit to a cell are also benchmarked at 1, 2 and 3 rescue UMIs while retaining the existing pair threshold. The current three-UMI KD threshold is not itself an established optimal assignment threshold.

## 3. Implementation work

| Component | Planned change |
|---|---|
| `scripts/guide_design.py` and a shared recovery module | Reference-derived suffix index, candidate boundaries, ambiguity and evidence classes |
| `scripts/03_match_grna_ibar.py` | Replace the single-endpoint gate for recovery; share the implementation across legacy and design-aware layouts |
| `scripts/design_preflight.py`, `scripts/grna_preflight.py` | Independent structural discovery, family-level validation and frozen recovery policy |
| `scripts/02_extract_grna_barcodes.py` | Optional bounded pilot sampling during the existing pass; preserve barcode decoding behavior |
| `scripts/04_assign_grna_to_cells.py` | Joint molecule counting, rescue support, candidate construct handling and transition summaries |
| `scripts/05_build_anndata.py` and QC | Preserve recovery class, uncertainty, resolved policy and observed/inferred status |
| Optional Step 1c and Q1 | Reuse the shared parser with source-specific adapters; keep compatibility checks consistent |
| `analysis/` | Recovered-cell cohorts, knockdown models and biological benchmark figures |

Record the observed spacer, candidate guide set, possible bases lost, TSO endpoints/scores, scaffold support, iBAR state, source and resolution basis. At molecule/cell level, retain normal-only, short-only, mixed and uncertain-length evidence. Deduplicate shared reads and UMIs across sources before counting.

Do not alter the pilot's GEX expression matrix or GEX read-removal criteria in this change. Any future expansion of guide-contaminant removal needs its own false-removal benchmark, especially for shorter sequencing chemistry.

## 4. Technical benchmark

Freeze the effective v2.0 source, configuration, reference, barcode inputs and retained cell set. The workspace contains pre-existing edits, so a version label alone is insufficient; save source hashes and configuration snapshots. Compare these arms on identical data:

| Arm | Change from frozen v2.0 |
|---|---|
| Baseline | Original extraction, matching and assignment |
| Endpoint only | Joint endpoints; intact spacer; original TSO limit |
| Spacer only | Original endpoint; 1–3-base loss |
| Combined conservative | Joint endpoints and loss; at least 18 bases retained |
| Combined extended | Same, allowing 17 bases where reference-safe |
| Frozen automatic policy | Settings selected without expression outcomes |
| Experimental sensitivities | Scaffold indels, 4–5 TSO edits or partial TSO, tested separately |

Use all four pilot libraries for full molecule/cell comparisons. Keep the existing 100k probes as discovery only; assess performance on held-out families and the rest of the data. Use the saved Gavin reference as a dense-codebook and decoupling challenge, and synthetic/asymmetric fixtures for layouts lacking a comparable biological dataset. Synthetic compatibility is not biological validation of those layouts.

Report extraction and identity rates with explicit denominators, unique molecules, newly assigned cells, recovered partners, additional pairs, singlet/multiplet transitions, lost or changed calls, per-guide yield, ambiguity, runtime and memory. Separate evidence already counted by v2.0 from genuinely incremental molecules.

Challenge the implementation with held-out intact molecules whose 5′ ends are masked, TSO substitutions/indels, repeated-G junctions, scaffold cross-pairs, retained-spacer errors, iBAR swaps, absent/truncated iBARs, barcode errors, ambient low-count evidence and ordinary GEX sequences. Include same-spacer/different-construct swaps, which cannot be resolved from spacer identity alone.

## 5. Biological benchmark

Use the unchanged full-depth pilot expression counts. Freeze guide calls before examining target expression. Compare recovered cells with non-targeting controls and normal-guide cells matched by guide design, library, treatment, expression depth and partner target where relevant.

Report distinct cohorts:

- **Endpoint-rescued intact guides:** test whether newly assigned cells show the expected on-target knockdown.
- **Shortened-only guides:** exclude any detected normal guide targeting that gene, across both positions, all iBARs and sources, including subthreshold evidence. Also exclude reads/molecules with a plausible full-length interpretation under the new boundary model. This is the direct extension of the user's abnormal-only efficacy question.
- **Mixed or uncertain-length recovery:** useful for validating assignment, but not evidence that a shortened guide itself is functional.
- **Recovered missing partners:** test the newly observed target while accounting for the already detected partner; do not treat inferred guides as observed evidence.

Show target expression, knockdown with confidence intervals, sample sizes and adjusted significance by recovery class and length tier. Include non-targeting and shuffled-label controls. Assess changes in cell quality and target coverage so enrichment for easier-to-detect cells is visible.

Use a prespecified count model with library/treatment and technical-depth adjustment; report independent replicate support where it exists. Four sublibraries are not automatically four biological replicates. Estimate power from target expression and cohort size before interpreting a negative result. The earlier audit's three eligible observations cannot establish shortened-only efficacy.

## 6. Release decisions and deliverables

Suggested gates to finalize before implementation:

- `off` reproduces baseline assignments and molecule counts; no duplicated evidence, fabricated partner UMIs or forced construct identities in any supported layout.
- Each automatically enabled rescue class/layout meets a proposed 99.5% lower confidence bound for identity accuracy on held-out known-identity challenges, using a one-sided 95% interval and at least 1,000 accepted independent families, with separate negative-control reporting. This is a challenge-set criterion, not an empirical guarantee of whole-library purity. Insufficient support is an inconclusive result.
- Review all target switches and quantify additional-pair burden. Recovery alone is not a passing result.
- Demonstrate on-target KD in adequately powered recovered-cell cohorts using prespecified targets and multiplicity control. A proposed meaningful-effect threshold is 20% KD; report uncertainty rather than treating that cutoff as proof. Inconclusive classes remain optional/audit-only and are labeled unvalidated biologically.
- Measure computational cost and document the tradeoff; do not infer scalability from the small pilot reference.

Deliver in four stages: shared parser and fixtures; reference-aware preflight and frozen policy; full technical and biological benchmarks; release documentation and recommended presets. Publish the evidence document with a concise main narrative and linked methods/results.

Old configurations should retain baseline behavior when recovery is unspecified. After the gates pass, ship a new recommended v2.1 preset using validated automatic recovery, plus explicit `off` and `audit` modes. If only endpoint rescue passes, recommend that component while leaving shortened-guide or weak-TSO classes experimental. This permits useful improvement without overstating biological validation.

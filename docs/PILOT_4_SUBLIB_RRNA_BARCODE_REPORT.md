# C1/C2/T1/T2 GEX Census Report: rRNA and Barcode Recovery

> **Report scope (v2.0 documentation update).** Results below retain the
> dataset, pipeline version and methods stated in their provenance. For current
> installation and commands, use the [quick start](QUICKSTART.md),
> [configuration reference](CONFIGURATION.md) and [UGER resource guide](UGER_RESOURCES.md).
> Historical counts/timings are not measurements of the new kallisto backend.
> Experiment-specific audit workflows remain in Git history and the local
> UGER archive; their original commands are provenance, not release entry points.
> Raw barcode–UMI recovery is a molecule proxy; it is not automatically
> a gene-resolved expression UMI count under the Step 1a correction/assignment policy.


> **Full four-sublibrary census.** This classifies every one of the
> 1,680,846,154 raw PE150 GEX read pairs from
> C1, C2, T1, T2 — a complete enumeration, not a sample. It supersedes
> the interim 100k C1 smoke report, now retired to Git history and the local archive.

## 1. rRNA abundance

**547,811,208/1,680,846,154 pairs
(32.59%) were rRNA-positive.**

Every count below requires direct alignment of the 150-nt R1 cDNA mate to the
pinned human rRNA reference: at least 50 aligned bases at ≥95% identity. STAR
local alignment soft-clips TSO or adapter sequence. R2 barcode, linker, probe,
and adapter sequence is excluded from the rRNA call. The reference excludes
ribosomal-protein mRNAs.

| Sublibrary | Raw pairs | R1-supported rRNA | rRNA fraction | Wilson 95% CI |
|---|---:|---:|---:|---:|
| C1 | 452,146,811 | 149,325,804 | 33.03% | 33.02%–33.03% |
| C2 | 363,646,678 | 121,621,234 | 33.44% | 33.44%–33.45% |
| T1 | 440,532,021 | 143,593,650 | 32.60% | 32.59%–32.60% |
| T2 | 424,520,644 | 133,270,520 | 31.39% | 31.39%–31.40% |
| **Combined** | **1,680,846,154** | **547,811,208** | **32.59%** | — |

The unweighted four-library mean is
32.61% (95% Student-t CI
31.21%–34.02%),
which describes uncertainty in this pilot mean rather than a prediction for a
future library. Two looser calls bracket the headline count from above:
accepting an rRNA call on either mate gives
554,332,182 pairs
(32.98%), and relaxing the R1 threshold
to the sensitivity tier of at least 30 aligned bases at ≥90% identity gives
579,243,804 pairs
(34.46%).

| rRNA family | R1-supported pairs | Share of rRNA pairs |
|---|---:|---:|
| 18S | 366,535,589 | 66.91% |
| 28S | 112,212,465 | 20.48% |
| MT-RNR2 | 50,583,005 | 9.23% |
| MT-RNR1 | 17,498,570 | 3.19% |
| Other annotated rRNA regions | 981,579 | 0.18% |

## 2. Barcode extraction

Step 1a accepted **1,117,772,521 pairs
(66.50%)** and rejected
**563,073,633
(33.50%)**.
No read pair in any of the four libraries was too short for the
complete 58-base barcode amplicon, and only 33,035 pairs
(0.0020%) carried an `N` inside a
fixed barcode slice.

| Fixed R2 slice | Exact | Uniquely corrected | Ambiguous | No whitelist match ≤2 |
|---|---:|---:|---:|---:|
| BC1 | 1,084,428,541 (64.52%) | 209,698,255 (12.48%) | 95,605,900 (5.69%) | 291,112,764 (17.32%) |
| BC2 | 1,371,904,454 (81.62%) | 125,611,868 (7.47%) | 36,085,166 (2.15%) | 147,213,304 (8.76%) |
| BC3 | 1,528,913,130 (90.96%) | 62,620,284 (3.73%) | 33,662,209 (2.00%) | 55,649,552 (3.31%) |

Among rejected pairs, **413,239,715
(73.39%)** failed at
least one independent barcode round. The other
**149,833,918
rejected pairs
(26.61%)**
passed all three slices independently but failed the production
joint-combination support rule.

| Sublibrary | Rejected pairs | Failed ≥1 independent round | Passed all slices, failed combination support |
|---|---:|---:|---:|
| C1 | 145,126,058 | 107,050,120 (73.76%) | 38,075,938 (26.24%) |
| C2 | 119,122,192 | 88,677,257 (74.44%) | 30,444,935 (25.56%) |
| T1 | 144,874,193 | 103,633,353 (71.53%) | 41,240,840 (28.47%) |
| T2 | 153,951,190 | 113,878,985 (73.97%) | 40,072,205 (26.03%) |
| **Combined** | **563,073,633** | **413,239,715 (73.39%)** | **149,833,918 (26.61%)** |

At census scale the independent pass-mask is tabulated over all raw pairs
rather than per barcode outcome, so the per-round failure counts share that
denominator:
BC1 386,719,358 (23.01%), BC2 183,329,832 (10.91%), BC3 89,312,740 (5.31%).
Those counts are not exclusive to rejected pairs, but only
27,509,593 accepted pairs
(2.46% of accepted pairs)
carry any independent-round failure, which places a floor of
92.89% (BC1), 84.99% (BC2), 69.20% (BC3)
on the rejected share of each count.

| Sublibrary | Accepted pairs | rRNA in accepted | Rejected pairs | rRNA in rejected |
|---|---:|---:|---:|---:|
| C1 | 307,020,753 | 30.22% | 145,126,058 | 38.97% |
| C2 | 244,524,486 | 30.93% | 119,122,192 | 38.61% |
| T1 | 295,657,828 | 28.88% | 144,874,193 | 40.17% |
| T2 | 270,569,454 | 27.88% | 153,951,190 | 37.57% |
| **Combined** | **1,117,772,521** | **29.45%** | **563,073,633** | **38.82%** |

Across the census, 39.90% of rRNA pairs
failed barcode extraction, versus
30.40% of non-rRNA pairs (descriptive
risk ratio 1.31; four-library mean
1.31, 95% Student-t CI
1.22–1.40).
This is an association, not a mechanistic assignment.

## 3. rRNA captured by random-hexamer versus poly(T) Round-1 probes

All 1,117,772,521 barcode-valid pairs—including all
329,224,928 barcode-valid rRNA pairs—were assigned
unambiguously to an N6 or poly(T) Round-1 probe from the recorded BC1.
Probe-type ratios therefore use the barcode-valid cohort only; rejected pairs
with unresolved BC1 are excluded.

| Round-1 probe | Barcode-valid pairs | R1-supported rRNA | rRNA within probe | Share of valid rRNA |
|---|---:|---:|---:|---:|
| Random hexamer (N6) | 526,165,581 | 221,321,072 | 42.06% | 67.22% |
| Poly(T) | 591,606,940 | 107,903,856 | 18.24% | 32.78% |

Random-hexamer BC1 molecules therefore contributed **221,321,072/
329,224,928 (67.22%)**
of valid-barcode rRNA pairs, versus
**107,903,856/329,224,928
(32.78%)** from poly(T) BC1 molecules.
The same ordering holds in every sublibrary.

Here, apparent rRNA capture rate means the fraction of barcode-valid read
pairs assigned to a probe type whose R1 meets the primary rRNA criterion.
The N6/poly(T) rate ratio is the N6 fraction divided by the poly(T) fraction;
the parenthetical relative difference is `(rate ratio - 1) × 100%`. These are
read-pair metrics and can include RT, amplification, and sequencing biases.

| Sublibrary | Valid-pair N6:poly(T) | rRNA-pair N6:poly(T) | rRNA rate, N6 vs poly(T) | N6/poly(T) rate ratio (relative difference) |
|---|---:|---:|---:|---:|
| C1 | 0.85 | 1.94 | 43.34% vs 19.05% | 2.28× (127.54%) |
| C2 | 0.94 | 2.12 | 43.40% vs 19.21% | 2.26× (125.92%) |
| T1 | 0.91 | 2.11 | 41.17% vs 17.71% | 2.32× (132.40%) |
| T2 | 0.87 | 2.06 | 40.38% vs 17.03% | 2.37× (137.08%) |
| **Combined** | **0.89** | **2.05** | **42.06% vs 18.24%** | **2.31× (130.62%)** |

The probe × rRNA-family cross-tabulation was emitted for T1, T2 only,
so it is not reported here.

The unweighted four-library mean N6/poly(T) rRNA rate ratio is
2.31×
(95% Student-t CI
2.23–
2.39).

Compact source tables:
`output/pilot_4_sublib_audit/rrna_barcode_diagnostic/aggregate/` (local run output).
Reconciled `central_summary.csv` SHA-256:
`644b71d33b9667f84e46eef729c5dc80dcd02db2a3bec78ef23ab7a06a2c07a2`.
Shared rRNA reference-manifest SHA-256:
`af529502ede35d6a959e2250f1458374850df1a14947ff0d4ae63187b9b2627d`.
Per-library run-manifest SHA-256:

- C1 `509a16b936ffa9ee9cbe6173087cac4352ece9d92ebee6c4d48449b40c1bb71b`
- C2 `d8bab1a7a3cf8882c34dde29758713a1e805056c8505c56e04b6e3a07e6ec4ec`
- T1 `759042abab15364df91bd3e5e888a60067b77214a588d11521f07acd136b692c`
- T2 `c3b62159954e94d29cd61afd6a1e86c5463ed67d444107e8c27e56af17870c22`

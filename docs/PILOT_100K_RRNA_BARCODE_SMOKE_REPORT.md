# C1 GEX 100k-Pair Smoke Report: rRNA and Barcode Recovery

> **Report scope (v2.0 documentation update).** Results below retain the
> dataset, pipeline version and methods stated in their provenance. For current
> installation and commands, use the [quick start](QUICKSTART.md),
> [configuration reference](CONFIGURATION.md) and [UGER resource guide](UGER_RESOURCES.md).
> Historical counts/timings are not measurements of the new kallisto backend.
> The maintained read-recovery and rRNA entry points are
> `submit/submit_gex_grna_read_recovery.sh` and `submit/submit_gex_rrna_audit.sh`.
> Earlier one-use wrappers are [archived](../archive/audit_workflows/2026_pilot/README.md).
> Raw barcode–UMI recovery is a molecule proxy; it is not automatically
> a gene-resolved expression UMI count under the Step 1a correction/assignment policy.


> **Interim C1 100k result.** This covers the first 100,000 PE150 C1 GEX
> records. It is not a random sample and is superseded by the full
> C1/C2/T1/T2 census in
> [PILOT_4_SUBLIB_RRNA_BARCODE_REPORT.md](PILOT_4_SUBLIB_RRNA_BARCODE_REPORT.md).

## 1. rRNA abundance

**32,619/100,000 pairs
(32.62%) were rRNA-positive.**

Every count below requires direct alignment of the 150-nt R1 cDNA mate to the
pinned human rRNA reference: at least 50 aligned bases at ≥95% identity. STAR
local alignment soft-clips TSO or adapter sequence. R2 barcode, linker, probe,
and adapter sequence is excluded from the rRNA call. The reference excludes
ribosomal-protein mRNAs.

| rRNA family | R1-supported pairs | Share of rRNA pairs |
|---|---:|---:|
| 18S | 23,057 | 70.69% |
| 28S | 5,393 | 16.53% |
| MT-RNR2 | 3,137 | 9.62% |
| MT-RNR1 | 1,005 | 3.08% |
| Other annotated rRNA regions | 27 | 0.08% |

## 2. Barcode extraction

Step 1a accepted **67,784 pairs
(67.78%)** and rejected **32,216
(32.22%)**. All 100,000 R2 reads were long
enough for the complete 58-base barcode amplicon, and none of the fixed
barcode slices contained `N`.

| Fixed R2 slice | Exact | Uniquely corrected | Ambiguous | No whitelist match ≤2 |
|---|---:|---:|---:|---:|
| BC1 | 65,826 (65.83%) | 12,098 (12.10%) | 5,385 (5.38%) | 16,691 (16.69%) |
| BC2 | 81,143 (81.14%) | 7,605 (7.61%) | 2,098 (2.10%) | 9,154 (9.15%) |
| BC3 | 90,774 (90.77%) | 3,722 (3.72%) | 2,006 (2.01%) | 3,498 (3.50%) |

Among rejected pairs, **23,825
(73.95%)** failed at least one
independent barcode round. The overlapping failure counts were BC1
20,772
(64.48%), BC2
11,010
(34.18%), and BC3
4,990
(15.49%). The other
**8,391 rejected pairs
(26.05%)** passed all three slices
independently but failed the production joint-combination support rule.

| Barcode outcome | Pairs | R1-supported rRNA | rRNA fraction |
|---|---:|---:|---:|
| Valid | 67,784 | 20,227 | 29.84% |
| Rejected | 32,216 | 12,392 | 38.47% |

In this window, 37.99% of rRNA pairs failed barcode extraction,
versus 29.42% of non-rRNA pairs (descriptive risk ratio
1.29). This is an association, not a mechanistic assignment.

## 3. rRNA captured by random-hexamer versus poly(T) Round-1 probes

All 20,227 barcode-valid rRNA pairs had an unambiguous Round-1 probe
barcode.

| Round-1 probe | Barcode-valid pairs | R1-supported rRNA | rRNA within probe | Share of valid rRNA |
|---|---:|---:|---:|---:|
| Random hexamer (N6) | 31,146 | 13,297 | 42.69% | 65.74% |
| Poly(T) | 36,638 | 6,930 | 18.91% | 34.26% |

Random-hexamer BC1 molecules therefore contributed **13,297/
20,227 (65.74%)** of valid-barcode rRNA pairs,
versus **6,930/20,227
(34.26%)** from poly(T) BC1 molecules.

Compact source tables:
[`output/pilot_4_sublib_audit/rrna_barcode_diagnostic/smoke_C1_100k/`](../output/pilot_4_sublib_audit/rrna_barcode_diagnostic/smoke_C1_100k/).
Source run-manifest SHA-256: `54570e5d0d67e5f1e01383f01732318fc51b94db427696102208204ac6cedf2b`.

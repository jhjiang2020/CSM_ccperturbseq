# v2.1 opt-in guide recovery

**4 November 2026.** The implementation and four-library benchmarks recover more guide evidence. Downstream expression signatures support biological consistency in a subset of recovered cells, while the primary cis-knockdown test and shortened-only efficacy remain inconclusive. Recovery therefore stays **opt-in**; existing defaults and pilot expression counts are unchanged. This guide covers activation, parameter changes, rationale and measured recovery. The [original plan](GRNA_RECOVERY_V2_1_PLAN.md) preserves the design criteria; the [illustrated benchmark](../analysis/output/invalid_construct_kd_union_20261023/fiveprime_20261104/v2_1/report.html) provides the detailed figures.

## Choosing a mode for a new run

| Choice | When it is useful | Benefit | Tradeoff |
|---|---|---|---|
| **Standard (`off`)** | Strict comparability with an existing analysis, or no convincing structural evidence among rejected reads | Established behavior; no additional guide calls | Discards some informative reads and cells |
| **Audit** | First evaluation of a new library, construct or read length | Measures eligible recovery without changing assignments | Extra processing; no immediate cell gain; audit alone does not validate a setting |
| **Automatic recovery** | The audit shows unique spacer/scaffold evidence and held-out identity checks pass | More guide molecules, cells and observed partners; settings adapt to the reference and pilot | Some singlets become multiplets; real-read assignment error and biological benefit need review |
| **Constrained automatic recovery** | Prefer a smaller change while retaining validation | Set `spacer.min_retained: 18` or lower `max_5prime_loss`; auto searches within that limit | May lose useful evidence; the 18-nt floor lost only seven newly assigned cells in this pilot |
| **Manual recovery** | A fixed policy has already been evaluated for the intended design | Explicit, reproducible limits | Bypasses automatic validation; unsupported choices need their own benchmark |

Start with an audit for an unfamiliar library. Inspect the fraction recovered,
ambiguous candidates, length uncertainty and iBAR disagreement. If useful
recovery is present, run automatic preflight and review the saved policy before
using recovered cells. An inconclusive pilot is a reason to keep standard calls,
not to lower validation criteria just to enable recovery.

Before adopting the relaxed result, compare strict and recovered assignments,
net singlet gain, new multiplets and support from independent molecules. Check
cis-KD where expression permits; otherwise use independently validated trans
signatures with matched strict/control cells. Keep strict and recovered cells
identifiable so conclusions can be checked with and without recovery. A matching
trans signature supports biological consistency but is not a measure of guide
assignment purity or proof that a shortened RNA was active.

Use `ibar_corroboration: report` when spacer–iBAR coupling is uncertain or when
comparing libraries without iBARs. Choose `require` only when that coupling is
trusted and losing discordant reads is acceptable. Neither choice lets an iBAR
resolve an ambiguous spacer. Intentionally absent iBARs remain supported;
missing sequence from a truncated read is not treated as absence by design.

## A. Rationale

Two problems can discard useful evidence: a spacer can be shorter than its design, or the parser can choose the wrong TSO endpoint despite recognizable guide structure. v2.0's design-aware parser improved construct handling but still selected that endpoint first. The v2.1 change evaluates TSO, spacer and scaffold together, retaining the existing three-edit TSO limit.

**5′ damage is plausible, not established.** U6-driven guides lack the conventional mRNA cap, but RNA structure, end chemistry and Cas binding also affect stability. Terminal modifications can improve guide performance, whereas mature sgRNAs resisted a 5′ exonuclease in a plant study. Neither finding proves when the pilot's shortening occurred. The relevant comparison is with capped mRNA templates; cDNA is DNA. [Hendel et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4729442/), [Cody and Scholthof, 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7536693/).

**Short spacers can retain biochemical activity.** Li et al. measured Cas9 cleavage with 16-nt guides **in vitro**. Cleavage varied by more than 50-fold with the DNA sequence beyond the RNA–DNA hybrid. This supports investigating shorter guides, but does not demonstrate cellular CRISPRi knockdown or justify a universal 16-nt mapping cutoff. The study itself distinguishes its in-vitro findings from previous cellular limitations. [Li et al., 2023, *Biochemistry* 62:2541–2548](https://doi.org/10.1021/acs.biochem.3c00250), [open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC10578059/).

**The TSO is introduced during reverse transcription.** Template switching and non-templated nucleotide addition can create junction variation; preparation damage and PCR/sequencing errors are other possibilities. Repeated G bases can make several boundaries plausible. The parser should retain that uncertainty rather than infer a unique amount of damage. [Wulf et al., 2019](https://pubmed.ncbi.nlm.nih.gov/31640989/), [Zajac et al., 2013](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0085270).

In the earlier pilot audit, normal and short reads co-occurred in **92.2%** of well-supported cell/position/iBAR groups containing short reads; their median short fraction was **16.0%**. Only **4.3%** of short-containing UMI families also contained normal reads. This is consistent with a shared guide origin, but cannot distinguish intracellular shortening from preparation artifacts. Recovering identity does not establish that a shortened RNA caused the observed phenotype.

## B. Activation and changes from standard mode

**A code upgrade is required; editing a v2.0 YAML alone is insufficient.** The
new mapper evaluates TSO endpoints, spacer identity and scaffold together.
Step 2 can collect the pilot sample, Step 3 selects and applies the policy,
Step 4 counts rescued and standard reads together, and Step 5 records the
provenance. Step 1c reuses the dedicated-library policy for GEX-derived evidence. After upgrading
these scripts, switching modes and tuning limits requires only configuration.

Omitting `grna_recovery`, or setting `grna_recovery.mode: off`, preserves standard
v2.0 behavior. Recovery is an additional pass over failed/unmatched reads;
existing successful calls and cell-assignment thresholds are preserved.
Do not widen `grna_structure.grna_lengths` or the general guide edit distance
to imitate this feature: those settings cannot distinguish a unique 5′ suffix
from other edits or competing guide identities.

The [full opt-in overlay](../examples/grna_recovery_v2_1_enabled.yaml) enables
both relaxed matching and incorporation of GEX-derived guide evidence. The
[audit overlay](../examples/grna_recovery_v2_1.yaml) changes only the matching
audit. Merge individual keys into the complete run configuration, preserving
its existing paths, construct and filter settings.

### Parameters

Keys below are relative to `grna_recovery` unless explicitly qualified.
Standard-mode values describe the four-library baseline; other designs can
have different existing settings.

| Parameter or rule | Standard mode | Opt-in recovery preset |
|---|---|---|
| `mode` | `off`, including when omitted | `audit` to inspect; `auto` to enable after validation; `manual` for a fixed policy |
| `tso.endpoint_policy` | Choose TSO endpoint before spacer extraction | `joint`: consider all compatible endpoints with spacer and scaffold; `original` is an ablation |
| `tso.max_edit_distance` | `grna_structure.tso_max_edit_dist: 3` | `3`, unchanged; full TSO only (`partial_tso: false`) |
| `tso.search_end` | First 50 R1 bases | `50`, unchanged |
| `spacer.max_5prime_loss` | Designed/configured lengths (20/21 nt in this pilot); no explicit suffix rescue | Up to `3` bases from each guide's own 5′ end; auto tests limits 0–3 |
| `spacer.min_retained` | No separate rescue floor | `auto` compares 18 and 17 nt; an integer fixes the floor |
| `spacer.retained_edit_distance` | General matching uses `grna_matching.max_lev_dist` (3 in this pilot) | `0`: exact retained sequence; unrelated general matching settings are unchanged |
| `spacer.min_quality` | No added recovery-specific quality filter | Every retained spacer base ≥Q30 |
| Scaffold support | Existing configured anchors | Both stem and downstream anchor required, using the existing combined mismatch limit (2 in this pilot) |
| Identity ambiguity | Existing mapper | Reject competing identities across the full reference, including plausible suffix competitors outside the accepted loss limit |
| `ibar_corroboration` | Existing guide/iBAR rules | `report` annotates agreement without requiring it for spacer recovery; optional `require` rejects disagreement/unknown iBAR |
| `gex_decontamination.step1c.enabled` | `false` | `true` in the full opt-in preset: identify guides in GEX filter-positive evidence |
| `gex_decontamination.combine_with_grna_arm.enabled` | `false` | `true` in the full opt-in preset: include accepted GEX guide evidence in final calls |
| UMI and cell thresholds | Existing `guide_assignment` settings | Unchanged; rescued reads join the same UMI clustering before assignment |

For a 20-nt design, a three-base loss permits 17 nt; a 21-nt design still needs
18 nt. Intact reference guides are not discarded merely because they are
shorter than the rescue floor. One read with both intact and shortened
interpretations is counted once and labeled length-uncertain.

The tunable automatic controls are `pilot.sample_reads` (up to 100,000),
`seed`, `min_validation_families`, `min_accuracy` and `near_optimal_fraction`.
The sample size, seed and validation criteria below are the benchmarked values.
Auto selects the loss limit and retention floor; it does **not** learn the TSO
edit budget, base-quality threshold or iBAR policy. Reducing quality/validation
requirements or using a floor below 17 nt is outside this benchmark.

### Include guide reads from GEX

**GEX decontamination and GEX guide assignment are separate.** Standard future
runs remove recognized guide contaminants from GEX, but do not add them to
guide calls. The full v2.1 opt-in enables both identity recovery and source
merging:

```yaml
gex_decontamination:
  enabled: true
  step1c:
    enabled: true
  combine_with_grna_arm:
    enabled: true
    combine_before_umi_collapse: true
    retain_source_provenance: true
```

Step 1c consumes the corrected-barcode evidence collected during full GEX
preprocessing and reuses the dedicated library's frozen recovery policy.
The preset retains the conservative GEX acceptance rules: guide/iBAR edit
distance at most one, compatible resolved iBAR in dual-guide mode, and no
decoupled calls. Design-aware layouts allow genuine iBAR absence. Short reads
that cannot establish identity remain unresolved.

Step 4 combines accepted dedicated and GEX reads before UMI correction. A shared
molecule is counted once, and source contributions remain available. Step 5
uses the combined calls while retaining dedicated and GEX matrices for comparison.
The benefit is additional evidence from a second library; the tradeoffs are
extra processing and possible new multiplets, and many GEX molecules may already
be represented in the dedicated library.

These switches do not broaden the GEX contaminant filter: only its captured
candidates enter Step 1c. Matching a larger fraction of that evidence is different
from scanning all GEX reads for additional candidates. For an existing run,
reuse compatible Step 1a corrected-barcode evidence. If it is unavailable,
the historical run needs a separate evidence-extraction/preprocessing pass;
retain the original expression matrix. See [existing-run support](GEX_DECONTAMINATION_V1_4.md#existing-split-pipe-runs).

`mode: audit` only audits the relaxed threshold. It does not disable separately
enabled GEX source merging. To audit without adding a new evidence source,
use the audit overlay on the standard configuration; apply the full opt-in
overlay when ready to change final guide calls. To return fully to standard
calls, disable recovery and both GEX identity/merge switches.

### Configure and inspect

After installing v2.1, merge the [example overlay](../examples/grna_recovery_v2_1.yaml) into a complete pipeline configuration. It starts in `audit` mode, which reports candidates without changing assignments. To opt in to assignment-changing recovery, use `mode: auto` as below:

```yaml
grna_recovery:
  mode: auto                    # off, audit, auto, manual
  tso:
    endpoint_policy: joint      # original is an ablation
    max_edit_distance: 3
    partial_tso: false
  spacer:
    max_5prime_loss: 3           # relative to each guide's designed length
    min_retained: auto          # compare 18 and 17 bases
    min_quality: 30
    retained_edit_distance: 0
  ibar_corroboration: report     # require is an optional stricter subset
  pilot:
    sample_reads: 100000
    seed: 21
    split_by: cell_umi
    min_validation_families: 1000
    min_accuracy: 0.995
    near_optimal_fraction: 0.99
  resolved_contract: null
```

`auto` uses a deterministic sample from the complete barcode-decoded library.
New Step 2 runs collect it during extraction; existing inputs require a streaming
pass. Each cell/UMI family belongs to discovery or validation, never both.
Settings are selected without expression data. The least permissive setting
within 99% of the best passing recovery is selected. Validation requires at
least 1,000 accepted families per challenged loss class and layout, a one-sided
95% accuracy lower bound of 99.5%, and a shuffled-spacer negative-control upper
bound no greater than 0.5%. These are synthetic challenge criteria, not a
measured false-assignment rate for experimental reads.

Run preflight separately to inspect its decision:

```bash
python scripts/03_match_grna_ibar.py --config config.yaml --library C1 --preflight-only
```

The policy and validation evidence are saved as
`grna_matches/C1_recovery.json`. To reuse that exact choice, set
`resolved_contract` to its path. Reference or setting mismatches are rejected.
An inconclusive automatic pilot retains baseline calls. Raw FASTQ preflight,
without validated cell barcodes and UMIs, cannot enable automatic recovery.

`audit` reports eligible candidates in the failed-read companion file without
adding them to counts. It uses the specified retention floor, or 18 bases when
the floor is `auto`; use `auto --preflight-only` to inspect automatic selection. `manual` applies the specified loss limit and retention
floor but still rejects ambiguous identities. Its settings bypass automatic validation. `off` preserves the original output schema and matching behavior.

The spacer determines guide identity. `ibar_corroboration: report` retains
unambiguous spacer evidence and records exact iBAR agreement, disagreement or
absence by design. `require` accepts a rescue only with agreement or genuine
absence by design. Neither mode corrects a spacer using iBAR identity. Unresolved
or foreign iBARs cannot force a construct assignment. Shared guide sequences
may remain ambiguous even when they identify a target gene.

Both scaffold anchors must be observed. No-iBAR and asymmetric-iBAR constructs
use the existing [design-aware configuration](ASYMMETRIC_GUIDE_DESIGNS.md)
(`guide_design.enabled: true`) and their reference layouts. A truncated read is not evidence of intentional
iBAR absence. The established single-guide short-read route remains available
when the new recovery lacks a downstream anchor.

Step 4 combines baseline and rescue reads before UMI filtering and correction.
Recovery annotations describe those same molecule clusters. `baseline_only`
means counted by the previous mapper; it does not establish an intact physical
spacer. `shortened_only` and `length_uncertain` describe the candidate boundaries.
Cell summaries are written beside the assignment outputs and copied into
AnnData, with frozen policies in `uns['grna_recovery_policies']`.

Step 1c reuses the dedicated library's frozen policy. Run its preflight first.
This changes guide evidence recovery, not which GEX reads are removed or the
pilot expression counts. Q1 amplicon primers are not TSOs; its existing library-QC
parser is a separate assay and does not automatically inherit this RT-junction
policy.

### Apply to an existing run

Use a separate output namespace for the relaxed call, preserving the original
GEX expression products and strict assignments for comparison. Reuse compatible
Step 2 barcode files; an existing file only needs a streaming pilot pass, not
re-extraction. Run Step 3 preflight, then matching, Step 4 assignment and Step 5
integration with the same configuration. If Step 1c guide recovery is used,
run it after the dedicated policy is frozen and before the combined Step 4.
The main launcher orders these jobs automatically when recovery is enabled.
A Step-1-only submission requires an already frozen policy.
Neither a GEX realignment nor retrospective GEX decontamination is required
solely to change this guide-mapping policy. Recovery requires `edlib` in the
production environment; keep other dependencies from the existing installation.

## C. Expected recovery: technical benchmark

### Inputs and automatic choice

These are measured gains for this pilot, not guaranteed gains for another reference, construct or read length.

The comparison scanned all **84,099,473** barcode-decoded dedicated-library reads from C1, C2, T1 and T2. Cell comparisons use the same **37,086** GEX-retained cells, reference and assignment thresholds throughout. The baseline is the effective v2.0 code/configuration snapshot saved before implementation.

Each library used a uniform 100,000-read pilot, split by cell/UMI family. Settings were selected without expression outcomes. All four selected joint endpoints, up to **three bases lost relative to each guide's own design**, and at least **17 retained bases**. Retained spacer sequence must match exactly at Q30 or better; competing reference identities are rejected. Full-length and short interpretations of one identity may coexist and are labeled uncertain. iBAR agreement is reported separately and cannot determine an ambiguous spacer identity.

### Measured recovery

| Setting | Added reads, all cells | Newly assigned cells | Net singlet gain | Recovered missing partners | Additional cell–pair calls |
|---|---:|---:|---:|---:|---:|
| Endpoint only | 4,001,925 | 339 | 87 | 434 | 673 |
| Spacer only | 2,714,002 | 287 | 65 | 534 | 585 |
| Combined, ≥18 nt | 4,772,949 | 430 | 124 | 679 | 847 |
| Combined, ≥17 nt / automatic | **4,931,644** | **437** | **127** | **694** | **862** |
| Require iBAR agreement | 4,710,521 | 437 | 127 | 694 | 862 |

Automatic recovery adds **5.86 percentage points** of all input reads. Among retained cells, matched reads rise from 54,979,254 to 59,475,425 and jointly deduplicated guide molecules from 2,146,814 to 2,367,497. Assigned cells rise from 31,246 to 31,683. No baseline cell–pair call is lost.

The gain is smaller for analyses requiring singlets: 433 newly assigned singlets are offset by **306 existing singlets becoming multiplets**, leaving 127 net additional singlets. Additional pairs could reflect real multiple constructs or contamination; the benchmark does not establish their origin. Requiring at least two or three incremental rescue molecules reduces the newly assigned-cell count from 437 to **264** or **130**, respectively. These are sensitivity counts, not a newly validated assignment threshold.

Endpoint recovery supplies most of the gain. The 18-nt setting retains nearly all cell recovery; allowing 17 nt adds seven new cells and 15 missing partners. Requiring exact iBAR agreement reduces read evidence but leaves resolved pair calls unchanged in this pilot, whose construct resolution already requires compatible iBAR evidence. This does not assume spacer–iBAR linkage is reliable in every library.

### Accuracy, compatibility and cost

- **Baseline parity:** all 54,979,254 retained-cell match rows and 434,797 molecule-count groups agree with production. Current-code recounts reproduce all four automatic-arm count tables. Old and rescued reads enter the same UMI correction; shared-source reads are counted once.
- **Automatic challenges:** every enabled loss class/layout passes the held-out criteria. The smallest one-sided 95% accuracy lower bound is **99.9378%**. Shuffled-spacer negatives yield zero accepted matches, with upper bounds of 0.0147–0.0161%. These synthetic bounds are not estimates of experimental false assignment.
- **Reference challenges:** the pilot and dense 10,502-spacer Gavin reference produce no wrong identities in the tested truncation/TSO-error challenges. Ambiguous suffixes, altered retained spacers, cross-scaffold constructs and missing downstream anchors are rejected. No-iBAR and asymmetric layouts pass fixtures; their biological performance has not been established here.
- **Regression checks:** the benchmark passed 66 recovery/design, GEX-guide processing and preflight tests. Release verification passes 96 pipeline checks and seven biological-method checks, including the launcher, opt-in presets and audit count preservation. AnnData metadata round-trips without changing expression. `off` matches the frozen baseline on 400,000 sampled reads.
- **Cost:** the shared six-arm scans take 7.0–15.6 minutes per library using four slots. Complete benchmark tasks take 27.5–52.4 minutes; scheduler maximum virtual memory is 3.57–3.68 GiB. This includes repeated counting and parity checks and is not an isolated estimate of production overhead.

A separate check applies the frozen policies to **84,415 previously archived GEX guide reads**, recovering another **6,120** matches, of which 5,432 belong to retained cells. This selected archive is not the full GEX denominator, and these reads are not added to the dedicated-library gains above. It supplies additional normal-guide exclusions for the biological audit. The 6,120 are newly matched archived reads, before the full preset’s GEX strict filter and joint source counting; they are not a measured incremental cell or UMI gain from the combined preset. The cell and biological benchmarks above/below use dedicated-library recovery. Additional biological validation of the combined calls remains separate. The original 150-base paired-end GEX matrix is preserved; retrospective decontamination is unnecessary for this comparison.

Archived pilot outputs (paths below are relative to the project workspace): [technical summary](../output/grna_spacer_scan_20261023/fiveprime_20261104/v2_1/technical_combined.csv), [library/read results](../output/grna_spacer_scan_20261023/fiveprime_20261104/v2_1/read_recovery.csv), [cell transitions](../output/grna_spacer_scan_20261023/fiveprime_20261104/v2_1/cell_transitions.csv), [held-out tests](../output/grna_spacer_scan_20261023/fiveprime_20261104/v2_1/automatic_validation.csv), [run manifest](../output/grna_spacer_scan_20261023/fiveprime_20261104/v2_1/manifest.json).

## D. Expected recovery: biological validation

The primary test requires **three incremental rescue UMI clusters**, an observed target guide and a compatible singlet assignment. Models compare with 1,518 non-targeting-control cells and normal-guide cells, adjusting for expression depth, library, guide and partner. Four sublibraries are technical partitions, not independent biological replicates.

For shortened-only efficacy, the exclusion is gene-wide: **any detected normal guide targeting that gene**, across positions, guides, iBARs and dedicated/GEX sources, excludes the cell–target observation. Plausible intact boundaries, mixed molecules and uncertain lengths also exclude it. Normal-guide comparisons are controls, not evidence included in the shortened-only cohort.

| Primary newly assigned-cell contrast | Cells | Estimated KD | 95% interval | Adjusted P |
|---|---:|---:|---:|---:|
| CCM2 | 13 | 9.9% | −242.5% to 76.3% | 0.878 |
| TLNRD1 | 12 | −23.0% | −224.8% to 53.4% | 0.878 |

Neither contrast establishes on-target KD. Missing-partner groups have only 1–8 eligible cells per target. The shortened-only cohort has **two cell–target observations** at this threshold, one TLNRD1 and one ZEB2, so no shortened-only effect can be estimated. No strictly intact-endpoint recovery survives the three-molecule threshold; four cell–target observations survive at one molecule. Optimistic count-based power for a 20% effect is below 10% in the primary recovery cohorts. Shuffled-label comparisons are also inconclusive.

One- and two-molecule sensitivity analyses are reported separately. At one molecule, mixed/uncertain-length recovery shows exploratory KD for JUN (27 cells; 79.5%, 95% CI 37.6–93.2%), LAMB2 (24; 52.9%, 15.9–73.6%) and MEF2C (22; 79.9%, 32.8–94.0%). Each has adjusted P = 0.029 within this sensitivity family. These findings do not validate shortened-only activity or satisfy the primary release gate. Full estimates, matched normal-guide effects, cohort sizes, exclusions and model failures are retained, including negative results.

Data: [illustrated report](../analysis/output/invalid_construct_kd_union_20261023/fiveprime_20261104/v2_1/report.html), [all KD estimates](../analysis/output/invalid_construct_kd_union_20261023/fiveprime_20261104/v2_1/kd_results.csv), [cohort sizes/power](../analysis/output/invalid_construct_kd_union_20261023/fiveprime_20261104/v2_1/cohort_sizes.csv), [expression summaries](../analysis/output/invalid_construct_kd_union_20261023/fiveprime_20261104/v2_1/descriptive_expression.csv), [earlier normal/short audit](../analysis/output/invalid_construct_kd_union_20261023/fiveprime_20261104/report.html).

### Do downstream expression changes agree with strict cells?

**Yes for some recovered cohorts.** This additional analysis learns expression signatures from 15,123 strict complete singlets, holding out 25% of cells for validation. Recovered cells are excluded from training and validation. The target is compared with a non-targeting replacement in the same partner-guide background; any control with detected normal-guide evidence for that target is excluded. All 13 directly perturbed genes are removed from the expression signatures, along with mitochondrial and ribosomal protein genes. Models account for library, expression depth, guide design and partner.

Seven of 12 tested signatures reproduce in held-out strict cells: **CCM2, CEBPB, JUN, KLF2, MEF2C, NFKB1 and TLNRD1**. ADAMTS7, FDX1, KLF4, LAMB2 and ZEB2 are uninformative here; MAPK7 lacks sufficient matched strict cells.

| Newly recovered cells | Minimum rescue molecules | Cells | Downstream signature effect, control SD units | Adjusted P |
|---|---:|---:|---:|---:|
| CCM2 | 1 | 33 | 0.91 (95% CI 0.48–1.35) | 0.00033 |
| CCM2 | 2 | 18 | 0.94 (0.26–1.62) | 0.014 |
| TLNRD1 | 1 | 30 | 0.64 (0.16–1.11) | 0.037 |

At three rescue molecules, target-level tests remain inconclusive. An **exploratory pooled analysis**, added after inspecting those tests, combines 29 cells across CCM2, MEF2C and TLNRD1. Their signature effect is **0.50 control SD** (95% CI 0.13–0.87; adjusted P = 0.012), compared with **1.29 SD** in matched strict cells. Thus the direction agrees, but the average response is weaker. Pooling weights reference cells to the recovered target distribution and accounts for repeated observations from the same cell. These units measure expression-pattern strength, not percent knockdown or assignment purity.

Recovered cells have lower expression depth and gene detection than strict assigned cells, reinforcing the need for these adjustments. Trans agreement supports retained perturbation-associated biology, but shared responses or cell-state differences remain possible. No shortened-only cohort has enough observations for a target-level trans test. This analysis does not establish that shortening preceded guide activity, or that the shortened RNA itself was active.

Data: [trans-expression report and figure](../analysis/output/invalid_construct_kd_union_20261023/fiveprime_20261104/v2_1/trans/report.html), [strict validation](../analysis/output/invalid_construct_kd_union_20261023/fiveprime_20261104/v2_1/trans/strict_validation.csv), [target contrasts](../analysis/output/invalid_construct_kd_union_20261023/fiveprime_20261104/v2_1/trans/recovery_alignment.csv), [pooled contrasts](../analysis/output/invalid_construct_kd_union_20261023/fiveprime_20261104/v2_1/trans/pooled_alignment.csv), [cell quality](../analysis/output/invalid_construct_kd_union_20261023/fiveprime_20261104/v2_1/cell_quality.csv), [loss-length breakdowns](../analysis/output/invalid_construct_kd_union_20261023/fiveprime_20261104/v2_1/loss_tier_cohorts.csv).

## Recommendation and limits

The release keeps **recovery and GEX source merging off by default**, with an audit preset and an explicitly experimental full opt-in preset. Technical recovery and some downstream phenotype consistency are demonstrated; biological validation is insufficient to promote either endpoint or shortened-spacer recovery to a validated default. More recovered cells with adequate target expression, or a targeted experiment, are needed for that decision.

The 16-nt paper strengthens the rationale, not the release evidence. The preset does not introduce a 16-nt floor or increase the TSO edit limit. Partial TSO and scaffold indels are unsupported. Q1 amplicon QC retains its separate primer-based parser. Real-read false-assignment rates, biological validation of asymmetric/no-iBAR libraries, and production-scale overhead on a dense reference remain unresolved.

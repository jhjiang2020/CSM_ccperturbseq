# v2.1 optional gRNA recovery

v2.1 can recover guide evidence lost to a short spacer or an uncertain TSO boundary, and can include guide reads captured from GEX. Configure these options in the **same YAML used for the rest of the run**. They remain off by default: the four-library pilot gained guide molecules and cells, but some singlets became multiplets and biological support was mixed.

## Configure once

A **v2.1 code upgrade is required**. Changing a v2.0 YAML alone cannot add the joint TSO/spacer/scaffold matching or restart checks. After upgrading, edit the annotated entries already present in the standard config templates; no overlay file or separate audit run is required. If you used an earlier v2.1 draft config, replace its `grna_recovery` section with these entries; do not keep both.

```yaml
# Edit these entries inside the existing sections; retain your construct settings.
grna_structure:
  relaxed_structure: auto       # false = standard; auto = discover; true = use limits below
  relaxed_max_5prime_loss: 3     # at most 3 nt shorter than each designed spacer
  relaxed_min_retained: 17      # retained-length floor; auto may choose a higher floor
  relaxed_max_tso_shift: 3      # endpoint may move at most 3 nt from the standard TSO call
  relaxed_min_quality: 30       # every retained spacer base must meet this Phred score
  relaxed_require_ibar: false   # report iBAR agreement; true requires it where iBAR exists

gex_decontamination:
  # Keep the existing GEX filtering settings, including the top-level enabled flag.
  step1c:
    enabled: true               # identify guides in captured GEX evidence
  combine_with_grna_arm:
    enabled: true               # include those reads in final guide calls
    combine_before_umi_collapse: true
    retain_source_provenance: true
```

Use `auto` when the abnormal structures are unknown. It samples up to 100,000 decoded reads **per sublibrary**, reuses the existing construct/iBAR orientation diagnostics, and evaluates combinations of spacer loss and TSO endpoint shift. Multiple supported patterns can contribute; it does not assume that every abnormal read has one fixed shift or loss. It checks guide identity against the complete reference, including competing shortened spacers. Expression outcomes never select settings.

The pilot chooses **one shared parameter set for every sublibrary and for both dedicated and GEX guide reads**. Each sublibrary must pass held-out checks; within that shared set, auto chooses the least permissive setting within 1% of the best pooled recovery. Validated iBAR orientation corrections are shared too. Conflicting strand evidence stops the run; insufficient recovery-validation evidence keeps standard calls for all sublibraries. The decision, pattern counts and validation results are saved in `grna_matches/relaxed_structure.json`.

`true` uses the stated limits directly, for a setting already evaluated for the intended reference. It still requires unambiguous guide identity but skips automatic parameter selection. `false` retains standard matching. The other `relaxed_*` entries are limits for auto and fixed values for `true`; they do nothing when the option is false. For example, a 20-nt design may retain 17 nt after three-base loss, whereas a 21-nt design needs at least 18 nt. The benchmark does not validate a universal 16-nt cutoff.

| What changes | Standard workflow | With relaxed structure |
|---|---|---|
| Spacer length | Designed/configured lengths | Also a unique 5′-shortened suffix within the selected limits |
| TSO boundary | TSO endpoint chosen before spacer matching | Compatible endpoints evaluated jointly with spacer and scaffold |
| Retained spacer | Existing general error correction | Additional calls require an **exact** retained sequence at the configured quality floor |
| Scaffold | Existing construct parser | Both stem and downstream anchor must be observed |
| iBAR | Existing construct assignment rules | Spacer recovery reports agreement independently; optional corroboration cannot resolve an ambiguous spacer |
| Molecule/cell thresholds | Existing `guide_assignment` values | Unchanged |

Keep `grna_lengths`, the general guide edit distance and the existing TSO edit budget unchanged. The pilot used a three-edit TSO budget; relaxation changes endpoint selection, not that budget. Partial TSO rescue and scaffold indels are unsupported. Full and short interpretations of the same read can coexist and are recorded as uncertain rather than assigned a definite damage length.

No-iBAR and asymmetric-iBAR designs use their existing [design-aware reference settings](ASYMMETRIC_GUIDE_DESIGNS.md). Intentional absence is supported, including when corroboration is required. Failure to sequence an expected iBAR is not absence by design. With corroboration off, a unique spacer can be retained despite iBAR disagreement, but it cannot force an unsupported construct identity. Step 1c retains its existing, stricter GEX acceptance rules, so some such reads still do not enter combined calls.

The two GEX switches are independent of structural relaxation. **Decontamination alone does not add reads to guide calls.** Step 1c identifies guides in the evidence already captured by the GEX filter; `combine_with_grna_arm.enabled` adds accepted evidence. These switches do not broaden the contaminant filter or scan all GEX reads anew.

## Run or restart

For a new run, follow the [quick start](QUICKSTART.md), make the choices above while configuring the run, and use the usual launcher:

```bash
bash submit/submit_pipeline.sh --config /absolute/path/to/config.yaml
```

The launcher runs one shared pilot after Step 2, then maps all sublibraries with its saved policy. Step 1c uses that same policy. Step 4 combines read-level evidence from dedicated and GEX libraries **before UMI filtering and correction**, then assigns cells. Step 5 rebuilds AnnData from the existing GEX products and new guide results.

For a completed run, edit the same config and rerun the guide stages:

```bash
bash submit/submit_pipeline.sh --config /absolute/path/to/config.yaml --step 3,4,5
```

This also reruns Step 1c and its aggregation when enabled. It reuses Step 2 barcode files (include Step 2 in the restart if they are no longer available) and existing GEX expression products; it does not schedule GEX alignment or expression counting. An older Step 2 file without a cached pilot sample needs one streaming scan. On another scheduler, run the shared pilot once with `python scripts/03_match_grna_ibar.py --config config.yaml --preflight-only` before launching per-library matching. A per-library pilot alone cannot establish the shared policy.

**Check these restart conditions:**

- Captured, corrected-barcode GEX evidence must exist to use Step 1c. Otherwise use the documented [existing-run evidence extraction](GEX_DECONTAMINATION_V1_4.md#existing-split-pipe-runs); enabling a switch does not recreate missing reads. Preserve the original expression matrix.
- Step 1c hashes the full config. Rerun it after changing the YAML; prior Step 1c matches/manifests are archived before rebuilding. Do not edit hashes by hand. GEX completion validation permits proven changes only to the guide-only `step1c.enabled` and `combine_with_grna_arm.enabled` switches; GEX filter, alignment, counting, reference and sample changes remain incompatible.
- Step 3 binds the shared policy to the reference, settings, sublibrary list and decoded inputs. Step 4 rejects stale or incomplete relaxed matches. Step 5 checks that the counted products still match their inputs, settings and shared policy. Rerun all sublibraries together after changing settings.
- **Step 5 replaces `combined_adata.h5ad`; downstream annotations are not merged into the new file.** The launcher and Step 5 warn about this. Step 5 first writes the rebuilt file successfully, then preserves the complete previous AnnData and its annotation inventory under `grna_restarts/`. Reapply downstream analysis to the new result as appropriate.

Relaxed reads retain the observed `grna_seq` and record the designed sequence in `grna_seq_corrected`. Existing `grna_dist` records the sequence difference, including missing 5′ bases, and `match_type` identifies `relaxed_structure`. One `structure_correction` field records plausible losses, TSO endpoints/shifts, iBAR agreement and the shared policy. This distinguishes structural recovery from ordinary spacer-error correction without adding a column for every diagnostic.

Step 4 summarizes standard-only, relaxed-only and mixed evidence on the **same UMI clusters used for counting**, in a compact count-table field and cell-summary sidecars. Step 5 stores the shared policy and source-specific molecule summaries in `uns['grna_structure']`, without adding a series of recovery columns to `obs`. These records support comparisons excluding relaxed calls or requiring additional independent molecules. Standard mapper support does not itself prove an intact physical RNA.

## Why this may help

A short spacer can arise before or during library preparation. U6-driven guide RNAs lack the conventional mRNA cap, making 5′ damage a plausible explanation, but RNA structure and Cas binding also affect stability. The relevant comparison is capped RNA, not cDNA. The pilot does not establish when shortening occurred. [Hendel et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4729442/), [Cody and Scholthof, 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7536693/).

The TSO is introduced during reverse transcription. Junction variation, preparation damage, PCR errors and sequencing errors can complicate its boundary. Repeated G bases can support several endpoints. v2.0 improved construct handling but still chose the TSO endpoint first; this ambiguity remained recoverable. [Wulf et al., 2019](https://pubmed.ncbi.nlm.nih.gov/31640989/), [Zajac et al., 2013](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0085270).

Li et al. observed Cas9 cleavage with **16-nt guides in vitro**, with more than 50-fold sequence-dependent variation. This supports investigating shorter guides, but does not establish cellular CRISPRi knockdown or a mapping threshold. [Li et al., 2023, *Biochemistry*](https://doi.org/10.1021/acs.biochem.3c00250), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC10578059/).

Normal and short reads co-occurred in **92.2%** of well-supported cell/position/iBAR groups containing short reads; the median short fraction was **16.0%**. Only **4.3%** of short-containing UMI families also had normal reads. This is consistent with a common guide origin, but does not distinguish intracellular shortening from preparation artifacts or prove that a shortened RNA caused the phenotype.

## Expected gain and tradeoffs: four-library pilot

These are measured results for C1, C2, T1 and T2, not a promise for every design. The updated benchmark scanned **84,099,473** decoded dedicated-library reads, holding the **37,086** GEX-retained cells and assignment thresholds fixed. Auto selected one policy across all four sublibraries: up to **three bases of 5′ loss**, at least **17 retained bases**, and a TSO endpoint shift of at most **three bases**. iBAR corroboration was reported but not required. The shared pilot took about 24 minutes on one CPU; full matching and counting are additional work.

**Technical recovery**

| Metric | Standard | Relaxed benchmark | Change |
|---|---:|---:|---:|
| Matched reads, all decoded cells | 60,360,104 | 65,275,272 | **+8.14%** |
| Matched reads, retained cells | 54,979,254 | 59,459,931 | +8.15% |
| Guide UMIs, retained cells | 2,146,814 | 2,366,989 | **+10.26%** |
| Cells with guide assignments | 31,246 | 31,683 | **+437 / +1.40%** |
| Cell–pair assignments | 38,479 | 39,339 | **+860 / +2.23%** |
| Complete pairs | 25,922 | 26,952 | +1,030 |
| Singlets | 25,478 | 25,605 | **+127 / +0.50%** |

The additional 4,915,168 matched reads equal **5.84 percentage points of input reads**. Recovery found 694 missing partners and lost no baseline cell–pair calls. However, 433 newly assigned singlets were offset by **306 previous singlets becoming multiplets**. New pairs could reflect genuine multiple constructs or contamination; the benchmark cannot distinguish them. Requiring two or three incremental rescue molecules reduced newly assigned cells to 264 or 130.

Earlier ablations, before adding the explicit three-base endpoint-shift cap, help explain the gain: endpoint-only recovery added 4,001,925 reads and 339 assigned cells; spacer-only recovery added 2,714,002 reads and 287 cells. Their gains overlap. A combined 18-nt floor retained 430 of the 437 new cells. Requiring exact iBAR agreement reduced added reads to 4,710,521 but left resolved pair calls unchanged in this pilot. Optional corroboration therefore provides a useful sensitivity check without assuming linkage is reliable in every library.

Held-out synthetic identity tests had a minimum one-sided accuracy bound of **99.9378%**. Shuffled-spacer controls produced no accepted matches. These are not experimental false-assignment estimates. Tests also covered a dense 10,502-spacer reference and no/asymmetric-iBAR layouts; those layouts do not yet have biological validation here.

A separate selected archive of **84,415 GEX guide reads** yielded **6,120 additional matches**, before Step 1c's strict filter and joint source counting. It is not a full GEX denominator. **An additional cell/UMI gain from enabling GEX merging has not been established**, and must not be added to the table above. The pilot's 150-nt paired-end GEX matrix was preserved; retrospective decontamination was unnecessary.

**Biological validation**

Direct target expression was too sparse for a firm primary conclusion. At three incremental rescue molecules, newly assigned CCM2 cells (n=13) had estimated KD of 9.9%, and TLNRD1 cells (n=12) −23.0%; both intervals were wide and adjusted P=0.878. The shortened-only cohort contained just **two cell–target observations**. Any cell with a detected normal guide targeting the same gene was excluded from that cohort, across guides, positions, iBARs and evidence sources.

Trans-expression signatures provided partial support. Signatures were trained on strict cells and validated in separate strict cells, excluding recovered cells and all directly perturbed genes. Seven of 12 tested targets validated. Newly recovered CCM2 cells agreed with the strict signature at one- and two-molecule thresholds; TLNRD1 agreed at one molecule. Three-molecule target-level tests remained inconclusive. An exploratory pool of 29 CCM2/MEF2C/TLNRD1 cells had a **0.50 control-SD** signature effect (95% CI 0.13–0.87; adjusted P=0.012), compared with **1.29 SD** in matched strict cells: aligned direction, weaker response.

Recovery is therefore useful for retaining evidence and some biologically consistent cells, with the costs of extra processing, more multiplet calls and uncertain benefit for weakly supported assignments. These results do not prove shortened-only guide activity or establish assignment purity. Keep recovery identifiable in analyses; treat target-specific and shortened-only claims according to their available support.

Implementation checks reproduced the standard calls and counts in all four sublibraries. Recounting every retained-cell read with the compact production metadata reproduced the benchmark counts exactly. The pipeline suite passed 151 tests (three environment-dependent skips); all seven biological-method tests passed.

Detailed results: [updated technical tables](../output/grna_spacer_scan_20261023/fiveprime_20261104/v2_1/shared_structure_20261104/technical_combined.csv), [revision manifest](../output/grna_spacer_scan_20261023/fiveprime_20261104/v2_1/jobs/shared_structure_20261104/manifest.json), [updated cis-KD results](../analysis/output/invalid_construct_kd_union_20261023/fiveprime_20261104/v2_1/shared_structure_20261104/kd_results.csv), [updated trans-expression report](../analysis/output/invalid_construct_kd_union_20261023/fiveprime_20261104/v2_1/shared_structure_20261104/trans/report.html). The biological analyses were rerun with the revised assignments and counts; their conclusions were unchanged.

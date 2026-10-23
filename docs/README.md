# v2.0 documentation index

For a new run, begin with the [quick start](QUICKSTART.md). The pipeline does
not require access to the rgupta shared folder. Shared software/reference paths
are an optional Broad setup; the self-contained installation route uses your
own storage and Parse-provided software.

## Operating guides

| Guide | Purpose |
| --- | --- |
| [Quick start](QUICKSTART.md) | Install/select dependencies, configure an experiment and submit the first run. |
| [Installation and references](INSTALLATION.md) | Own production environment, optional shared prefix, Parse reference and kallisto reference. |
| [Configuration reference](CONFIGURATION.md) | Current defaults, stage ownership, GEX versus gRNA UMI settings and feature-level choices. |
| [Execution runbook](RUNNING.md) | Broad UGER, another Linux scheduler, dependencies, resources, outputs and restarts. |
| [UGER resources](UGER_RESOURCES.md) | Optional YAML CPU, memory, time, array concurrency and site settings for every maintained submission path. |
| [Architecture](../pipeline_overview.md) | Data flow and repository boundaries. |
| [Detailed workflow](DETAILED_WORKFLOW.md) | Stage inputs/outputs, reference/matching logic and AnnData integration. |
| [Counting choices](GEX_COUNTING.md) | STAR/kallisto × unique/fractional × three UMI modes; Parse defaults and Tony workflow comparison. |
| [GEX backend details](KALLISTO_BACKEND.md) | Native BUS counting, UMI policies, gene/target resolution and output contracts. |
| [GEX decontamination](GEX_DECONTAMINATION_V1_4.md) | Removal before either backend, optional Step 1c, source union and historical STAR BAM reuse. The filename preserves old links; the guide covers current behavior. |
| [Sequencing scenarios](ILLUMINA_100_CYCLE_SEQUENCING.md) | Read reach, short-GEX versus short-guide decisions and compatible identity models. |
| [Preflight and scaling](SCALABILITY_AND_PREFLIGHT.md) | iBAR diagnostics, CPU/memory controls and stage-specific constraints. |
| [Asymmetric guide designs](ASYMMETRIC_GUIDE_DESIGNS.md) | Optional position-specific iBARs, no-iBAR constructs, design capabilities and validated pilots. |
| [Migration](MIGRATING_TO_V2.md) | Preserve old results, choose new namespaces and rebuild affected stages. |
| [Release notes](../RELEASE_NOTES.md) | Versioned changes; historical entries describe their original releases. |

## QC guides

| Guide | Scope |
| --- | --- |
| [Q1–Q4 overview](../QC_overview.md) | Standalone library QC, recombination, expression/KD and fitness notebooks/modules. Uses the separate analysis environment. |
| [Recombination model](../QC_recombination_model.md) | Statistical model specification; its model version is independent of pipeline v2.0. |
| [GEX carry-over QC](GEX_GRNA_CONTAMINATION_QC.md) | Separate raw-read/controlled-STAR/BAM audit using production `splitpipe`. It does not enable production filtering by itself. |

Use gene-level products for gene-target QC. A primary transcript matrix is not a
drop-in replacement for every gene-based notebook. BAM-coordinate audits remain
STAR-specific; kallisto pseudoalignments do not supply those coordinates.

## Benchmarks and historical evidence

These reports retain their measured dataset, code/policy and resource scope.
They are evidence, not installation recipes. Current operating instructions
are the guides above. Old STAR experiments do not become kallisto benchmarks
when the production backend is changed.

| Report | Interpretation |
| --- | --- |
| [Asymmetric gRNA UGER validation](ASYMMETRIC_GUIDE_UGER_VALIDATION.md) | Small pooled pilot-derived 7/4-nt fixture, production-environment matching/count conservation, symmetric legacy agreement and regression tests. |
| [Measured v2 Step 1a benchmark](GEX_BENCHMARK_V2_REPORT.md) · [design](GEX_BENCHMARK_V2.md) | Completed 50M GEX pairs per pilot sublibrary: all 48 products, measured stage runtime/resources, 264 agreement comparisons and 72 fixed-cohort UMI-change comparisons; full checksum validation and scheduler accounting. |
| [Four-sublibrary technical audit](PILOT_4_SUBLIB_TECHNICAL_AUDIT.md) | Historical pilot processing/count validation and biological/QC context. |
| [GEX guide carry-over](GEX_GRNA_CARRYOVER_REPORT.md) | Recorded carry-over measurement for its specified input cohort. |
| [C2 misalignment case study](C2_GEX_GRNA_MISALIGNMENT_CASE_STUDY.md) | Controlled STAR evidence for short guide-derived false expression; retains its original experiment. |
| [Full rRNA/barcode audit](PILOT_4_SUBLIB_RRNA_BARCODE_REPORT.md) | Full pilot raw-read audit with its own barcode/rRNA definitions. |
| [rRNA protocol test](GEX_RRNA_PROTOCOL_TEST_REPORT.md) | Recorded protocol comparison, not a v2.0 backend comparison. |
| [Comparative rRNA audit](RRNA_AUDIT_COMPARATIVE_REPORT_TELOHAEC_THP1.md) | Cohort-specific rRNA/barcode comparison and provenance. |
| [UMI recovery/saturation](GEX_UMI_SATURATION_REPORT.md) | Barcode–UMI recovery/rarefaction; these molecule proxies are not automatically gene-resolved Step 1a UMIs. |

The current 50M-pair report replaces both earlier 10M-pair benchmark reports.
Superseded reports, assets and one-use audit tools remain in Git history and
an ignored local UGER archive. Historical full-audit evidence above remains
available with its measured definitions. The preliminary 100k-read rRNA smoke
report has been retired in favor of the full audit.

Reusable QC configuration schemas are under [`examples/qc/`](../examples/qc/README.md).
Copy and edit one for a run; production QC launchers require an explicit config.

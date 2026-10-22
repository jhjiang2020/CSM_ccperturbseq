# Migrate an existing run to v2.0

Use a complete v2.0 checkout and a copied run configuration. Parse remains
**1.6.4**; do not replace or patch its installation. The v2.0 changes live in
repository wrappers, native counting, validation and sparse combination.
There is no requirement to copy data or software from the rgupta folder.

## Preserve the baseline

Keep the old config, output, logs and reference identities read-only. Record
the code commit and environment/tool versions. Start the new configuration from
the matching public template, transfer your experiment-specific values, and
choose a new `paths.output_dir` and `paths.log_dir`.

Use a full release checkout rather than the old v1.4 manual script-copy list.
The new backend has additional runtime modules, manifest validators and a
parallel-counting module; copying only `run_kallisto_gex.py` is insufficient.
Follow [installation](INSTALLATION.md) for an independent environment, or keep
using an accessible validated production prefix without modifying it.

## Compatibility decisions

| Existing behavior | v2.0 action |
| --- | --- |
| Existing config omits `gex_alignment` | It still selects STAR and Parse GEX UMI correction. |
| Keep STAR/Parse counts and only use later guide/QC work | Retain that policy explicitly and validate the existing output before reuse. Legacy unmarked STAR output is supported only for Parse-default correction. |
| Select `directional` or `none` for STAR GEX | New output root; rerun molecule processing through the current Step 1 wrapper and rebuild Step 1b/5. Do not relabel the old matrix. |
| Select kallisto | Keep Parse preprocessing/reference; build a matching kallisto reference and run fresh Step 1a/1b. No historical STAR BAM is consumed by this backend. |
| Select fractional genes | Treat it as a distinct quantification policy; preserve floats and compare explicitly against unique-gene counts. |
| Need transcript output alongside standard analysis | Prefer `output_level: both`; gene counts stay primary. |
| Select primary transcript output | Check downstream feature assumptions. Gene-target QC and integer/gene-only methods may require the gene product instead. |
| Old guide arm already complete | Recreate or carefully validate/copy compatible guide products into the new namespace before Step 5; do not point new outputs over the historical run. |

A change only to count worker number is an execution change, but config identity
checks can still require a fresh validated count run. Changing reference,
backend, UMI correction or assignment semantics is a scientific change.

## Keep comparison inputs matched

For a STAR/kallisto comparison, hold raw FASTQ selection, Parse kit/chemistry,
barcode preprocessing, GEX decontamination/trimming and source annotation
constant. Kallisto uses transcript/target evidence; it need not give the same
counts as STAR. The optimization regression asks a different question: the
optimized kallisto counter must reproduce its earlier native-counter counts
for each chosen UMI/assignment mode.

Compare UMIs/genes per cell, called-cell overlap, a fixed cell cohort and gene
expression agreement. Keep runtime components and one-time reference building
separate. The [original four-library report](KALLISTO_PILOT_BENCHMARK_REPORT.md)
is the **serial-counter baseline**; it must not be cited as the optimized v2.0
runtime. Completed optimization results are published separately.

## Existing STAR BAM decontamination

The feature introduced in v1.4 remains available for compatible historical
STAR outputs. Set `gex_alignment.backend: star`, preserve the old source, and
use `gex_decontamination.existing_alignment` with distinct source/sandbox
paths. Raw FASTQs are still needed to reproduce filter decisions. The sandbox
removes positive QNAMEs, harmonizes alignment statistics and resumes Parse at
molecule counting without rerunning STAR.

This route cannot retroactively apply the newer survivor TSO trimming to reads
already aligned. Use fresh preprocessing/alignment when that trimming policy
is part of the experiment. It cannot be combined with kallisto.
[Full configuration and validation](GEX_DECONTAMINATION_V1_4.md#existing-split-pipe-runs).

## GEX filtering and optional guide recovery

Do not enable Step 1c/source merging solely because the GEX backend changed.
The required guide/iBAR bases must be sequenced. Keep those options disabled
for short dual-guide GEX R1 when identity is unobservable; the contaminant
removal filter can still run.

When enabling source union, set both Step 1c and
`combine_with_grna_arm.enabled`, retain `combine_before_umi_collapse: true`,
and provide complete dedicated Step 3 evidence in the new output root. Step 4
then unions reads before its own depth filter/directional collapse; Step 5
retains source-specific products and records which calls drive the result.
GEX transcriptomic `umi_correction` does not change this guide-UMI policy.

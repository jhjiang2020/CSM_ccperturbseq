# Guide counts, primary assignment and optional rescue

Steps 4–5 do not fit statistical guide models or require crispat. Step 4 writes
the deduplicated guide-unit ledger before the assignment cutoff and applies the
configured primary rule. Step 5 joins those products to the existing expression
result. Optional rescue belongs to the [analysis submodule](../analysis/docs/guide_rescue.md).

## Guide evidence and calls

The existing Step 4 filenames remain unchanged. For dedicated evidence the long
cache is `grna_assignment/grna_umi_counts.csv.gz`; for combined evidence it is
`grna_assignment/grna_umi_counts_by_source.csv.gz`. Combined evidence is merged
before UMI collapse. The cache retains UMI counts, observable unit identity,
source provenance where available, barcode state, and construct candidates.
Those columns describe different levels of evidence and should not be flattened
into independently sorted lists for analysis. Use the Step 4 assignment table
for the exact accepted construct rows and construct-arm counts.

In dual mode the primary observed matrix uses the resolved Step 4 pair calls.
For example, 3+3 passes a pair-sum threshold of five; a 5+0 accepted construct
calls only the five-UMI arm. Unresolved and guide-only rows remain in the Step 4
QC tables. In single mode UMIs for the same named guide are summed across units
before applying the fixed guide threshold and the configured `all`, `top`, or
`discard` policy. The returned rows keep each selected unit's evidence; sum
their `umi_count` values to recover the selected guide total.

## AnnData contract, version 1

| Location | Meaning |
| --- | --- |
| `.obsm['gRNA_counts']` | Unthresholded, cells × guides UMI counts for the chosen evidence source. |
| `.obsm['gRNA_assignment']` | Binary observed calls from Step 4; same axes as counts. Replaces the old secondary per-guide cutoff. |
| `.uns['gRNA_names']` | Guide column order, preserved for compatibility. |
| `.uns['grna_features']` | One row per guide column: guide ID, unambiguous spacer/target/position if known, identity status, barcode design state. |
| `.obs['grna_observed_umis']`, `.obs['grna_observed_n_guides']` | Cached UMI total and guide columns with observed molecules. |
| `.obs['grna_unindexed_umis']` | Cached UMIs outside the named matrix axis, if any. |
| `.obs['grna_unresolved_identity_umis']` | UMIs on non-resolvable guide identities, including unindexed evidence. |
| `.obs['grna_primary_n_guides']`, `.obs['grna_primary_guides']`, `.obs['grna_primary_targets']` | Count and JSON arrays of primary observed guide/target IDs. Empty sets use `[]`. |
| `.obs['grna_primary_n_resolved_constructs']` | Number of accepted construct identities (dual only). |
| `.uns['grna_inferred_partners']` | Cell, construct, position and expected unobserved guide, if any. These rows do not create guide calls. |
| `.uns['grna_construct_targets']` | Reference target identities per construct, for filtering potential carriers. |
| `.uns['grna_schema']` | Version, primary mode/threshold, evidence source, relative guide-cache and primary-call table paths, and the cache SHA-256. |

The guide table represents a legacy name mapping to multiple spacers, targets,
or observed positions as ambiguous. Its UMI evidence remains auditable; optional
rescue does not assign it to a particular identity. Symmetric and asymmetric
iBARs use the same contract; barcode applicability is per guide, not inferred
from a global length. For dual constructs without an iBAR, Step 4 must use
design-aware matching and have usable spacer/position evidence. A non-unique
spacer/position does not prove a construct. Single-guide runs do not require
construct linkage or a barcode to use the primary rule.

Older notebooks use fields such as `target_gene`, `combo_type`, `pos1_grna`,
`grna_names`, and `dedicated_*`. They are retained by default for compatibility.
Set `guide_assignment.legacy_anndata_metadata: false` for a compact AnnData with
the versioned summaries, distinct source matrices, and no duplicate combined
matrix or legacy guide observation aliases. The GEX observation columns and
Step 4 source files are retained with either setting. New analysis should
use the versioned matrices rather than parsing legacy aliases.

To upgrade from current Step 4 and GEX products without replacing the original
AnnData, rerun only Step 5 with `--destination path/to/new_primary.h5ad`. This
builds a pipeline snapshot; use the analysis rescue stage to retain existing
downstream annotations on a versioned primary object.

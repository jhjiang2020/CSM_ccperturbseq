# Archived 2026 pilot-audit workflows

These files preserve the exact single-use submission wrappers, validators, and
example extractors developed while diagnosing the four-sublibrary pilot. They
are retained for provenance and are not active production entry points; some
refer to their original paths and would need to be restored before rerunning.

Use the two maintained submission interfaces instead:

- `submit/submit_gex_grna_read_recovery.sh` — GEX/gRNA raw-read recovery,
  deterministic rarefaction, sequencing saturation, and integrated guide QC.
- `submit/submit_gex_rrna_audit.sh` — full-census GEX rRNA/barcode audits, with
  both the four-library `pilot` profile and a generic FASTQ-`manifest` profile.

The maintained analysis modules remain under `scripts/`:

- `pilot_screen_audit.py`
- `pilot_cell_guide_diagnostic.py`
- `rrna_barcode_diagnostic.py`
- `external_rrna_audit.py`
- `barcode_contract.py`
- `barcode_rejection_qc.py`
- `barcode_crossover_audit.py`

Archived project-specific wrappers include the already-submitted
`THP_1_sublib11` chain. UGER spools submitted job scripts, so archiving the
workspace copies does not alter jobs that were already queued or running.

# Standalone QC examples

These are configuration schemas for Q1–Q4 and the optional GEX/gRNA carry-over
audit. Copy the required example to a writable run folder and edit its paths,
construct sequences, thresholds and library settings. The original root-level
`qc_config*.yaml` files were pilot run configurations and are no longer shipped.
Local copies remain usable when passed explicitly.

| Example | Entry point |
| --- | --- |
| `qc_config.example.yaml` | Q2–Q4 modules and QC notebook |
| `qc_config_guideseq.example.yaml` | `submit/submit_q1_guideseq.sh` / Q1 CLI |
| `qc_config_gex_grna_contamination.example.yaml` | `submit/submit_qc_gex_grna_contamination_pipeline.sh` |

`base_config` and `sc_config` resolve relative to the QC YAML, so update them
after copying. Their checked-in values point to `config_dual.yaml` only to show
the inheritance schema. Q1 can run standalone after removing `sc_config` and
configuring its environment, paths and reference explicitly. See
[QC overview](../../QC_overview.md) and [GEX carry-over QC](../../docs/GEX_GRNA_CONTAMINATION_QC.md).

For GEX carry-over QC, pass `--config /absolute/path/to/my_qc.yaml` to the
launcher. Direct array/aggregate jobs require
`CSM_GEX_CONTAM_CONFIG=/absolute/path/to/my_qc.yaml`. There is no implicit
pilot-config fallback. Q1 similarly takes `CSM_QC_CONFIG`.

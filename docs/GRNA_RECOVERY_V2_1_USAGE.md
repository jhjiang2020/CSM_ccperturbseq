# v2.1 guide recovery

The consolidated [v2.1 guide](GRNA_RECOVERY_V2_1_EVIDENCE.md) contains the rationale,
standard-versus-relaxed parameter table, activation instructions, automatic
settings and four-library technical and biological benchmarks.

Recovery remains off when unspecified. Start with the
[audit overlay](../examples/grna_recovery_v2_1.yaml); change `mode` to `auto`
to opt in to assignment-changing recovery after automatic validation.

The [full opt-in overlay](../examples/grna_recovery_v2_1_enabled.yaml) also enables
GEX guide identity recovery and source merging into final calls.

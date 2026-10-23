# Small pilot-derived asymmetric fixture

This is one pooled **7/4-nt iBAR** validation dataset, not a new biological
screen. It uses 480 distinct R1 records from the C1, C2, T1 and T2 pilot
sublibraries. The builder scans at most the first 100,000 R1 records per
sublibrary and selects 20 complete constructs per sublibrary, with three source
reads per position. The pooled reference contains 23 distinct constructs.

- `asymmetric_R1_001.fastq.gz`: 1,920 derived reads. Position 1 is unchanged;
  position 2 loses the final three bases of its expected seven-base iBAR.
- `reference_7_4.csv`: matching synthetic reference, retaining iBAR2 bases 1–4.
- `control_R1_001.fastq.gz`, `reference_7_7.csv`: the unedited symmetric control.
- `*_barcodes.tsv.gz`: the same R1s with **artificial** cell/UMI labels for
  Steps 3–4. They are not outputs of Parse R2 barcode extraction.
- `read_truth.tsv.gz`, `cell_truth.csv`: independent boundary/count truth.
- `provenance.tsv.gz`, `manifest.json`: source filenames, source record numbers,
  hashed source identifiers/content, half-open deletion coordinates, and checksums.
- `config_template.yaml`: matching policy inherited from the pilot; the runner
  supplies fixture/reference/output paths and switches the symmetric control.

Expected coordinates are TSO offset 0 or 7, TSO length 30, spacer length 20 or
21, stem length 19, iBAR length 7, then the full configured tracr anchor. Exact
fixed sequences and reference identities select clean reads; the production
extractor is **not** used to build truth. Sequence and quality are shortened
together; all other bases are retained. The reference is normalized to the
dominant exact observed strand before truncating iBAR2.

Artificial cells contain both arms, only position 1, or only position 2. Each
observed arm has three distant 12-nt UMIs, each repeated twice. Expected totals
are 240 artificial cells and 960 molecules. These labels deliberately isolate
guide matching/counting from barcode demultiplexing.

In the production environment, from the repository root:

```bash
python -B scripts/validate_asymmetric_pilot.py validate \
  --fixture tests/fixtures/asymmetric_pilot_v2 \
  --output output/asymmetric_pilot_check --workers 2
```

Run on a compute node, or use `submit/validate_asymmetric_pilot.sh` on UGER.
The output directory must not already exist. To rebuild from the private pilot
FASTQs, use the `build` subcommand with `--config config.yaml --destination`
pointing to a new directory. No private raw FASTQs are required to validate the
checked-in fixture. See `docs/ASYMMETRIC_GUIDE_UGER_VALIDATION.md` for measured
results and limitations.

# v2.0 asymmetric gRNA arm: UGER validation

Date: 2026-10-23. Scope: a small production-environment correctness check, not a
new screen analysis or performance benchmark.

## Dataset

The [shipped fixture](../tests/fixtures/asymmetric_pilot_v2/README.md) pools real
gRNA R1 sequences from all four pilot sublibraries, C1, C2, T1 and T2. Their
original design uses seven-base iBARs at both positions. The only asymmetric
sequence edit is to **delete the final three expected iBAR2 bases**, yielding
one **7/4-nt** dataset. Position 1 is unchanged. Quality strings receive the
same deletion, and the synthetic reference retains the first four iBAR2 bases.
The original 7/7-nt records/reference provide a symmetric control.

The builder scans at most 100,000 reads per sublibrary, using exact fixed-offset
TSO/stem/tracr checks, 20/21-nt spacers, and exact reference identities to select
clean source reads. It does not call the new production extractor. Deletion
coordinates are `ibar_start + 4` through `ibar_start + 7`, zero-based, half-open;
all other sequence/quality characters remain unchanged.

| Fixture property | Value |
| --- | ---: |
| Source reads inspected | 100,000 per sublibrary; 400,000 total |
| Distinct source R1 records retained | 480; 120 per sublibrary |
| Construct examples | 20 per sublibrary; 23 distinct constructs after pooling |
| Distinct spacer sequences | 12, covering 20 and 21 nt |
| Derived reads | 1,920; 480 per source sublibrary |
| Position 1 / position 2 reads | 960 / 960 |
| Artificial cells | 240: 80 complete, 80 position-1-only, 80 position-2-only |
| Expected molecules | 960: three UMIs per observed arm, two read copies per UMI |

Cell identifiers and 12-nt UMIs are deliberately synthetic. The barcode TSVs
are controlled **Step 3 inputs**, not Parse R2 demultiplexing results. Source
filenames, record ordinals, hashed identifiers/content, edits, expected counts,
and fixture checksums are included in the fixture provenance files. Raw pilot
FASTQs are not needed to rerun validation.

## Production environment and execution

The unchanged production `splitpipe` environment was activated on UGER:

- Prefix: `/medpop/rgupta/Jiahao/conda/envs/splitpipe`.
- Python 3.12.12; NumPy 2.2.6; pandas 2.3.3; SciPy 1.17.0;
  edlib 1.3.9.post1; PyYAML 6.0.3.
- Worker: `submit/validate_asymmetric_pilot.sh`, RedHat8, two bound SMP slots,
  8 GB per slot, one-hour limit, numerical thread pools fixed at one.
- Fixture build: job **62581433**. Final validation: job **62581448**.
- Final accounting: exit status 0, scheduler failure 0, 35.900 seconds wall
  time, 1.165 GB maximum virtual memory on `uger-rh8-d004`.
- No dependency changes were made. Because pytest is not installed in this
  production environment, the existing unittest-based tests run with Python's
  standard-library test runner.

The runner executes raw-FASTQ preflight, two-worker Step 3 matching, and Step 4
assignment for the legacy 7/7 control, the opt-in 7/7 control, and the opt-in 7/4
dataset. It checks every read's extracted spacer/iBAR, design-aware positions
and candidate sets, per-cell molecule conservation, and construct compatibility
against independently generated truth. It also compares every legacy matching
output column between the two symmetric paths.

## Results

| Run | Matched reads | Cells | Molecules | Complete calls | Inferred calls |
| --- | ---: | ---: | ---: | ---: | ---: |
| Legacy symmetric 7/7 | 1,920 | 240 | 960 | 80 | 160 |
| Design-aware symmetric 7/7 | 1,920 | 240 | 960 | 80 | 160 |
| Design-aware asymmetric 7/4 | 1,920 | 240 | 960 | 80 | 160 |

- No lost/duplicated reads, boundary errors, per-cell count discrepancies,
  design-aware candidate/position errors, or unsupported construct calls.
- The raw-FASTQ pilot independently recovers 960 seven-base position-1 iBARs
  and 960 four-base position-2 iBARs. Both positions resolve to `as_reference`
  in discovery and held-out validation; pilot status is `ok`, with no warnings.
- The design-aware 7/7 output matches all legacy Step 3 columns exactly.
- All 46 observable units identify one construct in this particular subset,
  allowing 160 supported missing-partner calls. This is not a general guarantee
  that a four-base iBAR uniquely identifies a construct.
- The resolved 7/4 contract disables the unsupported legacy recombination
  likelihood, retains observable-unit UMI correction, and identifies position
  from the nonoverlapping iBAR lengths. The inherited scalar `ibar_length: 7`
  produces the expected warning and is superseded by the reference-derived
  per-position lengths; insufficiently separated iBARs receive reduced
  correction radii.
- **102 regression tests pass**, with zero failures, errors or skips. They
  include asymmetric, one-base and no-iBAR designs, ambiguous/shared units,
  truncated reads, scaffold conflicts, orientation pilots, Step 0 capabilities,
  source deduplication, linked compatibility QC and the packaged deletion check.

Machine-readable results, package versions, exact commands, source/test hashes,
raw pilot output, test logs and UGER accounting are preserved under
[`assets/asymmetric_pilot_v2_20261023/`](assets/asymmetric_pilot_v2_20261023/).
Earlier fixture-builder attempts 62581429/62581430 and a FASTQ-naming attempt
62581443 failed in the validation harness and were corrected. Job 62581444
passed the initial 101-test validation; the final run adds the explicit
whole-fixture sequence/quality deletion test. These attempts remain in the
accounting history, rather than being presented as successful runs.

## Reproduce

From the repository root, inside a compute allocation with `splitpipe` active:

```bash
python -B scripts/validate_asymmetric_pilot.py validate \
  --fixture tests/fixtures/asymmetric_pilot_v2 \
  --output output/asymmetric_pilot_recheck --workers 2
```

For UGER, create a fresh durable run directory, then submit the worker with
`PROJECT_ROOT`, `RUN_DIR`, `FIXTURE_DIR`, `CONDA_EXE` and `MODE=validate`
through `qsub -v`. Set `CONDA_EXE` to the configured Conda executable.
Use shared absolute paths; the worker writes job-specific success/failure
sentinels. Capture the `qsub -terse` job ID and log path. The validation output
directory must not already exist. `MODE=build` rebuilds from the private pilot
FASTQs using `config.yaml` and requires a new fixture directory.

## Limits

This deliberately clean, small, partly replicated fixture validates code paths,
boundaries and count accounting. It does **not** estimate full-screen extraction
sensitivity, sequencing-error tolerance, biological recombination, library
uniformity, cell assignment accuracy or throughput. Artificial cell co-detection
is not physical linkage. No-iBAR and one-base behavior is covered by regression
tests, not additional pilot-derived datasets. Parse R2 extraction, GEX alignment,
and full Step 5 H5AD serialization are outside this smoke test.

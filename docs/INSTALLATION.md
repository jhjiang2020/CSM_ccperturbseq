# Installation and references for v2.0

Use Linux with Bash 4+, a working Conda installation, and storage visible to
compute jobs. The production `splitpipe` environment and the pristine Parse
1.6.4 source distribution are required for **both** GEX backends. The optional
kallisto tools environment does not replace them.

Commands that install packages, build references or import/process scientific
data should run in a compute allocation. On Broad, initialize UGER and request
an allocation first, for example:

```bash
source /broad/software/scripts/useuse
reuse -q UGER
qrsh -q broad -pe smp 8 -binding linear:8 \
  -l os=RedHat8,h_vmem=8G,h_rt=12:00:00
```

Set paths again inside the allocated shell if its environment differs. On
another cluster, use its equivalent allocation command. Resource needs depend
on reference size; the example reserves 64 GB in total.

## Own production environment

Obtain **ParseBiosciences-Pipeline.1.6.4** from Parse Biosciences. The repository
does not redistribute it. Parse provides its download/installation route through
its [data-analysis support page](https://resources.parsebiosciences.com/support-suite/hc/en-us/articles/27066395947412-how-do-i-analyze-my-parse-biosciences-data).
Obtain the supported 1.6.4 archive rather than substituting the newest version:
the runtime overlays validate that version and its source checksums.

Extract the pristine tree at your chosen `paths.splitpipe_dir`. A conventional
layout is:

```text
CSM_ccperturbseq/
├── ParseBiosciences-Pipeline.1.6.4/   # obtained separately
│   ├── split-pipe
│   ├── splitpipe/
│   ├── setup.py
│   └── install_dependencies_conda.sh
├── scripts/
└── submit/
```

The tested production Python version is 3.12.8. In your own installation:

```bash
conda create -n splitpipe python=3.12.8
conda activate splitpipe
cd ParseBiosciences-Pipeline.1.6.4
bash install_dependencies_conda.sh --install --yes
python -m pip install .
cd ..
python -m pip install -r requirements_extra.txt
python -m pip check
python -c 'import edlib, pandas, pysam, scipy, anndata, splitpipe, yaml'
STAR --version
```

Keep the archive's installation instructions alongside it. Do not physically
apply the repository's `.patch` to Parse: production wrappers apply validated
changes in memory. Do not update a shared or previously validated environment
in place to follow this new-installation recipe.

Configure the Conda base and environment you actually installed:

```bash
conda info --base
```

Put that result in `environment.conda_prefix`; put `splitpipe` (or its absolute
prefix) in `environment.conda_env`. The base must contain `condabin/conda`.
`requirements_extra.txt` supplies PyYAML and the required C-backed `edlib`
matcher for production GEX decontamination. STAR, `pysam`, SciPy and AnnData
come from the Parse installation. Production/QC workers use this environment.

## Optional shared rgupta setup

This route is only for colleagues who can read the shared directories and are
authorized to use the Parse distribution. Other colleagues use their own
installation above; no shared folder is necessary.

The published standalone production prefix is:

```yaml
environment:
  conda_prefix: /medpop/rgupta/Jiahao/conda/envs/splitpipe-share
  conda_env: /medpop/rgupta/Jiahao/conda/envs/splitpipe-share
```

Both fields intentionally point to the same standalone prefix. It is a
read-only dependency, not your output or scratch directory. Check its completion
manifest and imports **in a compute allocation**:

```bash
SHARED_ENV=/medpop/rgupta/Jiahao/conda/envs/splitpipe-share
test -r "${SHARED_ENV}/csm_shared_environment.complete.json"
"${SHARED_ENV}/condabin/conda" run --prefix "${SHARED_ENV}" \
  python -c 'import edlib, pandas, pysam, scipy, anndata, splitpipe, yaml'
```

An existing Parse tree and Parse reference used by the pilot are:

```text
/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq/ParseBiosciences-Pipeline.1.6.4
/medpop/rgupta/Jiahao/reference_genome/parse_refdata-gex-GRCh38-2024-A
```

These are optional lab resources, not portable defaults. Configure your own
checkout, FASTQs, reference CSV, output and log paths. The shared environment's
publication manifest identifies its dependency snapshot; the pipeline code
comes from your v2.0 checkout. A separate kallisto executable/reference must
still be configured when selecting that backend.

### Shared kallisto companion (published 2026-10-22)

The optional tools prefix keeps a stable public name, separate from production Python:

```text
/medpop/rgupta/Jiahao/conda/envs/csm_kallisto-share
  -> csm_kallisto-share_20261022
```

Use the stable path in run configs. The symlink points to the validated dated
build; it can be switched atomically after a future build passes validation.
Keep dated builds referenced by existing runs. Record the resolved target and
completion manifest when reproducing an analysis. The production environment
name remains `splitpipe-share`.

It supplies kallisto **0.52.0**, bustools **0.45.1** and kb-python **0.29.5**.
It was validated with the existing `splitpipe-share` Python: 60 backend,
real-executable, UMI-correction and resource-interface checks passed. Package
inventories, binary/source hashes and a group-permission audit are recorded in
`csm_shared_kallisto.complete.json` inside the tools prefix. The production
publication manifest, Conda inventory and pip inventory checksums were unchanged.
A subsequent targeted resource/report/backend regression suite passed 144 tests
and 49 subtests; its test-only pytest runner was installed in job scratch.

Configure the executable literally in your run YAML:

```yaml
gex_alignment:
  backend: kallisto
  gene_assignment: unique
  umi_correction: parse
  output_level: gene
  kallisto:
    executable: /medpop/rgupta/Jiahao/conda/envs/csm_kallisto-share/bin/kallisto
    reference_dir: /medpop/rgupta/Jiahao/reference_genome/kallisto_refdata-gex-GRCh38-2024-A_nac_20261023
```

The reference above is the published lab reference described below. For building
another reference, use this prefix's `bin/kb` together with its explicit
`bin/kallisto`; kb-python also bundles an older kallisto, which is not the pinned
0.52.0 runtime. The repository builder forwards `--kallisto` explicitly.

Colleagues without access use the own-installation route. Lab administrators can
reproduce publication with `environment_kallisto.yml` and
`submit/submit_shared_kallisto_environment.sh`, selecting a **new dated prefix**
and a writable publication-record directory. `CSM_TOOLS_ALIAS` defaults to
`csm_kallisto-share` beside the dated build; the publisher updates only a symlink
and refuses to replace a directory or regular file at that public name. Submit it with the
`shared_kallisto_environment` resource profile. The worker refuses to overwrite
a sealed tools prefix, validates the production environment without modifying
it, and grants `rguptalab` read/execute access only in the new tools tree.
`CSM_VALIDATE_ONLY=true` resumes validation of a built, unsealed prefix; it is not
an in-place environment upgrade.

### Shared GRCh38 kallisto reference

The lab reference published on **2026-10-23** is:

```text
/medpop/rgupta/Jiahao/reference_genome/kallisto_refdata-gex-GRCh38-2024-A_nac_20261023
```

It is a fresh build with the shared kallisto **0.52.0** executable and
kb-python **0.29.5**, using workflow **`nac`** and Parse genome label **`hg38`**.
It contains **226,005 mature transcript targets**, **38,606 nascent targets**
and **38,606 genes**. It does not depend on a project output directory remaining
available. The earlier benchmark reference matches the source release but
failed index-loading checks with the shared Conda executable, so its index was
not reused. Both executables report 0.52.0; the version string alone did not
establish compatibility.

The matching Parse reference is
`/medpop/rgupta/Jiahao/reference_genome/parse_refdata-gex-GRCh38-2024-A`.
Its saved `process/mkref_def.json` and build log identify these exact inputs:

| Input | Shared path |
| --- | --- |
| FASTA | `/medpop/rgupta/Jiahao/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa` |
| GTF | `/medpop/rgupta/Jiahao/reference_genome/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz` |

This is the 10x **GRCh38 2024-A** reference, using its filtered **GENCODE v44**
annotation. The source and published index checksums are:

```text
FASTA SHA256  fb7421217e7058120cd60a5277445198e8deef2dff1edc46cd1e98b31fe64cbb
GTF SHA256    f589e2e8f31221c715548a84b1db12797624ce9f40c371d1544e5116962d5139
index.idx SHA256  6cf01de62ac69d29dbb73cdef567d20ef348faf28f2001bf8658533d157ac4b0
```

Validation checks the artifact and source checksums, matches all **194** source
chromosomes to the actual Parse STAR genome sequence, and compares the gene
universe, genome labels and GTF gene/exon coordinates with Parse, including
**483,352 retained exon-start entries**. It also reads the copied index with the
shared executable and maps **12** synthetic source exon/intron reads back to their
expected genes. `reference.json` records runtime-asset checksums;
`shared_reference.complete.json` records the source/Parse/tool identities,
validation results, build/publication jobs and ownership. The reference contains the
validation script, index-inspection log, executable-comparison logs and smoke-test
inputs/outputs.

The published tree is owned by **`jjiahao:rguptalab`**: directories are **0750**
and files **0640**. Lab members can read files and traverse directories; only
the owner can write. Group write permission is also removed from the parent
`reference_genome/`, preventing replacement or renaming of the protected tree.
Existing reference contents and benchmark outputs are unchanged. Future builds
should use a new destination rather than replace this dated reference.

An authorized rgupta user can use this complete dependency configuration in
their own run YAML:

```yaml
environment:
  conda_prefix: /medpop/rgupta/Jiahao/conda/envs/splitpipe-share
  conda_env: /medpop/rgupta/Jiahao/conda/envs/splitpipe-share
paths:
  splitpipe_dir: /medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq/ParseBiosciences-Pipeline.1.6.4
  genome_dir: /medpop/rgupta/Jiahao/reference_genome/parse_refdata-gex-GRCh38-2024-A
gex_alignment:
  backend: kallisto
  gene_assignment: unique       # or fractional
  umi_correction: parse         # or directional / none
  output_level: gene
  count_threads: 0
  kallisto:
    executable: /medpop/rgupta/Jiahao/conda/envs/csm_kallisto-share/bin/kallisto
    reference_dir: /medpop/rgupta/Jiahao/reference_genome/kallisto_refdata-gex-GRCh38-2024-A_nac_20261023
    strand: forward
```

Merge these fields into a copied project config; retain its other required
fields and set your own checkout, inputs, output and log directories. Pipeline
Python stays in `splitpipe-share`; the kallisto binary is invoked by absolute
path. No kallisto module or separate interactive activation is required. Users
with an existing production `splitpipe` installation may retain its environment
settings and use the same shared executable/reference paths.

The Broad `use .kallisto-0.45.0` and `use .kallisto-0.46.0` modules are older
than this validated stack. Kallisto 0.50 introduced an incompatible index format
and BUS strand-selection options used by the pipeline; those modules are not
substitutes for the 0.52.0 runtime and reference. See the
[upstream 0.50 release notes](https://github.com/pachterlab/kallisto/releases/tag/v0.50.0).

## Parse reference for both backends

`paths.genome_dir` must identify a complete Parse-compatible reference even
for kallisto. Parse initialization reads its gene annotations, genome labels,
kit/sample parameters and cell-calling context. A kallisto index cannot replace
that directory.

Reuse an accessible compatible reference, or build one from your chosen genome
FASTA and annotation GTF with the original Parse launcher in a compute job:

```bash
PROJECT_DIR="$(pwd -P)"
SPLITPIPE_DIR="${PROJECT_DIR}/ParseBiosciences-Pipeline.1.6.4"
FASTA="$(realpath references/source/genome.fa)"
GTF="$(realpath references/source/annotation.gtf)"
PARSE_REFERENCE="${PROJECT_DIR}/references/parse_GRCh38"

"${SPLITPIPE_DIR}/split-pipe" --mode mkref \
  --genome_name hg38 --fasta "${FASTA}" --genes "${GTF}" \
  --output_dir "${PARSE_REFERENCE}" --nthreads 8
```

Replace the source files and genome label with those for your experiment; the
label must be reused in the kallisto build. Retain the source release, source
files and build logs. Use the resulting reference directory as
`paths.genome_dir`. Reference construction is separate from the Steps 1–5
submission graph.

## Optional kallisto tools and reference

You can choose tools and references independently:

| Preference | Setup |
| --- | --- |
| Own kallisto, shared reference | Install the pinned tools below; set `kallisto.executable` to your absolute `bin/kallisto` and retain the shared `reference_dir`. |
| Shared kallisto, own reference | Build with the shared companion's `bin/kb` and explicit `bin/kallisto`; write the new reference in your own writable storage. |
| Own tools and reference | Use the installation and build recipe below, then configure both absolute paths. |

Before using a shared index with a different executable build, run
`"${KALLISTO}" inspect -t 8 "${KALLISTO_REFERENCE}/index.idx"` in a compute
allocation, with those variables set to your binary and the shared reference.
If it cannot read the index, build a personal reference with that executable.

When reusing the lab Parse reference, build from the exact FASTA/GTF paths and
`hg38` label listed above. A different assembly/annotation requires its matching
Parse reference too. A bare kallisto index or an arbitrary downloaded reference
does not supply this pipeline's required `reference.json`, `features.tsv` and
canonical gene mapping; use `scripts/build_kallisto_reference.py` to create the
complete contract. Shared reference files require no edits when changing
unique/fractional assignment or UMI correction.

Install the pinned optional tools separately, in your own writable prefix:

```bash
conda env create -f environment_kallisto.yml
conda run -n csm_kallisto kallisto version
conda run -n csm_kallisto python -c "from importlib.metadata import version; print(version('kb-python'))"
conda run -n csm_kallisto kb ref --help
```

The tracked environment pins kallisto 0.52.0, bustools 0.45.1 and kb-python
0.29.5. Production Python remains `splitpipe`; native counting uses the
repository's counter, not an unconditional `bustools umicorrect` invocation.
See the [kallisto release](https://github.com/pachterlab/kallisto/releases/tag/v0.52.0)
and [upstream reference-building documentation](https://www.kallistobus.tools/kb_usage/kb_ref/)
for the underlying tools; this repository's builder adds the compatibility
manifest and annotation checks required by Step 1a.

Resolve the executable paths inside the allocation:

```bash
KALLISTO="$(conda run -n csm_kallisto which kallisto)"
KB="$(conda run -n csm_kallisto which kb)"
KALLISTO_REFERENCE="${PROJECT_DIR}/references/kallisto_GRCh38"
conda activate splitpipe
python scripts/build_kallisto_reference.py \
  --fasta "${FASTA}" --gtf "${GTF}" --genome hg38 \
  --output "${KALLISTO_REFERENCE}" --kb "${KB}" \
  --kallisto "${KALLISTO}" --threads 8
```

Use the **same source FASTA/GTF release** as the Parse reference. By default the
builder uses a mature-plus-nascent (`nac`) reference, retaining intronic targets;
`--workflow standard` is an explicit mature-transcript-only alternative. Build
a new destination: the builder refuses an existing directory. Source hashes,
gene-ID alias reconciliation, tool versions and runtime-asset hashes are saved
in `reference.json`. Do not assemble or edit that manifest manually.

For a standalone UGER build after configuration, the equivalent worker is:

```bash
bash submit/submit_job.sh --config "${CONFIG}" --profile kallisto_reference -- \
  -v "CSM_CONFIG=${CONFIG},CSM_PROJECT_DIR=${PROJECT_DIR}" \
  -o "${PROJECT_DIR}/logs/kallisto-reference.\$JOB_ID.log" \
  "${PROJECT_DIR}/submit/submit_kallisto_reference.sh" \
  --fasta "${FASTA}" --gtf "${GTF}" --genome hg38 \
  --output "${KALLISTO_REFERENCE}" --kb "${KB}" --kallisto "${KALLISTO}"
```

Create the log directory first and use absolute paths passed to the worker.
Tune `uger.kallisto_reference` using [the resource guide](UGER_RESOURCES.md).
Configure the resolved kallisto binary and reference in your YAML afterward.
Wait for reference construction to finish before submitting Step 1.

## Separate analysis environment

`environment_qc.yml` creates `csm_analysis` for standalone Q1–Q4 notebooks and
analysis modules. It does not replace production `splitpipe`. The GEX carry-over
QC submission workflow uses the configured production environment because it
replays Parse preprocessing and controlled STAR processing. See
[QC overview](../QC_overview.md) and the [carry-over QC guide](GEX_GRNA_CONTAMINATION_QC.md).

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
  umi_correction: parse
  kallisto:
    executable: /medpop/rgupta/Jiahao/conda/envs/csm_kallisto-share/bin/kallisto
    reference_dir: /absolute/path/to/your/validated/kallisto_reference
    gene_assignment: unique
    output_level: gene
```

The reference path above is a placeholder to replace, not a supplied new lab
reference. Build one or reuse an explicitly validated, readable existing
reference. For building, use this prefix's `bin/kb` together with its explicit
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

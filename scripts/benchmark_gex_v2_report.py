"""All-pairs agreement, fixed-cohort UMI gains and stage costs for v2 Step 1a."""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

try:
    from scripts.benchmark_gex_v2 import ARMS, LIBRARIES, arm_options, read_json, validate_artifacts
    from scripts.benchmark_gex_common import load_product, correlation, sparse_cell_correlations, markdown_table
    from scripts.gex_alignment import atomic_json
except ModuleNotFoundError:
    from benchmark_gex_v2 import ARMS, LIBRARIES, arm_options, read_json, validate_artifacts
    from benchmark_gex_common import load_product, correlation, sparse_cell_correlations, markdown_table
    from gex_alignment import atomic_json

BASELINE = "star_unique_parse"


def assignment_cache_description(max_multimaps=3):
    return """## What the assignment cache contains

An assignment cache is reusable **read-to-gene annotation evidence between alignment and molecule counting**. It avoids repeating gene-annotation work for each UMI policy. The alignment cache is the original BAM/BUS output; the assignment cache adds gene information. Final corrected UMIs, molecule counts, fractional weights and cell calls are computed afterward.

Both STAR assignment caches read the **same complete, queryname-sorted BAM from one STAR alignment**. They retain different information:

| Aspect | STAR unique cache | STAR fractional cache |
| --- | --- | --- |
| Eligible alignments | Primary records passing native Parse's mapped-read/MAPQ gate; secondary alignments are excluded | Primary plus all retained secondary hits, requiring one primary per original read and an `NH` tag matching the complete hit set; the multimapping limit is MULTIMAP_LIMIT loci in this experiment |
| Gene annotation | Installed Parse `get_gene` selects one stranded, fully containing gene, or leaves the read unassigned. For overlapping genes, its ranking favors protein-coding status, gene-name preference and exonic overlap | Retains the union of all strand-compatible genes whose genomic spans contain each alignment, across every admitted locus. It preserves multiple candidates instead of selecting Parse's preferred overlapping gene; the unique-only MAPQ gate is not applied |
| UMI identities | Preserves native Parse UMI handling | Excludes empty or non-ACGT UMIs with diagnostics, matching the BUS/shared-counter representation |
| Saved representation | `cache/<library>/star_unique/assigned.bam`: primary reads with `GX` containing the chosen gene (empty if unassigned) and `RE` containing the binary exonic flag | `cache/<library>/star_fractional/reads.sqlite`: corrected barcode, observed UMI, candidate-class ID and read support; `evidence.pickle` maps those classes to gene sets and stores feature/barcode metadata and diagnostics |
| Decisions made later | Native Parse molecule grouping, the chosen UMI correction policy, count-matrix assembly and cell calling. A small adapter reads the saved `GX/RE` tags instead of rerunning `get_gene` | Exact-UMI evidence/collision resolution, refinement of candidate sets, the chosen UMI correction policy, and allocation of each retained molecule across its remaining candidate genes |

The fractional spool records each original read once, regardless of how many alignment hits it has. It rejects incomplete hit sets rather than treating a primary-only BAM as complete evidence. The unique cache has already removed alternative alignments and resolved overlapping-gene choices, so it cannot supply the missing evidence needed to reconstruct the fractional cache. STAR fractional changes both multimapping admission and overlapping-gene handling, including for reads with only one genomic alignment. It does not report native Parse exonic-fraction metadata.

For example, a high-confidence alignment contained in overlapping genes **A and B** may become `GX=A` in the unique cache, while the fractional cache retains **{A,B}**. If another read from the same cell and exact UMI supports **A** alone, the shared counter can resolve the compatible molecule to A. If a retained molecule remains ambiguous between A and B after resolution/correction, it contributes **0.5 to A and 0.5 to B**, conserving one whole molecule. A read with two admitted genomic alignments at A/B can also enter the fractional cache even when the native unique MAPQ gate excludes it.

Each STAR cache feeds the three `parse`, `directional` and `none` UMI policies. Kallisto instead caches BUS equivalence-class evidence and its gene annotations once, then reuses that evidence for all six assignment/UMI combinations. The measured `*.assign` phases include cache construction and serialization; `*.molecules` and `*.publish` measure the additional per-policy work. The serial STAR unique cache builder includes Parse's gene-ranking/exon calculations and BAM writing; the fractional builder uses gene-span compatibility sets and a SQLite spool. Their cache-construction times therefore describe different work and formats.

""".replace("MULTIMAP_LIMIT", str(max_multimaps))


def align_matrix(product, cells, genes):
    matrix, obs, features = product
    ids = pd.Index(features["gene_id"].astype(str))
    if set(ids) != set(genes):
        raise ValueError("Benchmark gene universes differ")
    rows = obs.index.get_indexer(cells)
    present = np.flatnonzero(rows >= 0)
    selected = matrix[rows[present]][:, ids.get_indexer(genes)].tocoo()
    return sparse.csr_matrix((selected.data, (present[selected.row], selected.col)), shape=(len(cells), len(genes)))


def relaxation_edges():
    edges = []
    for backend in ("star", "kallisto"):
        for umi in ("parse", "directional", "none"):
            edges.append((f"{backend}_unique_{umi}", f"{backend}_fractional_{umi}", "assignment"))
        for assignment in ("unique", "fractional"):
            for first, second in (("parse", "directional"), ("parse", "none"), ("directional", "none")):
                edges.append((f"{backend}_{assignment}_{first}", f"{backend}_{assignment}_{second}", "umi_policy"))
    return edges


def gain_metrics(first, second):
    first, second = np.asarray(first).ravel(), np.asarray(second).ravel()
    delta = second - first
    ratio = np.divide(second, first, out=np.full(len(first), np.nan), where=first > 0)
    return {"first_umis": float(first.sum()), "second_umis": float(second.sum()),
            "delta_umis": float(delta.sum()),
            "gain_percent": float(100 * delta.sum() / first.sum()) if first.sum() else float("nan"),
            "median_cell_delta": float(np.median(delta)), "median_cell_ratio": float(np.nanmedian(ratio)),
            "cells_gain": int((delta > 1e-6).sum()), "cells_lose": int((delta < -1e-6).sum()),
            "cells_equal": int((np.abs(delta) <= 1e-6).sum())}


def arm_stages(arm):
    b, a, _ = arm_options(arm)
    return (["preprocess", "star.align", "star.sort", f"star.{a}.assign"] if b == "star" else
            ["preprocess", "kallisto.format", "kallisto.align", "kallisto.assign"]) + [f"{arm}.molecules", f"{arm}.publish"]


def report(args):
    root = args.work_dir
    success = root / "BENCHMARK_SUCCESS.json"
    success.unlink(missing_ok=True)
    config = read_json(root / "experiment.json")["config"]
    output = root / "comparison"
    output.mkdir(exist_ok=True)
    summaries, agreements, gains, runtime, stages = [], [], [], [], []
    for library in LIBRARIES:
        sample = next(lib["sample_name"] for lib in config["libraries"] if lib["sublibrary_name"] == library)
        required = {stage for arm in ARMS for stage in arm_stages(arm)} | {"sample"}
        records = {}
        for stage in sorted(required):
            record = read_json(root / "manifests" / f"{library}.{stage}.json")
            validate_artifacts(record, full=args.verify_cache)
            for dependency, digest in record["dependencies"].items():
                try:
                    from scripts.gex_alignment import sha256_file
                except ModuleNotFoundError:
                    from gex_alignment import sha256_file
                if sha256_file(root / "manifests" / f"{library}.{dependency}.json") != digest:
                    raise ValueError(f"Stale dependency: {library}/{stage}")
            records[stage] = record
            metrics = record["metrics"]
            stages.append({"library": library, "stage": stage, **{k: v for k, v in metrics.items() if k not in ("command", "cpu_affinity")},
                           "output_bytes": sum(a["bytes"] for a in record["artifacts"])})
        called, allwell = {}, {}
        for arm in ARMS:
            source = root / arm / f"gex_{library}"
            called[arm] = load_product(source / sample / "DGE_filtered")
            allwell[arm] = load_product(source / "process", process=True)
        cohort = called[BASELINE][1].index
        genes = pd.Index(called[BASELINE][2]["gene_id"].astype(str))
        fixed = {arm: align_matrix(allwell[arm], cohort, genes) for arm in ARMS}
        totals = {arm: np.asarray(matrix.sum(axis=1)).ravel() for arm, matrix in fixed.items()}
        pd.DataFrame(totals, index=cohort).to_csv(output / f"{library}.fixed_cohort_umis.csv.gz")
        gene_rows = []
        for arm in ARMS:
            b, a, u = arm_options(arm)
            matrix, obs, _ = called[arm]
            observed = np.asarray(matrix.sum(axis=1)).ravel()
            summaries.append({"library": library, "arm": arm, "backend": b, "assignment": a, "umi_correction": u,
                              "called_cells": len(obs), "called_total_umis": float(observed.sum()),
                              "called_median_umis": float(np.median(observed)),
                              "fixed_cells": len(cohort), "fixed_total_umis": float(totals[arm].sum()),
                              "fixed_median_genes": float(np.median(np.asarray((fixed[arm] > 0).sum(axis=1)).ravel())),
                              "gain_vs_baseline_percent": gain_metrics(totals[BASELINE], totals[arm])["gain_percent"]})
            selected = arm_stages(arm)
            incremental = selected[-2:]
            runtime.append({"library": library, "arm": arm,
                            "reconstructed_step1a_seconds": sum(records[s]["metrics"]["wall_seconds"] for s in selected),
                            "cached_policy_seconds": sum(records[s]["metrics"]["wall_seconds"] for s in incremental),
                            "reconstructed_cpu_hours": sum(records[s]["metrics"]["cpu_hours"] for s in selected),
                            "reconstructed_reserved_core_hours": sum(records[s]["metrics"]["reserved_core_hours"] for s in selected),
                            "max_stage_process_rss_gb": max(records[s]["metrics"]["max_process_rss_kb"] / 1024**2 for s in selected)})
            # Baseline gene-level differences expose gains concentrated in a few features.
            before = np.asarray(fixed[BASELINE].sum(axis=0)).ravel()
            after = np.asarray(fixed[arm].sum(axis=0)).ravel()
            gene_rows.append(pd.DataFrame({"arm": arm, "gene_id": genes, "baseline_umis": before,
                                          "umis": after, "delta_umis": after - before}))
        pd.concat(gene_rows, ignore_index=True).to_csv(output / f"{library}.gene_gains.csv.gz", index=False)
        for first, second in itertools.combinations(ARMS, 2):
            shared = called[first][1].index.intersection(called[second][1].index)
            union = called[first][1].index.union(called[second][1].index)
            x, y = fixed[first], fixed[second]
            common_x = align_matrix(allwell[first], shared, genes)
            common_y = align_matrix(allwell[second], shared, genes)
            cell_corr = sparse_cell_correlations(x, y)
            agreements.append({"library": library, "first": first, "second": second,
                "fixed_cells": len(cohort), "shared_called_cells": len(shared), "union_called_cells": len(union),
                "cell_call_jaccard": len(shared) / len(union) if len(union) else np.nan,
                "fixed_umi_pearson": correlation(totals[first], totals[second]),
                "fixed_umi_spearman": correlation(totals[first], totals[second], "spearman"),
                "fixed_median_cell_expression_log_pearson": float(np.nanmedian(cell_corr)),
                "fixed_pseudobulk_log_pearson": correlation(np.log1p(np.asarray(x.sum(axis=0)).ravel()), np.log1p(np.asarray(y.sum(axis=0)).ravel())),
                "common_called_median_cell_expression_log_pearson": float(np.nanmedian(sparse_cell_correlations(common_x, common_y))) if len(shared) else np.nan})
        for first, second, factor in relaxation_edges():
            gains.append({"library": library, "first": first, "second": second, "factor": factor,
                          "cohort": BASELINE, "fixed_cells": len(cohort), **gain_metrics(totals[first], totals[second])})
    frames = {"library_summary": pd.DataFrame(summaries), "agreement": pd.DataFrame(agreements),
              "umi_gains": pd.DataFrame(gains), "runtime": pd.DataFrame(runtime), "stage_resources": pd.DataFrame(stages)}
    for name, frame in frames.items():
        frame.to_csv(output / f"{name}.csv", index=False)
    failed = []
    for path in sorted((root / "failed_attempts").glob("*/*/*.metrics.json")):
        metrics = read_json(path)
        failed.append({"library": path.parents[1].name, "attempt": path.parent.name,
                       **{k: v for k, v in metrics.items() if k not in ("command", "cpu_affinity")}})
    failed = pd.DataFrame(failed, columns=None if failed else ["library", "attempt", "wall_seconds", "cpu_hours", "reserved_core_hours"])
    failed.to_csv(output / "failed_attempts.csv", index=False)
    retired, equality = [], []
    baseline = root / "timing_baselines/shared_sqlite_v1"
    for path in sorted((baseline / "manifests").glob("*.json")):
        record = read_json(path)
        retired.append({"library": record["library"], "stage": record["stage"], "storage": "shared_sqlite_v1",
                        **{k: v for k, v in record["metrics"].items() if k not in ("command", "cpu_affinity")}})
        if record["stage"].endswith(".molecules"):
            arm = record["stage"].removesuffix(".molecules")
            proof = read_json(root / arm / f"gex_{record['library']}" / "process/count_staging_equality.json")
            if not all(proof.get(key) is True for key in ("complete", "exact_matrices_and_metadata", "identical_molecule_assignments")):
                raise ValueError("Local cache staging equality verification is missing or failed")
            equality.append(proof)
    retired = pd.DataFrame(retired, columns=None if retired else ["library", "stage", "storage", "wall_seconds", "cpu_hours", "reserved_core_hours"])
    retired.to_csv(output / "cache_storage_baseline.csv", index=False)
    pd.DataFrame(equality, columns=None if equality else ["library", "arm", "complete", "exact_matrices_and_metadata", "identical_molecule_assignments"]).to_csv(output / "count_staging_equality.csv", index=False)
    plot(frames, output)
    summary = frames["library_summary"].groupby("arm", sort=False).agg(
        called_cells=("called_cells", "sum"), fixed_umis=("fixed_total_umis", "sum"),
        median_gain_percent=("gain_vs_baseline_percent", "median")).reset_index()
    times = frames["runtime"].groupby("arm", sort=False).median(numeric_only=True).reset_index()
    summary = summary.merge(times[["arm", "reconstructed_step1a_seconds", "cached_policy_seconds"]], on="arm")
    actual = frames["stage_resources"]
    step1a = actual[actual.stage != "sample"]
    stage_summary = actual.groupby("stage", sort=True).agg(
        median_wall_seconds=("wall_seconds", "median"), cpu_hours=("cpu_hours", "sum"),
        reserved_core_hours=("reserved_core_hours", "sum"), peak_process_rss_kb=("max_process_rss_kb", "max")).reset_index()
    stage_summary["peak_process_rss_gb"] = stage_summary.pop("peak_process_rss_kb") / 1024**2
    stage_summary.to_csv(output / "stage_summary.csv", index=False)
    baseline_agreement = frames["agreement"][frames["agreement"]["first"] == BASELINE].groupby("second").agg(
        median_cell_call_jaccard=("cell_call_jaccard", "median"),
        median_fixed_umi_spearman=("fixed_umi_spearman", "median"),
        median_cell_expression_log_pearson=("fixed_median_cell_expression_log_pearson", "median"),
        minimum_library_expression_agreement=("fixed_median_cell_expression_log_pearson", "min")).reset_index()
    gains_summary = frames["umi_gains"].groupby(["factor", "first", "second"], sort=False).agg(
        total_delta_umis=("delta_umis", "sum"), median_library_gain_percent=("gain_percent", "median"),
        min_library_gain_percent=("gain_percent", "min"), max_library_gain_percent=("gain_percent", "max")).reset_index()
    failed_wall = failed.wall_seconds.sum() / 3600 if len(failed) else 0
    failed_cpu = failed.cpu_hours.sum() if len(failed) else 0
    failed_reserved = failed.reserved_core_hours.sum() if len(failed) else 0
    retired_wall = retired.wall_seconds.sum() / 3600 if len(retired) else 0
    retired_cpu = retired.cpu_hours.sum() if len(retired) else 0
    retired_reserved = retired.reserved_core_hours.sum() if len(retired) else 0
    text = f"""# v2 Step 1a benchmark: {args.pairs:,} GEX read pairs per sublibrary

Completed C1/C2/T1/T2 comparison of all twelve configurations. Each library uses one fixed Parse 1.6.4 barcode-extraction/decontamination result, one STAR alignment and one kallisto alignment. Gene assignment caches feed the three UMI policies. All {len(agreements)} library/configuration-pair comparisons use the same STAR unique + Parse called-cell cohort within each library. Fractional gene weights are retained.

{markdown_table(summary, list(summary.columns))}

![Fixed-cohort gains](umi_gains.png)
![Expression agreement](agreement.png)
![Reconstructed stage costs](runtime.png)

Successful Step 1a work (shared stages counted once): {step1a.wall_seconds.sum()/3600:.3f} summed task-hours, {step1a.cpu_hours.sum():.3f} CPU-hours, {step1a.reserved_core_hours.sum():.3f} reserved core-hours. Recorded failed command attempts add {failed_wall:.3f} task-hours, {failed_cpu:.3f} CPU-hours and {failed_reserved:.3f} reserved core-hours. These summed task-hours are not elapsed makespan. Sampling and scheduler/setup overhead are separate; the scheduler accounting files retain complete job costs, including dependency failures before a command starts.

The initial shared-filesystem SQLite counting measurements were retained as a storage baseline and replaced by consistent job-local staging timings. Retired successful command measurements add {retired_wall:.3f} task-hours, {retired_cpu:.3f} CPU-hours and {retired_reserved:.3f} reserved core-hours to actual experiment costs; they are excluded from the final per-configuration timings. All {len(equality)} completed baseline molecule products match the retimed products exactly, including matrices, ordered metadata, scientific counters and molecule-assignment bytes. Local cache-copy time and disk space are included in the final molecule phases. The STAR unique annotation cache uses one streaming process; ordinary Parse production annotates inside parallel molecule workers, so reconstructed totals describe this cache experiment rather than a direct production throughput test.

{assignment_cache_description(config.get('gex_alignment', {}).get('star', {}).get('max_multimaps', 3))}## Measured stages

Wall time is the median of four libraries; CPU/core-hours sum the four libraries, and RSS is the maximum. Sampling is shown separately from Step 1a.

{markdown_table(stage_summary, list(stage_summary.columns))}

## Agreement with STAR unique + Parse

These are library medians on the same baseline-called cells. The complete table contains every one of the 66 configuration pairs per library.

{markdown_table(baseline_agreement, list(baseline_agreement.columns))}

## Changes from assignment and UMI settings

Absolute deltas sum the four fixed cohorts; percentage columns summarize library-specific changes. Parse → directional is an alternative graph comparison, not a guaranteed relaxation.

{markdown_table(gains_summary, list(gains_summary.columns))}

Runtime paths are reconstructed sums of measured stages, including cache I/O, process startup and reference initialization. They are not independent end-to-end runs or cold-cache aligner speed tests. cached_policy_seconds measures only molecule correction/matrix assembly and cell calling/publication with existing assignment caches. wait4 max_process_rss_kb is the largest individual process RSS high-water mark. peak_tree_rss_gb is sampled every two seconds and sums worker RSS, counting shared pages more than once. Sampled peak_temp_bytes and peak_working_file_bytes describe apparent file sizes; the latter includes linked input caches and is not extra physical disk allocation. Resource requests are recorded; qacct maxvmem is virtual memory. Scheduler accounting should be retained alongside these measurements.

[Per-stage resources](stage_resources.csv), [stage summaries](stage_summary.csv), [per-configuration runtime](runtime.csv), [yield](library_summary.csv), [all 66 configuration pairs per library](agreement.csv), [factor-matched signed UMI changes](umi_gains.csv), [failed command attempts](failed_attempts.csv), [retired storage measurements](cache_storage_baseline.csv), [exact staging equality](count_staging_equality.csv).

UMI gains are changes in counted molecules, not proven biological recovery. Unique-to-fractional changes include assignment/collision-policy effects. Parse and directional correction are alternative graphs with no guaranteed yield ordering; negative deltas are retained. No correction still deduplicates exact UMIs. Called-cell totals and fixed-cohort totals are reported separately. One observation per stage/library provides no timing confidence interval; libraries are distinct datasets, not technical replicates. Reference construction, sampling and Step 1b are outside reconstructed Step 1a time.
"""
    (output / "REPORT.md").write_text(text)
    atomic_json(success, {"complete": True, "pairs_per_library": args.pairs, "libraries": list(LIBRARIES), "arms": list(ARMS),
                          "configuration_pairs_per_library": 66, "baseline_cohort": BASELINE,
                          "full_cache_checksums_verified": args.verify_cache})


def plot(frames, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    summary = frames["library_summary"]
    fig, ax = plt.subplots(figsize=(13, 6), constrained_layout=True)
    summary.pivot(index="arm", columns="library", values="gain_vs_baseline_percent").reindex(ARMS).plot.bar(ax=ax)
    ax.set(ylabel="UMI change vs STAR unique + Parse (%)", xlabel="", title="Same baseline called cells in every configuration")
    fig.savefig(output / "umi_gains.png", dpi=180)
    plt.close(fig)
    agreement = frames["agreement"]
    grid = np.eye(len(ARMS))
    for (first, second), group in agreement.groupby(["first", "second"]):
        i, j = ARMS.index(first), ARMS.index(second)
        grid[i, j] = grid[j, i] = group.fixed_median_cell_expression_log_pearson.median()
    fig, ax = plt.subplots(figsize=(12, 11), constrained_layout=True)
    view = ax.imshow(grid, vmin=0, vmax=1, cmap="viridis")
    ax.set_xticks(range(len(ARMS)), ARMS, rotation=90)
    ax.set_yticks(range(len(ARMS)), ARMS)
    ax.set_title("Median across libraries: fixed-cohort per-cell log-expression Pearson")
    fig.colorbar(view, ax=ax)
    fig.savefig(output / "agreement.png", dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(12, 7), constrained_layout=True)
    rows = []
    metrics = frames["stage_resources"].set_index(["library", "stage"])
    for arm in ARMS:
        components = arm_stages(arm)
        rows.append([np.median([metrics.loc[(library, s), "wall_seconds"] for library in LIBRARIES]) / 60 for s in components])
    # Both paths have six stages: pre, align/format, sort/align, assign, molecules, publish.
    data = pd.DataFrame(rows, index=ARMS, columns=["Parse pre", "STAR align / kallisto format", "STAR sort / kallisto align", "assignment cache", "UMI + matrix", "cell calls + output"])
    data.plot.barh(stacked=True, ax=ax)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2)
    ax.set(xlabel="Sum of library-median stage minutes", ylabel="", title="Reconstructed Step 1a cost; cache I/O included")
    fig.savefig(output / "runtime.png", dpi=180)
    plt.close(fig)

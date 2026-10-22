"""Gene/cell agreement and measured runtime for the matched pilot benchmark."""

import datetime
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as sio, sparse, stats

try:
    from scripts.gex_alignment import atomic_json
    from scripts.benchmark_gex_backends import ARMS, BASE_ARMS, UMI_FLAVORS, arm_options
except ModuleNotFoundError:
    from gex_alignment import atomic_json
    from benchmark_gex_backends import ARMS, BASE_ARMS, UMI_FLAVORS, arm_options


def correlation(first, second, method="pearson"):
    first, second = np.asarray(first), np.asarray(second)
    finite = np.isfinite(first) & np.isfinite(second)
    first, second = first[finite], second[finite]
    if len(first) < 3 or np.std(first) == 0 or np.std(second) == 0:
        return float("nan")
    return float(stats.spearmanr(first, second).statistic if method == "spearman"
                 else np.corrcoef(first, second)[0, 1])


def sparse_cell_correlations(first, second):
    """Per-cell Pearson across all genes, including shared implicit zeros."""
    first, second = first.copy().astype(float).tocsr(), second.copy().astype(float).tocsr()
    first.data = np.log1p(first.data)
    second.data = np.log1p(second.data)
    n = first.shape[1]
    sum1, sum2 = np.asarray(first.sum(axis=1)).ravel(), np.asarray(second.sum(axis=1)).ravel()
    cross = np.asarray(first.multiply(second).sum(axis=1)).ravel() - sum1 * sum2 / n
    var1 = np.asarray(first.multiply(first).sum(axis=1)).ravel() - sum1**2 / n
    var2 = np.asarray(second.multiply(second).sum(axis=1)).ravel() - sum2**2 / n
    denominator = np.sqrt(np.maximum(var1 * var2, 0))
    return np.divide(cross, denominator, out=np.full(len(cross), np.nan), where=denominator > 0)


def load_product(path, process=False):
    if process:
        matrix_name, obs_name = "aw_uf_count_matrix.mtx", "aw_uf_cell_metadata.csv"
    else:
        matrix_name, obs_name = "count_matrix.mtx", "cell_metadata.csv"
    matrix = sparse.csr_matrix(sio.mmread(path / matrix_name))
    obs = pd.read_csv(path / obs_name, index_col="bc_wells")
    genes = pd.read_csv(path / "all_genes.csv", dtype=str)
    obs.index = obs.index.astype(str)
    if matrix.shape != (len(obs), len(genes)) or not obs.index.is_unique or genes["gene_id"].duplicated().any():
        raise ValueError(f"Malformed benchmark product: {path}")
    return matrix, obs, genes


def star_alignment_seconds(path):
    times = {}
    for line in path.read_text().splitlines():
        if "|" not in line:
            continue
        name, value = (part.strip() for part in line.split("|", 1))
        if name in ("Started job on", "Started mapping on", "Finished on"):
            times[name] = datetime.datetime.strptime(value, "%b %d %H:%M:%S")
    if not {"Started job on", "Finished on"} <= set(times):
        return float("nan")
    elapsed = (times["Finished on"] - times["Started job on"]).total_seconds()
    return elapsed if elapsed >= 0 else elapsed + 86400


def markdown_table(frame, columns):
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in frame[columns].itertuples(index=False, name=None):
        lines.append("| " + " | ".join(f"{value:.3g}" if isinstance(value, (float, np.floating)) else str(value)
                                       for value in row) + " |")
    return "\n".join(lines)


def plot_runtime(runtimes, output_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 2, figsize=(13, 14), constrained_layout=True)
    for row, flavor in enumerate(UMI_FLAVORS):
        subset = runtimes[runtimes["umi_correction"] == flavor]
        for axis, metric, title in zip(axes[row], ("alignment_seconds", "core_seconds"),
                                        ("Alignment (includes index loading)", "Alignment + counting + DGE")):
            subset.pivot(index="library", columns="backend", values=metric).reindex(columns=BASE_ARMS).plot.bar(ax=axis)
            axis.set(title=f"{title}; UMI correction: {flavor}", ylabel="Seconds", xlabel="Library")
            axis.set_ylim(0, float(subset[metric].max()) * 1.4)
            axis.tick_params(axis="x", rotation=0)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def report(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import yaml

    config = yaml.safe_load(Path(args.config).read_text())
    root = args.work_dir
    arms = ARMS
    # Exercise the real Step 1b handoff as an integration check. Its validation
    # work is separate from the matched per-library Step 1a runtime comparison.
    try:
        from scripts.combine_kallisto_gex import combine
        from scripts.gex_alignment import validate_combined, validate_star_policies
    except ModuleNotFoundError:
        from combine_kallisto_gex import combine
        from gex_alignment import validate_combined, validate_star_policies
    for arm in arms:
        config_path = root / "configs" / f"{arm}.yaml"
        arm_config = yaml.safe_load(config_path.read_text())
        if arm_options(arm)[0] == "star":
            validate_star_policies(arm_config, config_path)
            continue
        combine(arm_config, config_path)
        validate_combined(arm_config, config_path)
    output = root / "comparison"
    output.mkdir(exist_ok=True)
    summary, agreement, runtimes = [], [], []
    figures = {flavor: plt.subplots(4, 2, figsize=(10, 16), constrained_layout=True)
               for flavor in UMI_FLAVORS}
    comparisons = []
    for flavor in UMI_FLAVORS:
        group = [base if flavor == "parse" else f"{base}_{flavor}" for base in BASE_ARMS]
        comparisons.extend(((group[0], group[1]), (group[0], group[2]), (group[1], group[2])))
    for base in BASE_ARMS:
        comparisons.extend((base, f"{base}_{flavor}") for flavor in UMI_FLAVORS[1:])
    for lib_index, library in enumerate(config["libraries"]):
        name, sample = library["sublibrary_name"], library["sample_name"]
        products, awuf, totals, manifests = {}, {}, {}, {}
        pseudoaligned_counts = []
        for arm in arms:
            base, flavor = arm_options(arm)
            completion = root / "manifests" / f"{name}.{arm}.json"
            manifests[arm] = json.loads(completion.read_text())
            if manifests[arm].get("complete") is not True:
                raise ValueError(f"Incomplete benchmark arm: {name}/{arm}")
            arm_root = root / arm / f"gex_{name}"
            products[arm] = load_product(arm_root / sample / "DGE_filtered")
            awuf[arm] = load_product(arm_root / "process", process=True)
            matrix, obs, _ = products[arm]
            values = np.asarray(matrix.sum(axis=1)).ravel()
            gene_counts = np.asarray((matrix > 0).sum(axis=1)).ravel()
            totals[arm] = pd.Series(values, index=obs.index)
            summary.append({"library": name, "arm": arm, "backend": base, "umi_correction": flavor,
                            "called_cells": len(obs),
                            "total_umis": float(values.sum()), "median_umis_per_cell": float(np.median(values)),
                            "mean_umis_per_cell": float(np.mean(values)), "median_genes_per_cell": float(np.median(gene_counts))})
            timing = manifests[arm]
            native_qc = (json.loads((arm_root / "process/kallisto_qc.json").read_text())
                         if base != "star" else {})
            native_timing = native_qc.get("timing", {})
            if native_qc:
                pseudoaligned_counts.append(native_qc["pseudoaligned_reads"])
            alignment_seconds = (star_alignment_seconds(arm_root / "process/barcode_headLog.final.out") if base == "star"
                                 else native_timing["alignment_seconds"])
            runtimes.append({"library": name, "arm": arm, "threads": timing["threads"],
                             "backend": base, "umi_correction": flavor,
                             "alignment_seconds": alignment_seconds, "core_seconds": timing["elapsed_seconds"],
                             "input_format_seconds": native_timing.get("input_format_seconds", float("nan")),
                             "umi_count_seconds": native_timing.get("umi_count_seconds", float("nan")),
                             "publish_seconds": native_timing.get("publish_seconds", float("nan")),
                             "metadata_republish_seconds": timing.get("metadata_roundoff_repair", {}).get("elapsed_seconds", 0),
                             "shared_preprocess_seconds": timing["shared_preprocess_seconds"],
                             "preprocess_plus_core_seconds": timing["elapsed_seconds"] + timing["shared_preprocess_seconds"],
                             "peak_process_rss_gb": timing.get("peak_rss_gb", float("nan"))})
        if len({row["survivor_sha256"] for row in manifests.values()}) != 1:
            raise ValueError("Backend input mismatch in benchmark")
        if len(set(pseudoaligned_counts)) != 1:
            raise ValueError("Kallisto mapping totals changed across assignment/correction modes")
        detailed = pd.concat(totals, axis=1).fillna(0)
        detailed.to_csv(output / f"{name}.called_union_umis.csv.gz")
        for first_arm, second_arm in comparisons:
            first, second = products[first_arm], products[second_arm]
            shared = first[1].index.intersection(second[1].index)
            union = first[1].index.union(second[1].index)
            gene_ids = pd.Index(first[2]["gene_id"])
            second_ids = pd.Index(second[2]["gene_id"])
            if set(gene_ids) != set(second_ids):
                raise ValueError("Benchmark gene universes differ")
            first_matrix = first[0][first[1].index.get_indexer(shared)]
            second_matrix = second[0][second[1].index.get_indexer(shared)][:, second_ids.get_indexer(gene_ids)]
            first_umis = np.asarray(first_matrix.sum(axis=1)).ravel()
            second_umis = np.asarray(second_matrix.sum(axis=1)).ravel()
            cell_correlations = sparse_cell_correlations(first_matrix, second_matrix)
            baseline_cohort = first[1].index
            second_awuf_matrix, second_awuf_obs, _ = awuf[second_arm]
            second_awuf_totals = pd.Series(np.asarray(second_awuf_matrix.sum(axis=1)).ravel(), index=second_awuf_obs.index)
            matched_totals = second_awuf_totals.reindex(baseline_cohort, fill_value=0)
            baseline_totals = totals[first_arm]
            ratio = matched_totals / baseline_totals
            agreement.append({
                "library": name, "first": first_arm, "second": second_arm,
                "shared_called_cells": len(shared), "union_called_cells": len(union),
                "cell_call_jaccard": len(shared) / len(union) if len(union) else float("nan"),
                "first_called_retention": len(shared) / len(first[1]) if len(first[1]) else float("nan"),
                "shared_umi_pearson": correlation(first_umis, second_umis),
                "shared_umi_spearman": correlation(first_umis, second_umis, "spearman"),
                "fixed_first_cohort_median_umi_ratio": float(np.median(ratio)),
                "fixed_first_cohort_log_umi_pearson": correlation(np.log1p(baseline_totals), np.log1p(matched_totals)),
                "median_cell_expression_log_pearson": float(np.nanmedian(cell_correlations)),
                "common_cell_pseudobulk_log_pearson": correlation(
                    np.log1p(np.asarray(first_matrix.sum(axis=0)).ravel()),
                    np.log1p(np.asarray(second_matrix.sum(axis=0)).ravel())),
            })
            first_base, first_flavor = arm_options(first_arm)
            second_base, second_flavor = arm_options(second_arm)
            if first_base == "star" and second_base != "star" and first_flavor == second_flavor:
                axis = figures[first_flavor][1][lib_index, BASE_ARMS.index(second_base) - 1]
                x, y = np.log1p(baseline_totals), np.log1p(matched_totals)
                axis.hexbin(x, y, gridsize=55, mincnt=1, bins="log", cmap="viridis")
                bound = max(x.max(), y.max())
                axis.plot([0, bound], [0, bound], color="tomato", linewidth=1)
                axis.set(title=f"{name}: {second_base}, {first_flavor}", xlabel="STAR log(1 + UMI)", ylabel="kallisto log(1 + UMI)")
    for flavor, (fig, _) in figures.items():
        suffix = "" if flavor == "parse" else "_" + flavor
        fig.savefig(output / f"umi_agreement{suffix}.png", dpi=180)
        plt.close(fig)
    summary, agreement, runtimes = pd.DataFrame(summary), pd.DataFrame(agreement), pd.DataFrame(runtimes)
    for name, frame in (("library_summary", summary), ("agreement", agreement), ("runtime", runtimes)):
        frame.to_csv(output / f"{name}.csv", index=False)
    baseline_runtime = runtimes[runtimes["arm"] == "star"].set_index("library")
    runtimes["core_speedup_vs_star"] = [baseline_runtime.loc[row.library, "core_seconds"] / row.core_seconds
                                         for row in runtimes.itertuples()]
    runtimes["alignment_speedup_vs_star"] = [baseline_runtime.loc[row.library, "alignment_seconds"] / row.alignment_seconds
                                              for row in runtimes.itertuples()]
    same_flavor_star = runtimes[runtimes["backend"] == "star"].set_index(["library", "umi_correction"])
    runtimes["core_speedup_vs_same_flavor_star"] = [
        same_flavor_star.loc[(row.library, row.umi_correction), "core_seconds"] / row.core_seconds
        for row in runtimes.itertuples()]
    runtimes.to_csv(output / "runtime.csv", index=False)
    plot_runtime(runtimes, output / "runtime.png")
    text = "# Four-library kallisto benchmark: 10 million paired reads per library\n\n"
    text += ("> For current installation and counting/resource options, see the [v2.0 quick start](QUICKSTART.md), "
             "[backend guide](KALLISTO_BACKEND.md) and [resource guide](UGER_RESOURCES.md).\n\n")
    text += ("STAR uses the existing V1.5 alignment, molecule and DGE implementations. "
             "Kallisto uses unique and equal-fractional gene assignment. All three paths are run with "
             "Parse-default, directional, and no neighboring-UMI correction (nine combinations). "
             "All arms receive byte-identical Parse-preprocessed reads and use the same gene annotation. "
             "Each arm has 16 allocated CPU slots. Sampling selects balanced deterministic prefixes across source FASTQ pairs.\n\n"
             "STAR retains its standard BAM intermediate route; kallisto uses native BUS intermediates in job-local scratch. "
             "These timings compare the implemented pipeline paths, including index loading and their respective output work.\n\n"
             "Kallisto uses a mature-plus-nascent reference to retain intronic sensitivity. "
             "Runtime values describe the counter implementation recorded for each run; the earlier serial-counter "
             "baseline and the v2.0 parallel-counter rerun are separate experiments. "
             "The measured core timings therefore describe these implementations; they do not isolate aligner throughput. "
             "The runtime CSV separately records kallisto input formatting, alignment, molecule counting and publishing.\n\n"
             "The timed core ends at DGE output for all arms: optional clustering and HTML reporting are excluded. "
             "Preprocessing is measured once per library and added separately to each core time. "
             "Each timed core independently runs alignment and counting; completed Parse-default runs are reused as the baseline. "
             "Reference construction, read sampling and Step 1b combination are excluded. All six kallisto combinations passed the actual "
             "Step 1b combine and completion checks before this report was generated. "
             "One run per arm provides observed wall times, not uncertainty estimates; "
             "node scheduling and filesystem cache can affect them. The original Parse-default arrays allowed two concurrent jobs each; "
             "the correction-mode extension allowed four. Each individual run retained 16 CPU slots and 8 GB per slot. "
             "RSS is the child-process high-water mark reported by getrusage, not summed process-tree memory.\n\n")
    overview = summary.groupby("arm").agg(
        called_cells_total=("called_cells", "sum"),
        median_library_umis_per_cell=("median_umis_per_cell", "median"))
    overview["median_core_seconds"] = runtimes.groupby("arm")["core_seconds"].median()
    overview = overview.reindex(arms).reset_index()
    text += ("## Overview across four libraries\n\n"
             "Cell counts are summed across libraries; UMI and runtime summaries are medians of the four library values.\n\n")
    text += markdown_table(overview, ["arm", "called_cells_total", "median_library_umis_per_cell", "median_core_seconds"]) + "\n\n"
    text += "## Cell and molecule yield\n\n" + markdown_table(summary, ["library", "arm", "called_cells", "total_umis", "median_umis_per_cell", "median_genes_per_cell"]) + "\n\n"
    text += "## Agreement\n\n"
    text += ("UMI correlations use the common called-cell cohort. The fixed-cohort UMI ratio uses every cell called by the first arm, "
             "retrieving the other arm's counts from its all-well matrix and filling missing barcodes with zero. "
             "Expression correlations use log(1 + UMI) across the same genes; per-cell correlations include shared zeros.\n\n")
    text += markdown_table(agreement, ["library", "first", "second", "cell_call_jaccard", "shared_umi_spearman",
                                       "fixed_first_cohort_median_umi_ratio", "median_cell_expression_log_pearson",
                                       "common_cell_pseudobulk_log_pearson"]) + "\n\n"
    for flavor in UMI_FLAVORS:
        suffix = "" if flavor == "parse" else "_" + flavor
        text += f"![UMIs in the fixed STAR-called cohort, {flavor} correction](umi_agreement{suffix}.png)\n\n"
    text += "## Runtime\n\n" + markdown_table(runtimes, ["library", "arm", "alignment_seconds", "core_seconds",
                                                              "preprocess_plus_core_seconds", "core_speedup_vs_same_flavor_star", "peak_process_rss_gb"]) + "\n\n"
    text += "![Measured wall times](runtime.png)\n\n"
    text += ("Gene resolution and cell-call policies can change counts; these measurements do not establish biological equivalence. "
             "Transcript output is disabled in this gene-level comparison. Full-depth scaling remains unmeasured.\n")
    if (runtimes["metadata_republish_seconds"] > 0).any():
        text += ("\nThe completed Parse/default and directional fractional products were republished to correct floating-point "
                 "roundoff at whole-molecule cell-calling cutoffs. Gene weights and UMI families were preserved; omitted wells "
                 "at the all-well minimum were restored from the retained molecule assignments. Original metadata was archived. "
                 "Core runtimes retain the original alignment/counting measurements; the separate republish pass is recorded as "
                 "`metadata_republish_seconds` in the runtime CSV and excluded from core timings. The no-correction fractional "
                 "runs used whole-molecule calling totals directly.\n")
    (output / "REPORT.md").write_text(text)
    if getattr(args, "publish_report", True):
        docs = Path(__file__).resolve().parents[1] / "docs"
        assets = docs / "assets/kallisto_benchmark_10m"
        assets.mkdir(parents=True, exist_ok=True)
        for name in ("library_summary.csv", "agreement.csv", "runtime.csv", "umi_agreement.png",
                     "umi_agreement_directional.png", "umi_agreement_none.png", "runtime.png"):
            shutil.copy2(output / name, assets / name)
        public = text.replace("](umi_agreement.png)", "](assets/kallisto_benchmark_10m/umi_agreement.png)")
        public = public.replace("](runtime.png)", "](assets/kallisto_benchmark_10m/runtime.png)")
        for flavor in UMI_FLAVORS[1:]:
            public = public.replace(f"](umi_agreement_{flavor}.png)", f"](assets/kallisto_benchmark_10m/umi_agreement_{flavor}.png)")
        (docs / "KALLISTO_PILOT_BENCHMARK_REPORT.md").write_text(public)
    atomic_json(root / "BENCHMARK_SUCCESS.json", {"complete": True, "libraries": list(config["libraries"]),
                                                "pairs_per_library": args.pairs, "arms": list(arms),
                                                "umi_correction_options": list(UMI_FLAVORS),
                                                "report": "comparison/REPORT.md"})

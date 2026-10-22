#!/usr/bin/env python3
"""Check exact count preservation and runtime against a saved backend benchmark.

Run on a compute node. Inputs are completed benchmark roots, never raw reads.
The baseline is read only. Matrix equality uses exact numeric comparisons,
including fractional weights; assignment CSV equality includes row order.
"""
from __future__ import annotations

import argparse
import csv
import sqlite3
from collections import Counter
import json
import shutil
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy import io as sio, sparse
import yaml

try:
    from scripts.benchmark_gex_backends import ARMS, arm_options
    from scripts.benchmark_gex_report import report, markdown_table
    from scripts.gex_alignment import atomic_json, sha256_file
except ModuleNotFoundError:
    from benchmark_gex_backends import ARMS, arm_options
    from benchmark_gex_report import report, markdown_table
    from gex_alignment import atomic_json, sha256_file


def exact_product(old, new, process=False, *, serial_reference=None):
    matrix_name = "aw_uf_count_matrix.mtx" if process else "count_matrix.mtx"
    metadata_name = "aw_uf_cell_metadata.csv" if process else "cell_metadata.csv"
    old_matrix = sparse.csr_matrix(sio.mmread(old / matrix_name))
    new_matrix = sparse.csr_matrix(sio.mmread(new / matrix_name))
    old_obs = pd.read_csv(old / metadata_name, index_col="bc_wells")
    new_obs = pd.read_csv(new / metadata_name, index_col="bc_wells")
    old_features = pd.read_csv(old / "all_genes.csv", dtype=str)
    new_features = pd.read_csv(new / "all_genes.csv", dtype=str)
    pd.testing.assert_frame_equal(old_obs, new_obs, check_exact=True)
    pd.testing.assert_frame_equal(old_features, new_features, check_exact=True)
    if old_matrix.shape != new_matrix.shape:
        raise AssertionError(f"Matrix shape changed: {old} / {new}")
    saved_changed = (old_matrix != new_matrix).nnz
    reference_details = {"reference_policy": "saved_baseline", "baseline_roundtrip_wells": 0,
                         "baseline_roundtrip_entries": 0, "baseline_roundtrip_max_abs": 0.}
    if saved_changed and serial_reference is not None:
        old_matrix, reference_details = serial_reference(old_matrix, old_obs, old_features)
    changed = (old_matrix != new_matrix).nnz
    if changed:
        raise AssertionError(f"{changed} matrix entries changed: {new}")
    # The deliverable H5AD and its raw count layer must match the matrix too.
    if not process:
        import anndata as ad
        for directory, expected, obs in ((old, old_matrix, old_obs), (new, new_matrix, new_obs)):
            data = ad.read_h5ad(directory / "anndata.h5ad")
            if list(data.obs_names) != list(obs.index.astype(str)) or list(data.var["gene_id"]) != list(new_features["gene_id"]):
                raise AssertionError(f"H5AD cell/feature ordering disagrees: {directory}")
            if (sparse.csr_matrix(data.X) != expected).nnz or (sparse.csr_matrix(data.layers["counts"]) != expected).nnz:
                raise AssertionError(f"H5AD raw counts disagree: {directory}")
    return dict(cells=new_matrix.shape[0], genes=new_matrix.shape[1], nonzero_entries=new_matrix.nnz,
                changed_entries=changed, exact_counts=True, exact_metadata=True,
                exact_saved_matrix=saved_changed == 0, **reference_details)


def serial_repaired_reference(matrix, obs, features, *, baseline, arm, library):
    """Recover original serial weights only for documented repaired AWUF wells.

    The previous benchmark restored omitted three-UMI wells by parsing decimal
    CSV weights with pandas, which can change their last bits. Read those saved
    weights with Python's round-trip float parser and reproduce the original
    per-barcode Counter additions followed by SQLite SUM. Existing baseline
    rows are immutable; no rounding or tolerance is used to accept new counts.
    """
    process = baseline / arm / f"gex_{library}" / "process"
    qc = json.loads((process / "kallisto_qc.json").read_text())
    repair = qc.get("metadata_roundoff_repair", {})
    if not repair.get("complete") or not repair.get("existing_gene_weights_unchanged"):
        raise AssertionError("No documented fractional metadata repair for serial reference")
    if qc["molecules_below_read_threshold"] != 0:
        raise AssertionError("Serial CSV reconstruction requires all molecule families above the read threshold")
    original = pd.read_csv(baseline / "fractional_metadata_originals" / arm / library /
                           "aw_uf_cell_metadata.csv", index_col="bc_wells")
    restored = set(obs.index) - set(original.index)
    if len(restored) != repair["restored_awuf_wells"]:
        raise AssertionError("Baseline restored-well archive disagrees with repair manifest")
    gene_index = dict(zip(features.gene_id, range(len(features))))
    per_barcode = {}
    with (process / "tscp_assignment.csv").open(newline="") as handle:
        reader = csv.reader(handle)
        columns = {name: i for i, name in enumerate(next(reader))}
        for row in reader:
            well = row[columns["bc_wells"]]
            if well not in restored:
                continue
            barcode = row[columns["cell_barcode"]]
            if barcode not in per_barcode:
                per_barcode[barcode] = well, Counter()
            previous_well, counts = per_barcode[barcode]
            if previous_well != well:
                raise AssertionError("Baseline barcode identity maps to multiple wells")
            counts[gene_index[row[columns["gene"]]]] += float(row[columns["assignment_weight"]])
    if {well for well, _ in per_barcode.values()} != restored:
        raise AssertionError("Baseline assignments do not cover every restored well")
    connection = sqlite3.connect(":memory:")
    try:
        connection.execute("CREATE TABLE gene_counts (well TEXT, feature INTEGER, n REAL)")
        for well, counts in per_barcode.values():
            connection.executemany("INSERT INTO gene_counts VALUES (?,?,?)",
                                   [(well, feature, n) for feature, n in counts.items()])
        values = list(connection.execute("SELECT well,feature,SUM(n) FROM gene_counts "
                                         "GROUP BY well,feature ORDER BY well,feature"))
    finally:
        connection.close()
    cell_index = {cell: i for i, cell in enumerate(obs.index)}
    original_coo = matrix.tocoo()
    restored_rows = [cell_index[well] for well in restored]
    keep = ~np.isin(original_coo.row, restored_rows)
    rows = np.concatenate((original_coo.row[keep], [cell_index[well] for well, _, _ in values]))
    columns = np.concatenate((original_coo.col[keep], [feature for _, feature, _ in values]))
    data = np.concatenate((original_coo.data[keep], [value for _, _, value in values]))
    reference = sparse.csr_matrix((data, (rows, columns)), shape=matrix.shape)
    delta = (reference - matrix).tocoo()
    changed_wells = set(obs.index[np.unique(delta.row)])
    if not changed_wells <= restored:
        raise AssertionError("Serial reference changed an originally serialized baseline well")
    return reference, dict(reference_policy="baseline_assignments_original_serial_aggregation",
                           baseline_restored_wells=len(restored),
                           baseline_roundtrip_wells=len(changed_wells),
                           baseline_roundtrip_entries=delta.nnz,
                           baseline_roundtrip_max_abs=float(np.max(np.abs(delta.data))) if delta.nnz else 0.)


def compare(args):
    config = yaml.safe_load(args.config.read_text())
    baseline, current = args.baseline, args.work_dir
    output = current / "optimization"
    output.mkdir(exist_ok=True)
    (current / "OPTIMIZATION_SUCCESS.json").unlink(missing_ok=True)
    (current / "OPTIMIZATION_FAILED.json").unlink(missing_ok=True)
    equality, timings = [], []
    for arm in ARMS:
        base, flavor = arm_options(arm)
        for library in config["libraries"]:
            name, sample = library["sublibrary_name"], library["sample_name"]
            old, new = (root / arm / f"gex_{name}" for root in (baseline, current))
            old_manifest, new_manifest = (json.loads((root / "manifests" / f"{name}.{arm}.json").read_text())
                                          for root in (baseline, current))
            for manifest in (old_manifest, new_manifest):
                if not manifest.get("complete") or manifest["reads"] != 10_000_000 or manifest["threads"] != 16:
                    raise AssertionError(f"Unmatched or incomplete input/resources: {arm}/{name}")
            if old_manifest["survivor_sha256"] != new_manifest["survivor_sha256"]:
                raise AssertionError(f"Survivor input changed: {arm}/{name}")
            old_qc, new_qc = ({}, {}) if base == "star" else (
                json.loads((root / "process/kallisto_qc.json").read_text()) for root in (old, new))
            timing = dict(arm=arm, backend=base, umi_correction=flavor, library=name,
                          baseline_core_seconds=old_manifest["elapsed_seconds"],
                          optimized_core_seconds=new_manifest["elapsed_seconds"])
            timing["core_speedup"] = timing["baseline_core_seconds"] / timing["optimized_core_seconds"]
            if base != "star":
                if sha256_file(old / "process/tscp_assignment.csv") != sha256_file(new / "process/tscp_assignment.csv"):
                    raise AssertionError(f"Molecule assignments changed: {arm}/{name}")
                for key in ("pseudoaligned_reads", "molecules", "assigned_read_support", "multigene_molecules",
                            "discarded_multigene_reads", "molecules_below_read_threshold", "transcript_unresolved_molecules", "samples"):
                    if old_qc[key] != new_qc[key]:
                        raise AssertionError(f"Scientific QC changed: {arm}/{name}/{key}")
                timing.update(baseline_count_seconds=old_qc["timing"]["umi_count_seconds"],
                              optimized_count_seconds=new_qc["timing"]["umi_count_seconds"],
                              **new_qc["counting_execution"])
                timing["count_speedup"] = timing["baseline_count_seconds"] / timing["optimized_count_seconds"]
            timings.append(timing)
            # STAR is a fresh timing control. Exact preservation is required for
            # every native count product and both gene-assignment policies.
            for label, old_dir, new_dir in (("all_wells", old / "process", new / "process"),
                                           ("unfiltered", old / sample / "DGE_unfiltered", new / sample / "DGE_unfiltered"),
                                           ("filtered", old / sample / "DGE_filtered", new / sample / "DGE_filtered")):
                if base == "star":
                    continue
                reference = None
                if label == "all_wells" and old_qc.get("metadata_roundoff_repair", {}).get("complete"):
                    reference = lambda matrix, obs, features: serial_repaired_reference(
                        matrix, obs, features, baseline=baseline, arm=arm, library=name)
                product = exact_product(old_dir, new_dir, process=label == "all_wells", serial_reference=reference)
                equality.append(dict(arm=arm, backend=base, umi_correction=flavor, library=name, product=label,
                                     exact_assignments=True, **product))
            pd.DataFrame(equality).to_csv(output / "exact_equality.csv", index=False)
            pd.DataFrame(timings).to_csv(output / "runtime_comparison.csv", index=False)
    if not args.skip_backend_report:
        report(SimpleNamespace(config=args.config, work_dir=current, pairs=10_000_000, publish_report=False))
    # Validate combined products too, after the real Step 1b run above.
    if not args.skip_backend_report:
        for arm in ARMS:
            if arm_options(arm)[0] == "star":
                continue
            for sample in sorted({library["sample_name"] for library in config["libraries"]}):
                for label in ("DGE_filtered", "DGE_unfiltered"):
                    old = baseline / arm / "gex_combined" / sample / label
                    new = current / arm / "gex_combined" / sample / label
                    product = exact_product(old, new)
                    equality.append(dict(arm=arm, backend=arm_options(arm)[0], umi_correction=arm_options(arm)[1],
                                         library="combined:" + sample, product=label, **product))
    equal = pd.DataFrame(equality)
    runtime = pd.DataFrame(timings)
    equal.to_csv(output / "exact_equality.csv", index=False)
    runtime.to_csv(output / "runtime_comparison.csv", index=False)
    summary = runtime.groupby(["backend", "umi_correction"], sort=False)[
        ["baseline_core_seconds", "optimized_core_seconds", "baseline_count_seconds", "optimized_count_seconds"]].median().reset_index()
    summary["core_speedup"] = summary.baseline_core_seconds / summary.optimized_core_seconds
    summary["count_speedup"] = summary.baseline_count_seconds / summary.optimized_count_seconds
    summary.to_csv(output / "runtime_summary.csv", index=False)
    if args.skip_backend_report:
        atomic_json(output / "diagnostic_success.json", dict(exact_native_library_products=len(equal)))
        return
    make_report(args, equal, runtime, summary)
    atomic_json(current / "OPTIMIZATION_SUCCESS.json", dict(complete=True, baseline=str(baseline),
                library_runs=36, native_runs=24, exact_native_products=len(equal),
                umi_correction=["parse", "directional", "none"], exact_counts=True,
                exact_assignment_csv=True,
                reference_policy="saved_baseline_or_documented_original_serial_reconstruction",
                exact_saved_matrix_products=int(equal.exact_saved_matrix.sum()),
                baseline_csv_roundtrip_products=int((equal.baseline_roundtrip_entries > 0).sum()),
                slots_per_job=16, pairs_per_library=10_000_000))


def make_report(args, equality, runtime, summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    output = args.work_dir / "optimization"
    native = summary[summary.backend != "star"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    labels = [f"{'Unique' if row.backend == 'kallisto_unique' else 'Fractional'}\n{row.umi_correction}"
              for row in native.itertuples()]
    x = np.arange(len(native))
    for axis, before, after, title in ((axes[0], "baseline_count_seconds", "optimized_count_seconds", "Native UMI counting"),
                                      (axes[1], "baseline_core_seconds", "optimized_core_seconds", "Alignment + counting + DGE")):
        axis.bar(x - .2, native[before] / 60, .4, label="Previous serial counter")
        axis.bar(x + .2, native[after] / 60, .4, label="Optimized parallel counter")
        axis.set(xticks=x, xticklabels=labels, ylabel="Median minutes across four libraries", title=title)
        axis.legend(fontsize=8)
        axis.set_ylim(0, float(native[before].max()) / 60 * 1.3)
    fig.savefig(output / "runtime_comparison.png", dpi=180)
    plt.close(fig)
    table = summary.copy()
    for column in ("baseline_core_seconds", "optimized_core_seconds", "baseline_count_seconds", "optimized_count_seconds"):
        table[column.replace("_seconds", "_minutes")] = table[column] / 60
    text = "# Parallel kallisto counting: matched four-library benchmark\n\n"
    text += ("The optimized counter uses complete-barcode worker shards, cached EC/gene annotations, "
             "a singleton shortcut for directional correction, and an ordered disk-backed merge. "
             "Parse UMI families retain the installed Parse kernel check. STAR remains the default; "
             "the Parse 1.6.4 installation is unchanged.\n\n")
    text += ("C1, C2, T1 and T2 each use the exact saved 10-million-pair sample and identical Parse survivors "
             "from the previous benchmark. All 36 library runs were repeated with 16 bound cores and four "
             "concurrent jobs. The same mature-plus-nascent reference and kallisto executable were reused. "
             "Sampling, shared preprocessing, reference construction, and Step 1b validation are outside core time.\n\n")
    text += "Median UMI counting time (minutes), with speedup relative to the serial counter:\n\n"
    text += markdown_table(table[table.backend != "star"], ["backend", "umi_correction", "baseline_count_minutes", "optimized_count_minutes", "count_speedup"]) + "\n\n"
    text += "Median core time (minutes), including fresh STAR controls:\n\n"
    text += markdown_table(table, ["backend", "umi_correction", "baseline_core_minutes", "optimized_core_minutes", "core_speedup"]) + "\n\n"
    text += "![Runtime comparison](assets/kallisto_parallel_benchmark_10m/runtime_comparison.png)\n\n"
    text += (f"All 24 native runs retained byte-identical molecule assignment CSVs and identical scientific QC. "
             f"All {len(equality)} tested count products exactly match their serial reference, ordered cells/features and metadata, "
             "including fractional weights. H5AD X and raw count layers agree with their matrices. "
             "The checks cover unique and fractional genes with parse, directional and none correction; "
             "per-library all-well, unfiltered and filtered matrices; and filtered/unfiltered Step 1b combinations. "
             "No tolerance is used for matrix equality. Optional transcript products are covered by serial/parallel regression fixtures.\n\n")
    roundtrip = equality[equality.baseline_roundtrip_entries > 0]
    if len(roundtrip):
        text += (f"The saved all-well matrices in {len(roundtrip)} baseline products contain CSV round-trip "
                 f"differences in {int(roundtrip.baseline_roundtrip_wells.sum())} previously restored three-UMI wells "
                 f"(maximum absolute difference {roundtrip.baseline_roundtrip_max_abs.max():.3g}). "
                 "These are confined to wells omitted before the baseline's metadata repair. "
                 "Their serial reference was reconstructed from the byte-identical baseline assignments using "
                 "Python's exact float round-trip parser, the original per-barcode additions, and SQLite SUM order. "
                 "The new matrices match that reference exactly. All saved baseline files remain untouched; "
                 "the sample and combined matrices match their saved baseline matrices directly. "
                 "The equality CSV records both the raw saved-matrix comparison and the strict serial-reference check.\n\n")
    text += ("The baseline fractional Parse/directional outputs include the documented whole-UMI metadata repair. "
             "Their original core timings are retained. New runs apply whole-UMI cell calling directly. "
             "Timing differences include host and shared-filesystem variability; speedups describe this implementation and reference.\n\n")
    text += ("[Per-run runtime and counting phases](assets/kallisto_parallel_benchmark_10m/runtime_comparison.csv), "
             "[exact equality checks](assets/kallisto_parallel_benchmark_10m/exact_equality.csv), "
             "[backend UMI/cell agreement](assets/kallisto_parallel_benchmark_10m/agreement.csv), "
             "[per-library yield](assets/kallisto_parallel_benchmark_10m/library_summary.csv).\n")
    (output / "REPORT.md").write_text(text.replace("assets/kallisto_parallel_benchmark_10m/", ""))
    docs = Path(__file__).resolve().parents[1] / "docs"
    assets = docs / "assets/kallisto_parallel_benchmark_10m"
    assets.mkdir(parents=True, exist_ok=True)
    for name in ("runtime_comparison.csv", "runtime_summary.csv", "exact_equality.csv", "runtime_comparison.png"):
        shutil.copy2(output / name, assets / name)
    if not args.skip_backend_report:
        for name in ("agreement.csv", "library_summary.csv", "runtime.csv"):
            shutil.copy2(args.work_dir / "comparison" / name, output / name)
            shutil.copy2(output / name, assets / name)
    (docs / "KALLISTO_PARALLEL_BENCHMARK_REPORT.md").write_text(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--skip-backend-report", action="store_true", help="Skip Step 1b and cross-backend analysis for a diagnostic pass")
    args = parser.parse_args()
    try:
        compare(args)
    except BaseException as exc:
        atomic_json(args.work_dir / "OPTIMIZATION_FAILED.json", {"error": str(exc)})
        raise


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Step 1a native kallisto counting and Parse-compatible sparse products.

The Parse DGE entry point supplies its initialized kit/sample/reference and
cell-calling policy. Only the in-process DGE callback is replaced; vendor
sources and the STAR execution path are untouched.
"""

from __future__ import annotations

import argparse
import json
import os
import runpy
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import io as sio, sparse

try:
    from scripts.gex_alignment import (SCHEMA, atomic_json, backend_options, load_reference,
                                       resolve_path, settings_identity, sha256_file, tool_executable, check_output_backend)
    from scripts.kallisto_count_parallel import barcode_shards, count_shards
    from scripts.kallisto_umi import (umi_policy, prepare_fastqs,
                                     read_equivalence_classes, spool_bus)
except ModuleNotFoundError:
    from gex_alignment import (SCHEMA, atomic_json, backend_options, load_reference,
                               resolve_path, settings_identity, sha256_file, tool_executable, check_output_backend)
    from kallisto_count_parallel import barcode_shards, count_shards
    from kallisto_umi import (umi_policy, prepare_fastqs,
                             read_equivalence_classes, spool_bus)


def read_features(reference, transcript_names):
    features = pd.read_csv(reference / "features.tsv", sep="\t", dtype=str).fillna("")
    required = {"transcript_id", "gene_id", "gene_name", "genome", "feature_type"}
    if not required <= set(features) or features["transcript_id"].duplicated().any():
        raise ValueError("Invalid kallisto features.tsv")
    if features[list(required)].eq("").any().any():
        raise ValueError("Empty reference feature identifier/annotation")
    if not features["feature_type"].isin(["transcript", "nascent"]).all():
        raise ValueError("Unknown reference feature_type")
    names = Path(transcript_names).read_text().splitlines()
    if len(names) != len(set(names)) or set(names) != set(features["transcript_id"]):
        raise ValueError("Kallisto target order/list disagrees with reference features")
    return features.set_index("transcript_id").loc[names].reset_index()


def matrix_from_database(connection, table, cells, n_features, integer=False):
    cell_index = {name: index for index, name in enumerate(cells)}
    rows, columns, values = [], [], []
    for well, feature, count in connection.execute(
        f"SELECT well,feature,SUM(n) FROM {table} GROUP BY well,feature ORDER BY well,feature"
    ):
        rows.append(cell_index[well])
        columns.append(feature)
        values.append(count)
    dtype = np.int64 if integer else np.float64
    if integer and any(float(value) != int(value) for value in values):
        raise ValueError("Fractional counts must not be cast to integers")
    return sparse.csr_matrix((np.array(values, dtype=dtype), (rows, columns)),
                             shape=(len(cells), n_features))


def whole_molecule_totals(matrix):
    """Gene weights conserve whole UMIs; keep float roundoff out of cell calls."""
    values = np.asarray(matrix.sum(axis=1)).ravel()
    rounded = np.rint(values)
    if not np.isfinite(values).all() or not np.allclose(values, rounded, rtol=0, atol=1e-6):
        raise ValueError("Gene weights do not conserve whole-molecule totals")
    return rounded.astype(np.int64)


def count_native(bus_path, transcript_path, ec_path, reference, gene_features,
                 lookup, bc_len, umi_len, database, assignment_path, options,
                 *, min_reads=1, kernel=None, mutation_function=None, threads=1, diagnostics=None):
    phase_started = time.perf_counter()
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        raise ValueError("Counting threads must be a positive integer")
    tx_features = read_features(reference, transcript_path)
    gene_features = gene_features.reset_index(drop=True).copy()
    gene_features = gene_features.loc[:, ["gene_id", "gene_name", "genome"]].copy()
    gene_features["gene_name"] = gene_features["gene_name"].fillna(gene_features["gene_id"])
    gene_features = gene_features.astype(str)
    # Parse's in-memory reference quotes symbols for its manual CSV writer.
    # The DGE CSV reader removes those quotes before constructing feature names.
    gene_features["gene_name"] = gene_features["gene_name"].str.strip('"')
    if gene_features["gene_id"].duplicated().any():
        raise ValueError("Parse gene IDs must be unique")
    gene_index = {gene: index for index, gene in enumerate(gene_features["gene_id"])}
    unknown = set(tx_features["gene_id"]) - set(gene_index)
    if unknown:
        raise ValueError(f"Kallisto genes absent from Parse reference: {sorted(unknown)[:10]}")
    native_annotation = {row.gene_id: (row.genome, row.gene_name)
                         for row in gene_features.itertuples()}
    for row in tx_features.itertuples():
        if native_annotation[row.gene_id][0] != row.genome:
            raise ValueError(f"Reference genome label mismatch for {row.gene_id}")
    classes = read_equivalence_classes(ec_path, len(tx_features))
    reference_seconds = time.perf_counter() - phase_started
    phase_started = time.perf_counter()
    connection, aligned_reads = spool_bus(bus_path, database, lookup, classes, bc_len, umi_len)
    spool_seconds = time.perf_counter() - phase_started
    phase_started = time.perf_counter()
    stats = Counter(pseudoaligned_reads=aligned_reads, molecules=0,
                    transcript_unresolved_molecules=0, multigene_molecules=0,
                    discarded_multigene_reads=0, assigned_read_support=0,
                    molecules_below_read_threshold=0)
    tx_genes = list(tx_features["gene_id"])
    if options["gene_assignment"] == "unique":
        multigene_ecs = {ec for ec, targets in classes.items()
                         if len({tx_genes[target] for target in targets}) > 1}
        stats["discarded_multigene_reads"] = sum(
            int(count) for ec, count in connection.execute("SELECT ec,SUM(n) FROM reads GROUP BY ec")
            if ec in multigene_ecs
        )
    include_transcripts = options["output_level"] in ("both", "transcript")
    class_evidence = {ec: (frozenset(tx_genes[target] for target in targets), frozenset(targets))
                      for ec, targets in classes.items()}
    ranges = barcode_shards(connection, threads)
    workers = min(threads, len(ranges)) if ranges else 1
    connection.commit()
    connection.close()
    context = dict(lookup=lookup, classes=classes, tx_genes=tx_genes, options=options,
                   annotation=native_annotation, gene_index=gene_index, min_reads=min_reads,
                   include_transcripts=include_transcripts, kernel=kernel,
                   mutation_function=mutation_function, class_evidence=class_evidence)
    connection, counters = count_shards(database, assignment_path, ranges, workers, context)
    stats.update(counters)
    correction_seconds = time.perf_counter() - phase_started
    phase_started = time.perf_counter()
    try:
        cells = [row[0] for row in connection.execute("SELECT DISTINCT well FROM gene_counts ORDER BY well")]
        if not cells:
            raise ValueError("No gene-resolved molecules; inspect reference/strand/assignment settings")
        gene_matrix = matrix_from_database(connection, "gene_counts", cells, len(gene_features),
                                           integer=options["gene_assignment"] == "unique")
        tx_matrix = (matrix_from_database(connection, "tx_counts", cells, len(tx_features),
                                          integer=options["transcript_assignment"] == "unique")
                     if include_transcripts else None)
        support = dict(connection.execute("SELECT well,SUM(n) FROM support GROUP BY well"))
    finally:
        connection.close()
    gene_features["feature_type"] = "gene"
    gene_features["feature_id"] = gene_features["gene_id"]
    tx_features["feature_id"] = tx_features["transcript_id"]
    obs = pd.DataFrame(index=pd.Index(cells, name="bc_wells"))
    obs["tscp_count"] = whole_molecule_totals(gene_matrix)
    obs["gene_count"] = np.asarray((gene_matrix > 0).sum(axis=1)).ravel()
    obs["mread_count"] = [support.get(cell, 0) for cell in cells]
    for index in range(3):
        obs[f"bc{index + 1}_wind"] = [cell.split("_")[index] for cell in cells]
    # The precise genomic exonic fraction is unavailable from native ECs.
    obs["tscp_exonic"] = np.nan
    if diagnostics is not None:
        diagnostics.update(reference_seconds=reference_seconds, spool_seconds=spool_seconds,
                           correction_merge_seconds=correction_seconds,
                           matrix_seconds=time.perf_counter() - phase_started,
                           workers=workers, shards=len(ranges), ordering="barcode_lexical_v1")
    return gene_matrix, gene_features, tx_matrix, tx_features, obs, dict(stats)


def write_product(directory, matrix, features, obs, provenance, level="gene"):
    import anndata as ad

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    matrix = matrix.tocsr()
    if matrix.shape != (len(obs), len(features)):
        raise ValueError("Matrix dimensions disagree with ordered feature/cell metadata")
    obs = obs.copy()
    obs["matrix_count"] = np.asarray(matrix.sum(axis=1)).ravel()
    obs["matrix_feature_count"] = np.asarray((matrix > 0).sum(axis=1)).ravel()
    (directory / ("all_transcripts.csv" if level == "gene" else "all_genes.csv")).unlink(missing_ok=True)
    sio.mmwrite(directory / "count_matrix.mtx", matrix)
    obs.to_csv(directory / "cell_metadata.csv")
    features.to_csv(directory / "features.csv", index=False)
    features.to_csv(directory / ("all_genes.csv" if level == "gene" else "all_transcripts.csv"),
                    index=(level == "gene"))
    var = features.copy()
    if level == "gene":
        var.index = pd.Index(var["gene_name"].fillna(var["gene_id"]) + "_" + var["genome"], name="gene")
    else:
        var.index = pd.Index(var["transcript_id"], name="transcript")
    # Parse enables pandas' future string inference. Write legacy H5AD string
    # encodings without changing global pandas/AnnData settings in production.
    for frame in (obs, var):
        frame.index = pd.Index(np.asarray(frame.index, dtype=object), dtype=object, name=frame.index.name)
        for column in frame:
            values = frame[column]
            if isinstance(values.dtype, pd.StringDtype):
                if values.isna().any():
                    categories = pd.Index(values.dropna().unique().astype(object), dtype=object)
                    frame[column] = pd.Categorical(values, categories=categories)
                else:
                    frame[column] = values.astype(object)
    adata = ad.AnnData(X=matrix, obs=obs.copy(), var=var)
    adata.var_names_make_unique()
    # Older AnnData index setters re-infer strings even from object indices.
    for frame in (adata.obs, adata.var):
        frame.index = pd.Index(np.asarray(frame.index, dtype=object), dtype=object, name=frame.index.name)
    adata.layers["counts"] = matrix.copy()
    adata.uns["gex_quantification"] = provenance
    adata.write_h5ad(directory / "anndata.h5ad", convert_strings_to_categoricals=False)


def publish_samples(spipe, root, matrices, options, provenance):
    """Use Parse's initialized sample masks and original cell-cutoff method."""
    gene_matrix, genes, tx_matrix, transcripts, obs, stats = matrices
    minimum = spipe.get_par_val("mol_awuf_tscp_thresh", as_int=True)
    obs["in_mat"] = (obs["tscp_count"] >= minimum).astype(int)
    process = root / "process"
    process.mkdir(parents=True, exist_ok=True)
    obs.to_csv(process / "aw_uf_bc_counts.csv")
    genes.to_csv(process / "all_genes.csv")
    awuf = obs["in_mat"].to_numpy().astype(bool)
    sio.mmwrite(process / "aw_uf_count_matrix.mtx", gene_matrix[awuf])
    obs.loc[awuf].to_csv(process / "aw_uf_cell_metadata.csv")
    sample_qc = {}
    # Multi-genome classification uses the same strict dominant-species threshold.
    purity = spipe.get_par_val("dge_species_purity_thresh", as_float=True)
    genomes = sorted(genes["genome"].unique())
    species_counts = np.column_stack([
        np.asarray(gene_matrix[:, genes["genome"].eq(genome).to_numpy()].sum(axis=1)).ravel()
        for genome in genomes
    ])
    fractions = species_counts / np.maximum(obs["tscp_count"].to_numpy()[:, None], 1)
    obs["species"] = [genomes[int(np.argmax(row))] if len(genomes) == 1 or max(row) >= purity
                      else "multiplet" for row in fractions]
    for sample in spipe.get_samples():
        if options["output_level"] != "both":
            for suffix in ("filtered", "unfiltered"):
                stale = root / sample.get_name() / f"DGE_transcript_{suffix}"
                if stale.is_dir():
                    shutil.rmtree(stale)
        sample_obs = sample.bc_df_for_sample(obs).copy()
        sample_obs = sample_obs[sample_obs["in_mat"] > 0]
        if sample_obs.empty:
            raise ValueError(f"No molecules in the allowed wells for {sample.get_name()}")
        threshold = options["cell_umi_threshold"]
        if threshold is None:
            ok, recipe = spipe.get_tscp_cutoff(sample_obs[["tscp_count"]], sample)
            if not ok:
                raise ValueError(f"Parse cell calling failed for {sample.get_name()}: {recipe}")
            threshold = recipe["cell_tscp_cutoff"]
        else:
            recipe = {"cell_call": "given_umi_threshold", "cell_tscp_cutoff": threshold}
        sample_obs["sample"] = sample.get_name()
        for index in range(3):
            barcode_data = spipe.get_bc_info().get_bc_data(index + 1)
            seq_to_wind = spipe.get_bc_info().get_seq_to_windzp(index + 1)
            well_map = {str(seq_to_wind[row.sequence]): str(row.well)
                        for row in barcode_data.itertuples()}
            sample_obs[f"bc{index + 1}_well"] = sample_obs[f"bc{index + 1}_wind"].map(well_map)
        for filtered, cutoff in ((False, spipe.get_par_val("dge_tscp_cut_umin", as_int=True)),
                                 (True, threshold)):
            selected = sample_obs[sample_obs["tscp_count"] >= cutoff].copy()
            rows = obs.index.get_indexer(selected.index)
            suffix = "filtered" if filtered else "unfiltered"
            primary_level = "transcript" if options["output_level"] == "transcript" else "gene"
            primary = tx_matrix if primary_level == "transcript" else gene_matrix
            features = transcripts if primary_level == "transcript" else genes
            write_product(root / sample.get_name() / f"DGE_{suffix}", primary[rows], features,
                          selected, dict(provenance, feature_level=primary_level), primary_level)
            if options["output_level"] == "both":
                write_product(root / sample.get_name() / f"DGE_transcript_{suffix}",
                              tx_matrix[rows], transcripts, selected,
                              dict(provenance, feature_level="transcript"), "transcript")
        sample_qc[sample.get_name()] = dict(recipe, called_cells=int((sample_obs["tscp_count"] >= threshold).sum()))
    atomic_json(process / "kallisto_qc.json", dict(stats, samples=sample_qc,
                                                  exonic_fraction="unavailable_from_native_equivalence_classes"))


def run_counting(spipe, config, config_path, library, threads):
    started = time.perf_counter()
    timings = {}
    options = backend_options(config)
    options["executable"] = tool_executable(options["executable"], config_path)
    reference = resolve_path(options["reference_dir"], config_path)
    reference_manifest = load_reference(reference)
    root = resolve_path(config["paths"]["output_dir"], config_path) / f"gex_{library}"
    process = root / "process"
    completion = process / "kallisto.complete.json"
    completion.unlink(missing_ok=True)
    source = process / "barcode_head.fastq"
    if not source.is_file():
        source = source.with_suffix(".fastq.gz")
    if not source.is_file():
        raise FileNotFoundError("Run Parse preprocessing before kallisto counting")
    from splitpipe import cy_utils
    from splitpipe.molinfo import collapse_polyN

    before = {name: sha256_file(process / name)
              for name in ("barcode_data.csv", "run_proc_def.json")}
    version = subprocess.check_output([options["executable"], "version"], text=True).strip()
    kallisto_dir = process / "kallisto"
    kallisto_dir.mkdir(exist_ok=True)
    # Never consume stale BUS files after a failed subprocess.
    with tempfile.TemporaryDirectory(prefix="csm-kallisto-") as tmp:
        tmp = Path(tmp)
        before_prepare = time.perf_counter()
        lookup, bc_len, umi_len, preparation = prepare_fastqs(
            source, tmp / "cdna.fastq", tmp / "tags.fastq"
        )
        timings["input_format_seconds"] = time.perf_counter() - before_prepare
        command = [options["executable"], "bus", "-i", str(reference / "index.idx"),
                   "-o", str(tmp / "bus"), "-x", f"1,0,{bc_len}:1,{bc_len},{bc_len + umi_len}:0,0,0",
                   "-t", str(threads)]
        strand_flag = {"forward": "--fr-stranded", "reverse": "--rf-stranded", "unstranded": "--unstranded"}
        command += [strand_flag[options["strand"]], str(tmp / "cdna.fastq"), str(tmp / "tags.fastq")]
        with (kallisto_dir / "alignment.log").open("w") as log:
            before_alignment = time.perf_counter()
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
        timings["alignment_seconds"] = time.perf_counter() - before_alignment
        native = tmp / "bus"
        before_count = time.perf_counter()
        count_execution = {}
        matrices = count_native(
            native / "output.bus", native / "transcripts.txt", native / "matrix.ec", reference,
            spipe.get_gene_info()["all_genes"], lookup, bc_len, umi_len,
            tmp / "counts.sqlite", process / "tscp_assignment.csv", options,
            min_reads=spipe.get_par_val("mol_awuf_read_thresh", as_int=True),
            kernel=collapse_polyN, mutation_function=cy_utils.single_mut_seqs,
            threads=min(threads, options.get("count_threads", 0) or threads), diagnostics=count_execution,
        )
        timings["umi_count_seconds"] = time.perf_counter() - before_count
        if matrices[-1]["pseudoaligned_reads"] > preparation["submitted_reads"]:
            raise ValueError("BUS read count exceeds submitted input count")
        # Small native artifacts retain EC membership and mapping diagnostics.
        for name in ("transcripts.txt", "matrix.ec", "run_info.json"):
            if (native / name).exists():
                (kallisto_dir / name).write_bytes((native / name).read_bytes())
    provenance = {
        "backend": "kallisto", "kallisto_version": version,
        "umi_correction": options["umi_correction"], "umi_policy": umi_policy(options["umi_correction"]),
        "gene_assignment": options["gene_assignment"], "transcript_assignment": options["transcript_assignment"],
        "output_level": options["output_level"], "reference_sha256": sha256_file(reference / "reference.json"),
        "reference_workflow": reference_manifest["workflow"], "matrix_units": "UMI_molecules",
        "cell_umi_total_policy": "whole_molecule_totals_v1",
        "fractional_gene_counts": options["gene_assignment"] == "fractional",
    }
    before_publish = time.perf_counter()
    publish_samples(spipe, root, matrices, options, provenance)
    timings["publish_seconds"] = time.perf_counter() - before_publish
    timings["count_callback_seconds"] = time.perf_counter() - started
    qc_path = process / "kallisto_qc.json"
    qc = json.loads(qc_path.read_text())
    qc["timing"] = timings
    qc["counting_execution"] = count_execution
    atomic_json(qc_path, qc)
    # Keep the common molecule counters usable without inventing STAR mapping
    # uniqueness or alignment-coordinate statistics for pseudoalignments.
    stat_path = process / "pre_align_stats.csv"
    pipeline_stats = pd.read_csv(stat_path, comment="#", index_col="statistic")
    updates = {
        "reads_map_transcriptome": matrices[-1]["assigned_read_support"],
        "number_of_tscp": matrices[-1]["molecules"],
        "reads_pseudoalign_input": preparation["submitted_reads"],
        "reads_pseudoaligned": matrices[-1]["pseudoaligned_reads"],
        "reads_multigene_excluded": matrices[-1]["discarded_multigene_reads"],
    }
    for key, value in updates.items():
        pipeline_stats.loc[key, "value"] = value
    pipeline_stats.to_csv(stat_path)
    atomic_json(kallisto_dir / "input.json", dict(preparation, survivor_fastq_sha256=sha256_file(source),
                                                 command=command, **provenance))
    # run_proc_def is rewritten by Parse after the callback, so finish seals the
    # contract only after the vendor entry point has returned successfully.
    atomic_json(process / "kallisto.pending.json", {
        "schema_version": SCHEMA, "backend": "kallisto", "library": library,
        "settings_sha256": settings_identity(config), "reference_sha256": provenance["reference_sha256"],
        "barcode_data_sha256": before["barcode_data.csv"], "provenance": provenance,
    })
    return True


def seal(config, config_path, library):
    root = resolve_path(config["paths"]["output_dir"], config_path)
    library_root = root / f"gex_{library}"
    process = library_root / "process"
    pending = process / "kallisto.pending.json"
    payload = json.loads(pending.read_text())
    if sha256_file(process / "barcode_data.csv") != payload.pop("barcode_data_sha256"):
        raise ValueError("Parse barcode contract changed during counting")
    payload["preprocessing"] = {name: sha256_file(process / name)
                                 for name in ("barcode_data.csv", "run_proc_def.json")}
    filter_path = root / "gex_grna_filter" / library / "rejected_pairs.complete.json"
    if config.get("gex_decontamination", {}).get("enabled", False):
        payload["filter_manifest_sha256"] = sha256_file(filter_path)
    artifacts = [process / name for name in (
        "tscp_assignment.csv", "aw_uf_count_matrix.mtx", "aw_uf_cell_metadata.csv",
        "aw_uf_bc_counts.csv", "all_genes.csv", "kallisto_qc.json", "kallisto/input.json")]
    for entry in config["libraries"]:
        if entry["sublibrary_name"] != library:
            continue
        sample_root = library_root / entry["sample_name"]
        for directory in sorted(sample_root.glob("DGE_*")):
            artifacts.extend(directory / name for name in
                             ("count_matrix.mtx", "cell_metadata.csv", "features.csv", "anndata.h5ad"))
            artifacts.extend(directory / name for name in ("all_genes.csv", "all_transcripts.csv")
                             if (directory / name).is_file())
    payload["artifacts"] = {path.relative_to(library_root).as_posix(): sha256_file(path) for path in artifacts}
    payload["complete"] = True
    atomic_json(process / "kallisto.complete.json", payload)
    pending.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--library", required=True)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    config_path = Path(args.config).resolve()
    config = yaml.safe_load(config_path.read_text())
    if backend_options(config)["backend"] != "kallisto":
        parser.error("gex_alignment.backend must be kallisto")
    if sum(lib["sublibrary_name"] == args.library for lib in config["libraries"]) != 1:
        parser.error("--library must identify exactly one configured library")
    splitpipe_dir = resolve_path(config["paths"]["splitpipe_dir"], config_path)
    output = resolve_path(config["paths"]["output_dir"], config_path) / f"gex_{args.library}"
    check_output_backend(backend_options(config), output / "process")
    sys.dont_write_bytecode = True
    try:
        from scripts.splitpipe_gex_wrapper import validate_splitpipe_tree
    except ModuleNotFoundError:
        from splitpipe_gex_wrapper import validate_splitpipe_tree
    validate_splitpipe_tree(splitpipe_dir)
    sys.path.insert(0, str(splitpipe_dir))
    from splitpipe import dge
    # A runtime callback, never a patch written into the Parse distribution.
    dge.run_dge = lambda spipe: run_counting(spipe, config, config_path, args.library, args.threads)
    binary = splitpipe_dir / "split-pipe"
    sys.argv = [str(binary), "--mode", "dge", "--one_step", "--output_dir", str(output),
                "--nthreads", str(args.threads), "--kit_score_skip"]
    try:
        runpy.run_path(str(binary), run_name="__main__")
    except SystemExit as exc:
        if exc.code not in (0, None):
            raise
    seal(config, config_path, args.library)


if __name__ == "__main__":
    main()

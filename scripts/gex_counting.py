"""Shared compatibility-set molecule counting, Parse cell calls and sparse products.

STAR and kallisto adapters supply evidence; assignment and UMI policies are
independent of the aligner. No vendor installation is modified.
"""
from __future__ import annotations

import shutil
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as sio, sparse

try:
    from scripts.gex_alignment import atomic_json
    from scripts.kallisto_count_parallel import barcode_shards, count_shards
except ModuleNotFoundError:
    from gex_alignment import atomic_json
    from kallisto_count_parallel import barcode_shards, count_shards


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


def count_evidence(connection, classes, tx_features, gene_features, lookup, database,
                   assignment_path, options, *, aligned_reads, min_reads=1, kernel=None,
                   mutation_function=None, threads=1, diagnostics=None,
                   reference_seconds=0.0, spool_seconds=0.0):
    """Count backend-neutral compatibility evidence from a disk-backed read spool."""
    try:
        return _count_evidence(
            connection, classes, tx_features, gene_features, lookup, database,
            assignment_path, options, aligned_reads=aligned_reads, min_reads=min_reads,
            kernel=kernel, mutation_function=mutation_function, threads=threads,
            diagnostics=diagnostics, reference_seconds=reference_seconds,
            spool_seconds=spool_seconds)
    finally:
        # Validation can fail before the spool is closed for worker handoff.
        # Release the handle on that path too, including on shared filesystems.
        connection.close()


def _count_evidence(connection, classes, tx_features, gene_features, lookup, database,
                    assignment_path, options, *, aligned_reads, min_reads, kernel,
                    mutation_function, threads, diagnostics, reference_seconds, spool_seconds):
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        raise ValueError("Counting threads must be a positive integer")
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
        raise ValueError(f"Assigned genes absent from Parse reference: {sorted(unknown)[:10]}")
    native_annotation = {row.gene_id: (row.genome, row.gene_name)
                         for row in gene_features.itertuples()}
    for row in tx_features.itertuples():
        if native_annotation[row.gene_id][0] != row.genome:
            raise ValueError(f"Reference genome label mismatch for {row.gene_id}")
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
    atomic_json(process / ("kallisto_qc.json" if options.get("backend", "kallisto") == "kallisto" else "star_qc.json"), dict(stats, samples=sample_qc,
                                                  exonic_fraction="unavailable_from_native_equivalence_classes"))

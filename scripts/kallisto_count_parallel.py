"""Deterministic, disk-backed barcode shards for native molecule counting.

A complete corrected barcode belongs to one shard. Fork shares immutable EC
annotations and the installed Parse kernel; workers never modify vendor files.
Shard summaries are streamed back in barcode order, preserving floating-point
addition and assignment CSV order from the original serial implementation.
"""
from __future__ import annotations

import csv
import itertools
import multiprocessing
import pickle
import shutil
import sqlite3
import tempfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

try:
    from scripts.kallisto_umi import count_barcode
except ModuleNotFoundError:
    from kallisto_umi import count_barcode


_CONTEXT = None
ASSIGNMENT_COLUMNS = ["bc_wells", "genome", "gene", "gene_name", "count", "exonic",
                      "rt_type", "bc_bcis", "cell_barcode", "polyN", "assignment_weight"]


def barcode_shards(connection, workers):
    """Contiguous ranges balanced by BUS record count, never splitting a cell."""
    sizes = list(connection.execute("SELECT barcode,COUNT(*) FROM reads GROUP BY barcode ORDER BY barcode"))
    if not sizes:
        return []
    target = max(1, (sum(n for _, n in sizes) + workers * 4 - 1) // (workers * 4))
    result, first, total = [], None, 0
    for barcode, n in sizes:
        if first is None:
            first = barcode
        total += n
        if total >= target:
            result.append((first, barcode))
            first, total = None, 0
    if first is not None:
        result.append((first, sizes[-1][0]))
    return result


def _count_shard(task):
    database, directory, index, first, last = task
    context = _CONTEXT
    connection = sqlite3.connect(Path(database).resolve().as_uri() + "?mode=ro", uri=True)
    csv_path = Path(directory) / f"{index}.csv"
    summary_path = Path(directory) / f"{index}.pickle"
    stats = Counter()
    try:
        cursor = connection.execute(
            "SELECT barcode,umi,ec,SUM(n) FROM reads WHERE barcode>=? AND barcode<=? "
            "GROUP BY barcode,umi,ec ORDER BY barcode,umi,ec", (first, last))
        with csv_path.open("w", newline="") as handle, summary_path.open("wb") as summaries:
            writer = csv.writer(handle)
            for barcode, rows in itertools.groupby(cursor, key=lambda row: row[0]):
                wells, bcis, primer, barcode_parts = context["lookup"][barcode]
                gene_counts, tx_counts = Counter(), Counter()
                read_support = 0
                kwargs = {"kernel": context["kernel"]}
                if context["mutation_function"] is not None:
                    kwargs["mutation_function"] = context["mutation_function"]
                options = context["options"]
                molecules = count_barcode(
                    ((umi, ec, reads) for _, umi, ec, reads in rows), context["classes"],
                    context["tx_genes"], fractional=options["gene_assignment"] == "fractional",
                    umi_correction=options.get("umi_correction", "parse"),
                    class_evidence=context["class_evidence"], **kwargs)
                for umi, genes, reads, targets in molecules:
                    stats["molecules"] += 1
                    stats["multigene_molecules"] += len(genes) > 1
                    stats["assigned_read_support"] += reads
                    weight = 1.0 / len(genes)
                    for gene in genes:
                        genome, name = context["annotation"][gene]
                        writer.writerow([wells, genome, gene, name, reads, "", primer, bcis,
                                         barcode_parts, umi, weight])
                    if reads < context["min_reads"]:
                        stats["molecules_below_read_threshold"] += 1
                        continue
                    read_support += reads
                    for gene in genes:
                        gene_counts[context["gene_index"][gene]] += weight
                    if context["include_transcripts"]:
                        if not targets or (options["transcript_assignment"] == "unique" and len(targets) != 1):
                            stats["transcript_unresolved_molecules"] += 1
                        else:
                            for target in targets:
                                tx_counts[target] += 1.0 / len(targets)
                # Preserve per-barcode insertion order, including distinct RT identities
                # mapping to the same well; never sum fractional weights in workers.
                pickle.dump((wells, list(gene_counts.items()), list(tx_counts.items()), read_support),
                            summaries, protocol=pickle.HIGHEST_PROTOCOL)
    finally:
        connection.close()
    return csv_path, summary_path, stats


def count_shards(database, assignment_path, ranges, workers, context):
    """Return an open result database and scientific counters; propagate failures."""
    global _CONTEXT
    if workers > 1 and "fork" not in multiprocessing.get_all_start_methods():
        raise RuntimeError("Parallel native counting requires fork; set count_threads: 1 on this platform")
    _CONTEXT = context
    connection = None
    try:
        with tempfile.TemporaryDirectory(prefix="umi-shards-", dir=Path(database).parent) as directory:
            tasks = [(str(database), directory, index, first, last)
                     for index, (first, last) in enumerate(ranges)]
            executor = (ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("fork"))
                        if workers > 1 and tasks else None)
            try:
                results = executor.map(_count_shard, tasks, chunksize=1) if executor else map(_count_shard, tasks)
                # Workers may still start new read transactions while earlier
                # shards merge. Never write the input spool: SQLite rollback
                # journals can lock late readers even across unrelated tables.
                connection = sqlite3.connect(str(database) + ".counts")
                connection.execute("PRAGMA journal_mode=OFF")
                connection.execute("PRAGMA synchronous=OFF")
                connection.execute("CREATE TABLE gene_counts (well TEXT, feature INTEGER, n REAL)")
                connection.execute("CREATE TABLE tx_counts (well TEXT, feature INTEGER, n REAL)")
                connection.execute("CREATE TABLE support (well TEXT, n INTEGER)")
                stats = Counter()
                with open(assignment_path, "w", newline="") as assignments:
                    csv.writer(assignments).writerow(ASSIGNMENT_COLUMNS)
                    # Executor.map yields submission order, independent of completion.
                    for csv_path, summary_path, counters in results:
                        stats.update(counters)
                        with csv_path.open(newline="") as source:
                            shutil.copyfileobj(source, assignments, length=1024 * 1024)
                        with summary_path.open("rb") as summaries:
                            while True:
                                try:
                                    well, genes, transcripts, support = pickle.load(summaries)
                                except EOFError:
                                    break
                                connection.executemany("INSERT INTO gene_counts VALUES (?,?,?)",
                                                       [(well, feature, n) for feature, n in genes])
                                connection.executemany("INSERT INTO tx_counts VALUES (?,?,?)",
                                                       [(well, feature, n) for feature, n in transcripts])
                                connection.execute("INSERT INTO support VALUES (?,?)", (well, support))
                        csv_path.unlink()
                        summary_path.unlink()
                connection.commit()
            finally:
                if executor:
                    executor.shutdown(wait=True, cancel_futures=True)
        return connection, stats
    except BaseException:
        if connection is not None:
            connection.close()
        raise
    finally:
        _CONTEXT = None

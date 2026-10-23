#!/usr/bin/env python3
"""STAR BAM compatibility evidence for the shared fractional molecule counter.

Input is Parse's name-sorted BAM before molinfo (which drops secondary hits).
Gene candidates use Parse's stranded, whole-read gene-span containment rule.
All overlapping candidates survive; Parse's best-gene ranking is not applied.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sqlite3
import tempfile
import time
from collections import Counter
from pathlib import Path

import pysam
import yaml
import pandas as pd

try:
    from scripts.gex_alignment import (SCHEMA, atomic_json, backend_options, library_output_dir,
                                      settings_identity, sha256_file, star_reference_identity)
    from scripts.gex_counting import count_evidence, publish_samples
    from scripts.gex_runtime import run_parse_counting
    from scripts.kallisto_umi import umi_policy
except ModuleNotFoundError:
    from gex_alignment import (SCHEMA, atomic_json, backend_options, library_output_dir,
                               settings_identity, sha256_file, star_reference_identity)
    from gex_counting import count_evidence, publish_samples
    from gex_runtime import run_parse_counting
    from kallisto_umi import umi_policy


def candidate_genes(read, gene_info):
    """Retain every gene satisfying Parse's strand and containment geometry."""
    positions = read.get_reference_positions()
    if not positions:
        return set()
    step = gene_info["gtf_dict_stepsize"]
    strand = "-" if read.is_reverse else "+"
    key = f"{read.reference_name}:{positions[0] // step * step}:{strand}"
    return {gene for gene in gene_info["gene_bins"].get(key, ())
            if gene_info["gene_starts"][gene] <= positions[0]
            and positions[-1] <= gene_info["gene_ends"][gene]}


def spool_bam(path, database, gene_info, *, max_multimaps=3):
    """One read contributes one compatibility record, regardless of hit count.

    NH and primary-record checks reject BAMs missing alternative alignments.
    Unique original read identifiers are required to distinguish PCR reads from
    multiple alignments of the same read. No BAM index or full-file RAM load.
    """
    connection = sqlite3.connect(database)
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute("CREATE TABLE reads (barcode TEXT, umi TEXT, ec INTEGER, n INTEGER)")
    lookup, classes, class_ids = {}, {}, {}
    features = gene_info["all_genes"].reset_index(drop=True).copy()
    gene_index = {gene: index for index, gene in enumerate(features["gene_id"])}
    stats = Counter(input_reads=0, mapped_reads=0, gene_compatible_reads=0,
                    invalid_umi_reads=0, nongenic_reads=0, over_multimap_limit_reads=0)
    batch = []
    try:
        with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
            if bam.header.to_dict().get("HD", {}).get("SO") != "queryname":
                raise ValueError("Fractional STAR requires a queryname-sorted BAM")
            for qname, records in itertools.groupby(bam, key=lambda read: read.query_name):
                records = list(records)
                if "__OH_" not in qname:
                    raise ValueError("Fractional STAR requires Parse pre_head_keep_original=True (__OH_ IDs)")
                if any(read.is_paired or read.is_supplementary for read in records):
                    raise ValueError("Expected single-end Parse cDNA alignments without supplementary records")
                primary = [read for read in records if not read.is_secondary]
                if len(primary) != 1:
                    raise ValueError("BAM lacks one primary record per original read; check duplicate QNAMEs")
                stats["input_reads"] += 1
                mapped = [read for read in records if not read.is_unmapped]
                if not mapped:
                    continue
                nh = {read.get_tag("NH") if read.has_tag("NH") else None for read in mapped}
                if nh != {len(mapped)}:
                    raise ValueError("Incomplete STAR alignment set (NH mismatch); use the pre-molinfo sorted BAM")
                stats["mapped_reads"] += 1
                if len(mapped) > max_multimaps:
                    stats["over_multimap_limit_reads"] += 1
                    continue
                # STAR MAPQ 3/1/0 denotes multi-mapping, not bad alignment.
                # NH completeness replaces the Parse unique-only MAPQ 255 gate.
                genes = set().union(*(candidate_genes(read, gene_info) for read in mapped))
                if not genes:
                    stats["nongenic_reads"] += 1
                    continue
                if genes.difference(gene_index):
                    raise ValueError("STAR annotation gene absent from Parse feature table")
                parts = qname.split("__")
                if len(parts) < 6:
                    raise ValueError("Malformed Parse BAM QNAME")
                wells, primer, bcis, barcode_parts, umi = parts[:5]
                if any(len(value.split("_")) != 3 for value in (wells, bcis, barcode_parts)):
                    raise ValueError("Malformed Parse barcode identity")
                # Match BUS representability for harmonized UMI evidence.
                if not umi or any(base not in "ACGT" for base in umi):
                    stats["invalid_umi_reads"] += 1
                    continue
                barcode = barcode_parts.replace("_", "")
                identity = (wells, bcis, primer, barcode_parts)
                if barcode in lookup and lookup[barcode] != identity:
                    raise ValueError("One corrected barcode maps to conflicting Parse identities")
                lookup[barcode] = identity
                targets = tuple(sorted(gene_index[gene] for gene in genes))
                if targets not in class_ids:
                    class_ids[targets] = len(classes)
                    classes[class_ids[targets]] = targets
                batch.append((barcode, umi, class_ids[targets], 1))
                stats["gene_compatible_reads"] += 1
                if len(batch) >= 100_000:
                    connection.executemany("INSERT INTO reads VALUES (?,?,?,?)", batch)
                    batch.clear()
        connection.executemany("INSERT INTO reads VALUES (?,?,?,?)", batch)
        connection.execute("CREATE INDEX reads_order ON reads(barcode,umi,ec)")
        connection.commit()
    except BaseException:
        connection.close()
        raise
    # Synthetic gene-level compatibility targets stay internal. They never
    # become transcript output; no isoform identity can be inferred here.
    features["transcript_id"] = features["gene_id"]
    features["feature_type"] = "gene"
    return connection, classes, features, lookup, dict(stats)


def run_counting(spipe, config, config_path, library, threads):
    options = backend_options(config)
    root = library_output_dir(config, config_path, library)
    process = root / "process"
    bam = process / "barcode_headAligned_sorted.bam"
    from splitpipe import cy_utils
    from splitpipe.molinfo import collapse_polyN

    before_barcode = sha256_file(process / "barcode_data.csv")
    before_bam = sha256_file(bam)
    reference_hash = star_reference_identity(config, config_path)
    started = time.perf_counter()
    execution = {}
    with tempfile.TemporaryDirectory(prefix="csm-star-count-") as temporary:
        database = Path(temporary) / "reads.sqlite"
        connection, classes, features, lookup, stats = spool_bam(
            bam, database, spipe.get_gene_info(), max_multimaps=options["max_multimaps"])
        spool_seconds = time.perf_counter() - started
        matrices = count_evidence(
            connection, classes, features, spipe.get_gene_info()["all_genes"], lookup,
            database, process / "tscp_assignment.csv", options,
            aligned_reads=stats["gene_compatible_reads"],
            min_reads=spipe.get_par_val("mol_awuf_read_thresh", as_int=True),
            kernel=collapse_polyN, mutation_function=cy_utils.single_mut_seqs,
            threads=min(threads, options["count_threads"] or threads), diagnostics=execution,
            spool_seconds=spool_seconds)
    matrices[-1].pop("pseudoaligned_reads")
    matrices[-1].update(stats)
    provenance = {"backend": "star", "gene_assignment": "fractional", "output_level": "gene",
                  "umi_correction": options["umi_correction"], "umi_policy": umi_policy(options["umi_correction"]),
                  "fractional_gene_counts": True, "matrix_units": "UMI_molecules",
                  "cell_umi_total_policy": "whole_molecule_totals_v1",
                  "annotation_policy": "parse_stranded_gene_span_all_candidates_v1",
                  "reference_sha256": reference_hash, "max_multimaps": options["max_multimaps"]}
    publish_samples(spipe, root, matrices, options, provenance)
    # Hash before and after counting: never seal an input changed mid-run.
    if before_bam != sha256_file(bam):
        raise ValueError("STAR BAM changed during counting")
    qc_path = process / "star_qc.json"
    qc = json.loads(qc_path.read_text())
    qc.update(counting_execution=execution, count_callback_seconds=time.perf_counter() - started)
    qc["exonic_fraction"] = "not_reported_for_gene_compatibility_sets"
    atomic_json(qc_path, qc)
    stat_path = process / "pre_align_stats.csv"
    pipeline_stats = pd.read_csv(stat_path, comment="#", index_col="statistic")
    for key, value in {"reads_map_transcriptome": matrices[-1]["assigned_read_support"],
                       "number_of_tscp": matrices[-1]["molecules"],
                       "reads_multigene_excluded": 0}.items():
        pipeline_stats.loc[key, "value"] = value
    pipeline_stats.to_csv(stat_path)
    atomic_json(process / "star.pending.json", {
        "schema_version": SCHEMA, "backend": "star", "library": library,
        "settings_sha256": settings_identity(config), "reference_sha256": reference_hash,
        "barcode_data_sha256": before_barcode, "provenance": provenance,
        "inputs": {bam.name: before_bam}})
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--library", required=True)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    path = Path(args.config).resolve()
    config = yaml.safe_load(path.read_text())
    options = backend_options(config)
    if options["backend"] != "star" or options["gene_assignment"] != "fractional":
        parser.error("requires STAR + fractional gene assignment")
    run_parse_counting(config, path, args.library, args.threads, run_counting)


if __name__ == "__main__":
    main()

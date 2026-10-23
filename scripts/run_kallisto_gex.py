#!/usr/bin/env python3
"""Step 1a native kallisto counting and Parse-compatible sparse products.

The Parse DGE entry point supplies its initialized kit/sample/reference and
cell-calling policy. Only the in-process DGE callback is replaced; vendor
sources are untouched; STAR fractional counting uses the same counting core.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import tempfile
import time
from pathlib import Path

import pandas as pd
import yaml

try:
    from scripts.gex_alignment import (SCHEMA, atomic_json, backend_options, load_reference,
                                       resolve_path, settings_identity, sha256_file, tool_executable)
    from scripts.kallisto_umi import (umi_policy, prepare_fastqs,
                                     read_equivalence_classes, spool_bus)
except ModuleNotFoundError:
    from gex_alignment import (SCHEMA, atomic_json, backend_options, load_reference,
                               resolve_path, settings_identity, sha256_file, tool_executable)
    from kallisto_umi import (umi_policy, prepare_fastqs,
                             read_equivalence_classes, spool_bus)


# Re-export the previous helper API for benchmark scripts and existing callers.
try:
    from scripts.gex_counting import (count_evidence, read_features, matrix_from_database,
                                      whole_molecule_totals, publish_samples, write_product)
except ModuleNotFoundError:
    from gex_counting import (count_evidence, read_features, matrix_from_database,
                              whole_molecule_totals, publish_samples, write_product)


def count_native(bus_path, transcript_path, ec_path, reference, gene_features,
                 lookup, bc_len, umi_len, database, assignment_path, options,
                 *, min_reads=1, kernel=None, mutation_function=None, threads=1, diagnostics=None):
    """Adapt native BUS evidence to the shared molecule counter."""
    started = time.perf_counter()
    tx_features = read_features(reference, transcript_path)
    classes = read_equivalence_classes(ec_path, len(tx_features))
    reference_seconds = time.perf_counter() - started
    started = time.perf_counter()
    connection, aligned_reads = spool_bus(bus_path, database, lookup, classes, bc_len, umi_len)
    return count_evidence(connection, classes, tx_features, gene_features, lookup,
                          database, assignment_path, options, aligned_reads=aligned_reads,
                          min_reads=min_reads, kernel=kernel, mutation_function=mutation_function,
                          threads=threads, diagnostics=diagnostics, reference_seconds=reference_seconds,
                          spool_seconds=time.perf_counter() - started)


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


try:
    from scripts.gex_runtime import seal, run_parse_counting
except ModuleNotFoundError:
    from gex_runtime import seal, run_parse_counting


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
    run_parse_counting(config, config_path, args.library, args.threads, run_counting)


if __name__ == "__main__":
    main()

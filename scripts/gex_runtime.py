"""Parse initialization and completion contracts for shared GEX counters."""
from __future__ import annotations

import json
import runpy
import sys

try:
    from scripts.gex_alignment import (atomic_json, backend_options, check_output_backend, check_counting_settings,
                                      completion_name, library_output_dir, resolve_path, sha256_file)
    from scripts.splitpipe_gex_wrapper import validate_splitpipe_tree
except ModuleNotFoundError:
    from gex_alignment import (atomic_json, backend_options, check_output_backend, check_counting_settings,
                               completion_name, library_output_dir, resolve_path, sha256_file)
    from splitpipe_gex_wrapper import validate_splitpipe_tree


def seal(config, config_path, library):
    options = backend_options(config)
    backend = options["backend"]
    root = resolve_path(config["paths"]["output_dir"], config_path)
    library_root = library_output_dir(config, config_path, library)
    process = library_root / "process"
    pending = process / f"{backend}.pending.json"
    payload = json.loads(pending.read_text())
    if sha256_file(process / "barcode_data.csv") != payload.pop("barcode_data_sha256"):
        raise ValueError("Parse barcode contract changed during counting")
    payload["preprocessing"] = {name: sha256_file(process / name)
                                for name in ("barcode_data.csv", "run_proc_def.json")}
    if config.get("gex_decontamination", {}).get("enabled", False):
        payload["filter_manifest_sha256"] = sha256_file(
            root / "gex_grna_filter" / library / "rejected_pairs.complete.json")
    names = ["tscp_assignment.csv", "aw_uf_count_matrix.mtx", "aw_uf_cell_metadata.csv",
             "aw_uf_bc_counts.csv", "all_genes.csv", f"{backend}_qc.json"]
    if backend == "kallisto":
        names.append("kallisto/input.json")
    else:
        names.append("pre_align_stats.csv")
    artifacts = [process / name for name in names]
    for entry in config["libraries"]:
        if entry["sublibrary_name"] != library:
            continue
        for directory in sorted((library_root / entry["sample_name"]).glob("DGE_*")):
            artifacts.extend(directory / name for name in
                             ("count_matrix.mtx", "cell_metadata.csv", "features.csv", "anndata.h5ad"))
            artifacts.extend(directory / name for name in ("all_genes.csv", "all_transcripts.csv")
                             if (directory / name).is_file())
    payload["artifacts"] = {path.relative_to(library_root).as_posix(): sha256_file(path) for path in artifacts}
    payload["complete"] = True
    atomic_json(process / completion_name(options), payload)
    pending.unlink()


def run_parse_counting(config, config_path, library, threads, callback):
    """Initialize Parse's reference, sample masks and cell calls for one callback."""
    if threads < 1:
        raise ValueError("threads must be positive")
    if sum(lib["sublibrary_name"] == library for lib in config["libraries"]) != 1:
        raise ValueError("library must identify exactly one configured library")
    options = backend_options(config)
    output = library_output_dir(config, config_path, library)
    check_output_backend(options, output / "process")
    check_counting_settings(config, output / "process")
    splitpipe_dir = resolve_path(config["paths"]["splitpipe_dir"], config_path)
    validate_splitpipe_tree(splitpipe_dir)
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(splitpipe_dir))
    from splitpipe import dge

    original = dge.run_dge
    called = False

    def count(spipe):
        nonlocal called
        called = True
        return callback(spipe, config, config_path, library, threads)

    dge.run_dge = count
    binary = splitpipe_dir / "split-pipe"
    previous_argv = sys.argv
    sys.argv = [str(binary), "--mode", "dge", "--one_step", "--output_dir", str(output),
                "--nthreads", str(threads), "--kit_score_skip"]
    process = output / "process"
    (process / completion_name(options)).unlink(missing_ok=True)
    (process / f'{options["backend"]}.pending.json').unlink(missing_ok=True)
    try:
        try:
            runpy.run_path(str(binary), run_name="__main__")
        except SystemExit as exc:
            if exc.code not in (0, None):
                raise
        if not called:
            raise RuntimeError("Parse did not invoke the GEX counting callback")
        seal(config, config_path, library)
    finally:
        dge.run_dge = original
        sys.argv = previous_argv

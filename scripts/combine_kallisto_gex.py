#!/usr/bin/env python3
"""Step 1b sparse combine for shared STAR/kallisto gene/transcript products.

Matrix/metadata contracts match the pipeline's DGE consumer. Floating point
counts and transcript identities are preserved throughout.
"""

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import io as sio, sparse

try:
    from scripts.gex_alignment import (SCHEMA, atomic_json, backend_options, resolve_path,
                                       settings_identity, sha256_file, validate_completion,
                                       uses_shared_counter, completion_name, library_output_dir)
    from scripts.gex_counting import write_product
except ModuleNotFoundError:
    from gex_alignment import (SCHEMA, atomic_json, backend_options, resolve_path,
                               settings_identity, sha256_file, validate_completion,
                               uses_shared_counter, completion_name, library_output_dir)
    from gex_counting import write_product


def combine(config, config_path):
    options = backend_options(config)
    if not uses_shared_counter(options):
        raise ValueError("Shared combine requires kallisto or fractional STAR")
    root = resolve_path(config["paths"]["output_dir"], config_path)
    destination = root / "gex_combined"
    destination.mkdir(parents=True, exist_ok=True)
    completion = destination / completion_name(options)
    completion.unlink(missing_ok=True)
    library_manifests = {}
    sample_libraries = {}
    names = [library["sublibrary_name"] for library in config["libraries"]]
    if len(names) != len(set(names)):
        raise ValueError("Duplicated sublibrary names")
    for library in config["libraries"]:
        name = library["sublibrary_name"]
        # Every library uses the same configured reference. Hash its large
        # index once per combine; still verify every library product in full.
        library_manifests[name] = validate_completion(
            config, config_path, name, verify_reference=not library_manifests)
        sample_libraries.setdefault(library["sample_name"], []).append(name)
    reference_hashes = {payload["reference_sha256"] for payload in library_manifests.values()}
    if len(reference_hashes) != 1:
        raise ValueError("Cannot combine GEX libraries built against different references")
    artifacts = []
    for sample, libraries in sample_libraries.items():
        if options["output_level"] != "both":
            for suffix in ("filtered", "unfiltered"):
                stale = destination / sample / f"DGE_transcript_{suffix}"
                if stale.is_dir():
                    shutil.rmtree(stale)
        directories = [(f"DGE_{suffix}", "transcript" if options["output_level"] == "transcript" else "gene")
                       for suffix in ("filtered", "unfiltered")]
        if options["output_level"] == "both":
            directories += [(f"DGE_transcript_{suffix}", "transcript")
                            for suffix in ("filtered", "unfiltered")]
        for dirname, level in directories:
            matrices, cells = [], []
            features = None
            provenance = None
            for suffix, library in enumerate(libraries, start=1):
                source = library_output_dir(config, config_path, library) / sample / dirname
                current_features = pd.read_csv(source / "features.csv", dtype=str).fillna("")
                if features is None:
                    features = current_features
                    provenance = library_manifests[library]["provenance"]
                elif not features.equals(current_features):
                    raise ValueError(f"Feature order/annotations differ in {library}/{dirname}")
                obs = pd.read_csv(source / "cell_metadata.csv", index_col="bc_wells")
                obs.index = obs.index.astype(str)
                matrix = sparse.csr_matrix(sio.mmread(source / "count_matrix.mtx"))
                if matrix.shape != (len(obs), len(features)) or not obs.index.is_unique:
                    raise ValueError(f"Malformed cell/feature matrix in {source}")
                if not np.isfinite(matrix.data).all() or (matrix.data < 0).any():
                    raise ValueError(f"Invalid counts in {source}")
                if not obs["sample"].eq(sample).all():
                    raise ValueError(f"Unexpected sample annotation in {source}")
                obs.index = pd.Index([f"{barcode}__s{suffix}" for barcode in obs.index], name="bc_wells")
                matrices.append(matrix)
                cells.append(obs)
            combined = sparse.vstack(matrices, format="csr")
            obs = pd.concat(cells)
            if not obs.index.is_unique:
                raise ValueError("Combined barcodes are not unique")
            target = destination / sample / dirname
            write_product(target, combined, features, obs,
                          dict(provenance, feature_level=level, sublibraries=libraries), level)
            artifacts += [target / name for name in ("count_matrix.mtx", "features.csv",
                                                     "cell_metadata.csv", "anndata.h5ad")]
            artifacts.append(target / ("all_genes.csv" if level == "gene" else "all_transcripts.csv"))
    payload = {
        "schema_version": SCHEMA, "complete": True, "backend": options["backend"],
        "settings_sha256": settings_identity(config), "reference_sha256": reference_hashes.pop(),
        "sources": {name: sha256_file(library_output_dir(config, config_path, name) / "process" / completion_name(options))
                    for name in names},
        "artifacts": {path.relative_to(destination).as_posix(): sha256_file(path) for path in artifacts},
    }
    atomic_json(completion, payload)
    return destination


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    config = yaml.safe_load(Path(args.config).read_text())
    print(combine(config, args.config))


if __name__ == "__main__":
    main()

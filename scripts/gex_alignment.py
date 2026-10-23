"""Backend configuration and completion contracts; importing requires no Parse."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
from pathlib import Path

import yaml


SCHEMA = 1


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def backend_options(config):
    section = config.get("gex_alignment", {}) or {}
    if not isinstance(section, dict):
        raise ValueError("gex_alignment must be a mapping")
    backend = section.get("backend", "star")
    if backend not in ("star", "kallisto"):
        raise ValueError("gex_alignment.backend must be star or kallisto")
    correction = section.get("umi_correction", "parse")
    if correction not in ("parse", "directional", "none"):
        raise ValueError("gex_alignment.umi_correction must be parse, directional, or none")
    options = (section.get("kallisto", {}) or {}) if backend == "kallisto" else {}
    if not isinstance(options, dict):
        raise ValueError("gex_alignment.kallisto must be a mapping")
    # Common policies live above the backend-specific tool settings. Accept old
    # kallisto configs, but refuse contradictory policies rather than guessing.
    def policy(key, default):
        if key in section and key in options and section[key] != options[key]:
            raise ValueError(f"Conflicting gex_alignment.{key} and kallisto.{key}")
        return section.get(key, options.get(key, default))

    level = policy("output_level", "gene")
    assignment = policy("gene_assignment", "unique")
    transcript_assignment = policy("transcript_assignment", "fractional")
    if level not in ("gene", "transcript", "both"):
        raise ValueError("gex_alignment.output_level must be gene, transcript, or both")
    if assignment not in ("unique", "fractional"):
        raise ValueError("gex_alignment.gene_assignment must be unique or fractional")
    if transcript_assignment not in ("unique", "fractional"):
        raise ValueError("gex_alignment.transcript_assignment must be unique or fractional")
    count_threads = policy("count_threads", 0)
    if isinstance(count_threads, bool) or not isinstance(count_threads, int) or count_threads < 0:
        raise ValueError("gex_alignment.count_threads must be a nonnegative integer")
    common = {"backend": backend, "umi_correction": correction,
              "output_level": level, "gene_assignment": assignment,
              "transcript_assignment": transcript_assignment, "count_threads": count_threads}
    if backend == "star":
        if level != "gene":
            raise ValueError("STAR supports gene output only; use kallisto for transcript/both")
        star = section.get("star", {}) or {}
        if not isinstance(star, dict):
            raise ValueError("gex_alignment.star must be a mapping")
        maximum = star.get("max_multimaps", 3)
        if isinstance(maximum, bool) or not isinstance(maximum, int) or maximum < 1:
            raise ValueError("star.max_multimaps must be a positive integer")
        return dict(common, max_multimaps=maximum, cell_umi_threshold=None)
    if not options.get("reference_dir"):
        raise ValueError("kallisto.reference_dir is required")
    strand = options.get("strand", "forward")
    if strand not in ("forward", "reverse", "unstranded"):
        raise ValueError("kallisto.strand must be forward, reverse, or unstranded")
    if config.get("gex_decontamination", {}).get("existing_alignment", {}).get("enabled", False):
        raise ValueError("kallisto cannot reuse an existing STAR/BAM alignment; use a separate output_dir")
    threshold = options.get("cell_umi_threshold")
    if threshold is not None and (isinstance(threshold, bool) or not isinstance(threshold, (int, float)) or not math.isfinite(threshold) or threshold <= 0):
        raise ValueError("kallisto.cell_umi_threshold must be positive or null (Parse automatic calling)")
    return {
        **common, "reference_dir": str(options["reference_dir"]),
        "executable": str(options.get("executable", "kallisto")),
        "output_level": level, "gene_assignment": assignment,
        "transcript_assignment": transcript_assignment, "strand": strand,
        "cell_umi_threshold": threshold, "count_threads": count_threads,
    }


def uses_shared_counter(options):
    return options["backend"] == "kallisto" or options.get("gene_assignment") == "fractional"


def completion_name(options):
    return "kallisto.complete.json" if options["backend"] == "kallisto" else "star.complete.json"


def resolve_path(value, config_path):
    path = Path(value).expanduser()
    return path if path.is_absolute() else Path(config_path).resolve().parent / path


def tool_executable(value, config_path):
    candidate = str(resolve_path(value, config_path)) if "/" in value else value
    resolved = shutil.which(candidate)
    if resolved is None:
        raise FileNotFoundError(f"Kallisto executable is unavailable: {candidate}")
    return resolved


def check_output_backend(options, process):
    process = Path(process)
    if options["backend"] == "star" and any((process / name).exists() for name in
            ("kallisto.complete.json", "kallisto.pending.json", "kallisto")):
        raise ValueError("Output contains kallisto products; choose a separate STAR output_dir")
    if options["backend"] == "star" and (process / "aw_uf_count_matrix.mtx").exists():
        shared = process / "star.complete.json"
        if (shared.exists() or (process / "star.pending.json").exists()) != uses_shared_counter(options):
            raise ValueError("Output contains another gene assignment policy; choose a separate output_dir")
        if shared.is_file():
            payload = json.loads(shared.read_text())
            provenance = payload.get("provenance", {})
            recorded = provenance.get("umi_correction")
            if recorded != options.get("umi_correction", "parse"):
                raise ValueError("Output contains another UMI correction policy; choose a separate output_dir")
            if provenance.get("max_multimaps") != options.get("max_multimaps", 3):
                raise ValueError("STAR mapping limit changed; choose a separate output_dir")
            return
        marker = process / "gex_umi_policy.json"
        policy = json.loads(marker.read_text()) if marker.is_file() else {}
        recorded = policy.get("umi_correction", "parse")
        if recorded != options.get("umi_correction", "parse"):
            raise ValueError("Output contains another UMI correction policy; choose a separate output_dir")
        if policy.get("max_multimaps", 3) != options.get("max_multimaps", 3):
            raise ValueError("STAR mapping limit changed; choose a separate output_dir")
    if options["backend"] == "kallisto":
        if any(process.glob("barcode_headAligned*.bam")):
            raise ValueError("Output contains STAR BAMs; choose a separate kallisto output_dir")
        if (process / "aw_uf_count_matrix.mtx").exists() and not (process / "kallisto").exists():
            raise ValueError("Output contains prior GEX counts; choose a separate kallisto output_dir")


def load_reference(directory, verify=True):
    directory = Path(directory)
    manifest_path = directory / "reference.json"
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema_version") != SCHEMA or manifest.get("complete") is not True:
        raise ValueError("Incomplete or unsupported kallisto reference manifest")
    required = {"index.idx", "t2g.txt", "features.tsv"}
    artifacts = manifest.get("artifacts", {})
    if not required <= set(artifacts):
        raise ValueError("Kallisto reference lacks required artifacts")
    for name, digest in artifacts.items():
        path = directory / name
        if path.parent != directory or not path.is_file():
            raise ValueError(f"Missing/invalid kallisto reference artifact: {name}")
        if verify and sha256_file(path) != digest:
            raise ValueError(f"Kallisto reference checksum mismatch: {name}")
    return manifest


def settings_identity(config):
    # Include sample ordering, barcode settings, guide policy and reference path.
    # Path-independent output fields are omitted so a run can be relocated.
    keys = ("gex_alignment", "splitpipe", "libraries", "gex_decontamination")
    payload = {key: config.get(key) for key in keys}
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def check_counting_settings(config, process):
    """Bind shared-counter runs before Parse can skip existing intermediates."""
    options = backend_options(config)
    if not uses_shared_counter(options):
        return
    process = Path(process)
    identity = settings_identity(config)
    paths = [process / "gex_counting.settings.json", process / completion_name(options),
             process / f'{options["backend"]}.pending.json']
    for path in paths:
        if path.is_file() and json.loads(path.read_text()).get("settings_sha256") != identity:
            raise ValueError("GEX counting settings changed; choose a separate output_dir")
    atomic_json(paths[0], {"schema_version": SCHEMA, "settings_sha256": identity})


def library_output_dir(config, config_path, library):
    root = resolve_path(config["paths"]["output_dir"], config_path)
    existing = (config.get("gex_decontamination") or {}).get("existing_alignment") or {}
    if existing.get("enabled", False):
        template = existing.get("sandbox_output_dir_template", "{output_dir}/gex_{library}_v1_4")
        return resolve_path(template.format(output_dir=root, library=library), config_path)
    return root / f"gex_{library}"


def _policy_payload(config):
    try:
        from scripts.kallisto_umi import umi_policy
    except ModuleNotFoundError:
        from kallisto_umi import umi_policy
    options = backend_options(config)
    if options["backend"] != "star" or uses_shared_counter(options):
        raise ValueError("STAR UMI policy requires the STAR backend")
    flavor = options["umi_correction"]
    return {"schema_version": SCHEMA, "complete": True, "backend": "star",
            "gene_assignment": "unique", "max_multimaps": options["max_multimaps"],
            "umi_correction": flavor, "umi_policy": umi_policy(flavor)}


def _check_star_policy(payload, config):
    for key, expected in _policy_payload(config).items():
        legacy_defaults = {"gene_assignment": "unique", "max_multimaps": 3}
        if payload.get(key, legacy_defaults.get(key)) != expected:
            raise ValueError("GEX UMI correction policy changed; rerun Step 1 with a separate output_dir")


def write_star_policy(config, directory):
    process = Path(directory) / "process"
    payload = _policy_payload(config)
    assignment = "tscp_assignment.csv" if (process / "tscp_assignment.csv").is_file() else "tscp_assignment.csv.gz"
    # The existing combine wrapper may rename samples in run_proc_def.json.
    # Bind immutable molecule/barcode products without preventing that handoff.
    payload["artifacts"] = {name: sha256_file(process / name) for name in
                            ("aw_uf_count_matrix.mtx", assignment, "barcode_data.csv")}
    atomic_json(process / "gex_umi_policy.json", payload)


def validate_star_policies(config, config_path):
    if uses_shared_counter(backend_options(config)):
        raise ValueError("Fractional STAR requires shared completion validation")
    sources = {}
    for library in config["libraries"]:
        name = library["sublibrary_name"]
        process = library_output_dir(config, config_path, name) / "process"
        marker = process / "gex_umi_policy.json"
        if not marker.is_file():
            options = backend_options(config)
            if options["umi_correction"] != "parse" or options["max_multimaps"] != 3:
                raise ValueError(f"Missing UMI correction provenance for {name}; rerun Step 1a")
            sources[name] = None  # Historical Parse-default output remains usable.
            continue
        payload = json.loads(marker.read_text())
        _check_star_policy(payload, config)
        if not payload.get("artifacts"):
            raise ValueError(f"Empty UMI correction provenance for {name}")
        for relative, digest in payload["artifacts"].items():
            if Path(relative).name != relative or sha256_file(process / relative) != digest:
                raise ValueError(f"GEX molecule output changed for {name}: {relative}")
        sources[name] = sha256_file(marker)
    return sources


def write_star_combined_policy(config, config_path):
    destination = resolve_path(config["paths"]["output_dir"], config_path) / "gex_combined"
    payload = _policy_payload(config)
    payload["sources"] = validate_star_policies(config, config_path)
    payload["artifacts"] = {}
    for sample in sorted({lib["sample_name"] for lib in config["libraries"]}):
        for name in ("anndata.h5ad", "count_matrix.mtx", "cell_metadata.csv"):
            relative = f"{sample}/DGE_filtered/{name}"
            payload["artifacts"][relative] = sha256_file(destination / relative)
    atomic_json(destination / "gex_umi_policy.json", payload)


def validate_star_combined_policy(config, config_path):
    destination = resolve_path(config["paths"]["output_dir"], config_path) / "gex_combined"
    marker = destination / "gex_umi_policy.json"
    sources = validate_star_policies(config, config_path)
    if not marker.is_file():
        if backend_options(config)["umi_correction"] != "parse" or any(sources.values()):
            raise ValueError("Missing combined UMI correction provenance; rerun Step 1b")
        return  # Legacy default-Parse outputs.
    payload = json.loads(marker.read_text())
    _check_star_policy(payload, config)
    if payload.get("sources") != sources or not payload.get("artifacts"):
        raise ValueError("GEX molecule inputs changed since Step 1b; rerun combination")
    for relative, digest in payload["artifacts"].items():
        if Path(relative).is_absolute() or ".." in Path(relative).parts or sha256_file(destination / relative) != digest:
            raise ValueError(f"Combined GEX output changed: {relative}")


def validate_completion(config, config_path, library, *, verify_outputs=True, verify_reference=True):
    options = backend_options(config)
    root = resolve_path(config["paths"]["output_dir"], config_path)
    process = library_output_dir(config, config_path, library) / "process"
    path = process / completion_name(options)
    payload = json.loads(path.read_text())
    if not uses_shared_counter(options) or payload.get("backend") != options["backend"]:
        raise ValueError("Backend mismatch in GEX completion manifest")
    if payload.get("schema_version") != SCHEMA or payload.get("complete") is not True:
        raise ValueError("Incomplete GEX Step 1a")
    if payload.get("library") != library or payload.get("settings_sha256") != settings_identity(config):
        raise ValueError("Stale GEX library/configuration; rerun Step 1a")
    if options["backend"] == "kallisto":
        reference = resolve_path(options["reference_dir"], config_path)
        if sha256_file(reference / "reference.json") != payload.get("reference_sha256"):
            raise ValueError("Stale kallisto reference; rerun Step 1a")
        load_reference(reference, verify=verify_outputs and verify_reference)
    else:
        if star_reference_identity(config, config_path) != payload.get("reference_sha256"):
            raise ValueError("Stale STAR annotation reference; rerun Step 1a")
        if not payload.get("inputs"):
            raise ValueError("Missing STAR alignment provenance")
        for name, digest in payload["inputs"].items():
            if Path(name).name != name or (verify_outputs and sha256_file(process / name) != digest):
                raise ValueError(f"Changed STAR alignment input: {name}")
    for name, digest in payload.get("preprocessing", {}).items():
        if sha256_file(process / name) != digest:
            raise ValueError(f"Changed Parse preprocessing artifact: {name}")
    configured_genome = config.get("paths", {}).get("genome_dir")
    if configured_genome:
        run_info = json.loads((process / "run_proc_def.json").read_text())
        used_genome = run_info.get("genome", {}).get("genome_dir")
        if not used_genome or Path(used_genome).resolve() != resolve_path(configured_genome, config_path).resolve():
            raise ValueError("Parse genome configuration differs from Step 1a preprocessing")
    decontam = config.get("gex_decontamination", {})
    if decontam.get("enabled", False):
        filter_path = root / "gex_grna_filter" / library / "rejected_pairs.complete.json"
        if sha256_file(filter_path) != payload.get("filter_manifest_sha256"):
            raise ValueError("Changed GEX decontamination evidence; rerun Step 1a")
    artifacts = payload.get("artifacts", {})
    if not artifacts:
        raise ValueError("Empty GEX completion artifact set")
    for relative, digest in artifacts.items():
        artifact = process.parent / relative
        if Path(relative).is_absolute() or ".." in Path(relative).parts or not artifact.is_file():
            raise ValueError(f"Missing GEX output: {relative}")
        if verify_outputs and sha256_file(artifact) != digest:
            raise ValueError(f"Changed GEX output: {relative}")
    return payload


def validate_combined(config, config_path, *, verify_outputs=True):
    options = backend_options(config)
    root = resolve_path(config["paths"]["output_dir"], config_path)
    destination = root / "gex_combined"
    payload = json.loads((destination / completion_name(options)).read_text())
    if payload.get("schema_version") != SCHEMA or payload.get("complete") is not True or payload.get("backend") != options["backend"]:
        raise ValueError("Incomplete GEX Step 1b")
    if payload.get("settings_sha256") != settings_identity(config):
        raise ValueError("Stale combined GEX configuration")
    expected = {lib["sublibrary_name"] for lib in config["libraries"]}
    if set(payload.get("sources", {})) != expected:
        raise ValueError("Combined GEX library set is stale")
    for name, digest in payload["sources"].items():
        path = library_output_dir(config, config_path, name) / "process" / completion_name(options)
        if sha256_file(path) != digest:
            raise ValueError(f"GEX Step 1a changed since combine: {name}")
        validate_completion(config, config_path, name, verify_outputs=False)
    if not payload.get("artifacts"):
        raise ValueError("Empty GEX Step 1b artifacts")
    for relative, digest in payload["artifacts"].items():
        path = destination / relative
        if Path(relative).is_absolute() or ".." in Path(relative).parts or not path.is_file():
            raise ValueError(f"Missing combined GEX artifact: {relative}")
        if verify_outputs and sha256_file(path) != digest:
            raise ValueError(f"Changed combined GEX artifact: {relative}")
    return payload


def star_reference_identity(config, config_path):
    reference = resolve_path(config["paths"]["genome_dir"], config_path)
    artifacts = {name: sha256_file(reference / name) for name in ("gene_info.json", "all_genes.csv")}
    return hashlib.sha256(json.dumps(artifacts, sort_keys=True).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--library")
    parser.add_argument("--field", choices=("backend", "gene_assignment", "max_multimaps", "shared_counter"))
    parser.add_argument("--check-output", action="store_true")
    parser.add_argument("--write-star-policy", metavar="LIBRARY_OUTPUT_DIR")
    parser.add_argument("--validate-star-policies", action="store_true")
    parser.add_argument("--write-star-combined-policy", action="store_true")
    args = parser.parse_args()
    config = yaml.safe_load(Path(args.config).read_text())
    options = backend_options(config)
    if args.field:
        print(str(uses_shared_counter(options)).lower() if args.field == "shared_counter" else options[args.field])
        return
    if args.write_star_policy or args.validate_star_policies or args.write_star_combined_policy:
        if options["backend"] != "star":
            parser.error("STAR UMI policy operations require the STAR backend")
        if args.write_star_policy:
            write_star_policy(config, resolve_path(args.write_star_policy, args.config))
        elif args.write_star_combined_policy:
            write_star_combined_policy(config, args.config)
        else:
            validate_star_policies(config, args.config)
    if args.check_output:
        if not args.library:
            parser.error("--check-output requires --library")
        process = library_output_dir(config, args.config, args.library) / "process"
        check_output_backend(options, process)
        check_counting_settings(config, process)
        if options["backend"] == "kallisto":
            load_reference(resolve_path(options["reference_dir"], args.config), verify=False)
            tool_executable(options["executable"], args.config)
    print(options["backend"])


if __name__ == "__main__":
    main()

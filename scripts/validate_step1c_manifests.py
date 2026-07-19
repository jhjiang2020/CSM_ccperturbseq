#!/usr/bin/env python3
"""Validate that every configured library completed Step 1c under this config."""

import argparse
import hashlib
import json
import os
from pathlib import Path

import yaml


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate(config_path):
    config_path = Path(config_path).expanduser().resolve()
    with config_path.open() as handle:
        config = yaml.safe_load(handle) or {}
    decontam = config.get("gex_decontamination", {}) or {}
    step1c = decontam.get("step1c", {}) or {}
    if not step1c.get("enabled", False):
        raise ValueError("gex_decontamination.step1c.enabled is false")

    output_root = Path(config["paths"]["output_dir"]).expanduser()
    if not output_root.is_absolute():
        output_root = config_path.parent / output_root
    namespace = Path(step1c.get("output_namespace", "gex_grna_identity"))
    if not namespace.is_absolute():
        namespace = output_root / namespace
    expected_config_hash = sha256_file(config_path)
    errors = []
    manifests = []
    for library in config.get("libraries", []):
        name = library["sublibrary_name"]
        path = namespace / "manifests" / f"{name}.step1c.complete.json"
        if not path.is_file():
            errors.append(f"{name}: missing {path}")
            continue
        try:
            payload = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            errors.append(f"{name}: cannot read completion manifest: {exc}")
            continue
        if payload.get("complete") is not True or payload.get("status") != "complete":
            errors.append(f"{name}: Step 1c manifest is incomplete")
        if payload.get("schema_version") != 2:
            errors.append(f"{name}: unsupported Step 1c manifest schema")
        if payload.get("source_library") != name:
            errors.append(f"{name}: manifest library mismatch")
        if payload.get("config_sha256") != expected_config_hash:
            errors.append(f"{name}: Step 1c was produced from a stale config")
        for key in ("barcode_output", "match_output"):
            raw = payload.get(key)
            if not raw or not Path(raw).is_file():
                errors.append(f"{name}: missing {key}")
                continue
            expected = payload.get(f"{key}_sha256")
            if not expected or expected != sha256_file(raw):
                errors.append(f"{name}: {key} SHA-256 mismatch")
        strict = payload.get("strict_matching", {}) or {}
        if step1c.get("strict_matching", {}).get("enabled", True):
            if not strict.get("enabled"):
                errors.append(f"{name}: strict Step 1c product was not enabled")
            raw_match = payload.get("unfiltered_match_output")
            if not raw_match or not Path(raw_match).is_file():
                errors.append(f"{name}: unfiltered Step 3 audit product is missing")
            elif payload.get("unfiltered_match_output_sha256") != sha256_file(raw_match):
                errors.append(f"{name}: unfiltered Step 3 SHA-256 mismatch")
        step1a_raw = payload.get("step1a_manifest")
        if not step1a_raw or not Path(step1a_raw).is_file():
            errors.append(f"{name}: missing Step 1a source manifest")
        elif payload.get("step1a_manifest_sha256") != sha256_file(step1a_raw):
            errors.append(f"{name}: Step 1a source manifest changed after Step 1c")
        manifests.append(path)
    if errors:
        raise ValueError("Step 1c validation failed:\n  " + "\n  ".join(errors))
    return manifests


def write_audit_summary(config_path, manifests):
    """Collate Step 1a removal and Step 1c identity diagnostics."""
    config_path = Path(config_path).expanduser().resolve()
    with config_path.open() as handle:
        config = yaml.safe_load(handle) or {}
    output_root = Path(config["paths"]["output_dir"]).expanduser()
    if not output_root.is_absolute():
        output_root = config_path.parent / output_root
    step1c = config["gex_decontamination"]["step1c"]
    namespace = Path(step1c.get("output_namespace", "gex_grna_identity"))
    if not namespace.is_absolute():
        namespace = output_root / namespace

    totals = {
        "filter_reads_examined": 0,
        "tso_candidate_reads": 0,
        "filter_positive_reads": 0,
        "valid_barcode_reads_removed_from_star": 0,
        "step1c_rejected_raw_pairs": 0,
    }
    libraries = {}
    key_map = {
        "filter_reads_examined": "csm_gex_filter_reads_examined",
        "tso_candidate_reads": "csm_gex_tso_candidate_reads",
        "filter_positive_reads": "csm_gex_filter_positive_reads",
        "valid_barcode_reads_removed_from_star": (
            "csm_gex_filter_positive_valid_barcode_reads"
        ),
        "step1c_rejected_raw_pairs": "csm_gex_step1c_rejected_raw_pairs",
    }
    for manifest_path in manifests:
        step1c_payload = json.loads(Path(manifest_path).read_text())
        step1a_path = Path(step1c_payload["step1a_manifest"])
        step1a_payload = json.loads(step1a_path.read_text())
        filter_stats = step1a_payload.get("statistics", {}) or {}
        for total_name, source_name in key_map.items():
            totals[total_name] += int(filter_stats.get(source_name, 0))
        libraries[step1c_payload["source_library"]] = {
            "step1a_filter_statistics": filter_stats,
            "step1a_filter_rates": step1a_payload.get("rates", {}) or {},
            "step1c_identity_statistics": step1c_payload.get("statistics", {}) or {},
            "step1c_strict_matching": step1c_payload.get("strict_matching", {}) or {},
            "step1a_manifest": str(step1a_path),
            "step1c_manifest": str(manifest_path),
        }
    totals["filter_positive_fraction_of_examined_reads"] = (
        totals["filter_positive_reads"] / totals["filter_reads_examined"]
        if totals["filter_reads_examined"] else 0.0
    )
    totals["filter_positive_fraction_of_tso_candidates"] = (
        totals["filter_positive_reads"] / totals["tso_candidate_reads"]
        if totals["tso_candidate_reads"] else 0.0
    )
    payload = {
        "pipeline_version": "1.4",
        "totals": totals,
        "libraries": libraries,
    }
    destination = namespace / "gex_grna_decontamination_summary.yaml"
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp.{os.getpid()}")
    with temporary.open("w") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, destination)
    return destination


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    try:
        manifests = validate(args.config)
        summary = write_audit_summary(args.config, manifests)
    except (OSError, KeyError, TypeError, ValueError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")
    print(f"Validated {len(manifests)} Step 1c library manifests")
    print(f"Step 1a/1c audit summary: {summary}")


if __name__ == "__main__":
    main()

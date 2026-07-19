#!/usr/bin/env python3
"""Validate Step 1a GEX decontamination manifests before combine/resume."""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import yaml

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from gex_grna_filter import filter_settings_sha256, load_filter_from_config  # noqa: E402

EXPECTED_FILTER_SCHEMA = 4
EXPECTED_SURVIVOR_TRIM_POLICY = "full_tso_endpoint_then_parse_exact_17mer"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate(config_path, only_library=None, verify_bam_hashes=True):
    config_path = Path(config_path).expanduser().resolve()
    with config_path.open() as handle:
        config = yaml.safe_load(handle) or {}
    decontamination = config.get("gex_decontamination", {}) or {}
    if not bool(decontamination.get("enabled", False)):
        return []

    classifier = load_filter_from_config(config_path)
    expected_hash = classifier.index.stats.reference_sha256
    expected_settings_hash = filter_settings_sha256(classifier.settings)
    output_root = Path(config["paths"]["output_dir"]).expanduser()
    existing = decontamination.get("existing_alignment", {}) or {}
    errors = []
    manifests = []
    configured_libraries = config.get("libraries", [])
    if only_library is not None:
        configured_libraries = [
            library
            for library in configured_libraries
            if library.get("sublibrary_name") == only_library
        ]
        if len(configured_libraries) != 1:
            raise ValueError(
                f"expected one configured library named {only_library!r}, "
                f"found {len(configured_libraries)}"
            )
    for library in configured_libraries:
        name = library["sublibrary_name"]
        path = output_root / "gex_grna_filter" / name / "rejected_pairs.complete.json"
        if not path.is_file():
            errors.append(f"{name}: missing {path}")
            continue
        try:
            payload = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            errors.append(f"{name}: cannot read manifest: {exc}")
            continue
        if payload.get("complete") is not True or payload.get("status") != "complete":
            errors.append(f"{name}: manifest is not complete")
        if payload.get("library") != name:
            errors.append(f"{name}: manifest library is {payload.get('library')!r}")
        if payload.get("schema_version") != EXPECTED_FILTER_SCHEMA:
            errors.append(f"{name}: filter manifest schema is stale")
        if payload.get("survivor_tso_trim_policy") != EXPECTED_SURVIVOR_TRIM_POLICY:
            errors.append(f"{name}: survivor TSO trim policy is stale")
        if payload.get("reference_sha256") != expected_hash:
            errors.append(f"{name}: guide reference hash is stale")
        if payload.get("filter_backend") != "packed_inline":
            errors.append(f"{name}: unsupported filter backend in manifest")
        if payload.get("filter_settings_sha256") != expected_settings_hash:
            errors.append(f"{name}: removal-policy settings are stale")
        if payload.get("write_rejected_raw_pairs"):
            for key in ("rejected_r1", "rejected_r2"):
                candidate = payload.get(key)
                if not candidate or not Path(candidate).is_file():
                    errors.append(f"{name}: missing rejected FASTQ {key}")
        if existing.get("enabled", False):
            template = existing.get(
                "sandbox_output_dir_template", "{output_dir}/gex_{library}_v1_4"
            )
            run_root = Path(
                template.format(output_dir=output_root, library=name)
            )
        else:
            run_root = output_root / f"gex_{name}"
        harmonized_path = (
            run_root / "process" / "csm_harmonized_bam.complete.json"
        )
        if not harmonized_path.is_file():
            errors.append(f"{name}: missing harmonized BAM manifest {harmonized_path}")
        else:
            try:
                harmonized = json.loads(harmonized_path.read_text())
            except (OSError, json.JSONDecodeError) as exc:
                errors.append(f"{name}: cannot read harmonized BAM manifest: {exc}")
            else:
                if harmonized.get("complete") is not True:
                    errors.append(f"{name}: harmonized BAM manifest is incomplete")
                if harmonized.get("library") != name:
                    errors.append(f"{name}: harmonized BAM library mismatch")
                if harmonized.get("filter_reference_sha256") != expected_hash:
                    errors.append(f"{name}: harmonized BAM guide reference is stale")
                if (
                    harmonized.get("filter_settings_sha256")
                    != expected_settings_hash
                ):
                    errors.append(f"{name}: harmonized BAM removal policy is stale")
                expected_trim_policy = (
                    "not_applied_existing_alignment"
                    if existing.get("enabled", False)
                    else EXPECTED_SURVIVOR_TRIM_POLICY
                )
                if harmonized.get("survivor_tso_trim_policy") != expected_trim_policy:
                    errors.append(f"{name}: harmonized BAM survivor trim policy is stale")
                bam = harmonized.get("output_bam")
                if not bam or not Path(bam).is_file():
                    errors.append(f"{name}: harmonized BAM is missing")
                elif harmonized.get("alignment_reused") is True:
                    source_value = harmonized.get("source_bam")
                    source = Path(source_value) if source_value else None
                    output_bam = Path(bam)
                    if source is None or not source.is_file():
                        errors.append(f"{name}: historical source BAM is missing")
                    elif source.resolve() == output_bam.resolve():
                        errors.append(f"{name}: source and output BAM are in-place")
                    else:
                        expected_source_sha = harmonized.get("source_bam_sha256")
                        if verify_bam_hashes and sha256_file(source) != expected_source_sha:
                            errors.append(
                                f"{name}: production source BAM SHA-256 changed"
                            )
                    expected_output_sha = harmonized.get("output_bam_sha256")
                    if verify_bam_hashes and sha256_file(output_bam) != expected_output_sha:
                        errors.append(
                            f"{name}: harmonized output BAM SHA-256 changed after mol"
                        )
                    if harmonized.get("name_sort_validated") is not True:
                        errors.append(f"{name}: query-name validation is missing")
                    try:
                        records_before = int(harmonized["records_before"])
                        records_removed = int(harmonized["records_removed"])
                        records_after = int(harmonized["records_after"])
                        qnames_before = int(harmonized["qnames_before"])
                        qnames_removed = int(harmonized["qnames_removed"])
                        qnames_after = int(harmonized["qnames_after"])
                        original_headers = int(
                            harmonized["qnames_with_original_header"]
                        )
                        positives_seen = int(harmonized["positive_qnames_seen"])
                        positives_missing = int(
                            harmonized["positive_ids_missing_from_bam"]
                        )
                        valid_positives = int(
                            payload["statistics"][
                                "csm_gex_filter_positive_valid_barcode_reads"
                            ]
                        )
                    except (KeyError, TypeError, ValueError):
                        errors.append(
                            f"{name}: historical BAM count contract is incomplete"
                        )
                    else:
                        if records_before != records_removed + records_after:
                            errors.append(f"{name}: BAM record invariant failed")
                        if qnames_before != qnames_removed + qnames_after:
                            errors.append(f"{name}: BAM QNAME invariant failed")
                        if records_before and original_headers <= 0:
                            errors.append(
                                f"{name}: no original-header QNAMEs were validated"
                            )
                        if positives_missing != 0:
                            errors.append(
                                f"{name}: filter-positive IDs are missing from source BAM"
                            )
                        if not (
                            positives_seen == qnames_removed == valid_positives
                        ):
                            errors.append(
                                f"{name}: valid-barcode removal is not exact: "
                                f"seen={positives_seen}, removed={qnames_removed}, "
                                f"expected={valid_positives}"
                            )
        manifests.append(path)
    if errors:
        raise ValueError("Step 1a GEX filter validation failed:\n  " + "\n  ".join(errors))
    return manifests


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--library")
    args = parser.parse_args()
    try:
        manifests = validate(args.config, only_library=args.library)
    except (OSError, KeyError, TypeError, ValueError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")
    if manifests:
        print(f"Validated {len(manifests)} Step 1a GEX filter manifests")
    else:
        print("GEX decontamination disabled; no filter manifests required")


if __name__ == "__main__":
    main()

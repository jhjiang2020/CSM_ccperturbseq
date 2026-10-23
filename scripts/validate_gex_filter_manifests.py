#!/usr/bin/env python3
"""Validate Step 1a GEX decontamination manifests before combine/resume."""

import argparse
import csv
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
EXPECTED_HARMONIZED_SCHEMA = 2
EXPECTED_ALIGNMENT_STATS_SCHEMA = 1
EXPECTED_SURVIVOR_TRIM_POLICY = "full_tso_endpoint_then_parse_exact_17mer"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_stats_values(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or not {"statistic", "value"}.issubset(
            reader.fieldnames
        ):
            raise ValueError(f"invalid split-pipe statistics header: {path}")
        values = {}
        for row in reader:
            key = str(row.get("statistic", "")).strip()
            if not key or key in values:
                raise ValueError(f"invalid or duplicate statistic {key!r}: {path}")
            values[key] = str(row.get("value", "")).strip()
    return values


def required_int(mapping, key, context):
    try:
        return int(mapping[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"missing/invalid {key!r} in {context}") from exc


def validate_existing_alignment_contract(
    name, filter_payload, harmonized, run_root, *, verify_bam_hashes, errors
):
    if harmonized.get("schema_version") != EXPECTED_HARMONIZED_SCHEMA:
        errors.append(f"{name}: harmonized BAM manifest schema is stale")

    source_value = harmonized.get("source_bam")
    output_value = harmonized.get("output_bam")
    source = Path(source_value) if source_value else None
    output_bam = Path(output_value) if output_value else None
    if source is None or not source.is_file():
        errors.append(f"{name}: historical source BAM is missing")
    elif output_bam is not None and source.resolve() == output_bam.resolve():
        errors.append(f"{name}: source and output BAM are in-place")
    else:
        expected_source_sha = harmonized.get("source_bam_sha256")
        if verify_bam_hashes and sha256_file(source) != expected_source_sha:
            errors.append(f"{name}: production source BAM SHA-256 changed")
    if output_bam is None or not output_bam.is_file():
        errors.append(f"{name}: harmonized output BAM is missing")
    else:
        expected_output_sha = harmonized.get("output_bam_sha256")
        if verify_bam_hashes and sha256_file(output_bam) != expected_output_sha:
            errors.append(f"{name}: harmonized output BAM SHA-256 changed after mol")
    if harmonized.get("name_sort_validated") is not True:
        errors.append(f"{name}: query-name validation is missing")

    try:
        records_before = required_int(harmonized, "records_before", "BAM manifest")
        records_removed = required_int(harmonized, "records_removed", "BAM manifest")
        records_after = required_int(harmonized, "records_after", "BAM manifest")
        qnames_before = required_int(harmonized, "qnames_before", "BAM manifest")
        qnames_removed = required_int(harmonized, "qnames_removed", "BAM manifest")
        qnames_after = required_int(harmonized, "qnames_after", "BAM manifest")
        original_headers = required_int(
            harmonized, "qnames_with_original_header", "BAM manifest"
        )
        positive_records = required_int(
            harmonized, "positive_id_records", "BAM manifest"
        )
        positive_unique = required_int(
            harmonized, "positive_unique_ids", "BAM manifest"
        )
        positives_seen = required_int(
            harmonized, "positive_qnames_seen", "BAM manifest"
        )
        positives_missing = required_int(
            harmonized, "positive_ids_missing_from_bam", "BAM manifest"
        )
        filter_positives = required_int(
            filter_payload.get("statistics", {}),
            "csm_gex_filter_positive_reads",
            "filter manifest",
        )
        valid_positives = required_int(
            filter_payload.get("statistics", {}),
            "csm_gex_filter_positive_valid_barcode_reads",
            "filter manifest",
        )
    except ValueError as exc:
        errors.append(f"{name}: historical BAM count contract is incomplete: {exc}")
        return

    if records_before != records_removed + records_after:
        errors.append(f"{name}: BAM record invariant failed")
    if qnames_before != qnames_removed + qnames_after:
        errors.append(f"{name}: BAM QNAME invariant failed")
    if records_before and original_headers <= 0:
        errors.append(f"{name}: no original-header QNAMEs were validated")
    if positive_records != filter_positives or positive_unique != filter_positives:
        errors.append(
            f"{name}: positive evidence count differs from the filter manifest"
        )
    if positives_missing != positive_unique - positives_seen:
        errors.append(f"{name}: missing-positive QNAME arithmetic failed")
    if positives_seen != qnames_removed:
        errors.append(
            f"{name}: removed QNAMEs do not equal observed positive QNAMEs"
        )
    if positives_seen > valid_positives:
        errors.append(
            f"{name}: aligned positive QNAMEs exceed valid-barcode positives"
        )
    if harmonized.get("positive_valid_barcode_records") != valid_positives:
        errors.append(f"{name}: valid-positive provenance count is inconsistent")

    handoff = harmonized.get("alignment_stats_contract")
    if not isinstance(handoff, dict):
        errors.append(f"{name}: mol alignment-statistics handoff is missing")
        return
    if handoff.get("schema_version") != EXPECTED_ALIGNMENT_STATS_SCHEMA:
        errors.append(f"{name}: mol alignment-statistics schema is stale")
    try:
        source_unique = required_int(
            handoff, "source_reads_align_unique", "alignment-statistics handoff"
        )
        source_multimap = required_int(
            handoff, "source_reads_align_multimap", "alignment-statistics handoff"
        )
        removed_unique = required_int(
            handoff, "removed_unique_qnames", "alignment-statistics handoff"
        )
        removed_multimap = required_int(
            handoff, "removed_multimap_qnames", "alignment-statistics handoff"
        )
        filtered_unique = required_int(
            handoff, "filtered_reads_align_unique", "alignment-statistics handoff"
        )
        filtered_multimap = required_int(
            handoff, "filtered_reads_align_multimap", "alignment-statistics handoff"
        )
        align_input = required_int(
            handoff, "source_reads_align_input", "alignment-statistics handoff"
        )
        number_of_reads = required_int(
            handoff, "number_of_reads", "alignment-statistics handoff"
        )
        source_valid = required_int(
            handoff, "source_reads_valid_barcode", "alignment-statistics handoff"
        )
        sandbox_valid = required_int(
            handoff, "sandbox_reads_valid_barcode", "alignment-statistics handoff"
        )
        handoff_valid_positives = required_int(
            handoff, "valid_positive_records", "alignment-statistics handoff"
        )
    except ValueError as exc:
        errors.append(f"{name}: mol alignment-statistics contract is incomplete: {exc}")
        return

    if source_unique + source_multimap != qnames_before:
        errors.append(f"{name}: source alignment counts do not match source BAM")
    if removed_unique + removed_multimap != qnames_removed:
        errors.append(f"{name}: removed alignment counts do not match removed QNAMEs")
    if filtered_unique != source_unique - removed_unique:
        errors.append(f"{name}: filtered unique-alignment count is inconsistent")
    if filtered_multimap != source_multimap - removed_multimap:
        errors.append(f"{name}: filtered multimap count is inconsistent")
    if filtered_unique + filtered_multimap != qnames_after:
        errors.append(f"{name}: filtered alignment counts do not match output BAM")
    if handoff_valid_positives != valid_positives:
        errors.append(f"{name}: handoff valid-positive count is inconsistent")
    if source_valid - sandbox_valid != valid_positives:
        errors.append(f"{name}: handoff valid-barcode count difference is inconsistent")
    if align_input != source_valid:
        errors.append(f"{name}: source STAR input differs from valid-barcode reads")

    source_stats_value = handoff.get("source_pre_align_stats")
    source_stats = Path(source_stats_value) if source_stats_value else None
    if source_stats is None or not source_stats.is_file():
        errors.append(f"{name}: historical pre-alignment statistics are missing")
    else:
        expected_stats_sha = handoff.get("source_pre_align_stats_sha256")
        if sha256_file(source_stats) != expected_stats_sha:
            errors.append(f"{name}: historical pre-alignment statistics changed")
        if source is not None and source_stats.resolve() != (
            source.parent / "pre_align_stats.csv"
        ).resolve():
            errors.append(f"{name}: historical BAM/statistics paths are inconsistent")
        try:
            historical_stats = read_stats_values(source_stats)
            historical_reads = required_int(
                historical_stats, "number_of_reads", source_stats
            )
            historical_valid = required_int(
                historical_stats, "reads_valid_barcode", source_stats
            )
            historical_input = required_int(
                historical_stats, "reads_align_input", source_stats
            )
            historical_unique = required_int(
                historical_stats, "reads_align_unique", source_stats
            )
            historical_multimap = required_int(
                historical_stats, "reads_align_multimap", source_stats
            )
        except (OSError, ValueError) as exc:
            errors.append(f"{name}: cannot validate historical alignment statistics: {exc}")
        else:
            if (
                historical_reads != number_of_reads
                or historical_valid != source_valid
                or historical_input != align_input
                or historical_unique != source_unique
                or historical_multimap != source_multimap
            ):
                errors.append(f"{name}: historical alignment-statistics values changed")

    sandbox_stats = run_root / "process" / "pre_align_stats.csv"
    recorded_sandbox = handoff.get("sandbox_pre_align_stats")
    if not sandbox_stats.is_file() or not recorded_sandbox:
        errors.append(f"{name}: sandbox mol statistics are missing")
    elif Path(recorded_sandbox).resolve() != sandbox_stats.resolve():
        errors.append(f"{name}: sandbox mol-statistics path is inconsistent")
    else:
        try:
            current_stats = read_stats_values(sandbox_stats)
            current_input = required_int(
                current_stats, "reads_align_input", sandbox_stats
            )
            current_reads = required_int(
                current_stats, "number_of_reads", sandbox_stats
            )
            current_valid = required_int(
                current_stats, "reads_valid_barcode", sandbox_stats
            )
            current_unique = required_int(
                current_stats, "reads_align_unique", sandbox_stats
            )
            current_multimap = required_int(
                current_stats, "reads_align_multimap", sandbox_stats
            )
        except (OSError, ValueError) as exc:
            errors.append(f"{name}: cannot validate sandbox mol statistics: {exc}")
        else:
            if current_input != align_input:
                errors.append(f"{name}: sandbox STAR-input statistic changed")
            if current_reads != number_of_reads or current_valid != sandbox_valid:
                errors.append(f"{name}: sandbox preprocessing statistics changed")
            if current_unique != filtered_unique or current_multimap != filtered_multimap:
                errors.append(f"{name}: sandbox alignment counts changed after handoff")


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
        if config.get("gex_alignment", {}).get("backend", "star") == "kallisto":
            try:
                try:
                    from scripts.gex_alignment import validate_completion
                except ModuleNotFoundError:
                    from gex_alignment import validate_completion
                validate_completion(config, config_path, name, verify_outputs=False)
            except (OSError, ValueError, KeyError) as exc:
                errors.append(f"{name}: native kallisto handoff is invalid: {exc}")
            manifests.append(path)
            continue
        try:
            try:
                from scripts.gex_alignment import backend_options, uses_shared_counter, validate_completion
            except ModuleNotFoundError:
                from gex_alignment import backend_options, uses_shared_counter, validate_completion
            if uses_shared_counter(backend_options(config)):
                validate_completion(config, config_path, name, verify_outputs=False)
        except (OSError, ValueError, KeyError) as exc:
            errors.append(f"{name}: fractional STAR handoff is invalid: {exc}")
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
                if harmonized.get("alignment_reused") is True:
                    validate_existing_alignment_contract(
                        name,
                        payload,
                        harmonized,
                        run_root,
                        verify_bam_hashes=verify_bam_hashes,
                        errors=errors,
                    )
                elif not bam or not Path(bam).is_file():
                    errors.append(f"{name}: harmonized BAM is missing")
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

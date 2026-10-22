#!/usr/bin/env python3
"""Archived validator for the completed barcode-decoder crossover audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import yaml

try:
    from scripts.barcode_crossover_audit import (
        CROSSOVER_PRODUCTS,
        EVALUATION_READS,
        FORBIDDEN_PUBLIC_COLUMNS,
    )
except ImportError:
    from barcode_crossover_audit import (
        CROSSOVER_PRODUCTS,
        EVALUATION_READS,
        FORBIDDEN_PUBLIC_COLUMNS,
    )


DECODERS = {
    "step1a_native",
    "step2_historical",
    "step2_future_contract_aligned",
}
DECODER_PAIRS = {
    ("step1a_native", "step2_historical"),
    ("step1a_native", "step2_future_contract_aligned"),
    ("step2_historical", "step2_future_contract_aligned"),
}


def _as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().map({"true": True, "false": False})


def validate(
    config_path: Path,
    diagnostic_dir: Path,
    asset_dir: Path,
    evaluation_reads: int = EVALUATION_READS,
    maximum_asset_table_bytes: int = 1_000_000,
) -> dict:
    cfg = yaml.safe_load(config_path.read_text())
    libraries = [item["sublibrary_name"] for item in cfg["libraries"]]
    expected_library_arms = {
        (library, arm)
        for library in libraries
        for arm in ("GEX", "gRNA")
    }
    frames = {}
    for name in CROSSOVER_PRODUCTS:
        path = diagnostic_dir / name
        if not path.exists():
            raise FileNotFoundError(path)
        frame = pd.read_csv(path)
        forbidden = FORBIDDEN_PUBLIC_COLUMNS & set(frame.columns)
        if forbidden:
            raise AssertionError(f"{name} exposes forbidden columns: {forbidden}")
        frames[name] = frame

    decoder = frames["decoder_summary.csv"]
    expected_decoder_rows = {
        (library, arm, decoder_name)
        for library, arm in expected_library_arms
        for decoder_name in DECODERS
    }
    observed_decoder_rows = set(
        zip(decoder["sublibrary"], decoder["input_arm"], decoder["decoder"])
    )
    if observed_decoder_rows != expected_decoder_rows:
        raise AssertionError("Decoder summary does not cover every profile")
    if not (decoder["sampled_reads"] == evaluation_reads).all():
        raise AssertionError("Decoder profiles did not use the exact sample size")
    if not (
        decoder["accepted_reads"] <= decoder["valid_unrestricted_reads"]
    ).all():
        raise AssertionError("Accepted calls exceed unrestricted valid calls")
    warnings = _as_bool(decoder["replay_warning_over_3pp"]).fillna(False)
    expected_warnings = (
        decoder["sample_minus_full_percentage_points"].abs() > 3
    ).fillna(False)
    if not (warnings == expected_warnings).all():
        raise AssertionError("Three-percentage-point replay flags are inconsistent")

    decoder_round = frames["decoder_round_summary.csv"]
    decoder_round_partition = decoder_round.groupby(
        ["sublibrary", "input_arm", "decoder", "bc_round"]
    )["read_count"].sum()
    expected_decoder_round = {
        (library, arm, decoder_name, round_number)
        for library, arm in expected_library_arms
        for decoder_name in DECODERS
        for round_number in (1, 2, 3)
    }
    if set(decoder_round_partition.index) != expected_decoder_round:
        raise AssertionError("Decoder-specific per-round products are incomplete")
    if not (decoder_round_partition == evaluation_reads).all():
        raise AssertionError(
            "Decoder-specific per-round categories do not partition reads"
        )

    round_summary = frames["round_mapping_summary.csv"]
    round_partition = round_summary.groupby(
        ["sublibrary", "input_arm", "bc_round"]
    )["read_count"].sum()
    expected_round_rows = {
        (library, arm, round_number)
        for library, arm in expected_library_arms
        for round_number in (1, 2, 3)
    }
    if set(round_partition.index) != expected_round_rows:
        raise AssertionError("Per-round products are incomplete")
    if not (round_partition == evaluation_reads).all():
        raise AssertionError("Per-round mapping categories do not partition reads")

    for name, group_columns in (
        ("bc1_reference_summary.csv", ["sublibrary", "input_arm"]),
        ("failure_mask_summary.csv", ["sublibrary", "input_arm"]),
        ("linker_summary.csv", ["sublibrary", "input_arm"]),
        ("linker_rescue_summary.csv", ["sublibrary", "input_arm"]),
    ):
        partition = frames[name].groupby(group_columns)["read_count"].sum()
        if set(partition.index) != expected_library_arms:
            raise AssertionError(f"{name} does not cover every library/arm")
        if not (partition == evaluation_reads).all():
            raise AssertionError(f"{name} does not partition evaluation reads")
    masks = frames["failure_mask_summary.csv"]
    expected_masks = {
        (library, arm, f"{value:03b}")
        for library, arm in expected_library_arms
        for value in range(8)
    }
    observed_masks = set(
        zip(
            masks["sublibrary"],
            masks["input_arm"],
            masks["bc1_bc2_bc3_unique_mask"].astype(str).str.zfill(3),
        )
    )
    if observed_masks != expected_masks:
        raise AssertionError("All eight BC1/BC2/BC3 masks were not emitted")

    linker_quality = frames["linker_quality_by_status.csv"]
    linker_quality_partition = linker_quality.groupby(
        ["sublibrary", "input_arm", "segment"]
    )["reads"].sum()
    expected_linker_quality = {
        (library, arm, linker)
        for library, arm in expected_library_arms
        for linker in ("L23", "L12")
    }
    if set(linker_quality_partition.index) != expected_linker_quality:
        raise AssertionError("Linker quality summaries are incomplete")
    if not (linker_quality_partition == evaluation_reads).all():
        raise AssertionError("Linker quality statuses do not partition reads")

    concordance = frames["decoder_concordance.csv"]
    concordance_partition = concordance.groupby(
        ["sublibrary", "input_arm", "left_decoder", "right_decoder"]
    )["read_count"].sum()
    expected_concordance = {
        (library, arm, left, right)
        for library, arm in expected_library_arms
        for left, right in DECODER_PAIRS
    }
    if set(concordance_partition.index) != expected_concordance:
        raise AssertionError("Decoder concordance pairs are incomplete")
    if not (concordance_partition == evaluation_reads).all():
        raise AssertionError("Decoder concordance does not partition reads")

    policy = frames["policy_grid.csv"]
    if not (policy["sampled_reads"] == evaluation_reads).all():
        raise AssertionError("Policy grid sample sizes differ")
    post = policy[_as_bool(policy["namespace_enforced"])].copy()
    for key, frame in post.groupby(
        ["sublibrary", "input_arm", "calibration_arm"]
    ):
        ordered = frame.sort_values("maximum_edit_distance")[
            "accepted_reads"
        ].tolist()
        if ordered != sorted(ordered):
            raise AssertionError(f"Non-monotonic policy acceptance for {key}")
    wide = policy.pivot_table(
        index=[
            "sublibrary",
            "input_arm",
            "calibration_arm",
            "maximum_edit_distance",
        ],
        columns="namespace_enforced",
        values="accepted_reads",
        aggfunc="first",
    )
    if True in wide.columns and False in wide.columns:
        if not (wide[True] <= wide[False]).all():
            raise AssertionError("Namespace enforcement increased acceptance")

    rescue = frames["linker_rescue_summary.csv"]
    rescued = rescue[rescue["rescue_state"] == "rescued_to_intended"]
    if not rescued.empty:
        if not (rescued["linker_status"] == "unique_pair").all():
            raise AssertionError("Linker rescue lacks a primary unique linker pair")
        if rescued["fixed_bc1_reference_category"].str.startswith(
            "intended_a1_b12"
        ).any():
            raise AssertionError("Already-intended fixed calls counted as rescue")
        if not rescued["anchored_bc1_category"].str.startswith("unique_").all():
            raise AssertionError("Linker rescue lacks a unique intended anchored call")

    offsets = frames["linker_offset_summary.csv"]
    expected_l23_delta = (
        offsets["l23_end_offset"] - offsets["l23_start_offset"]
    )
    expected_l12_delta = (
        offsets["l12_end_offset"] - offsets["l12_start_offset"]
    )
    if not (offsets["l23_length_delta"] == expected_l23_delta).all():
        raise AssertionError("L23 length deltas do not reconcile with offsets")
    if not (offsets["l12_length_delta"] == expected_l12_delta).all():
        raise AssertionError("L12 length deltas do not reconcile with offsets")
    if (
        offsets[["l23_length_delta", "l12_length_delta"]].abs() > 2
    ).any(axis=None):
        raise AssertionError("Primary linker indel length exceeds edit-distance bound")

    manifest_path = diagnostic_dir / "manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("contains_read_identifiers_or_sequences") is not False:
        raise AssertionError("Manifest privacy assertion is missing")
    if len(manifest.get("libraries", [])) != len(libraries):
        raise AssertionError("Aggregate manifest does not cover every library")
    sample_digests = set()
    for item in manifest["libraries"]:
        if set(item.get("decoder_profiles", {})) != DECODERS:
            raise AssertionError("Run manifest decoder profiles are incomplete")
        input_r2 = item.get("input_r2", {})
        if set(input_r2) != {"GEX", "gRNA"}:
            raise AssertionError("Run manifest input R2 identities are incomplete")
        for identities in input_r2.values():
            if not identities:
                raise AssertionError("Run manifest has an empty R2 input list")
            for identity in identities:
                path = Path(identity["path"])
                if not path.exists():
                    raise FileNotFoundError(path)
                stat = path.stat()
                if (
                    int(identity["size_bytes"]) != stat.st_size
                    or int(identity["mtime_ns"]) != stat.st_mtime_ns
                ):
                    raise AssertionError(f"R2 input identity changed: {path}")
        for key in ("step1a_calibration", "step2_calibration"):
            sample = item[key]
            calibration_start, calibration_end = sample["calibration_window"]
            pool_start, pool_end = sample["evaluation_pool"]
            if max(calibration_start, pool_start) < min(calibration_end, pool_end):
                raise AssertionError("Calibration overlaps evaluation pool")
            if sample["evaluation_reads"] != evaluation_reads:
                raise AssertionError("Manifest evaluation sample size differs")
            sample_digests.add(sample["selected_index_digest_sha256"])
    if len(sample_digests) != 2 * len(libraries):
        raise AssertionError("Evaluation sample digests are not arm/library-specific")

    for name in CROSSOVER_PRODUCTS:
        if (asset_dir / name).exists():
            raise AssertionError(
                f"Detailed crossover table is in docs assets: {name}"
            )
    oversized = [
        str(path)
        for path in asset_dir.glob("*.csv")
        if path.stat().st_size > maximum_asset_table_bytes
    ]
    if oversized:
        raise AssertionError(f"Oversized docs asset tables: {oversized}")

    result = {
        "status": "passed",
        "libraries": libraries,
        "evaluation_reads_per_library_arm": evaluation_reads,
        "total_evaluation_reads": evaluation_reads * len(expected_library_arms),
        "native_replay_warnings_over_3pp": int(warnings.sum()),
        "product_rows": {
            name: int(len(frame)) for name, frame in frames.items()
        },
        "published_identifiers_or_sequences": False,
        "detailed_tables_in_docs_assets": False,
    }
    (diagnostic_dir / "validation.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument(
        "--diagnostic-dir",
        default="output/pilot_4_sublib_audit/barcode_crossover",
    )
    parser.add_argument(
        "--asset-dir", default="docs/assets/pilot_4_sublib_audit"
    )
    parser.add_argument("--evaluation-reads", type=int, default=EVALUATION_READS)
    parser.add_argument("--maximum-asset-table-bytes", type=int, default=1_000_000)
    args = parser.parse_args()
    print(
        validate(
            Path(args.config).resolve(),
            Path(args.diagnostic_dir).resolve(),
            Path(args.asset_dir).resolve(),
            args.evaluation_reads,
            args.maximum_asset_table_bytes,
        )
    )


if __name__ == "__main__":
    main()

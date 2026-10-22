#!/usr/bin/env python3
"""Archived validator for the completed barcode-rejection audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import yaml


PRODUCTS = (
    "barcode_contract_parity.csv",
    "barcode_rejection_summary.csv",
    "barcode_failure_by_round.csv",
    "barcode_quality_by_status.csv",
    "barcode_policy_grid.csv",
    "barcode_calibration_transfer.csv",
    "barcode_well_balance.csv",
    "barcode_error_spectrum.csv",
)
FORBIDDEN_COLUMNS = {
    "read_id",
    "sequence",
    "r2_seq",
    "barcode_sequence",
    "bc_wells",
    "cell_barcode",
}


def _as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().map({"true": True, "false": False})


def validate(
    config_path: Path,
    diagnostic_dir: Path,
    asset_dir: Path,
    maximum_asset_table_bytes: int = 1_000_000,
) -> dict:
    cfg = yaml.safe_load(config_path.read_text())
    libraries = [item["sublibrary_name"] for item in cfg["libraries"]]
    frames = {}
    for name in PRODUCTS:
        path = diagnostic_dir / name
        if not path.exists():
            raise FileNotFoundError(path)
        frame = pd.read_csv(path)
        forbidden = FORBIDDEN_COLUMNS & set(frame.columns)
        if forbidden:
            raise AssertionError(f"{name} exposes forbidden columns: {forbidden}")
        frames[name] = frame

    summary = frames["barcode_rejection_summary.csv"]
    observed_libraries = set(summary["sublibrary"])
    if observed_libraries != set(libraries):
        raise AssertionError(
            f"Summary libraries {sorted(observed_libraries)} != {sorted(libraries)}"
        )
    partition = summary.groupby(
        ["sublibrary", "arm", "window"], as_index=False
    ).agg(
        category_reads=("read_count", "sum"),
        stated_reads=("window_read_count", "first"),
        stated_values=("window_read_count", "nunique"),
    )
    if not (partition["stated_values"] == 1).all():
        raise AssertionError("Window denominators are inconsistent across categories")
    if not (partition["category_reads"] == partition["stated_reads"]).all():
        raise AssertionError("Rejection categories do not partition sampled reads")

    policy = frames["barcode_policy_grid.csv"]
    post_namespace = policy[_as_bool(policy["namespace_enforced"])].copy()
    monotonic_groups = [
        "sublibrary",
        "input_arm",
        "whitelist_source",
        "calibration_arm",
    ]
    for key, frame in post_namespace.groupby(monotonic_groups):
        ordered = frame.sort_values("maximum_edit_distance")["accepted_reads"].tolist()
        if ordered != sorted(ordered):
            raise AssertionError(f"Non-monotonic edit-distance acceptance for {key}")
    policy_wide = policy.pivot_table(
        index=monotonic_groups + ["maximum_edit_distance"],
        columns="namespace_enforced",
        values="accepted_reads",
        aggfunc="first",
    )
    if True in policy_wide.columns and False in policy_wide.columns:
        if not (policy_wide[True] <= policy_wide[False]).all():
            raise AssertionError("Namespace enforcement increased accepted reads")

    parity = frames["barcode_contract_parity.csv"]
    expected_contract_rows = {
        (library, round_number)
        for library in libraries
        for round_number in (1, 2, 3)
    }
    observed_contract_rows = set(zip(parity["sublibrary"], parity["bc_round"]))
    if observed_contract_rows != expected_contract_rows:
        raise AssertionError("Contract parity rows do not cover every library/round")
    metadata_columns = [
        "canonical_subset_match",
        "kit_match",
        "chemistry_match",
        "amplicon_coordinate_match",
    ]
    if not all(_as_bool(parity[column]).all() for column in metadata_columns):
        raise AssertionError("Canonical Step-1a subset or runtime metadata mismatch")

    for name in PRODUCTS:
        if (asset_dir / name).exists():
            raise AssertionError(
                f"Detailed product must remain under output, not docs assets: {name}"
            )
    oversized_asset_tables = [
        str(path)
        for path in asset_dir.glob("*.csv")
        if path.stat().st_size > maximum_asset_table_bytes
    ]
    if oversized_asset_tables:
        raise AssertionError(
            f"Oversized docs asset tables: {oversized_asset_tables}"
        )

    result = {
        "status": "passed",
        "libraries": libraries,
        "sample_partitions": len(partition),
        "sampled_reads_by_library_arm": {
            f"{row.sublibrary}:{row.arm}": int(
                partition.loc[
                    (partition["sublibrary"] == row.sublibrary)
                    & (partition["arm"] == row.arm),
                    "category_reads",
                ].sum()
            )
            for row in partition[["sublibrary", "arm"]]
            .drop_duplicates()
            .itertuples(index=False)
        },
        "product_rows": {
            name: int(len(frame)) for name, frame in frames.items()
        },
        "published_identifiers_or_sequences": False,
        "detailed_tables_in_docs_assets": False,
        "maximum_asset_table_bytes": maximum_asset_table_bytes,
    }
    destination = diagnostic_dir / "validation.json"
    destination.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument(
        "--diagnostic-dir",
        default="output/pilot_4_sublib_audit/barcode_rejection",
    )
    parser.add_argument(
        "--asset-dir", default="docs/assets/pilot_4_sublib_audit"
    )
    parser.add_argument("--maximum-asset-table-bytes", type=int, default=1_000_000)
    args = parser.parse_args()
    print(
        validate(
            Path(args.config).resolve(),
            Path(args.diagnostic_dir).resolve(),
            Path(args.asset_dir).resolve(),
            args.maximum_asset_table_bytes,
        )
    )


if __name__ == "__main__":
    main()

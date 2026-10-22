#!/usr/bin/env python3
"""Archived: select anonymized poly(T)-barcode/rRNA smoke-run examples."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import yaml

try:
    from scripts.qc_fastq_resolver import resolve_library_pairs
    from scripts.rrna_barcode_diagnostic import (
        FAMILY_NAME,
        PROBE_CODE,
        alignment_candidate,
        iter_fastq_pairs,
        synthetic_serial,
    )
except ImportError:
    from qc_fastq_resolver import resolve_library_pairs
    from rrna_barcode_diagnostic import (
        FAMILY_NAME,
        PROBE_CODE,
        alignment_candidate,
        iter_fastq_pairs,
        synthetic_serial,
    )


def annotate_r2(sequence: str) -> dict[str, str]:
    return {
        "UMI": sequence[0:10],
        "BC3": sequence[10:18],
        "L23": sequence[18:30],
        "BC2": sequence[30:38],
        "L12": sequence[38:50],
        "BC1": sequence[50:58],
        "polyT15": sequence[58:73],
        "VN": sequence[73:75],
        "cDNA_after_primer": sequence[75:],
    }


def annotated_string(parts: dict[str, str]) -> str:
    return "|".join(f"{name}:{value}" for name, value in parts.items())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--library", default="C1")
    parser.add_argument("--smoke-root", required=True)
    parser.add_argument("--reference-dir", required=True)
    parser.add_argument("--examples", type=int, default=4)
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    cfg = yaml.safe_load(config_path.read_text())
    output_dir = Path(args.smoke_root).resolve() / "output" / args.library
    manifest = json.loads((output_dir / "run_manifest.json").read_text())
    scratch = Path(manifest["scratch"])
    total = int(manifest["summary"]["raw"])
    annotation = np.memmap(
        scratch / "barcode_annotations.uint8.mmap",
        mode="r",
        dtype=np.uint8,
        shape=(total, 12),
    )
    r1_calls = np.memmap(
        scratch / "r1_rrna.uint8.mmap",
        mode="r",
        dtype=np.uint8,
        shape=(total, 6),
    )
    eligible = (
        (annotation[:, 0] == 1)
        & (annotation[:, 5] == PROBE_CODE["polyT"])
        & (r1_calls[:, 0] == 1)
    )
    candidates = {}
    for family_code in np.unique(r1_calls[eligible, 2]):
        family = FAMILY_NAME[int(family_code)]
        serials = np.flatnonzero(eligible & (r1_calls[:, 2] == family_code))
        candidates[family] = set(int(value) for value in serials[:200])
    target_serials = set().union(*candidates.values())

    entry = next(
        item
        for item in cfg["libraries"]
        if item["sublibrary_name"] == args.library
    )
    pairs = resolve_library_pairs(cfg, str(config_path), entry)
    records = {}
    for serial, (r1, r2) in enumerate(iter_fastq_pairs(pairs)):
        if serial >= total:
            break
        if serial in target_serials:
            parts = annotate_r2(r2.sequence)
            records[serial] = {
                "serial": serial,
                "family": FAMILY_NAME[int(r1_calls[serial, 2])],
                "r1_full_sequence": r1.sequence,
                "r2_full_sequence": r2.sequence,
                "r2_segments": parts,
                "r2_annotated": annotated_string(parts),
                "polyT15_T_count": parts["polyT15"].count("T"),
            }

    reference_manifest = json.loads(
        (Path(args.reference_dir) / "reference_manifest.json").read_text()
    )
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("pysam is required") from exc
    with pysam.AlignmentFile(
        str(scratch / "r1_Aligned.out.bam"), "rb", check_sq=False
    ) as handle:
        for alignment in handle.fetch(until_eof=True):
            serial = synthetic_serial(alignment.query_name)
            if serial not in records or "alignment" in records[serial]:
                continue
            call = alignment_candidate(alignment, reference_manifest)
            records[serial]["alignment"] = {
                key: call[key]
                for key in (
                    "family",
                    "reference_name",
                    "reference_start",
                    "strand",
                    "query_start",
                    "query_end",
                    "aligned",
                    "nm",
                    "identity",
                )
            }

    family_priority = ("18S", "28S", "MT-RNR2", "MT-RNR1")
    selected = []
    for family in family_priority:
        options = [
            value
            for value in records.values()
            if value["family"] == family and "alignment" in value
        ]
        if options:
            selected.append(
                max(
                    options,
                    key=lambda value: (
                        value["polyT15_T_count"],
                        value["alignment"]["aligned"],
                    ),
                )
            )
        if len(selected) >= args.examples:
            break
    if len(selected) < args.examples:
        remaining = [
            value
            for value in records.values()
            if value not in selected and "alignment" in value
        ]
        selected.extend(
            sorted(
                remaining,
                key=lambda value: (
                    value["polyT15_T_count"],
                    value["alignment"]["aligned"],
                ),
                reverse=True,
            )[: args.examples - len(selected)]
        )
    payload = {
        "selection": (
            "production-valid polyT BC1, strict R1 rRNA; no raw read identifiers"
        ),
        "population_counts": {
            family: int(
                (
                    eligible
                    & (r1_calls[:, 2] == family_code)
                ).sum()
            )
            for family_code, family in FAMILY_NAME.items()
            if family_code and (r1_calls[eligible, 2] == family_code).any()
        },
        "examples": selected,
    }
    result_path = Path(args.smoke_root) / "private_polyT_rrna_examples.json"
    result_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Archived: extract identifier-free GEX R2 examples with invalid BC1."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import heapq
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Mapping

import pandas as pd
import yaml

try:
    from scripts.barcode_contract import load_step1a_contract
    from scripts.barcode_crossover_audit import (
        AMPLICON_LENGTH,
        FixedLengthBarcodeEditIndex,
        classify_round,
        fixed_parts,
        linker_anchored_parts,
        linker_offset_pattern,
        locate_linker_pair,
    )
    from scripts.barcode_rejection_qc import resolve_r2_paths, setup_runtime
except ImportError:
    from barcode_contract import load_step1a_contract
    from barcode_crossover_audit import (
        AMPLICON_LENGTH,
        FixedLengthBarcodeEditIndex,
        classify_round,
        fixed_parts,
        linker_anchored_parts,
        linker_offset_pattern,
        locate_linker_pair,
    )
    from barcode_rejection_qc import resolve_r2_paths, setup_runtime


DEFAULT_INTERVAL = (4_000_000, 6_000_000)
DEFAULT_CONTEXT_END = 100
DEFAULT_POOL_PER_CATEGORY = 12
DEFAULT_TOP_BC1_PER_LIBRARY = 1_000
MIN_ARTIFACT_EXACT_MATCH = 12

# Split TSO is the project construct sequence. Illumina adapter-trimming
# sequences follow Illumina's published TruSeq/UD-index documentation.
KNOWN_SEQUENCE_ARTIFACTS = {
    "Split_TSO": "AAGCAGTGGTATCAACGCAGAGTGAATGGG",
    "Illumina_Read1_adapter": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    "Illumina_Read2_adapter": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
}


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def normalize_fastq_id(header: str) -> str:
    token = header.split(maxsplit=1)[0].lstrip("@")
    return token.removesuffix("/1").removesuffix("/2")


def paired_gex_paths(cfg: Mapping, library: str) -> list[tuple[Path, Path]]:
    pairs = []
    for r2_path in resolve_r2_paths(cfg, library, "GEX"):
        r1_name = r2_path.name.replace("_R2_", "_R1_", 1)
        if r1_name == r2_path.name:
            raise ValueError(f"Cannot derive R1 path from {r2_path}")
        r1_path = r2_path.with_name(r1_name)
        if not r1_path.exists():
            raise FileNotFoundError(f"Paired R1 is missing: {r1_path}")
        pairs.append((r1_path, r2_path))
    return pairs


def iter_paired_records(
    pairs: list[tuple[Path, Path]],
    stop: int,
):
    """Yield synchronized R1/R2 sequence and quality without retaining read IDs."""
    logical_index = 0
    for r1_path, r2_path in pairs:
        if logical_index >= stop:
            break
        opener1 = gzip.open if str(r1_path).endswith(".gz") else open
        opener2 = gzip.open if str(r2_path).endswith(".gz") else open
        with opener1(r1_path, "rt") as r1_handle, opener2(
            r2_path, "rt"
        ) as r2_handle:
            while logical_index < stop:
                r1_header = r1_handle.readline()
                r2_header = r2_handle.readline()
                if not r1_header and not r2_header:
                    break
                if not r1_header or not r2_header:
                    raise ValueError(
                        f"Paired FASTQ length mismatch: {r1_path} versus {r2_path}"
                    )
                r1_sequence = r1_handle.readline().rstrip("\r\n").upper()
                r2_sequence = r2_handle.readline().rstrip("\r\n").upper()
                r1_plus = r1_handle.readline()
                r2_plus = r2_handle.readline()
                r1_quality = r1_handle.readline().rstrip("\r\n")
                r2_quality = r2_handle.readline().rstrip("\r\n")
                if (
                    not r1_sequence
                    or not r2_sequence
                    or not r1_plus
                    or not r2_plus
                    or not r1_quality
                    or not r2_quality
                ):
                    raise ValueError(
                        f"Truncated paired FASTQ record: {r1_path} / {r2_path}"
                    )
                if normalize_fastq_id(r1_header) != normalize_fastq_id(r2_header):
                    raise ValueError(
                        f"Paired FASTQ identifier mismatch at record {logical_index}"
                    )
                if len(r1_sequence) != len(r1_quality) or len(r2_sequence) != len(
                    r2_quality
                ):
                    raise ValueError(
                        f"FASTQ sequence/quality mismatch at record {logical_index}"
                    )
                yield (
                    logical_index,
                    r1_sequence,
                    r1_quality,
                    r2_sequence,
                    r2_quality,
                )
                logical_index += 1
    if logical_index < stop:
        raise ValueError(f"Only recovered {logical_index:,} pairs; need {stop:,}")


def longest_exact_overlap(query: str, reference: str) -> dict:
    """Return the longest exact common substring and zero-based coordinates."""
    best_length = 0
    best_query_start = -1
    best_reference_start = -1
    previous = [0] * (len(reference) + 1)
    for query_index, query_base in enumerate(query):
        current = [0] * (len(reference) + 1)
        for reference_index, reference_base in enumerate(reference):
            if query_base != reference_base:
                continue
            current[reference_index + 1] = previous[reference_index] + 1
            length = current[reference_index + 1]
            if length > best_length:
                best_length = length
                best_query_start = query_index - length + 1
                best_reference_start = reference_index - length + 1
        previous = current
    return {
        "length": best_length,
        "sequence": (
            query[best_query_start : best_query_start + best_length]
            if best_length
            else ""
        ),
        "query_start": best_query_start,
        "reference_start": best_reference_start,
    }


def sequence_after_l12_prefix(sequence: str) -> tuple[int | None, str]:
    """Return sequence after the BC2-coupled TCCAAC half of L12."""
    expected_start = 38
    if sequence[expected_start : expected_start + 6] == "TCCAAC":
        return expected_start, sequence[expected_start + 6 :]
    start = sequence.find("TCCAAC", expected_start - 4, expected_start + 10)
    if start < 0:
        return None, ""
    return start, sequence[start + 6 :]


def screen_sequence_artifacts(sequence: str) -> dict:
    """Screen post-TCCAAC sequence for exact TSO/Illumina adapter fragments."""
    prefix_start, query = sequence_after_l12_prefix(sequence)
    comparisons = []
    for artifact, reference in KNOWN_SEQUENCE_ARTIFACTS.items():
        for orientation, oriented_reference in (
            ("forward", reference),
            ("reverse_complement", reverse_complement(reference)),
        ):
            overlap = longest_exact_overlap(query, oriented_reference)
            comparisons.append(
                {
                    "artifact": artifact,
                    "orientation": orientation,
                    **overlap,
                }
            )
    comparisons.sort(
        key=lambda row: (
            -row["length"],
            row["artifact"],
            row["orientation"],
        )
    )
    return {
        "l12_prefix_start": prefix_start,
        "sequence_after_tccaac": query,
        "best_artifact_overlap": comparisons[0] if comparisons else None,
        "supported_artifact_hits": [
            row
            for row in comparisons
            if row["length"] >= MIN_ARTIFACT_EXACT_MATCH
        ],
        "minimum_supported_exact_match": MIN_ARTIFACT_EXACT_MATCH,
    }


def exact_adapter_start(
    sequence: str,
    adapter: str,
    minimum_anchor: int = MIN_ARTIFACT_EXACT_MATCH,
) -> int | None:
    """Find the earliest read-through adapter prefix supported by an exact anchor."""
    for start in range(len(sequence)):
        overlap = min(len(sequence) - start, len(adapter))
        if overlap < minimum_anchor:
            break
        if sequence[start : start + overlap] == adapter[:overlap]:
            return start
    return None


def reconstruct_paired_fragment(
    r1_sequence: str,
    r2_sequence: str,
    minimum_overlap: int = 20,
    maximum_mismatch_fraction: float = 0.10,
) -> dict:
    """Trim exact adapter read-through and reconstruct an ungapped PE insert."""
    r1_adapter_start = exact_adapter_start(
        r1_sequence, KNOWN_SEQUENCE_ARTIFACTS["Illumina_Read1_adapter"]
    )
    r2_adapter_start = exact_adapter_start(
        r2_sequence, KNOWN_SEQUENCE_ARTIFACTS["Illumina_Read2_adapter"]
    )
    r1_insert = (
        r1_sequence[:r1_adapter_start]
        if r1_adapter_start is not None
        else r1_sequence
    )
    r2_insert = (
        r2_sequence[:r2_adapter_start]
        if r2_adapter_start is not None
        else r2_sequence
    )
    r2_insert_rc = reverse_complement(r2_insert)
    candidates = []
    maximum_overlap = min(len(r1_insert), len(r2_insert_rc))
    for overlap in range(minimum_overlap, maximum_overlap + 1):
        left = r1_insert[-overlap:]
        right = r2_insert_rc[:overlap]
        mismatches = sum(a != b for a, b in zip(left, right))
        mismatch_fraction = mismatches / overlap
        if mismatch_fraction <= maximum_mismatch_fraction:
            candidates.append((overlap, mismatches, mismatch_fraction))
    if not candidates:
        return {
            "status": "no_supported_overlap",
            "r1_adapter_start": r1_adapter_start,
            "r2_adapter_start": r2_adapter_start,
            "r1_insert_sequence": r1_insert,
            "r2_insert_sequence": r2_insert,
            "overlap_length": None,
            "overlap_mismatches": None,
            "overlap_mismatch_fraction": None,
            "estimated_fragment_length": None,
            "reconstructed_fragment_r1_orientation": "",
            "reconstructed_fragment_r2_orientation": "",
        }
    overlap, mismatches, mismatch_fraction = max(
        candidates,
        key=lambda value: (value[0], -value[1]),
    )
    reconstructed_r1 = r1_insert + r2_insert_rc[overlap:]
    return {
        "status": "supported_overlap",
        "r1_adapter_start": r1_adapter_start,
        "r2_adapter_start": r2_adapter_start,
        "r1_insert_sequence": r1_insert,
        "r2_insert_sequence": r2_insert,
        "overlap_length": overlap,
        "overlap_mismatches": mismatches,
        "overlap_mismatch_fraction": mismatch_fraction,
        "estimated_fragment_length": len(reconstructed_r1),
        "reconstructed_fragment_r1_orientation": reconstructed_r1,
        "reconstructed_fragment_r2_orientation": reverse_complement(
            reconstructed_r1
        ),
    }


def mean_phred(quality: str) -> float | None:
    if not quality:
        return None
    return sum(ord(value) - 33 for value in quality) / len(quality)


def longest_run(sequence: str, base: str) -> int:
    longest = current = 0
    for value in sequence:
        current = current + 1 if value == base else 0
        longest = max(longest, current)
    return longest


def stable_rank(library: str, sequence: str, quality: str) -> int:
    payload = f"bc1_examples_v1\x1f{library}\x1f{sequence}\x1f{quality}".encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")


def classify_example(
    bc1_category: str,
    linker_status: str,
    offset_pattern: str,
    anchored_bc1_category: str,
) -> str:
    if linker_status == "unique_pair":
        if offset_pattern == "intact_position":
            suffix = (
                "no_neighbor"
                if bc1_category == "no_neighbor_d2"
                else "ambiguous"
            )
            return f"intact_linkers_bc1_{suffix}"
        if offset_pattern == "intact_position_linker_damage":
            return "intact_position_linker_damage"
        if anchored_bc1_category.startswith("unique_"):
            return "linker_shift_rescues_bc1"
        return "linker_shift_without_bc1_rescue"
    return "unresolved_or_partial_linkers"


def add_bounded_candidate(
    heaps: dict[str, list],
    category: str,
    rank: int,
    record: dict,
    maximum_per_category: int,
) -> None:
    item = (-rank, record["sequence_first_100"], record)
    heap = heaps[category]
    if len(heap) < maximum_per_category:
        heapq.heappush(heap, item)
    elif item > heap[0]:
        heapq.heapreplace(heap, item)


def scan_library(
    config_path: Path,
    library: str,
    destination: Path,
    start: int = DEFAULT_INTERVAL[0],
    end: int = DEFAULT_INTERVAL[1],
    context_end: int = DEFAULT_CONTEXT_END,
    pool_per_category: int = DEFAULT_POOL_PER_CATEGORY,
    top_bc1_per_library: int = DEFAULT_TOP_BC1_PER_LIBRARY,
) -> dict:
    cfg = yaml.safe_load(config_path.read_text())
    bcinfo, _kits, _kit, _chemistry, amplicon = setup_runtime(cfg)
    contract = load_step1a_contract(cfg, config_path, library, amplicon)
    indices = {
        round_number: FixedLengthBarcodeEditIndex(
            contract.barcode_table.loc[
                contract.barcode_table["bc_round"] == round_number,
                "sequence",
            ],
            max_dist=2,
        )
        for round_number in (1, 2, 3)
    }
    intended_bc1 = set(
        contract.barcode_table.loc[
            (contract.barcode_table["bc_round"] == 1)
            & contract.barcode_table["well"].isin(contract.sample_wells),
            "sequence",
        ]
    )
    intended_bc1_index = FixedLengthBarcodeEditIndex(
        intended_bc1, max_dist=2
    )
    counts = Counter()
    bc1_sequence_counts = Counter()
    bc1_outcome_counts: dict[str, Counter] = defaultdict(Counter)
    bc1_structure_counts: dict[str, Counter] = defaultdict(Counter)
    representative_by_bc1: dict[str, tuple[int, dict]] = {}
    heaps: dict[str, list] = defaultdict(list)
    examined = 0
    for (
        record_index,
        r1_sequence,
        _r1_quality,
        sequence,
        quality,
    ) in iter_paired_records(paired_gex_paths(cfg, library), end):
        if record_index < start:
            continue
        examined += 1
        parts = fixed_parts(sequence, contract)
        bc2_category, _bc2_distance, bc2_matches = classify_round(
            parts[2], indices[2]
        )
        bc3_category, _bc3_distance, bc3_matches = classify_round(
            parts[3], indices[3]
        )
        if bc2_category != "unique_d0" or bc3_category != "unique_d0":
            continue
        bc1_category, _bc1_distance, _bc1_matches = classify_round(
            parts[1], indices[1]
        )
        if bc1_category.startswith("unique_"):
            continue
        pair = locate_linker_pair(sequence)
        pattern = linker_offset_pattern(pair)
        anchored = linker_anchored_parts(sequence, pair)
        anchored_category, _anchored_distance, _anchored_matches = (
            classify_round(anchored[1], intended_bc1_index)
            if anchored[1]
            else ("not_evaluable", None, set())
        )
        category = classify_example(
            bc1_category,
            pair.status,
            pattern,
            anchored_category,
        )
        counts[category] += 1
        bc1_sequence_counts[parts[1]] += 1
        bc1_outcome_counts[parts[1]][bc1_category] += 1
        bc1_structure_counts[parts[1]][category] += 1
        downstream = sequence[AMPLICON_LENGTH:]
        record = {
            "sublibrary": library,
            "example_category": category,
            "bc1_category": bc1_category,
            "linker_status": pair.status,
            "offset_pattern": pattern,
            "anchored_bc1_category": anchored_category,
            "full_r1_sequence": r1_sequence,
            "full_r2_sequence": sequence,
            "sequence_first_100": sequence[:context_end],
            "read_length": len(sequence),
            "umi": sequence[0:10],
            "bc3": parts[3],
            "l23_fixed": sequence[18:30],
            "bc2": parts[2],
            "l12_fixed": sequence[38:50],
            "bc1_fixed": parts[1],
            "post_bc1_context": downstream,
            "post_bc1_t_fraction": (
                downstream.count("T") / len(downstream)
                if downstream
                else None
            ),
            "post_bc1_longest_t_run": longest_run(downstream, "T"),
            "umi_mean_phred": mean_phred(quality[0:10]),
            "bc3_mean_phred": mean_phred(quality[10:18]),
            "l23_mean_phred": mean_phred(quality[18:30]),
            "bc2_mean_phred": mean_phred(quality[30:38]),
            "l12_mean_phred": mean_phred(quality[38:50]),
            "bc1_mean_phred": mean_phred(quality[50:58]),
            "post_bc1_mean_phred": mean_phred(
                quality[AMPLICON_LENGTH:context_end]
            ),
            "bc2_reference": next(iter(bc2_matches)),
            "bc3_reference": next(iter(bc3_matches)),
            "anchored_bc1": anchored[1],
            "l23_start": pair.l23.start if pair.l23 else None,
            "l23_end": pair.l23.end if pair.l23 else None,
            "l23_edit_distance": pair.l23.distance if pair.l23 else None,
            "l12_start": pair.l12.start if pair.l12 else None,
            "l12_end": pair.l12.end if pair.l12 else None,
            "l12_edit_distance": pair.l12.distance if pair.l12 else None,
            "selection_rank": stable_rank(library, sequence, quality),
        }
        add_bounded_candidate(
            heaps,
            category,
            record["selection_rank"],
            record,
            pool_per_category,
        )
        prior = representative_by_bc1.get(parts[1])
        if prior is None or record["selection_rank"] < prior[0]:
            representative_by_bc1[parts[1]] = (
                record["selection_rank"],
                record,
            )
    if examined != end - start:
        raise AssertionError(
            f"{library}: examined {examined:,}, expected {end - start:,}"
        )
    candidates = [
        item[2]
        for heap in heaps.values()
        for item in sorted(heap, reverse=True)
    ]
    top_bc1_examples = []
    for sequence_value, count in sorted(
        bc1_sequence_counts.items(), key=lambda item: (-item[1], item[0])
    )[:top_bc1_per_library]:
        record = dict(representative_by_bc1[sequence_value][1])
        record["library_bc1_occurrences"] = count
        record["library_bc1_outcome_counts"] = dict(
            sorted(bc1_outcome_counts[sequence_value].items())
        )
        record["library_bc1_structure_counts"] = dict(
            sorted(bc1_structure_counts[sequence_value].items())
        )
        top_bc1_examples.append(record)
    payload = {
        "schema_version": 1,
        "sublibrary": library,
        "arm": "GEX",
        "interval_zero_based_half_open": [start, end],
        "examined_reads": examined,
        "criteria": {
            "bc2": "exact canonical whitelist match",
            "bc3": "exact canonical whitelist match",
            "bc1": "not uniquely mappable through edit distance 2",
        },
        "category_counts": dict(sorted(counts.items())),
        "qualifying_reads": int(sum(bc1_sequence_counts.values())),
        "bc1_sequence_counts": dict(sorted(bc1_sequence_counts.items())),
        "bc1_outcome_counts": {
            sequence_value: dict(sorted(values.items()))
            for sequence_value, values in sorted(bc1_outcome_counts.items())
        },
        "bc1_structure_counts": {
            sequence_value: dict(sorted(values.items()))
            for sequence_value, values in sorted(
                bc1_structure_counts.items()
            )
        },
        "top_bc1_examples": top_bc1_examples,
        "candidate_pool": candidates,
        "contains_read_ids": False,
        "contains_r2_sequences": True,
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def choose_examples(candidate_pool: list[dict], per_category: int) -> list[dict]:
    selected = []
    frame = pd.DataFrame(candidate_pool)
    if frame.empty:
        return selected
    for category, subset in frame.groupby("example_category", sort=True):
        subset = subset.sort_values("selection_rank")
        used_libraries = set()
        rows = []
        for row in subset.to_dict("records"):
            if row["sublibrary"] not in used_libraries:
                rows.append(row)
                used_libraries.add(row["sublibrary"])
            if len(rows) == per_category:
                break
        if len(rows) < per_category:
            selected_ranks = {row["selection_rank"] for row in rows}
            for row in subset.to_dict("records"):
                if row["selection_rank"] in selected_ranks:
                    continue
                rows.append(row)
                if len(rows) == per_category:
                    break
        selected.extend(rows)
    return sorted(
        selected,
        key=lambda row: (row["example_category"], row["selection_rank"]),
    )


def annotated_line(row: Mapping) -> str:
    return (
        f"🟪 **UMI** `{row['umi']}` · "
        f"🟦 **BC3 ✓** `{row['bc3']}` · "
        f"🟨 **L23** `{row['l23_fixed']}` · "
        f"🟦 **BC2 ✓** `{row['bc2']}` · "
        f"🟨 **L12** `{row['l12_fixed']}` · "
        f"🟥 **BC1 ✗** `{row['bc1_fixed']}` · "
        f"⬜ **after BC1** `{row['post_bc1_context']}`"
    )


def artifact_hit_text(row: Mapping) -> str:
    hits = row.get("supported_artifact_hits", [])
    if not hits:
        best = row.get("best_artifact_overlap") or {}
        return (
            "No supported exact TSO/Illumina match "
            f"(threshold {row.get('minimum_supported_exact_match', MIN_ARTIFACT_EXACT_MATCH)} nt; "
            f"longest observed overlap {best.get('length', 0)} nt)."
        )
    return "; ".join(
        (
            f"{hit['artifact']} ({hit['orientation']}), "
            f"{hit['length']}-nt exact match `{hit['sequence']}` at "
            f"post-`TCCAAC` offset {hit['query_start']}"
        )
        for hit in hits
    )


def paired_fragment_text(row: Mapping) -> str:
    paired = row["paired_fragment"]
    adapter_text = (
        f"R1 adapter start {paired['r1_adapter_start']}; "
        f"R2 adapter start {paired['r2_adapter_start']}"
    )
    if paired["status"] != "supported_overlap":
        return f"No supported ≥20-nt mate overlap; {adapter_text}."
    return (
        f"{paired['overlap_length']}-nt overlap with "
        f"{paired['overlap_mismatches']} mismatch(es) "
        f"({100 * paired['overlap_mismatch_fraction']:.2f}%); "
        f"estimated insert length {paired['estimated_fragment_length']} nt; "
        f"{adapter_text}."
    )


def compact_markdown_table(frame: pd.DataFrame) -> str:
    columns = list(frame.columns)
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in frame.itertuples(index=False, name=None):
        values = [
            (
                f"{value:.6f}"
                if isinstance(value, float)
                else str(value)
            )
            for value in row
        ]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def write_markdown(
    destination: Path,
    selected: list[dict],
    counts: pd.DataFrame,
    top_sequences: pd.DataFrame,
) -> None:
    pair_rows = []
    for index, row in enumerate(selected, start=1):
        paired = row["paired_fragment"]
        hits = row.get("supported_artifact_hits", [])
        pair_rows.append(
            {
                "Example": index,
                "Library": row["sublibrary"],
                "Fixed BC1": row["bc1_fixed"],
                "BC1-side L12 half": (
                    "CACCTC"
                    if row["full_r2_sequence"][44:50] == "CACCTC"
                    else row["full_r2_sequence"][44:50]
                ),
                "TSO/adapter evidence": (
                    ", ".join(
                        f"{hit['artifact']} {hit['orientation']} {hit['length']}nt"
                        for hit in hits
                    )
                    if hits
                    else "none ≥12nt"
                ),
                "Mate overlap": (
                    f"{paired['overlap_length']}nt/"
                    f"{paired['overlap_mismatches']}mm"
                    if paired["status"] == "supported_overlap"
                    else "not supported"
                ),
                "Insert length": (
                    paired["estimated_fragment_length"]
                    if paired["estimated_fragment_length"] is not None
                    else "not bounded"
                ),
            }
        )
    pair_summary = pd.DataFrame(pair_rows)
    lines = [
        "# Representative GEX R2 reads with exact BC2/BC3 and invalid BC1",
        "",
        "Expected fixed-cycle structure:",
        "",
        "```text",
        "cycles    1-10   11-18  19-30  31-38  39-50  51-58  59-150",
        "field     UMI    BC3    L23    BC2    L12    BC1    downstream context",
        "```",
        "",
        "L12 is physically composite: `TCCAAC` is coupled to BC2, while",
        "`CACCTC` is supplied by the Round-1 probe immediately before BC1.",
        "Consequently, reads with `TCCAAC` but without the following `CACCTC`",
        "have lost the BC1-side half-linker; their nominal fixed-cycle BC1 can",
        "instead be sequence from the underlying insert.",
        "",
        "BC2 and BC3 are exact canonical whitelist matches. BC1 is not uniquely",
        "mappable through edit distance 2. Read identifiers are intentionally omitted.",
        "Each example below includes the complete paired R1/R2 sequences. Exact Split-TSO and",
        "Illumina adapter fragments are screened in both orientations after the",
        f"`TCCAAC` half-linker; matches shorter than {MIN_ARTIFACT_EXACT_MATCH} nt",
        "are treated as chance-compatible rather than positive evidence.",
        "",
        "## Paired-read and artifact summary",
        "",
        compact_markdown_table(pair_summary),
        "",
        "Complete R1/R2, post-`TCCAAC`, and reconstructed-insert FASTA files",
        "are stored beside this report for manual BLASTN or another aligner.",
        "",
        "## Most frequent rejected fixed-slice BC1 sequences",
        "",
        compact_markdown_table(top_sequences),
        "",
    ]
    for index, row in enumerate(selected, start=1):
        lines.extend(
            [
                f"## Example {index}: {row['sublibrary']} / {row['example_category']}",
                "",
                f"**Complete {len(row['full_r1_sequence'])}-nt R1 sequence**",
                "",
                "```text",
                row["full_r1_sequence"],
                "```",
                "",
                f"**Complete {row['read_length']}-nt R2 sequence**",
                "",
                "```text",
                row["full_r2_sequence"],
                "```",
                "",
                annotated_line(row),
                "",
                (
                    f"- Fixed BC1 rank: {row.get('global_bc1_rank', 'NA')}; "
                    f"{row.get('global_bc1_occurrences', 'NA')} occurrences "
                    f"({100 * row.get('fraction_of_qualifying_reads', 0):.3f}% "
                    "of qualifying rejected reads); selection role: "
                    f"`{row.get('selection_role', 'unspecified')}`."
                ),
                (
                    f"- BC1 outcome: `{row['bc1_category']}`; linker outcome: "
                    f"`{row['linker_status']}` / `{row['offset_pattern']}`."
                ),
                (
                    f"- Segment mean Phred: UMI {row['umi_mean_phred']:.1f}, "
                    f"BC3 {row['bc3_mean_phred']:.1f}, "
                    f"L23 {row['l23_mean_phred']:.1f}, "
                    f"BC2 {row['bc2_mean_phred']:.1f}, "
                    f"L12 {row['l12_mean_phred']:.1f}, "
                    f"BC1 {row['bc1_mean_phred']:.1f}."
                ),
                (
                    f"- Post-BC1 context: {100 * row['post_bc1_t_fraction']:.1f}% T; "
                    f"longest T run {row['post_bc1_longest_t_run']} nt."
                ),
                f"- Artifact screen: {artifact_hit_text(row)}",
                f"- Paired-read reconstruction: {paired_fragment_text(row)}",
            ]
        )
        paired = row["paired_fragment"]
        if paired["status"] == "supported_overlap":
            lines.extend(
                [
                    "",
                    "**Reconstructed insert in R1 orientation**",
                    "",
                    "```text",
                    paired["reconstructed_fragment_r1_orientation"],
                    "```",
                    "",
                    "**Reconstructed insert in R2/barcode orientation**",
                    "",
                    "```text",
                    paired["reconstructed_fragment_r2_orientation"],
                    "```",
                    "",
                    (
                        "- Reconstructed sequence following `TCCAAC`: "
                        f"`{paired['reconstructed_sequence_after_tccaac']}`"
                    ),
                ]
            )
        if row["anchored_bc1_category"].startswith("unique_"):
            lines.append(
                f"- Linker-anchored BC1 rescue: `{row['anchored_bc1']}` "
                f"(`{row['anchored_bc1_category']}`)."
            )
        lines.append("")
    lines.extend(
        [
            "## Candidate counts in the audited interval",
            "",
            compact_markdown_table(counts),
            "",
        ]
    )
    destination.write_text("\n".join(lines))


def write_fasta_products(destination: Path, selected: list[dict]) -> None:
    r1_lines = []
    full_lines = []
    suffix_lines = []
    reconstructed_lines = []
    for index, row in enumerate(selected, start=1):
        label = (
            f"example_{index:02d}_{row['sublibrary']}_"
            f"bc1_{row['bc1_fixed']}"
        )
        r1_lines.extend([f">{label}_full_R1", row["full_r1_sequence"]])
        full_lines.extend([f">{label}_full_R2", row["full_r2_sequence"]])
        suffix_lines.extend(
            [
                f">{label}_after_TCCAAC",
                row.get("sequence_after_tccaac", ""),
            ]
        )
        paired = row["paired_fragment"]
        if paired["status"] == "supported_overlap":
            reconstructed_lines.extend(
                [
                    f">{label}_insert_R1_orientation",
                    paired["reconstructed_fragment_r1_orientation"],
                    f">{label}_insert_R2_barcode_orientation",
                    paired["reconstructed_fragment_r2_orientation"],
                ]
            )
    (destination / "representative_full_r1.fasta").write_text(
        "\n".join(r1_lines) + "\n"
    )
    (destination / "representative_full_r2.fasta").write_text(
        "\n".join(full_lines) + "\n"
    )
    (destination / "representative_after_tccaac.fasta").write_text(
        "\n".join(suffix_lines) + "\n"
    )
    (destination / "representative_reconstructed_inserts.fasta").write_text(
        ("\n".join(reconstructed_lines) + "\n") if reconstructed_lines else ""
    )


def aggregate(
    config_path: Path,
    work_dir: Path,
    output_dir: Path,
    top_sequences: int,
    top_examples: int,
) -> dict:
    cfg = yaml.safe_load(config_path.read_text())
    libraries = [item["sublibrary_name"] for item in cfg["libraries"]]
    payloads = [
        json.loads((work_dir / f"{library}.json").read_text())
        for library in libraries
    ]
    counts = pd.DataFrame(
        [
            {
                "sublibrary": payload["sublibrary"],
                "example_category": category,
                "read_count": count,
                "examined_reads": payload["examined_reads"],
                "fraction_of_examined_gex_r2": (
                    count / payload["examined_reads"]
                ),
            }
            for payload in payloads
            for category, count in payload["category_counts"].items()
        ]
    )
    candidates = [
        row for payload in payloads for row in payload["candidate_pool"]
    ]
    combined_bc1_counts = Counter()
    combined_outcome_counts: dict[str, Counter] = defaultdict(Counter)
    combined_structure_counts: dict[str, Counter] = defaultdict(Counter)
    for payload in payloads:
        combined_bc1_counts.update(payload["bc1_sequence_counts"])
        for sequence_value, values in payload["bc1_outcome_counts"].items():
            combined_outcome_counts[sequence_value].update(values)
        for sequence_value, values in payload["bc1_structure_counts"].items():
            combined_structure_counts[sequence_value].update(values)
    qualifying_reads = sum(combined_bc1_counts.values())
    ranked_sequences = sorted(
        combined_bc1_counts.items(), key=lambda item: (-item[1], item[0])
    )
    top_rows = []
    for rank, (sequence_value, count) in enumerate(
        ranked_sequences[:top_sequences], start=1
    ):
        outcome = combined_outcome_counts[sequence_value]
        structure = combined_structure_counts[sequence_value]
        top_rows.append(
            {
                "rank": rank,
                "bc1_fixed_sequence": sequence_value,
                "occurrences": count,
                "fraction_of_qualifying_reads": count / qualifying_reads,
                **{
                    f"{library}_occurrences": int(
                        next(
                            payload["bc1_sequence_counts"].get(
                                sequence_value, 0
                            )
                            for payload in payloads
                            if payload["sublibrary"] == library
                        )
                    )
                    for library in libraries
                },
                "dominant_bc1_outcome": (
                    outcome.most_common(1)[0][0] if outcome else ""
                ),
                "dominant_structure": (
                    structure.most_common(1)[0][0] if structure else ""
                ),
            }
        )
    top_frame = pd.DataFrame(top_rows)
    example_lookup: dict[str, list[dict]] = defaultdict(list)
    for payload in payloads:
        for row in payload["top_bc1_examples"]:
            example_lookup[row["bc1_fixed"]].append(row)
    selected = []
    for global_rank, (sequence_value, count) in enumerate(
        ranked_sequences[:top_examples], start=1
    ):
        available = example_lookup.get(sequence_value, [])
        if not available:
            continue
        record = dict(
            min(available, key=lambda row: row["selection_rank"])
        )
        record["global_bc1_rank"] = global_rank
        record["global_bc1_occurrences"] = count
        record["fraction_of_qualifying_reads"] = count / qualifying_reads
        record["selection_role"] = "top_frequency"
        selected.append(record)
    for global_rank, (sequence_value, count) in enumerate(
        ranked_sequences, start=1
    ):
        dominant_structure = combined_structure_counts[
            sequence_value
        ].most_common(1)[0][0]
        if not dominant_structure.startswith("intact_linkers_"):
            continue
        available = [
            row
            for row in example_lookup.get(sequence_value, [])
            if row["example_category"].startswith("intact_linkers_")
        ]
        if not available:
            continue
        if any(row["bc1_fixed"] == sequence_value for row in selected):
            break
        record = dict(
            min(available, key=lambda row: row["selection_rank"])
        )
        record["global_bc1_rank"] = global_rank
        record["global_bc1_occurrences"] = count
        record["fraction_of_qualifying_reads"] = count / qualifying_reads
        record["selection_role"] = "intact_linker_control"
        selected.append(record)
        break
    for record in selected:
        record.update(screen_sequence_artifacts(record["full_r2_sequence"]))
        paired = reconstruct_paired_fragment(
            record["full_r1_sequence"],
            record["full_r2_sequence"],
        )
        if paired["status"] == "supported_overlap":
            reconstructed_r2 = paired[
                "reconstructed_fragment_r2_orientation"
            ]
            _prefix_start, reconstructed_suffix = sequence_after_l12_prefix(
                reconstructed_r2
            )
            paired["reconstructed_sequence_after_tccaac"] = (
                reconstructed_suffix
            )
            paired["reconstructed_sequence_after_fixed_bc1"] = (
                reconstructed_r2[AMPLICON_LENGTH:]
                if len(reconstructed_r2) > AMPLICON_LENGTH
                else ""
            )
        else:
            paired["reconstructed_sequence_after_tccaac"] = ""
            paired["reconstructed_sequence_after_fixed_bc1"] = ""
        record["paired_fragment"] = paired
    output_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(output_dir / "category_counts.csv", index=False)
    top_frame.to_csv(output_dir / "top_bc1_sequences.csv", index=False)
    (output_dir / "representative_reads.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "selected_examples": selected,
                "contains_read_ids": False,
                "contains_r2_sequences": True,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    write_markdown(
        output_dir / "REPRESENTATIVE_BC1_FAILURES.md",
        selected,
        counts,
        top_frame,
    )
    write_fasta_products(output_dir, selected)
    summary = {
        "libraries": libraries,
        "examined_reads": sum(payload["examined_reads"] for payload in payloads),
        "selected_examples": len(selected),
        "qualifying_reads": qualifying_reads,
        "unique_rejected_bc1_sequences": len(combined_bc1_counts),
        "contains_read_ids": False,
        "contains_r2_sequences": True,
        "sequence_extent": "complete_paired_R1_R2",
        "artifact_screen_minimum_exact_match": MIN_ARTIFACT_EXACT_MATCH,
        "fasta_products": [
            "representative_full_r1.fasta",
            "representative_full_r2.fasta",
            "representative_after_tccaac.fasta",
            "representative_reconstructed_inserts.fasta",
        ],
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    scan = subparsers.add_parser("scan")
    scan.add_argument("--config", default="config.yaml")
    scan.add_argument("--library", required=True)
    scan.add_argument(
        "--destination",
        required=True,
        help="Private JSON output; contains selected R2 sequences but no read IDs",
    )
    scan.add_argument("--start", type=int, default=DEFAULT_INTERVAL[0])
    scan.add_argument("--end", type=int, default=DEFAULT_INTERVAL[1])
    scan.add_argument("--context-end", type=int, default=DEFAULT_CONTEXT_END)
    scan.add_argument(
        "--pool-per-category", type=int, default=DEFAULT_POOL_PER_CATEGORY
    )
    scan.add_argument(
        "--top-bc1-per-library",
        type=int,
        default=DEFAULT_TOP_BC1_PER_LIBRARY,
    )
    combine = subparsers.add_parser("aggregate")
    combine.add_argument("--config", default="config.yaml")
    combine.add_argument("--work-dir", required=True)
    combine.add_argument("--output-dir", required=True)
    combine.add_argument("--top-sequences", type=int, default=25)
    combine.add_argument("--top-examples", type=int, default=8)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    config_path = Path(args.config).resolve()
    if args.command == "scan":
        result = scan_library(
            config_path,
            args.library,
            Path(args.destination).resolve(),
            args.start,
            args.end,
            args.context_end,
            args.pool_per_category,
            args.top_bc1_per_library,
        )
        print(
            {
                "sublibrary": result["sublibrary"],
                "examined_reads": result["examined_reads"],
                "category_counts": result["category_counts"],
            }
        )
    else:
        result = aggregate(
            config_path,
            Path(args.work_dir).resolve(),
            Path(args.output_dir).resolve(),
            args.top_sequences,
            args.top_examples,
        )
        print(result)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""QC and filtering benchmark for gRNA carry-over in Parse GEX libraries.

The module uses full-length pilot R1 reads to define a conservative structural
surrogate truth, truncates those same reads in memory to a requested length,
and evaluates filters that remain possible at that shorter length.  It can
also inspect split-pipe's annotated BAM to measure how often contaminating
construct reads were uniquely aligned and assigned to genes.

Modes
-----
benchmark
    Scan raw paired FASTQs, classify full-length R1, benchmark short-read
    filters, and optionally write a paired synthetic short-read FASTQ subset.
bam
    Scan barcode_headAligned_anno.bam and quantify mapped/assigned construct
    reads.  Reverse-strand BAM sequences are restored to sequencing orientation.
all
    Run benchmark followed by BAM audit.
filter
    Apply one configured short-read filter to all input pairs and write paired
    clean/rejected FASTQs suitable for a later split-pipe comparison.
alignment
    Audit one controlled STAR arm (current full R1 or raw R1 truncated to the
    configured short length) using ground-truth metadata retained in QNAME.
guide
    Compare Step 1c GEX-derived guide molecules with the dedicated guide arm,
    including exact overlap and combined-matrix deltas.
aggregate
    Combine completed per-library outputs.

The full-length label is deliberately structural: fuzzy Split_TSO followed by
the expected spacer-length interval, a complete CSM stem, an iBAR-sized block,
and the correctly paired partial tracr.  Guide/iBAR reference matches provide
additional confidence tiers but are not required for the primary structure
label, so unknown or recombined contaminants are not silently treated as GEX.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import itertools
import json
import math
import os
import random
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

import pandas as pd
import yaml

try:
    from .qc_fastq_resolver import (
        discover_token_bounded_pairs,
        resolve_library_pairs,
    )
except ImportError:  # Direct execution from scripts/ on the HPC.
    from qc_fastq_resolver import (
        discover_token_bounded_pairs,
        resolve_library_pairs,
    )

try:
    import edlib

    HAS_EDLIB = True
except ImportError:
    edlib = None
    HAS_EDLIB = False


DNA = "ACGT"
REVCOMP_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(REVCOMP_TABLE)[::-1]


def hamming_distance(a: str, b: str, max_dist: Optional[int] = None) -> int:
    """Hamming distance, counting a length difference as mismatching bases."""
    distance = abs(len(a) - len(b))
    if max_dist is not None and distance > max_dist:
        return distance
    for x, y in zip(a, b):
        if x != y:
            distance += 1
            if max_dist is not None and distance > max_dist:
                break
    return distance


def normalize_read_id(header: str) -> str:
    token = header.strip().split()[0]
    if token.startswith("@"):
        token = token[1:]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def stable_fraction(read_id: str) -> float:
    digest = hashlib.blake2b(read_id.encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "big") / float(2**64)


def open_text(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def file_sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_stat_identity(path: str) -> dict:
    stat = os.stat(path)
    return {
        "path": str(path),
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }


def benchmark_config_sha256(cfg: dict, library_name: str) -> str:
    qc = cfg.get("gex_grna_contamination_qc", {})
    qc_relevant = {
        key: value
        for key, value in qc.items()
        if key not in {"output_dir", "plots"}
    }
    payload = {
        "library": get_library(cfg, library_name),
        "qc": qc_relevant,
        "grna_structure": cfg.get("grna_structure", {}),
        "paths": {
            key: cfg.get("paths", {}).get(key, "")
            for key in ("fastq_dir", "output_dir", "grna_reference", "genome_dir")
        },
    }
    serialized = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(serialized).hexdigest()


@dataclass
class FastqRecord:
    header: str
    sequence: str
    plus: str
    quality: str


def read_fastq_record(handle) -> Optional[FastqRecord]:
    header = handle.readline()
    if not header:
        return None
    sequence = handle.readline()
    plus = handle.readline()
    quality = handle.readline()
    if not sequence or not plus or not quality:
        raise ValueError("Truncated FASTQ record")
    record = FastqRecord(
        header.rstrip("\r\n"),
        sequence.rstrip("\r\n").upper(),
        plus.rstrip("\r\n"),
        quality.rstrip("\r\n"),
    )
    if not record.header.startswith("@"):
        raise ValueError(f"Invalid FASTQ header: {record.header[:80]}")
    if not record.plus.startswith("+"):
        raise ValueError(f"Invalid FASTQ plus line for {record.header[:80]}")
    if len(record.sequence) != len(record.quality):
        raise ValueError(
            f"FASTQ sequence/quality length mismatch for {record.header[:80]}: "
            f"{len(record.sequence)} != {len(record.quality)}"
        )
    return record


def iter_fastq_pairs(r1_path: str, r2_path: Optional[str] = None):
    """Yield synchronized FASTQ records, validating IDs when R2 is present."""
    with open_text(r1_path, "rt") as r1_handle:
        if r2_path is None:
            while True:
                r1 = read_fastq_record(r1_handle)
                if r1 is None:
                    return
                yield r1, None
        else:
            with open_text(r2_path, "rt") as r2_handle:
                record_number = 0
                while True:
                    r1 = read_fastq_record(r1_handle)
                    r2 = read_fastq_record(r2_handle)
                    if r1 is None and r2 is None:
                        return
                    if r1 is None or r2 is None:
                        raise ValueError(
                            f"FASTQ length mismatch: {r1_path} versus {r2_path}"
                        )
                    record_number += 1
                    id1 = normalize_read_id(r1.header)
                    id2 = normalize_read_id(r2.header)
                    if id1 != id2:
                        raise ValueError(
                            f"FASTQ ID mismatch at record {record_number}: {id1} != {id2}"
                        )
                    yield r1, r2


def iter_fastq_records(path: str):
    with open_text(path, "rt") as handle:
        while True:
            record = read_fastq_record(handle)
            if record is None:
                return
            yield record


def normalize_fastq_paths(paths) -> List[str]:
    """Normalize one path or an ordered path collection without splitting strings."""
    if isinstance(paths, (str, os.PathLike)):
        normalized = [str(paths)]
    else:
        normalized = [str(path) for path in paths]
    if not normalized:
        raise ValueError("at least one FASTQ path is required")
    return normalized


def iter_fastq_records_across(paths):
    """Treat multiple FASTQ lanes as one ordered logical record stream."""
    for path in normalize_fastq_paths(paths):
        yield from iter_fastq_records(path)


def _open_fastq_binary(path: str):
    return gzip.open(path, "rb") if str(path).endswith(".gz") else open(path, "rb")


def _read_fastq_record_binary(handle, path: str, record_number: int):
    """Read and validate one FASTQ record while decoding only its sequence."""

    header = handle.readline()
    if not header:
        return None
    sequence = handle.readline()
    plus = handle.readline()
    quality = handle.readline()
    if not sequence or not plus or not quality:
        raise ValueError(f"Truncated FASTQ record {record_number} in {path}")
    header = header.rstrip(b"\r\n")
    sequence = sequence.rstrip(b"\r\n").upper()
    quality = quality.rstrip(b"\r\n")
    if not header.startswith(b"@"):
        raise ValueError(f"Invalid FASTQ header in {path}: {header[:80]!r}")
    if not plus.startswith(b"+"):
        raise ValueError(f"Invalid FASTQ plus line in {path}: {header[:80]!r}")
    if len(sequence) != len(quality):
        raise ValueError(
            f"FASTQ sequence/quality mismatch at record {record_number} in {path}"
        )
    read_id = header.split(maxsplit=1)[0].lstrip(b"@")
    if read_id.endswith((b"/1", b"/2")):
        read_id = read_id[:-2]
    return read_id, sequence


@dataclass(frozen=True, slots=True)
class BinaryFastqRecord:
    header: bytes
    sequence: bytes
    plus: bytes
    quality: bytes

    def as_text(self) -> FastqRecord:
        return FastqRecord(
            self.header.decode("ascii"),
            self.sequence.decode("ascii"),
            self.plus.decode("ascii"),
            self.quality.decode("ascii"),
        )


def _read_full_fastq_record_binary(handle, path: str, record_number: int):
    header = handle.readline()
    if not header:
        return None
    sequence = handle.readline()
    plus = handle.readline()
    quality = handle.readline()
    if not sequence or not plus or not quality:
        raise ValueError(f"Truncated FASTQ record {record_number} in {path}")
    record = BinaryFastqRecord(
        header.rstrip(b"\r\n"),
        sequence.rstrip(b"\r\n").upper(),
        plus.rstrip(b"\r\n"),
        quality.rstrip(b"\r\n"),
    )
    if not record.header.startswith(b"@"):
        raise ValueError(f"Invalid FASTQ header in {path}: {record.header[:80]!r}")
    if not record.plus.startswith(b"+"):
        raise ValueError(f"Invalid FASTQ plus line in {path}: {record.header[:80]!r}")
    if len(record.sequence) != len(record.quality):
        raise ValueError(
            f"FASTQ sequence/quality mismatch at record {record_number} in {path}"
        )
    read_id = record.header.split(maxsplit=1)[0].lstrip(b"@")
    if read_id.endswith((b"/1", b"/2")):
        read_id = read_id[:-2]
    return read_id, record


def iter_fastq_pairs_binary_raw(r1_path: str, r2_path: str):
    """Yield paired raw bytes while validating the complete FASTQ contract."""

    with _open_fastq_binary(r1_path) as r1_handle, _open_fastq_binary(
        r2_path
    ) as r2_handle:
        record_number = 0
        while True:
            next_number = record_number + 1
            r1 = _read_full_fastq_record_binary(r1_handle, r1_path, next_number)
            r2 = _read_full_fastq_record_binary(r2_handle, r2_path, next_number)
            if r1 is None and r2 is None:
                return
            if r1 is None or r2 is None:
                raise ValueError(f"FASTQ length mismatch: {r1_path} versus {r2_path}")
            record_number = next_number
            if r1[0] != r2[0]:
                raise ValueError(
                    f"FASTQ ID mismatch at record {record_number}: "
                    f"{r1[0]!r} != {r2[0]!r}"
                )
            yield r1[0], r1[1], r2[1]


def iter_fastq_sequence_pairs(r1_path: str, r2_path: str):
    """Yield validated uppercase sequence pairs with minimal allocations."""

    with _open_fastq_binary(r1_path) as r1_handle, _open_fastq_binary(
        r2_path
    ) as r2_handle:
        record_number = 0
        while True:
            next_number = record_number + 1
            r1 = _read_fastq_record_binary(r1_handle, r1_path, next_number)
            r2 = _read_fastq_record_binary(r2_handle, r2_path, next_number)
            if r1 is None and r2 is None:
                return
            if r1 is None or r2 is None:
                raise ValueError(f"FASTQ length mismatch: {r1_path} versus {r2_path}")
            record_number = next_number
            if r1[0] != r2[0]:
                raise ValueError(
                    f"FASTQ ID mismatch at record {record_number}: "
                    f"{r1[0]!r} != {r2[0]!r}"
                )
            yield r1[1].decode("ascii"), r2[1].decode("ascii")


def write_fastq_record(handle, rec: FastqRecord, length: Optional[int] = None):
    seq = rec.sequence if length is None else rec.sequence[:length]
    qual = rec.quality if length is None else rec.quality[:length]
    handle.write(f"{rec.header}\n{seq}\n{rec.plus}\n{qual}\n")


def splitpipe_processed_record(
    rec: FastqRecord,
    tso_sequence: str = "AACGCAGAGTGAATGGG",
    tso_start_limit: int = 15,
    raw_length: Optional[int] = None,
) -> Tuple[str, str, int]:
    """Emulate preprocessing of raw R1 immediately before STAR."""
    sequence = rec.sequence if raw_length is None else rec.sequence[:raw_length]
    quality = rec.quality if raw_length is None else rec.quality[:raw_length]
    index = sequence.find(
        tso_sequence,
        0,
        tso_start_limit + len(tso_sequence),
    )
    if index >= 0:
        trim_end = index + len(tso_sequence)
        sequence = sequence[trim_end:]
        quality = quality[trim_end:]
    return sequence, quality, index


def sanitize_qname_field(value: str) -> str:
    return str(value).replace("|", "_").replace(" ", "_").replace("\t", "_")


def write_alignment_truth_record(
    handle,
    truth_index: int,
    sequence: str,
    quality: str,
    barcode: BarcodeCall,
    evidence: FullEvidence,
    tso_trim_index: int,
    production_removed: bool,
):
    targets = ",".join(evidence.target_genes) if evidence.target_genes else "UNKNOWN"
    guide_query = evidence.guide_query or "UNKNOWN"
    header = (
        f"@GQ{truth_index:09d}"
        f"|C={sanitize_qname_field(barcode.cell)}"
        f"|B={sanitize_qname_field(barcode.barcode_indices)}"
        f"|U={sanitize_qname_field(barcode.umi)}"
        f"|T={sanitize_qname_field(targets)}"
        f"|Q={sanitize_qname_field(guide_query)}"
        f"|X={int(tso_trim_index >= 0)}"
        f"|O={evidence.tso_start}"
        f"|R={int(production_removed)}"
    )
    if len(header) > 250:
        raise ValueError(f"Alignment truth FASTQ header too long: {header}")
    handle.write(f"{header}\n{sequence}\n+\n{quality}\n")


def find_tso(
    sequence: str,
    tso: str,
    max_dist: int,
    search_end: int = 50,
) -> Tuple[int, int, int]:
    """Return (start, exclusive end, edit distance), or (-1, -1, -1)."""
    sequence = sequence.upper()
    tso = tso.upper()
    target = sequence[: min(len(sequence), search_end)]
    exact_start = target.find(tso)
    if 0 <= exact_start <= 20:
        return exact_start, exact_start + len(tso), 0
    if HAS_EDLIB:
        result = edlib.align(tso, target, mode="HW", k=max_dist, task="locations")
        if result["editDistance"] < 0 or not result["locations"]:
            return -1, -1, -1
        # Prefer the leftmost location; guide contaminants are expected to start
        # with the TSO, while still permitting small library-prep offsets.
        start, end = sorted(result["locations"], key=lambda loc: (loc[0], loc[1]))[0]
        return start, end + 1, int(result["editDistance"])

    # Dependency-free fallback. Partitioning the query into k+1 exact seeds is
    # lossless for edit distance <=k (pigeonhole principle), and avoids a Python
    # base-by-base scan for virtually every read in a >10^9-read pilot. Candidate
    # starts are expanded by k to accommodate insertions/deletions before a seed.
    def levenshtein(a: str, b: str, cutoff: int) -> int:
        if abs(len(a) - len(b)) > cutoff:
            return cutoff + 1
        previous = list(range(len(b) + 1))
        for row_index, a_base in enumerate(a, start=1):
            current = [row_index]
            row_min = row_index
            for col_index, b_base in enumerate(b, start=1):
                value = min(
                    current[-1] + 1,
                    previous[col_index] + 1,
                    previous[col_index - 1] + int(a_base != b_base),
                )
                current.append(value)
                row_min = min(row_min, value)
            if row_min > cutoff:
                return cutoff + 1
            previous = current
        return previous[-1]

    n_seeds = max_dist + 1
    candidate_starts: Set[int] = set()
    for seed_number in range(n_seeds):
        seed_start = len(tso) * seed_number // n_seeds
        seed_end = len(tso) * (seed_number + 1) // n_seeds
        seed = tso[seed_start:seed_end]
        observed_at = target.find(seed)
        while observed_at >= 0:
            nominal_start = observed_at - seed_start
            for shift in range(-max_dist, max_dist + 1):
                start = nominal_start + shift
                if 0 <= start <= 20:
                    candidate_starts.add(start)
            observed_at = target.find(seed, observed_at + 1)

    substitution_best = None
    for start in sorted(candidate_starts):
        observed = target[start : start + len(tso)]
        if len(observed) != len(tso):
            continue
        distance = hamming_distance(observed, tso, max_dist)
        candidate = (distance, start)
        if distance <= max_dist and (
            substitution_best is None or candidate < substitution_best
        ):
            substitution_best = candidate
    if substitution_best is not None:
        return (
            substitution_best[1],
            substitution_best[1] + len(tso),
            substitution_best[0],
        )

    best = None
    for start in sorted(candidate_starts):
        for length_delta in range(-max_dist, max_dist + 1):
            observed_length = len(tso) + length_delta
            if observed_length <= 0 or start + observed_length > len(target):
                continue
            observed = target[start : start + observed_length]
            distance = levenshtein(tso, observed, max_dist)
            candidate = (distance, start, start + observed_length)
            if distance <= max_dist and (best is None or candidate < best):
                best = candidate
    if best is None:
        return -1, -1, -1
    return best[1], best[2], best[0]


BASE_TO_BITS = {"A": 0, "C": 1, "G": 2, "T": 3}


def _encode_dna(sequence: str) -> int:
    packed = 0
    for base in sequence:
        try:
            code = BASE_TO_BITS[base]
        except KeyError as exc:
            raise ValueError(
                f"packed QC guide index requires A/C/G/T only: {sequence!r}"
            ) from exc
        packed = (packed << 2) | code
    return packed


def _packed_hamming(left: int, right: int, length: int) -> int:
    different = left ^ right
    different |= different >> 1
    low_bits = ((1 << (2 * length)) - 1) // 3
    return (different & low_bits).bit_count()


def _position_mask(length: int, start: int, end: int) -> int:
    mask = 0
    for position in range(start, end):
        mask |= 0b11 << (2 * (length - position - 1))
    return mask


@dataclass(frozen=True)
class PackedIndexStats:
    sequence_count: int
    lengths: Tuple[int, ...]
    seed_key_counts: Dict[int, Tuple[int, int, int]]
    maximum_seed_bucket_sizes: Dict[int, Tuple[int, int, int]]
    total_seed_memberships: int


class HammingGuideIndex:
    """Exact plus three-seed full-verification index for Hamming distance <=2.

    At most two substituted positions leave at least one of three disjoint
    seeds unchanged.  Memory therefore grows linearly with the reference
    instead of materializing every distance-two neighbor.
    """

    def __init__(self, sequences: Iterable[str], max_dist: int = 2):
        if max_dist not in (0, 1, 2):
            raise ValueError("HammingGuideIndex supports max_dist 0, 1, or 2")
        self.max_dist = max_dist
        canonical = sorted({str(value).upper() for value in sequences if str(value)})
        if not canonical:
            raise ValueError("packed guide index requires at least one sequence")
        if any(any(base not in BASE_TO_BITS for base in sequence) for sequence in canonical):
            raise ValueError("packed QC guide index requires A/C/G/T only")
        self.sequences = tuple(canonical)
        self._packed = tuple(_encode_dna(sequence) for sequence in self.sequences)
        self._indices_by_length: Dict[int, List[int]] = defaultdict(list)
        for index, sequence in enumerate(self.sequences):
            self._indices_by_length[len(sequence)].append(index)
        self._exact: Dict[int, Dict[int, int]] = {}
        self._masks: Dict[int, Tuple[int, int, int]] = {}
        self._seeds: Dict[int, Tuple[Dict[int, Tuple[int, ...]], ...]] = {}
        seed_key_counts: Dict[int, Tuple[int, int, int]] = {}
        maximums: Dict[int, Tuple[int, int, int]] = {}
        memberships = 0
        for length, indices in self._indices_by_length.items():
            boundaries = tuple((part * length) // 3 for part in range(4))
            masks = tuple(
                _position_mask(length, boundaries[part], boundaries[part + 1])
                for part in range(3)
            )
            mutable = tuple(defaultdict(list) for _ in range(3))
            exact: Dict[int, int] = {}
            for index in indices:
                packed = self._packed[index]
                exact[packed] = index
                for part, mask in enumerate(masks):
                    mutable[part][packed & mask].append(index)
                    memberships += 1
            frozen = tuple(
                {key: tuple(value) for key, value in table.items()}
                for table in mutable
            )
            self._exact[length] = exact
            self._masks[length] = masks  # type: ignore[assignment]
            self._seeds[length] = frozen
            seed_key_counts[length] = tuple(len(table) for table in frozen)  # type: ignore[assignment]
            maximums[length] = tuple(
                max((len(value) for value in table.values()), default=0)
                for table in frozen
            )  # type: ignore[assignment]
        self.stats = PackedIndexStats(
            sequence_count=len(self.sequences),
            lengths=tuple(sorted(self._indices_by_length)),
            seed_key_counts=seed_key_counts,
            maximum_seed_bucket_sizes=maximums,
            total_seed_memberships=memberships,
        )

    def match(
        self, query: str, max_dist: Optional[int] = None
    ) -> Tuple[Optional[int], Set[str]]:
        max_dist = self.max_dist if max_dist is None else int(max_dist)
        if max_dist < 0 or max_dist > self.max_dist:
            raise ValueError(
                f"query max_dist must be between 0 and initialized {self.max_dist}"
            )
        query = str(query).upper()
        length = len(query)
        if length not in self._exact or any(base not in BASE_TO_BITS for base in query):
            return None, set()
        packed = _encode_dna(query)
        exact_index = self._exact[length].get(packed)
        if exact_index is not None:
            return 0, {self.sequences[exact_index]}
        if max_dist == 0:
            return None, set()
        candidates: Set[int] = set()
        for mask, table in zip(self._masks[length], self._seeds[length]):
            candidates.update(table.get(packed & mask, ()))
        best = max_dist + 1
        matches: Set[str] = set()
        for index in candidates:
            distance = _packed_hamming(packed, self._packed[index], length)
            if distance < best:
                best = distance
                matches = {self.sequences[index]}
            elif distance == best:
                matches.add(self.sequences[index])
        if best > max_dist:
            return None, set()
        return best, matches


class FixedLengthBarcodeEditIndex:
    """Exact Levenshtein index for fixed-length DNA barcodes at distance <=2.

    split-pipe corrects barcodes with edit distance, not Hamming distance.  For
    equal-length reference and observed barcodes, every distance-two result is
    either two substitutions or one deletion plus one insertion.  Enumerating
    those bounded neighborhoods once avoids an edlib call for every read while
    preserving split-pipe's insertion/deletion rescues.
    """

    def __init__(self, sequences: Iterable[str], max_dist: int = 2):
        if max_dist not in (0, 1, 2):
            raise ValueError(
                "FixedLengthBarcodeEditIndex supports max_dist 0, 1, or 2"
            )
        canonical = sorted({str(value).upper() for value in sequences if str(value)})
        if not canonical:
            raise ValueError("barcode edit index requires at least one sequence")
        lengths = {len(sequence) for sequence in canonical}
        if len(lengths) != 1:
            raise ValueError("barcode edit index requires one fixed barcode length")
        if any(any(base not in BASE_TO_BITS for base in sequence) for sequence in canonical):
            raise ValueError("barcode edit index requires A/C/G/T only")
        self.max_dist = int(max_dist)
        self.length = lengths.pop()
        neighborhoods: Dict[str, Tuple[int, Set[str]]] = {}

        def add(query: str, distance: int, reference: str) -> None:
            prior = neighborhoods.get(query)
            if prior is None or distance < prior[0]:
                neighborhoods[query] = (distance, {reference})
            elif distance == prior[0]:
                prior[1].add(reference)

        for reference in canonical:
            add(reference, 0, reference)
            if self.max_dist == 0:
                continue
            for first in range(self.length):
                for base in DNA:
                    if base != reference[first]:
                        add(
                            reference[:first] + base + reference[first + 1 :],
                            1,
                            reference,
                        )
            if self.max_dist < 2:
                continue
            for first in range(self.length):
                for second in range(first + 1, self.length):
                    for first_base in DNA:
                        if first_base == reference[first]:
                            continue
                        for second_base in DNA:
                            if second_base == reference[second]:
                                continue
                            query = list(reference)
                            query[first] = first_base
                            query[second] = second_base
                            add("".join(query), 2, reference)
            for deleted_position in range(self.length):
                shortened = (
                    reference[:deleted_position] + reference[deleted_position + 1 :]
                )
                for inserted_position in range(self.length):
                    for base in DNA:
                        query = (
                            shortened[:inserted_position]
                            + base
                            + shortened[inserted_position:]
                        )
                        add(query, 2, reference)
        self._neighborhoods = neighborhoods

    def match(
        self, query: str, max_dist: Optional[int] = None
    ) -> Tuple[Optional[int], Set[str]]:
        max_dist = self.max_dist if max_dist is None else int(max_dist)
        if max_dist < 0 or max_dist > self.max_dist:
            raise ValueError(
                f"query max_dist must be between 0 and initialized {self.max_dist}"
            )
        query = str(query).upper()
        if len(query) != self.length or any(base not in BASE_TO_BITS for base in query):
            return None, set()
        result = self._neighborhoods.get(query)
        if result is None or result[0] > max_dist:
            return None, set()
        return result[0], set(result[1])


class HyperLogLog:
    """Small dependency-free cardinality sketch for streaming molecule counts."""

    def __init__(self, precision: int = 18):
        if not 10 <= precision <= 20:
            raise ValueError("HyperLogLog precision must be between 10 and 20")
        self.precision = precision
        self.n_registers = 1 << precision
        self.registers = bytearray(self.n_registers)

    @property
    def relative_standard_error(self) -> float:
        return 1.04 / math.sqrt(self.n_registers)

    def add(self, *parts):
        payload = "\x1f".join(str(part) for part in parts).encode("utf-8")
        value = int.from_bytes(
            hashlib.blake2b(payload, digest_size=8).digest(), "big"
        )
        suffix_bits = 64 - self.precision
        register = value >> suffix_bits
        suffix = value & ((1 << suffix_bits) - 1)
        rank = suffix_bits - suffix.bit_length() + 1 if suffix else suffix_bits + 1
        if rank > self.registers[register]:
            self.registers[register] = rank

    def estimate(self) -> float:
        m = self.n_registers
        if m == 16:
            alpha = 0.673
        elif m == 32:
            alpha = 0.697
        elif m == 64:
            alpha = 0.709
        else:
            alpha = 0.7213 / (1 + 1.079 / m)
        raw = alpha * m * m / sum(2.0 ** (-register) for register in self.registers)
        zero_registers = self.registers.count(0)
        if raw <= 2.5 * m and zero_registers:
            return m * math.log(m / zero_registers)
        return raw


@dataclass
class BarcodeCall:
    valid: bool
    cell: str = ""
    barcode_indices: str = ""
    umi: str = ""
    max_distance: int = -1


class PilotBarcodeCorrector:
    """Reproduce v3 split-pipe barcode correction from retained barcode data."""

    def __init__(
        self,
        barcode_df: pd.DataFrame,
        amplicon: str,
        max_mismatches: int = 2,
        perfect_cell_fraction: float = 0.92,
    ):
        required = {"bci", "sequence", "wind", "bc_round"}
        missing = required - set(barcode_df.columns)
        if missing:
            raise ValueError(f"barcode_data.csv missing columns: {sorted(missing)}")
        self.amplicon = amplicon
        self.max_mismatches = max_mismatches
        self.perfect_cell_fraction = perfect_cell_fraction
        self.spans = {}
        for symbol in ("N", "1", "2", "3"):
            positions = [i for i, value in enumerate(amplicon) if value == symbol]
            if not positions:
                raise ValueError(f"Amplicon has no '{symbol}' segment: {amplicon}")
            self.spans[symbol] = (min(positions), max(positions) + 1)

        self.sequences_by_round: Dict[int, Set[str]] = {}
        self.wind_by_round_sequence: Dict[Tuple[int, str], str] = {}
        self.bci_by_round_sequence: Dict[Tuple[int, str], str] = {}
        self.index_by_round: Dict[int, FixedLengthBarcodeEditIndex] = {}
        for bc_round in (1, 2, 3):
            frame = barcode_df[barcode_df["bc_round"].astype(int) == bc_round]
            sequences = set(frame["sequence"].astype(str).str.upper())
            if not sequences:
                raise ValueError(f"No round-{bc_round} barcodes in barcode_data.csv")
            self.sequences_by_round[bc_round] = sequences
            self.index_by_round[bc_round] = FixedLengthBarcodeEditIndex(
                sequences, max_dist=max_mismatches
            )
            for row in frame.to_dict("records"):
                sequence = str(row["sequence"]).upper()
                wind = str(int(row["wind"])).zfill(2)
                prior = self.wind_by_round_sequence.get((bc_round, sequence))
                if prior is not None and prior != wind:
                    raise ValueError(
                        f"Barcode {sequence} maps to multiple round-{bc_round} wells"
                    )
                self.wind_by_round_sequence[(bc_round, sequence)] = wind
                self.bci_by_round_sequence[(bc_round, sequence)] = str(
                    int(row["bci"])
                )

        self.perfect_combination_counts = Counter()
        self.perfect_count_threshold = 5

    def raw_parts(self, sequence: str) -> Tuple[str, str, str, str]:
        sequence = sequence.upper()
        if len(sequence) < len(self.amplicon):
            return "", "", "", ""
        umi = sequence[slice(*self.spans["N"])]
        bc1 = sequence[slice(*self.spans["1"])]
        bc2 = sequence[slice(*self.spans["2"])]
        bc3 = sequence[slice(*self.spans["3"])]
        return umi, bc1, bc2, bc3

    def fit_perfect_counts(
        self,
        r2_paths,
        sample_start: int = 3_000_000,
        sample_end: int = 4_000_000,
    ):
        counts = Counter()
        for index, record in enumerate(iter_fastq_records_across(r2_paths)):
            if index < sample_start:
                continue
            if index >= sample_end:
                break
            umi, bc1, bc2, bc3 = self.raw_parts(record.sequence)
            if (
                umi
                and bc1 in self.sequences_by_round[1]
                and bc2 in self.sequences_by_round[2]
                and bc3 in self.sequences_by_round[3]
            ):
                counts[(bc1, bc2, bc3)] += 1
        self.perfect_combination_counts = counts
        if counts:
            ordered = sorted(counts.values(), reverse=True)
            total = sum(ordered)
            cumulative = 0
            threshold = ordered[-1]
            best_delta = float("inf")
            for value in ordered:
                cumulative += value
                delta = abs(cumulative / total - self.perfect_cell_fraction)
                if delta < best_delta:
                    best_delta = delta
                    threshold = value
            self.perfect_count_threshold = max(5, int(threshold))

    def call(self, sequence: str) -> BarcodeCall:
        umi, raw1, raw2, raw3 = self.raw_parts(sequence)
        if not umi or "N" in umi + raw1 + raw2 + raw3:
            return BarcodeCall(False, umi=umi)
        matches = []
        distances = []
        for bc_round, raw in ((1, raw1), (2, raw2), (3, raw3)):
            distance, candidates = self.index_by_round[bc_round].match(
                raw, max_dist=self.max_mismatches
            )
            if distance is None or not candidates:
                return BarcodeCall(False, umi=umi)
            matches.append(sorted(candidates))
            distances.append(int(distance))

        if max(distances) == 0:
            corrected = [values[0] for values in matches]
        else:
            supported = []
            for combination in itertools.product(*matches):
                if (
                    self.perfect_combination_counts.get(combination, 0)
                    > self.perfect_count_threshold
                ):
                    supported.append(combination)
                    if len(supported) > 1:
                        break
            if len(supported) != 1:
                return BarcodeCall(False, umi=umi, max_distance=max(distances))
            corrected = list(supported[0])

        winds = [
            self.wind_by_round_sequence[(bc_round, sequence)]
            for bc_round, sequence in zip((1, 2, 3), corrected)
        ]
        barcode_indices = [
            self.bci_by_round_sequence[(bc_round, sequence)]
            for bc_round, sequence in zip((1, 2, 3), corrected)
        ]
        return BarcodeCall(
            True,
            cell="_".join(winds),
            barcode_indices="_".join(barcode_indices),
            umi=umi,
            max_distance=max(distances),
        )


@dataclass
class GuideMatch:
    length: int
    distance: int
    references: Set[str]
    query: str


@dataclass
class TemplatePrefixMatch:
    length: int
    distance: int
    prefixes: Set[str]
    guide_sequences: Set[str]
    guide_names: Set[str]
    target_genes: Set[str]
    variants: Set[int]
    minimum_boundary_bases: int
    identity_ambiguous: bool


@dataclass
class FullEvidence:
    tso_found: bool = False
    tso_start: int = -1
    tso_end: int = -1
    tso_distance: int = -1
    structure_truth: bool = False
    known_structure_truth: bool = False
    complete_reference_truth: bool = False
    tso_known_guide_truth: bool = False
    tier: str = "no_tso"
    guide_query: str = ""
    guide_length: int = 0
    guide_distance: int = -1
    guide_references: Tuple[str, ...] = ()
    guide_names: Tuple[str, ...] = ()
    target_genes: Tuple[str, ...] = ()
    variant: int = 0
    stem_mismatches: int = -1
    tracr_mismatches: int = -1
    ibar_query: str = ""
    ibar_distance: int = -1


class GuideModel:
    def __init__(self, ref_df: pd.DataFrame, structure_cfg: dict, qc_cfg: dict):
        required = {"grna_name", "grna_seq", "ibar_seq", "position", "target_gene"}
        missing = sorted(required - set(ref_df.columns))
        if missing:
            raise ValueError(f"Guide reference missing columns: {missing}")

        self.ref_df = ref_df.copy()
        for col in ("grna_seq", "ibar_seq"):
            self.ref_df[col] = self.ref_df[col].fillna("").astype(str).str.upper()
        self.ref_df["position"] = self.ref_df["position"].astype(int)

        self.tso = str(structure_cfg["fixed_seq_1"]).upper()
        self.stems = {
            1: str(structure_cfg["fixed_seq_2"]).upper(),
            2: str(structure_cfg.get("fixed_seq_3", "")).upper(),
        }
        self.tracrs = {
            1: str(structure_cfg.get("fixed_seq_4", "")).upper(),
            2: str(structure_cfg.get("fixed_seq_5", "")).upper(),
        }
        self.grna_lengths = sorted({int(x) for x in structure_cfg["grna_lengths"]})
        self.ibar_length = int(structure_cfg["ibar_length"])
        self.tso_max_dist = int(qc_cfg.get("tso_max_edit_dist", 3))
        self.anchor_max_mismatches = int(qc_cfg.get("anchor_max_mismatches", 2))
        self.guide_max_mismatches = int(qc_cfg.get("guide_max_mismatches", 2))
        self.ibar_max_mismatches = int(qc_cfg.get("ibar_max_mismatches", 2))
        self.tso_search_end = int(qc_cfg.get("tso_search_end", 50))
        self.splitpipe_tso = str(
            qc_cfg.get("splitpipe_tso_seq", "AACGCAGAGTGAATGGG")
        ).upper()
        self.splitpipe_tso_start_limit = int(qc_cfg.get("splitpipe_tso_start_limit", 15))

        guide_sequences = self.ref_df["grna_seq"].drop_duplicates().tolist()
        self.guide_index = HammingGuideIndex(guide_sequences, max_dist=2)
        self.rows_by_sequence: Dict[str, List[dict]] = defaultdict(list)
        for row in self.ref_df.to_dict("records"):
            self.rows_by_sequence[row["grna_seq"]].append(row)
        self.guide_prefixes = {
            k: {seq[:k] for seq in guide_sequences if len(seq) >= k}
            for k in qc_cfg.get("guide_prefix_lengths", [12, 15, 18])
        }

        # At 54 cycles, a TSO offset can leave fewer than 20 spacer bases.  Match
        # the entire observed suffix against expected guide+stem prefixes so
        # 20/21-nt guides and their scaffold boundary are evaluated jointly.
        self.template_prefix_min_length = int(
            qc_cfg.get("template_prefix_min_length", 15)
        )
        self.template_prefix_max_length = int(
            qc_cfg.get("template_prefix_max_length", 24)
        )
        self.template_prefix_max_mismatches = int(
            qc_cfg.get("template_prefix_max_mismatches", 2)
        )
        if self.template_prefix_max_mismatches not in (0, 1, 2):
            raise ValueError("template_prefix_max_mismatches must be 0, 1, or 2")
        self.template_metadata_by_prefix: Dict[str, List[dict]] = defaultdict(list)
        prefixes: Set[str] = set()
        unique_templates = set()
        for row in self.ref_df.to_dict("records"):
            variant = int(row["position"])
            stem = self.stems.get(variant, "")
            if not stem:
                continue
            guide = str(row["grna_seq"])
            template = guide + stem
            template_key = (template, guide, str(row["grna_name"]), variant)
            if template_key in unique_templates:
                continue
            unique_templates.add(template_key)
            max_length = min(len(template), self.template_prefix_max_length)
            for length in range(self.template_prefix_min_length, max_length + 1):
                prefix = template[:length]
                prefixes.add(prefix)
                self.template_metadata_by_prefix[prefix].append(
                    {
                        "grna_seq": guide,
                        "grna_name": str(row["grna_name"]),
                        "target_gene": str(row["target_gene"]),
                        "position": variant,
                        "guide_length": len(guide),
                    }
                )
        self.template_prefix_index = HammingGuideIndex(
            prefixes, max_dist=self.template_prefix_max_mismatches
        )

    def best_guide_match(self, sequence_after_tso: str) -> Optional[GuideMatch]:
        hits: List[GuideMatch] = []
        for length in self.grna_lengths:
            if len(sequence_after_tso) < length:
                continue
            query = sequence_after_tso[:length]
            distance, references = self.guide_index.match(
                query, max_dist=self.guide_max_mismatches
            )
            if distance is not None:
                hits.append(GuideMatch(length, distance, references, query))
        if not hits:
            return None
        best_distance = min(hit.distance for hit in hits)
        best = [hit for hit in hits if hit.distance == best_distance]
        # Prefer the hit with more observed stem bases, i.e. shorter guide when
        # edit distances tie; downstream stem consistency still guards length.
        best.sort(key=lambda hit: hit.length)
        return best[0]

    def best_template_prefix_match(
        self, sequence_after_tso: str
    ) -> Optional[TemplatePrefixMatch]:
        observed_length = min(
            len(sequence_after_tso), self.template_prefix_max_length
        )
        if observed_length < self.template_prefix_min_length:
            return None
        query = sequence_after_tso[:observed_length]
        distance, prefixes = self.template_prefix_index.match(
            query, max_dist=self.template_prefix_max_mismatches
        )
        if distance is None:
            return None
        metadata = [
            row
            for prefix in prefixes
            for row in self.template_metadata_by_prefix.get(prefix, [])
        ]
        guide_sequences = {row["grna_seq"] for row in metadata}
        guide_names = {row["grna_name"] for row in metadata}
        targets = {row["target_gene"] for row in metadata}
        variants = {int(row["position"]) for row in metadata}
        boundary_bases = [
            max(0, observed_length - int(row["guide_length"])) for row in metadata
        ]
        minimum_boundary_bases = min(boundary_bases) if boundary_bases else 0
        identities = {
            (row["grna_seq"], row["grna_name"], int(row["position"]))
            for row in metadata
        }
        return TemplatePrefixMatch(
            length=observed_length,
            distance=int(distance),
            prefixes=set(prefixes),
            guide_sequences=guide_sequences,
            guide_names=guide_names,
            target_genes=targets,
            variants=variants,
            minimum_boundary_bases=minimum_boundary_bases,
            identity_ambiguous=len(identities) > 1,
        )

    def candidate_rows(self, references: Iterable[str], variant: Optional[int] = None):
        rows = []
        for reference in references:
            for row in self.rows_by_sequence.get(reference, []):
                if variant is None or int(row["position"]) == int(variant):
                    rows.append(row)
        return rows

    def best_ibar_distance(self, query: str, rows: Sequence[dict]) -> Optional[int]:
        distances = []
        for row in rows:
            ibar = row["ibar_seq"]
            if ibar and len(ibar) == len(query):
                distances.append(hamming_distance(query, ibar, self.ibar_max_mismatches))
        if not distances:
            return None
        best = min(distances)
        return best if best <= self.ibar_max_mismatches else None

    def classify_full(self, sequence: str) -> FullEvidence:
        sequence = sequence.upper()
        evidence = FullEvidence()
        tso_start, tso_end, tso_distance = find_tso(
            sequence, self.tso, self.tso_max_dist, self.tso_search_end
        )
        if tso_start < 0:
            return evidence
        evidence.tso_found = True
        evidence.tso_start = tso_start
        evidence.tso_end = tso_end
        evidence.tso_distance = tso_distance
        evidence.tier = "tso_only"

        # Inclusive reference-supported tier.  This is deliberately not the
        # primary truth because a short genuine cDNA may also start with TSO;
        # complete downstream structure remains the central surrogate truth.
        possible_guide = self.best_guide_match(sequence[tso_end:])
        if possible_guide is not None:
            evidence.tso_known_guide_truth = True
            evidence.guide_query = possible_guide.query
            evidence.guide_length = possible_guide.length
            evidence.guide_distance = possible_guide.distance
            evidence.guide_references = tuple(sorted(possible_guide.references))
            possible_rows = self.candidate_rows(possible_guide.references)
            evidence.guide_names = tuple(
                sorted({str(row["grna_name"]) for row in possible_rows})
            )
            evidence.target_genes = tuple(
                sorted({str(row["target_gene"]) for row in possible_rows})
            )
            evidence.tier = "tso_known_guide_no_structure"

        candidates = []
        for guide_length in self.grna_lengths:
            guide_start = tso_end
            guide_end = guide_start + guide_length
            for variant in (1, 2):
                stem = self.stems.get(variant, "")
                tracr = self.tracrs.get(variant, "")
                if not stem or not tracr:
                    continue
                stem_start = guide_end
                stem_end = stem_start + len(stem)
                ibar_start = stem_end
                ibar_end = ibar_start + self.ibar_length
                tracr_end = ibar_end + len(tracr)
                if tracr_end > len(sequence):
                    continue
                stem_mm = hamming_distance(sequence[stem_start:stem_end], stem)
                tracr_mm = hamming_distance(sequence[ibar_end:tracr_end], tracr)
                anchor_mm = stem_mm + tracr_mm
                if anchor_mm > self.anchor_max_mismatches:
                    continue
                guide_query = sequence[guide_start:guide_end]
                guide_dist, guide_refs = self.guide_index.match(
                    guide_query, self.guide_max_mismatches
                )
                matching_rows = self.candidate_rows(guide_refs, variant) if guide_refs else []
                ibar_query = sequence[ibar_start:ibar_end]
                ibar_dist = self.best_ibar_distance(ibar_query, matching_rows)
                candidates.append(
                    (
                        anchor_mm,
                        guide_dist if guide_dist is not None else 99,
                        ibar_dist if ibar_dist is not None else 99,
                        guide_length,
                        variant,
                        stem_mm,
                        tracr_mm,
                        guide_query,
                        guide_refs,
                        matching_rows,
                        ibar_query,
                        ibar_dist,
                    )
                )

        if not candidates:
            return evidence

        candidates.sort(key=lambda item: item[:4])
        best = candidates[0]
        (
            _, guide_dist_sort, _, guide_length, variant, stem_mm, tracr_mm,
            guide_query, guide_refs, matching_rows, ibar_query, ibar_dist,
        ) = best
        evidence.structure_truth = True
        evidence.tier = "structure_unknown_guide"
        evidence.guide_query = guide_query
        evidence.guide_length = guide_length
        evidence.variant = variant
        evidence.stem_mismatches = stem_mm
        evidence.tracr_mismatches = tracr_mm
        evidence.ibar_query = ibar_query

        if guide_refs:
            evidence.known_structure_truth = True
            evidence.tso_known_guide_truth = True
            evidence.guide_distance = int(guide_dist_sort)
            evidence.guide_references = tuple(sorted(guide_refs))
            all_rows = self.candidate_rows(guide_refs)
            evidence.guide_names = tuple(sorted({str(row["grna_name"]) for row in all_rows}))
            evidence.target_genes = tuple(sorted({str(row["target_gene"]) for row in all_rows}))
            evidence.tier = "structure_known_guide"
            if matching_rows and ibar_dist is not None:
                evidence.complete_reference_truth = True
                evidence.ibar_distance = int(ibar_dist)
                evidence.tier = "complete_known_construct"
        return evidence

    def short_features(self, sequence: str, read_length: int, cfg: dict) -> dict:
        short = sequence[:read_length].upper()
        tso_start, tso_end, tso_dist = find_tso(
            short, self.tso, self.tso_max_dist, self.tso_search_end
        )
        tso_found = tso_start >= 0
        max_filter_tso_start = int(cfg.get("max_tso_start", 2))
        tso_at_expected_start = bool(
            tso_found and 0 <= tso_start <= max_filter_tso_start
        )
        post_tso = short[tso_end:] if tso_found else ""
        guide_hit = self.best_guide_match(post_tso) if tso_found else None
        template_prefix_hit = (
            self.best_template_prefix_match(post_tso) if tso_found else None
        )

        tso_start_mm = None
        if len(short) >= len(self.tso):
            tso_start_mm = hamming_distance(short[: len(self.tso)], self.tso)

        partial_stem_mm = None
        partial_stem_length = 0
        if guide_hit is not None:
            stem_observed = post_tso[guide_hit.length :]
            for reference in guide_hit.references:
                for row in self.rows_by_sequence.get(reference, []):
                    variant = int(row["position"])
                    stem = self.stems.get(variant, "")
                    observed = stem_observed[: len(stem)]
                    if not observed:
                        continue
                    mm = hamming_distance(observed, stem[: len(observed)])
                    if partial_stem_mm is None or (mm, -len(observed)) < (
                        partial_stem_mm, -partial_stem_length
                    ):
                        partial_stem_mm = mm
                        partial_stem_length = len(observed)

        template_mm = None
        if tso_start_mm is not None and guide_hit is not None and partial_stem_mm is not None:
            template_mm = tso_start_mm + guide_hit.distance + partial_stem_mm

        splitpipe_index = short.find(
            self.splitpipe_tso,
            0,
            self.splitpipe_tso_start_limit + len(self.splitpipe_tso),
        )
        splitpipe_post = (
            short[splitpipe_index + len(self.splitpipe_tso) :]
            if splitpipe_index >= 0
            else ""
        )
        splitpipe_hit = self.best_guide_match(splitpipe_post) if splitpipe_post else None

        calls: Dict[str, bool] = {
            "tso_only_d3": tso_at_expected_start and tso_dist <= 3,
            "tso_exact30": short.startswith(self.tso),
        }
        guide_thresholds = [int(x) for x in cfg.get("guide_mismatch_thresholds", [0, 1, 2])]
        min_partial = int(cfg.get("min_partial_stem_bases", 3))
        min_boundary = int(cfg.get("min_boundary_bases", 1))
        for distance in guide_thresholds:
            guide_pass = guide_hit is not None and guide_hit.distance <= distance
            calls[f"tso_guide_d{distance}"] = bool(
                tso_at_expected_start and guide_pass
            )
            calls[f"tso_guide_d{distance}_partial_stem"] = bool(
                tso_at_expected_start
                and guide_pass
                and partial_stem_length >= min_partial
                and partial_stem_mm == 0
            )
            split_pass = splitpipe_hit is not None and splitpipe_hit.distance <= distance
            calls[f"splitpipe_trim_guide_d{distance}"] = bool(splitpipe_index >= 0 and split_pass)
            template_pass = (
                template_prefix_hit is not None
                and template_prefix_hit.distance <= distance
            )
            calls[f"tso_template_prefix_d{distance}"] = bool(
                tso_at_expected_start and template_pass
            )
            calls[f"tso_template_boundary_d{distance}"] = bool(
                tso_at_expected_start
                and template_pass
                and template_prefix_hit.minimum_boundary_bases >= min_boundary
            )
            calls[f"tso_anywhere_template_prefix_d{distance}"] = bool(
                tso_found and template_pass
            )

        for threshold in [int(x) for x in cfg.get("template_mismatch_thresholds", [0, 1, 2, 3])]:
            calls[f"template_hamming_d{threshold}"] = bool(
                template_mm is not None and template_mm <= threshold
            )

        for prefix_length in [int(x) for x in cfg.get("guide_prefix_lengths", [12, 15, 18])]:
            calls[f"tso_guide_prefix_{prefix_length}"] = bool(
                tso_at_expected_start
                and len(post_tso) >= prefix_length
                and post_tso[:prefix_length] in self.guide_prefixes.get(prefix_length, set())
            )

        return {
            "short_sequence": short,
            "tso_found": tso_found,
            "tso_at_expected_start": tso_at_expected_start,
            "tso_start": tso_start,
            "tso_end": tso_end,
            "tso_distance": tso_dist,
            "tso_start_hamming": tso_start_mm,
            "guide_hit": guide_hit,
            "template_prefix_hit": template_prefix_hit,
            "partial_stem_mismatches": partial_stem_mm,
            "partial_stem_length": partial_stem_length,
            "template_hamming": template_mm,
            "splitpipe_tso_index": splitpipe_index,
            "splitpipe_guide_hit": splitpipe_hit,
            "calls": calls,
        }

    def classify_post_splitpipe(self, sequence: str) -> FullEvidence:
        """Classify a BAM query sequence, which usually has its TSO removed."""
        sequence = sequence.upper()
        raw_evidence = self.classify_full(sequence)
        if raw_evidence.structure_truth:
            return raw_evidence

        evidence = FullEvidence(tier="post_splitpipe_no_structure")
        candidates = []
        for guide_length in self.grna_lengths:
            for variant in (1, 2):
                stem = self.stems.get(variant, "")
                tracr = self.tracrs.get(variant, "")
                stem_start = guide_length
                stem_end = stem_start + len(stem)
                ibar_start = stem_end
                ibar_end = ibar_start + self.ibar_length
                tracr_end = ibar_end + len(tracr)
                if not stem or not tracr or tracr_end > len(sequence):
                    continue
                stem_mm = hamming_distance(sequence[stem_start:stem_end], stem)
                tracr_mm = hamming_distance(sequence[ibar_end:tracr_end], tracr)
                if stem_mm + tracr_mm > self.anchor_max_mismatches:
                    continue
                guide_query = sequence[:guide_length]
                guide_dist, guide_refs = self.guide_index.match(
                    guide_query, self.guide_max_mismatches
                )
                rows = self.candidate_rows(guide_refs, variant) if guide_refs else []
                ibar_query = sequence[ibar_start:ibar_end]
                ibar_dist = self.best_ibar_distance(ibar_query, rows)
                candidates.append(
                    (
                        stem_mm + tracr_mm,
                        guide_dist if guide_dist is not None else 99,
                        ibar_dist if ibar_dist is not None else 99,
                        guide_length,
                        variant,
                        stem_mm,
                        tracr_mm,
                        guide_query,
                        guide_refs,
                        rows,
                        ibar_query,
                        ibar_dist,
                    )
                )
        if not candidates:
            return evidence
        candidates.sort(key=lambda item: item[:4])
        best = candidates[0]
        evidence.structure_truth = True
        evidence.tier = "post_splitpipe_structure_unknown"
        evidence.guide_length = best[3]
        evidence.variant = best[4]
        evidence.stem_mismatches = best[5]
        evidence.tracr_mismatches = best[6]
        evidence.guide_query = best[7]
        evidence.ibar_query = best[10]
        guide_refs = best[8]
        rows = self.candidate_rows(guide_refs) if guide_refs else []
        if guide_refs:
            evidence.known_structure_truth = True
            evidence.guide_distance = int(best[1])
            evidence.guide_references = tuple(sorted(guide_refs))
            evidence.guide_names = tuple(sorted({str(row["grna_name"]) for row in rows}))
            evidence.target_genes = tuple(sorted({str(row["target_gene"]) for row in rows}))
            evidence.tier = "post_splitpipe_structure_known"
            if best[9] and best[11] is not None:
                evidence.complete_reference_truth = True
                evidence.ibar_distance = int(best[11])
                evidence.tier = "post_splitpipe_complete_known"
        return evidence


class Confusion:
    def __init__(self):
        self.tp = self.fp = self.tn = self.fn = 0
        self.hard_negative_n = self.hard_negative_fp = 0

    def update(self, truth: bool, call: bool, hard_negative: bool = False):
        if truth and call:
            self.tp += 1
        elif truth and not call:
            self.fn += 1
        elif not truth and call:
            self.fp += 1
        else:
            self.tn += 1
        if hard_negative:
            self.hard_negative_n += 1
            self.hard_negative_fp += int(call)

    def as_dict(self) -> dict:
        positive_n = self.tp + self.fn
        negative_n = self.fp + self.tn
        sensitivity = self.tp / positive_n if positive_n else float("nan")
        specificity = self.tn / negative_n if negative_n else float("nan")
        precision = self.tp / (self.tp + self.fp) if self.tp + self.fp else float("nan")
        fpr = self.fp / negative_n if negative_n else float("nan")
        f1 = (
            2 * precision * sensitivity / (precision + sensitivity)
            if precision + sensitivity > 0
            else float("nan")
        )
        sensitivity_low, sensitivity_high = wilson_interval(self.tp, positive_n)
        fpr_low, fpr_high = wilson_interval(self.fp, negative_n)
        hard_fpr_low, hard_fpr_high = wilson_interval(
            self.hard_negative_fp, self.hard_negative_n
        )
        return {
            "tp": self.tp,
            "fp": self.fp,
            "tn": self.tn,
            "fn": self.fn,
            "sensitivity": sensitivity,
            "sensitivity_95_low": sensitivity_low,
            "sensitivity_95_high": sensitivity_high,
            "specificity": specificity,
            "false_positive_rate": fpr,
            "false_positive_rate_95_low": fpr_low,
            "false_positive_rate_95_high": fpr_high,
            "precision": precision,
            "f1": f1,
            "hard_negative_n": self.hard_negative_n,
            "hard_negative_fp": self.hard_negative_fp,
            "hard_negative_fpr": (
                self.hard_negative_fp / self.hard_negative_n
                if self.hard_negative_n
                else float("nan")
            ),
            "hard_negative_fpr_95_low": hard_fpr_low,
            "hard_negative_fpr_95_high": hard_fpr_high,
        }


def wilson_interval(successes: int, total: int, z: float = 1.96) -> Tuple[float, float]:
    if total == 0:
        return float("nan"), float("nan")
    p = successes / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def percentile(values: Sequence[float], probability: float) -> float:
    values = sorted(float(value) for value in values if math.isfinite(value))
    if not values:
        return float("nan")
    position = min(1.0, max(0.0, probability)) * (len(values) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return values[lower]
    weight = position - lower
    return values[lower] * (1 - weight) + values[upper] * weight


def cluster_bootstrap_ratio_interval(
    cluster_counts: Dict[Tuple[str, str], Counter],
    numerator_key: str,
    denominator_key: str,
    replicates: int,
    seed_text: str,
) -> Tuple[float, float]:
    """Percentile interval for a ratio under corrected-CB/pB cluster sampling."""
    values = [
        (float(counts[numerator_key]), float(counts[denominator_key]))
        for counts in cluster_counts.values()
        if counts[denominator_key] > 0
    ]
    if len(values) < 2 or replicates <= 0:
        return float("nan"), float("nan")
    seed = int.from_bytes(
        hashlib.blake2b(seed_text.encode("utf-8"), digest_size=8).digest(), "big"
    )
    rng = random.Random(seed)
    n_clusters = len(values)
    ratios = []
    for _ in range(replicates):
        numerator = denominator = 0.0
        for _ in range(n_clusters):
            value_numerator, value_denominator = values[rng.randrange(n_clusters)]
            numerator += value_numerator
            denominator += value_denominator
        if denominator > 0:
            ratios.append(numerator / denominator)
    return percentile(ratios, 0.025), percentile(ratios, 0.975)


def paired_cluster_bootstrap_change(
    full_frame: pd.DataFrame,
    short_frame: pd.DataFrame,
    numerator_key: str,
    denominator_key: str,
    replicates: int,
    seed_text: str,
) -> dict:
    """Paired CB/pB-cluster bootstrap for short-minus-full ratio change."""
    keys = ["cell", "barcode_indices"]
    needed = keys + [numerator_key, denominator_key]
    if any(column not in full_frame.columns for column in needed) or any(
        column not in short_frame.columns for column in needed
    ):
        return {}
    merged = full_frame[needed].merge(
        short_frame[needed],
        on=keys,
        how="inner",
        suffixes=("_full", "_short"),
    )
    if len(merged) < 2 or replicates <= 0:
        return {"paired_sampling_units": int(len(merged))}
    values = list(
        zip(
            merged[f"{numerator_key}_full"].astype(float),
            merged[f"{denominator_key}_full"].astype(float),
            merged[f"{numerator_key}_short"].astype(float),
            merged[f"{denominator_key}_short"].astype(float),
        )
    )
    seed = int.from_bytes(
        hashlib.blake2b(seed_text.encode("utf-8"), digest_size=8).digest(), "big"
    )
    rng = random.Random(seed)
    changes = []
    folds = []
    n_clusters = len(values)
    for _ in range(replicates):
        full_numerator = full_denominator = 0.0
        short_numerator = short_denominator = 0.0
        for _ in range(n_clusters):
            fn, fd, sn, sd = values[rng.randrange(n_clusters)]
            full_numerator += fn
            full_denominator += fd
            short_numerator += sn
            short_denominator += sd
        if full_denominator <= 0 or short_denominator <= 0:
            continue
        full_ratio = full_numerator / full_denominator
        short_ratio = short_numerator / short_denominator
        changes.append(short_ratio - full_ratio)
        if full_ratio > 0:
            folds.append(short_ratio / full_ratio)
    return {
        "paired_sampling_units": int(n_clusters),
        "absolute_delta_cluster_bootstrap_95_low": percentile(changes, 0.025),
        "absolute_delta_cluster_bootstrap_95_high": percentile(changes, 0.975),
        "fold_change_cluster_bootstrap_95_low": percentile(folds, 0.025),
        "fold_change_cluster_bootstrap_95_high": percentile(folds, 0.975),
    }


class MoleculeEstimator:
    """Estimate exact-key molecule cardinalities without retaining every UMI."""

    CATEGORIES = (
        "all",
        "tso",
        "structure_contaminant",
        "known_structure_contaminant",
        "inclusive_contaminant",
    )
    KEY_TYPES = ("cell_pb_umi", "cell_pb_umi_insert")

    def __init__(self, precision: int = 18):
        self.precision = precision
        self.sketches = {
            (category, key_type): HyperLogLog(precision)
            for category in self.CATEGORIES
            for key_type in self.KEY_TYPES
        }
        self.valid_read_counts = Counter()
        self.barcode_edit_counts = Counter()
        self.cell_read_counts: Dict[str, Counter] = defaultdict(Counter)

    def add(self, barcode: BarcodeCall, sequence: str, evidence: FullEvidence):
        if not barcode.valid:
            self.valid_read_counts["invalid_barcode"] += 1
            return
        self.valid_read_counts["all"] += 1
        self.barcode_edit_counts[barcode.max_distance] += 1
        categories = ["all"]
        if evidence.tso_found:
            categories.append("tso")
            self.valid_read_counts["tso"] += 1
        if evidence.structure_truth:
            categories.append("structure_contaminant")
            self.valid_read_counts["structure_contaminant"] += 1
        if evidence.known_structure_truth:
            categories.append("known_structure_contaminant")
            self.valid_read_counts["known_structure_contaminant"] += 1
        if evidence.structure_truth or evidence.tso_known_guide_truth:
            categories.append("inclusive_contaminant")
            self.valid_read_counts["inclusive_contaminant"] += 1

        cell_counts = self.cell_read_counts[barcode.cell]
        cell_counts["all"] += 1
        cell_counts["tso"] += int(evidence.tso_found)
        cell_counts["structure_contaminant"] += int(evidence.structure_truth)
        for category in categories:
            self.sketches[(category, "cell_pb_umi")].add(
                barcode.cell, barcode.barcode_indices, barcode.umi
            )
            self.sketches[(category, "cell_pb_umi_insert")].add(
                barcode.cell, barcode.barcode_indices, barcode.umi, sequence
            )

    def _ratio(self, numerator: float, denominator: float) -> dict:
        if denominator <= 0:
            return {
                "fraction": float("nan"),
                "nominal_95_low": float("nan"),
                "nominal_95_high": float("nan"),
            }
        fraction = numerator / denominator
        # Conservative propagation treating nested HLL estimates as independent.
        ratio_rse = math.sqrt(2) * 1.04 / math.sqrt(1 << self.precision)
        half = 1.96 * ratio_rse * fraction
        return {
            "fraction": fraction,
            "nominal_95_low": max(0.0, fraction - half),
            "nominal_95_high": min(1.0, fraction + half),
        }

    def summary(self) -> dict:
        estimates = {
            key_type: {
                category: self.sketches[(category, key_type)].estimate()
                for category in self.CATEGORIES
            }
            for key_type in self.KEY_TYPES
        }
        ratios = {}
        for key_type, values in estimates.items():
            ratios[key_type] = {
                "contaminant_of_all": self._ratio(
                    values["structure_contaminant"], values["all"]
                ),
                "contaminant_of_tso": self._ratio(
                    values["structure_contaminant"], values["tso"]
                ),
                "inclusive_contaminant_of_all": self._ratio(
                    values["inclusive_contaminant"], values["all"]
                ),
                "inclusive_contaminant_of_tso": self._ratio(
                    values["inclusive_contaminant"], values["tso"]
                ),
            }
        return {
            "method": "HyperLogLog over exact corrected-cell/UMI keys",
            "precision": self.precision,
            "registers_per_sketch": 1 << self.precision,
            "nominal_relative_standard_error": 1.04
            / math.sqrt(1 << self.precision),
            "key_definitions": {
                "cell_pb_umi": (
                    "Distinct (corrected CB/well triple, corrected pB/barcode-index "
                    "triple, exact 10-nt UMI); may merge UMI collisions across genes "
                    "and does not Hamming-collapse UMI errors."
                ),
                "cell_pb_umi_insert": (
                    "Distinct (corrected CB, corrected pB, exact UMI, full raw R1 "
                    "sequence); separates UMI collisions but sequencing errors can "
                    "split duplicates."
                ),
            },
            "estimated_distinct_keys": estimates,
            "fractions": ratios,
            "valid_barcode_read_counts": {
                key: int(value) for key, value in self.valid_read_counts.items()
            },
            "barcode_max_edit_distance_counts": {
                str(key): int(value) for key, value in self.barcode_edit_counts.items()
            },
        }


def load_merged_config(config_path: str) -> dict:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)
    from qc_config_loader import load_qc_config

    return load_qc_config(config_path)


def _load_step4_module():
    path = Path(__file__).resolve().with_name("04_assign_grna_to_cells.py")
    spec = importlib.util.spec_from_file_location("csm_step4_assignment", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def resolve_path(path: str, config_path: str) -> str:
    if not path:
        return path
    if os.path.isabs(path):
        return path
    return os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(config_path)), path))


def discover_fastq_pairs(fastq_dir: str, prefix: str) -> List[Tuple[str, str]]:
    """Backward-compatible token-bounded discovery helper."""
    return discover_token_bounded_pairs(fastq_dir, prefix)


def get_library(cfg: dict, name: str) -> dict:
    matches = [lib for lib in cfg.get("libraries", []) if lib["sublibrary_name"] == name]
    if len(matches) != 1:
        raise ValueError(f"Expected one library named '{name}', found {len(matches)}")
    return matches[0]


def make_model(cfg: dict, config_path: str) -> GuideModel:
    ref_path = resolve_path(cfg["paths"]["grna_reference"], config_path)
    ref_df = pd.read_csv(ref_path, comment="#")
    qc = cfg["gex_grna_contamination_qc"]
    return GuideModel(ref_df, cfg["grna_structure"], qc["full_length_truth"])


def _fit_perfect_counts_fast(
    corrector: PilotBarcodeCorrector,
    r2_paths,
    sample_start: int,
    sample_end: int,
) -> None:
    """Fit the split-pipe barcode threshold using the binary sequence parser."""

    if sample_start < 0 or sample_end <= sample_start:
        raise ValueError(
            "barcode calibration interval must satisfy 0 <= start < end"
        )
    counts = Counter()
    logical_index = 0
    for r2_path in normalize_fastq_paths(r2_paths):
        if logical_index >= sample_end:
            break
        with _open_fastq_binary(r2_path) as handle:
            record_number = 0
            while logical_index < sample_end:
                record_number += 1
                record = _read_fastq_record_binary(handle, r2_path, record_number)
                if record is None:
                    break
                if logical_index >= sample_start:
                    sequence = record[1].decode("ascii")
                    umi, bc1, bc2, bc3 = corrector.raw_parts(sequence)
                    if (
                        umi
                        and bc1 in corrector.sequences_by_round[1]
                        and bc2 in corrector.sequences_by_round[2]
                        and bc3 in corrector.sequences_by_round[3]
                    ):
                        counts[(bc1, bc2, bc3)] += 1
                logical_index += 1
    corrector.perfect_combination_counts = counts
    if counts:
        ordered = sorted(counts.values(), reverse=True)
        total = sum(ordered)
        cumulative = 0
        threshold = ordered[-1]
        best_delta = float("inf")
        for value in ordered:
            cumulative += value
            delta = abs(cumulative / total - corrector.perfect_cell_fraction)
            if delta < best_delta:
                best_delta = delta
                threshold = value
        corrector.perfect_count_threshold = max(5, int(threshold))


def make_barcode_corrector(
    cfg: dict,
    config_path: str,
    library_name: str,
    r2_paths,
) -> PilotBarcodeCorrector:
    molecule_cfg = cfg["gex_grna_contamination_qc"]["molecule_estimation"]
    template = molecule_cfg.get(
        "barcode_data_template",
        "{output_dir}/gex_{library}/process/barcode_data.csv",
    )
    barcode_path = template.format(
        output_dir=cfg["paths"]["output_dir"], library=library_name
    )
    barcode_path = resolve_path(barcode_path, config_path)
    if not os.path.exists(barcode_path):
        raise FileNotFoundError(
            f"Molecule estimation requires split-pipe barcode data: {barcode_path}"
        )
    barcode_df = pd.read_csv(barcode_path, comment="#")
    corrector = PilotBarcodeCorrector(
        barcode_df,
        str(
            molecule_cfg.get(
                "amplicon_sequence",
                "NNNNNNNNNN33333333ATGAGGGGTCAG22222222TCCAACCACCTC11111111",
            )
        ),
        max_mismatches=int(molecule_cfg.get("barcode_max_mismatches", 2)),
        perfect_cell_fraction=float(
            molecule_cfg.get("perfect_barcode_cell_fraction", 0.92)
        ),
    )
    _fit_perfect_counts_fast(
        corrector,
        r2_paths,
        int(molecule_cfg.get("barcode_sample_start", 3_000_000)),
        int(molecule_cfg.get("barcode_sample_end", 4_000_000)),
    )
    return corrector


def initialize_base_composition(max_length: int):
    return {
        group: [Counter() for _ in range(max_length)]
        for group in ("all", "structure_truth", "tso_positive_gex_like")
    }


def update_base_composition(base_comp, group: str, sequence: str):
    if group not in base_comp:
        return
    for index, base in enumerate(sequence[: len(base_comp[group])]):
        base_comp[group][index][base] += 1


def write_base_composition(base_comp, output_path: str):
    rows = []
    for group, positions in base_comp.items():
        for index, counts in enumerate(positions, start=1):
            total = sum(counts.values())
            if not total:
                continue
            for base in "ACGTN":
                rows.append(
                    {
                        "group": group,
                        "position": index,
                        "base": base,
                        "count": counts.get(base, 0),
                        "fraction": counts.get(base, 0) / total,
                    }
                )
    pd.DataFrame(
        rows, columns=["group", "position", "base", "count", "fraction"]
    ).to_csv(output_path, index=False)


def write_rows_csv(
    rows: Sequence[dict],
    columns: Sequence[str],
    output_path: str,
    sort_by: Optional[Sequence[str]] = None,
    ascending=True,
):
    """Write a schema-stable CSV, including when there are zero rows."""
    frame = pd.DataFrame(rows, columns=list(columns))
    if sort_by and not frame.empty:
        frame = frame.sort_values(list(sort_by), ascending=ascending)
    frame.to_csv(output_path, index=False)


def single_substitution_neighbors(sequence: str) -> Iterator[str]:
    for index, base in enumerate(sequence):
        for replacement in DNA:
            if replacement != base:
                yield sequence[:index] + replacement + sequence[index + 1 :]


def collapse_umi_payloads(payloads: Dict[str, Counter]) -> List[Counter]:
    """Apply split-pipe's abundance-directed Hamming-1 UMI collapse."""
    # split-pipe processes a query-name-sorted BAM, so equal-abundance UMIs are
    # encountered lexically within a barcode. Make that tie behavior explicit
    # for controlled STAR BAMs, which are intentionally left unsorted.
    ordered = sorted(payloads, key=lambda umi: (-payloads[umi]["total"], umi))
    neighbors: Dict[str, List[str]] = defaultdict(list)
    if len(ordered) > 50:
        umi_set = set(ordered)
        for umi in ordered:
            for mutant in single_substitution_neighbors(umi):
                if mutant in umi_set:
                    neighbors[mutant].append(umi)
    else:
        for first, second in itertools.combinations(ordered, 2):
            if hamming_distance(first, second, 1) <= 1:
                neighbors[first].append(second)
                neighbors[second].append(first)

    remaining = {umi: Counter(payloads[umi]) for umi in ordered}
    for umi in ordered:
        if umi not in remaining:
            continue
        current = remaining[umi]
        for neighbor in neighbors.get(umi, []):
            other = remaining.get(neighbor)
            if other is not None and other["total"] < current["total"]:
                current.update(other)
                del remaining[neighbor]
    return list(remaining.values())


def _open_synthetic_outputs(
    output_dir: str,
    library: str,
    enabled: bool,
    write_r2: bool,
    read_length: int,
):
    if not enabled:
        return {}, {}
    syn_dir = os.path.join(output_dir, library, f"synthetic_{read_length}")
    os.makedirs(syn_dir, exist_ok=True)
    handles = {}
    paths = {}
    for category in ("all", "clean", "rejected"):
        r1_path = os.path.join(
            syn_dir, f"{library}_{category}_R1_{read_length}.fastq.gz"
        )
        handles[(category, 1)] = gzip.open(r1_path, "wt")
        paths[f"{category}_r1"] = r1_path
        if write_r2:
            r2_path = os.path.join(syn_dir, f"{library}_{category}_R2.fastq.gz")
            handles[(category, 2)] = gzip.open(r2_path, "wt")
            paths[f"{category}_r2"] = r2_path
    return handles, paths


def _close_handles(handles: dict):
    for handle in handles.values():
        handle.close()


def benchmark_library(
    cfg: dict,
    config_path: str,
    library_name: str,
    threads: int = 1,
    *,
    execution_backend: str = "auto",
    batch_size: int = 20_000,
):
    benchmark_started = time.perf_counter()
    if threads < 1:
        raise ValueError("threads must be positive")
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    if execution_backend not in {"auto", "serial", "parallel"}:
        raise ValueError("execution_backend must be auto, serial, or parallel")
    use_parallel = execution_backend == "parallel" or (
        execution_backend == "auto" and threads > 1
    )
    qc = cfg["gex_grna_contamination_qc"]
    raw_cfg = qc["raw_benchmark"]
    short_cfg = qc["short_read_filters"]
    synthetic_cfg = qc.get("synthetic_fastq", {})
    molecule_cfg = qc.get("molecule_estimation", {})
    alignment_cfg = qc.get("alignment_benchmark", {})
    library = get_library(cfg, library_name)
    pairs = resolve_library_pairs(cfg, config_path, library)
    model = make_model(cfg, config_path)
    try:
        from .gex_grna_filter import GexGrnaFilter, PackedGuideIndex, settings_from_config
    except ImportError:
        from gex_grna_filter import GexGrnaFilter, PackedGuideIndex, settings_from_config
    base_config_path = cfg.get("_base_config_path") or config_path
    production_settings = settings_from_config(cfg, base_config_path)
    production_filter = GexGrnaFilter(
        PackedGuideIndex.from_reference_csv(
            production_settings.reference_path,
            guide_lengths=production_settings.guide_lengths,
            max_seed_bucket_size=production_settings.max_seed_bucket_size,
        ),
        production_settings,
    )
    molecule_enabled = bool(molecule_cfg.get("enabled", True))
    alignment_enabled = bool(alignment_cfg.get("enabled", True))
    barcode_corrector = None
    molecule_estimator = None
    if molecule_enabled or alignment_enabled:
        if float(raw_cfg.get("sample_fraction", 1.0)) != 1.0:
            raise ValueError(
                "molecule_estimation/alignment_benchmark requires "
                "raw_benchmark.sample_fraction=1.0"
            )
        print(f"  {library_name}: fitting split-pipe barcode correction sample")
        barcode_corrector = make_barcode_corrector(
            cfg,
            config_path,
            library_name,
            [r2_path for _, r2_path in pairs],
        )
        if molecule_enabled:
            molecule_estimator = MoleculeEstimator(
                int(molecule_cfg.get("hll_precision", 18))
            )

    output_root = resolve_path(qc["output_dir"], config_path)
    library_dir = os.path.join(output_root, library_name)
    os.makedirs(library_dir, exist_ok=True)
    read_length = int(qc.get("short_read_length", 54))
    max_reads = int(raw_cfg.get("max_reads", 0))
    sample_fraction = float(raw_cfg.get("sample_fraction", 1.0))
    read_table_max = int(raw_cfg.get("read_table_max", 200_000))
    selected_method = str(short_cfg.get("selected_method", "tso_guide_d1"))

    synthetic_enabled = bool(synthetic_cfg.get("enabled", True))
    synthetic_max = int(synthetic_cfg.get("max_reads", 250_000))
    synthetic_fraction = float(synthetic_cfg.get("sample_fraction", 1.0))
    if not 0 < synthetic_fraction <= 1:
        raise ValueError("synthetic_fastq.sample_fraction must be in (0, 1]")
    write_r2 = bool(synthetic_cfg.get("write_r2", True))
    syn_handles, syn_paths = _open_synthetic_outputs(
        output_root, library_name, synthetic_enabled, write_r2, read_length
    )

    alignment_handles = {}
    alignment_paths = {}
    alignment_manifest_handle = None
    alignment_manifest_writer = None
    alignment_max = int(alignment_cfg.get("max_truth_reads", 2_000_000))
    alignment_fraction = float(alignment_cfg.get("sample_fraction", 1.0))
    if not 0 < alignment_fraction <= 1:
        raise ValueError("alignment_benchmark.sample_fraction must be in (0, 1]")
    alignment_molecule_sketch = HyperLogLog(
        int(alignment_cfg.get("hll_precision", 18))
    )
    alignment_trim_counts = Counter()
    qc_module_digest = file_sha256(os.path.abspath(__file__))
    benchmark_config_digest = benchmark_config_sha256(cfg, library_name)
    alignment_panel_hasher = hashlib.blake2b(digest_size=16)
    alignment_panel_hasher.update(
        json.dumps(
            {
                "library": library_name,
                "read_length": read_length,
                "sample_fraction": alignment_fraction,
                "max_truth_reads": alignment_max,
                "splitpipe_tso": model.splitpipe_tso,
                "splitpipe_tso_start_limit": model.splitpipe_tso_start_limit,
                "genome_dir": cfg["paths"].get("genome_dir", ""),
                "star_clip3p_adapter": alignment_cfg.get(
                    "star_clip3p_adapter", "GAGGTGGTTGGA"
                ),
                "star_multimap_nmax": 3,
                "qc_module_sha256": qc_module_digest,
                "benchmark_config_sha256": benchmark_config_digest,
            },
            sort_keys=True,
        ).encode("utf-8")
    )
    if alignment_enabled:
        alignment_root = os.path.join(output_root, library_name, "alignment_input")
        for variant in ("full", f"short{read_length}"):
            variant_dir = os.path.join(alignment_root, variant)
            os.makedirs(variant_dir, exist_ok=True)
            suffix = "full" if variant == "full" else str(read_length)
            path = os.path.join(
                variant_dir,
                f"{library_name}_ground_truth_R1_{suffix}.fastq.gz",
            )
            alignment_handles[variant] = gzip.open(path, "wt")
            alignment_paths[variant] = path
        manifest_path = os.path.join(
            alignment_root, "ground_truth_manifest.tsv.gz"
        )
        alignment_manifest_handle = gzip.open(manifest_path, "wt", newline="")
        alignment_paths["manifest"] = manifest_path

    read_table_path = os.path.join(library_dir, "read_evidence.tsv.gz")
    read_table_handle = gzip.open(read_table_path, "wt", newline="")
    read_writer = None
    truth_counts = Counter()
    guide_counts = Counter()
    target_counts = Counter()
    tso_offsets = Counter()
    confusion_by_truth: Dict[str, Dict[str, Confusion]] = {
        "structure_truth": defaultdict(Confusion),
        "known_structure_truth": defaultdict(Confusion),
        "inclusive_truth": defaultdict(Confusion),
    }
    positive_stratum_totals = Counter()
    positive_stratum_calls = Counter()
    hard_negative_stratum_totals = Counter()
    hard_negative_stratum_calls = Counter()
    base_comp = initialize_base_composition(int(raw_cfg.get("base_composition_length", 100)))

    records_seen = records_evaluated = synthetic_written = table_written = 0
    alignment_written = 0
    alignment_candidates_after_hash = 0
    stop = False
    scan_started = time.perf_counter()
    if use_parallel:
        try:
            from . import qc_gex_grna_parallel as parallel_qc
        except ImportError:
            import qc_gex_grna_parallel as parallel_qc

        parallel_qc.configure_runtime(
            model=model,
            barcode_corrector=barcode_corrector,
            read_length=read_length,
            short_config=short_cfg,
            selected_method=selected_method,
            molecule_enabled=molecule_enabled,
            hll_precision=int(molecule_cfg.get("hll_precision", 18)),
            base_composition_length=int(
                raw_cfg.get("base_composition_length", 100)
            ),
            alignment_enabled=alignment_enabled,
            alignment_fraction=alignment_fraction,
        )

        def stable_fraction_bytes(value: bytes) -> float:
            digest = hashlib.blake2b(value, digest_size=8).digest()
            return int.from_bytes(digest, "big") / float(2**64)

        synthetic_selected_count = 0

        def parallel_batches():
            nonlocal records_seen, records_evaluated, synthetic_selected_count
            batch_id = 0
            tasks = []
            raw_records = []
            finished = False
            for r1_path, r2_path in pairs:
                for read_id_bytes, raw_r1, raw_r2 in iter_fastq_pairs_binary_raw(
                    r1_path, r2_path
                ):
                    records_seen += 1
                    if (
                        sample_fraction < 1.0
                        and stable_fraction_bytes(read_id_bytes) >= sample_fraction
                    ):
                        continue
                    records_evaluated += 1
                    include_read_row = bool(
                        read_table_max == 0 or records_evaluated <= read_table_max
                    )
                    include_synthetic = bool(
                        synthetic_enabled
                        and stable_fraction_bytes(read_id_bytes + b":synthetic")
                        < synthetic_fraction
                        and (
                            synthetic_max == 0
                            or synthetic_selected_count < synthetic_max
                        )
                    )
                    synthetic_selected_count += int(include_synthetic)
                    tasks.append(
                        (
                            raw_r1.sequence.decode("ascii"),
                            raw_r2.sequence.decode("ascii"),
                            include_read_row,
                            include_synthetic,
                        )
                    )
                    raw_records.append((read_id_bytes, raw_r1, raw_r2))
                    if len(tasks) >= batch_size:
                        yield batch_id, tasks, raw_records
                        batch_id += 1
                        tasks = []
                        raw_records = []
                    if max_reads and records_evaluated >= max_reads:
                        finished = True
                        break
                if finished:
                    break
            if tasks:
                yield batch_id, tasks, raw_records

        def consume_parallel_batch(batch_id, result, raw_records):
            del batch_id
            nonlocal read_writer, alignment_manifest_writer
            nonlocal table_written, synthetic_written, alignment_written
            nonlocal alignment_candidates_after_hash
            for local_index, row in result["read_rows"]:
                read_id_bytes, raw_r1, _ = raw_records[local_index]
                row = {"read_id": read_id_bytes.decode("ascii"), **row}
                if read_writer is None:
                    read_writer = csv.DictWriter(
                        read_table_handle,
                        fieldnames=list(row.keys()),
                        delimiter="\t",
                    )
                    read_writer.writeheader()
                read_writer.writerow(row)
                table_written += 1

            for local_index, selected_call in result["synthetic"]:
                _, raw_r1, raw_r2 = raw_records[local_index]
                r1 = raw_r1.as_text()
                r2 = raw_r2.as_text()
                write_fastq_record(syn_handles[("all", 1)], r1, read_length)
                if write_r2:
                    write_fastq_record(syn_handles[("all", 2)], r2)
                category = "rejected" if selected_call else "clean"
                write_fastq_record(syn_handles[(category, 1)], r1, read_length)
                if write_r2:
                    write_fastq_record(syn_handles[(category, 2)], r2)
                synthetic_written += 1

            for local_index, full, barcode_call in result["alignment"]:
                alignment_candidates_after_hash += 1
                if alignment_max and alignment_written >= alignment_max:
                    continue
                read_id_bytes, raw_r1, _ = raw_records[local_index]
                r1 = raw_r1.as_text()
                read_id = read_id_bytes.decode("ascii")
                production_removed = bool(
                    production_filter.route(r1.sequence).filter_positive
                )
                full_sequence, full_quality, full_trim_index = (
                    splitpipe_processed_record(
                        r1,
                        model.splitpipe_tso,
                        model.splitpipe_tso_start_limit,
                    )
                )
                short_sequence, short_quality, short_trim_index = (
                    splitpipe_processed_record(
                        r1,
                        model.splitpipe_tso,
                        model.splitpipe_tso_start_limit,
                        raw_length=read_length,
                    )
                )
                alignment_written += 1
                write_alignment_truth_record(
                    alignment_handles["full"],
                    alignment_written,
                    full_sequence,
                    full_quality,
                    barcode_call,
                    full,
                    full_trim_index,
                    production_removed,
                )
                write_alignment_truth_record(
                    alignment_handles[f"short{read_length}"],
                    alignment_written,
                    short_sequence,
                    short_quality,
                    barcode_call,
                    full,
                    short_trim_index,
                    production_removed,
                )
                alignment_molecule_sketch.add(
                    barcode_call.cell,
                    barcode_call.barcode_indices,
                    barcode_call.umi,
                    full.guide_query,
                )
                alignment_trim_counts["full_trimmed"] += int(
                    full_trim_index >= 0
                )
                alignment_trim_counts["short_trimmed"] += int(
                    short_trim_index >= 0
                )
                alignment_trim_counts["production_removed"] += int(
                    production_removed
                )
                alignment_trim_counts["production_retained"] += int(
                    not production_removed
                )
                alignment_panel_hasher.update(
                    "\x1f".join(
                        (
                            read_id,
                            barcode_call.cell,
                            barcode_call.barcode_indices,
                            barcode_call.umi,
                            full.guide_query,
                            full_sequence,
                            short_sequence,
                            str(int(production_removed)),
                        )
                    ).encode("utf-8")
                )
                alignment_panel_hasher.update(b"\x1e")
                manifest_row = {
                    "truth_id": f"GQ{alignment_written:09d}",
                    "original_read_id": read_id,
                    "cell": barcode_call.cell,
                    "barcode_indices": barcode_call.barcode_indices,
                    "umi": barcode_call.umi,
                    "truth_tier": full.tier,
                    "guide_query": full.guide_query,
                    "guide_names": ",".join(full.guide_names),
                    "target_genes": ",".join(full.target_genes),
                    "tso_offset": full.tso_start,
                    "raw_full_length": len(r1.sequence),
                    "processed_full_length": len(full_sequence),
                    "processed_short_length": len(short_sequence),
                    "full_tso_trimmed": full_trim_index >= 0,
                    "short_tso_trimmed": short_trim_index >= 0,
                    "production_filter_removed": production_removed,
                }
                if alignment_manifest_writer is None:
                    alignment_manifest_writer = csv.DictWriter(
                        alignment_manifest_handle,
                        fieldnames=list(manifest_row.keys()),
                        delimiter="\t",
                    )
                    alignment_manifest_writer.writeheader()
                alignment_manifest_writer.writerow(manifest_row)

        parallel_accumulator = parallel_qc.aggregate_artifact_batches(
            parallel_batches(),
            workers=threads,
            consume_batch=consume_parallel_batch,
        )
        if parallel_accumulator.records_evaluated != records_evaluated:
            raise RuntimeError(
                "parallel QC record count differs from the parent FASTQ stream"
            )
        truth_counts = parallel_accumulator.truth_counts
        guide_counts = parallel_accumulator.guide_counts
        target_counts = parallel_accumulator.target_counts
        tso_offsets = parallel_accumulator.tso_offsets
        confusion_by_truth = parallel_accumulator.confusion_by_truth
        positive_stratum_totals = parallel_accumulator.positive_stratum_totals
        positive_stratum_calls = parallel_accumulator.positive_stratum_calls
        hard_negative_stratum_totals = (
            parallel_accumulator.hard_negative_stratum_totals
        )
        hard_negative_stratum_calls = (
            parallel_accumulator.hard_negative_stratum_calls
        )
        base_comp = parallel_accumulator.base_composition
        molecule_estimator = parallel_accumulator.molecule_estimator

    scan_pairs = () if use_parallel else pairs
    try:
        for r1_path, r2_path in scan_pairs:
            for r1, r2 in iter_fastq_pairs(r1_path, r2_path):
                records_seen += 1
                read_id = normalize_read_id(r1.header)
                if sample_fraction < 1.0 and stable_fraction(read_id) >= sample_fraction:
                    continue
                records_evaluated += 1
                full = model.classify_full(r1.sequence)
                short = model.short_features(r1.sequence, read_length, short_cfg)
                barcode_call = (
                    barcode_corrector.call(r2.sequence)
                    if barcode_corrector is not None
                    else None
                )
                if molecule_estimator is not None:
                    molecule_estimator.add(barcode_call, r1.sequence, full)
                truth_counts[full.tier] += 1
                truth_counts["structure_truth"] += int(full.structure_truth)
                truth_counts["known_structure_truth"] += int(full.known_structure_truth)
                truth_counts["complete_reference_truth"] += int(full.complete_reference_truth)
                truth_counts["tso_known_guide_truth"] += int(full.tso_known_guide_truth)
                truth_counts["tso_found"] += int(full.tso_found)
                inclusive_truth = bool(
                    full.structure_truth or full.tso_known_guide_truth
                )
                truth_counts["inclusive_truth"] += int(inclusive_truth)
                if full.structure_truth:
                    tso_offsets[full.tso_start] += 1
                    if full.guide_names:
                        for guide_name in full.guide_names:
                            guide_counts[guide_name] += 1 / len(full.guide_names)
                    if full.target_genes:
                        for target in full.target_genes:
                            target_counts[target] += 1 / len(full.target_genes)

                update_base_composition(base_comp, "all", r1.sequence)
                if full.structure_truth:
                    update_base_composition(base_comp, "structure_truth", r1.sequence)
                elif full.tso_found:
                    update_base_composition(base_comp, "tso_positive_gex_like", r1.sequence)

                hard_negative = full.tso_found and not full.structure_truth
                for truth_name, truth in (
                    ("structure_truth", full.structure_truth),
                    ("known_structure_truth", full.known_structure_truth),
                    ("inclusive_truth", inclusive_truth),
                ):
                    for method, call in short["calls"].items():
                        confusion_by_truth[truth_name][method].update(
                            truth, bool(call), hard_negative=hard_negative and not truth
                        )

                if full.structure_truth:
                    stratum = (full.tso_start, full.guide_length, full.variant)
                    positive_stratum_totals[stratum] += 1
                    for method, call in short["calls"].items():
                        positive_stratum_calls[stratum + (method,)] += int(call)
                elif full.tso_found:
                    stratum = (full.tso_start, 0, 0)
                    hard_negative_stratum_totals[stratum] += 1
                    for method, call in short["calls"].items():
                        hard_negative_stratum_calls[stratum + (method,)] += int(call)

                if selected_method not in short["calls"]:
                    raise KeyError(
                        f"selected_method '{selected_method}' is not available. "
                        f"Choices: {sorted(short['calls'])}"
                    )
                selected_call = bool(short["calls"][selected_method])
                if full.structure_truth:
                    truth_counts["structure_selected_filter_call"] += int(selected_call)
                    truth_counts["structure_predicted_splitpipe_tso_trim"] += int(
                        short["splitpipe_tso_index"] >= 0
                    )
                    truth_counts["structure_splitpipe_trim_known_guide_d2"] += int(
                        short["calls"].get("splitpipe_trim_guide_d2", False)
                    )
                elif hard_negative:
                    truth_counts["hard_negative_selected_filter_call"] += int(
                        selected_call
                    )

                alignment_hash_selected = (
                    alignment_enabled
                    and full.structure_truth
                    and barcode_call is not None
                    and barcode_call.valid
                    and stable_fraction(
                        barcode_call.cell
                        + "|"
                        + barcode_call.barcode_indices
                        + ":alignment"
                    )
                    < alignment_fraction
                )
                alignment_candidates_after_hash += int(alignment_hash_selected)
                alignment_selected = bool(
                    alignment_hash_selected
                    and (alignment_max == 0 or alignment_written < alignment_max)
                )
                if alignment_selected:
                    production_removed = bool(
                        production_filter.route(r1.sequence).filter_positive
                    )
                    full_sequence, full_quality, full_trim_index = (
                        splitpipe_processed_record(
                            r1,
                            model.splitpipe_tso,
                            model.splitpipe_tso_start_limit,
                        )
                    )
                    short_sequence, short_quality, short_trim_index = (
                        splitpipe_processed_record(
                            r1,
                            model.splitpipe_tso,
                            model.splitpipe_tso_start_limit,
                            raw_length=read_length,
                        )
                    )
                    alignment_written += 1
                    write_alignment_truth_record(
                        alignment_handles["full"],
                        alignment_written,
                        full_sequence,
                        full_quality,
                        barcode_call,
                        full,
                        full_trim_index,
                        production_removed,
                    )
                    write_alignment_truth_record(
                        alignment_handles[f"short{read_length}"],
                        alignment_written,
                        short_sequence,
                        short_quality,
                        barcode_call,
                        full,
                        short_trim_index,
                        production_removed,
                    )
                    alignment_molecule_sketch.add(
                        barcode_call.cell,
                        barcode_call.barcode_indices,
                        barcode_call.umi,
                        full.guide_query,
                    )
                    alignment_trim_counts["full_trimmed"] += int(
                        full_trim_index >= 0
                    )
                    alignment_trim_counts["short_trimmed"] += int(
                        short_trim_index >= 0
                    )
                    alignment_trim_counts["production_removed"] += int(
                        production_removed
                    )
                    alignment_trim_counts["production_retained"] += int(
                        not production_removed
                    )
                    alignment_panel_hasher.update(
                        "\x1f".join(
                            (
                                read_id,
                                barcode_call.cell,
                                barcode_call.barcode_indices,
                                barcode_call.umi,
                                full.guide_query,
                                full_sequence,
                                short_sequence,
                                str(int(production_removed)),
                            )
                        ).encode("utf-8")
                    )
                    alignment_panel_hasher.update(b"\x1e")
                    manifest_row = {
                        "truth_id": f"GQ{alignment_written:09d}",
                        "original_read_id": read_id,
                        "cell": barcode_call.cell,
                        "barcode_indices": barcode_call.barcode_indices,
                        "umi": barcode_call.umi,
                        "truth_tier": full.tier,
                        "guide_query": full.guide_query,
                        "guide_names": ",".join(full.guide_names),
                        "target_genes": ",".join(full.target_genes),
                        "tso_offset": full.tso_start,
                        "raw_full_length": len(r1.sequence),
                        "processed_full_length": len(full_sequence),
                        "processed_short_length": len(short_sequence),
                        "full_tso_trimmed": full_trim_index >= 0,
                        "short_tso_trimmed": short_trim_index >= 0,
                        "production_filter_removed": production_removed,
                    }
                    if alignment_manifest_writer is None:
                        alignment_manifest_writer = csv.DictWriter(
                            alignment_manifest_handle,
                            fieldnames=list(manifest_row.keys()),
                            delimiter="\t",
                        )
                        alignment_manifest_writer.writeheader()
                    alignment_manifest_writer.writerow(manifest_row)

                synthetic_selected = (
                    synthetic_enabled
                    and stable_fraction(read_id + ":synthetic") < synthetic_fraction
                    and (synthetic_max == 0 or synthetic_written < synthetic_max)
                )
                if synthetic_selected:
                    write_fastq_record(syn_handles[("all", 1)], r1, read_length)
                    if write_r2 and r2 is not None:
                        write_fastq_record(syn_handles[("all", 2)], r2)
                    category = "rejected" if selected_call else "clean"
                    write_fastq_record(syn_handles[(category, 1)], r1, read_length)
                    if write_r2 and r2 is not None:
                        write_fastq_record(syn_handles[(category, 2)], r2)
                    synthetic_written += 1

                if read_table_max == 0 or table_written < read_table_max:
                    guide_hit = short["guide_hit"]
                    template_hit = short["template_prefix_hit"]
                    row = {
                        "read_id": read_id,
                        "full_length": len(r1.sequence),
                        "truth_tier": full.tier,
                        "structure_truth": full.structure_truth,
                        "known_structure_truth": full.known_structure_truth,
                        "complete_reference_truth": full.complete_reference_truth,
                        "tso_known_guide_truth": full.tso_known_guide_truth,
                        "inclusive_truth": inclusive_truth,
                        "full_tso_offset": full.tso_start,
                        "full_tso_distance": full.tso_distance,
                        "full_guide_query": full.guide_query,
                        "full_guide_names": ",".join(full.guide_names),
                        "full_target_genes": ",".join(full.target_genes),
                        "full_variant": full.variant,
                        "short_tso_distance": short["tso_distance"],
                        "short_tso_at_expected_start": short[
                            "tso_at_expected_start"
                        ],
                        "predicted_splitpipe_tso_index": short[
                            "splitpipe_tso_index"
                        ],
                        "short_guide_distance": guide_hit.distance if guide_hit else -1,
                        "short_guide_length": guide_hit.length if guide_hit else 0,
                        "short_template_prefix_length": (
                            template_hit.length if template_hit else 0
                        ),
                        "short_template_prefix_distance": (
                            template_hit.distance if template_hit else -1
                        ),
                        "short_template_boundary_bases": (
                            template_hit.minimum_boundary_bases if template_hit else 0
                        ),
                        "short_template_identity_ambiguous": (
                            template_hit.identity_ambiguous if template_hit else False
                        ),
                        "short_template_candidate_guides": (
                            ",".join(sorted(template_hit.guide_names))
                            if template_hit
                            else ""
                        ),
                        "short_partial_stem_length": short["partial_stem_length"],
                        "short_partial_stem_mismatches": (
                            short["partial_stem_mismatches"]
                            if short["partial_stem_mismatches"] is not None
                            else -1
                        ),
                        "short_template_hamming": (
                            short["template_hamming"]
                            if short["template_hamming"] is not None
                            else -1
                        ),
                    }
                    row.update({f"call_{k}": bool(v) for k, v in short["calls"].items()})
                    if read_writer is None:
                        read_writer = csv.DictWriter(
                            read_table_handle, fieldnames=list(row.keys()), delimiter="\t"
                        )
                        read_writer.writeheader()
                    read_writer.writerow(row)
                    table_written += 1

                if records_evaluated % 1_000_000 == 0:
                    print(
                        f"  {library_name}: evaluated {records_evaluated:,} reads; "
                        f"structure contaminants {truth_counts['structure_truth']:,}",
                        flush=True,
                    )
                if max_reads and records_evaluated >= max_reads:
                    stop = True
                    break
            if stop:
                break
    finally:
        read_table_handle.close()
        _close_handles(syn_handles)
        _close_handles(alignment_handles)
        if alignment_manifest_handle is not None:
            alignment_manifest_handle.close()
    scan_wall_seconds = time.perf_counter() - scan_started

    metric_rows = []
    for truth_name, methods in confusion_by_truth.items():
        for method, confusion in sorted(methods.items()):
            row = {"library": library_name, "truth_definition": truth_name, "method": method}
            row.update(confusion.as_dict())
            metric_rows.append(row)
    metric_columns = [
        "library",
        "truth_definition",
        "method",
        "tp",
        "fp",
        "tn",
        "fn",
        "sensitivity",
        "sensitivity_95_low",
        "sensitivity_95_high",
        "specificity",
        "false_positive_rate",
        "false_positive_rate_95_low",
        "false_positive_rate_95_high",
        "precision",
        "f1",
        "hard_negative_n",
        "hard_negative_fp",
        "hard_negative_fpr",
        "hard_negative_fpr_95_low",
        "hard_negative_fpr_95_high",
    ]
    metrics_df = pd.DataFrame(metric_rows, columns=metric_columns)
    metrics_path = os.path.join(
        library_dir, f"short{read_length}_filter_metrics.csv"
    )
    metrics_df.to_csv(metrics_path, index=False)

    stratum_rows = []
    for group, totals, calls in (
        ("structure_positive", positive_stratum_totals, positive_stratum_calls),
        (
            "tso_positive_structure_negative",
            hard_negative_stratum_totals,
            hard_negative_stratum_calls,
        ),
    ):
        for (offset, guide_length, variant), total in sorted(totals.items()):
            methods = sorted(
                key[-1]
                for key in calls
                if key[:3] == (offset, guide_length, variant)
            )
            for method in methods:
                called = int(calls[(offset, guide_length, variant, method)])
                stratum_rows.append(
                    {
                        "library": library_name,
                        "group": group,
                        "tso_offset": offset,
                        "guide_length": guide_length,
                        "variant": variant,
                        "method": method,
                        "n": int(total),
                        "called": called,
                        "rate": called / total if total else float("nan"),
                    }
                )
    stratum_path = os.path.join(
        library_dir, f"short{read_length}_metrics_by_tso_offset.csv"
    )
    pd.DataFrame(
        stratum_rows,
        columns=[
            "library",
            "group",
            "tso_offset",
            "guide_length",
            "variant",
            "method",
            "n",
            "called",
            "rate",
        ],
    ).to_csv(stratum_path, index=False)

    primary = metrics_df[
        metrics_df["truth_definition"] == "structure_truth"
    ].copy()
    max_hard_fpr = float(short_cfg.get("recommendation_max_hard_negative_fpr", 0.001))
    hard_fpr_gate = (
        "hard_negative_fpr_95_high"
        if short_cfg.get("recommendation_use_upper_ci", True)
        else "hard_negative_fpr"
    )
    eligible = primary[
        (primary[hard_fpr_gate].fillna(1.0) <= max_hard_fpr)
        & ~primary["method"].str.startswith("tso_only")
        & ~primary["method"].str.startswith("tso_anywhere")
        & (primary["method"] != "tso_exact30")
    ]
    recommended_method = ""
    if not eligible.empty:
        eligible = eligible.sort_values(
            ["sensitivity", "precision", "false_positive_rate"],
            ascending=[False, False, True],
        )
        recommended_method = str(eligible.iloc[0]["method"])

    contamination_count = int(truth_counts["structure_truth"])
    ci_low, ci_high = wilson_interval(contamination_count, records_evaluated)
    tso_count = int(truth_counts["tso_found"])
    tso_contam_low, tso_contam_high = wilson_interval(
        contamination_count, tso_count
    )
    molecule_summary = (
        molecule_estimator.summary() if molecule_estimator is not None else None
    )
    alignment_cap_reached = bool(
        alignment_max and alignment_candidates_after_hash > alignment_written
    )
    confidence_tiers = {}
    for tier_name, count_key, interpretation in (
        (
            "definitive_reference_construct",
            "complete_reference_truth",
            "TSO + structure + known spacer + compatible reference iBAR",
        ),
        (
            "reference_supported_structure",
            "known_structure_truth",
            "TSO + complete structure + known spacer",
        ),
        (
            "primary_structural_surrogate",
            "structure_truth",
            "TSO + spacer-sized segment + paired stem/iBAR-sized block/tracr",
        ),
        (
            "inclusive_possible_contaminant",
            "inclusive_truth",
            "Primary structural call or TSO followed by a known spacer",
        ),
    ):
        tier_count = int(truth_counts[count_key])
        tier_low, tier_high = wilson_interval(tier_count, records_evaluated)
        confidence_tiers[tier_name] = {
            "count": tier_count,
            "fraction": tier_count / records_evaluated if records_evaluated else 0.0,
            "wilson_95_low": tier_low,
            "wilson_95_high": tier_high,
            "interpretation": interpretation,
        }
    summary = {
        "library": library_name,
        "read_length_simulated": read_length,
        "records_seen": records_seen,
        "records_evaluated": records_evaluated,
        "raw_scan_complete": max_reads == 0,
        "raw_max_reads": max_reads,
        "sample_fraction": sample_fraction,
        "input_fastq_pairs": [{"r1": a, "r2": b} for a, b in pairs],
        "input_fastq_metadata": [
            {"r1": file_stat_identity(r1), "r2": file_stat_identity(r2)}
            for r1, r2 in pairs
        ],
        "raw_read_prevalence": {
            "all_reads": records_evaluated,
            "tso_reads": tso_count,
            "structure_contaminant_reads": contamination_count,
            "contaminant_fraction_of_all_reads": (
                contamination_count / records_evaluated if records_evaluated else 0.0
            ),
            "contaminant_fraction_of_all_reads_wilson_95": [ci_low, ci_high],
            "contaminant_fraction_of_tso_reads": (
                contamination_count / tso_count if tso_count else 0.0
            ),
            "contaminant_fraction_of_tso_reads_wilson_95": [
                tso_contam_low,
                tso_contam_high,
            ],
        },
        "prealignment_molecule_prevalence": molecule_summary,
        "alignment_ground_truth_input": {
            "enabled": alignment_enabled,
            "reads_written_per_variant": alignment_written,
            "truth_panel_token": alignment_panel_hasher.hexdigest(),
            "qc_module_sha256": qc_module_digest,
            "benchmark_config_sha256": benchmark_config_digest,
            "truth_candidates_after_hash_sampling": alignment_candidates_after_hash,
            "safety_cap_reached": alignment_cap_reached,
            "estimated_distinct_cell_pb_guide_umi_molecules": (
                alignment_molecule_sketch.estimate() if alignment_enabled else 0.0
            ),
            "sample_fraction": alignment_fraction,
            "sampling_unit": (
                "corrected (CB, pB) group; all reads/UMIs in a selected group are "
                "retained so PCR families and Hamming-1 UMI collapse remain auditable"
            ),
            "max_truth_reads": alignment_max,
            "predicted_tso_trim_counts": {
                key: int(value) for key, value in alignment_trim_counts.items()
            },
            "paths": alignment_paths,
            "note": (
                "These FASTQs are already transformed to the exact sequence passed "
                "to STAR: full raw pilot or raw R1 truncated to the requested short "
                "length, followed by split-pipe's exact 17-mer TSO trim. QNAME keeps "
                "truth, cell, and UMI fields for direct outcome joins."
            ),
        },
        "full_length_surrogate_truth": {
            "structure_contaminants": contamination_count,
            "fraction": contamination_count / records_evaluated if records_evaluated else 0.0,
            "wilson_95_low": ci_low,
            "wilson_95_high": ci_high,
            "known_structure_contaminants": int(truth_counts["known_structure_truth"]),
            "complete_reference_contaminants": int(truth_counts["complete_reference_truth"]),
            "tso_known_guide_contaminants": int(
                truth_counts["tso_known_guide_truth"]
            ),
            "inclusive_contaminants": int(truth_counts["inclusive_truth"]),
            "tso_positive_reads": int(truth_counts["tso_found"]),
            "confidence_tiers": confidence_tiers,
            "tier_counts": {
                key: int(value)
                for key, value in truth_counts.items()
                if key
                not in {
                    "structure_truth",
                    "known_structure_truth",
                    "complete_reference_truth",
                    "tso_known_guide_truth",
                    "inclusive_truth",
                    "tso_found",
                    "structure_selected_filter_call",
                    "structure_predicted_splitpipe_tso_trim",
                    "structure_splitpipe_trim_known_guide_d2",
                    "hard_negative_selected_filter_call",
                }
            },
        },
        "short_read_hazard": {
            "structure_truth_with_exact_splitpipe_tso_trim": int(
                truth_counts["structure_predicted_splitpipe_tso_trim"]
            ),
            "structure_truth_with_splitpipe_trim_and_known_guide_d2": int(
                truth_counts["structure_splitpipe_trim_known_guide_d2"]
            ),
            "structure_truth_called_by_selected_filter": int(
                truth_counts["structure_selected_filter_call"]
            ),
            "tso_positive_structure_negative_called_by_selected_filter": int(
                truth_counts["hard_negative_selected_filter_call"]
            ),
        },
        "selected_filter": selected_method,
        "execution": {
            "backend": "bounded_multiprocessing" if use_parallel else "serial",
            "workers": threads if use_parallel else 1,
            "batch_size": batch_size,
            "maximum_pending_batches": 2 * threads if use_parallel else 1,
            "maximum_pending_records": (
                2 * threads * batch_size if use_parallel else batch_size
            ),
            "scan_wall_seconds": scan_wall_seconds,
            "scan_read_pairs_per_second": (
                records_evaluated / scan_wall_seconds
                if scan_wall_seconds
                else 0.0
            ),
            "wall_seconds_before_summary_write": (
                time.perf_counter() - benchmark_started
            ),
        },
        "recommended_filter_from_benchmark": recommended_method,
        "recommendation_max_hard_negative_fpr": max_hard_fpr,
        "recommendation_fpr_gate": hard_fpr_gate,
        "expected_tso_start_max": int(short_cfg.get("max_tso_start", 2)),
        "synthetic_reads_written": synthetic_written,
        "synthetic_sample_fraction": synthetic_fraction,
        "synthetic_fastqs": syn_paths,
        "outputs": {
            "read_evidence": read_table_path,
            "filter_metrics": metrics_path,
            "filter_metrics_by_tso_offset": stratum_path,
        },
        "edlib_available": HAS_EDLIB,
    }
    if barcode_corrector is not None:
        summary["barcode_correction"] = {
            "perfect_combination_sample_count": int(
                sum(barcode_corrector.perfect_combination_counts.values())
            ),
            "perfect_combination_count_threshold": int(
                barcode_corrector.perfect_count_threshold
            ),
            "maximum_barcode_mismatches": barcode_corrector.max_mismatches,
        }
    summary_path = os.path.join(library_dir, "benchmark_summary.yaml")
    with open(summary_path, "w") as handle:
        yaml.safe_dump(summary, handle, sort_keys=False)

    write_rows_csv(
        [{"guide_name": key, "count": value} for key, value in guide_counts.items()],
        ["guide_name", "count"],
        os.path.join(library_dir, "contaminant_guide_counts.csv"),
        sort_by=["count", "guide_name"],
        ascending=[False, True],
    )
    write_rows_csv(
        [{"target_gene": key, "count": value} for key, value in target_counts.items()],
        ["target_gene", "count"],
        os.path.join(library_dir, "contaminant_target_counts.csv"),
        sort_by=["count", "target_gene"],
        ascending=[False, True],
    )
    write_rows_csv(
        [{"tso_offset": key, "count": value} for key, value in tso_offsets.items()],
        ["tso_offset", "count"],
        os.path.join(library_dir, "contaminant_tso_offsets.csv"),
        sort_by=["tso_offset"],
    )
    write_base_composition(
        base_comp, os.path.join(library_dir, "per_cycle_base_composition.csv")
    )
    if molecule_estimator is not None:
        cell_rows = []
        for cell, counts in molecule_estimator.cell_read_counts.items():
            all_reads = int(counts["all"])
            tso_reads = int(counts["tso"])
            contaminants = int(counts["structure_contaminant"])
            cell_rows.append(
                {
                    "cell": cell,
                    "valid_barcode_reads": all_reads,
                    "tso_reads": tso_reads,
                    "structure_contaminant_reads": contaminants,
                    "contaminant_fraction_of_all": (
                        contaminants / all_reads if all_reads else 0.0
                    ),
                    "contaminant_fraction_of_tso": (
                        contaminants / tso_reads if tso_reads else 0.0
                    ),
                }
            )
        write_rows_csv(
            cell_rows,
            [
                "cell",
                "valid_barcode_reads",
                "tso_reads",
                "structure_contaminant_reads",
                "contaminant_fraction_of_all",
                "contaminant_fraction_of_tso",
            ],
            os.path.join(library_dir, "cell_read_contamination.csv"),
            sort_by=[
                "structure_contaminant_reads",
                "valid_barcode_reads",
                "cell",
            ],
            ascending=[False, False, True],
        )
    print(f"  Wrote benchmark summary: {summary_path}")
    if alignment_cap_reached and alignment_cfg.get(
        "fail_if_safety_cap_reached", True
    ):
        raise RuntimeError(
            "alignment_benchmark.max_truth_reads was reached after corrected "
            "(CB,pB)-group sampling. The partial panel would truncate PCR/UMI "
            "families; lower alignment_benchmark.sample_fraction and rerun. "
            f"Raw prevalence outputs remain available in {library_dir}."
        )
    return summary


def audit_bam(cfg: dict, config_path: str, library_name: str, threads: int = 1):
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("pysam is required for BAM audit mode") from exc

    qc = cfg["gex_grna_contamination_qc"]
    bam_cfg = qc.get("bam_audit", {})
    output_root = resolve_path(qc["output_dir"], config_path)
    library_dir = os.path.join(output_root, library_name)
    os.makedirs(library_dir, exist_ok=True)
    model = make_model(cfg, config_path)
    bam_path = bam_cfg.get("bam_path_template", "").format(
        output_dir=cfg["paths"]["output_dir"], library=library_name
    )
    if not bam_path:
        bam_path = os.path.join(
            cfg["paths"]["output_dir"],
            f"gex_{library_name}",
            "process",
            "barcode_headAligned_anno.bam",
        )
    bam_path = resolve_path(bam_path, config_path)
    if not os.path.exists(bam_path):
        raise FileNotFoundError(f"Annotated BAM not found: {bam_path}")

    max_reads = int(bam_cfg.get("max_primary_reads", 0))
    stats = Counter()
    gene_counts = Counter()
    target_gene_counts = Counter()
    historical_gene_molecule_counts = Counter()
    gene_names = {}
    contaminant_mapq_counts = Counter()
    track_exact_umi_groups = bool(
        bam_cfg.get("track_exact_contaminant_umis", True)
    )
    current_cell = None
    cell_umi_payloads: Dict[Tuple[str, str, str], Counter] = {}
    cell_molecule_rows = []

    def flush_cell_molecules():
        nonlocal cell_umi_payloads, current_cell
        if current_cell is None:
            return
        by_gene_pb: Dict[Tuple[str, str], Dict[str, Counter]] = defaultdict(dict)
        for (pb, gene, umi), payload in cell_umi_payloads.items():
            by_gene_pb[(pb, gene)][umi] = payload
            if track_exact_umi_groups:
                stats["exact_cell_pb_gene_umi_groups"] += int(
                    payload["guide"] > 0
                )
        cell_counts = Counter()
        for (_pb, gene), payload_by_umi in by_gene_pb.items():
            for cluster in collapse_umi_payloads(payload_by_umi):
                cell_counts["all"] += 1
                historical_gene_molecule_counts[gene] += 1
                cell_counts["contaminated_any"] += int(cluster["guide"] > 0)
                cell_counts["contaminated_all"] += int(
                    cluster["guide"] == cluster["total"]
                )
                cell_counts["contaminated_majority"] += int(
                    2 * cluster["guide"] > cluster["total"]
                )
                cell_counts["guide_target_any"] += int(cluster["guide_target"] > 0)
        for key, value in cell_counts.items():
            stats[f"counted_molecules_{key}"] += value
        cell_molecule_rows.append(
            {
                "cell": current_cell,
                "counted_molecules": int(cell_counts["all"]),
                "contaminated_molecules_any_guide_read": int(
                    cell_counts["contaminated_any"]
                ),
                "contaminated_molecules_all_guide_reads": int(
                    cell_counts["contaminated_all"]
                ),
                "contaminated_molecules_majority_guide_reads": int(
                    cell_counts["contaminated_majority"]
                ),
                "guide_target_molecules_any": int(cell_counts["guide_target_any"]),
            }
        )
        cell_umi_payloads = {}

    with pysam.AlignmentFile(bam_path, "rb", threads=max(1, threads)) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            stats["primary_reads"] += 1
            stats["mapped_primary_reads"] += int(not read.is_unmapped)
            stats["unmapped_primary_reads"] += int(read.is_unmapped)
            sequence = read.query_sequence or ""
            if read.is_reverse:
                sequence = reverse_complement(sequence)
            evidence = model.classify_post_splitpipe(sequence)
            is_contaminant = evidence.structure_truth
            gx = read.get_tag("GX") if read.has_tag("GX") else ""
            gn = (
                normalize_gene_name(read.get_tag("GN"))
                if read.has_tag("GN")
                else ""
            )
            if gx:
                gene_names[gx] = gn
            assigned = bool(gx)
            expected_targets = set(evidence.target_genes)
            expected_match = bool(is_contaminant and gn and gn in expected_targets)
            unique = (not read.is_unmapped) and read.mapping_quality == 255
            stats["all_unique_reads"] += int(unique)
            stats["all_gene_assigned_reads"] += int(assigned)

            parts = read.query_name.split("__")
            cell = (
                read.get_tag("CB")
                if read.has_tag("CB")
                else (parts[0] if parts else "")
            )
            pb = (
                read.get_tag("pB")
                if read.has_tag("pB")
                else (parts[2] if len(parts) > 2 else "")
            )
            umi = (
                read.get_tag("pN")
                if read.has_tag("pN")
                else (parts[4] if len(parts) > 4 else "")
            )
            if assigned:
                if current_cell is None:
                    current_cell = cell
                elif cell != current_cell:
                    flush_cell_molecules()
                    current_cell = cell
                key = (pb, gx, umi)
                payload = cell_umi_payloads.setdefault(key, Counter())
                payload["total"] += 1
                payload["guide"] += int(is_contaminant)
                payload["guide_target"] += int(expected_match)

            if is_contaminant:
                stats["structure_contaminants"] += 1
                stats["mapped_structure_contaminants"] += int(not read.is_unmapped)
                stats["unmapped_structure_contaminants"] += int(read.is_unmapped)
                stats["known_structure_contaminants"] += int(
                    evidence.known_structure_truth
                )
                contaminant_mapq_counts[int(read.mapping_quality)] += 1
                stats["unique_contaminants"] += int(unique)
                stats["gene_assigned_contaminants"] += int(assigned)
                stats["assigned_to_guide_target"] += int(expected_match)
            if assigned:
                if is_contaminant:
                    gene_counts[(gx, gn, expected_match)] += 1
                    for target in evidence.target_genes or ("UNKNOWN",):
                        target_gene_counts[(target, gn, expected_match)] += 1
            if stats["primary_reads"] % 1_000_000 == 0:
                print(
                    f"  {library_name} BAM: {stats['primary_reads']:,} primary; "
                    f"{stats['structure_contaminants']:,} contaminants",
                    flush=True,
                )
            if max_reads and stats["primary_reads"] >= max_reads:
                break

    flush_cell_molecules()

    if not track_exact_umi_groups:
        stats["exact_cell_pb_gene_umi_groups"] = -1
    contaminants = stats["structure_contaminants"]
    counted_molecules = stats["counted_molecules_all"]
    any_contaminant_molecules = stats["counted_molecules_contaminated_any"]
    all_contaminant_molecules = stats["counted_molecules_contaminated_all"]
    majority_contaminant_molecules = stats[
        "counted_molecules_contaminated_majority"
    ]
    pre_align_path = os.path.join(
        os.path.dirname(bam_path), "pre_align_stats.csv"
    )
    vendor_tscp = -1
    if os.path.exists(pre_align_path):
        pre_align = pd.read_csv(pre_align_path)
        matches = pre_align.loc[
            pre_align["statistic"] == "number_of_tscp", "value"
        ]
        if not matches.empty:
            vendor_tscp = int(float(matches.iloc[0]))
    summary = {
        "library": library_name,
        "bam": bam_path,
        "qc_module_sha256": file_sha256(os.path.abspath(__file__)),
        "qc_config_sha256": benchmark_config_sha256(cfg, library_name),
        "bam_size_bytes": os.path.getsize(bam_path),
        "bam_mtime_ns": os.stat(bam_path).st_mtime_ns,
        "max_primary_reads": max_reads,
        "bam_scan_complete": max_reads == 0,
        "primary_reads": int(stats["primary_reads"]),
        "mapped_primary_reads": int(stats["mapped_primary_reads"]),
        "unmapped_primary_reads": int(stats["unmapped_primary_reads"]),
        "structure_contaminants": int(contaminants),
        "fraction_of_primary_records": (
            contaminants / stats["primary_reads"] if stats["primary_reads"] else 0.0
        ),
        "mapped_structure_contaminants": int(stats["mapped_structure_contaminants"]),
        "unmapped_structure_contaminants": int(
            stats["unmapped_structure_contaminants"]
        ),
        "fraction_of_mapped_primary_reads": (
            stats["mapped_structure_contaminants"] / stats["mapped_primary_reads"]
            if stats["mapped_primary_reads"]
            else 0.0
        ),
        "known_structure_contaminants": int(stats["known_structure_contaminants"]),
        "unique_contaminants": int(stats["unique_contaminants"]),
        "gene_assigned_contaminants": int(stats["gene_assigned_contaminants"]),
        "assigned_to_guide_target": int(stats["assigned_to_guide_target"]),
        "gene_assigned_contaminant_fraction_of_all_gene_assigned_reads": (
            stats["gene_assigned_contaminants"] / stats["all_gene_assigned_reads"]
            if stats["all_gene_assigned_reads"]
            else 0.0
        ),
        "counted_molecule_impact": {
            "reconstructed_all_molecules": int(counted_molecules),
            "vendor_number_of_tscp": int(vendor_tscp),
            "vendor_reconstruction_validation_applicable": bool(
                vendor_tscp >= 0 and max_reads == 0
            ),
            "reconstruction_difference": (
                int(counted_molecules - vendor_tscp)
                if vendor_tscp >= 0 and max_reads == 0
                else None
            ),
            "molecules_with_any_guide_read_upper_bound": int(
                any_contaminant_molecules
            ),
            "molecules_with_all_guide_reads_lower_bound": int(
                all_contaminant_molecules
            ),
            "molecules_with_majority_guide_reads": int(
                majority_contaminant_molecules
            ),
            "guide_target_molecules_with_any_guide_read": int(
                stats["counted_molecules_guide_target_any"]
            ),
            "fraction_with_any_guide_read": (
                any_contaminant_molecules / counted_molecules
                if counted_molecules
                else 0.0
            ),
            "fraction_with_all_guide_reads": (
                all_contaminant_molecules / counted_molecules
                if counted_molecules
                else 0.0
            ),
            "fraction_with_majority_guide_reads": (
                majority_contaminant_molecules / counted_molecules
                if counted_molecules
                else 0.0
            ),
        },
        "exact_cell_pb_gene_umi_groups_upper_bound": int(
            stats["exact_cell_pb_gene_umi_groups"]
        ),
        "notes": [
            "split-pipe v1.6.4 defaults to align_star_keep_no_map=false, so this "
            "annotated BAM normally excludes raw unmapped input reads; BAM prevalence "
            "is therefore not a raw-library contamination estimate.",
            "Only MAPQ 255 alignments receive gene assignments in molinfo; STAR is "
            "configured to emit up to three mappings.",
            "Counted molecule groups include corrected pB and use split-pipe's "
            "abundance-directed Hamming-distance-1 collapse. Any-guide and all-guide "
            "clusters bracket collisions between GEX and contaminant reads.",
            "A tiny reconstruction difference can be expected from split-pipe "
            "v1.6.4's molinfo chunk-EOF handling: output_one_bc_files slices "
            "lists as [:-1] even for the final chunk/cell, potentially omitting "
            "the last assigned read of each processing chunk. The audit counts "
            "every annotated-BAM record and reports the difference explicitly.",
        ],
    }
    summary_path = os.path.join(library_dir, "bam_audit_summary.yaml")
    with open(summary_path, "w") as handle:
        yaml.safe_dump(summary, handle, sort_keys=False)

    write_rows_csv(
        [
            {
                "gene_id": gx,
                "gene_name": gn,
                "is_guide_target": expected,
                "read_count": count,
            }
            for (gx, gn, expected), count in gene_counts.items()
        ],
        ["gene_id", "gene_name", "is_guide_target", "read_count"],
        os.path.join(library_dir, "bam_contaminant_gene_assignments.csv"),
        sort_by=["read_count"],
        ascending=False,
    )
    write_rows_csv(
        [
            {
                "guide_target": target,
                "aligned_gene": aligned,
                "is_guide_target": expected,
                "read_count": count,
            }
            for (target, aligned, expected), count in target_gene_counts.items()
        ],
        ["guide_target", "aligned_gene", "is_guide_target", "read_count"],
        os.path.join(library_dir, "bam_guide_target_concordance.csv"),
        sort_by=["read_count"],
        ascending=False,
    )
    write_rows_csv(
        [
            {
                "gene_id": gene_id,
                "gene_name": gene_names.get(gene_id, ""),
                "historical_gene_molecules": count,
            }
            for gene_id, count in historical_gene_molecule_counts.items()
        ],
        ["gene_id", "gene_name", "historical_gene_molecules"],
        os.path.join(library_dir, "bam_historical_gene_molecule_totals.csv"),
        sort_by=["historical_gene_molecules"],
        ascending=False,
    )
    write_rows_csv(
        [
            {"mapq": mapq, "read_count": count}
            for mapq, count in contaminant_mapq_counts.items()
        ],
        ["mapq", "read_count"],
        os.path.join(library_dir, "bam_contaminant_mapq_counts.csv"),
        sort_by=["mapq"],
    )
    write_rows_csv(
        cell_molecule_rows,
        [
            "cell",
            "counted_molecules",
            "contaminated_molecules_any_guide_read",
            "contaminated_molecules_all_guide_reads",
            "contaminated_molecules_majority_guide_reads",
            "guide_target_molecules_any",
        ],
        os.path.join(library_dir, "bam_cell_molecule_contamination.csv"),
        sort_by=["contaminated_molecules_any_guide_read", "counted_molecules"],
        ascending=False,
    )
    print(f"  Wrote BAM audit summary: {summary_path}")
    return summary


def parse_alignment_truth_qname(query_name: str) -> dict:
    fields = query_name.split("|")
    parsed = {"truth_id": fields[0]}
    for field in fields[1:]:
        if "=" in field:
            key, value = field.split("=", 1)
            parsed[key] = value
    parsed["targets"] = tuple(
        target
        for target in parsed.get("T", "UNKNOWN").split(",")
        if target and target != "UNKNOWN"
    )
    parsed["tso_trimmed"] = parsed.get("X", "0") == "1"
    parsed["production_removed"] = parsed.get("R", "0") == "1"
    return parsed


def merge_intervals(intervals: Dict[int, int]) -> List[Tuple[int, int]]:
    merged = []
    for start, end in sorted((int(start), int(end)) for start, end in intervals.items()):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(int(start), int(end)) for start, end in merged]


def normalize_gene_name(value) -> str:
    """Remove GTF serialization quotes without changing the gene symbol."""
    name = str(value).strip()
    while (
        len(name) >= 2
        and name[0] == name[-1]
        and name[0] in {"'", '"'}
    ):
        name = name[1:-1].strip()
    return name


class GeneAnnotator:
    """Self-contained implementation of split-pipe molinfo gene assignment."""

    def __init__(self, gene_info_path: str):
        with open(gene_info_path) as handle:
            info = json.load(handle)
        self.gene_bins = info["gene_bins"]
        self.bin_size = int(info["gtf_dict_stepsize"])
        self.starts = {key: int(value) for key, value in info["gene_starts"].items()}
        self.ends = {key: int(value) for key, value in info["gene_ends"].items()}
        self.names = {
            key: normalize_gene_name(value)
            for key, value in info["gene_id_to_name"].items()
        }
        self.biotypes = {
            key: str(value) for key, value in info["gene_id_to_biotype"].items()
        }
        self.exons = {
            gene: merge_intervals(intervals)
            for gene, intervals in info["genes_to_exons"].items()
        }

    def exonic_fraction(self, positions: Sequence[int], gene_id: str) -> float:
        if not positions or gene_id not in self.exons:
            return 0.0
        hits = 0
        for position in positions:
            if any(start <= position <= end for start, end in self.exons[gene_id]):
                hits += 1
        return hits / len(positions)

    def assign(self, read) -> Tuple[str, str, float]:
        positions = read.get_reference_positions(full_length=False)
        if read.is_unmapped or not positions:
            return "", "", 0.0
        strand = "-" if read.is_reverse else "+"
        bin_start = int(math.floor(positions[0] / self.bin_size) * self.bin_size)
        key = f"{read.reference_name}:{bin_start}:{strand}"
        candidates = self.gene_bins.get(key, [])
        matching = [
            gene
            for gene in candidates
            if self.starts[gene] <= positions[0] and positions[-1] <= self.ends[gene]
        ]
        if not matching:
            return "", "", 0.0
        if len(matching) == 1:
            gene = matching[0]
            return gene, self.names.get(gene, ""), self.exonic_fraction(positions, gene)

        best = None
        for gene in matching:
            name = self.names.get(gene, "")
            exonic_fraction = self.exonic_fraction(positions, gene)
            score = 4 * (self.biotypes.get(gene) == "protein_coding")
            score += 2 * ("." not in name and "-" not in name)
            score += int(exonic_fraction > 0.5)
            candidate = (score, gene, name, exonic_fraction)
            if best is None or candidate[0] > best[0]:
                best = candidate
        return best[1], best[2], best[3]


def parse_star_final_log(path: str) -> dict:
    values = {}
    if not os.path.exists(path):
        return values
    with open(path) as handle:
        for line in handle:
            if "|" not in line:
                continue
            key, value = [part.strip() for part in line.split("|", 1)]
            normalized = (
                key.lower()
                .replace("%", "pct")
                .replace(":", "")
                .replace(",", "")
                .replace(" ", "_")
            )
            values[normalized] = value
    return values


def audit_alignment_benchmark(
    cfg: dict,
    config_path: str,
    library_name: str,
    variant: str,
    threads: int = 1,
):
    del threads
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("pysam is required for alignment audit mode") from exc

    qc = cfg["gex_grna_contamination_qc"]
    alignment_cfg = qc["alignment_benchmark"]
    read_length = int(qc.get("short_read_length", 54))
    valid_variants = {"full", f"short{read_length}"}
    if variant not in valid_variants:
        raise ValueError(f"Alignment variant must be one of {sorted(valid_variants)}")
    output_root = resolve_path(qc["output_dir"], config_path)
    library_dir = os.path.join(output_root, library_name)
    input_suffix = "full" if variant == "full" else str(read_length)
    input_fastq = os.path.join(
        library_dir,
        "alignment_input",
        variant,
        f"{library_name}_ground_truth_R1_{input_suffix}.fastq.gz",
    )
    star_dir = os.path.join(
        library_dir,
        alignment_cfg.get("output_subdir", "alignment_star"),
        variant,
    )
    bam_path = os.path.join(star_dir, f"{library_name}_Aligned.out.bam")
    star_log_path = os.path.join(star_dir, f"{library_name}_Log.final.out")
    benchmark_summary_path = os.path.join(library_dir, "benchmark_summary.yaml")
    truth_panel_token = ""
    expected_truth_reads = None
    benchmark_module_digest = ""
    benchmark_config_digest = ""
    if os.path.exists(benchmark_summary_path):
        with open(benchmark_summary_path) as handle:
            benchmark_summary = yaml.safe_load(handle) or {}
        current_module_digest = file_sha256(os.path.abspath(__file__))
        benchmark_current, benchmark_reason = benchmark_summary_status(
            benchmark_summary_path,
            cfg,
            library_name,
            current_module_digest,
        )
        if not benchmark_current:
            raise RuntimeError(
                f"Refusing alignment audit because benchmark is stale: "
                f"{benchmark_reason}."
            )
        truth_input = benchmark_summary.get("alignment_ground_truth_input", {})
        truth_panel_token = str(truth_input.get("truth_panel_token", ""))
        expected_truth_reads = int(
            truth_input.get("reads_written_per_variant", 0)
        )
        benchmark_module_digest = str(truth_input.get("qc_module_sha256", ""))
        benchmark_config_digest = str(
            truth_input.get("benchmark_config_sha256", "")
        )
        if truth_input.get("safety_cap_reached", False):
            raise RuntimeError(
                "Refusing alignment audit because the truth-input safety cap "
                "truncated selected corrected CB/pB groups."
            )
        if benchmark_module_digest and benchmark_module_digest != current_module_digest:
            raise RuntimeError(
                "Refusing alignment audit because the QC module changed after "
                "the truth panel was generated; rerun the raw benchmark."
            )
        current_config_digest = benchmark_config_sha256(cfg, library_name)
        if benchmark_config_digest and benchmark_config_digest != current_config_digest:
            raise RuntimeError(
                "Refusing alignment audit because QC configuration changed after "
                "the truth panel was generated; rerun the raw benchmark."
            )
    if not os.path.exists(input_fastq):
        raise FileNotFoundError(f"Ground-truth STAR input not found: {input_fastq}")
    if not os.path.exists(bam_path):
        raise FileNotFoundError(f"Ground-truth STAR BAM not found: {bam_path}")

    gene_info_path = os.path.join(
        resolve_path(cfg["paths"]["genome_dir"], config_path), "gene_info.json"
    )
    annotator = GeneAnnotator(gene_info_path)
    precision = int(alignment_cfg.get("hll_precision", 18))
    sketches = {
        key: HyperLogLog(precision)
        for key in (
            "input",
            "mapped_primary",
            "unique",
            "gene_assigned",
            "guide_target_assigned",
        )
    }
    exact_collapse = bool(alignment_cfg.get("exact_molecule_collapse", True))
    input_umi_payloads: Dict[Tuple[str, str, str], Dict[str, Counter]] = (
        defaultdict(dict) if exact_collapse else {}
    )
    assigned_umi_payloads: Dict[Tuple[str, str, str], Dict[str, Counter]] = (
        defaultdict(dict) if exact_collapse else {}
    )
    cluster_counts: Dict[Tuple[str, str], Counter] = defaultdict(Counter)
    stats = Counter()
    input_length_counts = Counter()
    trim_input_counts = Counter()
    for record in iter_fastq_records(input_fastq):
        truth = parse_alignment_truth_qname(normalize_read_id(record.header))
        cluster_key = (truth.get("C", ""), truth.get("B", ""))
        stats["input_reads"] += 1
        cluster_counts[cluster_key]["input_reads"] += 1
        input_length_counts[len(record.sequence)] += 1
        trim_input_counts[truth["tso_trimmed"]] += 1
        stats[
            "input_production_removed"
            if truth["production_removed"]
            else "input_production_retained"
        ] += 1
        sketches["input"].add(
            truth.get("C", ""),
            truth.get("B", ""),
            truth.get("Q", ""),
            truth.get("U", ""),
        )
        if exact_collapse:
            group = (
                truth.get("C", ""),
                truth.get("B", ""),
                truth.get("Q", ""),
            )
            umi = truth.get("U", "")
            payload = input_umi_payloads[group].setdefault(umi, Counter())
            payload["total"] += 1
            payload["tso_trimmed"] += int(truth["tso_trimmed"])
    if expected_truth_reads is not None and stats["input_reads"] != expected_truth_reads:
        raise RuntimeError(
            f"Ground-truth FASTQ has {stats['input_reads']} records but current "
            f"benchmark summary expects {expected_truth_reads}; refuse stale or "
            "partial alignment input."
        )

    gene_counts = Counter()
    trim_outcomes: Dict[bool, Counter] = defaultdict(Counter)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            truth = parse_alignment_truth_qname(read.query_name)
            cluster_key = (truth.get("C", ""), truth.get("B", ""))
            trimmed = truth["tso_trimmed"]
            input_molecule = (
                truth.get("C", ""),
                truth.get("B", ""),
                truth.get("Q", ""),
                truth.get("U", ""),
            )
            stats["bam_primary_records"] += 1
            if read.is_unmapped:
                stats["bam_unmapped_primary_records"] += 1
                cluster_counts[cluster_key]["bam_unmapped_primary_records"] += 1
                trim_outcomes[trimmed]["bam_unmapped_primary_records"] += 1
                continue
            stats["mapped_primary_reads"] += 1
            cluster_counts[cluster_key]["mapped_primary_reads"] += 1
            trim_outcomes[trimmed]["mapped_primary_reads"] += 1
            sketches["mapped_primary"].add(*input_molecule)
            if read.mapping_quality != 255:
                stats["multimapping_primary_reads"] += 1
                cluster_counts[cluster_key]["multimapping_primary_reads"] += 1
                trim_outcomes[trimmed]["multimapping_primary_reads"] += 1
                continue
            stats["unique_reads"] += 1
            cluster_counts[cluster_key]["unique_reads"] += 1
            trim_outcomes[trimmed]["unique_reads"] += 1
            sketches["unique"].add(*input_molecule)
            gene_id, gene_name, exonic_fraction = annotator.assign(read)
            if not gene_id:
                continue
            stats["false_gene_assigned_reads"] += 1
            stats[
                "false_gene_assigned_production_removed"
                if truth["production_removed"]
                else "false_gene_assigned_production_retained"
            ] += 1
            cluster_counts[cluster_key]["false_gene_assigned_reads"] += 1
            trim_outcomes[trimmed]["false_gene_assigned_reads"] += 1
            output_molecule = (
                truth.get("C", ""),
                truth.get("B", ""),
                gene_id,
                truth.get("U", ""),
            )
            sketches["gene_assigned"].add(*output_molecule)
            expected = bool(gene_name and gene_name in set(truth["targets"]))
            stats["guide_target_assigned_reads"] += int(expected)
            cluster_counts[cluster_key]["guide_target_assigned_reads"] += int(
                expected
            )
            trim_outcomes[trimmed]["guide_target_assigned_reads"] += int(expected)
            if expected:
                sketches["guide_target_assigned"].add(*output_molecule)
            if exact_collapse:
                group = (
                    truth.get("C", ""),
                    truth.get("B", ""),
                    gene_id,
                )
                umi = truth.get("U", "")
                payload = assigned_umi_payloads[group].setdefault(umi, Counter())
                payload["total"] += 1
                payload["guide_target"] += int(expected)
                payload["production_removed"] += int(
                    truth["production_removed"]
                )
            gene_counts[(gene_id, gene_name, expected, exonic_fraction > 0.5)] += 1

    molecule_estimates = {key: sketch.estimate() for key, sketch in sketches.items()}
    input_molecules = molecule_estimates["input"]
    exact_counts = Counter()
    controlled_gene_molecule_counts = Counter()
    cell_exact_counts: Dict[str, Counter] = defaultdict(Counter)
    if exact_collapse:
        for (cell, pb, _guide), payloads in input_umi_payloads.items():
            for cluster in collapse_umi_payloads(payloads):
                exact_counts["input"] += 1
                exact_counts["input_tso_trimmed_any"] += int(
                    cluster["tso_trimmed"] > 0
                )
                exact_counts["input_tso_trimmed_all"] += int(
                    cluster["tso_trimmed"] == cluster["total"]
                )
                cell_exact_counts[cell]["input"] += 1
                cluster_counts[(cell, pb)]["input_molecules"] += 1
        for (cell, pb, gene), payloads in assigned_umi_payloads.items():
            for cluster in collapse_umi_payloads(payloads):
                exact_counts["gene_assigned"] += 1
                exact_counts["guide_target_any"] += int(
                    cluster["guide_target"] > 0
                )
                exact_counts["guide_target_all"] += int(
                    cluster["guide_target"] == cluster["total"]
                )
                exact_counts["guide_target_majority"] += int(
                    2 * cluster["guide_target"] > cluster["total"]
                )
                cell_exact_counts[cell]["gene_assigned"] += 1
                cluster_counts[(cell, pb)]["false_gene_assigned_molecules"] += 1
                cell_exact_counts[cell]["guide_target_any"] += int(
                    cluster["guide_target"] > 0
                )
                cluster_counts[(cell, pb)]["guide_target_molecules_any"] += int(
                    cluster["guide_target"] > 0
                )
                stratum = (
                    "removed" if cluster["production_removed"] > 0 else "retained"
                )
                controlled_gene_molecule_counts[(gene, stratum)] += 1
    exact_input_molecules = int(exact_counts["input"])
    exact_assigned_molecules = int(exact_counts["gene_assigned"])
    nominal_false_read_ci = wilson_interval(
        int(stats["false_gene_assigned_reads"]), int(stats["input_reads"])
    )
    nominal_target_read_ci = wilson_interval(
        int(stats["guide_target_assigned_reads"]), int(stats["input_reads"])
    )
    bootstrap_replicates = int(
        alignment_cfg.get("cluster_bootstrap_replicates", 1000)
    )
    false_read_cluster_ci = cluster_bootstrap_ratio_interval(
        cluster_counts,
        "false_gene_assigned_reads",
        "input_reads",
        bootstrap_replicates,
        f"{library_name}|{variant}|false_reads",
    )
    target_read_cluster_ci = cluster_bootstrap_ratio_interval(
        cluster_counts,
        "guide_target_assigned_reads",
        "input_reads",
        bootstrap_replicates,
        f"{library_name}|{variant}|target_reads",
    )
    false_molecule_cluster_ci = cluster_bootstrap_ratio_interval(
        cluster_counts,
        "false_gene_assigned_molecules",
        "input_molecules",
        bootstrap_replicates,
        f"{library_name}|{variant}|false_molecules",
    )
    summary = {
        "library": library_name,
        "variant": variant,
        "star_input_fastq": input_fastq,
        "star_bam": bam_path,
        "star_final_log": star_log_path,
        "truth_panel_token": truth_panel_token,
        "benchmark_expected_input_reads": expected_truth_reads,
        "benchmark_qc_module_sha256": benchmark_module_digest,
        "benchmark_config_sha256": benchmark_config_digest,
        "alignment_audit_qc_module_sha256": file_sha256(
            os.path.abspath(__file__)
        ),
        "input_reads": int(stats["input_reads"]),
        "bam_primary_records": int(stats["bam_primary_records"]),
        "bam_unmapped_primary_records": int(stats["bam_unmapped_primary_records"]),
        "mapped_primary_reads": int(stats["mapped_primary_reads"]),
        "unique_reads": int(stats["unique_reads"]),
        "multimapping_primary_reads": int(stats["multimapping_primary_reads"]),
        "false_gene_assigned_reads": int(stats["false_gene_assigned_reads"]),
        "guide_target_assigned_reads": int(stats["guide_target_assigned_reads"]),
        "production_filter_stratification": {
            "removed": {
                "input_reads": int(stats["input_production_removed"]),
                "false_gene_assigned_reads": int(
                    stats["false_gene_assigned_production_removed"]
                ),
            },
            "retained": {
                "input_reads": int(stats["input_production_retained"]),
                "false_gene_assigned_reads": int(
                    stats["false_gene_assigned_production_retained"]
                ),
            },
        },
        "false_gene_assignment_fraction_of_input_reads": (
            stats["false_gene_assigned_reads"] / stats["input_reads"]
            if stats["input_reads"]
            else 0.0
        ),
        "false_gene_assignment_fraction_of_input_reads_cluster_bootstrap_95": list(
            false_read_cluster_ci
        ),
        "false_gene_assignment_fraction_of_input_reads_naive_wilson_95": list(
            nominal_false_read_ci
        ),
        "guide_target_assignment_fraction_of_input_reads": (
            stats["guide_target_assigned_reads"] / stats["input_reads"]
            if stats["input_reads"]
            else 0.0
        ),
        "guide_target_assignment_fraction_of_input_reads_cluster_bootstrap_95": list(
            target_read_cluster_ci
        ),
        "guide_target_assignment_fraction_of_input_reads_naive_wilson_95": list(
            nominal_target_read_ci
        ),
        "cluster_sampling_uncertainty": {
            "sampling_units": len(cluster_counts),
            "bootstrap_replicates": bootstrap_replicates,
            "primary_interval_method": (
                "percentile bootstrap over corrected (CB,pB) sampling units"
            ),
            "naive_wilson_note": (
                "Read-level Wilson intervals ignore CB/pB clustering and PCR "
                "dependence; retained only as a nominal diagnostic."
            ),
        },
        "estimated_exact_key_molecules_hll": molecule_estimates,
        "hll_key_definitions": {
            "input_mapped_primary_unique": "(CB, pB, ground-truth guide, exact UMI)",
            "gene_assigned_guide_target_assigned": "(CB, pB, assigned gene, exact UMI)",
        },
        "false_gene_assignment_fraction_of_input_molecules": (
            molecule_estimates["gene_assigned"] / input_molecules
            if input_molecules
            else 0.0
        ),
        "guide_target_assignment_fraction_of_input_molecules": (
            molecule_estimates["guide_target_assigned"] / input_molecules
            if input_molecules
            else 0.0
        ),
        "splitpipe_equivalent_hamming1_molecule_impact": {
            "enabled": exact_collapse,
            "ground_truth_input_molecules_by_guide": exact_input_molecules,
            "ground_truth_input_molecules_tso_trimmed_any": int(
                exact_counts["input_tso_trimmed_any"]
            ),
            "ground_truth_input_molecules_tso_trimmed_all": int(
                exact_counts["input_tso_trimmed_all"]
            ),
            "false_gene_assigned_molecules": exact_assigned_molecules,
            "guide_target_molecules_any_read": int(
                exact_counts["guide_target_any"]
            ),
            "guide_target_molecules_all_reads": int(
                exact_counts["guide_target_all"]
            ),
            "guide_target_molecules_majority_reads": int(
                exact_counts["guide_target_majority"]
            ),
            "false_gene_assignment_fraction_of_input_molecules": (
                exact_assigned_molecules / exact_input_molecules
                if exact_input_molecules
                else 0.0
            ),
            "false_gene_assignment_fraction_of_input_molecules_cluster_bootstrap_95": list(
                false_molecule_cluster_ci
            ),
            "guide_target_assignment_fraction_of_input_molecules": (
                exact_counts["guide_target_any"] / exact_input_molecules
                if exact_input_molecules
                else 0.0
            ),
            "input_grouping": (
                "(CB, pB, ground-truth guide), followed by split-pipe's "
                "abundance-directed Hamming-distance-1 UMI collapse"
            ),
            "assigned_grouping": (
                "(CB, pB, assigned gene), followed by split-pipe's "
                "abundance-directed Hamming-distance-1 UMI collapse"
            ),
        },
        "input_processed_length_counts": {
            str(key): int(value) for key, value in input_length_counts.items()
        },
        "predicted_tso_trim_input_counts": {
            str(key): int(value) for key, value in trim_input_counts.items()
        },
        "outcomes_by_predicted_tso_trim": {
            str(trimmed): {key: int(value) for key, value in counts.items()}
            for trimmed, counts in trim_outcomes.items()
        },
        "star_log_metrics": parse_star_final_log(star_log_path),
        "interpretation": (
            "Every input record is a full-length structural gRNA contaminant. Any "
            "split-pipe-equivalent gene assignment is therefore a false GEX count; "
            "guide-target assignments are the subset most likely to mimic the "
            "intended perturbation's expression signal."
        ),
    }
    summary_path = os.path.join(star_dir, "alignment_impact_summary.yaml")
    with open(summary_path, "w") as handle:
        yaml.safe_dump(summary, handle, sort_keys=False)
    write_rows_csv(
        [
            {
                "gene_id": gene_id,
                "gene_name": gene_name,
                "is_guide_target": expected,
                "region": "exonic" if exonic else "intronic",
                "read_count": count,
            }
            for (gene_id, gene_name, expected, exonic), count in gene_counts.items()
        ],
        ["gene_id", "gene_name", "is_guide_target", "region", "read_count"],
        os.path.join(star_dir, "false_gene_assignments.csv"),
        sort_by=["read_count"],
        ascending=False,
    )
    write_rows_csv(
        [
            {
                "cell": cell,
                "ground_truth_input_molecules": int(counts["input"]),
                "false_gene_assigned_molecules": int(counts["gene_assigned"]),
                "guide_target_molecules_any_read": int(
                    counts["guide_target_any"]
                ),
                "false_gene_assignment_fraction": (
                    counts["gene_assigned"] / counts["input"]
                    if counts["input"]
                    else 0.0
                ),
            }
            for cell, counts in cell_exact_counts.items()
        ],
        [
            "cell",
            "ground_truth_input_molecules",
            "false_gene_assigned_molecules",
            "guide_target_molecules_any_read",
            "false_gene_assignment_fraction",
        ],
        os.path.join(star_dir, "cell_false_molecule_impact.csv"),
        sort_by=["false_gene_assigned_molecules", "ground_truth_input_molecules"],
        ascending=False,
    )
    gene_name_by_id = {
        gene_id: gene_name
        for gene_id, gene_name, _expected, _exonic in gene_counts
    }
    write_rows_csv(
        [
            {
                "gene_id": gene_id,
                "gene_name": gene_name_by_id.get(gene_id, ""),
                "production_filter_stratum": stratum,
                "false_gene_assigned_molecules": count,
            }
            for (gene_id, stratum), count in controlled_gene_molecule_counts.items()
        ],
        [
            "gene_id", "gene_name", "production_filter_stratum",
            "false_gene_assigned_molecules",
        ],
        os.path.join(star_dir, "false_gene_assigned_molecules_by_gene.csv"),
        sort_by=["false_gene_assigned_molecules"],
        ascending=False,
    )
    write_rows_csv(
        [
            {
                "cell": cell,
                "barcode_indices": pb,
                "input_reads": int(counts["input_reads"]),
                "mapped_primary_reads": int(counts["mapped_primary_reads"]),
                "multimapping_primary_reads": int(
                    counts["multimapping_primary_reads"]
                ),
                "unique_reads": int(counts["unique_reads"]),
                "false_gene_assigned_reads": int(
                    counts["false_gene_assigned_reads"]
                ),
                "guide_target_assigned_reads": int(
                    counts["guide_target_assigned_reads"]
                ),
                "input_molecules": int(counts["input_molecules"]),
                "false_gene_assigned_molecules": int(
                    counts["false_gene_assigned_molecules"]
                ),
                "guide_target_molecules_any": int(
                    counts["guide_target_molecules_any"]
                ),
            }
            for (cell, pb), counts in cluster_counts.items()
        ],
        [
            "cell",
            "barcode_indices",
            "input_reads",
            "mapped_primary_reads",
            "multimapping_primary_reads",
            "unique_reads",
            "false_gene_assigned_reads",
            "guide_target_assigned_reads",
            "input_molecules",
            "false_gene_assigned_molecules",
            "guide_target_molecules_any",
        ],
        os.path.join(star_dir, "cb_pb_cluster_alignment_outcomes.csv"),
        sort_by=["false_gene_assigned_reads", "input_reads"],
        ascending=False,
    )
    print(f"  Wrote controlled alignment impact: {summary_path}")
    return summary


def filter_library(cfg: dict, config_path: str, library_name: str):
    qc = cfg["gex_grna_contamination_qc"]
    short_cfg = qc["short_read_filters"]
    selected_method = str(short_cfg["selected_method"])
    read_length = int(qc.get("short_read_length", 54))
    library = get_library(cfg, library_name)
    pairs = resolve_library_pairs(cfg, config_path, library)
    model = make_model(cfg, config_path)
    output_root = resolve_path(qc["output_dir"], config_path)
    filter_dir = os.path.join(output_root, library_name, f"filtered_{read_length}")
    os.makedirs(filter_dir, exist_ok=True)
    clean_r1 = os.path.join(
        filter_dir, f"{library_name}_clean_R1_{read_length}.fastq.gz"
    )
    clean_r2 = os.path.join(filter_dir, f"{library_name}_clean_R2.fastq.gz")
    reject_r1 = os.path.join(
        filter_dir, f"{library_name}_rejected_R1_{read_length}.fastq.gz"
    )
    reject_r2 = os.path.join(filter_dir, f"{library_name}_rejected_R2.fastq.gz")
    stats = Counter()
    with gzip.open(clean_r1, "wt") as c1, gzip.open(clean_r2, "wt") as c2, gzip.open(
        reject_r1, "wt"
    ) as x1, gzip.open(reject_r2, "wt") as x2:
        for r1_path, r2_path in pairs:
            for r1, r2 in iter_fastq_pairs(r1_path, r2_path):
                stats["total"] += 1
                short = model.short_features(r1.sequence, read_length, short_cfg)
                if selected_method not in short["calls"]:
                    raise KeyError(
                        f"selected_method '{selected_method}' is not available: "
                        f"{sorted(short['calls'])}"
                    )
                rejected = bool(short["calls"][selected_method])
                stats["rejected"] += int(rejected)
                out1, out2 = (x1, x2) if rejected else (c1, c2)
                write_fastq_record(out1, r1, read_length)
                write_fastq_record(out2, r2)
                if stats["total"] % 1_000_000 == 0:
                    print(
                        f"  {library_name} filter: {stats['total']:,} reads; "
                        f"{stats['rejected']:,} rejected",
                        flush=True,
                    )
    summary = {
        "library": library_name,
        "method": selected_method,
        "short_read_length": read_length,
        "total": int(stats["total"]),
        "rejected": int(stats["rejected"]),
        "rejected_fraction": stats["rejected"] / stats["total"] if stats["total"] else 0.0,
        "clean_r1": clean_r1,
        "clean_r2": clean_r2,
        "rejected_r1": reject_r1,
        "rejected_r2": reject_r2,
    }
    with open(os.path.join(filter_dir, "filter_summary.yaml"), "w") as handle:
        yaml.safe_dump(summary, handle, sort_keys=False)
    return summary


def make_aggregate_plots(
    summary_df: pd.DataFrame,
    metrics_df: pd.DataFrame,
    alignment_df: pd.DataFrame,
    output_root: str,
):
    """Create compact pilot-level plots; CSV/YAML remain the source of truth."""
    plot_dir = os.path.join(output_root, "plots")
    os.makedirs(plot_dir, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", plot_dir)
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  WARNING: matplotlib unavailable; skipping aggregate plots")
        return

    if not summary_df.empty:
        frame = summary_df.sort_values("library")
        values = frame["contamination_fraction"].astype(float)
        lower = values - frame["wilson_95_low"].astype(float)
        upper = frame["wilson_95_high"].astype(float) - values
        fig, ax = plt.subplots(figsize=(max(5, len(frame) * 1.2), 4))
        ax.bar(frame["library"], values, color="#35618f")
        ax.errorbar(
            frame["library"], values, yerr=[lower, upper], fmt="none", color="black", capsize=4
        )
        ax.set_ylabel("Full-length structural contamination fraction")
        ax.set_xlabel("GEX sublibrary")
        ax.set_title("gRNA carry-over in pilot GEX R1")
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
        fig.tight_layout()
        fig.savefig(os.path.join(plot_dir, "contamination_fraction.png"), dpi=180)
        plt.close(fig)

    if not metrics_df.empty:
        frame = metrics_df[metrics_df["truth_definition"] == "structure_truth"].copy()
        if not frame.empty:
            grouped = (
                frame.groupby("method", as_index=False)[
                    ["sensitivity", "hard_negative_fpr"]
                ]
                .mean(numeric_only=True)
                .dropna()
            )
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.scatter(
                grouped["hard_negative_fpr"],
                grouped["sensitivity"],
                color="#b44b3f",
                alpha=0.8,
            )
            for row in grouped.itertuples(index=False):
                ax.annotate(
                    row.method,
                    (row.hard_negative_fpr, row.sensitivity),
                    fontsize=6,
                    xytext=(3, 2),
                    textcoords="offset points",
                )
            ax.set_xscale("symlog", linthresh=1e-6)
            ax.set_xlim(left=0)
            ax.set_ylim(0, 1.03)
            ax.set_xlabel("FPR among TSO-positive, structure-negative reads")
            ax.set_ylabel("Sensitivity vs full-length structural truth")
            ax.set_title("54-cycle filter trade-off (mean across libraries)")
            fig.tight_layout()
            fig.savefig(os.path.join(plot_dir, "short_read_filter_tradeoff.png"), dpi=180)
            plt.close(fig)

    if not alignment_df.empty:
        value_column = "false_gene_assignment_fraction_of_input_reads"
        if value_column in alignment_df.columns:
            frame = alignment_df.copy()
            variants = list(dict.fromkeys(frame["variant"].astype(str)))
            libraries = sorted(frame["library"].astype(str).unique())
            width = 0.8 / max(1, len(variants))
            positions = list(range(len(libraries)))
            fig, ax = plt.subplots(figsize=(max(6, len(libraries) * 1.4), 4.5))
            colors = ["#35618f", "#b44b3f", "#6f8f3d"]
            for variant_index, variant in enumerate(variants):
                by_library = (
                    frame[frame["variant"] == variant]
                    .set_index("library")[value_column]
                    .to_dict()
                )
                offsets = [
                    position
                    + (variant_index - (len(variants) - 1) / 2) * width
                    for position in positions
                ]
                ax.bar(
                    offsets,
                    [float(by_library.get(library, 0.0)) for library in libraries],
                    width=width,
                    label=variant,
                    color=colors[variant_index % len(colors)],
                )
            ax.set_xticks(positions, libraries)
            ax.set_ylim(bottom=0)
            ax.set_ylabel("False gene assignments / truth reads")
            ax.set_xlabel("GEX sublibrary")
            ax.set_title("Ground-truth gRNA reads under current and short R1")
            ax.legend(title="STAR input arm")
            fig.tight_layout()
            fig.savefig(
                os.path.join(plot_dir, "controlled_alignment_false_assignment.png"),
                dpi=180,
            )
            plt.close(fig)


def _report_value(value, column: str) -> str:
    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return "NA"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, (int, float)):
        if "fraction" in column or column.endswith("_fpr"):
            return f"{float(value):.3%}"
        if float(value).is_integer():
            return f"{int(value):,}"
        return f"{float(value):.4g}"
    return str(value).replace("|", "\\|")


def _markdown_table(frame: pd.DataFrame, columns: Sequence[str]) -> str:
    if frame.empty:
        return "_No completed outputs were found._"
    columns = [column for column in columns if column in frame.columns]
    header = "| " + " | ".join(columns) + " |"
    rule = "| " + " | ".join("---" for _ in columns) + " |"
    rows = [header, rule]
    for record in frame[columns].to_dict("records"):
        rows.append(
            "| "
            + " | ".join(_report_value(record[column], column) for column in columns)
            + " |"
        )
    return "\n".join(rows)


def write_aggregate_audit_report(
    summary_df: pd.DataFrame,
    alignment_df: pd.DataFrame,
    delta_df: pd.DataFrame,
    completeness_df: pd.DataFrame,
    output_root: str,
    read_length: int,
):
    sections = [
        "# GEX gRNA Cross-Contamination Audit",
        "",
        "This report is generated from the pilot full-length structural truth scan, "
        f"the in-silico raw-R1={read_length} benchmark, and controlled STAR arms. "
        "Every controlled-alignment input is a ground-truth structural guide read; "
        "therefore every gene assignment in that arm is a false GEX assignment.",
        "",
        "## Raw-read prevalence",
        "",
        _markdown_table(
            summary_df,
            [
                "library",
                "raw_reads",
                "raw_tso_reads",
                "raw_structure_contaminant_reads",
                "raw_contaminant_fraction_of_all",
                "raw_contaminant_fraction_of_tso",
            ],
        ),
        "",
        "## Pre-alignment UMI-family prevalence",
        "",
        "The `cell_pb_umi` estimate is the collision-prone lower molecular-diversity "
        "bound. The `cell_pb_umi_insert` estimate separates distinct inserts but can "
        "split a PCR family at sequencing errors. Both use corrected CB/pB and exact "
        "UMIs; they are HLL estimates, not split-pipe Hamming-1 molecule counts.",
        "",
        _markdown_table(
            summary_df,
            [
                "library",
                "valid_barcode_reads",
                "cell_pb_umi_all",
                "cell_pb_umi_structure_contaminants",
                "cell_pb_umi_contaminant_fraction_of_all",
                "cell_pb_umi_contaminant_fraction_of_tso",
                "cell_pb_umi_insert_contaminant_fraction_of_all",
                "cell_pb_umi_insert_contaminant_fraction_of_tso",
            ],
        ),
        "",
        "## Controlled alignment impact",
        "",
        _markdown_table(
            alignment_df,
            [
                "library",
                "variant",
                "input_reads",
                "mapped_primary_reads",
                "unique_reads",
                "false_gene_assigned_reads",
                "false_gene_assignment_fraction_of_input_reads",
                "false_gene_assignment_fraction_of_input_reads_cluster_bootstrap_95_low",
                "false_gene_assignment_fraction_of_input_reads_cluster_bootstrap_95_high",
                "hamming1_false_gene_assigned_molecules",
                "hamming1_false_gene_assignment_fraction_of_input_molecules",
                "projected_false_gene_assigned_fraction_of_valid_reads",
            ],
        ),
        "",
        f"## Full versus short{read_length}",
        "",
        _markdown_table(
            delta_df,
            [
                "library",
                "full_read_false_assignment_fraction",
                f"short{read_length}_read_false_assignment_fraction",
                "read_false_assignment_absolute_delta",
                "read_absolute_delta_cluster_bootstrap_95_low",
                "read_absolute_delta_cluster_bootstrap_95_high",
                "read_false_assignment_fold_change",
                "full_hamming1_molecule_false_assignment_fraction",
                f"short{read_length}_hamming1_molecule_false_assignment_fraction",
                "hamming1_molecule_absolute_delta",
                "hamming1_molecule_absolute_delta_cluster_bootstrap_95_low",
                "hamming1_molecule_absolute_delta_cluster_bootstrap_95_high",
                "hamming1_molecule_fold_change",
            ],
        ),
        "",
        "## Interpretation constraints",
        "",
        "- The primary contamination label requires fuzzy full Split_TSO plus a "
        "spacer-sized segment, complete paired CSM stem, 7-nt iBAR-sized block, and "
        "paired partial tracr. TSO alone is never contaminant truth.",
        "- Raw-read fractions use every sequenced pair. UMI-family estimates require "
        "a valid split-pipe-equivalent R2 barcode call and therefore use a different "
        "denominator.",
        "- The controlled arms truncate raw R1 first and then reproduce split-pipe's "
        "exact 17-mer TSO trim before the unchanged STAR command. This preserves the "
        "ground-truth label even when only spacer/stem bases remain after trimming.",
        "- Projected whole-library burdens multiply contamination prevalence by the "
        "conditional ground-truth alignment rate. They are explicit projections; "
        "the controlled conditional rates are direct measurements, while the "
        "existing-BAM audit describes only mapped survivors when STAR-unmapped "
        "records were not retained.",
        "- Controlled false molecules are contaminant-only `(CB,pB,gene,UMI)` "
        "clusters. They are not guaranteed to be incremental DGE molecules after "
        "mixing with genuine GEX reads: a real-GEX collision can absorb or alter a "
        "Hamming-1 cluster. Historical-BAM any/all bounds measure that competition "
        "only among records retained in that BAM; short-read burden remains an "
        "isolated-arm projection.",
        "- A production filter should be selected only after reviewing sensitivity "
        "by TSO offset and the Wilson upper bound on false positives among "
        "TSO-positive, structure-negative GEX reads.",
        "",
        "## Completion status",
        "",
        _markdown_table(
            completeness_df,
            [
                "library",
                "benchmark_complete",
                "truth_safety_cap_reached",
                "bam_audit_expected",
                "bam_audit_required",
                "bam_audit_complete",
                "alignment_full_complete",
                f"alignment_short{read_length}_complete",
                "audit_complete",
            ],
        ),
        "",
    ]
    report_path = os.path.join(output_root, "GEX_GRNA_CONTAMINATION_AUDIT_REPORT.md")
    with open(report_path, "w") as handle:
        handle.write("\n".join(sections))
    return report_path


def benchmark_summary_status(
    summary_path: str,
    cfg: dict,
    library_name: str,
    current_module_digest: str,
) -> Tuple[bool, str]:
    if not os.path.exists(summary_path):
        return False, "missing"
    with open(summary_path) as handle:
        summary = yaml.safe_load(handle) or {}
    if not summary.get("raw_scan_complete", False):
        return False, "partial_raw_scan"
    truth_input = summary.get("alignment_ground_truth_input", {})
    if truth_input.get("qc_module_sha256", "") != current_module_digest:
        return False, "qc_module_mismatch"
    if truth_input.get("benchmark_config_sha256", "") != benchmark_config_sha256(
        cfg, library_name
    ):
        return False, "configuration_mismatch"
    metadata = summary.get("input_fastq_metadata", [])
    if not metadata:
        return False, "input_fastq_metadata_missing"
    for pair in metadata:
        for read_name in ("r1", "r2"):
            identity = pair.get(read_name, {})
            path = str(identity.get("path", ""))
            if not path or not os.path.exists(path):
                return False, f"input_{read_name}_missing"
            stat = os.stat(path)
            if int(identity.get("size_bytes", -1)) != stat.st_size:
                return False, f"input_{read_name}_size_changed"
            if int(identity.get("mtime_ns", -1)) != stat.st_mtime_ns:
                return False, f"input_{read_name}_mtime_changed"
    return True, "current"


def alignment_impact_status(
    impact_path: str,
    benchmark_summary: dict,
    current_module_digest: str,
    current_config_digest: str,
) -> Tuple[bool, str]:
    if not os.path.exists(impact_path):
        return False, "missing"
    truth_input = (benchmark_summary or {}).get("alignment_ground_truth_input", {})
    if truth_input.get("safety_cap_reached", False):
        return False, "truth_safety_cap_reached"
    expected_token = str(truth_input.get("truth_panel_token", ""))
    expected_reads = int(truth_input.get("reads_written_per_variant", 0))
    expected_module = str(truth_input.get("qc_module_sha256", ""))
    expected_config = str(truth_input.get("benchmark_config_sha256", ""))
    with open(impact_path) as handle:
        impact = yaml.safe_load(handle) or {}
    if not expected_token or impact.get("truth_panel_token", "") != expected_token:
        return False, "truth_panel_token_mismatch"
    if int(impact.get("input_reads", -1)) != expected_reads:
        return False, "input_read_count_mismatch"
    if not expected_module or expected_module != current_module_digest:
        return False, "benchmark_module_mismatch"
    if not expected_config or expected_config != current_config_digest:
        return False, "benchmark_config_mismatch"
    if impact.get("benchmark_qc_module_sha256", "") != expected_module:
        return False, "impact_benchmark_module_mismatch"
    if impact.get("alignment_audit_qc_module_sha256", "") != current_module_digest:
        return False, "alignment_audit_module_mismatch"
    if impact.get("benchmark_config_sha256", "") != current_config_digest:
        return False, "impact_config_mismatch"
    return True, "current"


def bam_audit_status(
    summary_path: str,
    current_module_digest: str,
    current_config_digest: str,
) -> Tuple[bool, str]:
    if not os.path.exists(summary_path):
        return False, "missing"
    with open(summary_path) as handle:
        summary = yaml.safe_load(handle) or {}
    if not summary.get("bam_scan_complete", False):
        return False, "partial_bam_scan"
    if summary.get("qc_module_sha256", "") != current_module_digest:
        return False, "qc_module_mismatch"
    if summary.get("qc_config_sha256", "") != current_config_digest:
        return False, "configuration_mismatch"
    bam_path = str(summary.get("bam", ""))
    if not bam_path or not os.path.exists(bam_path):
        return False, "source_bam_missing"
    stat = os.stat(bam_path)
    if int(summary.get("bam_size_bytes", -1)) != stat.st_size:
        return False, "source_bam_size_changed"
    if int(summary.get("bam_mtime_ns", -1)) != stat.st_mtime_ns:
        return False, "source_bam_mtime_changed"
    return True, "current"


def _base_output_root(cfg: dict, qc_config_path: str) -> str:
    qc_path = Path(qc_config_path).expanduser().resolve()
    with qc_path.open() as handle:
        qc_only = yaml.safe_load(handle) or {}
    base_raw = qc_only.get("base_config")
    base_path = (
        (qc_path.parent / base_raw).resolve()
        if base_raw and not Path(base_raw).is_absolute()
        else Path(base_raw).resolve() if base_raw else qc_path
    )
    raw = str(cfg["paths"]["output_dir"])
    return raw if os.path.isabs(raw) else str((base_path.parent / raw).resolve())


def _format_guide_path(template: str, cfg: dict, config_path: str, library: str) -> Path:
    output_dir = _base_output_root(cfg, config_path)
    rendered = str(template).format(output_dir=output_dir, library=library)
    path = Path(rendered).expanduser()
    return path.resolve() if path.is_absolute() else (Path(config_path).resolve().parent / path).resolve()


def _production_assignment_sets(frame: pd.DataFrame, dual_mode: bool) -> dict:
    """Return the identities emitted by production Step 4 for each cell."""
    assignments = {}
    if frame.empty:
        return assignments
    identity_column = "pair_id" if dual_mode else "grna_name"
    if identity_column not in frame.columns:
        raise ValueError(
            f"production assignment output lacks {identity_column!r}"
        )
    for cell, cell_frame in frame.groupby("bc_wells", sort=False):
        identities = {
            str(value)
            for value in cell_frame[identity_column].dropna().astype(str)
            if str(value)
        }
        if identities:
            assignments[str(cell)] = identities
    return assignments


def _read_match_header(path: Path) -> list[str]:
    return pd.read_csv(path, sep="\t", nrows=0).columns.tolist()


def _cross_source_overlap(gex_path: Path, dedicated_path: Path, library: str):
    required = {"bc_wells", "bc_bci", "umi", "grna_seq", "ibar_seq_corrected"}
    for label, path in (("GEX", gex_path), ("dedicated", dedicated_path)):
        missing = required - set(_read_match_header(path))
        if missing:
            raise ValueError(f"{label} match file lacks overlap columns: {sorted(missing)}")

    usecols = sorted(required)
    primary_gex = set()
    strict_gex = set()
    guide_only_gex = set()
    gex_rows = resolved_rows = 0
    for chunk in pd.read_csv(gex_path, sep="\t", usecols=usecols, chunksize=250_000):
        chunk = chunk.fillna("")
        for row in chunk.itertuples(index=False):
            gex_rows += 1
            ibar = str(row.ibar_seq_corrected)
            common = (library, str(row.bc_wells), str(row.grna_seq), str(row.umi))
            if ibar:
                resolved_rows += 1
                primary_gex.add((*common[:3], ibar, common[3]))
                strict_gex.add((*common[:3], ibar, common[3], str(row.bc_bci)))
            else:
                guide_only_gex.add(common)

    exact_overlap = set()
    strict_overlap = set()
    guide_only_overlap = set()
    for chunk in pd.read_csv(
        dedicated_path, sep="\t", usecols=usecols, chunksize=500_000
    ):
        chunk = chunk.fillna("")
        for row in chunk.itertuples(index=False):
            common = (library, str(row.bc_wells), str(row.grna_seq), str(row.umi))
            ibar = str(row.ibar_seq_corrected)
            if common in guide_only_gex:
                guide_only_overlap.add(common)
            if ibar:
                primary = (*common[:3], ibar, common[3])
                if primary in primary_gex:
                    exact_overlap.add(primary)
                    strict = (*common[:3], ibar, common[3], str(row.bc_bci))
                    if strict in strict_gex:
                        strict_overlap.add(strict)

    rows = []
    for key in sorted(exact_overlap):
        rows.append({
            "overlap_type": "guide_ibar_umi",
            "sublibrary": key[0],
            "cell_barcode": key[1],
            "guide_sequence": key[2],
            "corrected_ibar": key[3],
            "umi": key[4],
            "pb": "",
        })
    for key in sorted(strict_overlap):
        rows.append({
            "overlap_type": "guide_ibar_umi_pb",
            "sublibrary": key[0],
            "cell_barcode": key[1],
            "guide_sequence": key[2],
            "corrected_ibar": key[3],
            "umi": key[4],
            "pb": key[5],
        })
    for key in sorted(guide_only_overlap):
        rows.append({
            "overlap_type": "guide_only_no_ibar",
            "sublibrary": key[0],
            "cell_barcode": key[1],
            "guide_sequence": key[2],
            "corrected_ibar": "",
            "umi": key[3],
            "pb": "",
        })
    summary = {
        "gex_matched_reads": int(gex_rows),
        "gex_reads_with_resolved_ibar": int(resolved_rows),
        "gex_distinct_guide_ibar_umi_keys": len(primary_gex),
        "exact_cross_source_guide_ibar_umi_overlap": len(exact_overlap),
        "supporting_overlap_also_matching_pb": len(strict_overlap),
        "gex_distinct_guide_only_keys_without_ibar": len(guide_only_gex),
        "guide_only_cross_source_overlap_without_ibar": len(guide_only_overlap),
        "interpretation": (
            "Exact cross-source overlap is evidence, not proof, that the same RT "
            "molecule entered both size-selected libraries; barcode/UMI collisions "
            "remain possible. pB-matching overlap is stricter supporting evidence."
        ),
    }
    return pd.DataFrame(rows), summary


def validate_step1a_barcode_evidence_with_qc_decoder(
    cfg: dict,
    config_path: str,
    library_name: str,
    step1a_manifest_path,
    step1a: Optional[dict] = None,
    *,
    barcode_corrector=None,
):
    """Check Step 1a evidence against the QC decoder on rejected pairs only.

    The decoder fits the ordinary full-GEX 3–4M calibration interval, but it
    then decodes only the small synchronized Step 1a positive FASTQs. This is
    an independent release check, not a second full-library Step 1c scan.
    """

    step1a_manifest_path = Path(step1a_manifest_path).expanduser().resolve()
    step1a = step1a or json.loads(step1a_manifest_path.read_text())
    if step1a.get("complete") is not True or step1a.get("library") != library_name:
        raise ValueError("Step 1a manifest is incomplete or library-mismatched")

    def artifact(key):
        raw = step1a.get(key)
        if not raw:
            raise ValueError(f"Step 1a manifest is missing {key}")
        path = Path(raw).expanduser()
        if not path.is_absolute():
            path = step1a_manifest_path.parent / path
        path = path.resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Step 1a artifact is missing: {path}")
        expected = step1a.get(f"{key}_sha256")
        if not expected or file_sha256(str(path)) != expected:
            raise ValueError(f"Step 1a {key} SHA-256 mismatch")
        return path

    rejected_r1 = artifact("rejected_r1")
    rejected_r2 = artifact("rejected_r2")
    evidence = artifact("barcode_evidence")
    if barcode_corrector is None:
        library = get_library(cfg, library_name)
        full_pairs = resolve_library_pairs(cfg, config_path, library)
        barcode_corrector = make_barcode_corrector(
            cfg,
            config_path,
            library_name,
            [r2_path for _, r2_path in full_pairs],
        )

    expected_header = [
        "read_id", "bc_wells", "bc_bci", "umi", "r1_seq", "r1_qual"
    ]
    digest = hashlib.sha256()
    raw_pairs = valid_rows = 0
    with gzip.open(evidence, "rt", newline="") as evidence_handle:
        reader = csv.DictReader(evidence_handle, delimiter="\t")
        if list(reader.fieldnames or ()) != expected_header:
            raise ValueError("Step 1a barcode-evidence header is incompatible")
        expected_row = next(reader, None)
        for r1, r2 in iter_fastq_pairs(str(rejected_r1), str(rejected_r2)):
            raw_pairs += 1
            call = barcode_corrector.call(r2.sequence)
            if call is None or not call.valid:
                continue
            if expected_row is None:
                raise ValueError("QC decoder produced more valid rows than Step 1a")
            observed = {
                "read_id": normalize_read_id(r1.header),
                "bc_wells": str(call.cell),
                "bc_bci": str(call.barcode_indices),
                "umi": str(call.umi),
                "r1_seq": str(r1.sequence),
                "r1_qual": str(r1.quality),
            }
            if observed != expected_row:
                differing = [
                    key for key in expected_header
                    if observed.get(key) != expected_row.get(key)
                ]
                raise ValueError(
                    "Step 1a/QC barcode evidence mismatch at valid row "
                    f"{valid_rows + 1}: {differing}"
                )
            digest.update("\t".join(observed[key] for key in expected_header).encode())
            digest.update(b"\n")
            valid_rows += 1
            expected_row = next(reader, None)
        if expected_row is not None:
            raise ValueError("Step 1a evidence has rows absent from the QC decoder")

    expected_pairs = int(step1a.get("rejected_pair_count", -1))
    expected_rows = int(step1a.get("barcode_evidence_rows", -1))
    if raw_pairs != expected_pairs or valid_rows != expected_rows:
        raise ValueError(
            "Step 1a/QC evidence counts differ: "
            f"pairs {raw_pairs}/{expected_pairs}, valid {valid_rows}/{expected_rows}"
        )
    return {
        "complete": True,
        "scope": "Step 1a rejected pairs only after full-GEX calibration",
        "rejected_pairs_decoded": raw_pairs,
        "valid_rows_compared": valid_rows,
        "ordered_canonical_sha256": digest.hexdigest(),
        "all_fields_identical": True,
    }


def audit_guide_assignment_impact(
    cfg: dict, config_path: str, library_name: str, threads: int = 1
):
    del threads  # Step 4 streams one library; UGER binding still limits the job.
    qc = cfg["gex_grna_contamination_qc"]
    guide_cfg = qc.get("guide_assignment_impact", {}) or {}
    if not bool(guide_cfg.get("enabled", False)):
        raise ValueError("guide_assignment_impact.enabled is false")
    step1c_template = guide_cfg.get(
        "step1c_manifest_template",
        "{output_dir}/gex_grna_identity/manifests/{library}.step1c.complete.json",
    )
    dedicated_template = guide_cfg.get(
        "dedicated_match_template",
        "{output_dir}/grna_matches/{library}_grna_matches.tsv.gz",
    )
    step1c_manifest_path = _format_guide_path(
        step1c_template, cfg, config_path, library_name
    )
    dedicated_path = _format_guide_path(
        dedicated_template, cfg, config_path, library_name
    )
    if not step1c_manifest_path.is_file() or not dedicated_path.is_file():
        raise FileNotFoundError(
            f"guide-impact inputs missing: {step1c_manifest_path}, {dedicated_path}"
        )
    step1c_manifest = json.loads(step1c_manifest_path.read_text())
    if step1c_manifest.get("complete") is not True:
        raise ValueError("Step 1c completion manifest is incomplete")
    gex_path = Path(step1c_manifest["match_output"]).expanduser().resolve()
    if not gex_path.is_file():
        raise FileNotFoundError(f"Step 1c match output is missing: {gex_path}")
    expected_gex_hash = step1c_manifest.get("match_output_sha256")
    observed_gex_hash = file_sha256(str(gex_path))
    if not expected_gex_hash or expected_gex_hash != observed_gex_hash:
        raise ValueError("Step 1c match-output SHA-256 mismatch")

    step4 = _load_step4_module()
    ga = cfg["guide_assignment"]
    dedup_ibar = bool(ga.get("umi_dedup_ibar", True))
    group_cols = step4._determine_group_cols(str(dedicated_path), dedup_ibar)
    dual_mode = bool(ga.get("dual_guide_counting_mode", False))
    library = get_library(cfg, library_name)
    sublib_index = next(
        index for index, row in enumerate(cfg["libraries"])
        if row["sublibrary_name"] == library_name
    )
    combined, merge_stats = step4.stream_library_umi_dedup_sources(
        [("dedicated", dedicated_path), ("gex", gex_path)],
        library_name,
        library.get("sample_name", library_name),
        sublib_index,
        group_cols=group_cols,
        max_umi_edits=int(ga.get("max_umi_edits", 1)),
        min_reads_per_umi=int(ga.get("min_reads_per_umi", 1)),
        require_resolved_ibar_sources={"gex"} if dual_mode else set(),
    )

    matrix_columns = [
        "bc_wells", "grna_name", "target_gene", "sublibrary",
        "exact_umi_count_dedicated", "umi_count_dedicated",
        "exact_umi_count_gex", "umi_count_gex", "umi_count", "read_count",
        "exact_umi_count_overlap", "umi_count_both_sources",
    ]
    for column in matrix_columns:
        if column not in combined:
            combined[column] = 0 if "count" in column else ""
    matrix = combined[matrix_columns].copy()
    matrix["combined_minus_dedicated_collapsed_umis"] = (
        matrix["umi_count"] - matrix["umi_count_dedicated"]
    )
    per_cell_guide = matrix.groupby(
        ["bc_wells", "grna_name", "target_gene"], as_index=False
    ).sum(numeric_only=True)
    per_guide = per_cell_guide.groupby(
        ["grna_name", "target_gene"], as_index=False
    ).sum(numeric_only=True)

    reference_path = resolve_path(
        cfg["paths"]["grna_reference"],
        cfg.get("_base_config_path") or config_path,
    )
    ref_df = pd.read_csv(reference_path) if os.path.exists(reference_path) else None
    annotation_columns = [
        "grna_name", "position", "stem_variant", "pair_id",
        "target_combo", "target_version", "gi_pair", "gi_version",
    ]
    available_columns = (
        [column for column in annotation_columns if column in ref_df.columns]
        if ref_df is not None
        else []
    )
    dedicated_for_assignment = combined.copy()
    dedicated_for_assignment["umi_count"] = dedicated_for_assignment[
        "umi_count_dedicated"
    ]
    dedicated_for_assignment["umi_count_raw"] = dedicated_for_assignment[
        "exact_umi_count_dedicated"
    ]
    combined_for_assignment = combined.copy()
    dedicated_for_assignment, dedicated_assignment_frame = (
        step4._annotate_and_assign(
            dedicated_for_assignment,
            ref_df,
            available_columns,
            group_cols,
            ga,
        )
    )
    combined_for_assignment, combined_assignment_frame = step4._annotate_and_assign(
        combined_for_assignment,
        ref_df,
        available_columns,
        group_cols,
        ga,
    )
    dedicated_assignments = _production_assignment_sets(
        dedicated_assignment_frame, dual_mode
    )
    combined_assignments = _production_assignment_sets(
        combined_assignment_frame, dual_mode
    )
    assignment_rows = []
    assignment_status = Counter()
    for cell in sorted(set(dedicated_assignments) | set(combined_assignments)):
        before = dedicated_assignments.get(cell, set())
        after = combined_assignments.get(cell, set())
        added = after - before
        removed = before - after
        if not before and after:
            status = "added"
        elif before and not after:
            status = "removed"
        elif before != after:
            status = "changed"
        else:
            status = "unchanged"
        assignment_status[status] += 1
        assignment_rows.append({
            "cell_barcode": cell,
            "dedicated_assignment": ",".join(sorted(before)),
            "combined_assignment": ",".join(sorted(after)),
            "guides_added": ",".join(sorted(added)),
            "guides_removed": ",".join(sorted(removed)),
            "status": status,
            "assignment_identity": "pair_id" if dual_mode else "grna_name",
        })

    overlap_rows, overlap_summary = _cross_source_overlap(
        gex_path, dedicated_path, library_name
    )
    output_root = Path(resolve_path(qc["output_dir"], config_path))
    output_dir = output_root / library_name / guide_cfg.get(
        "output_subdir", "guide_assignment_impact"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(output_dir / "combined_umi_counts_by_source.csv.gz", index=False)
    per_cell_guide.to_csv(output_dir / "per_cell_guide_umi_delta.csv", index=False)
    per_guide.to_csv(output_dir / "per_guide_umi_delta.csv", index=False)
    pd.DataFrame(assignment_rows).to_csv(
        output_dir / "cell_assignment_changes.csv", index=False
    )
    overlap_rows.to_csv(output_dir / "cross_source_molecule_overlap.csv", index=False)

    step1a_manifest = Path(step1c_manifest["step1a_manifest"])
    step1a = json.loads(step1a_manifest.read_text())
    barcode_evidence_validation = validate_step1a_barcode_evidence_with_qc_decoder(
        cfg, config_path, library_name, step1a_manifest, step1a
    )
    summary = {
        "schema_version": 1,
        "complete": True,
        "library": library_name,
        "qc_module_sha256": file_sha256(os.path.abspath(__file__)),
        "benchmark_config_sha256": benchmark_config_sha256(cfg, library_name),
        "inputs": {
            "step1c_manifest": str(step1c_manifest_path),
            "step1c_manifest_sha256": file_sha256(str(step1c_manifest_path)),
            "gex_match": str(gex_path),
            "gex_match_sha256": observed_gex_hash,
            "dedicated_match": str(dedicated_path),
            "dedicated_match_sha256": file_sha256(str(dedicated_path)),
        },
        "gex_positive_reads_with_valid_barcodes": int(
            step1a.get("barcode_evidence_rows", 0)
        ),
        "barcode_evidence_qc_validation": barcode_evidence_validation,
        "gex_matched_guide_reads": overlap_summary["gex_matched_reads"],
        "gex_reads_with_resolved_ibar": overlap_summary[
            "gex_reads_with_resolved_ibar"
        ],
        "molecule_counts": {
            "dedicated_exact": int(matrix["exact_umi_count_dedicated"].sum()),
            "dedicated_collapsed": int(matrix["umi_count_dedicated"].sum()),
            "gex_exact": int(matrix["exact_umi_count_gex"].sum()),
            "gex_collapsed": int(matrix["umi_count_gex"].sum()),
            "combined_collapsed": int(matrix["umi_count"].sum()),
            "combined_increment": int(
                matrix["combined_minus_dedicated_collapsed_umis"].sum()
            ),
        },
        "assignment_changes": {
            key: int(value) for key, value in assignment_status.items()
        },
        "assignment_semantics": {
            "production_step4": True,
            "identity": "pair_id" if dual_mode else "grna_name",
            "min_umi_threshold": int(ga.get("min_umi_threshold", 5)),
            "multiplet_strategy": ga.get("multiplet_strategy"),
            "min_dominance_fraction": float(
                ga.get("min_dominance_fraction", 0.0)
            ),
        },
        "cross_source_overlap": overlap_summary,
        "production_merge_statistics": merge_stats,
        "semantics": {
            "combine_order": "reads_before_depth_filter_and_directional_hamming1_collapse",
            "primary_overlap_key": "(sublibrary, cell, guide sequence, corrected iBAR, exact UMI)",
            "supporting_overlap_key": "primary key plus corrected pB",
        },
    }
    with (output_dir / "guide_assignment_impact_summary.yaml").open("w") as handle:
        yaml.safe_dump(summary, handle, sort_keys=False)
    print(f"  Guide-assignment impact written under {output_dir}")
    return summary


def guide_impact_status(path, module_digest, config_digest):
    if not os.path.exists(path):
        return False, "missing"
    with open(path) as handle:
        payload = yaml.safe_load(handle) or {}
    if payload.get("complete") is not True:
        return False, "incomplete"
    if payload.get("qc_module_sha256") != module_digest:
        return False, "stale QC module hash"
    if payload.get("benchmark_config_sha256") != config_digest:
        return False, "stale configuration hash"
    return True, "complete"


def aggregate(cfg: dict, config_path: str):
    qc = cfg["gex_grna_contamination_qc"]
    output_root = resolve_path(qc["output_dir"], config_path)
    os.makedirs(output_root, exist_ok=True)
    read_length = int(qc.get("short_read_length", 54))
    summary_rows = []
    benchmark_by_library = {}
    metric_frames = []
    stratum_frames = []
    bam_rows = []
    gene_frames = []
    alignment_rows = []
    guide_rows = []
    burden_frames = []
    current_module_digest = file_sha256(os.path.abspath(__file__))
    alignment_cfg = qc.get("alignment_benchmark", {})
    variants = [
        f"short{read_length}" if str(value) == "short54" else str(value)
        for value in alignment_cfg.get("variants", ["full", f"short{read_length}"])
    ]
    for library in cfg.get("libraries", []):
        name = library["sublibrary_name"]
        lib_dir = os.path.join(output_root, name)
        guide_summary_path = os.path.join(
            lib_dir,
            qc.get("guide_assignment_impact", {}).get(
                "output_subdir", "guide_assignment_impact"
            ),
            "guide_assignment_impact_summary.yaml",
        )
        guide_current, _ = guide_impact_status(
            guide_summary_path,
            current_module_digest,
            benchmark_config_sha256(cfg, name),
        )
        if guide_current:
            with open(guide_summary_path) as handle:
                guide_summary = yaml.safe_load(handle) or {}
            guide_rows.append({
                "library": name,
                "gex_positive_reads_with_valid_barcodes": guide_summary.get(
                    "gex_positive_reads_with_valid_barcodes", 0
                ),
                "gex_matched_guide_reads": guide_summary.get(
                    "gex_matched_guide_reads", 0
                ),
                "gex_reads_with_resolved_ibar": guide_summary.get(
                    "gex_reads_with_resolved_ibar", 0
                ),
                **{
                    f"molecules_{key}": value
                    for key, value in (guide_summary.get("molecule_counts", {}) or {}).items()
                },
                **{
                    f"assignments_{key}": value
                    for key, value in (guide_summary.get("assignment_changes", {}) or {}).items()
                },
                **{
                    f"overlap_{key}": value
                    for key, value in (guide_summary.get("cross_source_overlap", {}) or {}).items()
                    if isinstance(value, (int, float))
                },
            })
        benchmark_path = os.path.join(lib_dir, "benchmark_summary.yaml")
        benchmark_current, benchmark_reason = benchmark_summary_status(
            benchmark_path, cfg, name, current_module_digest
        )
        if benchmark_current:
            with open(benchmark_path) as handle:
                summary = yaml.safe_load(handle)
            benchmark_by_library[name] = summary
            truth = summary["full_length_surrogate_truth"]
            raw = summary.get("raw_read_prevalence", {})
            molecule = summary.get("prealignment_molecule_prevalence") or {}
            valid_counts = molecule.get("valid_barcode_read_counts", {})
            estimates = molecule.get("estimated_distinct_keys", {})
            fractions = molecule.get("fractions", {})
            lower = estimates.get("cell_pb_umi", {})
            upper = estimates.get("cell_pb_umi_insert", {})
            lower_fractions = fractions.get("cell_pb_umi", {})
            upper_fractions = fractions.get("cell_pb_umi_insert", {})
            valid_all = int(valid_counts.get("all", 0))
            valid_tso = int(valid_counts.get("tso", 0))
            valid_contaminants = int(
                valid_counts.get("structure_contaminant", 0)
            )
            summary_rows.append(
                {
                    "library": name,
                    "records_evaluated": summary["records_evaluated"],
                    "raw_reads": int(raw.get("all_reads", summary["records_evaluated"])),
                    "raw_tso_reads": int(raw.get("tso_reads", truth["tso_positive_reads"])),
                    "raw_structure_contaminant_reads": int(
                        raw.get("structure_contaminant_reads", truth["structure_contaminants"])
                    ),
                    "structure_contaminants": truth["structure_contaminants"],
                    "contamination_fraction": truth["fraction"],
                    "raw_contaminant_fraction_of_all": raw.get(
                        "contaminant_fraction_of_all_reads", truth["fraction"]
                    ),
                    "raw_contaminant_fraction_of_tso": raw.get(
                        "contaminant_fraction_of_tso_reads", float("nan")
                    ),
                    "wilson_95_low": truth["wilson_95_low"],
                    "wilson_95_high": truth["wilson_95_high"],
                    "known_structure_contaminants": truth["known_structure_contaminants"],
                    "complete_reference_contaminants": truth["complete_reference_contaminants"],
                    "inclusive_possible_contaminants": truth.get(
                        "inclusive_contaminants", float("nan")
                    ),
                    "inclusive_possible_contaminant_fraction_of_all": (
                        truth.get("inclusive_contaminants", 0)
                        / summary["records_evaluated"]
                        if summary["records_evaluated"]
                        else float("nan")
                    ),
                    "valid_barcode_reads": valid_all,
                    "valid_barcode_tso_reads": valid_tso,
                    "valid_barcode_structure_contaminant_reads": valid_contaminants,
                    "valid_barcode_read_contaminant_fraction_of_all": (
                        valid_contaminants / valid_all if valid_all else float("nan")
                    ),
                    "valid_barcode_read_contaminant_fraction_of_tso": (
                        valid_contaminants / valid_tso if valid_tso else float("nan")
                    ),
                    "cell_pb_umi_all": lower.get("all", float("nan")),
                    "cell_pb_umi_tso": lower.get("tso", float("nan")),
                    "cell_pb_umi_structure_contaminants": lower.get(
                        "structure_contaminant", float("nan")
                    ),
                    "cell_pb_umi_contaminant_fraction_of_all": lower_fractions.get(
                        "contaminant_of_all", {}
                    ).get("fraction", float("nan")),
                    "cell_pb_umi_contaminant_fraction_of_tso": lower_fractions.get(
                        "contaminant_of_tso", {}
                    ).get("fraction", float("nan")),
                    "cell_pb_umi_inclusive_contaminant_fraction_of_all": lower_fractions.get(
                        "inclusive_contaminant_of_all", {}
                    ).get("fraction", float("nan")),
                    "cell_pb_umi_inclusive_contaminant_fraction_of_tso": lower_fractions.get(
                        "inclusive_contaminant_of_tso", {}
                    ).get("fraction", float("nan")),
                    "cell_pb_umi_insert_all": upper.get("all", float("nan")),
                    "cell_pb_umi_insert_tso": upper.get("tso", float("nan")),
                    "cell_pb_umi_insert_structure_contaminants": upper.get(
                        "structure_contaminant", float("nan")
                    ),
                    "cell_pb_umi_insert_contaminant_fraction_of_all": upper_fractions.get(
                        "contaminant_of_all", {}
                    ).get("fraction", float("nan")),
                    "cell_pb_umi_insert_contaminant_fraction_of_tso": upper_fractions.get(
                        "contaminant_of_tso", {}
                    ).get("fraction", float("nan")),
                    "cell_pb_umi_insert_inclusive_contaminant_fraction_of_all": upper_fractions.get(
                        "inclusive_contaminant_of_all", {}
                    ).get("fraction", float("nan")),
                    "cell_pb_umi_insert_inclusive_contaminant_fraction_of_tso": upper_fractions.get(
                        "inclusive_contaminant_of_tso", {}
                    ).get("fraction", float("nan")),
                    "selected_filter": summary["selected_filter"],
                    "recommended_filter": summary["recommended_filter_from_benchmark"],
                }
            )
        elif os.path.exists(benchmark_path):
            print(
                f"  WARNING: ignoring {name} benchmark summary: "
                f"{benchmark_reason}"
            )
        metric_path = os.path.join(lib_dir, f"short{read_length}_filter_metrics.csv")
        if benchmark_current and os.path.exists(metric_path):
            metric_frames.append(pd.read_csv(metric_path))
        stratum_path = os.path.join(
            lib_dir, f"short{read_length}_metrics_by_tso_offset.csv"
        )
        if benchmark_current and os.path.exists(stratum_path):
            stratum_frames.append(pd.read_csv(stratum_path))
        bam_path = os.path.join(lib_dir, "bam_audit_summary.yaml")
        bam_current, bam_status_reason = bam_audit_status(
            bam_path,
            current_module_digest,
            benchmark_config_sha256(cfg, name),
        )
        if bam_current:
            with open(bam_path) as handle:
                bam_summary = yaml.safe_load(handle)
            counted = bam_summary.get("counted_molecule_impact", {})
            bam_row = {
                key: bam_summary.get(key)
                for key in (
                    "library",
                    "primary_reads",
                    "mapped_primary_reads",
                    "structure_contaminants",
                    "mapped_structure_contaminants",
                    "unique_contaminants",
                    "gene_assigned_contaminants",
                    "assigned_to_guide_target",
                    "gene_assigned_contaminant_fraction_of_all_gene_assigned_reads",
                )
            }
            bam_row.update(
                {f"molecule_{key}": value for key, value in counted.items()}
            )
            bam_rows.append(bam_row)
        elif os.path.exists(bam_path):
            print(
                f"  WARNING: ignoring {name} BAM audit: {bam_status_reason}"
            )
        gene_path = os.path.join(lib_dir, "bam_contaminant_gene_assignments.csv")
        if bam_current and os.path.exists(gene_path):
            frame = pd.read_csv(gene_path)
            frame.insert(0, "library", name)
            gene_frames.append(frame)

        for variant in variants:
            impact_path = os.path.join(
                lib_dir,
                alignment_cfg.get("output_subdir", "alignment_star"),
                variant,
                "alignment_impact_summary.yaml",
            )
            if not os.path.exists(impact_path):
                continue
            impact_current, stale_reason = alignment_impact_status(
                impact_path,
                benchmark_by_library.get(name, {}),
                current_module_digest,
                benchmark_config_sha256(cfg, name),
            )
            if not impact_current:
                print(
                    f"  WARNING: ignoring {name}/{variant} alignment impact: "
                    f"{stale_reason}"
                )
                continue
            with open(impact_path) as handle:
                impact = yaml.safe_load(handle)
            hll = impact.get("estimated_exact_key_molecules_hll", {})
            exact = impact.get("splitpipe_equivalent_hamming1_molecule_impact", {})
            raw_row = next(
                (row for row in summary_rows if row["library"] == name), {}
            )
            conditional_molecule = (
                exact.get(
                    "false_gene_assignment_fraction_of_input_molecules",
                    float("nan"),
                )
                if exact.get("enabled", False)
                else impact.get(
                    "false_gene_assignment_fraction_of_input_molecules",
                    float("nan"),
                )
            )
            trimmed_counts = impact.get("predicted_tso_trim_input_counts", {})
            trim_outcomes = impact.get("outcomes_by_predicted_tso_trim", {})
            trimmed_input = int(trimmed_counts.get("True", 0))
            untrimmed_input = int(trimmed_counts.get("False", 0))
            trimmed_false = int(
                trim_outcomes.get("True", {}).get("false_gene_assigned_reads", 0)
            )
            untrimmed_false = int(
                trim_outcomes.get("False", {}).get("false_gene_assigned_reads", 0)
            )
            false_read_ci = impact.get(
                "false_gene_assignment_fraction_of_input_reads_cluster_bootstrap_95",
                [float("nan"), float("nan")],
            )
            naive_read_ci = impact.get(
                "false_gene_assignment_fraction_of_input_reads_naive_wilson_95",
                [float("nan"), float("nan")],
            )
            molecule_cluster_ci = exact.get(
                "false_gene_assignment_fraction_of_input_molecules_cluster_bootstrap_95",
                [float("nan"), float("nan")],
            )
            alignment_rows.append(
                {
                    "library": name,
                    "variant": variant,
                    "input_reads": impact.get("input_reads", 0),
                    "mapped_primary_reads": impact.get("mapped_primary_reads", 0),
                    "unique_reads": impact.get("unique_reads", 0),
                    "multimapping_primary_reads": impact.get(
                        "multimapping_primary_reads", 0
                    ),
                    "false_gene_assigned_reads": impact.get(
                        "false_gene_assigned_reads", 0
                    ),
                    "guide_target_assigned_reads": impact.get(
                        "guide_target_assigned_reads", 0
                    ),
                    "false_gene_assignment_fraction_of_input_reads": impact.get(
                        "false_gene_assignment_fraction_of_input_reads", float("nan")
                    ),
                    "false_gene_assignment_fraction_of_input_reads_cluster_bootstrap_95_low": false_read_ci[0],
                    "false_gene_assignment_fraction_of_input_reads_cluster_bootstrap_95_high": false_read_ci[1],
                    "false_gene_assignment_fraction_of_input_reads_naive_wilson_95_low": naive_read_ci[0],
                    "false_gene_assignment_fraction_of_input_reads_naive_wilson_95_high": naive_read_ci[1],
                    "guide_target_assignment_fraction_of_input_reads": impact.get(
                        "guide_target_assignment_fraction_of_input_reads", float("nan")
                    ),
                    "tso_trimmed_input_reads": trimmed_input,
                    "tso_trimmed_false_gene_assignment_fraction": (
                        trimmed_false / trimmed_input
                        if trimmed_input
                        else float("nan")
                    ),
                    "tso_untrimmed_input_reads": untrimmed_input,
                    "tso_untrimmed_false_gene_assignment_fraction": (
                        untrimmed_false / untrimmed_input
                        if untrimmed_input
                        else float("nan")
                    ),
                    "hll_input_molecules": hll.get("input", float("nan")),
                    "hll_false_gene_assigned_molecules": hll.get(
                        "gene_assigned", float("nan")
                    ),
                    "hamming1_input_molecules": exact.get(
                        "ground_truth_input_molecules_by_guide", float("nan")
                    ),
                    "hamming1_false_gene_assigned_molecules": exact.get(
                        "false_gene_assigned_molecules", float("nan")
                    ),
                    "hamming1_false_gene_assignment_fraction_of_input_molecules": conditional_molecule,
                    "hamming1_false_gene_assignment_fraction_of_input_molecules_cluster_bootstrap_95_low": molecule_cluster_ci[0],
                    "hamming1_false_gene_assignment_fraction_of_input_molecules_cluster_bootstrap_95_high": molecule_cluster_ci[1],
                    "hamming1_guide_target_molecules_any_read": exact.get(
                        "guide_target_molecules_any_read", float("nan")
                    ),
                    "projected_false_gene_assigned_fraction_of_valid_reads": (
                        raw_row.get(
                            "valid_barcode_read_contaminant_fraction_of_all",
                            float("nan"),
                        )
                        * impact.get(
                            "false_gene_assignment_fraction_of_input_reads",
                            float("nan"),
                        )
                    ),
                    "projected_cell_pb_umi_false_molecule_fraction": (
                        raw_row.get(
                            "cell_pb_umi_contaminant_fraction_of_all", float("nan")
                        )
                        * conditional_molecule
                    ),
                    "projected_cell_pb_umi_insert_false_molecule_fraction": (
                        raw_row.get(
                            "cell_pb_umi_insert_contaminant_fraction_of_all",
                            float("nan"),
                        )
                        * conditional_molecule
                    ),
                }
            )

    for library in cfg.get("libraries", []):
        name = library["sublibrary_name"]
        lib_dir = os.path.join(output_root, name)
        historical_path = os.path.join(
            lib_dir, "bam_historical_gene_molecule_totals.csv"
        )
        if not os.path.exists(historical_path):
            continue
        historical = pd.read_csv(historical_path)
        library_burden = []
        for variant in variants:
            controlled_path = os.path.join(
                lib_dir,
                alignment_cfg.get("output_subdir", "alignment_star"),
                variant,
                "false_gene_assigned_molecules_by_gene.csv",
            )
            if not os.path.exists(controlled_path):
                continue
            controlled = pd.read_csv(controlled_path)
            if controlled.empty:
                continue
            controlled = controlled.groupby(
                ["gene_id", "gene_name", "production_filter_stratum"],
                as_index=False,
            )["false_gene_assigned_molecules"].sum()
            joined = controlled.merge(
                historical,
                on="gene_id",
                how="left",
                suffixes=("_controlled", "_historical"),
            )
            joined.insert(0, "variant", variant)
            joined.insert(0, "library", name)
            joined["historical_gene_molecules"] = joined[
                "historical_gene_molecules"
            ].fillna(0).astype(int)
            denominator = joined["historical_gene_molecules"]
            joined["controlled_false_fraction_of_historical_gene_molecules"] = (
                joined["false_gene_assigned_molecules"]
                / denominator.where(denominator > 0)
            )
            joined["interpretation"] = (
                "projected_counterfactual_burden"
                if variant == f"short{read_length}"
                else "controlled_full_length_potential_burden"
            )
            library_burden.append(joined)
        if library_burden:
            burden = pd.concat(library_burden, ignore_index=True)
            burden.to_csv(
                os.path.join(lib_dir, "gene_expression_burden.csv"), index=False
            )
            burden_frames.append(burden)

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(
        os.path.join(output_root, "all_libraries_contamination_summary.csv"),
        index=False,
    )
    metrics_df = pd.concat(metric_frames, ignore_index=True) if metric_frames else pd.DataFrame()
    if metric_frames:
        metrics_df.to_csv(
            os.path.join(
                output_root, f"all_libraries_short{read_length}_filter_metrics.csv"
            ),
            index=False,
        )
    if stratum_frames:
        pd.concat(stratum_frames, ignore_index=True).to_csv(
            os.path.join(
                output_root,
                f"all_libraries_short{read_length}_metrics_by_tso_offset.csv",
            ),
            index=False,
        )
    if bam_rows:
        pd.DataFrame(bam_rows).to_csv(
            os.path.join(output_root, "all_libraries_bam_audit_summary.csv"), index=False
        )
    if gene_frames:
        pd.concat(gene_frames, ignore_index=True).to_csv(
            os.path.join(output_root, "all_libraries_bam_gene_assignments.csv"), index=False
        )
    if guide_rows:
        pd.DataFrame(guide_rows).to_csv(
            os.path.join(output_root, "all_libraries_guide_assignment_impact.csv"),
            index=False,
        )
    if burden_frames:
        pd.concat(burden_frames, ignore_index=True).to_csv(
            os.path.join(output_root, "all_libraries_gene_expression_burden.csv"),
            index=False,
        )

    alignment_df = pd.DataFrame(alignment_rows)
    if not alignment_df.empty:
        alignment_df.to_csv(
            os.path.join(output_root, "all_libraries_alignment_impact.csv"),
            index=False,
        )
    delta_rows = []
    short_variant = f"short{read_length}"
    if not alignment_df.empty:
        for name, frame in alignment_df.groupby("library"):
            full = frame[frame["variant"] == "full"]
            short = frame[frame["variant"] == short_variant]
            if full.empty or short.empty:
                continue
            full_row = full.iloc[0]
            short_row = short.iloc[0]
            full_read = float(full_row["false_gene_assignment_fraction_of_input_reads"])
            short_read = float(short_row["false_gene_assignment_fraction_of_input_reads"])
            full_molecule = float(
                full_row["hamming1_false_gene_assignment_fraction_of_input_molecules"]
            )
            short_molecule = float(
                short_row["hamming1_false_gene_assignment_fraction_of_input_molecules"]
            )
            read_change = {}
            molecule_change = {}
            full_cluster_path = os.path.join(
                output_root,
                name,
                alignment_cfg.get("output_subdir", "alignment_star"),
                "full",
                "cb_pb_cluster_alignment_outcomes.csv",
            )
            short_cluster_path = os.path.join(
                output_root,
                name,
                alignment_cfg.get("output_subdir", "alignment_star"),
                short_variant,
                "cb_pb_cluster_alignment_outcomes.csv",
            )
            if os.path.exists(full_cluster_path) and os.path.exists(short_cluster_path):
                full_clusters = pd.read_csv(full_cluster_path)
                short_clusters = pd.read_csv(short_cluster_path)
                replicates = int(
                    alignment_cfg.get("cluster_bootstrap_replicates", 1000)
                )
                read_change = paired_cluster_bootstrap_change(
                    full_clusters,
                    short_clusters,
                    "false_gene_assigned_reads",
                    "input_reads",
                    replicates,
                    f"{name}|full_vs_{short_variant}|reads",
                )
                molecule_change = paired_cluster_bootstrap_change(
                    full_clusters,
                    short_clusters,
                    "false_gene_assigned_molecules",
                    "input_molecules",
                    replicates,
                    f"{name}|full_vs_{short_variant}|molecules",
                )
            delta_row = {
                    "library": name,
                    "full_read_false_assignment_fraction": full_read,
                    f"{short_variant}_read_false_assignment_fraction": short_read,
                    "read_false_assignment_absolute_delta": short_read - full_read,
                    "read_false_assignment_fold_change": (
                        short_read / full_read if full_read > 0 else float("nan")
                    ),
                    "full_hamming1_molecule_false_assignment_fraction": full_molecule,
                    f"{short_variant}_hamming1_molecule_false_assignment_fraction": short_molecule,
                    "hamming1_molecule_absolute_delta": short_molecule - full_molecule,
                    "hamming1_molecule_fold_change": (
                        short_molecule / full_molecule
                        if full_molecule > 0
                        else float("nan")
                    ),
                    "full_projected_false_read_fraction": float(
                        full_row["projected_false_gene_assigned_fraction_of_valid_reads"]
                    ),
                    f"{short_variant}_projected_false_read_fraction": float(
                        short_row["projected_false_gene_assigned_fraction_of_valid_reads"]
                    ),
                }
            delta_row.update(
                {f"read_{key}": value for key, value in read_change.items()}
            )
            delta_row.update(
                {f"hamming1_molecule_{key}": value for key, value in molecule_change.items()}
            )
            delta_rows.append(delta_row)
    delta_df = pd.DataFrame(delta_rows)
    if not delta_df.empty:
        delta_df.to_csv(
            os.path.join(
                output_root,
                f"all_libraries_full_vs_{short_variant}_alignment_delta.csv",
            ),
            index=False,
        )
    completeness_rows = []
    bam_cfg = qc.get("bam_audit", {}) or {}
    bam_expected = bool(bam_cfg.get("enabled", True))
    bam_required = bool(
        bam_cfg.get("required_for_audit", bam_expected)
    )
    alignment_expected = bool(alignment_cfg.get("enabled", True))
    guide_cfg = qc.get("guide_assignment_impact", {}) or {}
    guide_expected = bool(guide_cfg.get("enabled", False))
    for library in cfg.get("libraries", []):
        name = library["sublibrary_name"]
        lib_dir = os.path.join(output_root, name)
        benchmark_path = os.path.join(lib_dir, "benchmark_summary.yaml")
        benchmark_complete, benchmark_reason = benchmark_summary_status(
            benchmark_path, cfg, name, current_module_digest
        )
        cap_reached = False
        truth_reads = 0
        benchmark_summary = {}
        if benchmark_complete:
            with open(benchmark_path) as handle:
                benchmark_summary = yaml.safe_load(handle) or {}
            truth_input = benchmark_summary.get("alignment_ground_truth_input", {})
            cap_reached = bool(truth_input.get("safety_cap_reached", False))
            truth_reads = int(truth_input.get("reads_written_per_variant", 0))
        bam_summary_path = os.path.join(lib_dir, "bam_audit_summary.yaml")
        bam_complete, bam_reason = bam_audit_status(
            bam_summary_path,
            current_module_digest,
            benchmark_config_sha256(cfg, name),
        )
        row = {
            "library": name,
            "benchmark_complete": benchmark_complete,
            "benchmark_status": benchmark_reason,
            "truth_safety_cap_reached": cap_reached,
            "alignment_truth_reads_per_variant": truth_reads,
            "bam_audit_expected": bam_expected,
            "bam_audit_required": bam_required,
            "bam_audit_complete": bam_complete,
            "bam_audit_status": bam_reason,
            "controlled_alignment_expected": alignment_expected,
        }
        guide_summary_path = os.path.join(
            lib_dir,
            guide_cfg.get("output_subdir", "guide_assignment_impact"),
            "guide_assignment_impact_summary.yaml",
        )
        guide_complete, guide_reason = guide_impact_status(
            guide_summary_path,
            current_module_digest,
            benchmark_config_sha256(cfg, name),
        )
        row["guide_assignment_impact_expected"] = guide_expected
        row["guide_assignment_impact_complete"] = guide_complete
        row["guide_assignment_impact_status"] = guide_reason
        alignment_complete = True
        for variant in variants:
            impact_path = os.path.join(
                lib_dir,
                alignment_cfg.get("output_subdir", "alignment_star"),
                variant,
                "alignment_impact_summary.yaml",
            )
            complete, reason = alignment_impact_status(
                impact_path,
                benchmark_summary,
                current_module_digest,
                benchmark_config_sha256(cfg, name),
            )
            row[f"alignment_{variant}_complete"] = complete
            row[f"alignment_{variant}_status"] = reason
            alignment_complete = alignment_complete and complete
        row["audit_complete"] = bool(
            benchmark_complete
            and not cap_reached
            and (row["bam_audit_complete"] or not bam_required)
            and (alignment_complete or not alignment_expected)
            and (guide_complete or not guide_expected)
        )
        completeness_rows.append(row)
    completeness_df = pd.DataFrame(completeness_rows)
    completeness_df.to_csv(
        os.path.join(output_root, "audit_completeness.csv"), index=False
    )
    report_path = write_aggregate_audit_report(
        summary_df,
        alignment_df,
        delta_df,
        completeness_df,
        output_root,
        read_length,
    )
    if qc.get("plots", {}).get("enabled", True):
        make_aggregate_plots(summary_df, metrics_df, alignment_df, output_root)
    print(f"  Aggregated outputs under {output_root}")
    print(f"  Wrote audit report: {report_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="QC YAML extending config.yaml")
    parser.add_argument("--library", help="Sublibrary name (required except aggregate)")
    parser.add_argument(
        "--mode",
        choices=["benchmark", "bam", "all", "filter", "alignment", "guide", "aggregate"],
        default="all",
    )
    parser.add_argument(
        "--alignment-variant",
        help="Controlled STAR arm (full or short<length>), required for --mode alignment",
    )
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument(
        "--execution-backend",
        choices=["auto", "serial", "parallel"],
        default="auto",
        help="Raw FASTQ benchmark execution backend (auto uses workers when --threads >1)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        help="Records per bounded worker batch (default: performance.batch_size or 20000)",
    )
    args = parser.parse_args()

    cfg = load_merged_config(args.config)
    if args.mode == "aggregate":
        aggregate(cfg, args.config)
        return
    if not args.library:
        parser.error("--library is required unless --mode aggregate")

    performance_cfg = cfg["gex_grna_contamination_qc"].get("performance", {}) or {}
    batch_size = (
        args.batch_size
        if args.batch_size is not None
        else int(performance_cfg.get("batch_size", 20_000))
    )

    if args.mode in ("benchmark", "all"):
        benchmark_library(
            cfg,
            args.config,
            args.library,
            args.threads,
            execution_backend=args.execution_backend,
            batch_size=batch_size,
        )
    if args.mode in ("bam", "all"):
        if cfg["gex_grna_contamination_qc"].get("bam_audit", {}).get("enabled", True):
            audit_bam(cfg, args.config, args.library, args.threads)
    if args.mode in ("guide", "all"):
        if cfg["gex_grna_contamination_qc"].get(
            "guide_assignment_impact", {}
        ).get("enabled", False):
            audit_guide_assignment_impact(
                cfg, args.config, args.library, args.threads
            )
    if args.mode == "filter":
        filter_library(cfg, args.config, args.library)
    if args.mode == "alignment":
        if not args.alignment_variant:
            parser.error("--alignment-variant is required for --mode alignment")
        audit_alignment_benchmark(
            cfg,
            args.config,
            args.library,
            args.alignment_variant,
            args.threads,
        )


if __name__ == "__main__":
    main()

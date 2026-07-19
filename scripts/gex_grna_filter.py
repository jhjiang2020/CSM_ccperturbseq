#!/usr/bin/env python3
"""Removal-only gRNA classifier for GEX reads.

This module implements the packed inline matcher used by the optional Step 1a
GEX decontamination hook.  It deliberately answers only one biological
question: is a TSO-bearing read supported by a guide sequence in the supplied
reference under the removal policy?  It does *not* assign guide names, iBARs,
constructs, pairs, targets, cells, or UMIs.  Those calls belong to Step 1c and
the ordinary Steps 2--4.

The guide index uses two bits per base, an exact packed lookup, and two
complementary non-overlapping seed masks for each guide length.  If a query and
reference differ by at most one substitution, at least one seed is unchanged;
the union of the two seed buckets is therefore a lossless candidate set for
Hamming distance <= 1.  Full packed-guide verification prevents seed hits from
becoming false-positive matches.

The production TSO route uses edlib's C-backed bit-parallel matcher.  The
dependency-free Myers implementation remains the executable correctness
reference used by tests and by the rare constrained-location fallback.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple

import yaml

try:
    import edlib

    HAS_EDLIB = True
except ImportError:  # pragma: no cover - exercised by dependency-gate tests
    edlib = None
    HAS_EDLIB = False


DNA_BASES = frozenset("ACGT")
BASE_TO_BITS = {"A": 0, "C": 1, "G": 2, "T": 3}


def _validate_dna(sequence: str, label: str, *, allow_empty: bool = False) -> str:
    """Validate and return an uppercase DNA sequence.

    Reference and configuration sequences are intentionally strict.  Read
    queries are handled separately and an N or other non-ACGT symbol simply
    produces a negative diagnostic rather than terminating preprocessing.
    """

    sequence = str(sequence).strip()
    if not sequence and allow_empty:
        return sequence
    if not sequence:
        raise ValueError(f"{label} must not be empty")
    if sequence != sequence.upper() or any(base not in DNA_BASES for base in sequence):
        raise ValueError(f"{label} must contain uppercase A/C/G/T only: {sequence!r}")
    return sequence


def encode_dna(sequence: str) -> int:
    """Pack an uppercase A/C/G/T sequence into two bits per base."""

    packed = 0
    for base in sequence:
        try:
            code = BASE_TO_BITS[base]
        except KeyError as exc:
            raise ValueError(
                f"DNA sequence must contain uppercase A/C/G/T only: {sequence!r}"
            ) from exc
        packed = (packed << 2) | code
    return packed


def packed_hamming_distance(left: int, right: int, length: int) -> int:
    """Return base-wise Hamming distance between equal-length packed guides."""

    different = left ^ right
    # Fold each two-bit base into its low bit, then count non-zero base pairs.
    different |= different >> 1
    # (2 ** (2L) - 1) / 3 is binary 0101...01 (L one-bits).
    pair_low_bits = ((1 << (2 * length)) - 1) // 3
    return (different & pair_low_bits).bit_count()


def _bounded_levenshtein(left: str, right: str, max_distance: int) -> int:
    """Dependency-free Levenshtein distance, capped at max_distance + 1."""

    if abs(len(left) - len(right)) > max_distance:
        return max_distance + 1
    if len(left) < len(right):
        left, right = right, left
    previous = list(range(len(right) + 1))
    for row_index, left_base in enumerate(left, start=1):
        current = [row_index]
        row_minimum = row_index
        for column_index, right_base in enumerate(right, start=1):
            value = min(
                previous[column_index] + 1,
                current[column_index - 1] + 1,
                previous[column_index - 1] + (left_base != right_base),
            )
            current.append(value)
            row_minimum = min(row_minimum, value)
        if row_minimum > max_distance:
            return max_distance + 1
        previous = current
    return min(previous[-1], max_distance + 1)


def _pattern_masks(pattern: str) -> Dict[str, int]:
    """Return the per-base bit masks used by the single-word Myers core."""

    masks: Dict[str, int] = {}
    for index, base in enumerate(pattern):
        masks[base] = masks.get(base, 0) | (1 << index)
    return masks


def _myers_prefix_distances(
    pattern: str,
    text: str,
    masks: Optional[Dict[str, int]] = None,
) -> List[int]:
    """Return global edit distance for every non-empty prefix of ``text``.

    The configured TSO is short enough for one machine word.  Traversing the
    Myers state once per possible start replaces a separate dynamic-programming
    matrix for every overlapping candidate window.
    """

    if not pattern:
        return list(range(1, len(text) + 1))
    if len(pattern) > 63:
        raise ValueError("single-word Myers implementation supports <=63 bases")
    masks = _pattern_masks(pattern) if masks is None else masks
    width_mask = (1 << len(pattern)) - 1
    highest_bit = 1 << (len(pattern) - 1)
    positive = width_mask
    negative = 0
    score = len(pattern)
    distances: List[int] = []
    for base in text:
        equal = masks.get(base, 0)
        xv = equal | negative
        xh = (((equal & positive) + positive) ^ positive) | equal
        ph = negative | ~(xh | positive)
        mh = positive & xh
        if ph & highest_bit:
            score += 1
        elif mh & highest_bit:
            score -= 1
        ph = ((ph << 1) | 1) & width_mask
        mh = (mh << 1) & width_mask
        positive = (mh | ~(xv | ph)) & width_mask
        negative = ph & xv
        distances.append(score)
    return distances


@dataclass(frozen=True)
class TsoMatch:
    """A full-motif TSO match at an allowed read-start offset."""

    found: bool
    start: int = -1
    end: int = -1
    distance: int = -1
    observed_length: int = 0
    match_type: str = "none"


def find_tso(
    r1_sequence: str,
    tso_sequence: str,
    max_edit_distance: int = 3,
    max_start: int = 2,
) -> TsoMatch:
    """Find a full TSO motif starting within ``0..max_start``.

    Exact matches are selected before fuzzy matches.  Fuzzy matching aligns
    the complete configured TSO against windows that are up to
    ``max_edit_distance`` bases shorter or longer, allowing substitutions and
    small indels while preserving the correct post-TSO extraction boundary.
    The TSO is never searched later in the read than ``max_start``.
    """

    if max_edit_distance < 0:
        raise ValueError("max_edit_distance must be >= 0")
    if max_start < 0:
        raise ValueError("max_start must be >= 0")
    tso_sequence = _validate_dna(tso_sequence, "TSO sequence")
    read = str(r1_sequence).strip().upper()
    tso_length = len(tso_sequence)

    # Search a small guard band beyond the accepted boundary.  Without this,
    # an exact TSO at max_start + 1 could be reinterpreted as an allowed fuzzy
    # match with one leading insertion.  We locate the best motif first and
    # apply the configured start boundary to that location.
    search_max_start = max_start + max_edit_distance

    # Preserve exact-first semantics even if a fuzzy candidate starts earlier.
    for start in range(search_max_start + 1):
        end = start + tso_length
        if end <= len(read) and read[start:end] == tso_sequence:
            if start <= max_start:
                return TsoMatch(True, start, end, 0, tso_length, "exact")
            return TsoMatch(False)

    if max_edit_distance == 0:
        return TsoMatch(False)

    # Pigeonhole prefilter: an alignment with at most k edits must preserve at
    # least one of k+1 non-overlapping query blocks as an exact target
    # substring.  Most non-TSO GEX reads fail these C-level substring searches,
    # avoiding dozens of small dynamic-programming matrices per read.
    block_count = max_edit_distance + 1
    target_limit = min(
        len(read), search_max_start + tso_length + max_edit_distance
    )
    target_prefix = read[:target_limit]
    has_exact_block = False
    for block_number in range(block_count):
        block_start = (block_number * tso_length) // block_count
        block_end = ((block_number + 1) * tso_length) // block_count
        if tso_sequence[block_start:block_end] in target_prefix:
            has_exact_block = True
            break
    if not has_exact_block:
        return TsoMatch(False)

    masks = _pattern_masks(tso_sequence)
    best_candidate: Optional[Tuple[int, int, int, int]] = None
    minimum_window = max(1, tso_length - max_edit_distance)
    maximum_window = tso_length + max_edit_distance
    for start in range(search_max_start + 1):
        prefix_distances = _myers_prefix_distances(
            tso_sequence,
            read[start : start + maximum_window],
            masks,
        )
        for observed_length in range(minimum_window, maximum_window + 1):
            end = start + observed_length
            if observed_length > len(prefix_distances):
                continue
            distance = prefix_distances[observed_length - 1]
            if distance <= max_edit_distance:
                candidate = (
                    distance,
                    start,
                    abs(observed_length - tso_length),
                    end,
                )
                if best_candidate is None or candidate < best_candidate:
                    best_candidate = candidate

    if best_candidate is None:
        return TsoMatch(False)
    distance, start, _, end = best_candidate
    if start > max_start:
        return TsoMatch(False)
    return TsoMatch(True, start, end, distance, end - start, "fuzzy")


class TsoMatcher:
    """Precomputed production TSO matcher with exact Myers-equivalent semantics."""

    def __init__(
        self,
        tso_sequence: str,
        max_edit_distance: int = 3,
        max_start: int = 2,
        backend: str = "edlib",
    ) -> None:
        if max_edit_distance < 0:
            raise ValueError("max_edit_distance must be >= 0")
        if max_start < 0:
            raise ValueError("max_start must be >= 0")
        backend = str(backend).strip().lower()
        if backend not in {"edlib", "myers"}:
            raise ValueError("TSO backend must be edlib or myers")
        if backend == "edlib" and not HAS_EDLIB:
            raise RuntimeError(
                "gex_decontamination TSO backend edlib requires edlib>=1.3; "
                "install requirements_extra.txt in the splitpipe environment"
            )

        self.tso_sequence = _validate_dna(tso_sequence, "TSO sequence")
        self.max_edit_distance = int(max_edit_distance)
        self.max_start = int(max_start)
        self.backend = backend
        self.tso_length = len(self.tso_sequence)
        self.search_max_start = self.max_start + self.max_edit_distance
        self.minimum_window = max(1, self.tso_length - self.max_edit_distance)
        self.maximum_window = self.tso_length + self.max_edit_distance
        self.target_limit = self.search_max_start + self.maximum_window

        block_count = self.max_edit_distance + 1
        self.blocks = tuple(
            self.tso_sequence[
                (block_number * self.tso_length) // block_count :
                ((block_number + 1) * self.tso_length) // block_count
            ]
            for block_number in range(block_count)
        )

    def find(self, r1_sequence: str) -> TsoMatch:
        """Normalize and match one R1 sequence."""

        return self.find_normalized(str(r1_sequence).strip().upper())

    def find_normalized(self, read: str) -> TsoMatch:
        """Match an already normalized R1 sequence without per-read setup."""

        exact_stop = min(
            len(read), self.search_max_start + self.tso_length
        )
        exact_start = read.find(self.tso_sequence, 0, exact_stop)
        if exact_start >= 0:
            if exact_start <= self.max_start:
                return TsoMatch(
                    True,
                    exact_start,
                    exact_start + self.tso_length,
                    0,
                    self.tso_length,
                    "exact",
                )
            return TsoMatch(False)

        if self.max_edit_distance == 0:
            return TsoMatch(False)

        target_stop = min(len(read), self.target_limit)
        if not any(read.find(block, 0, target_stop) >= 0 for block in self.blocks):
            return TsoMatch(False)

        if self.backend == "myers":
            return find_tso(
                read,
                self.tso_sequence,
                self.max_edit_distance,
                self.max_start,
            )

        target = read[:target_stop]
        result = edlib.align(  # type: ignore[union-attr]
            self.tso_sequence,
            target,
            mode="HW",
            k=self.max_edit_distance,
            task="locations",
        )
        distance = int(result.get("editDistance", -1))
        if distance < 0:
            return TsoMatch(False)

        candidates = []
        for raw_start, raw_end in result.get("locations", ()):
            if raw_start is None or raw_end is None:
                continue
            start = int(raw_start)
            end = int(raw_end) + 1
            observed_length = end - start
            if not 0 <= start <= self.search_max_start:
                continue
            if not self.minimum_window <= observed_length <= self.maximum_window:
                continue
            candidates.append(
                (distance, start, abs(observed_length - self.tso_length), end)
            )

        if not candidates:
            # HW can find a better target substring just outside the permitted
            # start range.  In that rare case it must not hide a valid but
            # higher-distance constrained candidate.
            return find_tso(
                read,
                self.tso_sequence,
                self.max_edit_distance,
                self.max_start,
            )

        distance, start, _, end = min(candidates)
        if start > self.max_start:
            return TsoMatch(False)
        return TsoMatch(
            True,
            start,
            end,
            distance,
            end - start,
            "fuzzy",
        )


@dataclass(frozen=True)
class SequenceGroup:
    """Anonymous reference-sequence group; contains no biological identity."""

    group_id: str
    length: int
    sequence: str
    packed_sequence: int
    reference_row_count: int


@dataclass(frozen=True)
class LengthIndexStats:
    """Immutable construction statistics for one guide length."""

    guide_length: int
    sequence_groups: int
    reference_rows: int
    seed_masks: Tuple[int, int]
    seed_position_counts: Tuple[int, int]
    seed_key_counts: Tuple[int, int]
    seed_max_bucket_sizes: Tuple[int, int]
    seed_bucket_histograms: Tuple[Tuple[Tuple[int, int], ...], Tuple[Tuple[int, int], ...]]


@dataclass(frozen=True)
class IndexStats:
    """Immutable manifest-ready summary of a packed reference index."""

    guide_lengths: Tuple[int, ...]
    reference_rows: int
    unique_sequence_groups: int
    duplicate_reference_rows: int
    reference_sha256: str
    by_length: Tuple[LengthIndexStats, ...]


@dataclass(frozen=True)
class GuideMatch:
    """Removal evidence for one complete TSO-relative guide-length query."""

    query_length: int
    supported: bool
    distance: Optional[int]
    sequence_group_ids: Tuple[str, ...]
    reference_row_count: int
    sequence_ambiguous: bool
    shared_reference_sequence: bool
    candidate_count: int
    reason: str


@dataclass(frozen=True)
class FilterSettings:
    """Validated removal settings inherited from the main config."""

    tso_sequence: str
    tso_max_edit_distance: int
    max_tso_start: int
    guide_lengths: Tuple[int, ...]
    max_guide_mismatches: int
    remove_ambiguous_hits: bool
    max_seed_bucket_size: int
    reference_path: str
    tso_backend: str = "myers"


def filter_settings_sha256(settings: FilterSettings) -> str:
    """Return a stable fingerprint for every removal-policy setting.

    The guide-reference content has its own checksum in :class:`IndexStats`;
    this fingerprint prevents split-pipe from reusing a preprocessed FASTQ
    after changing TSO tolerance, guide mismatch policy, ambiguity handling,
    or any other classifier setting.
    """

    payload = json.dumps(asdict(settings), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class FilterDecision:
    """Removal-only decision emitted by Step 1a.

    ``sequence_group_ids`` are stable hashes of reference sequences.  They are
    useful for ambiguity audits but are not guide or construct assignments.
    """

    read_id: str
    raw_pair_location: str
    tso_found: bool
    tso_start: int
    tso_end: int
    tso_distance: int
    tso_match_type: str
    reference_supported: bool
    filter_positive: bool
    filter_reason: str
    confidence_tier: str
    guide_lengths: Tuple[int, ...]
    guide_distance: Optional[int]
    sequence_group_ids: Tuple[str, ...]
    sequence_ambiguous: bool
    shared_reference_sequence: bool
    candidate_count: int
    length_matches: Tuple[GuideMatch, ...]

    @property
    def guide_length(self) -> Optional[int]:
        """Return the sole supported guide length, or None for no hit/a tie."""

        return self.guide_lengths[0] if len(self.guide_lengths) == 1 else None

    def to_dict(self) -> dict:
        """Return a JSON/TSV-friendly diagnostic representation."""

        return asdict(self)


@dataclass(frozen=True, slots=True)
class CompactDecision:
    """Minimum immutable decision needed by the split-pipe routing hook."""

    tso_found: bool
    tso_start: int
    tso_end: int
    tso_distance: int
    tso_match_type: str
    filter_positive: bool
    sequence_ambiguous: bool
    confidence_tier: str


NO_TSO_DECISION = CompactDecision(
    False, -1, -1, -1, "none", False, False, "none"
)


class PackedGuideIndex:
    """Two-bit exact + complementary-seed guide index.

    The object is completely constructed before it is exposed.  Callers should
    preload one instance before preprocessing workers fork and treat it as
    immutable.  Public sequence groups and all index statistics are frozen
    dataclasses/tuples; lookup tables are private and are never modified during
    queries.
    """

    def __init__(
        self,
        guides: Iterable[str],
        *,
        guide_lengths: Sequence[int] = (20, 21),
        max_seed_bucket_size: int = 10_000,
    ) -> None:
        lengths = tuple(sorted({int(length) for length in guide_lengths}))
        if not lengths or any(length < 2 for length in lengths):
            raise ValueError("guide_lengths must contain lengths >= 2")
        if max_seed_bucket_size < 1:
            raise ValueError("max_seed_bucket_size must be >= 1")

        row_counts: Counter[str] = Counter()
        reference_rows = 0
        for row_index, raw_sequence in enumerate(guides, start=1):
            sequence = _validate_dna(raw_sequence, f"reference guide row {row_index}")
            if len(sequence) not in lengths:
                raise ValueError(
                    f"reference guide row {row_index} has length {len(sequence)}, "
                    f"not configured guide_lengths {lengths}"
                )
            row_counts[sequence] += 1
            reference_rows += 1
        if reference_rows == 0:
            raise ValueError("guide reference contains no rows")

        sequence_groups: List[SequenceGroup] = []
        groups_by_length: Dict[int, List[int]] = defaultdict(list)
        canonical_hash = hashlib.sha256()
        for sequence in sorted(row_counts, key=lambda item: (len(item), item)):
            length = len(sequence)
            group_digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()[:16]
            group = SequenceGroup(
                group_id=f"sg{length}_{group_digest}",
                length=length,
                sequence=sequence,
                packed_sequence=encode_dna(sequence),
                reference_row_count=row_counts[sequence],
            )
            group_index = len(sequence_groups)
            sequence_groups.append(group)
            groups_by_length[length].append(group_index)
            canonical_hash.update(
                f"{length}\t{sequence}\t{row_counts[sequence]}\n".encode("ascii")
            )

        self.sequence_groups: Tuple[SequenceGroup, ...] = tuple(sequence_groups)
        self.guide_lengths: Tuple[int, ...] = lengths
        self._exact_by_length: Dict[int, Dict[int, int]] = {}
        self._seed_by_length: Dict[int, Tuple[Dict[int, Tuple[int, ...]], ...]] = {}
        self._seed_masks_by_length: Dict[int, Tuple[int, int]] = {}
        self._exact_string_by_length: Dict[int, Dict[str, int]] = {}
        self._seed_string_by_length: Dict[
            int, Tuple[Dict[str, Tuple[int, ...]], Dict[str, Tuple[int, ...]]]
        ] = {}
        self._seed_split_by_length: Dict[int, int] = {}
        length_stats: List[LengthIndexStats] = []

        for length in lengths:
            group_indices = groups_by_length.get(length, [])
            exact: Dict[int, int] = {}
            midpoint = length // 2
            first_mask = _position_mask(length, range(0, midpoint))
            second_mask = _position_mask(length, range(midpoint, length))
            masks = (first_mask, second_mask)
            mutable_seed_tables: Tuple[Dict[int, List[int]], Dict[int, List[int]]] = (
                defaultdict(list),
                defaultdict(list),
            )
            for group_index in group_indices:
                group = self.sequence_groups[group_index]
                exact[group.packed_sequence] = group_index
                for seed_number, mask in enumerate(masks):
                    mutable_seed_tables[seed_number][group.packed_sequence & mask].append(
                        group_index
                    )

            frozen_seed_tables: List[Dict[int, Tuple[int, ...]]] = []
            histograms = []
            maximums = []
            for seed_number, mutable_table in enumerate(mutable_seed_tables):
                frozen_table = {
                    key: tuple(values) for key, values in mutable_table.items()
                }
                bucket_histogram = Counter(len(values) for values in frozen_table.values())
                maximum = max(bucket_histogram, default=0)
                if maximum > max_seed_bucket_size:
                    raise ValueError(
                        f"guide length {length} seed {seed_number} maximum bucket "
                        f"size {maximum} exceeds configured limit {max_seed_bucket_size}"
                    )
                frozen_seed_tables.append(frozen_table)
                histograms.append(tuple(sorted(bucket_histogram.items())))
                maximums.append(maximum)

            self._exact_by_length[length] = exact
            self._seed_by_length[length] = tuple(frozen_seed_tables)  # type: ignore[assignment]
            self._seed_masks_by_length[length] = masks
            self._exact_string_by_length[length] = {
                self.sequence_groups[index].sequence: index
                for index in group_indices
            }
            string_seed_tables = []
            for seed_number in range(2):
                mutable_string_table: Dict[str, List[int]] = defaultdict(list)
                for group_index in group_indices:
                    sequence = self.sequence_groups[group_index].sequence
                    key = sequence[:midpoint] if seed_number == 0 else sequence[midpoint:]
                    mutable_string_table[key].append(group_index)
                string_seed_tables.append(
                    {key: tuple(values) for key, values in mutable_string_table.items()}
                )
            self._seed_string_by_length[length] = (
                string_seed_tables[0],
                string_seed_tables[1],
            )
            self._seed_split_by_length[length] = midpoint
            length_stats.append(
                LengthIndexStats(
                    guide_length=length,
                    sequence_groups=len(group_indices),
                    reference_rows=sum(
                        self.sequence_groups[index].reference_row_count
                        for index in group_indices
                    ),
                    seed_masks=masks,
                    seed_position_counts=(midpoint, length - midpoint),
                    seed_key_counts=tuple(
                        len(table) for table in frozen_seed_tables
                    ),  # type: ignore[arg-type]
                    seed_max_bucket_sizes=tuple(maximums),  # type: ignore[arg-type]
                    seed_bucket_histograms=tuple(histograms),  # type: ignore[arg-type]
                )
            )

        self.stats = IndexStats(
            guide_lengths=lengths,
            reference_rows=reference_rows,
            unique_sequence_groups=len(sequence_groups),
            duplicate_reference_rows=reference_rows - len(sequence_groups),
            reference_sha256=canonical_hash.hexdigest(),
            by_length=tuple(length_stats),
        )

    def match_compact(
        self, query: str, max_mismatches: int = 1
    ) -> Optional[Tuple[int, Tuple[int, ...]]]:
        """Return ``(distance, group indexes)`` without negative-path objects.

        String exact and complementary-half lookups keep the overwhelmingly
        common random cDNA query in C-level dictionary operations.  Packed
        Hamming verification remains authoritative whenever a seed hits.
        """

        if max_mismatches not in (0, 1):
            raise ValueError("packed guide index supports max_mismatches 0 or 1")
        length = len(query)
        exact_table = self._exact_string_by_length.get(length)
        if exact_table is None:
            return None
        exact_index = exact_table.get(query)
        if exact_index is not None:
            return 0, (exact_index,)
        if max_mismatches == 0:
            return None

        midpoint = self._seed_split_by_length[length]
        first_table, second_table = self._seed_string_by_length[length]
        first = first_table.get(query[:midpoint], ())
        second = second_table.get(query[midpoint:], ())
        if not first and not second:
            return None
        if any(base not in DNA_BASES for base in query):
            return None

        if not first:
            candidate_indices = second
        elif not second:
            candidate_indices = first
        else:
            candidate_indices = tuple(set(first).union(second))
        packed_query = encode_dna(query)
        best_distance = max_mismatches + 1
        best_indices: List[int] = []
        for group_index in candidate_indices:
            group = self.sequence_groups[group_index]
            distance = packed_hamming_distance(
                packed_query, group.packed_sequence, length
            )
            if distance < best_distance:
                best_distance = distance
                best_indices = [group_index]
            elif distance == best_distance:
                best_indices.append(group_index)
        if best_distance > max_mismatches:
            return None
        return best_distance, tuple(
            sorted(best_indices, key=lambda index: self.sequence_groups[index].group_id)
        )

    @classmethod
    def from_reference_csv(
        cls,
        reference_path: str | Path,
        *,
        guide_lengths: Sequence[int] = (20, 21),
        max_seed_bucket_size: int = 10_000,
        sequence_column: str = "grna_seq",
    ) -> "PackedGuideIndex":
        """Load a comment-tolerant Step 0 reference CSV."""

        path = Path(reference_path)
        with path.open("r", newline="") as handle:
            reader = csv.DictReader(
                line for line in handle if not line.lstrip().startswith("#")
            )
            if not reader.fieldnames or sequence_column not in reader.fieldnames:
                raise ValueError(
                    f"guide reference {path} is missing required column "
                    f"{sequence_column!r}"
                )
            guides = [row[sequence_column] for row in reader]
        return cls(
            guides,
            guide_lengths=guide_lengths,
            max_seed_bucket_size=max_seed_bucket_size,
        )

    def match(self, query: str, max_mismatches: int = 1) -> GuideMatch:
        """Match one complete guide query using exact-first Hamming semantics."""

        if max_mismatches not in (0, 1):
            raise ValueError("packed guide index supports max_mismatches 0 or 1")
        query = str(query).strip().upper()
        length = len(query)
        if length not in self._exact_by_length:
            return GuideMatch(
                length, False, None, (), 0, False, False, 0, "unsupported_length"
            )
        if any(base not in DNA_BASES for base in query):
            return GuideMatch(
                length, False, None, (), 0, False, False, 0, "invalid_query_bases"
            )
        packed_query = encode_dna(query)

        # Exact matches never probe or depend on fuzzy seed buckets.
        exact_index = self._exact_by_length[length].get(packed_query)
        if exact_index is not None:
            group = self.sequence_groups[exact_index]
            return GuideMatch(
                query_length=length,
                supported=True,
                distance=0,
                sequence_group_ids=(group.group_id,),
                reference_row_count=group.reference_row_count,
                sequence_ambiguous=False,
                shared_reference_sequence=group.reference_row_count > 1,
                candidate_count=1,
                reason=(
                    "shared_sequence_exact"
                    if group.reference_row_count > 1
                    else "unique_sequence_exact"
                ),
            )
        if max_mismatches == 0:
            return GuideMatch(length, False, None, (), 0, False, False, 0, "no_match")

        masks = self._seed_masks_by_length[length]
        seed_tables = self._seed_by_length[length]
        candidate_indices = set()
        for seed_number, mask in enumerate(masks):
            candidate_indices.update(
                seed_tables[seed_number].get(packed_query & mask, ())
            )

        best_distance = max_mismatches + 1
        best_indices: List[int] = []
        for group_index in candidate_indices:
            group = self.sequence_groups[group_index]
            distance = packed_hamming_distance(
                packed_query, group.packed_sequence, length
            )
            if distance < best_distance:
                best_distance = distance
                best_indices = [group_index]
            elif distance == best_distance:
                best_indices.append(group_index)

        if best_distance > max_mismatches:
            return GuideMatch(
                length,
                False,
                None,
                (),
                0,
                False,
                False,
                len(candidate_indices),
                "no_match",
            )

        best_groups = sorted(
            (self.sequence_groups[index] for index in best_indices),
            key=lambda group: group.group_id,
        )
        ambiguous = len(best_groups) > 1
        shared = any(group.reference_row_count > 1 for group in best_groups)
        return GuideMatch(
            query_length=length,
            supported=True,
            distance=best_distance,
            sequence_group_ids=tuple(group.group_id for group in best_groups),
            reference_row_count=sum(
                group.reference_row_count for group in best_groups
            ),
            sequence_ambiguous=ambiguous,
            shared_reference_sequence=shared,
            candidate_count=len(candidate_indices),
            reason=("tied_sequences_hamming_1" if ambiguous else "unique_sequence_hamming_1"),
        )


def _position_mask(length: int, positions: Iterable[int]) -> int:
    """Return a two-bit packed mask for zero-based sequence positions."""

    mask = 0
    for position in positions:
        shift = 2 * (length - position - 1)
        mask |= 0b11 << shift
    return mask


class GexGrnaFilter:
    """TSO-aware wrapper around :class:`PackedGuideIndex`."""

    def __init__(self, index: PackedGuideIndex, settings: FilterSettings) -> None:
        if index.guide_lengths != settings.guide_lengths:
            raise ValueError(
                "index guide lengths do not match filter settings: "
                f"{index.guide_lengths} != {settings.guide_lengths}"
            )
        self.index = index
        self.settings = settings
        self.tso_matcher = TsoMatcher(
            settings.tso_sequence,
            settings.tso_max_edit_distance,
            settings.max_tso_start,
            settings.tso_backend,
        )

    @classmethod
    def from_config(cls, config_path: str | Path) -> "GexGrnaFilter":
        """Load, validate, and construct the Step 1a filter from main YAML."""

        config_file = Path(config_path).resolve()
        with config_file.open("r") as handle:
            config = yaml.safe_load(handle) or {}
        settings = settings_from_config(config, config_file)
        index = PackedGuideIndex.from_reference_csv(
            settings.reference_path,
            guide_lengths=settings.guide_lengths,
            max_seed_bucket_size=settings.max_seed_bucket_size,
        )
        return cls(index, settings)

    def classify(
        self,
        read_id: str,
        r1_sequence: str,
        raw_pair_location: Optional[str] = None,
    ) -> FilterDecision:
        """Classify one raw R1 without making a biological identity call."""

        read = str(r1_sequence).strip().upper()
        tso = self.tso_matcher.find_normalized(read)
        location = raw_pair_location or ""
        if not tso.found:
            return _negative_decision(read_id, location, tso, "no_tso")

        post_tso = read[tso.end :]
        available_lengths = tuple(
            length
            for length in self.settings.guide_lengths
            if len(post_tso) >= length
        )
        if not available_lengths:
            return _negative_decision(
                read_id, location, tso, "insufficient_full_guide"
            )

        length_matches = tuple(
            self.index.match(
                post_tso[:length], self.settings.max_guide_mismatches
            )
            for length in available_lengths
        )
        supported_matches = tuple(match for match in length_matches if match.supported)
        if not supported_matches:
            reason = (
                "invalid_guide_bases"
                if all(match.reason == "invalid_query_bases" for match in length_matches)
                else "no_reference_match"
            )
            return _negative_decision(
                read_id, location, tso, reason, length_matches=length_matches
            )

        # Evaluate every full configured length, then retain the global minimum
        # distance.  This makes mixed 20/21-mer behavior independent of query
        # order and preserves cross-length ties as ambiguity diagnostics.
        best_distance = min(
            int(match.distance) for match in supported_matches if match.distance is not None
        )
        best_matches = tuple(
            match for match in supported_matches if match.distance == best_distance
        )
        group_ids = tuple(
            sorted(
                {
                    group_id
                    for match in best_matches
                    for group_id in match.sequence_group_ids
                }
            )
        )
        ambiguous = len(group_ids) > 1
        shared = any(match.shared_reference_sequence for match in best_matches)
        reference_supported = True
        filter_positive = not ambiguous or self.settings.remove_ambiguous_hits
        confidence = "exact" if best_distance == 0 else "hamming_1"
        if ambiguous:
            reason = (
                f"reference_supported_ambiguous_{confidence}"
                if filter_positive
                else f"reference_supported_ambiguous_{confidence}_not_removed"
            )
        elif shared:
            reason = f"reference_supported_shared_sequence_{confidence}"
        else:
            reason = f"reference_supported_unique_{confidence}"

        return FilterDecision(
            read_id=str(read_id),
            raw_pair_location=location,
            tso_found=True,
            tso_start=tso.start,
            tso_end=tso.end,
            tso_distance=tso.distance,
            tso_match_type=tso.match_type,
            reference_supported=reference_supported,
            filter_positive=filter_positive,
            filter_reason=reason,
            confidence_tier=confidence,
            guide_lengths=tuple(sorted({match.query_length for match in best_matches})),
            guide_distance=best_distance,
            sequence_group_ids=group_ids,
            sequence_ambiguous=ambiguous,
            shared_reference_sequence=shared,
            candidate_count=sum(match.candidate_count for match in best_matches),
            length_matches=length_matches,
        )

    def route(self, r1_sequence: str) -> CompactDecision:
        """Return only the fields consumed by the default preprocessing path.

        Rich per-read diagnostics remain available through :meth:`classify`.
        The production hook uses this allocation-minimal route unless
        ``write_read_decisions`` is explicitly enabled.
        """

        read = str(r1_sequence).strip().upper()
        tso = self.tso_matcher.find_normalized(read)
        return self._route_normalized_after_tso(read, tso)

    def _route_normalized_after_tso(
        self, read: str, tso: TsoMatch
    ) -> CompactDecision:
        """Finish compact routing after a timed/precomputed TSO decision.

        This internal split keeps the production hot path unchanged while
        allowing the remediation profiler to measure TSO and guide routing
        independently on identical records.
        """

        if not tso.found:
            return NO_TSO_DECISION
        post_tso = read[tso.end :]
        if not any(len(post_tso) >= length for length in self.settings.guide_lengths):
            return CompactDecision(
                True,
                tso.start,
                tso.end,
                tso.distance,
                tso.match_type,
                False,
                False,
                "none",
            )
        best_distance: Optional[int] = None
        best_group_indexes = set()
        for length in self.settings.guide_lengths:
            if len(post_tso) < length:
                continue
            distance_limit = self.settings.max_guide_mismatches
            if best_distance is not None:
                distance_limit = min(distance_limit, best_distance)
            match = self.index.match_compact(post_tso[:length], distance_limit)
            if match is None:
                continue
            distance, group_indexes = match
            if best_distance is None or distance < best_distance:
                best_distance = distance
                best_group_indexes = set(group_indexes)
            elif distance == best_distance:
                best_group_indexes.update(group_indexes)
        if best_distance is None:
            return CompactDecision(
                True,
                tso.start,
                tso.end,
                tso.distance,
                tso.match_type,
                False,
                False,
                "none",
            )
        ambiguous = len(best_group_indexes) > 1
        return CompactDecision(
            True,
            tso.start,
            tso.end,
            tso.distance,
            tso.match_type,
            not ambiguous or self.settings.remove_ambiguous_hits,
            ambiguous,
            "exact" if best_distance == 0 else "hamming_1",
        )

    def classify_many(
        self, records: Iterable[Tuple[str, str]]
    ) -> Iterator[FilterDecision]:
        """Stream decisions for ``(read_id, raw_r1_sequence)`` records."""

        for read_id, sequence in records:
            yield self.classify(read_id, sequence)


def _negative_decision(
    read_id: str,
    raw_pair_location: str,
    tso: TsoMatch,
    reason: str,
    *,
    length_matches: Tuple[GuideMatch, ...] = (),
) -> FilterDecision:
    return FilterDecision(
        read_id=str(read_id),
        raw_pair_location=raw_pair_location,
        tso_found=tso.found,
        tso_start=tso.start,
        tso_end=tso.end,
        tso_distance=tso.distance,
        tso_match_type=tso.match_type,
        reference_supported=False,
        filter_positive=False,
        filter_reason=reason,
        confidence_tier="none",
        guide_lengths=(),
        guide_distance=None,
        sequence_group_ids=(),
        sequence_ambiguous=False,
        shared_reference_sequence=False,
        candidate_count=sum(match.candidate_count for match in length_matches),
        length_matches=length_matches,
    )


def settings_from_config(config: Mapping, config_path: str | Path) -> FilterSettings:
    """Validate the v1.4 Step 1a config and resolve its reference path."""

    try:
        structure = config["grna_structure"]
        paths = config["paths"]
    except (KeyError, TypeError) as exc:
        raise ValueError("config requires grna_structure and paths mappings") from exc
    decontamination = config.get("gex_decontamination", {}) or {}
    filter_matching = decontamination.get("filter_matching", {}) or {}
    resources = decontamination.get("resource_limits", {}) or {}

    guide_distance = filter_matching.get("guide_distance", "hamming")
    if guide_distance != "hamming":
        raise ValueError("packed Step 1a requires guide_distance: hamming")
    if bool(filter_matching.get("allow_guide_indels", False)):
        raise ValueError("packed Step 1a does not allow guide indels")
    if not bool(filter_matching.get("full_guide_required", True)):
        raise ValueError("packed Step 1a requires full_guide_required: true")
    backend = filter_matching.get("backend", "packed_inline")
    if backend != "packed_inline":
        raise ValueError(
            f"GexGrnaFilter supports backend packed_inline, not {backend!r}"
        )
    index_name = filter_matching.get("index", "two_bit_two_seed")
    if index_name != "two_bit_two_seed":
        raise ValueError(
            f"GexGrnaFilter supports index two_bit_two_seed, not {index_name!r}"
        )

    guide_lengths = tuple(sorted({int(value) for value in structure["grna_lengths"]}))
    max_mismatches = int(filter_matching.get("max_guide_mismatches", 1))
    if max_mismatches not in (0, 1):
        raise ValueError("packed Step 1a supports max_guide_mismatches 0 or 1")
    tso_sequence = _validate_dna(structure["fixed_seq_1"], "grna_structure.fixed_seq_1")

    raw_reference_path = Path(str(paths["grna_reference"])).expanduser()
    if not raw_reference_path.is_absolute():
        raw_reference_path = Path(config_path).resolve().parent / raw_reference_path
    reference_path = str(raw_reference_path.resolve())

    return FilterSettings(
        tso_sequence=tso_sequence,
        tso_max_edit_distance=int(
            decontamination.get(
                "tso_max_edit_dist", structure.get("tso_max_edit_dist", 3)
            )
        ),
        max_tso_start=int(decontamination.get("max_tso_start", 2)),
        guide_lengths=guide_lengths,
        max_guide_mismatches=max_mismatches,
        remove_ambiguous_hits=bool(
            filter_matching.get("remove_reference_supported_ambiguous_hits", True)
        ),
        max_seed_bucket_size=int(resources.get("max_seed_bucket_size", 10_000)),
        reference_path=reference_path,
        tso_backend=str(
            (decontamination.get("tso_matching", {}) or {}).get(
                "backend", "edlib"
            )
        ).strip().lower(),
    )


def load_filter_from_config(config_path: str | Path) -> GexGrnaFilter:
    """Convenience integration API used by the preprocessing hook."""

    return GexGrnaFilter.from_config(config_path)


def _main() -> None:
    parser = argparse.ArgumentParser(
        description="Classify one raw GEX R1 with the packed Step 1a filter"
    )
    parser.add_argument("--config", required=True, help="Main pipeline YAML")
    parser.add_argument("--sequence", help="Raw R1 sequence to classify")
    parser.add_argument("--read-id", default="query")
    parser.add_argument(
        "--stats", action="store_true", help="Print immutable index statistics"
    )
    args = parser.parse_args()

    classifier = load_filter_from_config(args.config)
    payload = {}
    if args.stats:
        payload["index_stats"] = asdict(classifier.index.stats)
    if args.sequence is not None:
        payload["decision"] = classifier.classify(
            args.read_id, args.sequence
        ).to_dict()
    if not payload:
        parser.error("provide --sequence and/or --stats")
    print(json.dumps(payload, sort_keys=True))


if __name__ == "__main__":
    _main()

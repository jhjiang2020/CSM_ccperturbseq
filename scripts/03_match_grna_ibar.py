#!/usr/bin/env python3
"""
Step 3: Match gRNA and iBAR sequences from R1 reads.

Pipeline:
  1. Locate Split_TSO in R1 via semi-global alignment (handles variable offset)
  2. Extract gRNA + iBAR at computed positions relative to TSO end
  3. Match gRNA/iBAR against reference using fast fuzzy matching

Alignment backend:
  - Uses `edlib` (bit-parallel Myers algorithm) if available — ~50x faster
  - Falls back to pure-Python sliding-window Hamming if edlib is not installed
  - Install via: pip install edlib

R1 Read Structure (CROPseq-multi / Parse Biosciences):
  ┌───────────────────────────────────────────────────────────────────┐
  │ N7 random prefix  │ Split_TSO (30bp) │ gRNA │ stem │ iBAR │ tracr│
  └───────────────────────────────────────────────────────────────────┘
  In intact dedicated-gRNA reads the TSO starts at offset 7 after the N7 random
  prefix. i7/i5 are separate index reads; no additional N3 R1 prefix is assumed.
  GEX-derived guide reads start at offset 0. We use semi-global alignment to
  tolerate damage/error and find the TSO endpoint before extracting downstream
  sequence.

Usage:
  python scripts/03_match_grna_ibar.py --config config.yaml
  python scripts/03_match_grna_ibar.py --config config.yaml --library C1 --n_processes 4
"""

import argparse
import os
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

try:
    from scripts.grna_preflight import check_ibar_reference, preflight_options, print_preflight
    from scripts.guide_design import DESIGN_COLUMNS, DesignMatcher, compile_design, emit_warnings, enabled, extract_unit
    from scripts.design_preflight import check_design
except ModuleNotFoundError:
    from grna_preflight import check_ibar_reference, preflight_options, print_preflight
    from guide_design import DESIGN_COLUMNS, DesignMatcher, compile_design, emit_warnings, enabled, extract_unit
    from design_preflight import check_design


def resolve_output_root(output_dir, output_namespace=None):
    """Resolve an optional source namespace, preserving historical defaults."""
    if not output_namespace:
        return Path(output_dir)
    namespace = Path(output_namespace).expanduser()
    return namespace if namespace.is_absolute() else Path(output_dir) / namespace

# Try to import edlib for fast alignment
try:
    import edlib
    HAS_EDLIB = True
except ImportError:
    HAS_EDLIB = False


# ---------------------------------------------------------------------------
# Alignment: TSO finding
# ---------------------------------------------------------------------------

def find_tso_edlib(r1_seq, tso_seq, max_dist=3, search_end=50):
    """Find TSO in R1 using edlib semi-global alignment (infix/HW mode).

    edlib's HW mode finds the best alignment of the query (TSO) within
    the target (R1), allowing the query to start and end at any position.
    This naturally handles the variable offset.

    Args:
        r1_seq: R1 read sequence
        tso_seq: Split_TSO sequence (query)
        max_dist: Max edit distance for TSO match
        search_end: Only search within first N bases of R1

    Returns:
        (tso_end_position, edit_distance) or (-1, -1) if not found
    """
    target = r1_seq[:search_end].upper()
    tso_seq = tso_seq.upper()
    result = edlib.align(tso_seq, target, mode="HW", k=max_dist, task="locations")

    if result['editDistance'] < 0:
        return -1, -1

    # Get the best (rightmost end) location — this is where TSO ends
    locations = result['locations']
    if not locations:
        return -1, -1

    # Pick the location whose end matches expectations best
    # (closest to expected TSO length from start)
    best_loc = locations[0]
    tso_end = best_loc[1] + 1  # edlib locations are inclusive, +1 for exclusive end
    return tso_end, result['editDistance']


def find_tso_sliding(r1_seq, tso_seq, max_dist=3, search_end=50):
    """Find TSO in R1 using sliding-window Hamming distance (pure Python fallback).

    Slides the TSO across the first `search_end` bases of R1, computing
    Hamming distance at each position. Returns the position with minimum
    mismatches.

    Args:
        r1_seq: R1 read sequence
        tso_seq: Split_TSO sequence
        max_dist: Max allowed Hamming distance
        search_end: Search limit from start of R1

    Returns:
        (tso_end_position, hamming_distance) or (-1, -1) if not found
    """
    tso_len = len(tso_seq)
    max_start = min(search_end - tso_len, 20)  # TSO unlikely to start >20bp in

    if max_start < 0:
        return -1, -1

    best_mm = max_dist + 1
    best_end = -1

    for start in range(0, max_start + 1):
        observed = r1_seq[start:start + tso_len]
        if len(observed) < tso_len:
            break

        mm = 0
        for c1, c2 in zip(observed, tso_seq):
            if c1 != c2:
                mm += 1
                if mm > best_mm:
                    break  # Early termination within Hamming
        else:
            if mm < best_mm:
                best_mm = mm
                best_end = start + tso_len

    if best_mm <= max_dist:
        return best_end, best_mm
    return -1, -1


def find_tso(r1_seq, tso_seq, max_dist=3, search_end=50):
    """Find TSO in R1 — dispatches to edlib (fast) or sliding (fallback)."""
    if HAS_EDLIB:
        return find_tso_edlib(r1_seq, tso_seq, max_dist, search_end)
    else:
        return find_tso_sliding(r1_seq, tso_seq, max_dist, search_end)


# ---------------------------------------------------------------------------
# Levenshtein distance with early termination
# ---------------------------------------------------------------------------

def levenshtein_distance(s1, s2, max_dist=None):
    """Compute Levenshtein edit distance with optional early termination."""
    if HAS_EDLIB:
        k = max_dist if max_dist is not None else -1
        result = edlib.align(s1, s2, mode="NW", k=k)
        d = result['editDistance']
        if d < 0:
            return (max_dist + 1) if max_dist is not None else len(s1) + len(s2)
        return d

    # Pure Python fallback
    len1, len2 = len(s1), len(s2)
    if max_dist is not None and abs(len1 - len2) > max_dist:
        return max_dist + 1
    if len1 < len2:
        return levenshtein_distance(s2, s1, max_dist)
    if len2 == 0:
        return len1

    prev_row = list(range(len2 + 1))
    for i, c1 in enumerate(s1):
        curr_row = [i + 1]
        row_min = i + 1
        for j, c2 in enumerate(s2):
            cost = 0 if c1 == c2 else 1
            val = min(prev_row[j + 1] + 1, curr_row[j] + 1, prev_row[j] + cost)
            curr_row.append(val)
            if val < row_min:
                row_min = val
        if max_dist is not None and row_min > max_dist:
            return max_dist + 1
        prev_row = curr_row
    return prev_row[-1]


def hamming_distance(s1, s2, max_mm=None):
    """Fast Hamming distance for equal-length strings."""
    d = 0
    for c1, c2 in zip(s1, s2):
        if c1 != c2:
            d += 1
            if max_mm is not None and d > max_mm:
                return d
    return d


# ---------------------------------------------------------------------------
# NeighborhoodMatcher: fast fuzzy lookup via precomputed neighborhoods
# ---------------------------------------------------------------------------

class NeighborhoodMatcher:
    """Fast approximate string matcher.

    Precomputes deletion and substitution neighborhoods for O(1) lookup
    at Levenshtein distance ≤ 1, with progressively wider search for
    larger distances using edlib when available.

    For dist 0: exact dict lookup (1 lookup)
    For dist 1: Hamming-1 substitutions (3L lookups) + 1-deletion (2L lookups)
    For dist 2+: edlib global alignment against candidates filtered by length
    """

    def __init__(self, seqs, names, max_dist=3):
        self.seqs = list(seqs)
        self.names = list(names)
        self.max_dist = max_dist
        self.n_refs = len(self.seqs)

        # Exact match dict: seq -> list of indices (multi-index for shared seqs)
        self.exact = {}
        for i, seq in enumerate(self.seqs):
            self.exact.setdefault(seq, []).append(i)

        # 1-deletion neighborhood: deletion_variant -> [ref_indices]
        self.del1 = defaultdict(list)
        for i, seq in enumerate(self.seqs):
            for j in range(len(seq)):
                variant = seq[:j] + seq[j+1:]
                self.del1[variant].append(i)

        # Length-grouped refs for brute-force fallback
        self.by_length = defaultdict(list)
        for i, seq in enumerate(self.seqs):
            self.by_length[len(seq)].append(i)

    def match(self, query):
        """Find best match(es) for query.

        Returns: (name_or_names, distance, is_unique)
          - Unique match: (str, int, True)
          - Position-ambiguous (same seq, different names): (list[str], int, False)
          - Sequence-ambiguous (different seqs at same dist): (None, int, False)
          - No match: (None, None, False)
        """
        # --- Distance 0: exact ---
        indices = self.exact.get(query)
        if indices is not None:
            if len(indices) == 1:
                return self.names[indices[0]], 0, True
            else:
                # Same sequence, multiple names (shared gRNA) — caller must disambiguate
                return [self.names[i] for i in indices], 0, False

        if self.max_dist == 0:
            return None, None, False

        # --- Distance 1: substitution + deletion neighborhood ---
        candidates_d1 = set()

        # Hamming-1 substitutions of query → exact dict
        for j in range(len(query)):
            for base in 'ACGT':
                if base != query[j]:
                    variant = query[:j] + base + query[j+1:]
                    for idx in self.exact.get(variant, []):
                        candidates_d1.add(idx)

        # 1-deletion of query → exact dict (catches 1-insertion in ref)
        for j in range(len(query)):
            qdel = query[:j] + query[j+1:]
            for idx in self.exact.get(qdel, []):
                candidates_d1.add(idx)

        # Query → ref 1-deletion dict (catches 1-deletion in ref)
        for idx in self.del1.get(query, []):
            candidates_d1.add(idx)

        # 1-insertion to query at each position → exact dict
        for j in range(len(query) + 1):
            for base in 'ACGT':
                variant = query[:j] + base + query[j:]
                for idx in self.exact.get(variant, []):
                    candidates_d1.add(idx)

        if candidates_d1:
            # Verify with exact Levenshtein
            best_dist = self.max_dist + 1
            best_indices = []
            for idx in candidates_d1:
                d = levenshtein_distance(query, self.seqs[idx], max_dist=1)
                if d < best_dist:
                    best_dist = d
                    best_indices = [idx]
                elif d == best_dist:
                    best_indices.append(idx)

            if best_dist <= 1 and len(best_indices) >= 1:
                # Deduplicate by sequence (multiple indices may share the same seq)
                unique_seqs = {}
                for idx in best_indices:
                    unique_seqs.setdefault(self.seqs[idx], []).append(idx)
                if len(unique_seqs) == 1:
                    # All indices share the same sequence
                    all_idx = best_indices
                    if len(all_idx) == 1:
                        return self.names[all_idx[0]], best_dist, True
                    else:
                        # Same seq, multiple names — position-ambiguous
                        return [self.names[i] for i in all_idx], best_dist, False
                else:
                    # Different sequences at same distance — sequence-ambiguous
                    return None, best_dist, False

        if self.max_dist < 2:
            return None, None, False

        # --- Distance 2+: length-filtered brute-force with edlib ---
        candidates = set()
        for offset in range(-self.max_dist, self.max_dist + 1):
            candidates.update(self.by_length.get(len(query) + offset, []))
        candidates -= candidates_d1  # Don't re-check

        best_dist = self.max_dist + 1
        best_indices = []

        for idx in candidates:
            d = levenshtein_distance(query, self.seqs[idx], max_dist=self.max_dist)
            if d < best_dist:
                best_dist = d
                best_indices = [idx]
            elif d == best_dist:
                # Distance <=1 was exhausted. Two distinct distance-2 hits
                # prove an irreducible tie; remaining references cannot fix it.
                if d == 2 and self.seqs[idx] != self.seqs[best_indices[0]]:
                    return None, 2, False
                best_indices.append(idx)

        if best_dist <= self.max_dist and len(best_indices) >= 1:
            unique_seqs = {}
            for idx in best_indices:
                unique_seqs.setdefault(self.seqs[idx], []).append(idx)
            if len(unique_seqs) == 1:
                all_idx = best_indices
                if len(all_idx) == 1:
                    return self.names[all_idx[0]], best_dist, True
                else:
                    return [self.names[i] for i in all_idx], best_dist, False
            else:
                return None, best_dist, False

        return None, None, False


# ---------------------------------------------------------------------------
# gRNA/iBAR extraction from R1
# ---------------------------------------------------------------------------

def extract_grna_ibar(
    r1_seq, tso_seq, stem1, stem2, grna_lengths, ibar_length,
    tracr1=None, tracr2=None, anchor_max_mm=2, tso_max_dist=3
):
    """Extract gRNA and iBAR from R1 read.

    Step 1: Find TSO anywhere in R1 via semi-global alignment
    Step 2: Try each gRNA length, score combined stem+tracr anchors
    Step 3: Pick the variant with lowest combined score

    Stem+tracr combined scoring:
      - Variant A: stem1 + tracr1 (correct pair for position 1)
      - Variant B: stem2 + tracr2 (correct pair for position 2)
      - anchor_max_mm applies to the COMBINED stem_mm + tracr_mm total

    Returns:
        dict with grna_seq, ibar_seq, variant, anchor_combo, etc.
        or None if extraction fails.
    """
    # Normalize user-provided reference sequences to uppercase
    r1_seq = r1_seq.upper() if r1_seq else ""
    tso_seq = tso_seq.upper() if tso_seq else ""
    stem1 = stem1.upper() if stem1 else ""
    stem2 = stem2.upper() if stem2 else ""
    tracr1 = tracr1.upper() if tracr1 else ""
    tracr2 = tracr2.upper() if tracr2 else ""

    # Step 1: Find TSO
    tso_end, tso_dist = find_tso(r1_seq, tso_seq, max_dist=tso_max_dist)
    if tso_end < 0:
        return None

    stem_len = len(stem1)
    tracr_len = len(tracr1) if tracr1 else 0

    # Step 2: Try each gRNA length
    for grna_len in sorted(grna_lengths, reverse=True):
        grna_start = tso_end
        grna_end = grna_start + grna_len
        stem_start = grna_end
        stem_end = stem_start + stem_len
        ibar_start = stem_end
        ibar_end = ibar_start + ibar_length

        if ibar_end > len(r1_seq):
            continue

        stem_observed = r1_seq[stem_start:stem_end]
        mm_stem1 = hamming_distance(stem_observed, stem1, anchor_max_mm + 1)
        if stem2:
            mm_stem2 = hamming_distance(stem_observed, stem2, anchor_max_mm + 1)
        else:
            mm_stem2 = anchor_max_mm + 2

        # Compute tracr mismatches if available
        mm_tracr1, mm_tracr2 = 0, 0
        has_tracr = False
        if tracr1 and ibar_end + tracr_len <= len(r1_seq):
            tracr_observed = r1_seq[ibar_end:ibar_end + tracr_len]
            mm_tracr1 = hamming_distance(tracr_observed, tracr1, anchor_max_mm + 1)
            if tracr2:
                mm_tracr2 = hamming_distance(tracr_observed, tracr2, anchor_max_mm + 1)
            else:
                mm_tracr2 = anchor_max_mm + 2
            has_tracr = True

        # Combined scores for all 4 combinations
        score_A = mm_stem1 + mm_tracr1   # stem1 + tracr1 (correct for pos1)
        score_B = mm_stem2 + mm_tracr2   # stem2 + tracr2 (correct for pos2)
        score_cross_12 = mm_stem1 + mm_tracr2  # stem1 + tracr2 (chimeric)
        score_cross_21 = mm_stem2 + mm_tracr1  # stem2 + tracr1 (chimeric)

        # Pick best pair (A or B only — cross-pairs are never "correct")
        if score_A < score_B and score_A <= anchor_max_mm:
            variant = 1
            mm_stem = mm_stem1
            mm_tracr = mm_tracr1
        elif score_B < score_A and score_B <= anchor_max_mm:
            variant = 2
            mm_stem = mm_stem2
            mm_tracr = mm_tracr2
        elif score_A == score_B and score_A <= anchor_max_mm:
            variant = 'AMBIGUOUS'
            mm_stem = mm_stem1
            mm_tracr = mm_tracr1
        else:
            continue

        grna_seq = r1_seq[grna_start:grna_end]
        ibar_seq = r1_seq[ibar_start:ibar_end]

        # Classify anchor combination for stats
        if variant == 1:
            anchor_combo = 'A'  # stem1 + tracr1
        elif variant == 2:
            anchor_combo = 'B'  # stem2 + tracr2
        else:
            anchor_combo = 'ambiguous'

        # Flag cross-pair if a chimeric combo scores better than the chosen one
        chosen_score = score_A if variant == 1 else score_B
        if has_tracr:
            if score_cross_12 < chosen_score or score_cross_21 < chosen_score:
                anchor_combo = 'cross'

        return {
            'grna_seq': grna_seq,
            'ibar_seq': ibar_seq,
            'variant': variant,
            'grna_length': grna_len,
            'tso_offset': tso_end - len(tso_seq),
            'tso_dist': tso_dist,
            'stem_mm': mm_stem,
            'tracr_mm': mm_tracr,
            'anchor_combo': anchor_combo,
        }

    return None


def extract_grna_single_short(
    r1_seq, tso_seq, stem1, grna_lengths, anchor_max_mm=2, tso_max_dist=3,
):
    """Extract gRNA from R1 in single-guide no-iBAR mode with short-read support.

    Handles four scenarios based on available read length after TSO:
      Scenario 1: Full gRNA + full stem  → normal extraction (best confidence)
      Scenario 2: Full gRNA + partial stem → partial-stem boundary check
      Scenario 3: Truncated gRNA (no stem) → rescue whatever follows TSO
      Scenario 4: No TSO or nothing after TSO → discard

    Partial stem mismatch rules:
      - partial_stem_len <= 2: 0 mismatches (strict, overrides user config)
      - partial_stem_len > 2:  min(1, anchor_max_mm) mismatches

    Returns:
        dict with grna_seq, scenario, stem_truncated, grna_truncated, etc.
        or None if extraction fails (Scenario 4).
    """
    r1_seq = r1_seq.upper() if r1_seq else ""
    tso_seq = tso_seq.upper() if tso_seq else ""
    stem1 = stem1.upper() if stem1 else ""

    # Step 1: Find TSO
    tso_end, tso_dist = find_tso(r1_seq, tso_seq, max_dist=tso_max_dist)
    if tso_end < 0:
        return None  # Scenario 4: no TSO

    remaining = len(r1_seq) - tso_end
    if remaining <= 0:
        return None  # Scenario 4: nothing after TSO

    stem_len = len(stem1)
    min_grna_len = min(grna_lengths)
    max_grna_len = max(grna_lengths)

    # --- Scenario 1 & 2: Try each gRNA length (longest first) ---
    # We attempt full-stem and partial-stem evaluation for each candidate length.
    best_result = None
    best_score = None  # (scenario, -partial_stem_len_or_full, stem_mm)

    for grna_len in sorted(grna_lengths, reverse=True):
        grna_end = tso_end + grna_len
        if grna_end > len(r1_seq):
            continue  # Not enough room for even this gRNA length

        avail_after_grna = len(r1_seq) - grna_end

        if avail_after_grna >= stem_len:
            # --- Scenario 1: Full stem available ---
            stem_observed = r1_seq[grna_end:grna_end + stem_len]
            mm = hamming_distance(stem_observed, stem1, anchor_max_mm + 1)
            if mm <= anchor_max_mm:
                result = {
                    'grna_seq': r1_seq[tso_end:grna_end],
                    'ibar_seq': '',
                    'variant': 1,
                    'grna_length': grna_len,
                    'tso_offset': tso_end - len(tso_seq),
                    'tso_dist': tso_dist,
                    'stem_mm': mm,
                    'tracr_mm': 0,
                    'anchor_combo': 'A',
                    'scenario': 1,
                    'stem_truncated': False,
                    'grna_truncated': False,
                    'partial_stem_len': stem_len,
                }
                # Scenario 1 with full stem is highest confidence; return immediately
                return result

        elif avail_after_grna > 0:
            # --- Scenario 2: Partial stem available ---
            partial_stem_len = avail_after_grna
            stem_observed = r1_seq[grna_end:grna_end + partial_stem_len]
            stem_expected = stem1[:partial_stem_len]

            # Strict mismatch rules for partial stems
            if partial_stem_len <= 2:
                allowed_mm = 0
            else:
                allowed_mm = min(1, anchor_max_mm)

            mm = hamming_distance(stem_observed, stem_expected, allowed_mm + 1)
            if mm <= allowed_mm:
                # Valid partial stem match — track as candidate
                # Prefer: lower scenario number, then longer partial stem, then fewer mm
                score = (2, -partial_stem_len, mm)
                if best_score is None or score < best_score:
                    best_score = score
                    best_result = {
                        'grna_seq': r1_seq[tso_end:grna_end],
                        'ibar_seq': '',
                        'variant': 1,
                        'grna_length': grna_len,
                        'tso_offset': tso_end - len(tso_seq),
                        'tso_dist': tso_dist,
                        'stem_mm': mm,
                        'tracr_mm': 0,
                        'anchor_combo': 'partial_stem',
                        'scenario': 2,
                        'stem_truncated': True,
                        'grna_truncated': False,
                        'partial_stem_len': partial_stem_len,
                    }
        else:
            # avail_after_grna == 0: gRNA fills exactly to end of read
            # This is a zero-partial-stem case; trivially matches 0mm
            # but we need to resolve ambiguity: prefer shortest gRNA length
            score = (2, 0, 0)  # partial_stem_len=0 is worst within scenario 2
            if best_score is None or score < best_score:
                best_score = score
                best_result = {
                    'grna_seq': r1_seq[tso_end:grna_end],
                    'ibar_seq': '',
                    'variant': 1,
                    'grna_length': grna_len,
                    'tso_offset': tso_end - len(tso_seq),
                    'tso_dist': tso_dist,
                    'stem_mm': 0,
                    'tracr_mm': 0,
                    'anchor_combo': 'no_anchor',
                    'scenario': 2,
                    'stem_truncated': True,
                    'grna_truncated': False,
                    'partial_stem_len': 0,
                }

    if best_result is not None:
        return best_result

    # --- Scenario 3: Truncated gRNA (remaining < min gRNA length) ---
    if remaining > 0:
        grna_seq = r1_seq[tso_end:]
        return {
            'grna_seq': grna_seq,
            'ibar_seq': '',
            'variant': 1,
            'grna_length': len(grna_seq),
            'tso_offset': tso_end - len(tso_seq),
            'tso_dist': tso_dist,
            'stem_mm': -1,
            'tracr_mm': -1,
            'anchor_combo': 'truncated',
            'scenario': 3,
            'stem_truncated': True,
            'grna_truncated': True,
            'partial_stem_len': 0,
        }

    return None  # Should not reach here, but defensive


# ---------------------------------------------------------------------------
# Short-read sampling: detect truncated gRNA lengths from first N reads
# ---------------------------------------------------------------------------

def sample_truncated_lengths(barcode_tsv_path, tso_seq, stem1, grna_lengths,
                             anchor_max_mm=2, tso_max_dist=3,
                             n_sample=100_000):
    """Sample first N reads to detect all naturally-occurring gRNA lengths.

    Used in single-guide no-iBAR mode to pre-build truncated matchers.
    Returns a sorted list of observed gRNA lengths (including truncated).
    """
    observed_lengths = Counter()
    scenario_counts = Counter()
    n_read = 0

    for chunk_df in pd.read_csv(barcode_tsv_path, sep='\t', chunksize=10_000):
        for r1_seq in chunk_df['r1_seq'].values:
            result = extract_grna_single_short(
                r1_seq, tso_seq, stem1, grna_lengths,
                anchor_max_mm, tso_max_dist,
            )
            if result is not None:
                observed_lengths[result['grna_length']] += 1
                scenario_counts[result['scenario']] += 1
            else:
                scenario_counts['discard'] += 1

            n_read += 1
            if n_read >= n_sample:
                break
        if n_read >= n_sample:
            break

    return sorted(observed_lengths.keys()), dict(observed_lengths), dict(scenario_counts)


# ---------------------------------------------------------------------------
# Chunk processing for multiprocessing
# ---------------------------------------------------------------------------

# Module-level worker state — set once per worker via its initializer
_worker_state = {}


def _match_ibar_within_guide(query_ibar, guide_ibars, max_dist):
    """Match iBAR against only the iBARs paired with a specific gRNA.

    Brute-force Levenshtein against a small set (~10-20 iBARs per gRNA).
    Returns (matched_ibar_seq, distance) or (None, None).
    """
    best_seq = None
    best_dist = max_dist + 1
    n_best = 0

    for ref_ibar in guide_ibars:
        d = levenshtein_distance(query_ibar, ref_ibar, max_dist=max_dist)
        if d < best_dist:
            best_dist = d
            best_seq = ref_ibar
            n_best = 1
        elif d == best_dist:
            n_best += 1

    if best_seq is not None and best_dist <= max_dist and n_best == 1:
        return best_seq, best_dist
    return None, None


def _init_worker(config_dict):
    """Initialize per-worker state (called once per process, not per chunk).

    This avoids pickle-serializing the matchers and lookup dicts with every
    chunk, and avoids rebuilding NeighborhoodMatcher ~300 times.
    """
    _worker_state['design_contract'] = config_dict.get('design_contract')
    if config_dict.get('design_contract'):
        _worker_state['design_matcher'] = DesignMatcher(config_dict['design_contract'])
    _worker_state['tso_seq'] = config_dict['tso_seq']
    _worker_state['stem1'] = config_dict['stem1']
    _worker_state['stem2'] = config_dict['stem2']
    _worker_state['tracr1'] = config_dict['tracr1']
    _worker_state['tracr2'] = config_dict['tracr2']
    _worker_state['grna_lengths'] = config_dict['grna_lengths']
    _worker_state['ibar_length'] = config_dict['ibar_length']
    _worker_state['anchor_max_mm'] = config_dict['anchor_max_mm']
    _worker_state['tso_max_dist'] = config_dict['tso_max_dist']
    _worker_state['combined_lookup'] = config_dict['combined_lookup']
    _worker_state['ref_grna_to_target'] = config_dict['ref_grna_to_target']
    _worker_state['grna_to_ibars'] = config_dict['grna_to_ibars']
    _worker_state['max_lev_ibar'] = config_dict['max_lev_ibar']
    _worker_state['grna_name_to_position'] = config_dict['grna_name_to_position']
    _worker_state['single_guide_mode'] = config_dict.get('single_guide_mode', False)
    _worker_state['seq_id_to_names'] = config_dict['seq_id_to_names']
    _worker_state['short_read_mode'] = config_dict.get('short_read_mode', False)
    _worker_state['save_failed_reads'] = config_dict.get('save_failed_reads', False)

    # Build matchers ONCE per worker
    _worker_state['grna_matcher'] = NeighborhoodMatcher(
        config_dict['grna_seqs'], config_dict['grna_names'],
        max_dist=config_dict['max_lev_grna']
    )
    if not config_dict.get('design_contract') and config_dict['ibar_length'] > 0:
        _worker_state['ibar_matcher'] = NeighborhoodMatcher(
            config_dict['ibar_seqs'], config_dict['ibar_names'],
            max_dist=config_dict['max_lev_ibar']
        )
    else:
        _worker_state['ibar_matcher'] = None

    # Build truncated matchers for short-read single-guide mode
    _worker_state['truncated_matchers'] = {}
    if config_dict.get('short_read_mode', False):
        truncated_lengths = config_dict.get('truncated_lengths', [])
        grna_seqs_full = config_dict['grna_seqs']
        grna_names_full = config_dict['grna_names']
        max_lev = config_dict['max_lev_grna']

        for L in truncated_lengths:
            # Skip lengths that match a full gRNA length (already handled by primary matcher)
            if L in config_dict['grna_lengths']:
                continue
            # Truncate references to length L
            trunc_seqs = [seq[:L] for seq in grna_seqs_full]
            _worker_state['truncated_matchers'][L] = NeighborhoodMatcher(
                trunc_seqs, grna_names_full, max_dist=max_lev
            )


def _resolve_shared_grna(candidate_names, ibar_seq, variant, grna_to_ibars,
                         grna_name_to_position, max_lev_ibar):
    """Resolve a shared gRNA sequence to the correct grna_name.

    When multiple grna_names share the same sequence (e.g., TLNRD1_A_1 at pos1
    and TLNRD1_B_1 at pos2), use iBAR and stem variant to disambiguate.

    Strategy:
      1. Try within-guide iBAR match for each candidate — if exactly one
         candidate has an iBAR match, use it (iBARs don't overlap across positions)
      2. Otherwise, filter by stem variant (position == variant)

    Returns:
        (grna_name, resolution_method, ibar_guide, ibar_dist_guide)
        or (None, 'unresolved', None, None) if cannot resolve
    """
    # 1. Try iBAR-based disambiguation
    ibar_matches = []
    for name in candidate_names:
        guide_ibars = grna_to_ibars.get(name, [])
        ibar_match, ibar_dist = _match_ibar_within_guide(
            ibar_seq, guide_ibars, max_lev_ibar
        )
        if ibar_match is not None:
            ibar_matches.append((name, ibar_match, ibar_dist))

    if len(ibar_matches) == 1:
        name, ibar_match, ibar_dist = ibar_matches[0]
        return name, 'ibar', ibar_match, ibar_dist

    # 2. Fall back to stem variant (only if unambiguous)
    if variant != 'AMBIGUOUS':
        variant_matches = [
            n for n in candidate_names
            if grna_name_to_position.get(n) == variant
        ]
        if len(variant_matches) == 1:
            name = variant_matches[0]
            # Also try iBAR match for the resolved name
            guide_ibars = grna_to_ibars.get(name, [])
            ibar_match, ibar_dist = _match_ibar_within_guide(
                ibar_seq, guide_ibars, max_lev_ibar
            )
            return name, 'variant', ibar_match, ibar_dist

    return None, 'unresolved', None, None


def _process_design_rows(rows):
    state = _worker_state
    contract = state['design_contract']
    resolver = state['design_matcher']
    stats = Counter()
    results, discarded = [], []
    spacers = list(dict.fromkeys(row['grna_seq'] for row in contract['rows']))
    for read_id, bc_wells, bc_bci, umi, sequence in rows:
        sequence = sequence.upper()
        stats['total_reads'] += 1
        start, _ = find_tso(sequence, state['tso_seq'], state['tso_max_dist'])
        observation = extract_unit(sequence, start, contract, state['anchor_max_mm']) if start >= 0 else None
        if observation is None:
            stats['extraction_failed'] += 1
            minimum = min(min(profile['grna_lengths']) + len(profile['stem']) + min(profile['ibar_lengths']) + len(profile['tracr']) for profile in contract['profiles'].values())
            reason = 'unobserved_truncated' if start >= 0 and len(sequence) - start < minimum else 'extraction_failed_or_ambiguous_boundary'
            stats[f'failure_{reason}'] += 1
            if state['save_failed_reads']:
                discarded.append((read_id, bc_wells, bc_bci, umi, sequence, '', '', reason))
            continue
        stats['extracted'] += 1
        spacer_id, distance, unique = state['grna_matcher'].match(observation['grna_seq'])
        resolved = resolver.resolve(observation, spacers[spacer_id]) if unique and spacer_id is not None else None
        if resolved is None:
            stats['grna_no_match'] += 1
            if state['save_failed_reads']:
                discarded.append((read_id, bc_wells, bc_bci, umi, sequence, observation['grna_seq'], observation['ibar_seq'], 'guide_unresolved'))
            continue
        match_type = 'exact' if distance == 0 and resolved['ibar_dist'] == 0 else f"fuzzy_g{distance}_i{resolved['ibar_dist']}"
        stats['exact_match' if match_type == 'exact' else 'fuzzy_match'] += 1
        stats[f"ibar_status_{resolved['ibar_status']}"] += 1
        stats[f"structure_{resolved['structure_status']}"] += 1
        results.append((read_id, bc_wells, bc_bci, umi, resolved['grna_name'], observation['grna_seq'],
                        observation['ibar_seq'], resolved['target_gene'], distance, resolved['ibar_dist'],
                        match_type, observation['variant'], observation['grna_length'], resolved['pair_id'],
                        resolved['target_combo'], resolved['ibar_seq_corrected'], resolved['ibar_decoupled'],
                        *(resolved[column] for column in DESIGN_COLUMNS)))
    return results, discarded, stats


def _process_rows(rows):
    """Process a list of (read_id, bc_wells, bc_bci, umi, r1_seq) tuples.

    Uses worker-local state initialized by _init_worker.
    Returns (results_list, discarded_list, stats_counter).
    """
    ws = _worker_state
    if ws.get('design_contract'):
        return _process_design_rows(rows)
    tso_seq = ws['tso_seq']
    stem1 = ws['stem1']
    stem2 = ws['stem2']
    tracr1 = ws['tracr1']
    tracr2 = ws['tracr2']
    grna_lengths = ws['grna_lengths']
    ibar_length = ws['ibar_length']
    anchor_max_mm = ws['anchor_max_mm']
    tso_max_dist = ws['tso_max_dist']
    grna_matcher = ws['grna_matcher']
    ibar_matcher = ws['ibar_matcher']
    combined_lookup = ws['combined_lookup']
    ref_grna_to_target = ws['ref_grna_to_target']
    grna_to_ibars = ws['grna_to_ibars']
    max_lev_ibar = ws['max_lev_ibar']
    grna_name_to_position = ws['grna_name_to_position']
    single_guide_mode = ws['single_guide_mode']
    seq_id_to_names = ws['seq_id_to_names']
    short_read_mode = ws.get('short_read_mode', False)
    truncated_matchers = ws.get('truncated_matchers', {})
    save_failed = ws.get('save_failed_reads', False)

    stats = Counter()
    results = []
    discarded = []

    for read_id, bc_wells, bc_bci, umi, r1_seq in rows:
        stats['total_reads'] += 1

        # --- Short-read path: single-guide no-iBAR mode ---
        if short_read_mode:
            result = extract_grna_single_short(
                r1_seq, tso_seq, stem1, grna_lengths,
                anchor_max_mm, tso_max_dist,
            )

            if result is None:
                stats['extraction_failed'] += 1
                stats['short_read_discarded'] += 1
                if save_failed:
                    discarded.append((read_id, bc_wells, bc_bci, umi, r1_seq, '', '', 'extraction_failed'))
                continue

            stats['extracted'] += 1
            scenario = result['scenario']
            stats[f'short_read_scenario_{scenario}'] += 1

            grna_seq = result['grna_seq']
            grna_len = result['grna_length']
            stem_truncated = result['stem_truncated']
            grna_truncated = result['grna_truncated']
            tso_offset = result['tso_offset']
            stats[f'tso_offset_{tso_offset}'] += 1
            anchor_combo = result.get('anchor_combo', '?')
            stats[f'anchor_combo_{anchor_combo}'] += 1

            # Select the appropriate matcher based on gRNA length
            if grna_len in [len(s) for s in ws['grna_matcher'].seqs]:
                # Full-length gRNA — use primary matcher
                active_matcher = grna_matcher
            elif grna_len in truncated_matchers:
                # Truncated gRNA — use pre-built truncated matcher
                active_matcher = truncated_matchers[grna_len]
                stats['truncated_matcher_used'] += 1
            else:
                # No matcher for this length (extremely rare/short)
                stats['grna_no_matcher_for_length'] += 1
                if save_failed:
                    discarded.append((read_id, bc_wells, bc_bci, umi, r1_seq, grna_seq, '', 'grna_no_matcher_for_length'))
                continue

            # Fuzzy match gRNA
            grna_match_raw, grna_dist, grna_seq_unique = active_matcher.match(grna_seq)

            if grna_match_raw is None:
                if grna_dist is not None:
                    stats['grna_ambiguous'] += 1
                    if grna_truncated:
                        stats['grna_ambiguous_truncated'] += 1
                    if save_failed:
                        discarded.append((read_id, bc_wells, bc_bci, umi, r1_seq, grna_seq, '', 'grna_ambiguous'))
                else:
                    stats['grna_no_match'] += 1
                    if save_failed:
                        discarded.append((read_id, bc_wells, bc_bci, umi, r1_seq, grna_seq, '', 'grna_no_match'))
                continue

            if isinstance(grna_match_raw, list):
                # Truncated sequences lost uniqueness and collides with multiple full sequences
                stats['grna_ambiguous'] += 1
                if grna_truncated:
                    stats['grna_ambiguous_truncated'] += 1
                if save_failed:
                    discarded.append((read_id, bc_wells, bc_bci, umi, r1_seq, grna_seq, '', 'grna_ambiguous'))
                continue

            # Map seq_id back to grna_names
            candidate_names = seq_id_to_names[grna_match_raw]
            if len(candidate_names) == 1:
                grna_name = candidate_names[0]
            else:
                # In single-guide mode, shared gRNA seqs are ambiguous
                stats['grna_ambiguous'] += 1
                if save_failed:
                    discarded.append((read_id, bc_wells, bc_bci, umi, r1_seq, grna_seq, '', 'grna_ambiguous_shared'))
                continue

            target_gene = ref_grna_to_target.get(grna_name)
            if target_gene is None:
                stats['grna_name_not_in_ref'] += 1
                if save_failed:
                    discarded.append((read_id, bc_wells, bc_bci, umi, r1_seq, grna_seq, '', 'grna_name_not_in_ref'))
                continue

            match_type = f"grna_only_d{grna_dist}"
            if grna_truncated:
                match_type = f"grna_truncated_d{grna_dist}"
            elif stem_truncated:
                match_type = f"grna_partial_stem_d{grna_dist}"
            stats['grna_only_match'] += 1

            results.append((
                read_id, bc_wells, bc_bci, umi, grna_name, grna_seq,
                '', target_gene, grna_dist, -1, match_type,
                1, grna_len, '', '',
                '', False,
                stem_truncated, grna_truncated, grna_len,
            ))
            continue

        # --- Standard path (dual-guide or single-guide with full reads) ---

        # Extract gRNA + iBAR (includes TSO finding + combined stem+tracr scoring)
        result = extract_grna_ibar(
            r1_seq, tso_seq, stem1, stem2,
            grna_lengths, ibar_length,
            tracr1, tracr2, anchor_max_mm, tso_max_dist
        )

        if result is None:
            stats['extraction_failed'] += 1
            if save_failed:
                discarded.append((read_id, bc_wells, bc_bci, umi, r1_seq, '', '', 'extraction_failed'))
            continue

        stats['extracted'] += 1
        grna_seq = result['grna_seq']
        ibar_seq = result['ibar_seq']
        variant = result['variant']
        grna_len = result['grna_length']
        anchor_combo = result.get('anchor_combo', '?')

        # Track anchor combo stats
        stats[f'anchor_combo_{anchor_combo}'] += 1

        # Track TSO offset distribution
        tso_offset = result['tso_offset']
        stats[f'tso_offset_{tso_offset}'] += 1

        # --- Matching ---

        # 1. Exact combined (gRNA + iBAR) lookup
        combined_key = (grna_seq, ibar_seq)
        if combined_key in combined_lookup:
            grna_name, target_gene, pair_id, target_combo = combined_lookup[combined_key]
            stats['exact_match'] += 1
            results.append((
                read_id, bc_wells, bc_bci, umi, grna_name, grna_seq,
                ibar_seq, target_gene, 0, 0, 'exact', variant, grna_len,
                pair_id, target_combo,
                ibar_seq, False,  # ibar_seq_corrected, ibar_decoupled
            ))
            continue

        # 2. Fuzzy match gRNA
        grna_match_raw, grna_dist, grna_seq_unique = grna_matcher.match(grna_seq)

        if grna_match_raw is None:
            if grna_dist is not None:
                stats['grna_ambiguous'] += 1
                if save_failed:
                    discarded.append((read_id, bc_wells, bc_bci, umi, r1_seq, grna_seq, ibar_seq, 'grna_ambiguous'))
            else:
                stats['grna_no_match'] += 1
                if save_failed:
                    discarded.append((read_id, bc_wells, bc_bci, umi, r1_seq, grna_seq, ibar_seq, 'grna_no_match'))
            continue

        # Map seq_id back to grna_names
        candidate_names = seq_id_to_names[grna_match_raw]
        if len(candidate_names) == 1:
            grna_match_result = candidate_names[0]
            grna_unique = True
        else:
            grna_match_result = candidate_names
            grna_unique = False

        # 3. Handle shared gRNA sequences (position-ambiguous)
        if isinstance(grna_match_result, list):
            if single_guide_mode:
                stats['grna_ambiguous'] += 1
                if save_failed:
                    discarded.append((read_id, bc_wells, bc_bci, umi, r1_seq, grna_seq, ibar_seq, 'grna_ambiguous_shared'))
                continue
            
            # Multiple grna_names share this sequence — disambiguate (dual mode)
            stats['shared_grna_hit'] += 1
            grna_name, resolution, ibar_guide, ibar_dist_guide = _resolve_shared_grna(
                grna_match_result, ibar_seq, variant,
                grna_to_ibars, grna_name_to_position, max_lev_ibar
            )
            if grna_name is None:
                stats['shared_grna_unresolved'] += 1
                if save_failed:
                    discarded.append((read_id, bc_wells, bc_bci, umi, r1_seq, grna_seq, ibar_seq, 'shared_grna_unresolved'))
                continue
            stats[f'shared_grna_{resolution}_resolved'] += 1
            grna_unique = True  # Now resolved
        else:
            grna_name = grna_match_result
            ibar_guide = None
            ibar_dist_guide = None

        if not grna_unique:
            # Sequence-ambiguous (different seqs at same distance)
            stats['grna_ambiguous'] += 1
            if save_failed:
                discarded.append((read_id, bc_wells, bc_bci, umi, r1_seq, grna_seq, ibar_seq, 'grna_ambiguous_seq'))
            continue

        # 4. Fuzzy match iBAR — dual matching (global + within-guide)
        if ibar_length == 0:
            target_gene = ref_grna_to_target.get(grna_name)
            if target_gene is None:
                stats['grna_name_not_in_ref'] += 1
                if save_failed:
                    discarded.append((read_id, bc_wells, bc_bci, umi, r1_seq, grna_seq, '', 'grna_name_not_in_ref'))
                continue
            
            match_type = f"grna_only_d{grna_dist}"
            stats['grna_only_match'] += 1
            results.append((
                read_id, bc_wells, bc_bci, umi, grna_name, grna_seq,
                '', target_gene, grna_dist, -1, match_type,
                variant, grna_len, '', '',
                '', False,
            ))
            continue
        ibar_global, ibar_dist_global, ibar_unique_global = ibar_matcher.match(ibar_seq)

        target_gene = ref_grna_to_target.get(grna_name)
        if target_gene is None:
            stats['grna_name_not_in_ref'] += 1
            if save_failed:
                discarded.append((read_id, bc_wells, bc_bci, umi, r1_seq, grna_seq, ibar_seq, 'grna_name_not_in_ref'))
            continue

        # Within-guide iBAR match (skip if already done during shared-gRNA resolution)
        if ibar_guide is None:
            guide_ibars = grna_to_ibars.get(grna_name, [])
            ibar_guide, ibar_dist_guide = _match_ibar_within_guide(
                ibar_seq, guide_ibars, max_lev_ibar
            )

        # Determine corrected iBAR and decoupling status
        ibar_decoupled = False
        if ibar_guide is not None:
            # Within-guide match found — this is the definitive correction
            ibar_seq_corrected = ibar_guide
            ibar_dist = ibar_dist_guide
            if ibar_global is not None and ibar_global != ibar_guide:
                # Global match disagrees with within-guide → decoupling
                ibar_decoupled = True
                stats['ibar_decoupled'] += 1
        elif ibar_global is not None and ibar_unique_global:
            # No within-guide match but global match found → potential decoupling
            ibar_seq_corrected = ibar_global
            ibar_dist = ibar_dist_global
            ibar_decoupled = True
            stats['ibar_decoupled'] += 1
        else:
            ibar_seq_corrected = ''
            ibar_dist = None

        if ibar_dist is not None:
            match_type = f"fuzzy_g{grna_dist}_i{ibar_dist}"
            stats['fuzzy_match'] += 1
        else:
            match_type = f"grna_only_d{grna_dist}"
            ibar_dist = -1
            stats['grna_only_match'] += 1

        results.append((
            read_id, bc_wells, bc_bci, umi, grna_name, grna_seq,
            ibar_seq, target_gene, grna_dist, ibar_dist, match_type,
            variant, grna_len, '', '',
            ibar_seq_corrected, ibar_decoupled,
        ))

    return results, discarded, stats


# Keep backward-compatible wrapper for single-process mode
def process_chunk(args):
    """Process a chunk of barcode TSV data (single-process mode)."""
    chunk_df, config_dict = args

    # Initialize worker state if not already done
    if not _worker_state:
        _init_worker(config_dict)

    rows = list(zip(
        chunk_df['read_id'], chunk_df['bc_wells'], chunk_df['bc_bci'],
        chunk_df['umi'], chunk_df['r1_seq'],
    ))
    results, discarded, stats = _process_rows(rows)
    return results, discarded, stats


def _chunk_generator(barcode_tsv_path, chunk_size):
    """Yield lightweight row-tuples from a (possibly gzipped) TSV, in chunks."""
    for chunk_df in pd.read_csv(
        barcode_tsv_path, sep='\t', chunksize=chunk_size,
        usecols=['read_id', 'bc_wells', 'bc_bci', 'umi', 'r1_seq'],
        dtype=str, keep_default_na=False,
    ):
        rows = list(zip(
            chunk_df['read_id'].values,
            chunk_df['bc_wells'].values,
            chunk_df['bc_bci'].values,
            chunk_df['umi'].values,
            chunk_df['r1_seq'].values,
        ))
        yield rows


def _bounded_results(executor, function, items, max_pending):
    """Bound queued work AND completed results while the parent writes output."""
    if max_pending < 1:
        raise ValueError('max_pending must be positive')
    iterator = iter(items)
    sentinel = object()
    pending = set()
    for _ in range(max_pending):
        item = next(iterator, sentinel)
        if item is sentinel:
            break
        pending.add(executor.submit(function, item))
    while pending:
        done, pending = wait(pending, return_when=FIRST_COMPLETED)
        while done:
            future = done.pop()
            yield future.result()
            item = next(iterator, sentinel)
            if item is not sentinel:
                pending.add(executor.submit(function, item))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def load_reference(ref_path):
    """Load gRNA/iBAR reference lookup table."""
    df = pd.read_csv(ref_path, comment='#')
    required_cols = ['grna_name', 'grna_seq', 'ibar_seq', 'target_gene']
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in reference: {missing}")
    print(f"  Loaded reference: {len(df)} entries, "
          f"{df['target_gene'].nunique()} target genes, "
          f"{df['grna_seq'].nunique()} unique gRNA seqs")
    return df


OUTPUT_COLS = [
    'read_id', 'bc_wells', 'bc_bci', 'umi', 'grna_name', 'grna_seq',
    'ibar_seq', 'target_gene', 'grna_dist', 'ibar_dist',
    'match_type', 'variant', 'grna_length', 'pair_id', 'target_combo',
    'ibar_seq_corrected', 'ibar_decoupled',
]

# Extended columns for single-guide short-read mode
OUTPUT_COLS_SHORT = OUTPUT_COLS + [
    'stem_truncated', 'grna_truncated', 'observed_grna_length',
]


FAILED_COLS = [
    'read_id', 'bc_wells', 'bc_bci', 'umi', 'r1_seq',
    'grna_seq', 'ibar_seq', 'failure_reason',
]


def process_library(barcode_tsv_path, output_path, config_dict, n_processes=1,
                    chunk_size=100_000):
    """Process one library with optional multiprocessing.

    Optimizations over the naive approach:
      - Worker state (matchers, lookup dicts) initialized ONCE per process
        via the worker initializer, not pickled per chunk.
      - Chunks are streamed as lightweight tuples, not full DataFrames.
      - Results are written incrementally to avoid accumulating everything
        in memory.
      - Progress output in both single and multi-process modes.
    """
    t_start = time.time()
    if n_processes < 1 or chunk_size < 1:
        raise ValueError('n_processes and chunk_size must be positive')
    print(f"  Matching: {n_processes} worker(s), {chunk_size:,} reads/chunk; "
          f"at most {2 * n_processes} chunks in flight", flush=True)

    total_stats = Counter()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    n_results_total = 0
    n_discarded_total = 0
    short_read_mode = config_dict.get('short_read_mode', False)
    save_failed = config_dict.get('save_failed_reads', False)
    output_cols = OUTPUT_COLS_SHORT if short_read_mode else OUTPUT_COLS
    if config_dict.get('design_contract'):
        output_cols = OUTPUT_COLS + DESIGN_COLUMNS

    # Failed reads output path (companion to main output)
    discarded_path = None
    if save_failed:
        discarded_path = output_path.replace('_grna_matches.tsv', '_grna_failed.tsv')
        if discarded_path == output_path:
            # Fallback if naming pattern didn't match
            discarded_path = output_path.replace('.tsv', '_failed.tsv')

    # Helper to open output file (gzip or plain)
    import gzip as _gzip
    def _open_out(path):
        if path.endswith('.gz'):
            return _gzip.open(path, 'wt')
        return open(path, 'w')

    if n_processes == 1:
        # Single-process: use _init_worker once, then stream chunks
        _worker_state.clear()
        _init_worker(config_dict)

        fout = _open_out(output_path)
        fout_failed = _open_out(discarded_path) if save_failed else None
        try:
            fout.write('\t'.join(output_cols) + '\n')
            if fout_failed:
                fout_failed.write('\t'.join(FAILED_COLS) + '\n')

            for i, rows in enumerate(_chunk_generator(barcode_tsv_path, chunk_size)):
                results, discarded, stats = _process_rows(rows)
                total_stats += stats

                # Write results immediately
                for r in results:
                    fout.write('\t'.join(str(x) for x in r) + '\n')
                n_results_total += len(results)

                # Write failed reads
                if fout_failed and discarded:
                    for d in discarded:
                        fout_failed.write('\t'.join(str(x) for x in d) + '\n')
                    n_discarded_total += len(discarded)

                processed = total_stats['total_reads']
                elapsed = time.time() - t_start
                rate = processed / elapsed if elapsed > 0 else 0
                matched = total_stats['exact_match'] + total_stats['fuzzy_match'] + total_stats['grna_only_match']
                match_pct = matched / processed * 100 if processed > 0 else 0
                print(f"    Chunk {i+1}: {processed:,} reads "
                      f"({rate:,.0f}/s, {match_pct:.1f}% matched)", flush=True)
        finally:
            fout.close()
            if fout_failed:
                fout_failed.close()
    else:
        # Bound pending chunks so a fast reader/worker cannot outrun gzip output.
        pool = ProcessPoolExecutor(max_workers=n_processes,
                                   initializer=_init_worker,
                                   initargs=(config_dict,))
        fout = _open_out(output_path)
        fout_failed = _open_out(discarded_path) if save_failed else None
        try:
            fout.write('\t'.join(output_cols) + '\n')
            if fout_failed:
                fout_failed.write('\t'.join(FAILED_COLS) + '\n')

            chunks_done = 0
            chunk_iter = _chunk_generator(barcode_tsv_path, chunk_size)

            for results, discarded, stats in _bounded_results(
                pool, _process_rows, chunk_iter, max_pending=2 * n_processes
            ):
                total_stats += stats
                chunks_done += 1

                for r in results:
                    fout.write('\t'.join(str(x) for x in r) + '\n')
                n_results_total += len(results)

                # Write failed reads
                if fout_failed and discarded:
                    for d in discarded:
                        fout_failed.write('\t'.join(str(x) for x in d) + '\n')
                    n_discarded_total += len(discarded)

                # Report the first completed chunk, then every 10 chunks.
                if chunks_done == 1 or chunks_done % 10 == 0:
                    processed = total_stats['total_reads']
                    elapsed = time.time() - t_start
                    rate = processed / elapsed if elapsed > 0 else 0
                    matched = (total_stats['exact_match'] +
                               total_stats['fuzzy_match'] +
                               total_stats['grna_only_match'])
                    match_pct = matched / processed * 100 if processed > 0 else 0
                    print(f"    Chunks {chunks_done}: "
                          f"{processed:,} reads ({rate:,.0f}/s, "
                          f"{match_pct:.1f}% matched)", flush=True)
        finally:
            fout.close()
            if fout_failed:
                fout_failed.close()
            pool.shutdown(wait=True, cancel_futures=True)

    elapsed = time.time() - t_start
    print(f"  Completed in {elapsed:.1f}s "
          f"({total_stats['total_reads']/elapsed:,.0f} reads/s)")
    print(f"  Wrote {n_results_total:,} matched reads to {output_path}")
    if save_failed:
        print(f"  Wrote {n_discarded_total:,} failed reads to {discarded_path}")

    return total_stats


def main():
    parser = argparse.ArgumentParser(
        description="Match gRNA and iBAR sequences from R1 reads"
    )
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument("--library", type=str, default=None)
    parser.add_argument("--n_processes", type=int, default=None)
    parser.add_argument("--chunk_size", type=int, default=100_000)
    parser.add_argument(
        '--preflight-only', action='store_true',
        help='Check a bounded read sample and write diagnostics without production matching.',
    )
    parser.add_argument(
        '--preflight-fastq', action='store_true',
        help='With --preflight-only, sample raw gRNA R1 FASTQs before Step 2 exists.',
    )
    parser.add_argument(
        "--barcode-input", default=None,
        help=("Explicit barcode-evidence TSV for one library. Production Step 1c "
              "passes the compact corrected-barcode evidence emitted by Step 1a."),
    )
    parser.add_argument(
        "--output-namespace", default=None,
        help=("Output root for this source. Relative paths are below paths.output_dir; "
              "default keeps the historical output layout."),
    )
    args = parser.parse_args()
    if args.preflight_fastq and (not args.preflight_only or args.barcode_input):
        parser.error('--preflight-fastq requires --preflight-only and cannot use --barcode-input')

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    ref_df = load_reference(cfg['paths']['grna_reference'])
    if enabled(cfg):
        ref_df = ref_df.fillna('')
    design_contract = compile_design(ref_df, cfg) if enabled(cfg) else None
    frozen_contract = None
    if design_contract and cfg.get('guide_design', {}).get('resolved_contract'):
        frozen_path = Path(cfg['guide_design']['resolved_contract']).expanduser()
        if not frozen_path.is_absolute():
            frozen_path = Path(args.config).resolve().parent / frozen_path
        with frozen_path.open() as handle:
            frozen_contract = yaml.safe_load(handle)
        if frozen_contract.get('input_fingerprint', frozen_contract.get('fingerprint')) != design_contract['fingerprint']:
            raise ValueError('resolved_contract is stale for this reference/settings; rerun the raw pilot')
        if compile_design(frozen_contract['rows'], cfg)['fingerprint'] != frozen_contract.get('fingerprint'):
            raise ValueError('resolved_contract content does not match its fingerprint')
        frozen_contract['capabilities']['infer_missing_partner'] = design_contract['capabilities']['infer_missing_partner']
    if design_contract:
        emit_warnings(design_contract)
        ref_df = pd.DataFrame(design_contract['rows'])

    gs = cfg['grna_structure']
    gm = cfg['grna_matching']
    preflight_enabled, preflight_sample_reads = preflight_options(gm)
    n_processes = (args.n_processes if args.n_processes is not None
                   else cfg['barcode_extraction'].get('n_processes', 1))
    if n_processes < 1 or args.chunk_size < 1:
        parser.error('--n_processes and --chunk_size must be positive')

    single_guide_mode = gs.get('single_guide_mode', False)
    tso_seq = gs['fixed_seq_1']
    stem1 = gs['fixed_seq_2']
    stem2 = gs.get('fixed_seq_3', None) if not single_guide_mode else None
    tracr1 = gs.get('fixed_seq_4', None)
    tracr2 = gs.get('fixed_seq_5', None) if not single_guide_mode else None
    grna_lengths = gs.get('grna_lengths', design_contract['grna_lengths'] if design_contract else [20, 21])
    ibar_length = gs.get('ibar_length', 7)
    if design_contract:
        grna_lengths = design_contract['grna_lengths']
        ibar_length = max(length for profile in design_contract['profiles'].values() for length in profile['ibar_lengths'])
    anchor_max_mm = gs['fixed_seq_max_mismatches']
    tso_max_dist = gs.get('tso_max_edit_dist', 3)

    print(f"\n{'='*60}")
    print(f"  Step 3: gRNA + iBAR Matching")
    print(f"{'='*60}")
    print(f"  Alignment backend: {'edlib (C/bit-parallel)' if HAS_EDLIB else 'pure Python (slower — install edlib for ~50x speedup)'}")
    print(f"  R1 structure: variable_offset + TSO({len(tso_seq)}) + gRNA({grna_lengths}) "
          f"+ stem({len(stem1)}) + iBAR({ibar_length}) + tracr({len(tracr1)})")
    print(f"  TSO max edit dist: {tso_max_dist}")
    print(f"  Stem max mismatches: {anchor_max_mm}")
    print(f"  Max Levenshtein (gRNA): {gm['max_lev_dist']}")
    print(f"  Max Levenshtein (iBAR): {gm['max_lev_dist_ibar']}")
    print(f"  Processes: {n_processes}")
    print(f"  Save failed reads: {gm.get('save_failed_reads', False)}")

    # Detect short-read mode: single-guide + no iBAR
    short_read_mode = single_guide_mode and ibar_length == 0 and not design_contract
    if short_read_mode:
        print(f"  ⚡ Short-read mode: ENABLED (single-guide, no iBAR)")
        # One-time warning about partial stem mismatch overrides
        stem_len = len(stem1)
        effective_mm_short = 0 if stem_len <= 2 else min(1, anchor_max_mm)
        if effective_mm_short != anchor_max_mm:
            print(f"  ⚠ WARNING: Partial stem mismatch thresholds differ from config")
            print(f"    Config fixed_seq_max_mismatches = {anchor_max_mm}")
            print(f"    Partial stem ≤2bp: 0 mismatches allowed (overridden)")
            if anchor_max_mm >= 1:
                print(f"    Partial stem >2bp: 1 mismatch allowed")
            else:
                print(f"    Partial stem >2bp: 0 mismatches allowed")
    else:
        print(f"  Short-read mode: disabled")
    print()

    # Build lookups
    combined_lookup = {}
    ref_grna_to_target = {}
    legacy_reference = ref_df.iloc[:0] if design_contract else ref_df
    for _, row in legacy_reference.iterrows():
        combined_lookup[(row['grna_seq'], row['ibar_seq'])] = (
            row['grna_name'], row['target_gene'],
            row.get('pair_id', ''), row.get('target_combo', ''),
        )
        ref_grna_to_target[row['grna_name']] = row['target_gene']

    # Build per-gRNA iBAR lookup for within-guide matching
    grna_to_ibars = defaultdict(list)
    for _, row in legacy_reference.iterrows():
        grna_to_ibars[row['grna_name']].append(row['ibar_seq'])
    grna_to_ibars = dict(grna_to_ibars)  # convert for pickle efficiency

    # Build grna_name → position mapping (for stem-variant disambiguation)
    grna_name_to_position = {}
    for _, row in legacy_reference.iterrows():
        grna_name_to_position[row['grna_name']] = row['position']

    # Deduplicate gRNA sequences for NeighborhoodMatcher efficiency
    # Map each unique seq to an internal seq_id, and track all grna_names it belongs to.
    grna_dedup = ref_df.drop_duplicates('grna_seq') if design_contract else ref_df.drop_duplicates('grna_name')
    grna_seq_to_names = defaultdict(list)
    for _, row in grna_dedup.iterrows():
        grna_seq_to_names[row['grna_seq']].append(row['grna_name'])
        
    unique_grna_seqs = list(grna_seq_to_names.keys())
    unique_seq_ids = list(range(len(unique_grna_seqs)))
    seq_id_to_names = {i: names for i, names in zip(unique_seq_ids, grna_seq_to_names.values())}

    ibar_seqs = [] if design_contract else ref_df['ibar_seq'].unique().tolist()

    # Detect shared gRNA sequences
    grna_seq_counts = grna_dedup.groupby('grna_seq')['grna_name'].count()
    shared_grna_seqs = grna_seq_counts[grna_seq_counts > 1]

    print("  Reference lookup sizes (matchers initialize after preflight)...")
    # Actual matchers are built once per worker, after preflight. Avoid a spare
    # parent copy, which becomes substantial for large reference libraries.
    n_unique_grna = len(unique_grna_seqs)
    print(f"    gRNA: {n_unique_grna} unique sequences ({len(grna_dedup)} initial entries)")
    if len(shared_grna_seqs) > 0:
        n_shared = len(shared_grna_seqs)
        n_affected = ref_df[ref_df['grna_seq'].isin(shared_grna_seqs.index)]['grna_name'].nunique()
        print(f"    ⚠ {n_shared} gRNA sequences shared across positions "
              f"({n_affected} grna_names affected — using iBAR+variant disambiguation)")
    if design_contract:
        capabilities = design_contract['capabilities']
        print(f"    Design: {capabilities['units_total']} observable units; "
              f"{capabilities['units_identifying_pair']} identify a pair independently")
    else:
        print(f"    iBAR: {len(ibar_seqs)} unique sequences")
    if grna_to_ibars:
        print(f"    Per-gRNA iBAR groups: {len(grna_to_ibars)} "
              f"(avg {sum(len(v) for v in grna_to_ibars.values())/len(grna_to_ibars):.1f} iBARs/gRNA)")

    config_dict = {
        'design_contract': design_contract,
        'tso_seq': tso_seq, 'stem1': stem1, 'stem2': stem2,
        'tracr1': tracr1, 'tracr2': tracr2,
        'grna_lengths': grna_lengths, 'ibar_length': ibar_length,
        'anchor_max_mm': anchor_max_mm, 'tso_max_dist': tso_max_dist,
        'max_lev_grna': gm['max_lev_dist'],
        'max_lev_ibar': gm['max_lev_dist_ibar'],
        'grna_seqs': unique_grna_seqs, 'grna_names': unique_seq_ids,
        'seq_id_to_names': seq_id_to_names,
        'ibar_seqs': ibar_seqs, 'ibar_names': ibar_seqs,
        'combined_lookup': combined_lookup,
        'ref_grna_to_target': ref_grna_to_target,
        'grna_to_ibars': grna_to_ibars,
        'grna_name_to_position': grna_name_to_position,
        'single_guide_mode': single_guide_mode,
        'short_read_mode': short_read_mode,
        'truncated_lengths': [],  # populated per-library via sampling
        'save_failed_reads': gm.get('save_failed_reads', False),
    }

    libraries = cfg['libraries']
    if args.barcode_input and not args.library and len(libraries) != 1:
        print("ERROR: --barcode-input requires --library when config has multiple libraries")
        sys.exit(1)
    if args.library:
        libraries = [lib for lib in libraries if lib['sublibrary_name'] == args.library]
        if not libraries:
            print(f"ERROR: Library '{args.library}' not found"); sys.exit(1)

    output_dir = str(resolve_output_root(
        cfg['paths']['output_dir'], args.output_namespace
    ))
    all_stats = {}

    for lib in libraries:
        lib_name = lib['sublibrary_name']
        # Support both .tsv.gz (preferred) and .tsv (legacy) input
        if args.barcode_input:
            barcode_tsv = str(Path(args.barcode_input).expanduser().resolve())
        else:
            barcode_tsv = os.path.join(output_dir, 'grna_barcodes', f"{lib_name}_barcodes.tsv.gz")
            if not os.path.exists(barcode_tsv):
                barcode_tsv = os.path.join(output_dir, 'grna_barcodes', f"{lib_name}_barcodes.tsv")
        if not args.preflight_fastq and not os.path.exists(barcode_tsv):
            if args.preflight_only:
                raise FileNotFoundError(f'{barcode_tsv} not found; use --preflight-fastq to check raw R1 before Step 2.')
            print(f"\n  WARNING: {barcode_tsv} not found — run Step 2 first!")
            continue

        output_path = os.path.join(output_dir, 'grna_matches', f"{lib_name}_grna_matches.tsv.gz")
        print(f"\n--- Library: {lib_name} ---")

        if preflight_enabled or args.preflight_only:
            fastq_paths = None
            if args.preflight_fastq:
                fastq_paths = sorted(Path(cfg['paths']['fastq_dir']).glob(
                    f"*{lib['grna_prefix']}*_R1_*.fastq.gz"
                ))
                if not fastq_paths:
                    raise FileNotFoundError(f"No gRNA R1 FASTQs found for {lib_name}: {lib['grna_prefix']}")
            if design_contract:
                preflight, effective_contract = check_design(
                    barcode_tsv, design_contract, cfg, find_tso,
                    sample_reads=preflight_sample_reads, fastq_paths=fastq_paths,
                )
                config_dict['design_contract'] = effective_contract
                if frozen_contract:
                    for position, orientation in preflight['orientation'].items():
                        previous = frozen_contract.get('resolved_orientation', {}).get(position, 'as_reference')
                        if orientation in ('as_reference', 'reverse_complement') and previous in ('as_reference', 'reverse_complement') and previous != orientation:
                            raise ValueError('Pilot orientation conflicts with the frozen design contract')
                    config_dict['design_contract'] = frozen_contract
            else:
                preflight = check_ibar_reference(
                    barcode_tsv, ref_df[['grna_seq', 'ibar_seq']].itertuples(index=False, name=None),
                    config_dict, extract_grna_ibar, sample_reads=preflight_sample_reads,
                    fastq_paths=fastq_paths,
                )
            preflight.update(library=lib_name, reference_path=cfg['paths']['grna_reference'],
                             configured_ibar_length=ibar_length)
            print_preflight(preflight, lib_name)
            preflight_path = Path(output_dir) / 'grna_matches' / f'{lib_name}_preflight.yaml'
            preflight_path.parent.mkdir(parents=True, exist_ok=True)
            with preflight_path.open('w') as handle:
                yaml.safe_dump(preflight, handle, sort_keys=False)
        else:
            print(f"  iBAR preflight [{lib_name}]: disabled by config", flush=True)
            config_dict['design_contract'] = frozen_contract or design_contract
            if design_contract and not frozen_contract and cfg.get('guide_design', {}).get('ibar_orientation') == 'auto':
                raise ValueError('ibar_orientation: auto requires enabled preflight')
        if design_contract:
            contract_path = Path(output_dir) / 'grna_matches' / f'{lib_name}_design.yaml'
            contract_path.parent.mkdir(parents=True, exist_ok=True)
            with contract_path.open('w') as handle:
                yaml.safe_dump(config_dict['design_contract'], handle, sort_keys=False)
        if args.preflight_only:
            continue

        # Run 100k-read sampling to detect truncated gRNA lengths (short-read mode)
        if short_read_mode:
            print(f"  Sampling first 100,000 reads to detect gRNA lengths...")
            detected_lengths, length_counts, scenario_counts = sample_truncated_lengths(
                barcode_tsv, tso_seq, stem1, grna_lengths,
                anchor_max_mm, tso_max_dist, n_sample=100_000,
            )
            config_dict['truncated_lengths'] = detected_lengths

            # Log detected vs expected lengths
            print(f"    Expected gRNA lengths (config): {grna_lengths}")
            print(f"    Detected gRNA lengths (sampled): {detected_lengths}")
            print(f"    Length distribution:")
            for L in sorted(length_counts.keys()):
                label = f"{L}bp"
                if L in grna_lengths:
                    label += " (expected)"
                elif L < min(grna_lengths):
                    label += " (truncated)"
                print(f"      {label}: {length_counts[L]:,} reads")
            print(f"    Scenario distribution (sampled):")
            for s in sorted(scenario_counts.keys(), key=str):
                print(f"      Scenario {s}: {scenario_counts[s]:,}")
            # Clear worker state so it rebuilds with new truncated_lengths
            _worker_state.clear()

        stats = process_library(
            barcode_tsv, output_path, config_dict,
            n_processes=n_processes, chunk_size=args.chunk_size
        )
        all_stats[lib_name] = dict(stats)

        total = stats['total_reads']
        exact = stats.get('exact_match', 0)
        fuzzy = stats.get('fuzzy_match', 0)
        grna_only = stats.get('grna_only_match', 0)
        total_matched = exact + fuzzy + grna_only

        print(f"\n  Summary for {lib_name}:")
        print(f"    Total reads:        {total:>12,}")
        print(f"    Extraction failed:  {stats.get('extraction_failed', 0):>12,}")
        print(f"    Extracted:          {stats.get('extracted', 0):>12,} "
              f"({stats.get('extracted', 0)/total*100:.1f}%)" if total else "")
        print(f"    Exact match:        {exact:>12,}")
        print(f"    Fuzzy match:        {fuzzy:>12,}")
        print(f"    gRNA-only match:    {grna_only:>12,}")
        print(f"    Total matched:      {total_matched:>12,} "
              f"({total_matched/total*100:.1f}%)" if total else "")
        print(f"    gRNA ambiguous:     {stats.get('grna_ambiguous', 0):>12,}")
        print(f"    gRNA no match:      {stats.get('grna_no_match', 0):>12,}")
        if stats.get('ibar_decoupled', 0) > 0:
            print(f"    iBAR decoupled:     {stats.get('ibar_decoupled', 0):>12,}")

        # Short-read mode specific stats
        if short_read_mode:
            print(f"    --- Short-read stats ---")
            for s in [1, 2, 3]:
                key = f'short_read_scenario_{s}'
                count = stats.get(key, 0)
                labels = {1: 'Full stem', 2: 'Partial stem', 3: 'Truncated gRNA'}
                print(f"    Scenario {s} ({labels[s]}):{count:>8,}")
            if stats.get('short_read_discarded', 0) > 0:
                print(f"    Short-read discarded:{stats['short_read_discarded']:>8,}")
            if stats.get('truncated_matcher_used', 0) > 0:
                print(f"    Truncated matcher:  {stats['truncated_matcher_used']:>8,}")
            if stats.get('grna_no_matcher_for_length', 0) > 0:
                print(f"    No matcher (rare L):{stats['grna_no_matcher_for_length']:>8,}")
            if stats.get('grna_ambiguous_truncated', 0) > 0:
                print(f"    Ambiguous (trunc):  {stats['grna_ambiguous_truncated']:>8,}")

        # TSO offset distribution
        offsets = {k: v for k, v in stats.items() if k.startswith('tso_offset_')}
        if offsets:
            print(f"    TSO offset distribution:")
            for key in sorted(offsets.keys(), key=lambda k: int(k.split('_')[-1])):
                offset = key.split('_')[-1]
                count = offsets[key]
                print(f"      offset {offset:>3s}bp: {count:>10,} "
                      f"({count/stats.get('extracted', 1)*100:.1f}%)")
    if args.preflight_only:
        print('\n  Preflight complete; production matching was not started.', flush=True)
        return
    lib_label = args.library if getattr(args, 'library', None) else "all"
    stats_path = os.path.join(output_dir, 'grna_matches', f'{lib_label}_matching_stats.yaml')
    os.makedirs(os.path.dirname(stats_path), exist_ok=True)
    with open(stats_path, 'w') as f:
        yaml.dump(all_stats, f, default_flow_style=False)

    print(f"\n{'='*60}")
    print(f"  Step 3 complete!")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()

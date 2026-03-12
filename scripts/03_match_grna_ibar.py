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
  │ variable (0-15bp) │ Split_TSO (30bp) │ gRNA │ stem │ iBAR │ tracr│
  └───────────────────────────────────────────────────────────────────┘
  The TSO can start at a variable offset (most commonly ~7bp from R1 start).
  We use semi-global alignment to find TSO position, then extract downstream.

Usage:
  python scripts/03_match_grna_ibar.py --config config.yaml
  python scripts/03_match_grna_ibar.py --config config.yaml --library C1 --n_processes 4
"""

import argparse
import os
import sys
import time
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

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
    r1_seq = r1_seq.upper()
    tso_seq = tso_seq.upper()
    stem1 = stem1.upper()
    stem2 = stem2.upper()
    if tracr1:
        tracr1 = tracr1.upper()
    if tracr2:
        tracr2 = tracr2.upper()

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
        mm_stem2 = hamming_distance(stem_observed, stem2, anchor_max_mm + 1)

        # Compute tracr mismatches if available
        mm_tracr1, mm_tracr2 = 0, 0
        has_tracr = False
        if tracr1 and tracr2 and ibar_end + tracr_len <= len(r1_seq):
            tracr_observed = r1_seq[ibar_end:ibar_end + tracr_len]
            mm_tracr1 = hamming_distance(tracr_observed, tracr1, anchor_max_mm + 1)
            mm_tracr2 = hamming_distance(tracr_observed, tracr2, anchor_max_mm + 1)
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


# ---------------------------------------------------------------------------
# Chunk processing for multiprocessing
# ---------------------------------------------------------------------------

# Module-level worker state — set once per worker via Pool(initializer=)
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

    # Build matchers ONCE per worker
    _worker_state['grna_matcher'] = NeighborhoodMatcher(
        config_dict['grna_seqs'], config_dict['grna_names'],
        max_dist=config_dict['max_lev_grna']
    )
    _worker_state['ibar_matcher'] = NeighborhoodMatcher(
        config_dict['ibar_seqs'], config_dict['ibar_names'],
        max_dist=config_dict['max_lev_ibar']
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


def _process_rows(rows):
    """Process a list of (read_id, bc_wells, bc_bci, umi, r1_seq) tuples.

    Uses worker-local state initialized by _init_worker.
    Returns (results_list, stats_counter).
    """
    ws = _worker_state
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

    stats = Counter()
    results = []

    for read_id, bc_wells, bc_bci, umi, r1_seq in rows:
        stats['total_reads'] += 1

        # Extract gRNA + iBAR (includes TSO finding + combined stem+tracr scoring)
        result = extract_grna_ibar(
            r1_seq, tso_seq, stem1, stem2,
            grna_lengths, ibar_length,
            tracr1, tracr2, anchor_max_mm, tso_max_dist
        )

        if result is None:
            stats['extraction_failed'] += 1
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
        grna_match_result, grna_dist, grna_unique = grna_matcher.match(grna_seq)

        if grna_match_result is None:
            stats['grna_no_match'] += 1
            continue

        # 3. Handle shared gRNA sequences (position-ambiguous)
        if isinstance(grna_match_result, list):
            # Multiple grna_names share this sequence — disambiguate
            stats['shared_grna_hit'] += 1
            grna_name, resolution, ibar_guide, ibar_dist_guide = _resolve_shared_grna(
                grna_match_result, ibar_seq, variant,
                grna_to_ibars, grna_name_to_position, max_lev_ibar
            )
            if grna_name is None:
                stats['shared_grna_unresolved'] += 1
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
            continue

        # 4. Fuzzy match iBAR — dual matching (global + within-guide)
        ibar_global, ibar_dist_global, ibar_unique_global = ibar_matcher.match(ibar_seq)

        target_gene = ref_grna_to_target.get(grna_name)
        if target_gene is None:
            stats['grna_name_not_in_ref'] += 1
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

    return results, stats


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
    return _process_rows(rows)


def _chunk_generator(barcode_tsv_path, chunk_size):
    """Yield lightweight row-tuples from a (possibly gzipped) TSV, in chunks."""
    for chunk_df in pd.read_csv(barcode_tsv_path, sep='\t', chunksize=chunk_size):
        rows = list(zip(
            chunk_df['read_id'].values,
            chunk_df['bc_wells'].values,
            chunk_df['bc_bci'].values,
            chunk_df['umi'].values,
            chunk_df['r1_seq'].values,
        ))
        yield rows


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


def process_library(barcode_tsv_path, output_path, config_dict, n_processes=1,
                    chunk_size=100_000):
    """Process one library with optional multiprocessing.

    Optimizations over the naive approach:
      - Worker state (matchers, lookup dicts) initialized ONCE per process
        via Pool(initializer=), not pickled per chunk.
      - Chunks are streamed as lightweight tuples, not full DataFrames.
      - Results are written incrementally to avoid accumulating everything
        in memory.
      - Progress output in both single and multi-process modes.
    """
    t_start = time.time()

    total_stats = Counter()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    n_results_total = 0

    # Helper to open output file (gzip or plain)
    import gzip as _gzip
    def _open_out():
        if output_path.endswith('.gz'):
            return _gzip.open(output_path, 'wt')
        return open(output_path, 'w')

    if n_processes == 1:
        # Single-process: use _init_worker once, then stream chunks
        _worker_state.clear()
        _init_worker(config_dict)

        with _open_out() as fout:
            fout.write('\t'.join(OUTPUT_COLS) + '\n')

            for i, rows in enumerate(_chunk_generator(barcode_tsv_path, chunk_size)):
                results, stats = _process_rows(rows)
                total_stats += stats

                # Write results immediately
                for r in results:
                    fout.write('\t'.join(str(x) for x in r) + '\n')
                n_results_total += len(results)

                processed = total_stats['total_reads']
                elapsed = time.time() - t_start
                rate = processed / elapsed if elapsed > 0 else 0
                matched = total_stats['exact_match'] + total_stats['fuzzy_match'] + total_stats['grna_only_match']
                match_pct = matched / processed * 100 if processed > 0 else 0
                print(f"    Chunk {i+1}: {processed:,} reads "
                      f"({rate:,.0f}/s, {match_pct:.1f}% matched)")
    else:
        # Multi-process: initialize workers via Pool(initializer=)
        pool = Pool(processes=n_processes,
                    initializer=_init_worker,
                    initargs=(config_dict,))
        fout = _open_out()
        try:
            fout.write('\t'.join(OUTPUT_COLS) + '\n')

            chunks_done = 0
            chunk_iter = _chunk_generator(barcode_tsv_path, chunk_size)

            for results, stats in pool.imap_unordered(_process_rows, chunk_iter):
                total_stats += stats
                chunks_done += 1

                for r in results:
                    fout.write('\t'.join(str(x) for x in r) + '\n')
                n_results_total += len(results)

                # Progress every 10 chunks
                if chunks_done % 10 == 0:
                    processed = total_stats['total_reads']
                    elapsed = time.time() - t_start
                    rate = processed / elapsed if elapsed > 0 else 0
                    matched = (total_stats['exact_match'] +
                               total_stats['fuzzy_match'] +
                               total_stats['grna_only_match'])
                    match_pct = matched / processed * 100 if processed > 0 else 0
                    print(f"    Chunks {chunks_done}: "
                          f"{processed:,} reads ({rate:,.0f}/s, "
                          f"{match_pct:.1f}% matched)")
        finally:
            fout.close()
            pool.close()
            pool.join()

    elapsed = time.time() - t_start
    print(f"  Completed in {elapsed:.1f}s "
          f"({total_stats['total_reads']/elapsed:,.0f} reads/s)")
    print(f"  Wrote {n_results_total:,} matched reads to {output_path}")

    return total_stats


def main():
    parser = argparse.ArgumentParser(
        description="Match gRNA and iBAR sequences from R1 reads"
    )
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument("--library", type=str, default=None)
    parser.add_argument("--n_processes", type=int, default=None)
    parser.add_argument("--chunk_size", type=int, default=100_000)
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    ref_df = load_reference(cfg['paths']['grna_reference'])

    gs = cfg['grna_structure']
    gm = cfg['grna_matching']
    n_processes = args.n_processes or cfg['barcode_extraction'].get('n_processes', 1)

    tso_seq = gs['fixed_seq_1']
    stem1 = gs['fixed_seq_2']
    stem2 = gs['fixed_seq_3']
    tracr1 = gs['fixed_seq_4']
    tracr2 = gs['fixed_seq_5']
    grna_lengths = gs['grna_lengths']
    ibar_length = gs['ibar_length']
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
    print()

    # Build lookups
    combined_lookup = {}
    ref_grna_to_target = {}
    for _, row in ref_df.iterrows():
        combined_lookup[(row['grna_seq'], row['ibar_seq'])] = (
            row['grna_name'], row['target_gene'],
            row.get('pair_id', ''), row.get('target_combo', ''),
        )
        ref_grna_to_target[row['grna_name']] = row['target_gene']

    # Build per-gRNA iBAR lookup for within-guide matching
    grna_to_ibars = defaultdict(list)
    for _, row in ref_df.iterrows():
        grna_to_ibars[row['grna_name']].append(row['ibar_seq'])
    grna_to_ibars = dict(grna_to_ibars)  # convert for pickle efficiency

    # Build grna_name → position mapping (for stem-variant disambiguation)
    grna_name_to_position = {}
    for _, row in ref_df.iterrows():
        grna_name_to_position[row['grna_name']] = row['position']

    # Pass entries deduplicated by grna_name (not grna_seq) so duplicate sequences
    # at different positions are both present and returned as position-ambiguous.
    # Without this, drop_duplicates('grna_seq') silently drops one position.
    grna_dedup = ref_df.drop_duplicates('grna_name')
    grna_seqs = grna_dedup['grna_seq'].tolist()
    grna_names = grna_dedup['grna_name'].tolist()
    ibar_seqs = ref_df['ibar_seq'].unique().tolist()

    # Detect shared gRNA sequences
    grna_seq_counts = ref_df.drop_duplicates('grna_name').groupby('grna_seq')['grna_name'].count()
    shared_grna_seqs = grna_seq_counts[grna_seq_counts > 1]

    print("  Building matchers...")
    grna_matcher = NeighborhoodMatcher(grna_seqs, grna_names, max_dist=gm['max_lev_dist'])
    ibar_matcher = NeighborhoodMatcher(ibar_seqs, ibar_seqs, max_dist=gm['max_lev_dist_ibar'])
    n_unique_grna = ref_df['grna_seq'].nunique()
    print(f"    gRNA: {n_unique_grna} unique sequences ({len(grna_seqs)} entries)")
    if len(shared_grna_seqs) > 0:
        n_shared = len(shared_grna_seqs)
        n_affected = ref_df[ref_df['grna_seq'].isin(shared_grna_seqs.index)]['grna_name'].nunique()
        print(f"    ⚠ {n_shared} gRNA sequences shared across positions "
              f"({n_affected} grna_names affected — using iBAR+variant disambiguation)")
    print(f"    iBAR: {len(ibar_seqs)} unique sequences")
    print(f"    Per-gRNA iBAR groups: {len(grna_to_ibars)} "
          f"(avg {sum(len(v) for v in grna_to_ibars.values())/len(grna_to_ibars):.1f} iBARs/gRNA)")

    config_dict = {
        'tso_seq': tso_seq, 'stem1': stem1, 'stem2': stem2,
        'tracr1': tracr1, 'tracr2': tracr2,
        'grna_lengths': grna_lengths, 'ibar_length': ibar_length,
        'anchor_max_mm': anchor_max_mm, 'tso_max_dist': tso_max_dist,
        'max_lev_grna': gm['max_lev_dist'],
        'max_lev_ibar': gm['max_lev_dist_ibar'],
        'grna_seqs': grna_seqs, 'grna_names': grna_names,
        'ibar_seqs': ibar_seqs, 'ibar_names': ibar_seqs,
        'combined_lookup': combined_lookup,
        'ref_grna_to_target': ref_grna_to_target,
        'grna_to_ibars': grna_to_ibars,
        'grna_name_to_position': grna_name_to_position,
    }

    libraries = cfg['libraries']
    if args.library:
        libraries = [lib for lib in libraries if lib['sublibrary_name'] == args.library]
        if not libraries:
            print(f"ERROR: Library '{args.library}' not found"); sys.exit(1)

    output_dir = cfg['paths']['output_dir']
    all_stats = {}

    for lib in libraries:
        lib_name = lib['sublibrary_name']
        # Support both .tsv.gz (preferred) and .tsv (legacy) input
        barcode_tsv = os.path.join(output_dir, 'grna_barcodes', f"{lib_name}_barcodes.tsv.gz")
        if not os.path.exists(barcode_tsv):
            barcode_tsv = os.path.join(output_dir, 'grna_barcodes', f"{lib_name}_barcodes.tsv")
        if not os.path.exists(barcode_tsv):
            print(f"\n  WARNING: {barcode_tsv} not found — run Step 2 first!")
            continue

        output_path = os.path.join(output_dir, 'grna_matches', f"{lib_name}_grna_matches.tsv.gz")
        print(f"\n--- Library: {lib_name} ---")

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

        # TSO offset distribution
        offsets = {k: v for k, v in stats.items() if k.startswith('tso_offset_')}
        if offsets:
            print(f"    TSO offset distribution:")
            for key in sorted(offsets.keys(), key=lambda k: int(k.split('_')[-1])):
                offset = key.split('_')[-1]
                count = offsets[key]
                print(f"      offset {offset:>3s}bp: {count:>10,} "
                      f"({count/stats.get('extracted', 1)*100:.1f}%)")

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

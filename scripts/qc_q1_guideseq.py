#!/usr/bin/env python3
"""
QC Q1: CSM Guide-Seq Bulk Library Analysis

Processes pre-lentiviral bulk paired-end guide-seq FASTQs to estimate:
  - gRNA-iBAR decoupling rate (c)
  - Pre-transduction recombination rates (r1, r2, r12, del1, del2)
  - Within-target vs across-target swap rates
  - Per-guide and per-combo frequency tables (intact constructs only)

This module is designed to run INDEPENDENTLY of the single-cell pipeline.
Use qc_config_guideseq.yaml as the standalone config.

Read structure:
  R1 → pos1: index_1 | gRNA | CSM_stem_1 | iBAR | CSM_partial_tracr_1
  R2 (rev-comp) → pos2: tRNA_incomplete | gRNA | CSM_stem_2 | iBAR | CSM_partial_tracr_2

Guide matching follows the single-cell pipeline (03_match_grna_ibar.py) exactly:
  1. Exact combined (gRNA, iBAR) lookup
  2. Fuzzy gRNA via NeighborhoodMatcher
  3. Fuzzy iBAR (global + within-guide)
  4. Decoupling detection

Usage (from notebook):
  import scripts.qc_q1_guideseq as q1
  results = q1.run_guideseq_qc(fastq_pairs, ref_df, guideseq_cfg, out_dir)
"""

import gzip
import os
from collections import Counter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

# Try to import edlib for fast alignment
try:
    import edlib
    HAS_EDLIB = True
except ImportError:
    HAS_EDLIB = False

# Import matching utilities from step 3
try:
    from scripts import match_grna_ibar_03 as _step3
except ImportError:
    try:
        # Fallback: try importing as sibling module
        import importlib.util
        import sys
        _step3_path = os.path.join(os.path.dirname(__file__), '03_match_grna_ibar.py')
        _spec = importlib.util.spec_from_file_location("match_grna_ibar_03", _step3_path)
        _step3 = importlib.util.module_from_spec(_spec)
        _spec.loader.exec_module(_step3)
    except Exception:
        raise ImportError(
            "Cannot import 03_match_grna_ibar.py. "
            "Make sure to run from the project root directory."
        )

NeighborhoodMatcher = _step3.NeighborhoodMatcher
hamming_distance = _step3.hamming_distance
_match_ibar_within_guide = _step3._match_ibar_within_guide


# ============================================================================
# Sequence Utilities
# ============================================================================

_COMP = str.maketrans('ACGTacgt', 'TGCAtgca')


def revcomp(seq):
    """Reverse complement a DNA sequence."""
    return seq.translate(_COMP)[::-1]


def find_anchor(seq, anchor, max_mm=2):
    """Find best position of anchor in sequence.

    Uses edlib semi-global alignment (HW mode) if available (C-backed, ~10-50x faster),
    falls back to sliding Hamming distance otherwise.

    Both sequences are uppercased internally for case-insensitive comparison.

    Args:
        seq: read sequence to search in
        anchor: anchor sequence to find
        max_mm: maximum allowed edit distance / mismatches

    Returns:
        (start_position, end_position, edit_distance) or (-1, -1, -1) if not found
    """
    seq = seq.upper()
    anchor = anchor.upper()
    if HAS_EDLIB:
        result = edlib.align(anchor, seq, mode="HW", k=max_mm, task="locations")
        if result['editDistance'] < 0 or not result['locations']:
            return -1, -1, -1
        loc = result['locations'][0]
        return loc[0], loc[1] + 1, result['editDistance']  # +1 for exclusive end

    # Pure Python fallback: sliding Hamming
    anchor_len = len(anchor)
    if len(seq) < anchor_len:
        return -1, -1, -1

    best_pos = -1
    best_mm = max_mm + 1

    for i in range(len(seq) - anchor_len + 1):
        mm = hamming_distance(seq[i:i + anchor_len], anchor, max_mm)
        if mm < best_mm:
            best_mm = mm
            best_pos = i
            if mm == 0:
                break

    if best_mm <= max_mm:
        return best_pos, best_pos + anchor_len, best_mm
    return -1, -1, -1


# ============================================================================
# Read Merging (tRNA_leader-anchored)
# ============================================================================

def merge_reads(r1_seq, r2_seq, cfg, min_overlap=10, max_mm_frac=0.15):
    """Merge R1 + R2rc into a single amplicon using tRNA_leader as anchor.

    The 6bp tRNA_leader is expected in both the 3' tail of R1 and the 5'
    head of R2rc for intact constructs.  Finding it in both gives a
    biologically-informed overlap offset.

    Falls back to generic edlib overlap if tRNA_leader anchoring fails.

    Args:
        r1_seq: R1 read sequence (forward)
        r2_seq: R2 read sequence (original orientation, NOT yet revcomp'd)
        cfg: guideseq config dict (needs tRNA_leader)
        min_overlap: minimum required overlap length
        max_mm_frac: max mismatch fraction in overlap region

    Returns:
        (merged_seq, overlap_len) or (None, 0) if merge fails
    """
    r2rc = revcomp(r2_seq)
    trna_leader = cfg.get('tRNA_leader', '')

    if trna_leader and HAS_EDLIB:
        # -- Anchor-based merge --
        # tRNA_leader should be near the 3' end of R1 (position ~130-140)
        # and near the 5' end of R2rc (position ~0-15)
        leader_mm = cfg.get('tRNA_leader_max_mm', 1)

        # Find in R1 (search last 40bp)
        r1_search_start = max(0, len(r1_seq) - 40)
        r1_tail = r1_seq[r1_search_start:]
        l1_pos, l1_end, l1_mm = find_anchor(r1_tail, trna_leader, leader_mm)

        # Find in R2rc (search first 20bp)
        r2rc_head = r2rc[:20]
        l2_pos, l2_end, l2_mm = find_anchor(r2rc_head, trna_leader, leader_mm)

        if l1_pos >= 0 and l2_pos >= 0:
            # tRNA_leader absolute positions
            r1_leader_pos = r1_search_start + l1_pos
            r2rc_leader_pos = l2_pos

            # Overlap: the leader should be at the same position in the merged construct
            # overlap_start_in_r1 = r1_leader_pos - r2rc_leader_pos
            overlap_start = r1_leader_pos - r2rc_leader_pos
            if overlap_start >= 0:
                overlap_len = len(r1_seq) - overlap_start
                if overlap_len >= min_overlap:
                    # Validate overlap quality
                    r1_tail_ov = r1_seq[overlap_start:]
                    r2rc_head_ov = r2rc[:overlap_len]
                    check_len = min(len(r1_tail_ov), len(r2rc_head_ov))
                    if check_len > 0:
                        mm = sum(a != b for a, b in zip(
                            r1_tail_ov[:check_len], r2rc_head_ov[:check_len]))
                        if mm / check_len <= max_mm_frac:
                            merged = r1_seq[:overlap_start] + r2rc
                            return merged, overlap_len

    # -- Fallback: generic edlib overlap detection --
    if HAS_EDLIB:
        query_len = min(50, len(r2rc))
        query = r2rc[:query_len]
        max_k = int(query_len * max_mm_frac)

        result = edlib.align(query, r1_seq, mode="HW", k=max_k, task="locations")
        if result['editDistance'] >= 0 and result['locations']:
            loc = result['locations'][0]
            overlap_start = loc[0]
            overlap_len = len(r1_seq) - overlap_start
            if overlap_len >= min_overlap:
                merged = r1_seq[:overlap_start] + r2rc
                return merged, overlap_len
    else:
        # Pure Python fallback
        for overlap in range(min(len(r1_seq), len(r2rc)), min_overlap - 1, -1):
            r1_tail = r1_seq[-overlap:]
            r2_head = r2rc[:overlap]
            mm = sum(a != b for a, b in zip(r1_tail, r2_head))
            if mm / overlap <= max_mm_frac:
                merged = r1_seq + r2rc[overlap:]
                return merged, overlap

    return None, 0


NOT_FOUND = "Not found"


# ============================================================================
# Per-Read Extraction Functions
# ============================================================================

AMBIGUOUS = "ambiguous"


def _trim_anchor(anchor_seq, max_anchor_length, side):
    """Trim a long anchor sequence to max_anchor_length.

    Long anchors (tracr ~68bp, tRNA ~72bp) may have sequencing errors
    far from the gRNA. Trimming keeps the high-quality portion.

    Args:
        anchor_seq: full anchor sequence
        max_anchor_length: max length (0 = keep full length)
        side: 'right' → remove 3' end (keep 5', for tracr in R1)
              'left'  → remove 5' end (keep 3', for tRNA in R2rc)

    Returns:
        trimmed anchor
    """
    if max_anchor_length <= 0 or len(anchor_seq) <= max_anchor_length:
        return anchor_seq
    if side == 'right':
        return anchor_seq[:max_anchor_length]   # keep 5' end
    else:  # side == 'left'
        return anchor_seq[-max_anchor_length:]  # keep 3' end


def _extract_guide_unit(seq, left_anchor, stem_1, stem_2,
                        right_anchor_1, right_anchor_2,
                        grna_lengths, ibar_length,
                        left_anchor_max_mm, stem_max_mm, right_anchor_max_mm):
    """Extract one guide unit using left-anchor + combined right-scoring.

    Mirrors Step 3's extract_grna_ibar logic exactly:
      1. Find left_anchor via semi-global alignment → anchor_end
      2. For each gRNA length (longest first):
         - grna_start = anchor_end
         - stem immediately after gRNA → Hamming score stem_1 / stem_2
         - ibar: 7bp after stem
         - right_anchor after ibar → Hamming score
         - Combined score = stem_mm + right_anchor_mm
         - Accept if combined ≤ threshold
      3. Return best match with Step-3-consistent metadata

    Args:
        seq: read sequence (R1 raw or R2rc)
        left_anchor: anchor seq immediately upstream of gRNA (R1_primer or tRNA)
        stem_1, stem_2: stem sequences for variant scoring
        right_anchor_1: first right anchor (tracr_1 for R1, R2_primer_rc for R2)
        right_anchor_2: second right anchor (tracr_2 for R1, None for R2)
        grna_lengths: list of candidate gRNA lengths [20, 21]
        ibar_length: iBAR length (7)
        left_anchor_max_mm: max edit distance for left anchor
        stem_max_mm: max combined stem+right_anchor mismatches
        right_anchor_max_mm: max edit distance when scoring right anchor alone

    Returns:
        dict with: grna_seq, ibar_seq, variant (1 or 2), grna_length,
                   anchor_offset, anchor_dist, stem_mm, right_anchor_mm,
                   anchor_combo
        or None if extraction fails.
    """
    # Normalize case for all sequence comparisons
    seq = seq.upper()
    stem_1 = stem_1.upper()
    stem_2 = stem_2.upper()
    if right_anchor_1:
        right_anchor_1 = right_anchor_1.upper()
    if right_anchor_2:
        right_anchor_2 = right_anchor_2.upper()

    # Step 1: Find left anchor (like TSO in Step 3)
    # (find_anchor already uppercases internally)
    if left_anchor:
        anchor_pos, anchor_end, anchor_dist = find_anchor(
            seq, left_anchor, left_anchor_max_mm)
        if anchor_pos < 0:
            return None
    else:
        # No left anchor — gRNA starts at position 0
        # (used by construct extraction where the split point IS the anchor)
        anchor_pos = 0
        anchor_end = 0
        anchor_dist = 0

    stem_len = len(stem_1)

    # Step 2: Try each gRNA length (longest first, same as Step 3)
    for grna_len in sorted(grna_lengths, reverse=True):
        grna_start = anchor_end
        grna_end = grna_start + grna_len
        stem_start = grna_end
        stem_end = stem_start + stem_len
        ibar_start = stem_end
        ibar_end = ibar_start + ibar_length

        if ibar_end > len(seq):
            continue

        # Score stem variants (Hamming at positionally-determined location)
        stem_observed = seq[stem_start:stem_end]
        mm_stem1 = hamming_distance(stem_observed, stem_1, stem_max_mm + 1)
        mm_stem2 = hamming_distance(stem_observed, stem_2, stem_max_mm + 1)

        # Score right anchors downstream of iBAR
        mm_ra1, mm_ra2 = 0, 0
        has_right_anchor = False

        if right_anchor_1:
            ra1_len = len(right_anchor_1)
            if ibar_end + ra1_len <= len(seq):
                ra_observed = seq[ibar_end:ibar_end + ra1_len]
                mm_ra1 = hamming_distance(
                    ra_observed, right_anchor_1, right_anchor_max_mm + 1)
                has_right_anchor = True

        if right_anchor_2:
            ra2_len = len(right_anchor_2)
            if ibar_end + ra2_len <= len(seq):
                ra_observed = seq[ibar_end:ibar_end + ra2_len]
                mm_ra2 = hamming_distance(
                    ra_observed, right_anchor_2, right_anchor_max_mm + 1)
                has_right_anchor = True

        # Combined scoring (stem + right_anchor), same as Step 3
        score_A = mm_stem1 + mm_ra1   # stem_1 + right_anchor_1
        score_B = mm_stem2 + mm_ra2   # stem_2 + right_anchor_2

        # For R2 where right_anchor_2 is None, only score_A and score_B
        # using the single right_anchor for both stem variants
        if right_anchor_2 is None and right_anchor_1:
            # R2 case: single right anchor, score is stem_mm + same ra_mm
            score_A = mm_stem1 + mm_ra1
            score_B = mm_stem2 + mm_ra1

        # Pick best (A or B), reject if combined exceeds threshold
        combined_threshold = stem_max_mm
        if score_A < score_B and score_A <= combined_threshold:
            variant = 1
            mm_stem = mm_stem1
            mm_ra = mm_ra1
        elif score_B < score_A and score_B <= combined_threshold:
            variant = 2
            mm_stem = mm_stem2
            mm_ra = mm_ra2 if right_anchor_2 else mm_ra1
        elif score_A == score_B and score_A <= combined_threshold:
            variant = 'AMBIGUOUS'
            mm_stem = mm_stem1
            mm_ra = mm_ra1
        else:
            continue  # Try next gRNA length

        grna_seq = seq[grna_start:grna_end]
        ibar_seq = seq[ibar_start:ibar_end]

        # Classify anchor combination for stats
        if variant == 1:
            anchor_combo = 'A'  # stem1 + right_anchor_1
        elif variant == 2:
            anchor_combo = 'B'  # stem2 + right_anchor_2
        else:
            anchor_combo = 'ambiguous'

        # Flag cross-pair if applicable (only for R1 where we have two right anchors)
        if has_right_anchor and right_anchor_2:
            score_cross_12 = mm_stem1 + mm_ra2  # stem1 + right2
            score_cross_21 = mm_stem2 + mm_ra1  # stem2 + right1
            chosen_score = score_A if variant == 1 else score_B
            if score_cross_12 < chosen_score or score_cross_21 < chosen_score:
                anchor_combo = 'cross'

        return {
            'grna_seq': grna_seq,
            'ibar_seq': ibar_seq,
            'variant': variant,
            'grna_length': grna_len,
            'anchor_offset': anchor_pos,
            'anchor_dist': anchor_dist,
            'stem_mm': mm_stem,
            'right_anchor_mm': mm_ra,
            'anchor_combo': anchor_combo,
        }

    return None


def extract_from_r1(r1_seq, cfg):
    """Extract pos1 guide unit from R1.

    R1 structure: [R1_primer(18)]─[gRNA₁(20/21)]─[stem(19)]─[iBAR₁(7)]─[tracr(~68)]

    Uses R1_primer as left anchor (like TSO in Step 3), then forward-computes
    positions with combined stem+tracr scoring to select gRNA length.

    Args:
        r1_seq: R1 read sequence (forward strand)
        cfg: guideseq config dict

    Returns:
        dict with Step-3-compatible metadata, or None if extraction fails.
        Keys: grna_seq, ibar_seq, variant, grna_length, anchor_offset,
              anchor_dist, stem_mm, right_anchor_mm, anchor_combo
    """
    max_anchor_length = cfg.get('max_anchor_length', 20)
    max_mm = cfg.get('anchor_max_mm', 5)
    stem_max_mm = cfg.get('stem_max_mm', 5)

    # Right anchors: tracr_1 / tracr_2, trimmed to keep 5' end (near iBAR)
    tracr_1 = _trim_anchor(cfg.get('tracr_1', ''), max_anchor_length, 'right')
    tracr_2 = _trim_anchor(cfg.get('tracr_2', ''), max_anchor_length, 'right')

    return _extract_guide_unit(
        seq=r1_seq,
        left_anchor=cfg.get('R1_primer', ''),
        stem_1=cfg['stem_1'], stem_2=cfg['stem_2'],
        right_anchor_1=tracr_1, right_anchor_2=tracr_2,
        grna_lengths=cfg['grna_lengths'],
        ibar_length=cfg['ibar_length'],
        left_anchor_max_mm=max_mm,
        stem_max_mm=stem_max_mm,
        right_anchor_max_mm=max_mm,
    )


def extract_from_r2(r2_seq, cfg):
    """Extract pos2 guide unit from R2.

    R2rc structure: [tRNA_leader(6)]─[tRNA(~72)]─[gRNA₂(20/21)]─[stem(19)]─[iBAR₂(7)]─[R2_primer_rc]

    Uses tRNA as left anchor, combined stem+R2_primer_rc scoring to select
    gRNA length.

    Args:
        r2_seq: R2 read sequence (original R2 orientation, NOT revcomp'd)
        cfg: guideseq config dict

    Returns:
        dict with Step-3-compatible metadata + 'tRNA' and 'tRNA_meta' keys,
        or None if extraction fails.
    """
    r2rc = revcomp(r2_seq)

    max_anchor_length = cfg.get('max_anchor_length', 20)
    max_mm = cfg.get('anchor_max_mm', 5)
    stem_max_mm = cfg.get('stem_max_mm', 5)

    # Right anchor: R2_primer is in original R2 orientation → revcomp for R2rc
    r2_primer_rc = revcomp(cfg.get('R2_primer', ''))

    # Left anchor: try each tRNA variant, pick best match
    trna_seqs = cfg.get('tRNA_seqs', {})
    tRNA_leader = cfg.get('tRNA_leader', '')

    best_result = None
    best_tRNA_name = NOT_FOUND
    best_tRNA_meta = NOT_FOUND

    for variant_name, trna_full_seq in trna_seqs.items():
        if not trna_full_seq:
            continue

        # Build left anchor: tRNA_leader + tRNA body (trimmed to keep 3' end near gRNA)
        tRNA_trimmed = _trim_anchor(trna_full_seq, max_anchor_length, 'left')
        left_anchor = tRNA_leader + tRNA_trimmed

        result = _extract_guide_unit(
            seq=r2rc,
            left_anchor=left_anchor,
            stem_1=cfg['stem_1'], stem_2=cfg['stem_2'],
            right_anchor_1=r2_primer_rc, right_anchor_2=None,
            grna_lengths=cfg['grna_lengths'],
            ibar_length=cfg['ibar_length'],
            left_anchor_max_mm=max_mm,
            stem_max_mm=stem_max_mm,
            right_anchor_max_mm=max_mm,
        )
        if result is not None:
            combined_score = result['stem_mm'] + result['right_anchor_mm']
            if best_result is None or combined_score < (
                    best_result['stem_mm'] + best_result['right_anchor_mm']):
                best_result = result
                best_tRNA_name = f"tRNA_{variant_name}"
                best_tRNA_meta = f"{variant_name}={result['anchor_dist']}"

    if best_result is not None:
        best_result['tRNA'] = best_tRNA_name
        best_result['tRNA_meta'] = best_tRNA_meta

    return best_result


def extract_from_construct(merged_seq, cfg):
    """Extract both guide units from a merged construct.

    Splits via tRNA_leader + tRNA as linker:
      - left half  → extract_from_r1 (pos1 guide unit)
      - right half → extract_from_r2 via revcomp (pos2 guide unit)

    Args:
        merged_seq: full merged amplicon (forward strand)
        cfg: guideseq config dict

    Returns:
        dict with pos1 + pos2 results and tRNA identity, or dict with 'fail' key
    """
    leader_max_mm = cfg.get('tRNA_leader_max_mm', 1)
    trna_leader = cfg.get('tRNA_leader', '')
    trna_seqs = cfg.get('tRNA_seqs', {})

    if not trna_leader:
        return {'fail': 'no_tRNA_leader_in_config'}

    # --- Find tRNA_leader + tRNA_X as linker ---
    linker_found = False
    linker_start = -1
    linker_end = -1
    tRNA_name = NOT_FOUND
    tRNA_meta = NOT_FOUND

    for variant_name, trna_seq in trna_seqs.items():
        if not trna_seq:
            continue
        full_linker = trna_leader + trna_seq
        linker_max_mm = max(2, int(len(full_linker) * 0.05))
        pos, end, mm = find_anchor(merged_seq, full_linker, linker_max_mm)
        if pos >= 0:
            linker_start = pos
            linker_end = end
            tRNA_name = f"tRNA_{variant_name}"
            tRNA_meta = f"{variant_name}={mm}"
            linker_found = True
            break

    # --- Fallback: tRNA_leader alone with bounded position ---
    if not linker_found:
        search_start = max(0, 100)
        search_end = min(len(merged_seq), 180)
        search_region = merged_seq[search_start:search_end]

        pos, end, mm = find_anchor(search_region, trna_leader, leader_max_mm)
        if pos >= 0:
            linker_start = search_start + pos
            linker_end = search_start + end
            linker_found = True

    if not linker_found:
        return {'fail': 'linker_not_found'}

    # --- Split: left half = pos1 (like R1), right half = pos2 ---
    left_half = merged_seq[:linker_start]
    right_half = merged_seq[linker_end:]

    # Pos1: extract like R1 (left half has R1_primer → gRNA₁ → stem → iBAR → tracr)
    pos1_result = extract_from_r1(left_half, cfg)

    # Pos2: right half structure is [gRNA₂]─[stem₂]─[iBAR₂]─[tracr₂]
    # This is the same as R1 after the primer, so use empty left_anchor
    # (gRNA starts at position 0) with tracr as right anchors.
    max_anchor_length = cfg.get('max_anchor_length', 20)
    max_mm = cfg.get('anchor_max_mm', 5)
    stem_max_mm_val = cfg.get('stem_max_mm', 5)
    tracr_1 = _trim_anchor(cfg.get('tracr_1', ''), max_anchor_length, 'right')
    tracr_2 = _trim_anchor(cfg.get('tracr_2', ''), max_anchor_length, 'right')

    pos2_result = _extract_guide_unit(
        seq=right_half,
        left_anchor='',  # gRNA₂ starts at position 0 (right after tRNA linker)
        stem_1=cfg['stem_1'], stem_2=cfg['stem_2'],
        right_anchor_1=tracr_1, right_anchor_2=tracr_2,
        grna_lengths=cfg['grna_lengths'],
        ibar_length=cfg['ibar_length'],
        left_anchor_max_mm=max_mm,
        stem_max_mm=stem_max_mm_val,
        right_anchor_max_mm=max_mm,
    )

    # If tRNA was resolved at the linker step, override
    if pos2_result is not None:
        if tRNA_name != NOT_FOUND:
            pos2_result['tRNA'] = tRNA_name
        if tRNA_meta != NOT_FOUND:
            pos2_result['tRNA_meta'] = tRNA_meta

    return {
        'pos1': pos1_result,
        'pos2': pos2_result,
        'tRNA': tRNA_name,
        'tRNA_meta': tRNA_meta,
    }



# ============================================================================
# Guide Matching (follows Step 3 pipeline)
# ============================================================================

def build_matchers(ref_df, max_lev_dist=3, max_lev_dist_ibar=3):
    """Build matchers and lookups from reference, matching step 3 logic.

    Args:
        ref_df: reference DataFrame
        max_lev_dist: max Levenshtein distance for gRNA
        max_lev_dist_ibar: max Levenshtein distance for iBAR

    Returns:
        dict with matchers and lookups
    """
    # gRNA matcher
    grna_matcher = NeighborhoodMatcher(
        ref_df['grna_seq'].tolist(),
        ref_df['grna_name'].tolist(),
        max_dist=max_lev_dist
    )

    # iBAR matcher (global)
    ibar_matcher = NeighborhoodMatcher(
        ref_df['ibar_seq'].tolist(),
        ref_df['ibar_seq'].tolist(),  # iBAR matched by sequence
        max_dist=max_lev_dist_ibar
    )

    # Exact combined lookup: (grna_seq, ibar_seq) → (grna_name, target_gene, pair_id, target_combo)
    combined_lookup = {}
    for _, row in ref_df.iterrows():
        key = (row['grna_seq'], row['ibar_seq'])
        combined_lookup[key] = (
            row['grna_name'], row['target_gene'],
            row['pair_id'], row['target_combo']
        )

    # grna_name → target_gene
    ref_grna_to_target = dict(zip(ref_df['grna_name'], ref_df['target_gene']))

    # grna_name → list of (ibar_seq,) for within-guide matching
    grna_to_ibars = {}
    for _, row in ref_df.iterrows():
        grna_to_ibars.setdefault(row['grna_name'], []).append(row['ibar_seq'])

    # grna_name → position (1 or 2)
    grna_to_position = dict(zip(ref_df['grna_name'], ref_df['position'].astype(int)))

    # (grna_name, ibar_seq) → pair_id
    grna_ibar_to_pairid = {}
    for _, row in ref_df.iterrows():
        grna_ibar_to_pairid[(row['grna_name'], row['ibar_seq'])] = row['pair_id']

    # grna_name → target_gene (for within/across classification)
    grna_to_target = dict(zip(ref_df['grna_name'], ref_df['target_gene']))

    # Detect shared gRNA sequences (appear at both pos1 and pos2)
    _seq_pos = ref_df.groupby('grna_seq')['position'].apply(set)
    shared_grna_seqs = set(_seq_pos[_seq_pos.apply(len) > 1].index)

    # grna_seq → {position: grna_name} for shared sequences
    grna_seq_to_pos_name = {}
    for seq in shared_grna_seqs:
        sub = ref_df[ref_df['grna_seq'] == seq].drop_duplicates('grna_name')
        grna_seq_to_pos_name[seq] = dict(
            zip(sub['position'].astype(int), sub['grna_name'])
        )

    if shared_grna_seqs:
        print(f"    NOTE: {len(shared_grna_seqs)} shared gRNA sequence(s) detected "
              f"(appear at both pos1 and pos2)")

    # pair_id → {position: grna_name} for swap classification (Stage 4)
    pairid_to_guides = {}
    for _, row in ref_df.iterrows():
        pid = row['pair_id']
        pos = int(row['position'])
        pairid_to_guides.setdefault(pid, {})[pos] = row['grna_name']

    return {
        'grna_matcher': grna_matcher,
        'ibar_matcher': ibar_matcher,
        'combined_lookup': combined_lookup,
        'ref_grna_to_target': ref_grna_to_target,
        'grna_to_ibars': grna_to_ibars,
        'grna_to_position': grna_to_position,
        'grna_ibar_to_pairid': grna_ibar_to_pairid,
        'grna_to_target': grna_to_target,
        'max_lev_ibar': max_lev_dist_ibar,
        'shared_grna_seqs': shared_grna_seqs,
        'grna_seq_to_pos_name': grna_seq_to_pos_name,
        'pairid_to_guides': pairid_to_guides,
    }


# ============================================================================
# Stage 2 Helper — gRNA → grna_name resolution
# ============================================================================

def _resolve_grna_name_for_pos(grna_seq, stem_name, grna_matcher,
                                shared_grna_seqs, grna_seq_to_pos_name):
    """Resolve grna_name for one position using gRNA seq + stem identity.

    Primary resolution: stem identity narrows shared gRNAs to a unique name.
    Fallback: gRNA-only match (if stem is ambiguous and gRNA is unique).

    Args:
        grna_seq: extracted gRNA sequence (or NOT_FOUND)
        stem_name: "stem_1", "stem_2", AMBIGUOUS, or NOT_FOUND
        grna_matcher: NeighborhoodMatcher for gRNAs
        shared_grna_seqs: set of shared gRNA sequences
        grna_seq_to_pos_name: {grna_seq: {position: grna_name}}

    Returns:
        (grna_name, candidates)
        grna_name: resolved str, AMBIGUOUS, or None (no match)
        candidates: list of candidate names when ambiguous, else None
    """
    if grna_seq == NOT_FOUND:
        return None, None

    # Determine expected position from stem
    if stem_name == 'stem_1':
        expected_position = 1
    elif stem_name == 'stem_2':
        expected_position = 2
    else:
        expected_position = None  # ambiguous or other

    # Fuzzy match gRNA
    grna_name, grna_dist, grna_unique = grna_matcher.match(grna_seq)

    if grna_name is None:
        return None, None

    if grna_unique:
        # Single match → resolved
        return grna_name, None

    # Multiple matches (shared gRNA sequence)
    candidates = grna_name if isinstance(grna_name, list) else [grna_name]

    # Try to disambiguate with stem
    if expected_position is not None and grna_seq in grna_seq_to_pos_name:
        pos_names = grna_seq_to_pos_name[grna_seq]
        if expected_position in pos_names:
            return pos_names[expected_position], candidates

    return AMBIGUOUS, candidates


# ============================================================================
# Stage 3 Helper — iBAR correction + ambiguous-guide resolution
# ============================================================================

def _correct_ibar_and_classify(grna_name, grna_candidates, ibar_seq,
                                grna_to_ibars, grna_ibar_to_pairid,
                                max_lev_ibar):
    """Correct iBAR and resolve pair_id.

    Path A (resolved grna_name): within-guide iBAR correction.
    Path B (ambiguous grna_name): joint gRNA-iBAR resolution across candidates.

    Returns:
        dict with: grna_name_final, ibar_corrected, pair_id, status
        status: "correct" | "grna_only" | "no_match" | "ambiguous"
    """
    result = {
        'grna_name_final': grna_name,
        'ibar_corrected': None,
        'pair_id': None,
        'status': 'no_match',
    }

    if grna_name is None:
        return result

    if ibar_seq == NOT_FOUND:
        result['status'] = 'grna_only'
        return result

    # Path B: ambiguous grna_name → joint gRNA-iBAR resolution
    if grna_name == AMBIGUOUS:
        if not grna_candidates:
            return result

        best_name = None
        best_ibar = None
        best_dist = max_lev_ibar + 1
        best_pairid = None
        n_best = 0

        for cand_name in grna_candidates:
            guide_ibars = grna_to_ibars.get(cand_name, [])
            matched_ibar, dist = _match_ibar_within_guide(
                ibar_seq, guide_ibars, max_lev_ibar)
            if matched_ibar is not None:
                if dist < best_dist:
                    best_dist = dist
                    best_name = cand_name
                    best_ibar = matched_ibar
                    best_pairid = grna_ibar_to_pairid.get(
                        (cand_name, matched_ibar))
                    n_best = 1
                elif dist == best_dist:
                    n_best += 1

        if best_name is not None and n_best == 1:
            result['grna_name_final'] = best_name
            result['ibar_corrected'] = best_ibar
            result['pair_id'] = best_pairid
            result['status'] = 'correct' if best_pairid else 'grna_only'
        else:
            result['status'] = 'ambiguous'
        return result

    # Path A: resolved grna_name → within-guide iBAR correction
    guide_ibars = grna_to_ibars.get(grna_name, [])
    matched_ibar, ibar_dist = _match_ibar_within_guide(
        ibar_seq, guide_ibars, max_lev_ibar)

    if matched_ibar is not None:
        result['ibar_corrected'] = matched_ibar
        pair_id = grna_ibar_to_pairid.get((grna_name, matched_ibar))
        result['pair_id'] = pair_id
        result['status'] = 'correct' if pair_id else 'grna_only'
    else:
        result['status'] = 'grna_only'

    return result


# ============================================================================
# Stage 4 Helper — 3-level swap taxonomy
# ============================================================================

def _classify_position_swap(actual_name, expected_name, grna_to_target,
                             counts, pos_label):
    """Classify a swap at one position into the 3-level taxonomy.

    Called only for non-intact reads where the actual guide and
    the expected guide (from the other position's pair_id) are known.

    Levels:
      within_guide  — same grna_name, different pair_id (iBAR swap)
      between_guide — different grna_name, same target_gene
      across_target — different target_gene
    """
    if actual_name == expected_name:
        counts[f'swap_within_guide_{pos_label}'] += 1
    else:
        actual_target = grna_to_target.get(actual_name, '')
        expected_target = grna_to_target.get(expected_name, '')
        if actual_target and expected_target and actual_target == expected_target:
            counts[f'swap_between_guide_{pos_label}'] += 1
        else:
            counts[f'swap_across_target_{pos_label}'] += 1


# ============================================================================
# FASTQ Auto-Detection
# ============================================================================

def detect_fastq_pairs(fastq_dir):
    """Auto-detect paired R1/R2 FASTQ files from a directory.

    Scans for files matching patterns like:
      *_R1*.fastq.gz / *_R2*.fastq.gz
      *_R1*.fq.gz    / *_R2*.fq.gz
      *_R1*.fastq    / *_R2*.fastq

    Pairs are matched by replacing '_R1' with '_R2' in the filename.
    Library name is derived from the common prefix before _R1.

    Args:
        fastq_dir: path to directory containing FASTQ files

    Returns:
        list of (name, r1_path, r2_path) tuples, sorted by name
    """
    import re

    fastq_dir = os.path.abspath(fastq_dir)
    if not os.path.isdir(fastq_dir):
        raise FileNotFoundError(f"Directory not found: {fastq_dir}")

    # Find all R1 files
    r1_pattern = re.compile(r'(.+)[_.]R1[_.](.+\.(fastq|fq)(\.gz)?)$', re.IGNORECASE)
    r1_files = []
    for fname in sorted(os.listdir(fastq_dir)):
        m = r1_pattern.match(fname)
        if m:
            r1_files.append(fname)

    pairs = []
    for r1_fname in r1_files:
        # Derive R2 filename by replacing R1 → R2
        r2_fname = re.sub(r'([_.])(R1)([_.])', r'\1R2\3', r1_fname, count=1)
        r2_path = os.path.join(fastq_dir, r2_fname)

        if not os.path.exists(r2_path):
            print(f"  Warning: R2 not found for {r1_fname}, skipping")
            continue

        # Derive library name from prefix before _R1
        name = re.split(r'[_.]R1[_.]', r1_fname, maxsplit=1)[0]
        r1_path = os.path.join(fastq_dir, r1_fname)

        pairs.append((name, r1_path, r2_path))

    if not pairs:
        print(f"  Warning: No paired FASTQ files found in {fastq_dir}")
    else:
        print(f"  Detected {len(pairs)} paired FASTQ libraries in {fastq_dir}:")
        for name, r1, r2 in pairs:
            print(f"    {name}: {os.path.basename(r1)} + {os.path.basename(r2)}")

    return pairs


# ============================================================================
# FASTQ Processing
# ============================================================================

def _read_fastq_pairs(r1_path, r2_path):
    """Yield paired (read_id, r1_seq, r2_seq) from FASTQ files.

    Assumes reads are synchronized (same order in R1 and R2).
    The read_id is extracted from the first word of the FASTQ header line
    (after stripping the '@' prefix).
    """
    open_fn = gzip.open if r1_path.endswith('.gz') else open
    mode = 'rt' if r1_path.endswith('.gz') else 'r'

    open_fn2 = gzip.open if r2_path.endswith('.gz') else open
    mode2 = 'rt' if r2_path.endswith('.gz') else 'r'

    with open_fn(r1_path, mode) as f1, open_fn2(r2_path, mode2) as f2:
        while True:
            # R1
            header1 = f1.readline().strip()
            seq1 = f1.readline().strip()
            plus1 = f1.readline().strip()
            qual1 = f1.readline().strip()

            # R2
            header2 = f2.readline().strip()
            seq2 = f2.readline().strip()
            plus2 = f2.readline().strip()
            qual2 = f2.readline().strip()

            if not header1 or not header2:
                break

            # Extract read ID (first word after '@')
            read_id = header1.split()[0].lstrip('@')
            yield read_id, seq1, seq2


# TSV output columns — one row per read pair, two guide units (pos1 + pos2)
_GU_FIELDS = [
    'grna_seq', 'ibar_seq', 'variant', 'grna_length',
    'anchor_offset', 'anchor_dist', 'stem_mm', 'right_anchor_mm', 'anchor_combo',
]

TSV_COLUMNS = ['read_id']
for _pref in ('pos1_', 'pos2_'):
    for _f in _GU_FIELDS:
        TSV_COLUMNS.append(_pref + _f)
TSV_COLUMNS.extend(['tRNA', 'tRNA_meta', 'source'])


def _gu_to_row(result, prefix):
    """Flatten a guide-unit result dict into prefixed row entries."""
    row = {}
    if result is None:
        for f in _GU_FIELDS:
            row[prefix + f] = NOT_FOUND
    else:
        for f in _GU_FIELDS:
            row[prefix + f] = result.get(f, NOT_FOUND)
    return row


def _extract_one_pair(read_id, r1_seq, r2_seq, guideseq_cfg, stats):
    """Extract guide units from one read pair using dual-path architecture.

    Path 1 (independent): extract_from_r1 + extract_from_r2 on raw reads.
    Path 2 (construct): merge → extract_from_construct on merged reads.

    Each read pair yields two guide units (pos1 from R1, pos2 from R2).

    Args:
        read_id: FASTQ read pair identifier
        r1_seq, r2_seq: raw read sequences
        guideseq_cfg: config dict
        stats: Counter updated in-place

    Returns:
        (indep_row, construct_row, fasta_entry)
        - indep_row: dict with TSV_COLUMNS or None (if filtered)
        - construct_row: dict with TSV_COLUMNS or None (if merge/split failed)
        - fasta_entry: (header, seq) tuple or None
    """
    # --- Pre-filter ---
    min_r1 = guideseq_cfg.get('min_r1_length', 65)
    min_r2 = guideseq_cfg.get('min_r2_length', 69)
    max_ns = guideseq_cfg.get('max_consecutive_ns', 10)
    ns_pattern = 'N' * max_ns

    if len(r1_seq) < min_r1 or len(r2_seq) < min_r2:
        stats['filtered_short'] += 1
        return None, None, None
    if ns_pattern in r1_seq or ns_pattern in r2_seq:
        stats['filtered_ns'] += 1
        return None, None, None

    # --- Path 1: Independent extraction (all reads) ---
    # R1 → pos1 guide unit (R1_primer left anchor handled inside extract_from_r1)
    pos1_result = extract_from_r1(r1_seq, guideseq_cfg)
    # R2 → pos2 guide unit (tRNA left anchor handled inside extract_from_r2)
    pos2_result = extract_from_r2(r2_seq, guideseq_cfg)

    indep_row = {'read_id': read_id, 'source': 'independent'}
    indep_row.update(_gu_to_row(pos1_result, 'pos1_'))
    indep_row.update(_gu_to_row(pos2_result, 'pos2_'))
    # tRNA from pos2 result
    indep_row['tRNA'] = pos2_result.get('tRNA', NOT_FOUND) if pos2_result else NOT_FOUND
    indep_row['tRNA_meta'] = pos2_result.get('tRNA_meta', NOT_FOUND) if pos2_result else NOT_FOUND

    stats['n_independent'] += 1
    if pos1_result is not None:
        stats['pos1_extracted'] += 1
        stats[f'pos1_anchor_combo_{pos1_result["anchor_combo"]}'] += 1
        stats[f'pos1_grna_length_{pos1_result["grna_length"]}'] += 1
    else:
        stats['pos1_failed'] += 1
    if pos2_result is not None:
        stats['pos2_extracted'] += 1
        stats[f'pos2_anchor_combo_{pos2_result["anchor_combo"]}'] += 1
        stats[f'pos2_grna_length_{pos2_result["grna_length"]}'] += 1
    else:
        stats['pos2_failed'] += 1

    # --- Path 2: Construct extraction (merged reads) ---
    construct_row = None
    fasta_entry = None

    # Merge uses untrimmed reads (full length needed for overlap detection)
    merged_seq, overlap_len = merge_reads(r1_seq, r2_seq, guideseq_cfg)

    if merged_seq is not None:
        stats['merge_success'] += 1
        fasta_entry = (f">{read_id} overlap={overlap_len}", merged_seq)

        construct_result = extract_from_construct(merged_seq, guideseq_cfg)

        if 'fail' not in construct_result:
            c_pos1 = construct_result.get('pos1')
            c_pos2 = construct_result.get('pos2')

            construct_row = {'read_id': read_id, 'source': 'construct'}
            construct_row.update(_gu_to_row(c_pos1, 'pos1_'))
            construct_row.update(_gu_to_row(c_pos2, 'pos2_'))
            construct_row['tRNA'] = construct_result.get('tRNA', NOT_FOUND)
            construct_row['tRNA_meta'] = construct_result.get('tRNA_meta', NOT_FOUND)
            stats['n_construct'] += 1
        else:
            stats[f'construct_fail_{construct_result["fail"]}'] += 1
    else:
        stats['merge_failed'] += 1

    return indep_row, construct_row, fasta_entry


# ============================================================================
# Library Processing — Extraction Stage (writes TSV + FASTA)
# ============================================================================

# Module-level state for multiprocessing workers
_extraction_worker_state = {}


def _init_extraction_worker(guideseq_cfg):
    """Initialize per-worker state for multiprocessing Pool."""
    _extraction_worker_state['cfg'] = guideseq_cfg


def _extract_batch(batch):
    """Process a batch of read pairs in a worker process.

    Args:
        batch: list of (read_id, r1_seq, r2_seq) tuples

    Returns:
        (indep_rows, construct_rows, fasta_entries, stats)
    """
    cfg = _extraction_worker_state['cfg']
    stats = Counter()
    indep_rows = []
    construct_rows = []
    fasta_entries = []

    for read_id, r1_seq, r2_seq in batch:
        indep_row, construct_row, fasta_entry = _extract_one_pair(
            read_id, r1_seq, r2_seq, cfg, stats)
        if indep_row is not None:
            indep_rows.append(indep_row)
        if construct_row is not None:
            construct_rows.append(construct_row)
        if fasta_entry is not None:
            fasta_entries.append(fasta_entry)

    return indep_rows, construct_rows, fasta_entries, stats

def process_guideseq_library(r1_path, r2_path, guideseq_cfg, out_dir,
                             max_reads=0, n_workers=1, verbose=True):
    """Extract elements from one guide-seq paired-end library.

    Runs dual-path extraction on every read pair and writes:
      - extraction_independent.tsv.gz (all reads)
      - extraction_construct.tsv.gz (successfully merged reads only)
      - constructs.fasta.gz (merged amplicon sequences)

    Supports multiprocessing via n_workers > 1.

    Args:
        r1_path, r2_path: paths to paired FASTQ files
        guideseq_cfg: dict with extraction config
        out_dir: output directory for this library
        max_reads: max reads to process (0 = all)
        n_workers: number of parallel workers (1 = serial)
        verbose: print progress

    Returns:
        dict with stats (Counter), indep_tsv_path, construct_tsv_path, fasta_path
    """
    os.makedirs(out_dir, exist_ok=True)

    indep_tsv_path = os.path.join(out_dir, 'extraction_independent.tsv.gz')
    construct_tsv_path = os.path.join(out_dir, 'extraction_construct.tsv.gz')
    fasta_path = os.path.join(out_dir, 'constructs.fasta.gz')

    stats = Counter()

    # Estimate total reads from file size
    try:
        file_size = os.path.getsize(r1_path)
        est_reads = file_size // 80 if r1_path.endswith('.gz') else file_size // 400
    except OSError:
        est_reads = None

    total = min(max_reads, est_reads) if (max_reads > 0 and est_reads) else est_reads

    header_line = '\t'.join(TSV_COLUMNS) + '\n'

    if n_workers > 1:
        # --- Parallel mode ---
        from multiprocessing import Pool

        # Collect read batches
        batch_size = 10_000
        batches = []
        current_batch = []
        n_total = 0

        for read_id, r1_seq, r2_seq in _read_fastq_pairs(r1_path, r2_path):
            n_total += 1
            if max_reads > 0 and n_total > max_reads:
                break
            current_batch.append((read_id, r1_seq, r2_seq))
            if len(current_batch) >= batch_size:
                batches.append(current_batch)
                current_batch = []
        if current_batch:
            batches.append(current_batch)

        if verbose:
            print(f"    {n_total:,} read pairs -> {len(batches)} batches (batch_size={batch_size:,})")

        with Pool(n_workers, initializer=_init_extraction_worker,
                  initargs=(guideseq_cfg,)) as pool:

            with (gzip.open(indep_tsv_path, 'wt') as f_indep,
                  gzip.open(construct_tsv_path, 'wt') as f_construct,
                  gzip.open(fasta_path, 'wt') as f_fasta):

                f_indep.write(header_line)
                f_construct.write(header_line)

                pbar = tqdm(
                    pool.imap_unordered(_extract_batch, batches),
                    total=len(batches),
                    desc="    Extracting",
                    unit=" batches",
                    disable=not verbose,
                )

                for batch_indep, batch_construct, batch_fasta, batch_stats in pbar:
                    stats += batch_stats
                    for row in batch_indep:
                        f_indep.write('\t'.join(str(row[c]) for c in TSV_COLUMNS) + '\n')
                    for row in batch_construct:
                        f_construct.write('\t'.join(str(row[c]) for c in TSV_COLUMNS) + '\n')
                    for header, seq in batch_fasta:
                        f_fasta.write(f"{header}\n{seq}\n")

        stats['n_total'] = n_total

    else:
        # --- Serial mode (n_workers == 1) ---
        with (gzip.open(indep_tsv_path, 'wt') as f_indep,
              gzip.open(construct_tsv_path, 'wt') as f_construct,
              gzip.open(fasta_path, 'wt') as f_fasta):

            f_indep.write(header_line)
            f_construct.write(header_line)

            pbar = tqdm(
                _read_fastq_pairs(r1_path, r2_path),
                total=total,
                desc="    Extracting",
                unit=" pairs",
                disable=not verbose,
            )

            n_total = 0
            for read_id, r1_seq, r2_seq in pbar:
                n_total += 1
                if max_reads > 0 and n_total > max_reads:
                    break

                indep_row, construct_row, fasta_entry = _extract_one_pair(
                    read_id, r1_seq, r2_seq, guideseq_cfg, stats)

                if indep_row is not None:
                    f_indep.write('\t'.join(str(indep_row[c]) for c in TSV_COLUMNS) + '\n')

                if construct_row is not None:
                    f_construct.write('\t'.join(str(construct_row[c]) for c in TSV_COLUMNS) + '\n')

                if fasta_entry is not None:
                    header, seq = fasta_entry
                    f_fasta.write(f"{header}\n{seq}\n")

        stats['n_total'] = n_total

    if verbose:
        print(f"    Total read pairs:  {stats.get('n_total', 0):,}")
        print(f"    Filtered (short):  {stats.get('filtered_short', 0):,}")
        print(f"    Filtered (Ns):     {stats.get('filtered_ns', 0):,}")
        print(f"    Pos1 extracted:    {stats.get('pos1_extracted', 0):,}")
        print(f"    Pos1 failed:       {stats.get('pos1_failed', 0):,}")
        print(f"    Pos2 extracted:    {stats.get('pos2_extracted', 0):,}")
        print(f"    Pos2 failed:       {stats.get('pos2_failed', 0):,}")
        print(f"    Merge success:     {stats.get('merge_success', 0):,}")
        print(f"    Construct rows:    {stats.get('n_construct', 0):,}")
        print(f"    Written: {os.path.basename(indep_tsv_path)}")
        print(f"    Written: {os.path.basename(construct_tsv_path)}")
        print(f"    Written: {os.path.basename(fasta_path)}")

    return {
        'stats': dict(stats),
        'indep_tsv_path': indep_tsv_path,
        'construct_tsv_path': construct_tsv_path,
        'fasta_path': fasta_path,
    }


# ============================================================================
# 4-Stage Matching Pipeline — from saved TSVs
# ============================================================================

def match_and_compute_from_tsv(tsv_path, ref_df, matchers, verbose=True):
    """Read extraction TSV and run the 4-stage matching + metrics pipeline.
    
    1. Match by gRNA
    2. Disambiguate by Stem
    3. Disambiguate by iBAR
    4. Verify iBAR & resolve pair ID
    """
    grna_matcher = matchers['grna_matcher']
    grna_seq_to_pos_name = matchers.get('grna_seq_to_pos_name', {})
    grna_to_ibars = matchers['grna_to_ibars']
    grna_ibar_to_pairid = matchers['grna_ibar_to_pairid']
    grna_to_target = matchers['grna_to_target']
    grna_to_position = matchers['grna_to_position']
    pairid_to_guides = matchers['pairid_to_guides']
    max_lev_ibar = matchers['max_lev_ibar']
    combined_lookup = matchers.get('combined_lookup', {})
    shared_grna_seqs = matchers.get('shared_grna_seqs', set())

    import pandas as pd
    from collections import Counter
    from qc_q1_guideseq import NOT_FOUND, AMBIGUOUS, _classify_position_swap, hamming_distance
    
    df = pd.read_csv(tsv_path, sep='\t', dtype=str, na_filter=False)
    N_total = len(df)

    if N_total == 0:
        empty_rates = {k: 0.0 for k in [
            'pos1_fail', 'pos2_fail',
            'grna_del1', 'grna_del2', 'decoupled_pos1', 'decoupled_pos2',
            'ambiguous_pos1', 'ambiguous_pos2',
            'swap_within_guide_pos1', 'swap_within_guide_pos2',
            'swap_between_guide_pos1', 'swap_between_guide_pos2',
            'swap_across_target_pos1', 'swap_across_target_pos2',
            'pos1_pos2_swap', 'intact_rate',
        ]}
        empty_rates['N_total'] = 0
        empty_rates['N_valid'] = 0
        return {'rates': empty_rates, 'stem_report': {}, 'counts': {}}

    # pos1/pos2 extraction failure = variant is NOT_FOUND
    pos1_fail_mask = df['pos1_variant'].astype(str) == NOT_FOUND
    pos2_fail_mask = df['pos2_variant'].astype(str) == NOT_FOUND
    N_pos1_fail = int(pos1_fail_mask.sum())
    N_pos2_fail = int(pos2_fail_mask.sum())
    N_both_fail = int((pos1_fail_mask & pos2_fail_mask).sum())

    valid_mask = (~pos1_fail_mask) & (~pos2_fail_mask)
    df_valid = df[valid_mask].reset_index(drop=True)
    N_valid = len(df_valid)
    
    stem_report = {
        'N_total': N_total, 'N_pos1_fail': N_pos1_fail, 'N_pos2_fail': N_pos2_fail,
        'N_both_fail': N_both_fail, 'N_valid': N_valid,
    }
    
    if verbose:
        print(f"    Stage 1: {N_total:,} reads -> {N_valid:,} with both positions extracted")

    if N_valid == 0:
        rates = {'pos1_fail': N_pos1_fail/N_total, 'pos2_fail': N_pos2_fail/N_total, 'N_total': N_total, 'N_valid': 0}
        return {'rates': rates, 'stem_report': stem_report, 'counts': {}}

    pos1_gseqs = df_valid['pos1_grna_seq'].values
    pos1_variants = df_valid['pos1_variant'].astype(str).values
    pos1_ibars = df_valid['pos1_ibar_seq'].values
    pos2_gseqs = df_valid['pos2_grna_seq'].values
    pos2_variants = df_valid['pos2_variant'].astype(str).values
    pos2_ibars = df_valid['pos2_ibar_seq'].values

    # Pre-cache gRNA mapping
    grna_cache = {}
    unique_gseqs = set(pos1_gseqs) | set(pos2_gseqs)
    unique_gseqs.discard(NOT_FOUND)
    for seq in unique_gseqs:
        grna_cache[seq] = grna_matcher.match(seq)

    def _resolve_guide_unit(grna_seq, variant_str, ibar_seq, expected_pos):
        """Resolve a guide unit to a pair_id.

        Follows Step 3 single-cell pipeline logic:
          Step 0: Exact combined (grna_seq, ibar_seq) lookup
          Step 1: Fuzzy gRNA match via NeighborhoodMatcher
          Step 2: For shared gRNAs, iBAR-based disambiguation
                  (closest unique-best per candidate, via _match_ibar_within_guide)
          Step 3: Fall back to stem variant (skipped if AMBIGUOUS)
          Step 4: iBAR verification & pair_id resolution

        variant_str: '1', '2', or 'AMBIGUOUS' (from extraction variant field)
        """
        if grna_seq == NOT_FOUND:
            return None, 'grna_deletion'

        # --- Step 0: Exact combined lookup ---
        if ibar_seq != NOT_FOUND:
            combined_key = (grna_seq, ibar_seq)
            if combined_key in combined_lookup:
                grna_name, target_gene, pair_id, target_combo = combined_lookup[combined_key]
                counts['exact_combined'] += 1
                return pair_id, 'ok'

        # --- Step 1: Fuzzy gRNA match ---
        grna_name, grna_dist, grna_unique = grna_cache[grna_seq]
        if grna_name is None:
            return None, 'grna_deletion'

        resolved_name = None
        candidates = grna_name if isinstance(grna_name, list) else [grna_name]
        is_shared = not grna_unique and len(candidates) > 1

        # Track shared gRNA hits
        if is_shared:
            counts['shared_grna_hit'] += 1

        # Map variant to stem identity
        try:
            variant_pos = int(variant_str)
        except (ValueError, TypeError):
            variant_pos = None

        # Unique gRNA → resolved immediately
        if grna_unique:
            resolved_name = candidates[0]

        # --- Step 2: iBAR-based disambiguation (shared gRNAs) ---
        if is_shared and ibar_seq != NOT_FOUND:
            ibar_matches = []
            for cand in candidates:
                guide_ibars = grna_to_ibars.get(cand, [])
                ibar_match, ibar_dist = _match_ibar_within_guide(
                    ibar_seq, guide_ibars, max_lev_ibar
                )
                if ibar_match is not None:
                    ibar_matches.append((cand, ibar_match, ibar_dist))

            if len(ibar_matches) == 1:
                resolved_name = ibar_matches[0][0]
                counts['shared_grna_ibar_resolved'] += 1
            elif len(ibar_matches) > 1:
                counts['shared_grna_ambig_ibar'] += 1
            else:
                counts['shared_grna_no_ibar_match'] += 1

        # --- Step 3: Fall back to stem variant (only if unambiguous) ---
        if resolved_name is None and is_shared and variant_str != 'AMBIGUOUS':
            if variant_pos in (1, 2):
                if grna_seq in grna_seq_to_pos_name:
                    pos_map = grna_seq_to_pos_name[grna_seq]
                    if variant_pos in pos_map:
                        resolved_name = pos_map[variant_pos]
                        counts['shared_grna_stem_resolved'] += 1

        if resolved_name is None and is_shared and variant_str == 'AMBIGUOUS':
            counts['shared_grna_ambig_stem'] += 1

        if resolved_name is None:
            return None, 'ambiguous'

        # --- Step 4: iBAR verification & pair_id resolution ---
        if ibar_seq == NOT_FOUND:
            return None, 'decoupled'

        # Use _match_ibar_within_guide for closest unique-best match
        guide_ibars = grna_to_ibars.get(resolved_name, [])
        ibar_match, ibar_dist = _match_ibar_within_guide(
            ibar_seq, guide_ibars, max_lev_ibar
        )
        if ibar_match is None:
            return None, 'decoupled'

        pid = grna_ibar_to_pairid.get((resolved_name, ibar_match))
        if pid is not None:
            return pid, 'ok'
        return None, 'decoupled'

    counts = Counter()
    pos1_pairs, pos2_pairs = [], []
    
    import functools
    @functools.lru_cache(maxsize=100000)
    def _cached_resolve(gseq, sname, iseq, pos):
        return _resolve_guide_unit(gseq, sname, iseq, pos)

    for i in range(N_valid):
        p1_id, p1_stat = _cached_resolve(pos1_gseqs[i], pos1_variants[i], pos1_ibars[i], 1)
        p2_id, p2_stat = _cached_resolve(pos2_gseqs[i], pos2_variants[i], pos2_ibars[i], 2)
        pos1_pairs.append(p1_id)
        pos2_pairs.append(p2_id)
        
        if p1_stat == 'grna_deletion': counts['grna_del1'] += 1
        elif p1_stat == 'ambiguous': counts['ambiguous_pos1'] += 1
        elif p1_stat == 'decoupled': counts['decoupled_pos1'] += 1
        if p2_stat == 'grna_deletion': counts['grna_del2'] += 1
        elif p2_stat == 'ambiguous': counts['ambiguous_pos2'] += 1
        elif p2_stat == 'decoupled': counts['decoupled_pos2'] += 1
    
    for i in range(N_valid):
        p1_id = pos1_pairs[i]
        p2_id = pos2_pairs[i]
        if p1_id is not None and p2_id is not None:
            if p1_id == p2_id:
                counts['intact'] += 1
                counts['intact_relaxed'] += 1  # Standard intact is also relaxed intact
            else:
                g1 = pairid_to_guides.get(p1_id, {}).get(1)
                exp_p1 = pairid_to_guides.get(p2_id, {}).get(1)
                g2 = pairid_to_guides.get(p2_id, {}).get(2)
                exp_p2 = pairid_to_guides.get(p1_id, {}).get(2)
                
                # Check relaxed intactness (same target, ignore stem version)
                p1_target = grna_to_target.get(g1) if g1 else None
                p2_target = grna_to_target.get(g2) if g2 else None
                if p1_target is not None and p2_target is not None and p1_target == p2_target:
                    counts['intact_relaxed'] += 1

                if g1 and exp_p1: _classify_position_swap(g1, exp_p1, grna_to_target, counts, 'pos1')
                if g2 and exp_p2: _classify_position_swap(g2, exp_p2, grna_to_target, counts, 'pos2')
                if grna_to_position.get(g1) == 2 and grna_to_position.get(g2) == 1:
                    counts['pos1_pos2_swap'] += 1

    denom, denom_total = max(N_valid, 1), max(N_total, 1)
    rates = {
        'pos1_fail': N_pos1_fail / denom_total,
        'pos2_fail': N_pos2_fail / denom_total,
        'grna_del1': counts.get('grna_del1', 0) / denom,
        'grna_del2': counts.get('grna_del2', 0) / denom,
        'decoupled_pos1': counts.get('decoupled_pos1', 0) / denom,
        'decoupled_pos2': counts.get('decoupled_pos2', 0) / denom,
        'ambiguous_pos1': counts.get('ambiguous_pos1', 0) / denom,
        'ambiguous_pos2': counts.get('ambiguous_pos2', 0) / denom,
        'swap_within_guide_pos1': counts.get('swap_within_guide_pos1', 0) / denom,
        'swap_within_guide_pos2': counts.get('swap_within_guide_pos2', 0) / denom,
        'swap_between_guide_pos1': counts.get('swap_between_guide_pos1', 0) / denom,
        'swap_between_guide_pos2': counts.get('swap_between_guide_pos2', 0) / denom,
        'swap_across_target_pos1': counts.get('swap_across_target_pos1', 0) / denom,
        'swap_across_target_pos2': counts.get('swap_across_target_pos2', 0) / denom,
        'pos1_pos2_swap': counts.get('pos1_pos2_swap', 0) / denom,
        'intact_rate': counts.get('intact', 0) / denom,
        'intact_with_ambiguous_stem': counts.get('intact_relaxed', 0) / denom,
        'N_total': N_total,
        'N_valid': N_valid,
        **{f'n_{k}': v for k, v in counts.items()},
    }
    
    if verbose:
        print(f"    Stage 4: intact={counts.get('intact',0):,} ({(counts.get('intact',0)/denom)*100:.2f}%)")
        print(f"             intact_relaxed={counts.get('intact_relaxed',0):,} ({(counts.get('intact_relaxed',0)/denom)*100:.2f}%)")
        print(f"             ambiguous1={counts.get('ambiguous_pos1', 0):,}  ambiguous2={counts.get('ambiguous_pos2', 0):,}")
        print(f"             decoupled1={counts.get('decoupled_pos1', 0):,}  decoupled2={counts.get('decoupled_pos2', 0):,}")
        # Shared gRNA diagnostics
        n_shared = counts.get('shared_grna_hit', 0)
        if n_shared > 0:
            n_exact = counts.get('exact_combined', 0)
            n_ibar = counts.get('shared_grna_ibar_resolved', 0)
            n_stem = counts.get('shared_grna_stem_resolved', 0)
            n_ambig_ibar = counts.get('shared_grna_ambig_ibar', 0)
            n_ambig_stem = counts.get('shared_grna_ambig_stem', 0)
            n_no_ibar = counts.get('shared_grna_no_ibar_match', 0)
            print(f"      --- Shared gRNA resolution (N={n_shared:,}) ---")
            print(f"      exact_combined: {n_exact:,}  ibar_resolved: {n_ibar:,}  stem_resolved: {n_stem:,}")
            print(f"      ambig_ibar: {n_ambig_ibar:,}  ambig_stem: {n_ambig_stem:,}  no_ibar_match: {n_no_ibar:,}")

    return {'rates': rates, 'stem_report': stem_report, 'counts': dict(counts),
            'pos1_pairs': pos1_pairs, 'pos2_pairs': pos2_pairs,
            'pos1_grna_names': [_cached_resolve(pos1_gseqs[i], pos1_variants[i], pos1_ibars[i], 1) for i in range(N_valid)],
            'pos2_grna_names': [_cached_resolve(pos2_gseqs[i], pos2_variants[i], pos2_ibars[i], 2) for i in range(N_valid)]}


# ============================================================================
# Plotting
# ============================================================================

def plot_guideseq_qc(per_lib_rates, merged_rates, out_dir, suffix=''):
    """Generate QC plots for guide-seq analysis.

    Args:
        per_lib_rates: list of (name, rates_dict) tuples
        merged_rates: merged rates dict
        out_dir: output directory
    """
    # --- Rate comparison bar chart ---
    rate_keys = ['intact_rate', 'grna_del1', 'grna_del2',
                 'swap_within_guide_pos1', 'swap_within_guide_pos2',
                 'swap_between_guide_pos1', 'swap_between_guide_pos2',
                 'swap_across_target_pos1', 'swap_across_target_pos2',
                 'pos1_pos2_swap', 'pos1_fail', 'pos2_fail']
    rate_labels = ['Intact', 'gRNA\ndel1', 'gRNA\ndel2',
                   'Within\nguide\npos1', 'Within\nguide\npos2',
                   'Between\nguide\npos1', 'Between\nguide\npos2',
                   'Across\ntarget\npos1', 'Across\ntarget\npos2',
                   'pos1\u2194pos2\nswap', 'Pos1\nfail', 'Pos2\nfail']

    n_groups = len(per_lib_rates) + 1  # per-lib + merged
    fig, ax = plt.subplots(figsize=(14, 5))

    x = np.arange(len(rate_keys))
    width = 0.8 / n_groups

    colors = plt.cm.Set2(np.linspace(0, 1, max(n_groups, 3)))

    for i, (name, rates) in enumerate(per_lib_rates):
        vals = [rates.get(k, 0) * 100 for k in rate_keys]
        ax.bar(x + i * width, vals, width, label=name, color=colors[i], alpha=0.8)

    # Merged
    vals = [merged_rates.get(k, 0) * 100 for k in rate_keys]
    ax.bar(x + len(per_lib_rates) * width, vals, width, label='Merged',
           color='#333333', alpha=0.9, edgecolor='white')

    ax.set_xticks(x + width * (n_groups - 1) / 2)
    ax.set_xticklabels(rate_labels, fontsize=8)
    ax.set_ylabel('Rate (%)', fontsize=12)
    ax.set_title('Guide-Seq QC: 4-Stage Matching Rates', fontsize=13)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.2, axis='y')

    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, f'q1_guideseq_rates_{suffix}.png'),
                dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved: q1_guideseq_rates_{suffix}.png")

    # --- 3-level swap breakdown per position ---
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    for idx, pos_label in enumerate(['pos1', 'pos2']):
        ax = axes[idx]
        k_within = f'swap_within_guide_{pos_label}'
        k_between = f'swap_between_guide_{pos_label}'
        k_across = f'swap_across_target_{pos_label}'

        within_vals = [r.get(k_within, 0) * 100 for _, r in per_lib_rates]
        between_vals = [r.get(k_between, 0) * 100 for _, r in per_lib_rates]
        across_vals = [r.get(k_across, 0) * 100 for _, r in per_lib_rates]
        within_vals.append(merged_rates.get(k_within, 0) * 100)
        between_vals.append(merged_rates.get(k_between, 0) * 100)
        across_vals.append(merged_rates.get(k_across, 0) * 100)

        names = [n for n, _ in per_lib_rates] + ['Merged']
        x = np.arange(len(names))

        ax.bar(x, within_vals, 0.35, label='Within guide', color='#4CAF50', alpha=0.8)
        bottom_1 = [w for w in within_vals]
        ax.bar(x, between_vals, 0.35, bottom=bottom_1,
               label='Between guide', color='#2196F3', alpha=0.8)
        bottom_2 = [b1 + b2 for b1, b2 in zip(bottom_1, between_vals)]
        ax.bar(x, across_vals, 0.35, bottom=bottom_2,
               label='Across target', color='#FF5722', alpha=0.8)

        ax.set_xticks(x)
        ax.set_xticklabels(names, fontsize=9, rotation=30, ha='right')
        ax.set_ylabel('Rate (%)', fontsize=10)
        ax.set_title(f'{pos_label}: 3-level swap breakdown', fontsize=11)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.2, axis='y')

    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, f'q1_guideseq_swap_breakdown_{suffix}.png'),
                dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved: q1_guideseq_swap_breakdown_{suffix}.png")

    # --- Construct-level pie chart (merged) ---
    fig, ax = plt.subplots(figsize=(7, 7))

    cat_labels = ['Intact', 'pos1 swap', 'pos2 swap', 'Cross swap',
                  'pos1 del', 'pos2 del', 'Both swapped', 'Stem swap']
    cat_keys = ['n_intact', 'n_pos1_swap', 'n_pos2_swap', 'n_cross_swap',
                'n_pos1_deleted', 'n_pos2_deleted', 'n_both_swapped', 'n_stem_swap']
    cat_vals = [merged_rates.get(k, 0) for k in cat_keys]

    # Filter out zeros
    nonzero = [(l, v) for l, v in zip(cat_labels, cat_vals) if v > 0]
    if nonzero:
        labels, vals = zip(*nonzero)
        colors_pie = plt.cm.Set3(np.linspace(0, 1, len(labels)))
        ax.pie(vals, labels=labels, autopct='%1.1f%%', colors=colors_pie,
               startangle=140, textprops={'fontsize': 9})
        ax.set_title(f'Construct Classification (N={merged_rates.get("n_valid", 0):,})',
                     fontsize=12)

    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, f'q1_guideseq_construct_pie_{suffix}.png'),
                dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved: q1_guideseq_construct_pie_{suffix}.png")


# ============================================================================
# Frequency Table Extraction (for Q4 fitness analysis)
# ============================================================================

def get_guide_frequencies(match_results, ref_df):
    """Extract per-guide frequencies from intact constructs only.

    Only constructs where both positions resolve to the same pair_id
    (i.e., intact) contribute to the frequency table.

    Args:
        match_results: dict returned by match_and_compute_from_tsv
        ref_df: reference DataFrame

    Returns:
        DataFrame with columns: [guide_name, target_gene, position, count, frequency]
    """
    pos1_pairs = match_results.get('pos1_pairs', [])
    pos2_pairs = match_results.get('pos2_pairs', [])

    # Build pair_id → guide info lookup
    pairid_to_pos_guide = {}  # pair_id → {position: grna_name}
    grna_to_info = {}  # grna_name → (target_gene, position)
    for _, row in ref_df.iterrows():
        pid = row['pair_id']
        pos = int(row['position'])
        pairid_to_pos_guide.setdefault(pid, {})[pos] = row['grna_name']
        grna_to_info[row['grna_name']] = (row['target_gene'], pos)

    # Count intact constructs per guide
    guide_counts = Counter()
    for i in range(len(pos1_pairs)):
        p1 = pos1_pairs[i]
        p2 = pos2_pairs[i]
        if p1 is not None and p2 is not None and p1 == p2:
            # Intact construct — count both guides
            for pos in [1, 2]:
                gname = pairid_to_pos_guide.get(p1, {}).get(pos)
                if gname:
                    guide_counts[gname] += 1

    total = sum(guide_counts.values())
    rows = []
    for gname, count in sorted(guide_counts.items()):
        tgene, pos = grna_to_info.get(gname, ('', 0))
        rows.append({
            'guide_name': gname,
            'target_gene': tgene,
            'position': pos,
            'count': count,
            'frequency': count / total if total > 0 else 0,
        })

    return pd.DataFrame(rows)


def get_combo_frequencies(match_results, ref_df):
    """Extract per-combo frequencies from intact constructs only.

    Only constructs where both positions resolve to the same pair_id
    (i.e., intact) contribute to the frequency table.

    Args:
        match_results: dict returned by match_and_compute_from_tsv
        ref_df: reference DataFrame

    Returns:
        DataFrame with columns: [pair_id, pos1_guide, pos2_guide, target_combo, count, frequency]
    """
    pos1_pairs = match_results.get('pos1_pairs', [])
    pos2_pairs = match_results.get('pos2_pairs', [])

    # Build pair_id → info lookup
    pairid_info = {}
    for pid, grp in ref_df.groupby('pair_id'):
        p1 = grp[grp['position'] == 1]
        p2 = grp[grp['position'] == 2]
        pairid_info[pid] = {
            'pos1_guide': p1.iloc[0]['grna_name'] if len(p1) else '',
            'pos2_guide': p2.iloc[0]['grna_name'] if len(p2) else '',
            'target_combo': p1.iloc[0].get('target_combo', '') if len(p1) else '',
        }

    # Count intact combos
    combo_counts = Counter()
    for i in range(len(pos1_pairs)):
        p1 = pos1_pairs[i]
        p2 = pos2_pairs[i]
        if p1 is not None and p2 is not None and p1 == p2:
            combo_counts[p1] += 1

    total = sum(combo_counts.values())
    rows = []
    for pid, count in sorted(combo_counts.items(), key=lambda x: -x[1]):
        info = pairid_info.get(pid, {})
        rows.append({
            'pair_id': pid,
            'pos1_guide': info.get('pos1_guide', ''),
            'pos2_guide': info.get('pos2_guide', ''),
            'target_combo': info.get('target_combo', ''),
            'count': count,
            'frequency': count / total if total > 0 else 0,
        })

    return pd.DataFrame(rows)


# ============================================================================
# Main Entry Point
# ============================================================================

def run_guideseq_qc(fastq_pairs, ref_df, guideseq_cfg, out_dir,
                    max_lev_dist=3, max_lev_dist_ibar=3,
                    max_reads=0, n_workers=1, verbose=True):
    """Run guide-seq QC analysis: extraction → matching → rates.

    Stage 1: Extract elements from each library, write TSV + FASTA files.
    Stage 2: Match extracted elements to reference, classify constructs.

    Args:
        fastq_pairs: list of (name, r1_path, r2_path) tuples
        ref_df: reference DataFrame (from grna_ibar_reference.csv)
        guideseq_cfg: dict with guide-seq extraction config
        out_dir: output directory
        max_lev_dist: max Levenshtein for gRNA matching
        max_lev_dist_ibar: max Levenshtein for iBAR matching
        max_reads: max reads per library (0 = all)
        n_workers: number of parallel workers for extraction
        verbose: print details

    Returns:
        dict with per_lib_rates, merged_rates
    """
    try:
        from scripts.guide_design import enabled
        from scripts.guide_design_qc import run_guideseq
    except ModuleNotFoundError:
        from guide_design import enabled
        from guide_design_qc import run_guideseq
    design_config = {'guide_design': guideseq_cfg.get('guide_design', {}), 'grna_structure': guideseq_cfg}
    if enabled(design_config):
        return run_guideseq(fastq_pairs, ref_df, guideseq_cfg, out_dir, find_anchor, _read_fastq_pairs,
                            max_lev_dist, max_lev_dist_ibar, max_reads, matcher_factory=NeighborhoodMatcher)
    os.makedirs(out_dir, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"  QC Q1: CSM Guide-Seq Bulk Library Analysis")
    print(f"{'='*60}")
    print(f"  Libraries: {len(fastq_pairs)}")
    print(f"  Matching: max_lev_dist={max_lev_dist}, max_lev_dist_ibar={max_lev_dist_ibar}")

    # Build matchers once (shared across libraries)
    print(f"  Building matchers...")
    matchers = build_matchers(ref_df, max_lev_dist, max_lev_dist_ibar)

    per_lib_extraction = []
    per_lib_rates = []

    for name, r1_path, r2_path in fastq_pairs:
        lib_out_dir = os.path.join(out_dir, name)

        print(f"\n  --- Library: {name} ---")
        print(f"    R1: {r1_path}")
        print(f"    R2: {r2_path}")

        # Stage 1: Extraction → TSV + FASTA
        extraction_result = process_guideseq_library(
            r1_path, r2_path, guideseq_cfg, lib_out_dir,
            max_reads=max_reads, n_workers=n_workers, verbose=verbose,
        )
        per_lib_extraction.append((name, extraction_result))

        # Stage 2: 4-stage matching from TSVs (both independent and construct)
        for tsv_label, tsv_path in [
            ('independent', extraction_result['indep_tsv_path']),
            ('construct', extraction_result['construct_tsv_path']),
        ]:
            print(f"\n    Matching ({tsv_label}):")
            result = match_and_compute_from_tsv(
                tsv_path, ref_df, matchers, verbose=verbose)

            rates = result['rates']
            rates['source'] = tsv_label

            # Add extraction stats
            ext_stats = extraction_result['stats']
            merge_total = ext_stats.get('merge_success', 0) + ext_stats.get('merge_failed', 0)
            rates['merge_rate'] = ext_stats.get('merge_success', 0) / max(merge_total, 1)
            rates['n_total_reads'] = ext_stats.get('n_total', 0)

            per_lib_rates.append((f"{name}_{tsv_label}", rates))

            if verbose:
                print(f"      Intact:            {rates.get('intact_rate', 0):.4f} ({rates.get('intact_rate', 0)*100:.2f}%)")
                print(f"      pos1_fail:         {rates.get('pos1_fail', 0):.4f}")
                print(f"      pos2_fail:         {rates.get('pos2_fail', 0):.4f}")
                print(f"      grna_del1:         {rates.get('grna_del1', 0):.4f}")
                print(f"      grna_del2:         {rates.get('grna_del2', 0):.4f}")
                print(f"      --- Swap rates ---")
                print(f"      within_guide:  pos1={rates.get('swap_within_guide_pos1', 0):.4f}  "
                      f"pos2={rates.get('swap_within_guide_pos2', 0):.4f}")
                print(f"      between_guide: pos1={rates.get('swap_between_guide_pos1', 0):.4f}  "
                      f"pos2={rates.get('swap_between_guide_pos2', 0):.4f}")
                print(f"      across_target: pos1={rates.get('swap_across_target_pos1', 0):.4f}  "
                      f"pos2={rates.get('swap_across_target_pos2', 0):.4f}")
                print(f"      pos1_pos2_swap:    {rates.get('pos1_pos2_swap', 0):.4f}")

            # Save per-lib per-source CSV
            rate_df = pd.DataFrame([rates])
            rate_df.to_csv(os.path.join(lib_out_dir, f'q1_rates_{tsv_label}.csv'), index=False)

    # Save combined summary CSV
    all_rates = []
    for label, rates in per_lib_rates:
        row = {'library': label, **rates}
        all_rates.append(row)

    summary_df = pd.DataFrame(all_rates)
    summary_path = os.path.join(out_dir, 'q1_guideseq_summary.csv')
    summary_df.to_csv(summary_path, index=False)
    print(f"\n  Saved summary: {summary_path}")

        # Plots
    independent_rates = [(l, r) for l, r in per_lib_rates if l.endswith('_independent')]
    construct_rates = [(l, r) for l, r in per_lib_rates if l.endswith('_construct')]

    merged_indep = independent_rates[0][1] if independent_rates else {}
    merged_const = construct_rates[0][1] if construct_rates else {}

    plot_guideseq_qc(independent_rates, merged_indep, out_dir, suffix='independent')
    plot_guideseq_qc(construct_rates, merged_const, out_dir, suffix='construct')


    # Compute frequency tables from construct-path matching results
    # (use the last construct result as the primary source)
    construct_match_results = []
    for name, r1_path, r2_path in fastq_pairs:
        lib_out_dir_freq = os.path.join(out_dir, name)
        construct_tsv = os.path.join(lib_out_dir_freq, 'extraction_construct.tsv.gz')
        if os.path.exists(construct_tsv):
            freq_result = match_and_compute_from_tsv(
                construct_tsv, ref_df, matchers, verbose=False)
            construct_match_results.append((name, freq_result))

    # Merge frequency tables across libraries
    all_guide_freq = []
    all_combo_freq = []
    for name, freq_result in construct_match_results:
        gf = get_guide_frequencies(freq_result, ref_df)
        cf = get_combo_frequencies(freq_result, ref_df)
        gf['library'] = name
        cf['library'] = name
        all_guide_freq.append(gf)
        all_combo_freq.append(cf)

    if all_guide_freq:
        guide_freq_df = pd.concat(all_guide_freq, ignore_index=True)
        guide_freq_path = os.path.join(out_dir, 'q1_guide_frequencies.csv')
        guide_freq_df.to_csv(guide_freq_path, index=False)
        print(f"  Saved guide frequencies: {guide_freq_path}")
    else:
        guide_freq_df = pd.DataFrame()

    if all_combo_freq:
        combo_freq_df = pd.concat(all_combo_freq, ignore_index=True)
        combo_freq_path = os.path.join(out_dir, 'q1_combo_frequencies.csv')
        combo_freq_df.to_csv(combo_freq_path, index=False)
        print(f"  Saved combo frequencies: {combo_freq_path}")
    else:
        combo_freq_df = pd.DataFrame()

    print(f"\n{'='*60}")
    print(f"  QC Q1 Complete")
    print(f"{'='*60}")

    return {
        'per_lib_rates': per_lib_rates,
        'per_lib_extraction': per_lib_extraction,
        'guide_frequencies': guide_freq_df,
        'combo_frequencies': combo_freq_df,
    }


# ============================================================================
# CLI Entry Point
# ============================================================================

if __name__ == '__main__':
    import argparse
    import sys as _sys

    # Add project root to sys.path so "from scripts.*" imports work
    _project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if _project_root not in _sys.path:
        _sys.path.insert(0, _project_root)

    parser = argparse.ArgumentParser(
        description='Q1: CSM Guide-Seq Bulk Library Analysis')
    parser.add_argument('--config', required=True,
                        help='Path to qc_config_guideseq.yaml (standalone or with sc_config cross-ref)')
    parser.add_argument('--n-workers', type=int, default=8,
                        help='Number of parallel workers (default: 8)')
    args = parser.parse_args()

    # Load config (standalone or cross-referenced)
    from scripts.qc_config_loader import load_guideseq_config
    cfg = load_guideseq_config(args.config)

    # Load reference
    ref_path = cfg['paths']['grna_reference']
    if not os.path.isabs(ref_path):
        ref_path = os.path.join(os.path.dirname(args.config), ref_path)
    ref_df = pd.read_csv(ref_path, comment='#')

    # FASTQ directory
    fastq_dir = cfg['paths'].get('guideseq_fastq_dir', '')
    if not fastq_dir:
        print('ERROR: guideseq_fastq_dir not set in config paths')
        _sys.exit(1)

    fastq_pairs = detect_fastq_pairs(fastq_dir)
    print(f'Found {len(fastq_pairs)} FASTQ pair(s)')

    # Guide-seq config
    guideseq_cfg = cfg.get('guideseq', {})

    # Output directory
    qc_dir = cfg['paths'].get('qc_dir', 'qc')
    if not os.path.isabs(qc_dir):
        qc_dir = os.path.join(os.path.dirname(args.config), qc_dir)
    out_dir = os.path.join(qc_dir, 'q1_guideseq')

    # Processing params
    qc_proc = cfg.get('qc_processing', {})

    results = run_guideseq_qc(
        fastq_pairs, ref_df, guideseq_cfg, out_dir,
        max_lev_dist=cfg.get('grna_matching', {}).get('max_lev_dist', 3),
        max_lev_dist_ibar=cfg.get('grna_matching', {}).get('max_lev_dist_ibar', 3),
        max_reads=qc_proc.get('max_reads', 0),
        n_workers=args.n_workers,
    )

    # Print summary
    if results.get('per_lib_rates'):
        # Use first construct rate as summary
        construct_rates = [(l, r) for l, r in results['per_lib_rates'] if l.endswith('_construct')]
        if construct_rates:
            rates = construct_rates[0][1]
            print(f'\nIntact rate: {rates.get("intact_rate", 0):.4f}')
            print(f'Pos1 fail: {rates.get("pos1_fail", 0):.4f}  Pos2 fail: {rates.get("pos2_fail", 0):.4f}')
            print(f'gRNA del1: {rates.get("grna_del1", 0):.4f}  del2: {rates.get("grna_del2", 0):.4f}')
            print(f'Merge rate: {rates.get("merge_rate", 0):.1%}')

    # Print frequency table summary
    gf = results.get('guide_frequencies')
    cf = results.get('combo_frequencies')
    if gf is not None and not gf.empty:
        print(f'\nGuide frequencies: {len(gf)} guides')
    if cf is not None and not cf.empty:
        print(f'Combo frequencies: {len(cf)} combos')


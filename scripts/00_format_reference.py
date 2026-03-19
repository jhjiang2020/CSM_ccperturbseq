#!/usr/bin/env python3
"""
Step 0: Format gRNA/iBAR reference from library design file.

Converts the library design file into the reference lookup table used by the
downstream pipeline. Supports both dual-guide mode (two entries per design row)
and single-guide mode (one entry per design row).

Mode is auto-resolved from config (single_guide_mode: true → single, else dual)
but can be overridden with --mode.

Key behaviors:
  - iBAR sequences in the design are on the oligo (sense) strand. In the
    gRNA read, iBARs appear as the reverse complement. By default this
    script reverse-complements iBARs; disable with --no-rev-comp.
  - When ibar_length is 0 (no iBAR), all iBAR-related checks are skipped.
    If duplicate grna_seq entries are detected, they are collapsed to the
    first grna_name and the user is notified.
  - Guide ID format: <spacer_X_target>_<target_posX>_<spacer_version_posX>
  - Pair-level annotation: from iBAR (and spacer+iBAR combination) we can
    infer <target>_<target_version> for the combination of target_1 + target_2.

Usage:
  python scripts/00_format_reference.py --config config.yaml
  python scripts/00_format_reference.py --config config.yaml --mode single
  python scripts/00_format_reference.py --design design.csv --output ref.csv
"""

import argparse
import os
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import pandas as pd
import yaml


# ---------------------------------------------------------------------------
# Reverse complement
# ---------------------------------------------------------------------------

COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(COMPLEMENT)[::-1]


# ---------------------------------------------------------------------------
# Core formatting
# ---------------------------------------------------------------------------

def format_reference(design_df, rev_comp_ibar=True, ibar_length=None):
    """Convert a library design DataFrame to gRNA/iBAR reference table.

    For each design row, creates TWO entries:
      1. spacer_1 + iBAR_1 (position 1, variant A / stem_1 + tracr_1)
      2. spacer_2 + iBAR_2 (position 2, variant B / stem_2 + tracr_2)

    Args:
        design_df: DataFrame with design columns (see REQUIRED_COLS)
        rev_comp_ibar: If True, reverse-complement iBAR from design file
                       (design iBARs are oligo-strand; reads are rev comp)

    Returns:
        DataFrame with reference columns
    """
    REQUIRED_COLS = [
        'spacer_1', 'spacer_2',
        'target_pos1', 'spacer_version_pos1',
        'target_pos2', 'spacer_version_pos2',
        'target', 'target_version',
        'spacer_1_target', 'spacer_2_target',
    ]
    missing = [c for c in REQUIRED_COLS if c not in design_df.columns]
    if missing:
        raise ValueError(
            f"Design file missing required columns: {missing}\n"
            f"Required: {REQUIRED_COLS}\n"
            f"Found: {list(design_df.columns)}"
        )
        
    has_ibar = 'iBAR_1' in design_df.columns and 'iBAR_2' in design_df.columns
    if not has_ibar:
        print("  WARNING: iBAR_1 or iBAR_2 not found in dual-guide design. Proceeding with empty iBAR sequences.")
        print("  → Make sure to set `ibar_length: 0` in your config.yaml.")

    records = []

    for idx, row in design_df.iterrows():
        # --- Shared pair-level info ---
        target = row['target']
        target_version = row['target_version']
        pair_id = f"{target}_{target_version}"

        # Optional columns (may or may not exist)
        gi_pair = row.get('GI_pair', '')
        gi_version = row.get('GI_version', '')
        gi_version_index = row.get('GI_version_index', '')
        trna_version = row.get('tRNA_version', '')
        dialout = row.get('dialout', '')

        # --- Position 1: spacer_1 + iBAR_1 + stem_1/tracr_1 ---
        ibar_1_raw = str(row.get('iBAR_1', '')).strip() if pd.notna(row.get('iBAR_1', '')) else ''
        ibar_1 = reverse_complement(ibar_1_raw) if rev_comp_ibar and ibar_1_raw else ibar_1_raw

        grna_name_1 = f"{row['spacer_1_target']}_{row['target_pos1']}_{row['spacer_version_pos1']}"
        records.append({
            'grna_name': grna_name_1,
            'grna_seq': str(row['spacer_1']).strip(),
            'ibar_seq': ibar_1,
            'ibar_design_seq': ibar_1_raw,
            'target_gene': str(row['spacer_1_target']).strip(),
            'position': 1,
            'stem_variant': 'A',  # stem_1 + tracr_1
            'pair_id': pair_id,
            'target_combo': target,
            'target_version': target_version,
            'gi_pair': gi_pair,
            'gi_version': gi_version,
            'gi_version_index': gi_version_index,
            'trna_version': trna_version,
            'dialout': dialout,
            'design_row_index': idx,
        })

        # --- Position 2: spacer_2 + iBAR_2 + stem_2/tracr_2 ---
        ibar_2_raw = str(row.get('iBAR_2', '')).strip() if pd.notna(row.get('iBAR_2', '')) else ''
        ibar_2 = reverse_complement(ibar_2_raw) if rev_comp_ibar and ibar_2_raw else ibar_2_raw

        grna_name_2 = f"{row['spacer_2_target']}_{row['target_pos2']}_{row['spacer_version_pos2']}"
        records.append({
            'grna_name': grna_name_2,
            'grna_seq': str(row['spacer_2']).strip(),
            'ibar_seq': ibar_2,
            'ibar_design_seq': ibar_2_raw,
            'target_gene': str(row['spacer_2_target']).strip(),
            'position': 2,
            'stem_variant': 'B',  # stem_2 + tracr_2
            'pair_id': pair_id,
            'target_combo': target,
            'target_version': target_version,
            'gi_pair': gi_pair,
            'gi_version': gi_version,
            'gi_version_index': gi_version_index,
            'trna_version': trna_version,
            'dialout': dialout,
            'design_row_index': idx,
        })

    ref_df = pd.DataFrame(records)

    # --- Validate and fix pair_id uniqueness ---
    # Each pair_id must correspond to exactly one pos1 row and one pos2 row.
    # If the design file gives two different clones the same target + target_version,
    # they will collide on pair_id and cause downstream errors in Step 4.
    # Fix: append a unique suffix (_d1, _d2, ...) to each duplicate pair_id.
    design_rows_per_pair = ref_df.groupby('pair_id')['design_row_index'].nunique()
    dup_pair_ids = design_rows_per_pair[design_rows_per_pair > 1].index.tolist()

    if dup_pair_ids:
        n_fixed = 0
        print(f"\n  ⚠ DUPLICATE pair_id: {len(dup_pair_ids)} pair_id(s) are shared by "
              f"multiple design rows. Appending unique suffixes:")
        print(f"    {'original pair_id':<30s}  {'target':<20s}  {'target_version':<16s}  → new pair_id")
        print(f"    {'─'*30}  {'─'*20}  {'─'*16}  {'─'*30}")
        for pid in sorted(dup_pair_ids):
            mask = ref_df['pair_id'] == pid
            # Get unique design rows for this pair_id, ordered by design_row_index
            design_rows = sorted(ref_df.loc[mask, 'design_row_index'].unique())
            for i, dri in enumerate(design_rows, start=1):
                new_pid = f"{pid}_d{i}"
                row_mask = mask & (ref_df['design_row_index'] == dri)
                ref_df.loc[row_mask, 'pair_id'] = new_pid
                # Get the target/target_version for display
                sample_row = ref_df.loc[row_mask].iloc[0]
                print(f"    {pid:<30s}  {sample_row['target_combo']:<20s}  "
                      f"{str(sample_row['target_version']):<16s}  → {new_pid}")
                n_fixed += 1
        print(f"\n    Fixed {n_fixed} pair_ids across {len(dup_pair_ids)} duplicate group(s).")
        print(f"    Note: target_combo is preserved; only pair_id was made unique.")
    else:
        print(f"  ✓ All pair_ids are unique (each maps to exactly one pos1 + one pos2 entry)")

    # --- Build iBAR → pair_id lookup ---
    if has_ibar:
        # From iBAR alone we should be able to infer the pair (target_target_version)
        ibar_to_pair = ref_df.groupby('ibar_seq')['pair_id'].apply(
            lambda x: x.unique()
        )
        ambiguous_ibars = ibar_to_pair[ibar_to_pair.apply(len) > 1]
        if len(ambiguous_ibars) > 0:
            print(f"  WARNING: {len(ambiguous_ibars)} iBAR sequences map to "
                  f"multiple pair_ids (iBAR alone is not sufficient for pair inference):")
            for ibar, pairs in ambiguous_ibars.head(5).items():
                print(f"    {ibar} → {list(pairs)}")
        else:
            print(f"  ✓ All iBAR sequences uniquely identify a pair_id")

        # --- Build (spacer, iBAR) → pair_id lookup ---
        combo_to_pair = ref_df.groupby(['grna_seq', 'ibar_seq'])['pair_id'].apply(
            lambda x: x.unique()
        )
        ambiguous_combos = combo_to_pair[combo_to_pair.apply(len) > 1]
        if len(ambiguous_combos) > 0:
            print(f"  WARNING: {len(ambiguous_combos)} (spacer, iBAR) combos map "
                  f"to multiple pair_ids")
        else:
            print(f"  ✓ All (spacer, iBAR) combos uniquely identify a pair_id")
    else:
        print("  ✓ Skipping iBAR overlap checks (no iBARs provided)")

    return ref_df


def format_reference_single(design_df, rev_comp_ibar=True, ibar_length=None):
    """Convert a library design DataFrame to gRNA/iBAR reference table for single-guide mode.

    For each design row, creates ONE entry:
      - spacer + iBAR (position 1, variant A / stem_1 + tracr_1)

    When ibar_length==0 and duplicate grna_seq are found, they are collapsed
    to the first grna_name and the user is notified.

    Args:
        design_df: DataFrame with design columns
        rev_comp_ibar: If True, reverse-complement iBAR from design file
        ibar_length: Expected iBAR length from config (0 = no iBAR)

    Returns:
        DataFrame with reference columns
    """
    REQUIRED_COLS = ['sgRNA_ID', 'spacer', 'target']
    missing = [c for c in REQUIRED_COLS if c not in design_df.columns]
    if missing:
        raise ValueError(
            f"Design file missing required columns: {missing}\n"
            f"Required: {REQUIRED_COLS}\n"
            f"Found: {list(design_df.columns)}"
        )

    has_ibar = 'iBAR' in design_df.columns
    if not has_ibar:
        print("  WARNING: iBAR not found in single-guide design. Proceeding with empty iBAR sequences.")
        print("  → Make sure to set `ibar_length: 0` in your config.yaml.")

    records = []
    for idx, row in design_df.iterrows():
        grna_name = str(row['sgRNA_ID']).strip()
        target = str(row['target']).strip()
        
        ibar_raw = str(row.get('iBAR', '')).strip() if pd.notna(row.get('iBAR', '')) else ''
        ibar = reverse_complement(ibar_raw) if rev_comp_ibar and ibar_raw else ibar_raw

        records.append({
            'grna_name': grna_name,
            'grna_seq': str(row['spacer']).strip(),
            'ibar_seq': ibar,
            'ibar_design_seq': ibar_raw,
            'target_gene': target,
            'position': 1,
            'stem_variant': 'A',
            'pair_id': '',
            'target_combo': target,
            'target_version': row.get('target_version', ''),
            'gi_pair': '',
            'gi_version': '',
            'gi_version_index': '',
            'trna_version': '',
            'dialout': '',
            'design_row_index': idx,
        })
        
    ref_df = pd.DataFrame(records)

    # --- No-iBAR dedup: collapse duplicate grna_seq ---
    if ibar_length is not None and ibar_length == 0:
        dup_mask = ref_df.duplicated(subset=['grna_seq'], keep=False)
        dup_seqs = ref_df[dup_mask]
        if len(dup_seqs) > 0:
            collapsed_info = []
            for seq, group in dup_seqs.groupby('grna_seq'):
                kept = group.iloc[0]['grna_name']
                dropped = group.iloc[1:]['grna_name'].tolist()
                collapsed_info.append((seq[:20], kept, dropped))
            print(f"\n  ⚠ {len(collapsed_info)} duplicated grna_seq detected (no iBAR to distinguish).")
            print(f"    Collapsing to first grna_name:")
            for seq_trunc, kept, dropped in collapsed_info[:10]:
                print(f"      {seq_trunc}... → keep '{kept}', drop {dropped}")
            if len(collapsed_info) > 10:
                print(f"      ... and {len(collapsed_info) - 10} more")
            ref_df = ref_df.drop_duplicates(subset=['grna_seq'], keep='first')
            print(f"    Reference reduced from {len(dup_seqs) + len(ref_df) - len(ref_df)} to {len(ref_df)} entries.")

    return ref_df


# ---------------------------------------------------------------------------
# Levenshtein distance (self-contained for Step 0)
# ---------------------------------------------------------------------------

try:
    import edlib
    HAS_EDLIB = True
except ImportError:
    HAS_EDLIB = False

def _levenshtein(s1, s2, max_dist=None):
    """Compute Levenshtein edit distance with optional early termination."""
    if HAS_EDLIB:
        d = edlib.align(s1, s2, mode="NW", task="distance", k=max_dist if max_dist is not None else -1)["editDistance"]
        if max_dist is not None and d == -1:
            return max_dist + 1
        return d

    if s1 == s2:
        return 0
    len1, len2 = len(s1), len(s2)
    if max_dist is not None and abs(len1 - len2) > max_dist:
        return max_dist + 1
    if len1 == 0:
        return len2
    if len2 == 0:
        return len1

    prev = list(range(len2 + 1))
    for i in range(1, len1 + 1):
        curr = [i] + [0] * len2
        for j in range(1, len2 + 1):
            cost = 0 if s1[i - 1] == s2[j - 1] else 1
            curr[j] = min(curr[j - 1] + 1, prev[j] + 1, prev[j - 1] + cost)
        if max_dist is not None and min(curr) > max_dist:
            return max_dist + 1
        prev = curr
    return prev[len2]


# ---------------------------------------------------------------------------
# Reference overlap / ambiguity check
# ---------------------------------------------------------------------------

def check_reference_overlap(ref_df, max_lev_dist_grna, max_lev_dist_ibar, ibar_length=None):
    """Check for potentially ambiguous gRNA and iBAR sequences at the given
    Levenshtein distance thresholds.

    Two reference sequences A and B are problematic if their pairwise
    Levenshtein distance d(A,B) is small relative to the matching tolerance:

      - CRITICAL (d ≤ max_lev_dist): A read from ref A with ≤max_lev_dist
        errors could match ref B equally well → guaranteed ambiguity zone.
        These reads will be DISCARDED as ambiguous (not mis-assigned).
      - WARNING  (max_lev_dist < d ≤ 2×max_lev_dist): Error balls of A and B
        overlap — some error-prone reads could fall in the overlap region.
        Lower-distance matches are preferred, so most reads still resolve
        uniquely, but a fraction may be lost.
      - SAFE     (d > 2×max_lev_dist): No overlap possible.

    Args:
        ref_df: Reference DataFrame from format_reference().
        max_lev_dist_grna: Max Levenshtein distance for gRNA matching.
        max_lev_dist_ibar: Max Levenshtein distance for iBAR matching.
        ibar_length: Expected iBAR length from config (0 = no iBAR, skip iBAR checks).

    Returns:
        True if there are CRITICAL overlaps, False otherwise.
    """
    has_critical = False
    has_ibar = ibar_length is None or ibar_length > 0

    # ------------------------------------------------------------------
    # 0. Shared gRNA sequences across positions
    # ------------------------------------------------------------------
    print(f"\n  Checking for shared gRNA sequences across positions...")

    # Find grna_seq that appears at multiple positions with different grna_names
    grna_pos = ref_df.drop_duplicates('grna_name').groupby('grna_seq').agg(
        names=('grna_name', list),
        positions=('position', lambda x: sorted(set(x))),
        n_names=('grna_name', 'count'),
    )
    shared = grna_pos[grna_pos['n_names'] > 1]

    if len(shared) > 0:
        n_shared = len(shared)
        affected_entries = ref_df[ref_df['grna_seq'].isin(shared.index)]
        n_affected = len(affected_entries)
        affected_combos = affected_entries['target_combo'].unique()

        print(f"\n    ⚠ {n_shared} gRNA sequence(s) used at BOTH pos1 and pos2:")
        print(f"      Affected entries: {n_affected}/{len(ref_df)} ({100*n_affected/len(ref_df):.1f}%)")
        print(f"      Affected combos: {len(affected_combos)}")
        for seq, row in list(shared.iterrows())[:5]:
            names = row['names']
            print(f"      {seq[:20]}... → {', '.join(names[:4])}"
                  f"{'...' if len(names) > 4 else ''}")
        if n_shared > 5:
            print(f"      ... and {n_shared - 5} more")

        if has_ibar:
            # Check iBAR uniqueness across positions for shared gRNAs
            ibar_overlap_count = 0
            for seq in shared.index:
                entries = ref_df[ref_df['grna_seq'] == seq]
                pos1_ibars = set(entries[entries['position'] == 1]['ibar_seq'])
                pos2_ibars = set(entries[entries['position'] == 2]['ibar_seq'])
                ibar_overlap_count += len(pos1_ibars & pos2_ibars)

            if ibar_overlap_count == 0:
                print(f"\n    ✓ iBARs do NOT overlap across positions for shared gRNAs")
                print(f"      → Step 3 can disambiguate via iBAR + stem variant")
            else:
                has_critical = True
                print(f"\n    ⚠ CRITICAL: {ibar_overlap_count} iBAR(s) overlap across positions!")
                print(f"      → These reads cannot be disambiguated even with iBAR + variant")
        else:
            has_critical = True
            print(f"\n    ⚠ CRITICAL: No iBAR to disambiguate shared gRNA sequences!")
            print(f"      → These reads cannot be assigned to a unique guide.")
    else:
        print(f"    ✓ No gRNA sequences shared across positions")

    # ------------------------------------------------------------------
    # 1. gRNA sequence overlap (Levenshtein distance)
    # ------------------------------------------------------------------
    print(f"\n  Checking gRNA sequence uniqueness (max_lev_dist={max_lev_dist_grna})...")
    if not HAS_EDLIB and len(ref_df) > 1000:
        print(f"    ⚠ WARNING: edlib is not installed. Computing {len(ref_df)**2//2} pairwise distances in pure Python will be slow!")

    unique_grnas = ref_df.drop_duplicates('grna_seq')[['grna_seq', 'grna_name']].values.tolist()
    n_grnas = len(unique_grnas)
    scan_dist = 2 * max_lev_dist_grna  # Check up to error-ball overlap range

    critical_pairs = []
    warning_pairs = []
    
    from collections import defaultdict
    grna_dist_counts = defaultdict(int)
    min_grna_d = float('inf')

    for i in range(n_grnas):
        seq_a, name_a = unique_grnas[i]
        for j in range(i + 1, n_grnas):
            seq_b, name_b = unique_grnas[j]
            d = _levenshtein(seq_a, seq_b) # compute full distance so we can track min and median
            grna_dist_counts[d] += 1
            if d < min_grna_d:
                min_grna_d = d
            
            if d <= max_lev_dist_grna:
                critical_pairs.append((name_a, name_b, d, seq_a, seq_b))
            elif d <= scan_dist:
                warning_pairs.append((name_a, name_b, d, seq_a, seq_b))

    if critical_pairs:
        has_critical = True
        print(f"\n    ⚠ CRITICAL: {len(critical_pairs)} gRNA pair(s) within max_lev_dist={max_lev_dist_grna}")
        print(f"      These reads will be DISCARDED as ambiguous (not mis-assigned).")
        for name_a, name_b, d, seq_a, seq_b in critical_pairs[:10]:
            print(f"      dist={d}: {name_a} ({seq_a}) ↔ {name_b} ({seq_b})")
        if len(critical_pairs) > 10:
            print(f"      ... and {len(critical_pairs) - 10} more")
    else:
        print(f"    ✓ No gRNA pairs within max_lev_dist={max_lev_dist_grna}")

    if warning_pairs:
        print(f"\n    ⚡ WARNING: {len(warning_pairs)} gRNA pair(s) within 2×max_lev_dist={scan_dist}")
        print(f"      Error balls partially overlap — some reads may be lost.")
        for name_a, name_b, d, seq_a, seq_b in warning_pairs[:5]:
            print(f"      dist={d}: {name_a} ↔ {name_b}")
        if len(warning_pairs) > 5:
            print(f"      ... and {len(warning_pairs) - 5} more")

    # Min pairwise distance among all gRNAs
    if n_grnas > 1:
        # compute median from histogram
        total_pairs = n_grnas * (n_grnas - 1) // 2
        median_d = 0
        cumulative = 0
        for md in sorted(grna_dist_counts.keys()):
            cumulative += grna_dist_counts[md]
            if cumulative >= total_pairs / 2:
                median_d = md
                break
                
        print(f"\n    gRNA pairwise distance stats (total {total_pairs} pairs):")
        print(f"      Minimum:  {min_grna_d}")
        print(f"      Median:   {median_d}")
        print(f"      Safe max_lev_dist for gRNA: ≤ {(min_grna_d - 1) // 2 if min_grna_d > 0 else 0} "
              f"(guarantees no overlap)")
        if min_grna_d <= max_lev_dist_grna:
            print(f"      → Current max_lev_dist={max_lev_dist_grna} EXCEEDS safe limit!")

    # ------------------------------------------------------------------
    # 2. iBAR overlap
    # ------------------------------------------------------------------
    if not has_ibar:
        print("\n  Skipping iBAR overlap checks (ibar_length=0)")
    else:
        ibar_scan_dist = 2 * max_lev_dist_ibar

        # --- 2a. Global iBAR distances (across all constructs) ---
        print(f"\n  Checking iBAR uniqueness globally (all constructs)...")

        unique_ibars = ref_df.drop_duplicates('ibar_seq')[['ibar_seq', 'pair_id']].values.tolist()
        n_ibars = len(unique_ibars)

        global_ibar_critical = []
        global_ibar_warning = []
        
        global_ibar_dist_counts = defaultdict(int)
        min_global_ibar_d = float('inf')

        for i in range(n_ibars):
            seq_a, pair_a = unique_ibars[i]
            for j in range(i + 1, n_ibars):
                seq_b, pair_b = unique_ibars[j]
                d = _levenshtein(seq_a, seq_b) # compute full
                global_ibar_dist_counts[d] += 1
                if d < min_global_ibar_d:
                    min_global_ibar_d = d
                
                if d <= max_lev_dist_ibar:
                    global_ibar_critical.append((pair_a, pair_b, d, seq_a, seq_b))
                elif d <= ibar_scan_dist:
                    global_ibar_warning.append((pair_a, pair_b, d, seq_a, seq_b))

        global_pairs_total = n_ibars * (n_ibars - 1) // 2
        if global_pairs_total > 0:
            median_d = 0
            cumulative = 0
            for md in sorted(global_ibar_dist_counts.keys()):
                cumulative += global_ibar_dist_counts[md]
                if cumulative >= global_pairs_total / 2:
                    median_d = md
                    break
                    
            print(f"\n    Global iBAR pairwise distance stats ({n_ibars} unique iBARs, {global_pairs_total} comparisons):")
            print(f"      Minimum:  {min_global_ibar_d}")
            print(f"      Median:   {median_d}")
            if global_ibar_critical:
                print(f"      {len(global_ibar_critical)} pair(s) within dist≤{max_lev_dist_ibar} globally")
                for pa, pb, d, sa, sb in global_ibar_critical[:5]:
                    print(f"        dist={d}: {pa} ({sa}) ↔ {pb} ({sb})")
                if len(global_ibar_critical) > 5:
                    print(f"        ... and {len(global_ibar_critical) - 5} more")
            if global_ibar_warning:
                print(f"      {len(global_ibar_warning)} pair(s) within dist≤{ibar_scan_dist} globally")
        else:
            print(f"    Only {n_ibars} unique iBAR(s) — nothing to compare")

        # --- 2b. Per-gRNA iBAR distances (what actually matters for matching) ---
        print(f"\n  Checking iBAR uniqueness per gRNA (max_lev_dist_ibar={max_lev_dist_ibar})...")
        print(f"    (This is the operationally relevant check — iBAR matching is gRNA-conditional)")

        # Group iBARs by gRNA
        grna_ibars = defaultdict(list)
        for _, row in ref_df.iterrows():
            grna_ibars[row['grna_name']].append((row['ibar_seq'], row['pair_id']))

        ibar_critical = []
        ibar_warning = []
        all_ibar_min_dist = float('inf')
        all_ibar_dist_counts = defaultdict(int)
        all_ibar_total_pairs = 0

        for grna_name, ibars in grna_ibars.items():
            if len(ibars) < 2:
                continue
            for (seq_a, pair_a), (seq_b, pair_b) in combinations(ibars, 2):
                d = _levenshtein(seq_a, seq_b)
                all_ibar_dist_counts[d] += 1
                all_ibar_total_pairs += 1
                if d < all_ibar_min_dist:
                    all_ibar_min_dist = d
                
                if d <= max_lev_dist_ibar:
                    ibar_critical.append((grna_name, pair_a, pair_b, d, seq_a, seq_b))
                elif d <= ibar_scan_dist:
                    ibar_warning.append((grna_name, pair_a, pair_b, d, seq_a, seq_b))

        if ibar_critical:
            has_critical = True
            print(f"\n    ⚠ CRITICAL: {len(ibar_critical)} iBAR pair(s) within max_lev_dist_ibar={max_lev_dist_ibar}")
            print(f"      These reads will be DISCARDED as ambiguous after gRNA is matched.")
            for grna, pa, pb, d, sa, sb in ibar_critical[:10]:
                print(f"      dist={d}: [{grna}] {pa} ({sa}) ↔ {pb} ({sb})")
            if len(ibar_critical) > 10:
                print(f"      ... and {len(ibar_critical) - 10} more")
        else:
            print(f"    ✓ No iBAR pairs within max_lev_dist_ibar={max_lev_dist_ibar} (per gRNA)")

        if ibar_warning:
            print(f"\n    ⚡ WARNING: {len(ibar_warning)} iBAR pair(s) within 2×max_lev_dist_ibar={ibar_scan_dist}")
            for grna, pa, pb, d, sa, sb in ibar_warning[:5]:
                print(f"      dist={d}: [{grna}] {pa} ↔ {pb}")
            if len(ibar_warning) > 5:
                print(f"      ... and {len(ibar_warning) - 5} more")

        if all_ibar_total_pairs > 0:
            median_d = 0
            cumulative = 0
            for md in sorted(all_ibar_dist_counts.keys()):
                cumulative += all_ibar_dist_counts[md]
                if cumulative >= all_ibar_total_pairs / 2:
                    median_d = md
                    break
                    
            print(f"\n    Per-gRNA iBAR pairwise distance stats ({all_ibar_total_pairs} comparisons):")
            print(f"      Minimum:  {all_ibar_min_dist}")
            print(f"      Median:   {median_d}")
            print(f"      Safe max_lev_dist_ibar: ≤ {(all_ibar_min_dist - 1) // 2 if all_ibar_min_dist > 0 else 0} "
                  f"(guarantees no overlap)")
            if all_ibar_min_dist <= max_lev_dist_ibar:
                print(f"      → Current max_lev_dist_ibar={max_lev_dist_ibar} EXCEEDS safe limit!")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    if has_critical:
        print(f"\n  {'='*58}")
        print(f"  ⚠ OVERLAP WARNING: Some gRNAs/iBARs are not uniquely")
        print(f"    resolvable at the configured Levenshtein distances.")
        print(f"    Reads matching these sequences will be marked AMBIGUOUS")
        print(f"    and discarded (no risk of wrong assignment, but data loss).")
        print(f"")
        print(f"    Options:")
        print(f"      1. Reduce max_lev_dist / max_lev_dist_ibar in config.yaml")
        print(f"         to the safe limits shown above.")
        print(f"      2. Proceed as-is (ambiguous reads are discarded,")
        print(f"         check Step 3 QC for '%% ambiguous').")
        print(f"  {'='*58}")
    else:
        print(f"\n    ✓ All sequences are uniquely resolvable at the configured distances.")

    return has_critical


# ---------------------------------------------------------------------------
# QC and reporting
# ---------------------------------------------------------------------------

def print_summary(ref_df):
    """Print a summary of the formatted reference."""
    print(f"\n  Reference Summary:")
    print(f"    Total entries:        {len(ref_df):>8,}")
    print(f"    Position 1 entries:   {(ref_df['position'] == 1).sum():>8,}")
    print(f"    Position 2 entries:   {(ref_df['position'] == 2).sum():>8,}")
    print(f"    Unique gRNA seqs:     {ref_df['grna_seq'].nunique():>8,}")
    print(f"    Unique iBAR seqs:     {ref_df['ibar_seq'].nunique():>8,}")
    print(f"    Unique target genes:  {ref_df['target_gene'].nunique():>8,}")
    print(f"    Unique pair_ids:      {ref_df['pair_id'].nunique():>8,}")
    print(f"    Unique target_combos: {ref_df['target_combo'].nunique():>8,}")

    # gRNA length distribution
    grna_lens = ref_df['grna_seq'].str.len().value_counts().sort_index()
    print(f"\n    gRNA length distribution:")
    for length, count in grna_lens.items():
        print(f"      {length}bp: {count:>8,}")

    # iBAR length distribution
    ibar_lens = ref_df['ibar_seq'].str.len().value_counts().sort_index()
    print(f"\n    iBAR length distribution:")
    for length, count in ibar_lens.items():
        print(f"      {length}bp: {count:>8,}")

    # Check for duplicate guide names
    dup_names = ref_df['grna_name'].value_counts()
    dup_names = dup_names[dup_names > 1]
    if len(dup_names) > 0:
        print(f"\n    WARNING: {len(dup_names)} duplicate guide names found!")
        for name, count in dup_names.head(5).items():
            print(f"      {name}: {count}x")
    else:
        print(f"\n    ✓ All guide names are unique")

    # Check for duplicate (grna_seq, ibar_seq)
    dup_combos = ref_df.duplicated(subset=['grna_seq', 'ibar_seq'], keep=False)
    if dup_combos.sum() > 0:
        print(f"    WARNING: {dup_combos.sum()} rows have duplicate (grna_seq, ibar_seq)")
    else:
        print(f"    ✓ All (grna_seq, ibar_seq) combinations are unique")

    # Target gene distribution (top 10)
    print(f"\n    Target gene distribution (top 10):")
    for gene, count in ref_df['target_gene'].value_counts().head(10).items():
        print(f"      {gene:>30s}: {count:>6,}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Format gRNA/iBAR reference from library design file"
    )
    parser.add_argument("--config", type=str, default=None,
                        help="Path to config.yaml (uses design_file, grna_reference, rev_comp_ibar)")
    parser.add_argument("--design", type=str, default=None,
                        help="Override: path to design CSV/TSV file")
    parser.add_argument("--output", type=str, default=None,
                        help="Override: path to output reference CSV")
    parser.add_argument("--no-rev-comp", action="store_true",
                        help="Do NOT reverse-complement iBAR sequences from design")
    parser.add_argument("--mode", type=str, choices=["dual", "single"], default=None,
                        help="Mode for reference formatting: dual or single. Auto-resolved from config single_guide_mode if not specified.")
    parser.add_argument("--sep", type=str, default=None,
                        help="Column separator for design file (auto-detected by default)")
    parser.add_argument("--check-only", type=str, default=None, metavar="REF_CSV",
                        help="Run overlap check on an existing reference CSV (skip reformatting)")
    parser.add_argument("--max-lev-dist", type=int, default=None,
                        help="Override max Levenshtein distance for gRNA (default: from config or 3)")
    parser.add_argument("--max-lev-dist-ibar", type=int, default=None,
                        help="Override max Levenshtein distance for iBAR (default: from config or 3)")
    parser.add_argument("--skip-check", action="store_true",
                        help="Bypass reference overlap distance check entirely (useful for very large libraries)")
    args = parser.parse_args()

    # --- Check-only mode: just run overlap analysis on existing reference ---
    if args.check_only:
        ref_path = args.check_only
        if not os.path.exists(ref_path):
            print(f"ERROR: Reference file not found: {ref_path}")
            sys.exit(1)

        # Read distance params from config or CLI
        max_lev_grna = args.max_lev_dist or 3
        max_lev_ibar = args.max_lev_dist_ibar or 3
        if args.config:
            with open(args.config) as f:
                cfg = yaml.safe_load(f)
            grna_cfg = cfg.get('grna_matching', {})
            max_lev_grna = args.max_lev_dist or grna_cfg.get('max_lev_dist', 3)
            max_lev_ibar = args.max_lev_dist_ibar or grna_cfg.get('max_lev_dist_ibar', 3)

        print(f"\n{'='*60}")
        print(f"  Reference Overlap Check")
        print(f"{'='*60}")
        print(f"  Reference: {ref_path}")
        ref_df = pd.read_csv(ref_path)
        print(f"  Loaded {len(ref_df)} entries")
        print_summary(ref_df)
        if not args.skip_check:
            has_critical = check_reference_overlap(ref_df, max_lev_grna, max_lev_ibar)
        else:
            print("\n  ⚠ Overlap check bypassed via --skip-check")
            has_critical = False
            
        print(f"\n{'='*60}")
        sys.exit(1 if has_critical else 0)

    # Resolve paths
    if args.config:
        with open(args.config) as f:
            cfg = yaml.safe_load(f)
        design_path = args.design or cfg['design']['design_file']
        output_path = args.output or cfg['paths']['grna_reference']
        rev_comp = not args.no_rev_comp and cfg['design'].get('rev_comp_ibar', True)
    elif args.design:
        design_path = args.design
        output_path = args.output or design_path.replace('.csv', '_reference.csv')
        rev_comp = not args.no_rev_comp
        cfg = None
    else:
        parser.error("Either --config or --design must be provided")
        return

    # Resolve mode: CLI > config > default (dual)
    if args.mode is not None:
        mode = args.mode
        mode_source = "CLI"
    elif args.config:
        single = cfg.get('grna_structure', {}).get('single_guide_mode', False)
        mode = "single" if single else "dual"
        mode_source = "config (single_guide_mode)"
    else:
        mode = "dual"
        mode_source = "default"

    # Resolve ibar_length from config (default 7)
    ibar_length = 7
    if args.config:
        ibar_length = cfg.get('grna_structure', {}).get('ibar_length', 7)

    print(f"\n{'='*60}")
    print(f"  Step 0: Format gRNA/iBAR Reference")
    print(f"{'='*60}")
    print(f"  Mode:          {mode} (from {mode_source})")
    print(f"  iBAR length:   {ibar_length}")
    print(f"  Design file:   {design_path}")
    print(f"  Output:        {output_path}")
    print(f"  Rev-comp iBAR: {rev_comp}")
    print()

    # Read design file
    if not os.path.exists(design_path):
        print(f"ERROR: Design file not found: {design_path}")
        sys.exit(1)

    # Auto-detect separator
    sep = args.sep
    if sep is None:
        with open(design_path) as f:
            first_line = f.readline()
        if '\t' in first_line:
            sep = '\t'
        else:
            sep = ','

    design_df = pd.read_csv(design_path, sep=sep)
    print(f"  Loaded design: {len(design_df)} rows, {len(design_df.columns)} columns")
    print(f"  Columns: {list(design_df.columns)}")

    # Format reference
    if mode == "single":
        ref_df = format_reference_single(design_df, rev_comp_ibar=rev_comp, ibar_length=ibar_length)
    else:
        ref_df = format_reference(design_df, rev_comp_ibar=rev_comp, ibar_length=ibar_length)

    # Print summary
    print_summary(ref_df)

    # Run overlap check if config provides distance parameters
    if args.config and not args.skip_check:
        grna_cfg = cfg.get('grna_matching', {})
        max_lev_grna = grna_cfg.get('max_lev_dist', 3)
        max_lev_ibar = grna_cfg.get('max_lev_dist_ibar', 3)
        check_reference_overlap(ref_df, max_lev_grna, max_lev_ibar, ibar_length=ibar_length)
    elif args.skip_check:
        print("\n  ⚠ Skipped reference overlap check (--skip-check)")

    # Save
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    ref_df.to_csv(output_path, index=False)
    print(f"\n  ✓ Saved reference to {output_path}")

    # Also save the iBAR → pair_id lookup as a separate file for convenience
    if ibar_length > 0:
        ibar_lookup = ref_df[['ibar_seq', 'pair_id', 'target_combo', 'target_version']].drop_duplicates()
        lookup_path = output_path.replace('.csv', '_ibar_lookup.csv')
        ibar_lookup.to_csv(lookup_path, index=False)
        print(f"  ✓ Saved iBAR→pair_id lookup to {lookup_path}")
    else:
        print(f"  ✓ Skipped iBAR→pair_id lookup (ibar_length=0)")

    print(f"\n{'='*60}")
    print(f"  Step 0 complete!")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()

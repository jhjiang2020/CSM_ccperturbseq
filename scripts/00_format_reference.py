#!/usr/bin/env python3
"""
Step 0: Format gRNA/iBAR reference from library design file.

Converts the dual-guide library design file into the reference lookup table
used by the downstream pipeline. Each design row contains a dual-guide
construct (spacer_1 + iBAR_1 paired with stem_1/tracr_1, spacer_2 + iBAR_2
paired with stem_2/tracr_2). This script produces two reference entries
per design row — one per guide position.

Key behaviors:
  - iBAR sequences in the design are on the oligo (sense) strand. In the
    gRNA read, iBARs appear as the reverse complement. By default this
    script reverse-complements iBARs; disable with --no-rev-comp.
  - Guide ID format: <spacer_X_target>_<target_posX>_<spacer_version_posX>
  - Pair-level annotation: from iBAR (and spacer+iBAR combination) we can
    infer <target>_<target_version> for the combination of target_1 + target_2.

Usage:
  python scripts/00_format_reference.py --config config.yaml
  python scripts/00_format_reference.py --design design.csv --output ref.csv
  python scripts/00_format_reference.py --design design.csv --no-rev-comp
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

def format_reference(design_df, rev_comp_ibar=True):
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
        'spacer_1', 'spacer_2', 'iBAR_1', 'iBAR_2',
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
        ibar_1_raw = str(row['iBAR_1']).strip()
        ibar_1 = reverse_complement(ibar_1_raw) if rev_comp_ibar else ibar_1_raw

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
        ibar_2_raw = str(row['iBAR_2']).strip()
        ibar_2 = reverse_complement(ibar_2_raw) if rev_comp_ibar else ibar_2_raw

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

    # --- Build iBAR → pair_id lookup ---
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

    return ref_df


# ---------------------------------------------------------------------------
# Levenshtein distance (self-contained for Step 0)
# ---------------------------------------------------------------------------

def _levenshtein(s1, s2, max_dist=None):
    """Compute Levenshtein edit distance with optional early termination."""
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

def check_reference_overlap(ref_df, max_lev_dist_grna, max_lev_dist_ibar):
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

    Returns:
        True if there are CRITICAL overlaps, False otherwise.
    """
    has_critical = False

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
        print(f"    ✓ No gRNA sequences shared across positions")

    # ------------------------------------------------------------------
    # 1. gRNA sequence overlap (Levenshtein distance)
    # ------------------------------------------------------------------
    print(f"\n  Checking gRNA sequence uniqueness (max_lev_dist={max_lev_dist_grna})...")

    unique_grnas = ref_df.drop_duplicates('grna_seq')[['grna_seq', 'grna_name']].values.tolist()
    n_grnas = len(unique_grnas)
    scan_dist = 2 * max_lev_dist_grna  # Check up to error-ball overlap range

    critical_pairs = []
    warning_pairs = []

    for i in range(n_grnas):
        for j in range(i + 1, n_grnas):
            seq_a, name_a = unique_grnas[i]
            seq_b, name_b = unique_grnas[j]
            d = _levenshtein(seq_a, seq_b, max_dist=scan_dist)
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
        all_grna_dists = []
        for i in range(n_grnas):
            for j in range(i + 1, n_grnas):
                d = _levenshtein(unique_grnas[i][0], unique_grnas[j][0])
                all_grna_dists.append(d)
        min_d = min(all_grna_dists)
        print(f"\n    gRNA pairwise distance stats:")
        print(f"      Minimum:  {min_d}")
        print(f"      Median:   {sorted(all_grna_dists)[len(all_grna_dists)//2]}")
        print(f"      Safe max_lev_dist for gRNA: ≤ {(min_d - 1) // 2} "
              f"(guarantees no overlap)")
        if min_d <= max_lev_dist_grna:
            print(f"      → Current max_lev_dist={max_lev_dist_grna} EXCEEDS safe limit!")

    # ------------------------------------------------------------------
    # 2. iBAR overlap
    # ------------------------------------------------------------------
    ibar_scan_dist = 2 * max_lev_dist_ibar

    # --- 2a. Global iBAR distances (across all constructs) ---
    print(f"\n  Checking iBAR uniqueness globally (all constructs)...")

    unique_ibars = ref_df.drop_duplicates('ibar_seq')[['ibar_seq', 'pair_id']].values.tolist()
    n_ibars = len(unique_ibars)

    global_ibar_dists = []
    global_ibar_critical = []
    global_ibar_warning = []

    for i in range(n_ibars):
        for j in range(i + 1, n_ibars):
            seq_a, pair_a = unique_ibars[i]
            seq_b, pair_b = unique_ibars[j]
            d = _levenshtein(seq_a, seq_b, max_dist=ibar_scan_dist)
            global_ibar_dists.append(d)
            if d <= max_lev_dist_ibar:
                global_ibar_critical.append((pair_a, pair_b, d, seq_a, seq_b))
            elif d <= ibar_scan_dist:
                global_ibar_warning.append((pair_a, pair_b, d, seq_a, seq_b))

    if global_ibar_dists:
        min_d = min(global_ibar_dists)
        print(f"\n    Global iBAR pairwise distance stats ({n_ibars} unique iBARs):")
        print(f"      Minimum:  {min_d}")
        print(f"      Median:   {sorted(global_ibar_dists)[len(global_ibar_dists)//2]}")
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
    all_ibar_min_dists = []

    for grna_name, ibars in grna_ibars.items():
        if len(ibars) < 2:
            continue
        for (seq_a, pair_a), (seq_b, pair_b) in combinations(ibars, 2):
            d = _levenshtein(seq_a, seq_b, max_dist=ibar_scan_dist)
            all_ibar_min_dists.append(d)
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

    if all_ibar_min_dists:
        min_d = min(all_ibar_min_dists)
        print(f"\n    Per-gRNA iBAR pairwise distance stats:")
        print(f"      Minimum:  {min_d}")
        print(f"      Median:   {sorted(all_ibar_min_dists)[len(all_ibar_min_dists)//2]}")
        print(f"      Safe max_lev_dist_ibar: ≤ {(min_d - 1) // 2} "
              f"(guarantees no overlap)")
        if min_d <= max_lev_dist_ibar:
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
    parser.add_argument("--sep", type=str, default=None,
                        help="Column separator for design file (auto-detected by default)")
    parser.add_argument("--check-only", type=str, default=None, metavar="REF_CSV",
                        help="Run overlap check on an existing reference CSV (skip reformatting)")
    parser.add_argument("--max-lev-dist", type=int, default=None,
                        help="Override max Levenshtein distance for gRNA (default: from config or 3)")
    parser.add_argument("--max-lev-dist-ibar", type=int, default=None,
                        help="Override max Levenshtein distance for iBAR (default: from config or 3)")
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
        has_critical = check_reference_overlap(ref_df, max_lev_grna, max_lev_ibar)
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
    else:
        parser.error("Either --config or --design must be provided")
        return

    print(f"\n{'='*60}")
    print(f"  Step 0: Format gRNA/iBAR Reference")
    print(f"{'='*60}")
    print(f"  Design file: {design_path}")
    print(f"  Output:      {output_path}")
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
    ref_df = format_reference(design_df, rev_comp_ibar=rev_comp)

    # Print summary
    print_summary(ref_df)

    # Run overlap check if config provides distance parameters
    if args.config:
        grna_cfg = cfg.get('grna_matching', {})
        max_lev_grna = grna_cfg.get('max_lev_dist', 3)
        max_lev_ibar = grna_cfg.get('max_lev_dist_ibar', 3)
        check_reference_overlap(ref_df, max_lev_grna, max_lev_ibar)

    # Save
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    ref_df.to_csv(output_path, index=False)
    print(f"\n  ✓ Saved reference to {output_path}")

    # Also save the iBAR → pair_id lookup as a separate file for convenience
    ibar_lookup = ref_df[['ibar_seq', 'pair_id', 'target_combo', 'target_version']].drop_duplicates()
    lookup_path = output_path.replace('.csv', '_ibar_lookup.csv')
    ibar_lookup.to_csv(lookup_path, index=False)
    print(f"  ✓ Saved iBAR→pair_id lookup to {lookup_path}")

    print(f"\n{'='*60}")
    print(f"  Step 0 complete!")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()

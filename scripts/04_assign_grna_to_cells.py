#!/usr/bin/env python3
"""
Step 4: Cell-level gRNA assignment.

Collapses read-level gRNA matches (from Step 3) to cell-level assignments
via UMI deduplication and threshold filtering.

Usage:
  python scripts/04_assign_grna_to_cells.py --config config.yaml
"""

import argparse
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy import sparse

try:
    from scripts.guide_design import assign_cells, compile_design, eligible_assignments, enabled
except ModuleNotFoundError:
    from guide_design import assign_cells, compile_design, eligible_assignments, enabled


def load_config(config_path):
    with open(config_path) as f:
        return yaml.safe_load(f)


def _resolve_match_path(output_dir, lib_name):
    """Find match file path, preferring .tsv.gz over .tsv."""
    for ext in ('.tsv.gz', '.tsv'):
        p = os.path.join(output_dir, 'grna_matches', f"{lib_name}_grna_matches{ext}")
        if os.path.exists(p):
            return p
    return None


def _resolve_gex_match_path(output_dir, namespace, lib_name):
    """Find a Step1c match file without changing the Step 3 row schema."""
    for ext in ('.tsv.gz', '.tsv'):
        path = os.path.join(
            output_dir,
            namespace,
            'grna_matches',
            f"{lib_name}_grna_matches{ext}",
        )
        if os.path.exists(path):
            return path
    return None


def stream_library_umi_dedup(match_path, lib_name, sample_name, sublib_index,
                             group_cols, max_umi_edits=1, min_reads_per_umi=1,
                             chunk_size=500_000):
    """Stream one library's match file and perform UMI dedup without loading all into memory.

    Reads the .tsv.gz in chunks, accumulates {group_key: {umi: read_count}},
    optionally filters UMIs below min_reads_per_umi, then applies directional
    UMI collapsing.

    Returns:
        (deduped_df, lib_stats) where lib_stats has read/cell counts and
        reads-per-UMI distribution data for QC plotting.
    """
    import time
    t0 = time.time()

    # Columns we actually need (skip the rest to save memory)
    needed_cols = list(set([c for c in group_cols if c != 'ibar_seq_unresolved'] + ['umi']))
    if 'ibar_seq_unresolved' in group_cols and 'ibar_seq' not in needed_cols:
        needed_cols.append('ibar_seq')

    # Accumulate UMI counts: {group_key_tuple: {umi_string: read_count}}
    umi_accum = defaultdict(Counter)
    total_reads = 0
    unique_cells = set()
    unique_guides = set()
    chunks_read = 0

    for chunk in pd.read_csv(match_path, sep='\t', usecols=needed_cols,
                             chunksize=chunk_size):
        if 'guide_unit_id' in group_cols:
            chunk = chunk.fillna('')
        total_reads += len(chunk)
        chunks_read += 1

        # Prefix the cell barcode with the sublibrary index to match split-pipe combine output
        if sublib_index is not None:
            chunk['bc_wells'] = f"{sublib_index}_" + chunk['bc_wells']

        unique_cells.update(chunk['bc_wells'].unique())
        unique_guides.update(chunk['grna_name'].unique())

        if 'ibar_seq_unresolved' in group_cols:
            chunk['ibar_seq_corrected'] = chunk['ibar_seq_corrected'].fillna('')
            chunk['ibar_seq'] = chunk['ibar_seq'].fillna('')
            unresolved_mask = chunk['ibar_seq_corrected'] == ''
            chunk['ibar_seq_unresolved'] = ''
            chunk.loc[unresolved_mask, 'ibar_seq_unresolved'] = chunk.loc[unresolved_mask, 'ibar_seq']

        # Accumulate UMI counts per group
        for row in chunk.itertuples(index=False):
            key = tuple(getattr(row, c) for c in group_cols)
            umi_accum[key][row.umi] += 1

        if chunks_read % 20 == 0:
            elapsed = time.time() - t0
            print(f"    {lib_name}: {total_reads:,} reads, "
                  f"{len(umi_accum):,} groups ({elapsed:.0f}s)")

    elapsed = time.time() - t0
    print(f"    {lib_name}: {total_reads:,} reads streamed in {elapsed:.1f}s, "
          f"{len(umi_accum):,} groups")

    # --- Collect reads-per-UMI distribution (before filtering) ---
    reads_per_umi_dist = Counter()  # {read_count_value: n_umis_with_that_count}
    total_umis_raw = 0
    for umi_counts in umi_accum.values():
        for rc in umi_counts.values():
            reads_per_umi_dist[rc] += 1
            total_umis_raw += 1

    # --- Filter UMIs below min_reads_per_umi ---
    umis_filtered = 0
    reads_filtered = 0
    if min_reads_per_umi > 1:
        for group_key in list(umi_accum.keys()):
            before_count = len(umi_accum[group_key])
            filtered = {u: c for u, c in umi_accum[group_key].items()
                        if c >= min_reads_per_umi}
            removed = before_count - len(filtered)
            umis_filtered += removed
            reads_filtered += sum(
                c for c in umi_accum[group_key].values() if c < min_reads_per_umi
            )
            if filtered:
                umi_accum[group_key] = Counter(filtered)
            else:
                del umi_accum[group_key]
        print(f"    {lib_name}: min_reads_per_umi={min_reads_per_umi} "
              f"→ {umis_filtered:,} UMIs removed ({reads_filtered:,} reads)")

    # Apply directional UMI collapsing on accumulated counts
    records = []
    total_umis_after_filter = 0
    total_umis_after_collapse = 0
    for group_key, umi_counts in umi_accum.items():
        n_umis_pre = len(umi_counts)
        total_umis_after_filter += n_umis_pre
        if max_umi_edits > 0:
            n_clusters = _directional_umi_collapse(dict(umi_counts),
                                                   max_edits=max_umi_edits)
        else:
            n_clusters = len(umi_counts)
        total_umis_after_collapse += n_clusters

        row = dict(zip(group_cols, group_key))
        row['umi_count'] = n_clusters
        row['umi_count_raw'] = n_umis_pre
        row['read_count'] = sum(umi_counts.values())
        row['sublibrary'] = lib_name
        row['sample_name'] = sample_name
        records.append(row)

    # Free the large accumulator
    del umi_accum

    deduped = pd.DataFrame(records)

    lib_stats = {
        'total_reads': total_reads,
        'unique_cells': len(unique_cells),
        'unique_guides': len(unique_guides),
        'groups': len(deduped),
        'elapsed_s': time.time() - t0,
        'total_umis_raw': total_umis_raw,
        'total_umis_after_filter': total_umis_after_filter,
        'total_umis_after_collapse': total_umis_after_collapse,
        'umis_filtered_low_reads': umis_filtered,
        'reads_filtered_low_reads': reads_filtered,
        'reads_per_umi_dist': dict(reads_per_umi_dist),
    }

    # Compute sequencing saturation: 1 - unique_umis / total_reads
    if total_reads > 0:
        lib_stats['sequencing_saturation'] = round(
            1.0 - total_umis_after_collapse / total_reads, 4
        )
    else:
        lib_stats['sequencing_saturation'] = 0.0

    if max_umi_edits > 0 and len(deduped) > 0:
        collapsed = total_umis_after_filter - total_umis_after_collapse
        lib_stats['umis_collapsed'] = int(collapsed)
        print(f"    {lib_name}: UMI collapsing done — {collapsed:,} UMIs collapsed")

    return deduped, lib_stats


def stream_library_umi_dedup_sources(
    source_paths,
    lib_name,
    sample_name,
    sublib_index,
    group_cols,
    max_umi_edits=1,
    min_reads_per_umi=1,
    require_resolved_ibar_sources=None,
    chunk_size=500_000,
):
    """Collapse read evidence from multiple sources in one UMI accumulator.

    ``source_paths`` is an ordered sequence of ``(source_name, match_path)``
    pairs. Read counts from all sources are combined *before*
    ``min_reads_per_umi`` and directional collapse. Source-specific standalone
    counts and union-cluster provenance are retained in the returned rows.

    In dual/iBAR mode callers may put ``"gex"`` in
    ``require_resolved_ibar_sources``. Those rows must carry a non-empty
    ``ibar_seq_corrected`` from the ordinary Step 3 parser before they are
    eligible for the construct-aware union.
    """
    import time

    t0 = time.time()
    require_resolved_ibar_sources = set(require_resolved_ibar_sources or ())
    source_paths = [(str(source), str(path)) for source, path in source_paths]
    source_names = [source for source, _ in source_paths]
    if len(source_names) != len(set(source_names)):
        raise ValueError(f"Duplicate evidence source names: {source_names}")

    needed_cols = list(
        set([c for c in group_cols if c != 'ibar_seq_unresolved'] + ['umi'])
    )
    if 'ibar_seq_unresolved' in group_cols and 'ibar_seq' not in needed_cols:
        needed_cols.append('ibar_seq')

    # group -> UMI -> total reads; and group -> UMI -> source -> reads
    umi_accum = defaultdict(Counter)
    source_accum = defaultdict(dict)
    source_stats = {source: Counter() for source in source_names}
    unique_cells = set()
    unique_guides = set()

    for source, match_path in source_paths:
        header_cols = pd.read_csv(match_path, sep='\t', nrows=0).columns.tolist()
        missing = [col for col in needed_cols if col not in header_cols]
        if missing:
            raise ValueError(
                f"{source} match file {match_path} is missing columns: {missing}"
            )
        read_cols = list(needed_cols)
        if source in require_resolved_ibar_sources:
            if 'ibar_seq_corrected' not in header_cols:
                raise ValueError(
                    f"{source} evidence requires resolved iBARs, but "
                    f"ibar_seq_corrected is absent from {match_path}"
                )
            if 'ibar_seq_corrected' not in read_cols:
                read_cols.append('ibar_seq_corrected')

        for chunk in pd.read_csv(
            match_path, sep='\t', usecols=read_cols, chunksize=chunk_size
        ):
            if 'guide_unit_id' in group_cols:
                chunk = chunk.fillna('')
            source_stats[source]['input_reads'] += len(chunk)
            if source in require_resolved_ibar_sources:
                corrected = chunk['ibar_seq_corrected'].fillna('').astype(str)
                if 'guide_unit_id' in group_cols:
                    eligible = chunk['ibar_status'].isin(['present_resolved', 'absent_by_design']) & chunk['structure_status'].eq('concordant')
                    chunk = chunk[eligible].copy()
                else:
                    chunk = chunk[corrected != ''].copy()
            source_stats[source]['eligible_reads'] += len(chunk)
            if len(chunk) == 0:
                continue

            if sublib_index is not None:
                chunk['bc_wells'] = f"{sublib_index}_" + chunk['bc_wells'].astype(str)

            if 'ibar_seq_unresolved' in group_cols:
                chunk['ibar_seq_corrected'] = chunk['ibar_seq_corrected'].fillna('')
                chunk['ibar_seq'] = chunk['ibar_seq'].fillna('')
                unresolved_mask = chunk['ibar_seq_corrected'] == ''
                chunk['ibar_seq_unresolved'] = ''
                chunk.loc[unresolved_mask, 'ibar_seq_unresolved'] = chunk.loc[
                    unresolved_mask, 'ibar_seq'
                ]

            unique_cells.update(chunk['bc_wells'].unique())
            unique_guides.update(chunk['grna_name'].unique())
            for row in chunk.itertuples(index=False):
                key = tuple(getattr(row, c) for c in group_cols)
                umi = row.umi
                umi_accum[key][umi] += 1
                per_umi = source_accum[key].setdefault(umi, Counter())
                per_umi[source] += 1

    reads_per_umi_dist = Counter()
    total_umis_raw = 0
    for umi_counts in umi_accum.values():
        for read_count in umi_counts.values():
            reads_per_umi_dist[read_count] += 1
            total_umis_raw += 1

    records = []
    total_reads = sum(s['eligible_reads'] for s in source_stats.values())
    total_umis_after_filter = 0
    total_umis_after_collapse = 0
    umis_filtered = 0
    reads_filtered = 0
    provenance_totals = Counter()

    for group_key, raw_union_counts in umi_accum.items():
        union_counts = Counter(
            {
                umi: count
                for umi, count in raw_union_counts.items()
                if count >= min_reads_per_umi
            }
        )
        umis_filtered += len(raw_union_counts) - len(union_counts)
        reads_filtered += sum(
            count
            for umi, count in raw_union_counts.items()
            if umi not in union_counts
        )
        if not union_counts:
            continue

        total_umis_after_filter += len(union_counts)
        union_map = _directional_umi_cluster_map(
            dict(union_counts), max_edits=max_umi_edits
        )
        union_roots = set(union_map.values())
        total_umis_after_collapse += len(union_roots)

        cluster_sources = defaultdict(set)
        for umi, root in union_map.items():
            cluster_sources[root].update(source_accum[group_key][umi].keys())

        cluster_provenance = Counter()
        for sources in cluster_sources.values():
            if len(sources) > 1:
                cluster_provenance['both'] += 1
            else:
                cluster_provenance[next(iter(sources))] += 1
        provenance_totals.update(cluster_provenance)

        row = dict(zip(group_cols, group_key))
        row['umi_count'] = len(union_roots)
        row['umi_count_raw'] = len(union_counts)
        row['read_count'] = sum(union_counts.values())
        row['sublibrary'] = lib_name
        row['sample_name'] = sample_name

        sources_present = set()
        for source in source_names:
            source_raw = Counter(
                {
                    umi: per_source[source]
                    for umi, per_source in source_accum[group_key].items()
                    if per_source.get(source, 0) >= min_reads_per_umi
                }
            )
            source_map = _directional_umi_cluster_map(
                dict(source_raw), max_edits=max_umi_edits
            )
            row[f'exact_umi_count_{source}'] = len(source_raw)
            row[f'umi_count_{source}'] = len(set(source_map.values()))
            row[f'read_count_{source}'] = sum(
                source_accum[group_key][umi].get(source, 0) for umi in union_counts
            )
            if row[f'read_count_{source}'] > 0:
                sources_present.add(source)

        row['exact_umi_count_overlap'] = sum(
            1
            for umi in union_counts
            if sum(
                source_accum[group_key][umi].get(source, 0) > 0
                for source in source_names
            )
            > 1
        )
        row['umi_count_both_sources'] = cluster_provenance.get('both', 0)
        for source in source_names:
            row[f'umi_count_{source}_only'] = cluster_provenance.get(source, 0)
        row['source'] = (
            'both' if len(sources_present) > 1 else next(iter(sources_present))
        )
        records.append(row)

    deduped = pd.DataFrame(records)
    lib_stats = {
        'total_reads': int(total_reads),
        'unique_cells': len(unique_cells),
        'unique_guides': len(unique_guides),
        'groups': len(deduped),
        'elapsed_s': time.time() - t0,
        'total_umis_raw': total_umis_raw,
        'total_umis_after_filter': total_umis_after_filter,
        'total_umis_after_collapse': total_umis_after_collapse,
        'umis_filtered_low_reads': umis_filtered,
        'reads_filtered_low_reads': reads_filtered,
        'reads_per_umi_dist': dict(reads_per_umi_dist),
        'source_stats': {
            source: {key: int(value) for key, value in stats.items()}
            for source, stats in source_stats.items()
        },
        'cluster_provenance': dict(provenance_totals),
    }
    lib_stats['sequencing_saturation'] = round(
        1.0 - total_umis_after_collapse / total_reads if total_reads else 0.0, 4
    )
    lib_stats['umis_collapsed'] = int(
        total_umis_after_filter - total_umis_after_collapse
    )
    print(
        f"    {lib_name}: {total_reads:,} eligible reads from "
        f"{', '.join(source_names)} -> {total_umis_after_collapse:,} UMI clusters"
    )
    return deduped, lib_stats


def _hamming_dist(a, b):
    """Fast Hamming distance for equal-length strings."""
    return sum(c1 != c2 for c1, c2 in zip(a, b))


def _directional_umi_collapse(umi_counts, max_edits=1):
    """Directional UMI collapsing (same algorithm as umi_tools).

    Args:
        umi_counts: dict of {umi_string: read_count}
        max_edits: max Hamming distance for collapsing

    Returns:
        Number of unique UMI clusters after collapsing.

    Algorithm (UMI-tools directional method):
        1. Sort UMIs by count descending
        2. Build adjacency graph: edge if Hamming dist ≤ max_edits
        3. For each edge (high→low): collapse if count_high >= 2*count_low - 1
           (i.e., the child's count is consistent with being a PCR error of the parent)
        4. Count remaining unmerged clusters
    """
    cluster_map = _directional_umi_cluster_map(umi_counts, max_edits=max_edits)
    return len(set(cluster_map.values()))


def _directional_umi_cluster_map(umi_counts, max_edits=1):
    """Return ``{observed_umi: cluster_representative}``.

    Retaining membership uses the historical directional-collapse traversal
    while allowing the optional dedicated+GEX merge to attach source
    provenance to each final molecule.
    """
    if len(umi_counts) <= 1 or max_edits <= 0:
        return {umi: umi for umi in umi_counts}

    # Sort by count descending
    umis = sorted(umi_counts.keys(), key=lambda u: umi_counts[u], reverse=True)
    counts = {u: umi_counts[u] for u in umis}

    # Track which UMIs are collapsed (merged into a parent)
    parent = {}  # umi → parent umi

    for i, umi_high in enumerate(umis):
        if umi_high in parent:
            continue  # already collapsed
        for j in range(i + 1, len(umis)):
            umi_low = umis[j]
            if umi_low in parent:
                continue  # already collapsed
            if _hamming_dist(umi_high, umi_low) <= max_edits:
                # Directional: collapse only if the parent count is ≥ 2× the child count - 1.
                # This ensures the child's read count is consistent with being a
                # PCR error offspring of the parent (introduced partway through amplification).
                if counts[umi_high] >= 2 * counts[umi_low] - 1:
                    parent[umi_low] = umi_high

    return {umi: parent.get(umi, umi) for umi in umis}


def assign_single_guides(deduped_df, min_umi_threshold=5,
                           multiplet_strategy='all',
                           min_dominance_fraction=0.0):
    """Assign individual guides to cells (Single Guide Mode).

    Args:
        deduped_df: UMI-deduplicated DataFrame
        min_umi_threshold: Minimum UMI count to consider a guide present
        multiplet_strategy: 'all', 'top', or 'discard'
        min_dominance_fraction: For 'top' strategy, minimum fraction of total UMIs

    Returns:
        DataFrame with per-cell guide assignments formatted for single-guide mode.
    """
    # Filter by UMI threshold
    above_thresh = deduped_df[deduped_df['umi_count'] >= min_umi_threshold].copy()

    print(f"\n  Single-guide assignment (threshold={min_umi_threshold}):")
    print(f"    Entries above threshold: {len(above_thresh):,} "
          f"(of {len(deduped_df):,})")
    print(f"    Cells with ≥1 guide:    {above_thresh['bc_wells'].nunique():,}")

    # Count guides per cell
    guides_per_cell = above_thresh.groupby('bc_wells')['grna_name'].nunique()
    singles = (guides_per_cell == 1).sum()
    multiplets = (guides_per_cell > 1).sum()

    print(f"    Single-guide cells:     {singles:,}")
    print(f"    Multi-guide cells:      {multiplets:,}")

    # Apply multiplet strategy
    if multiplet_strategy == 'all':
        # Keep all assignments
        result = above_thresh.copy()
        result['combo_type'] = 'singlet'
        result.loc[
            result['bc_wells'].isin(guides_per_cell[guides_per_cell > 1].index),
            'combo_type'
        ] = 'multiplet'

    elif multiplet_strategy == 'top':
        # For multi-guide cells, keep only the top guide
        result_parts = []

        for bc, group in above_thresh.groupby('bc_wells'):
            if len(group) == 1:
                row = group.iloc[0].copy()
                row['combo_type'] = 'singlet'
                result_parts.append(row)
            else:
                group_sorted = group.sort_values('umi_count', ascending=False)
                top = group_sorted.iloc[0]
                total_umi = group_sorted['umi_count'].sum()
                dominance = top['umi_count'] / total_umi

                if min_dominance_fraction > 0 and dominance < min_dominance_fraction:
                    continue

                row = top.copy()
                row['combo_type'] = 'multiplet'
                result_parts.append(row)

        result = pd.DataFrame(result_parts)

    elif multiplet_strategy == 'discard':
        # Discard all multi-guide cells
        single_cells = guides_per_cell[guides_per_cell == 1].index
        result = above_thresh[above_thresh['bc_wells'].isin(single_cells)].copy()
        result['combo_type'] = 'singlet'
        print(f"    Discarded {multiplets:,} multi-guide cells")

    else:
        raise ValueError(f"Unknown multiplet_strategy: {multiplet_strategy}")

    if len(result) > 0:
        result['n_guides'] = result['bc_wells'].map(guides_per_cell)

        print("    Combo types (per cell):")
        cell_level = result.drop_duplicates('bc_wells')
        for status, count in cell_level['combo_type'].value_counts().items():
            print(f"      {status:>15}: {count:,}")
            
    print(f"    Final assignments:      {len(result):,}")
    return result


def build_guide_umi_matrix(deduped_df, min_umi_threshold=0):
    """Build a cells × guides UMI count matrix (sparse).

    Args:
        deduped_df: UMI-deduplicated DataFrame
        min_umi_threshold: Only include cells with at least this many total UMIs

    Returns:
        sparse matrix, cell_barcodes list, guide_names list
    """
    # 1. First aggregate UMI counts for each cell-guide pair in long format
    # This prevents duplicate pairs if present and acts as a safe base.
    grouped_df = deduped_df.groupby(['bc_wells', 'grna_name'], as_index=False)['umi_count'].sum()

    # 2. Apply UMI threshold directly to cell sums
    if min_umi_threshold > 0:
        cell_totals = grouped_df.groupby('bc_wells')['umi_count'].sum()
        valid_cells = cell_totals[cell_totals >= min_umi_threshold].index
        grouped_df = grouped_df[grouped_df['bc_wells'].isin(valid_cells)]

    # 3. Create categorical variables for mapping barcodes/guides to integer indices
    bc_cat = pd.Categorical(grouped_df['bc_wells'])
    grna_cat = pd.Categorical(grouped_df['grna_name'])

    cell_barcodes = list(bc_cat.categories)
    guide_names = list(grna_cat.categories)
    
    data_vals = grouped_df['umi_count'].values
    row_ind = bc_cat.codes
    col_ind = grna_cat.codes

    # 4. Create the sparse matrix efficiently using COO format, then convert to CSR
    matrix = sparse.coo_matrix(
        (data_vals, (row_ind, col_ind)), 
        shape=(len(cell_barcodes), len(guide_names))
    ).tocsr()

    print(f"\n  Guide UMI matrix: {matrix.shape[0]:,} cells × {matrix.shape[1]:,} guides")
    if matrix.shape[0] > 0 and matrix.shape[1] > 0:
        print(f"    Non-zero entries: {matrix.nnz:,}")
        print(f"    Sparsity: {1 - matrix.nnz / (matrix.shape[0] * matrix.shape[1]):.4f}")
    else:
        print("    Empty matrix.")

    return matrix, cell_barcodes, guide_names


def assign_dual_guide_pairs(deduped_df, ref_df, min_umi_threshold=5, ibar_col='ibar_seq_corrected'):
    """Assign dual-guide pairs to cells, supporting iBAR-based inference for half-detected pairs."""
    if 'guide_unit_id' in deduped_df.columns:
        return assign_cells(deduped_df, ref_df, min_umi_threshold)
    print(f"\n  Dual-guide pair assignment (pair UMI threshold={min_umi_threshold}):")
    
    # 1. Pre-build lookup for inference
    pair_to_guides = {}
    if ref_df is not None:
        for pair_id, group in ref_df.groupby('pair_id'):
            pos1_row = group[group['position'] == 1]
            pos2_row = group[group['position'] == 2]
            
            pair_to_guides[pair_id] = {
                1: pos1_row.iloc[0].to_dict() if len(pos1_row) else {},
                2: pos2_row.iloc[0].to_dict() if len(pos2_row) else {}
            }
            
    cell_records = []
    
    # Sort for deterministic processing
    if 'pair_id' not in deduped_df.columns:
        print("    WARNING: 'pair_id' not found in deduped data. Ensure reference merge was successful.")
        return pd.DataFrame()
        
    sorted_df = deduped_df.dropna(subset=['pair_id']).sort_values(
        ['bc_wells', 'pair_id', 'umi_count'], ascending=[True, True, False]
    )
    
    n_processed = 0
    unique_cells = sorted_df['bc_wells'].unique()
    
    for bc, cell_df in sorted_df.groupby('bc_wells'):
        pairs_for_cell = []
        
        for pair_id, pair_df in cell_df.groupby('pair_id'):
            pos1_obs = pair_df[pair_df['position'] == 1]
            pos2_obs = pair_df[pair_df['position'] == 2]
            
            if len(pos1_obs) > 0 and len(pos2_obs) > 0:
                # Complete pair - take combinations
                for _, r1 in pos1_obs.iterrows():
                    for _, r2 in pos2_obs.iterrows():
                        combo_umi = r1['umi_count'] + r2['umi_count']
                        if combo_umi >= min_umi_threshold:
                            pairs_for_cell.append({
                                'bc_wells': bc,
                                'pair_id': pair_id,
                                'target_combo': r1.get('target_combo', ''),
                                'gi_pair': r1.get('gi_pair', ''),
                                'detection_status': 'complete',
                                
                                'pos1_grna': r1['grna_name'],
                                'pos1_ibar': r1.get(ibar_col, ''),
                                'pos1_target': r1.get('target_gene', ''),
                                'pos1_umi': r1['umi_count'],
                                'pos1_umi_raw': r1.get('umi_count_raw', r1['umi_count']),
                                
                                'pos2_grna': r2['grna_name'],
                                'pos2_ibar': r2.get(ibar_col, ''),
                                'pos2_target': r2.get('target_gene', ''),
                                'pos2_umi': r2['umi_count'],
                                'pos2_umi_raw': r2.get('umi_count_raw', r2['umi_count']),
                                
                                'inferred_pos1_grna': '',
                                'inferred_pos1_target': '',
                                'inferred_pos2_grna': '',
                                'inferred_pos2_target': '',
                                'inference_method': '',
                                
                                'pair_umi_total': combo_umi
                            })
                            
            elif len(pos1_obs) > 0:
                # Half pos1
                for _, r1 in pos1_obs.iterrows():
                    if r1['umi_count'] >= min_umi_threshold:
                        expected_pos2 = pair_to_guides.get(pair_id, {}).get(2, {})
                        pairs_for_cell.append({
                            'bc_wells': bc,
                            'pair_id': pair_id,
                            'target_combo': r1.get('target_combo', ''),
                            'gi_pair': r1.get('gi_pair', ''),
                            'detection_status': 'inferred',
                            
                            'pos1_grna': r1['grna_name'],
                            'pos1_ibar': r1.get(ibar_col, ''),
                            'pos1_target': r1.get('target_gene', ''),
                            'pos1_umi': r1['umi_count'],
                            'pos1_umi_raw': r1.get('umi_count_raw', r1['umi_count']),
                            
                            'pos2_grna': '',
                            'pos2_ibar': '',
                            'pos2_target': '',
                            'pos2_umi': 0,
                            'pos2_umi_raw': 0,
                            
                            'inferred_pos1_grna': '',
                            'inferred_pos1_target': '',
                            'inferred_pos2_grna': expected_pos2.get('grna_name', ''),
                            'inferred_pos2_target': expected_pos2.get('target_gene', ''),
                            'inference_method': 'known_pair_id',
                            
                            'pair_umi_total': r1['umi_count']
                        })
                        
            elif len(pos2_obs) > 0:
                # Half pos2
                for _, r2 in pos2_obs.iterrows():
                    if r2['umi_count'] >= min_umi_threshold:
                        expected_pos1 = pair_to_guides.get(pair_id, {}).get(1, {})
                        pairs_for_cell.append({
                            'bc_wells': bc,
                            'pair_id': pair_id,
                            'target_combo': r2.get('target_combo', ''),
                            'gi_pair': r2.get('gi_pair', ''),
                            'detection_status': 'inferred',
                            
                            'pos1_grna': '',
                            'pos1_ibar': '',
                            'pos1_target': '',
                            'pos1_umi': 0,
                            'pos1_umi_raw': 0,
                            
                            'pos2_grna': r2['grna_name'],
                            'pos2_ibar': r2.get(ibar_col, ''),
                            'pos2_target': r2.get('target_gene', ''),
                            'pos2_umi': r2['umi_count'],
                            'pos2_umi_raw': r2.get('umi_count_raw', r2['umi_count']),
                            
                            'inferred_pos1_grna': expected_pos1.get('grna_name', ''),
                            'inferred_pos1_target': expected_pos1.get('target_gene', ''),
                            'inferred_pos2_grna': '',
                            'inferred_pos2_target': '',
                            'inference_method': 'known_pair_id',
                            
                            'pair_umi_total': r2['umi_count']
                        })
                        
        n_combos = len(pairs_for_cell)
        if n_combos > 0:
            combo_type = 'singlet' if n_combos == 1 else 'multiplet'
            for p in pairs_for_cell:
                p['combo_type'] = combo_type
                p['n_combos'] = n_combos
            cell_records.extend(pairs_for_cell)
            
        n_processed += 1
        if n_processed % 100_000 == 0:
            print(f"    Processed {n_processed:,}/{len(unique_cells):,} cells...")

    result_df = pd.DataFrame(cell_records)
    n_assigned_cells = result_df['bc_wells'].nunique()
    print(f"    Total cells with pair assignment: {n_assigned_cells:,}")
    print(f"    Total pair-assignments (rows):    {len(result_df):,}")
    
    if len(result_df) > 0:
        # Cell-level counts (deduplicated by bc_wells)
        cell_level = result_df.drop_duplicates('bc_wells')
        print("    Combo types (per cell):")
        for status, count in cell_level['combo_type'].value_counts().items():
            print(f"      {status:>15}: {count:,}")
        print("    Detection status (per pair-assignment):")
        for status, count in result_df['detection_status'].value_counts().items():
            print(f"      {status:>15}: {count:,}")
            
    return result_df


def compute_qc_stats(lib_stats_all, deduped_df, assignment_df, guide_matrix,
                     cell_barcodes, guide_names, dual_guide_counting_mode=False):
    """Compute QC statistics.

    Args:
        lib_stats_all: dict of {lib_name: {total_reads, unique_cells, ...}}
    """
    stats = {}
    design_aware = 'pair_resolution_basis' in assignment_df.columns
    if design_aware:
        stats['guide_observation_records'] = len(assignment_df)
        stats['unresolved_observation_records'] = int(assignment_df['pair_id'].fillna('').eq('').sum())
        assignment_df = eligible_assignments(assignment_df)

    # Per-library stats from streaming summaries
    total_reads = 0
    total_unique_cells = 0
    all_guides = set()
    for lib_name, ls in lib_stats_all.items():
        stats[f'lib_{lib_name}_total_reads'] = ls['total_reads']
        stats[f'lib_{lib_name}_unique_cells'] = ls['unique_cells']
        total_reads += ls['total_reads']
        total_unique_cells += ls['unique_cells']

    # Overall stats
    stats['total_matched_reads'] = total_reads
    stats['total_unique_cells_with_reads'] = deduped_df['bc_wells'].nunique()
    stats['total_unique_guides_detected'] = deduped_df['grna_name'].nunique()

    # After dedup
    stats['unique_cell_guide_pairs'] = len(deduped_df)
    stats['median_umi_per_cell_guide'] = float(deduped_df['umi_count'].median())
    stats['mean_umi_per_cell_guide'] = float(deduped_df['umi_count'].mean())

    # UMI filtering and saturation stats
    total_umis_raw = sum(ls.get('total_umis_raw', 0) for ls in lib_stats_all.values())
    total_umis_after_filter = sum(ls.get('total_umis_after_filter', 0) for ls in lib_stats_all.values())
    total_umis_after_collapse = sum(ls.get('total_umis_after_collapse', 0) for ls in lib_stats_all.values())
    umis_filtered = sum(ls.get('umis_filtered_low_reads', 0) for ls in lib_stats_all.values())
    stats['total_umis_raw'] = int(total_umis_raw)
    stats['total_umis_after_filter'] = int(total_umis_after_filter)
    stats['total_umis_after_collapse'] = int(total_umis_after_collapse)
    stats['umis_filtered_low_reads'] = int(umis_filtered)
    stats['sequencing_saturation'] = round(
        1.0 - total_umis_after_collapse / total_reads if total_reads > 0 else 0.0, 4
    )
    if 'read_count' in deduped_df.columns and len(deduped_df) > 0:
        rpu = deduped_df['read_count'] / deduped_df['umi_count'].clip(lower=1)
        stats['median_reads_per_umi'] = round(float(rpu.median()), 2)
        stats['mean_reads_per_umi'] = round(float(rpu.mean()), 2)

    # After assignment
    stats['cells_with_assignment'] = assignment_df['bc_wells'].nunique()
    stats['assignments_total'] = len(assignment_df)

    if dual_guide_counting_mode:
        # --- Cell-level counts (each cell counted once) ---
        cell_level = assignment_df.drop_duplicates('bc_wells')
        n_singlet = int((cell_level['combo_type'] == 'singlet').sum())
        n_multiplet = int((cell_level['combo_type'] == 'multiplet').sum())
        stats['n_singlet_cells'] = n_singlet
        stats['n_multiplet_cells'] = n_multiplet

        # Classify multiplet cells by detection status composition
        if n_multiplet > 0:
            mult_cells = cell_level[cell_level['combo_type'] == 'multiplet']['bc_wells']
            mult_assignments = assignment_df[assignment_df['bc_wells'].isin(mult_cells)]
            # A multiplet cell has inferred pairs if ANY of its pair-assignments are 'inferred'
            cells_with_inferred = mult_assignments[
                mult_assignments['detection_status'] == 'inferred'
            ]['bc_wells'].unique()
            stats['n_multiplet_with_inferred_cells'] = len(cells_with_inferred)
            stats['n_multiplet_all_detected_cells'] = n_multiplet - len(cells_with_inferred)

        # --- Assignment-level counts (rows = individual pair-assignments) ---
        stats['n_pair_assignments_total'] = len(assignment_df)
        for dtype in assignment_df['detection_status'].unique():
            count = int((assignment_df['detection_status'] == dtype).sum())
            stats[f'pair_assignments_{dtype}'] = count

        pairs_per_cell = assignment_df.groupby('bc_wells')['pair_id'].nunique()
        if design_aware:
            pairs_per_cell = pairs_per_cell.reindex(deduped_df['bc_wells'].unique(), fill_value=0)
        stats['mean_pairs_per_cell'] = float(pairs_per_cell.mean()) if len(pairs_per_cell) else 0.0
        stats['median_pairs_per_cell'] = float(pairs_per_cell.median()) if len(pairs_per_cell) else 0.0
        stats['max_pairs_per_cell'] = int(pairs_per_cell.max()) if len(pairs_per_cell) else 0

        # Target combo counts (assignment-level is correct here)
        for tc in assignment_df['target_combo'].unique():
            count = int((assignment_df['target_combo'] == tc).sum())
            stats[f'target_combo_{tc}_assignments'] = count
    else:
        for atype in assignment_df['combo_type'].unique():
            count = int((assignment_df['combo_type'] == atype).sum())
            stats[f'combo_type_{atype}'] = count

        guides_per_cell = assignment_df.groupby('bc_wells')['grna_name'].nunique()
        stats['mean_guides_per_cell'] = float(guides_per_cell.mean())
        stats['median_guides_per_cell'] = float(guides_per_cell.median())
        stats['max_guides_per_cell'] = int(guides_per_cell.max())

        for gene in assignment_df['target_gene'].unique():
            count = int((assignment_df['target_gene'] == gene).sum())
            stats[f'target_{gene}_assignments'] = count

    return stats


def plot_umi_threshold_sensitivity(deduped_df, save_path=None):
    """Plot number of cells with valid gRNA at different UMI thresholds.

    Uses the sublibrary column from deduped_df (no merged_df needed).
    """
    thresholds = list(range(1, 21)) + [25, 30, 40, 50]

    sublibraries = sorted(deduped_df['sublibrary'].unique())

    fig, ax = plt.subplots(figsize=(10, 6))

    # For each threshold, count cells with at least one gRNA above threshold
    for lib in sublibraries:
        lib_df = deduped_df[deduped_df['sublibrary'] == lib]
        counts = []
        for t in thresholds:
            n = lib_df[lib_df['umi_count'] >= t]['bc_wells'].nunique()
            counts.append(n)
        ax.plot(thresholds, counts, marker='o', markersize=3, label=lib, linewidth=1.5)

    # Total line
    total_counts = []
    for t in thresholds:
        n = deduped_df[deduped_df['umi_count'] >= t]['bc_wells'].nunique()
        total_counts.append(n)
    ax.plot(thresholds, total_counts, marker='s', markersize=4, label='Total',
            linewidth=2.5, color='black', linestyle='--')

    ax.set_xlabel('UMI Threshold', fontsize=12)
    ax.set_ylabel('Number of Cells with Valid gRNA', fontsize=12)
    ax.set_title('UMI Threshold Sensitivity: Cells with Valid gRNA', fontsize=14)
    ax.legend(fontsize=10, loc='upper right')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    ax.grid(True, which='minor', alpha=0.1)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"\n  Saved UMI threshold sensitivity plot: {save_path}")
    plt.close(fig)


def _compute_saturation_curve(lib_stats_all, min_reads_per_umi=1,
                               n_points=10, n_reps=3):
    """Compute sequencing saturation downsampling curve from accumulated data.

    For each library and combined, estimate how many unique UMIs would be
    detected at various fractions of the actual sequencing depth. Uses
    Poisson approximation for speed.

    Args:
        lib_stats_all: dict of {lib_name: lib_stats} with reads_per_umi_dist
        min_reads_per_umi: UMI validity threshold
        n_points: Number of downsampling fractions
        n_reps: Number of repetitions for averaging

    Returns:
        dict of {lib_name: list_of_dicts, 'combined': list_of_dicts}
        Each dict: {fraction, reads_sampled, umis_detected, saturation}
    """
    rng = np.random.default_rng(42)
    fractions = np.linspace(0.1, 1.0, n_points)
    curve_data = {}

    # Flatten all libraries' reads-per-UMI distributions into (read_count, n_umis) arrays
    for lib_name, stats in lib_stats_all.items():
        dist = stats.get('reads_per_umi_dist', {})
        if not dist:
            continue

        total_reads = stats['total_reads']
        # dist maps read_count -> n_umis_with_that_count
        read_counts = np.array(list(dist.keys()), dtype=np.int64)
        n_umis = np.array(list(dist.values()), dtype=np.int64)

        lib_curve = []
        for frac in fractions:
            umis_detected_reps = []
            for _ in range(n_reps):
                # For each UMI with k reads, sample reads ~ Poisson(k * frac)
                sampled_reads = rng.poisson(read_counts * frac)
                # A UMI is detected if sampled_reads >= min_reads_per_umi
                detected = np.sum(n_umis[sampled_reads >= min_reads_per_umi])
                umis_detected_reps.append(int(detected))

            avg_umis = int(np.mean(umis_detected_reps))
            reads_sampled = int(total_reads * frac)
            sat = 1.0 - avg_umis / reads_sampled if reads_sampled > 0 else 0.0

            lib_curve.append({
                'fraction': round(float(frac), 2),
                'reads_sampled': reads_sampled,
                'umis_detected': avg_umis,
                'saturation': round(sat, 4),
            })
        curve_data[lib_name] = lib_curve

    # Combined curve: merge all libraries' distributions
    combined_dist = Counter()
    combined_total_reads = 0
    for stats in lib_stats_all.values():
        for rc_str, count in stats.get('reads_per_umi_dist', {}).items():
            combined_dist[int(rc_str)] += count
        combined_total_reads += stats['total_reads']

    if combined_dist:
        read_counts = np.array(list(combined_dist.keys()), dtype=np.int64)
        n_umis = np.array(list(combined_dist.values()), dtype=np.int64)

        combined_curve = []
        for frac in fractions:
            umis_detected_reps = []
            for _ in range(n_reps):
                sampled_reads = rng.poisson(read_counts * frac)
                detected = np.sum(n_umis[sampled_reads >= min_reads_per_umi])
                umis_detected_reps.append(int(detected))

            avg_umis = int(np.mean(umis_detected_reps))
            reads_sampled = int(combined_total_reads * frac)
            sat = 1.0 - avg_umis / reads_sampled if reads_sampled > 0 else 0.0

            combined_curve.append({
                'fraction': round(float(frac), 2),
                'reads_sampled': reads_sampled,
                'umis_detected': avg_umis,
                'saturation': round(sat, 4),
            })
        curve_data['combined'] = combined_curve

    return curve_data


def plot_reads_per_umi_distribution(lib_stats_all, min_reads_per_umi=1,
                                     save_path=None):
    """Plot reads-per-UMI distribution: histogram and reverse CDF.

    Left panel: Log-scale histogram with threshold line.
    Right panel: Reverse CDF (fraction of UMIs with >= x reads).
    """
    # Merge all libraries' distributions
    combined_dist = Counter()
    per_lib_dists = {}
    for lib_name, stats in lib_stats_all.items():
        dist = stats.get('reads_per_umi_dist', {})
        lib_dist = Counter()
        for rc, count in dist.items():
            combined_dist[int(rc)] += count
            lib_dist[int(rc)] += count
        per_lib_dists[lib_name] = lib_dist

    if not combined_dist:
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # --- Left panel: Histogram ---
    max_rc = max(combined_dist.keys())
    bins = np.logspace(0, np.log10(max(max_rc, 2)), min(50, max_rc))

    # Expand distribution into array of read counts for histogram
    all_values = []
    for rc, count in combined_dist.items():
        all_values.extend([rc] * count)
    all_values = np.array(all_values)

    ax1.hist(all_values, bins=bins, color='#5B8DEF', edgecolor='white',
             linewidth=0.2, alpha=0.8, label='All UMIs')

    if min_reads_per_umi > 1:
        ax1.axvline(min_reads_per_umi, color='#E74C3C', linewidth=2,
                    linestyle='--', label=f'Threshold = {min_reads_per_umi}')
        n_below = sum(c for rc, c in combined_dist.items() if rc < min_reads_per_umi)
        n_total = sum(combined_dist.values())
        pct_filtered = n_below / n_total * 100 if n_total > 0 else 0
        ax1.fill_betweenx([0, ax1.get_ylim()[1] * 10],
                          0.5, min_reads_per_umi,
                          color='#E74C3C', alpha=0.1)
        ax1.text(0.02, 0.95, f'{pct_filtered:.1f}% filtered',
                 transform=ax1.transAxes, fontsize=10, color='#E74C3C',
                 verticalalignment='top')

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('Reads per UMI', fontsize=12)
    ax1.set_ylabel('Number of UMIs', fontsize=12)
    ax1.set_title('Reads-per-UMI Distribution', fontsize=13)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)

    # --- Right panel: Reverse CDF per library ---
    for lib_name in sorted(per_lib_dists.keys()):
        lib_dist = per_lib_dists[lib_name]
        sorted_rcs = sorted(lib_dist.keys())
        total = sum(lib_dist.values())
        cumsum = 0
        x_vals, y_vals = [], []
        for rc in sorted_rcs:
            x_vals.append(rc)
            y_vals.append((total - cumsum) / total * 100)
            cumsum += lib_dist[rc]
        ax2.plot(x_vals, y_vals, marker='.', markersize=2, linewidth=1.5,
                 label=lib_name)

    # Combined line
    sorted_rcs = sorted(combined_dist.keys())
    total = sum(combined_dist.values())
    cumsum = 0
    x_vals, y_vals = [], []
    for rc in sorted_rcs:
        x_vals.append(rc)
        y_vals.append((total - cumsum) / total * 100)
        cumsum += combined_dist[rc]
    ax2.plot(x_vals, y_vals, marker='s', markersize=3, linewidth=2.5,
             color='black', linestyle='--', label='Combined')

    if min_reads_per_umi > 1:
        ax2.axvline(min_reads_per_umi, color='#E74C3C', linewidth=2,
                    linestyle='--', alpha=0.7)

    ax2.set_xscale('log')
    ax2.set_xlabel('Min reads per UMI (x)', fontsize=12)
    ax2.set_ylabel('% UMIs with \u2265 x reads', fontsize=12)
    ax2.set_title('UMI Survival (Reverse CDF)', fontsize=13)
    ax2.legend(fontsize=9, loc='upper right')
    ax2.set_ylim(0, 105)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"  Saved reads-per-UMI distribution plot: {save_path}")
    plt.close(fig)


def plot_saturation_curve(curve_data, save_path=None):
    """Plot sequencing saturation downsampling curve.

    Left panel: UMIs detected vs reads sampled (diminishing returns).
    Right panel: Saturation % vs reads sampled.
    """
    if not curve_data:
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    colors = plt.cm.Set2(np.linspace(0, 1, max(len(curve_data), 3)))
    lib_names = [k for k in sorted(curve_data.keys()) if k != 'combined']

    for i, lib_name in enumerate(lib_names):
        points = curve_data[lib_name]
        reads = [p['reads_sampled'] for p in points]
        umis = [p['umis_detected'] for p in points]
        sats = [p['saturation'] * 100 for p in points]

        ax1.plot(reads, umis, marker='o', markersize=4, linewidth=1.5,
                 color=colors[i], label=lib_name)
        ax2.plot(reads, sats, marker='o', markersize=4, linewidth=1.5,
                 color=colors[i], label=lib_name)

    # Combined line
    if 'combined' in curve_data:
        points = curve_data['combined']
        reads = [p['reads_sampled'] for p in points]
        umis = [p['umis_detected'] for p in points]
        sats = [p['saturation'] * 100 for p in points]

        ax1.plot(reads, umis, marker='s', markersize=5, linewidth=2.5,
                 color='black', linestyle='--', label='Combined')
        ax2.plot(reads, sats, marker='s', markersize=5, linewidth=2.5,
                 color='black', linestyle='--', label='Combined')

        # Mark current depth
        ax1.plot(reads[-1], umis[-1], '*', markersize=15, color='#E74C3C',
                 zorder=5, label='Current depth')
        ax2.plot(reads[-1], sats[-1], '*', markersize=15, color='#E74C3C',
                 zorder=5)

    # Left panel formatting
    max_reads = max(p['reads_sampled'] for pts in curve_data.values() for p in pts)
    ax1.plot([0, max_reads], [0, max_reads], 'k:', alpha=0.2, label='No saturation')
    ax1.set_xlabel('Reads sampled', fontsize=12)
    ax1.set_ylabel('Unique UMIs detected', fontsize=12)
    ax1.set_title('Sequencing Depth vs UMI Discovery', fontsize=13)
    ax1.legend(fontsize=8, loc='lower right')
    ax1.grid(True, alpha=0.3)
    ax1.ticklabel_format(style='sci', axis='both', scilimits=(0, 0))

    # Right panel formatting
    for pct, ls in [(50, ':'), (75, '--'), (90, '-.')]:
        ax2.axhline(pct, color='gray', linewidth=0.8, linestyle=ls, alpha=0.5)
        ax2.text(0.02, pct / 100 - 0.01, f'{pct}%', transform=ax2.get_yaxis_transform(),
                 fontsize=8, color='gray')
    ax2.set_xlabel('Reads sampled', fontsize=12)
    ax2.set_ylabel('Sequencing Saturation (%)', fontsize=12)
    ax2.set_title('Sequencing Saturation Curve', fontsize=13)
    ax2.legend(fontsize=8, loc='lower right')
    ax2.set_ylim(0, 105)
    ax2.grid(True, alpha=0.3)
    ax2.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"  Saved sequencing saturation curve: {save_path}")
    plt.close(fig)


def plot_median_reads_per_umi_per_cell(deduped_df, save_path=None):
    """Plot distribution of median reads-per-UMI across cells, split by sublibrary.

    Each cell's median reads-per-UMI is computed as median(read_count / umi_count)
    across that cell's groups. Higher values indicate deeper per-cell sequencing.
    """
    if len(deduped_df) == 0 or 'read_count' not in deduped_df.columns:
        return

    # Compute reads-per-UMI ratio per group, then median per cell
    df = deduped_df.copy()
    df['reads_per_umi'] = df['read_count'] / df['umi_count'].clip(lower=1)

    cell_medians = df.groupby(['bc_wells', 'sublibrary'])['reads_per_umi'].median().reset_index()

    sublibraries = sorted(cell_medians['sublibrary'].unique())
    n_libs = len(sublibraries)

    fig, ax = plt.subplots(figsize=(max(6, n_libs * 1.5 + 2), 5))

    # Violin plot
    data_by_lib = [cell_medians[cell_medians['sublibrary'] == lib]['reads_per_umi'].values
                   for lib in sublibraries]

    parts = ax.violinplot(data_by_lib, positions=range(len(sublibraries)),
                          showmedians=True, showextrema=False)

    for pc in parts['bodies']:
        pc.set_facecolor('#5B8DEF')
        pc.set_alpha(0.6)
    parts['cmedians'].set_color('#E74C3C')
    parts['cmedians'].set_linewidth(2)

    # Add quartile markers
    for i, lib_data in enumerate(data_by_lib):
        q25, q50, q75 = np.percentile(lib_data, [25, 50, 75])
        ax.vlines(i, q25, q75, color='black', linewidth=3, alpha=0.5)
        ax.text(i, q50, f' {q50:.1f}', fontsize=8, va='center', ha='left')

    ax.set_xticks(range(len(sublibraries)))
    ax.set_xticklabels(sublibraries, fontsize=11)
    ax.set_xlabel('Sublibrary', fontsize=12)
    ax.set_ylabel('Median reads per UMI (per cell)', fontsize=12)
    ax.set_title('Per-Cell Sequencing Depth', fontsize=13)
    ax.grid(True, axis='y', alpha=0.3)

    # Overall median line
    overall_median = cell_medians['reads_per_umi'].median()
    ax.axhline(overall_median, color='#E74C3C', linewidth=1.5, linestyle='--',
               alpha=0.7, label=f'Overall median = {overall_median:.1f}')
    ax.legend(fontsize=9)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"  Saved per-cell reads-per-UMI plot: {save_path}")
    plt.close(fig)


def _determine_group_cols(match_path, dedup_ibar):
    """Determine grouping columns by peeking at the header of a match file."""
    header_cols = pd.read_csv(match_path, sep='\t', nrows=0).columns.tolist()
    if 'guide_unit_id' in header_cols:
        columns = ['bc_wells', 'grna_name', 'target_gene', 'guide_unit_id', 'position',
                   'ibar_seq_corrected', 'ibar_seq_unresolved', 'ibar_status',
                   'pair_candidates', 'pair_candidate_count', 'structure_status',
                   'guide_design_fingerprint', 'resolution_status', 'legacy_model_eligible']
        return columns
    if dedup_ibar:
        if 'ibar_seq_corrected' in header_cols:
            return ['bc_wells', 'grna_name', 'target_gene', 'ibar_seq_corrected', 'ibar_seq_unresolved']
        elif 'ibar_seq' in header_cols:
            return ['bc_wells', 'grna_name', 'target_gene', 'ibar_seq']
        else:
            print("  WARNING: dedup_ibar is True, but no iBAR columns found. Disabling dedup_ibar.")
            return ['bc_wells', 'grna_name', 'target_gene']
    else:
        return ['bc_wells', 'grna_name', 'target_gene']


def _process_one_library(task, group_cols, max_umi_edits, min_reads_per_umi):
    """Module-level wrapper for multiprocessing (must be picklable)."""
    match_path, lib_name, sample_name, sublib_idx = task
    return stream_library_umi_dedup(
        match_path, lib_name, sample_name, sublib_idx,
        group_cols=group_cols,
        max_umi_edits=max_umi_edits,
        min_reads_per_umi=min_reads_per_umi
    )


def _source_mode_tasks(cfg, output_dir, source_mode, gex_namespace):
    """Resolve per-library inputs for one assignment product."""
    tasks = []
    missing = []
    for sublib_idx, lib in enumerate(cfg['libraries']):
        lib_name = lib['sublibrary_name']
        dedicated_path = _resolve_match_path(output_dir, lib_name)
        gex_path = _resolve_gex_match_path(output_dir, gex_namespace, lib_name)
        if source_mode == 'dedicated':
            if dedicated_path is None:
                print(f"  WARNING: Match file not found for {lib_name}")
                continue
            paths = [('dedicated', dedicated_path)]
        elif source_mode == 'gex':
            if gex_path is None:
                missing.append(f"{lib_name}: GEX Step1c match file")
                continue
            paths = [('gex', gex_path)]
        elif source_mode == 'combined':
            absent = []
            if dedicated_path is None:
                absent.append('dedicated')
            if gex_path is None:
                absent.append('GEX Step1c')
            if absent:
                missing.append(f"{lib_name}: {', '.join(absent)} match file")
                continue
            paths = [('dedicated', dedicated_path), ('gex', gex_path)]
        else:
            raise ValueError(f"Unknown source mode: {source_mode}")
        tasks.append((paths, lib_name, lib['sample_name'], sublib_idx))

    if missing:
        raise FileNotFoundError(
            f"Incomplete {source_mode} evidence; refusing a partial matrix:\n  "
            + "\n  ".join(missing)
        )
    if not tasks:
        if source_mode == 'dedicated':
            raise FileNotFoundError("No match files found. Run Step 3 first.")
        raise FileNotFoundError(f"No {source_mode} Step 3 match files found.")
    return tasks


def _annotate_and_assign(deduped_df, ref_df, available_cols, group_cols, ga):
    """Apply the unchanged design annotation and guide-calling logic."""
    dual_mode = ga.get('dual_guide_counting_mode', False)
    if 'guide_unit_id' in deduped_df.columns:
        if ref_df is None:
            raise ValueError('Design-aware assignment requires a reference')
        import json
        deduped_df = deduped_df.copy()
        candidate_lists = deduped_df['pair_candidates'].map(json.loads)
        deduped_df['pair_id'] = candidate_lists.map(lambda candidates: candidates[0] if len(candidates) == 1 else '')
        metadata = ref_df.drop_duplicates('pair_id').set_index('pair_id') if 'pair_id' in ref_df.columns else None
        if metadata is not None:
            for column in ('target_combo', 'target_version', 'gi_pair', 'gi_version'):
                if column in metadata:
                    deduped_df[column] = deduped_df['pair_id'].map(metadata[column]).fillna('')
        if dual_mode:
            assignment = assign_cells(deduped_df, ref_df, ga.get('min_umi_threshold', 5),
                                      infer_missing_partner=ga.get('infer_missing_partner', 'auto') is not False)
        else:
            assignment = assign_single_guides(deduped_df, ga.get('min_umi_threshold', 5),
                                               ga.get('multiplet_strategy', 'all'), ga.get('min_dominance_fraction', 0.0))
        return deduped_df, assignment
    if ref_df is not None:
        if dual_mode:
            ibar_cols = [
                c for c in group_cols
                if 'ibar' in c and c != 'ibar_seq_unresolved'
            ]
            if not ibar_cols:
                raise ValueError(
                    "dual_guide_counting_mode requires resolved iBAR grouping"
                )
            ibar_col = ibar_cols[0]
            ref_annot = ref_df[available_cols + ['ibar_seq']].drop_duplicates()
            deduped_df = deduped_df.merge(
                ref_annot,
                left_on=['grna_name', ibar_col],
                right_on=['grna_name', 'ibar_seq'],
                how='left',
            )
            print(
                "\n  Joined dual-guide design annotations on "
                f"grna_name + {ibar_col}"
            )
        else:
            ref_annot = ref_df[available_cols].drop_duplicates(
                subset=['grna_name']
            )
            deduped_df = deduped_df.merge(
                ref_annot, on='grna_name', how='left'
            )
            print("\n  Joined single-guide design annotations")

    if dual_mode:
        if ref_df is None:
            raise ValueError(
                "dual_guide_counting_mode requires paths.grna_reference"
            )
        ibar_cols = [
            c for c in group_cols
            if 'ibar' in c and c != 'ibar_seq_unresolved'
        ]
        if not ibar_cols:
            raise ValueError(
                "dual_guide_counting_mode requires resolved iBAR evidence"
            )
        assignment_df = assign_dual_guide_pairs(
            deduped_df,
            ref_df=ref_df,
            min_umi_threshold=ga['min_umi_threshold'],
            ibar_col=ibar_cols[0],
        )
    else:
        assignment_df = assign_single_guides(
            deduped_df,
            min_umi_threshold=ga['min_umi_threshold'],
            multiplet_strategy=ga['multiplet_strategy'],
            min_dominance_fraction=ga['min_dominance_fraction'],
        )
    return deduped_df, assignment_df


def _stream_assignment_product(
    cfg,
    output_dir,
    source_mode,
    gex_namespace,
    dedup_ibar,
    max_umi_edits,
    min_reads_per_umi,
    ncpus,
    ref_df,
    available_cols,
):
    """Build one dedicated, GEX-only, or union assignment product."""
    tasks = _source_mode_tasks(cfg, output_dir, source_mode, gex_namespace)
    group_cols = _determine_group_cols(tasks[0][0][0][1], dedup_ibar)
    if not enabled(cfg) and 'guide_unit_id' in group_cols:
        raise ValueError('Design-aware match outputs require guide_design.enabled; use a separate output directory for legacy runs')
    if enabled(cfg):
        expected_contract = compile_design(ref_df, cfg)
        effective_fingerprints = set()
        for paths, lib_name, _, _ in tasks:
            for _, match_path in paths:
                contract_path = Path(match_path).parent / f'{lib_name}_design.yaml'
                if not contract_path.is_file():
                    raise ValueError(f'Missing resolved design: {contract_path}; rerun Step 3')
                with contract_path.open() as handle:
                    saved = yaml.safe_load(handle)
                if saved.get('input_fingerprint', saved['fingerprint']) != expected_contract['fingerprint']:
                    raise ValueError('Reference or guide settings changed since Step 3; rerun matching')
                effective_fingerprints.add(saved['fingerprint'])
        if len(effective_fingerprints) != 1:
            raise ValueError('Libraries/evidence sources resolved incompatible guide designs; inspect orientation pilots')
        if 'guide_unit_id' not in group_cols:
            raise ValueError('guide_design.enabled requires design-aware Step 3 outputs')
    print(f"  Grouping by: {group_cols}")

    all_deduped = []
    lib_stats_all = {}
    ga = cfg['guide_assignment']
    dual_mode = ga.get('dual_guide_counting_mode', False)

    # Preserve the established multiprocessing path for the default arm.
    if source_mode == 'dedicated' and ncpus > 1 and len(tasks) > 1:
        from functools import partial
        from multiprocessing import Pool

        legacy_tasks = [
            (paths[0][1], lib_name, sample_name, sublib_idx)
            for paths, lib_name, sample_name, sublib_idx in tasks
        ]
        worker = partial(
            _process_one_library,
            group_cols=group_cols,
            max_umi_edits=max_umi_edits,
            min_reads_per_umi=min_reads_per_umi,
        )
        with Pool(min(ncpus, len(legacy_tasks))) as pool:
            results = pool.map(worker, legacy_tasks)
        for (deduped, stats), (_, lib_name, _, _) in zip(results, legacy_tasks):
            all_deduped.append(deduped)
            lib_stats_all[lib_name] = stats
    else:
        for paths, lib_name, sample_name, sublib_idx in tasks:
            if source_mode == 'dedicated':
                deduped, stats = stream_library_umi_dedup(
                    paths[0][1],
                    lib_name,
                    sample_name,
                    sublib_idx,
                    group_cols=group_cols,
                    max_umi_edits=max_umi_edits,
                    min_reads_per_umi=min_reads_per_umi,
                )
            else:
                deduped, stats = stream_library_umi_dedup_sources(
                    paths,
                    lib_name,
                    sample_name,
                    sublib_idx,
                    group_cols=group_cols,
                    max_umi_edits=max_umi_edits,
                    min_reads_per_umi=min_reads_per_umi,
                    require_resolved_ibar_sources=(
                        {'gex'} if dual_mode else set()
                    ),
                )
            all_deduped.append(deduped)
            lib_stats_all[lib_name] = stats

    nonempty = [frame for frame in all_deduped if len(frame) > 0]
    if not nonempty:
        if source_mode == 'gex':
            # Zero rejected pairs or zero ordinary Step-3 matches is a valid
            # Step 1c outcome, especially for short dual-guide reads with no
            # observable iBAR. Install a complete empty product instead of
            # turning a successful decontamination run into a pipeline error.
            dedup_columns = list(dict.fromkeys(
                group_cols + [
                    'umi_count', 'umi_count_raw', 'read_count',
                    'sublibrary', 'sample_name',
                    'exact_umi_count_gex', 'umi_count_gex',
                    'read_count_gex', 'exact_umi_count_overlap',
                    'umi_count_both_sources', 'umi_count_gex_only', 'source',
                ]
            ))
            assignment_columns = list(dict.fromkeys(
                ['bc_wells', 'grna_name', 'target_gene', 'umi_count']
            ))
            print(
                "  WARNING: no GEX-derived reads passed ordinary Step 3; "
                "writing a complete empty Step 1c identity product"
            )
            return {
                'mode': source_mode,
                'deduped': pd.DataFrame(columns=dedup_columns),
                'assignment': pd.DataFrame(columns=assignment_columns),
                'matrix': sparse.csr_matrix((0, 0), dtype=np.int64),
                'barcodes': [],
                'guides': [],
                'group_cols': group_cols,
                'lib_stats': lib_stats_all,
            }
        raise ValueError(
            f"No eligible {source_mode} guide evidence remained after filtering"
        )
    deduped_df = pd.concat(nonempty, ignore_index=True)
    if enabled(cfg) and set(deduped_df['guide_design_fingerprint']) != effective_fingerprints:
        raise ValueError('Match evidence does not match its resolved guide design')
    deduped_df, assignment_df = _annotate_and_assign(
        deduped_df, ref_df, available_cols, group_cols, ga
    )
    matrix, barcodes, guides = build_guide_umi_matrix(deduped_df)
    return {
        'mode': source_mode,
        'deduped': deduped_df,
        'assignment': assignment_df,
        'matrix': matrix,
        'barcodes': barcodes,
        'guides': guides,
        'group_cols': group_cols,
        'lib_stats': lib_stats_all,
    }


def _write_assignment_product(product, output_dir, gex_namespace):
    """Write a product using the documented source-specific filenames."""
    from scipy.io import mmwrite

    mode = product['mode']
    if mode == 'dedicated':
        assign_dir = os.path.join(output_dir, 'grna_assignment')
        names = {
            'assignment': 'grna_cell_assignment.csv',
            'counts': 'grna_umi_counts.csv.gz',
            'matrix': 'grna_umi_matrix.mtx',
            'barcodes': 'grna_matrix_barcodes.txt',
            'guides': 'grna_matrix_guides.txt',
        }
    elif mode == 'gex':
        assign_dir = os.path.join(
            output_dir, gex_namespace, 'grna_assignment'
        )
        names = {
            'assignment': 'gex_grna_assignments.csv.gz',
            'counts': 'gex_grna_umi_counts.csv.gz',
            'matrix': 'gex_grna_umi_matrix.mtx',
            'barcodes': 'gex_grna_matrix_barcodes.txt',
            'guides': 'gex_grna_matrix_guides.txt',
        }
    else:
        assign_dir = os.path.join(output_dir, 'grna_assignment')
        names = {
            'assignment': 'grna_cell_assignment_combined.csv',
            'counts': 'grna_umi_counts_by_source.csv.gz',
            'matrix': 'grna_umi_matrix_combined.mtx',
            'barcodes': 'grna_matrix_combined_barcodes.txt',
            'guides': 'grna_matrix_combined_guides.txt',
        }

    os.makedirs(assign_dir, exist_ok=True)
    product['assignment'].to_csv(
        os.path.join(assign_dir, names['assignment']),
        index=False,
        compression='gzip' if names['assignment'].endswith('.gz') else None,
    )
    product['deduped'].to_csv(
        os.path.join(assign_dir, names['counts']),
        index=False,
        compression='gzip',
    )
    mmwrite(os.path.join(assign_dir, names['matrix']), product['matrix'])
    pd.Series(product['barcodes']).to_csv(
        os.path.join(assign_dir, names['barcodes']), index=False, header=False
    )
    pd.Series(product['guides']).to_csv(
        os.path.join(assign_dir, names['guides']), index=False, header=False
    )

    summary = {
        'source_mode': mode,
        'combined_used_for_final_calls': mode == 'combined',
        'total_input_reads': int(sum(
            stats.get('total_reads', 0)
            for stats in product['lib_stats'].values()
        )),
        'total_umi_clusters': int(product['deduped']['umi_count'].sum()),
        'cells_with_evidence': int(product['deduped']['bc_wells'].nunique()),
        'cells_with_assignment': int(
            product['assignment']['bc_wells'].nunique()
            if 'bc_wells' in product['assignment'] else 0
        ),
        'libraries': product['lib_stats'],
    }
    if mode in ('gex', 'combined'):
        source_names = ('gex',) if mode == 'gex' else ('dedicated', 'gex')
        summary['source_counts'] = {}
        for source in source_names:
            summary['source_counts'][source] = {
                'reads_in_retained_union_umis': int(
                    product['deduped'].get(
                        f'read_count_{source}', pd.Series(dtype=int)
                    ).sum()
                ),
                'exact_umis_after_source_only_read_filter': int(
                    product['deduped'].get(
                        f'exact_umi_count_{source}', pd.Series(dtype=int)
                    ).sum()
                ),
                'collapsed_umis_after_source_only_read_filter': int(
                    product['deduped'].get(
                        f'umi_count_{source}', pd.Series(dtype=int)
                    ).sum()
                ),
                'union_clusters_from_source_only': int(
                    product['deduped'].get(
                        f'umi_count_{source}_only', pd.Series(dtype=int)
                    ).sum()
                ),
            }
        summary['overlap'] = {
            'exact_umis_observed_in_multiple_sources': int(
                product['deduped'].get(
                    'exact_umi_count_overlap', pd.Series(dtype=int)
                ).sum()
            ),
            'union_clusters_with_multiple_sources': int(
                product['deduped'].get(
                    'umi_count_both_sources', pd.Series(dtype=int)
                ).sum()
            ),
        }
        if mode == 'combined':
            summary['incremental_gex_union_clusters'] = int(
                product['deduped'].get(
                    'umi_count_gex_only', pd.Series(dtype=int)
                ).sum()
            )
    if mode == 'gex':
        summary_path = os.path.join(assign_dir, 'step1c_summary.yaml')
    elif mode == 'combined':
        summary_path = os.path.join(assign_dir, 'grna_gex_merge_summary.yaml')
    else:
        summary_path = os.path.join(assign_dir, 'grna_qc_stats.yaml')
    if mode != 'dedicated':
        with open(summary_path, 'w') as handle:
            yaml.safe_dump(summary, handle, sort_keys=False)
    print(f"  Saved {mode} assignment products under: {assign_dir}")
    return assign_dir


def main():
    parser = argparse.ArgumentParser(
        description="Assign gRNAs to cells via UMI deduplication"
    )
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument("--ncpus", type=int, default=1,
                        help="Number of CPUs for parallel library processing")
    parser.add_argument(
        "--source-mode", "--source_mode",
        choices=('auto', 'dedicated', 'gex', 'combined'),
        default='auto',
        help=(
            "Evidence product to build. auto preserves the dedicated-only "
            "workflow unless gex_decontamination Step1c/merge is enabled."
        ),
    )
    parser.add_argument(
        "--output-namespace", default=None,
        help="Override gex_decontamination.step1c.output_namespace",
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    output_dir = cfg['paths']['output_dir']
    ga = cfg['guide_assignment']
    dedup_ibar = ga.get('umi_dedup_ibar', True)
    max_umi_edits = ga.get('max_umi_edits', 1)
    min_reads_per_umi = ga.get('min_reads_per_umi', 1)
    compute_saturation = ga.get('compute_saturation_curve', True)
    decontam = cfg.get('gex_decontamination', {})
    step1c_cfg = decontam.get('step1c', {})
    combine_cfg = decontam.get('combine_with_grna_arm', {})
    gex_namespace = (
        args.output_namespace
        or step1c_cfg.get('output_namespace', 'gex_grna_identity')
    )

    # Auto-override: disable iBAR dedup when ibar_length == 0
    ibar_length = cfg.get('grna_structure', {}).get('ibar_length', 7)
    if ibar_length == 0 and dedup_ibar and not enabled(cfg):
        print("  NOTE: ibar_length=0 in config → overriding umi_dedup_ibar to False")
        dedup_ibar = False

    print(f"\n{'='*60}")
    print(f"  Step 4: Cell-Level gRNA Assignment")
    print(f"  iBAR length:             {ibar_length}")
    print(f"  UMI dedup iBAR grouping: {dedup_ibar}")
    print(f"  UMI max edits:           {max_umi_edits}")
    print(f"  Min reads per UMI:       {min_reads_per_umi}")
    print(f"  Saturation curve:        {compute_saturation}")
    print(f"  Parallel CPUs:           {args.ncpus}")
    print(f"  Requested source mode:   {args.source_mode}")
    print(f"{'='*60}")

    # Load reference table for design annotations
    ref_path = cfg['paths']['grna_reference']
    ref_df = None
    annotation_cols = ['grna_name', 'position', 'stem_variant', 'pair_id',
                       'target_combo', 'target_version', 'gi_pair', 'gi_version', 'guide_unit_id']
    if os.path.exists(ref_path):
        ref_df = pd.read_csv(ref_path).fillna('')
        if enabled(cfg):
            ref_df = pd.DataFrame(compile_design(ref_df, cfg)['rows'])
        available_cols = [c for c in annotation_cols if c in ref_df.columns]
        if available_cols:
            print(f"\n  Loaded reference annotations: {available_cols}")
        else:
            ref_df = None

    if args.source_mode == 'auto':
        modes = ['dedicated']
        step1c_enabled = bool(step1c_cfg.get('enabled', False))
        combine_enabled = bool(combine_cfg.get('enabled', False))
        if combine_enabled and not step1c_enabled:
            raise ValueError(
                "combine_with_grna_arm.enabled requires step1c.enabled"
            )
        if step1c_enabled:
            modes.append('gex')
        if combine_enabled:
            if not combine_cfg.get('combine_before_umi_collapse', True):
                raise ValueError(
                    "Only combine_before_umi_collapse=true is supported"
                )
            modes.append('combined')
    else:
        modes = [args.source_mode]

    products = {}
    for mode in modes:
        print(f"\n  Streaming {mode} UMI evidence across libraries...")
        products[mode] = _stream_assignment_product(
            cfg,
            output_dir,
            mode,
            gex_namespace,
            dedup_ibar,
            max_umi_edits,
            min_reads_per_umi,
            args.ncpus,
            ref_df,
            available_cols if ref_df is not None else [],
        )
        _write_assignment_product(products[mode], output_dir, gex_namespace)

    # The remaining historical QC products apply to the dedicated arm only.
    if 'dedicated' not in products:
        print(f"\n{'='*60}\n  Step 4 {modes[0]} product complete!\n{'='*60}")
        return

    dedicated = products['dedicated']
    deduped_df = dedicated['deduped']
    assignment_df = dedicated['assignment']
    guide_matrix = dedicated['matrix']
    cell_barcodes = dedicated['barcodes']
    guide_names = dedicated['guides']
    group_cols = dedicated['group_cols']
    lib_stats_all = dedicated['lib_stats']

    total_reads = sum(s['total_reads'] for s in lib_stats_all.values())
    total_umis_raw = sum(s.get('total_umis_raw', 0) for s in lib_stats_all.values())
    total_umis_after = sum(s.get('total_umis_after_collapse', 0) for s in lib_stats_all.values())
    overall_saturation = 1.0 - total_umis_after / total_reads if total_reads > 0 else 0.0
    print(f"\n  UMI dedup summary:")
    print(f"    Total reads streamed: {total_reads:>12,}")
    print(f"    Total UMIs (raw):     {total_umis_raw:>12,}")
    print(f"    UMIs after filter:    {sum(s.get('total_umis_after_filter', 0) for s in lib_stats_all.values()):>12,}")
    print(f"    UMIs after collapse:  {total_umis_after:>12,}")
    print(f"    Sequencing saturation:{overall_saturation:>11.1%}")
    print(f"    Unique groups:        {len(deduped_df):>12,}")
    print(f"    Unique cells:         {deduped_df['bc_wells'].nunique():>12,}")
    print(f"    Unique guides:        {deduped_df['grna_name'].nunique():>12,}")

    dual_guide_counting_mode = ga.get('dual_guide_counting_mode', False)

    # Compute QC stats
    qc_stats = compute_qc_stats(
        lib_stats_all, deduped_df, assignment_df,
        guide_matrix, cell_barcodes, guide_names,
        dual_guide_counting_mode=dual_guide_counting_mode
    )

    # --- Save outputs ---
    assign_dir = os.path.join(output_dir, 'grna_assignment')
    os.makedirs(assign_dir, exist_ok=True)

    # Guide annotation table (maps guide_name → design annotations)
    if ref_df is not None:
        annot_path = os.path.join(assign_dir, 'guide_annotations.csv')
        guide_annot = ref_df[available_cols].drop_duplicates()
        guide_annot.to_csv(annot_path, index=False)
        print(f"  Saved guide annotations: {annot_path}")

    # 5. QC stats
    qc_path = os.path.join(assign_dir, 'grna_qc_stats.yaml')
    with open(qc_path, 'w') as f:
        yaml.dump(qc_stats, f, default_flow_style=False)
    print(f"  Saved QC stats: {qc_path}")

    # 6. UMI threshold sensitivity plot
    plot_umi_threshold_sensitivity(
        deduped_df,
        save_path=os.path.join(assign_dir, 'umi_threshold_sensitivity.png')
    )

    # 7. Reads-per-UMI distribution plot
    plot_reads_per_umi_distribution(
        lib_stats_all,
        min_reads_per_umi=min_reads_per_umi,
        save_path=os.path.join(assign_dir, 'reads_per_umi_distribution.png')
    )

    # 8. Per-cell reads-per-UMI violin plot
    plot_median_reads_per_umi_per_cell(
        deduped_df,
        save_path=os.path.join(assign_dir, 'median_reads_per_umi_per_cell.png')
    )

    # 9. Sequencing saturation curve (optional, can be slow)
    if compute_saturation:
        print(f"\n  Computing sequencing saturation curve...")
        curve_data = _compute_saturation_curve(
            lib_stats_all,
            min_reads_per_umi=min_reads_per_umi
        )
        plot_saturation_curve(
            curve_data,
            save_path=os.path.join(assign_dir, 'saturation_curve.png')
        )
        # Add curve data to QC stats
        qc_stats['saturation_curve'] = {
            lib: [{k: v for k, v in pt.items()} for pt in pts]
            for lib, pts in curve_data.items()
        }

    # Print QC summary
    print(f"\n  QC Summary:")
    print(f"    Total matched reads:    {qc_stats['total_matched_reads']:>12,}")
    print(f"    Unique cells with reads:{qc_stats['total_unique_cells_with_reads']:>12,}")
    print(f"    Cells with assignment:  {qc_stats['cells_with_assignment']:>12,}")
    print(f"    Sequencing saturation:  {qc_stats.get('sequencing_saturation', 0):>11.1%}")
    print(f"    Median reads/UMI:       {qc_stats.get('median_reads_per_umi', 0):>12.1f}")
    
    if dual_guide_counting_mode:
        print(f"    Mean pairs/cell:        {qc_stats.get('mean_pairs_per_cell', 0):>12.2f}")
    else:
        print(f"    Mean guides/cell:       {qc_stats.get('mean_guides_per_cell', 0):>12.2f}")
        
    print(f"    Median UMI/cell-guide:  {qc_stats['median_umi_per_cell_guide']:>12.1f}")

    # Pair-level analysis (if dual guide mode)
    if dual_guide_counting_mode and len(assignment_df) > 0:
        print(f"\n  Pair-Level Analysis:")
        # assignment_df has one row per pair combination above threshold
        both_pos = (assignment_df['detection_status'] == 'complete').sum()
        inferred = (assignment_df['detection_status'] == 'inferred').sum()
        print(f"    Complete pairs:         {both_pos:>10,}")
        print(f"    Inferred pairs:         {inferred:>10,}")
        print(f"    Unique pair_ids:        {assignment_df['pair_id'].nunique():>10,}")

    print(f"\n{'='*60}")
    print(f"  Step 4 complete!")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()

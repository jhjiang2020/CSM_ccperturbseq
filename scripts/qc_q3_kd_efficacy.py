#!/usr/bin/env python3
"""
QC Q3: KD Efficacy + Inferred Source Classification

Assess per-guide knockdown efficacy using control cells as reference and
classify inferred singlets into Source 1 (dropout) vs Source 2 (recombinant).

Test sets:
  Set 1 (perfect_singlet): Both pos1 and pos2 detected, iBARs match a single
         correct combo in library design (combo_type=singlet, detection_status=complete).
  Set 2 (perfect_matched_multi): Both pos1 and pos2 detected, but iBARs indicate
         mismatched pos1+pos2 (recombination or undetected multiplet).

Reference: ctrl_ctrl singlet cells (both positions targeting nothing).

For X_ctrl cells → test KD of gene X (pos1 guide).
For ctrl_Y cells → test KD of gene Y (pos2 guide).

Q3c classifies inferred singlets, then repeats DE tests on Source 1 (true
dropout) inferred singlets to verify they show comparable KD to perfect singlets.

Metrics: log2FC, Wilcoxon p-value, Cohen's d, per guide and per sample.

Usage:
  python scripts/qc_q3_kd_efficacy.py --config config.yaml
"""

import argparse
import os
import sys
import warnings
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from scipy import stats


def load_config(config_path):
    with open(config_path) as f:
        return yaml.safe_load(f)


def cohens_d(group1, group2):
    """Compute Cohen's d effect size."""
    n1, n2 = len(group1), len(group2)
    if n1 < 2 or n2 < 2:
        return np.nan
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    if pooled_std == 0:
        return 0.0
    return (np.mean(group1) - np.mean(group2)) / pooled_std


def classify_cells_full(adata, verbose=True):
    """Comprehensive cell classification for dataset overview.

    Returns:
        dict mapping set_name → DataFrame of cell obs
        Also prints a summary table.

    Cell sets:
        Set 1 (perfect_singlet): singlet + complete
            → Both positions detected, exactly one guide combo.
        Set 2 (inferred_singlet): singlet + inferred
            → One position detected, other inferred from pair_id.
        Set 3 (perfect_multiplet): multiplet + at least some complete
            → Multiple guide combos, with both positions detected on ≥1 combo.
        Set 4 (inferred_multi_mismatch): multiplet + all complete, with mismatched pos1/pos2
            → Both positions detected but assigned to different pair IDs (recombination/multiplet).
        Set 5 (inferred_multi_single_pos): multiplet + all inferred
            → Multiple combos but only one position ever detected.
        Set 6 (unassigned): cells with no gRNA assignment.
    """
    obs = adata.obs.copy()
    assigned = obs[obs['has_grna'] == True].copy()
    unassigned = obs[obs['has_grna'] != True].copy()

    cell_sets = {}

    # --- Set 1: Perfect singlets (singlet + complete) ---
    set1_mask = (
        (assigned['combo_type'] == 'singlet') &
        (assigned['detection_status'] == 'complete')
    )
    cell_sets['perfect_singlet'] = assigned[set1_mask].copy()

    # --- Set 2: Inferred singlets (singlet + inferred) ---
    set2_mask = (
        (assigned['combo_type'] == 'singlet') &
        (assigned['detection_status'] == 'inferred')
    )
    cell_sets['inferred_singlet'] = assigned[set2_mask].copy()

    # --- Multiplets ---
    multi_mask = assigned['combo_type'] == 'multiplet'
    multiplets = assigned[multi_mask]

    # Set 3: Perfect multiplets — multiplet cells where ALL detection statuses are 'complete'
    # (both positions detected for every combo in the cell)
    set3_mask = multi_mask & (assigned['detection_status'] == 'complete')
    cell_sets['perfect_multiplet'] = assigned[set3_mask].copy()

    # Set 4: Multiplet mismatch — multiplet, contains both 'complete' and 'inferred' combos
    # This means both positions were detected, but the cell has multiple pair assignments
    # with mixed detection (some complete, some inferred → likely recombination)
    set4_mask = (
        multi_mask &
        assigned['detection_status'].str.contains('complete', na=False) &
        assigned['detection_status'].str.contains('inferred', na=False)
    )
    cell_sets['inferred_multi_mismatch'] = assigned[set4_mask].copy()

    # Set 5: Single-pos multiplets — multiplet where ALL detection statuses are 'inferred'
    # (only one position detected, the other always inferred)
    set5_mask = multi_mask & (assigned['detection_status'] == 'inferred')
    cell_sets['inferred_multi_single_pos'] = assigned[set5_mask].copy()

    # Catch remaining multiplets not classified above
    classified_multi = set3_mask | set4_mask | set5_mask
    remaining_multi = multi_mask & ~classified_multi
    if remaining_multi.sum() > 0:
        cell_sets['multiplet_other'] = assigned[remaining_multi].copy()

    # --- Set 6: Unassigned ---
    cell_sets['unassigned'] = unassigned.copy()

    # --- Summary ---
    if verbose:
        total = len(obs)
        print(f"\n  Comprehensive Cell Classification:")
        print(f"    {'Category':<30} {'Count':>8}  {'% of Total':>10}")
        print(f"    {'─' * 52}")
        for name, df in cell_sets.items():
            pct = len(df) / total * 100 if total > 0 else 0
            print(f"    {name:<30} {len(df):>8,}  {pct:>9.1f}%")
        print(f"    {'─' * 52}")
        print(f"    {'TOTAL':<30} {total:>8,}  {'100.0%':>10}")

        # Sanity check: classified assigned = all assigned
        classified = sum(len(df) for name, df in cell_sets.items() if name != 'unassigned')
        if classified != len(assigned):
            print(f"\n    ⚠ WARNING: {len(assigned) - classified:,} assigned cells not classified")

    return cell_sets


def classify_cells(adata, min_cells_per_set=None, group_by_gene_per_set=None):
    """Classify cells into the two primary test sets for KD analysis.

    Args:
        adata: AnnData object.
        min_cells_per_set: dict mapping set name to minimum cell count.
            Defaults to {'perfect_singlet': 10, 'inferred_singlet': 5}
        group_by_gene_per_set: dict mapping set name to bool for gene-level
            pooling. Defaults to {'perfect_singlet': False, 'inferred_singlet': True}

    Returns:
        dict of (DataFrame, min_cells, group_by_gene) tuples keyed by set name.
    """
    mc_defaults = {'perfect_singlet': 10, 'inferred_singlet': 5}
    gbg_defaults = {'perfect_singlet': False, 'inferred_singlet': True}
    if min_cells_per_set is not None:
        mc_defaults.update(min_cells_per_set)
    if group_by_gene_per_set is not None:
        gbg_defaults.update(group_by_gene_per_set)

    full = classify_cells_full(adata, verbose=True)
    return {
        name: (full[name], mc_defaults.get(name, 10), gbg_defaults.get(name, False))
        for name in ['perfect_singlet', 'inferred_singlet']
    }


def get_ctrl_ctrl_cells(cell_df):
    """Extract ctrl_ctrl cells from a classified set."""
    mask = cell_df['target_combo'] == 'ctrl_ctrl'
    ctrl = cell_df[mask]
    return ctrl


def get_grna_ctrl_cells(cell_df, group_by_gene=False):
    """Extract cells where one position is ctrl and the other is targeting.

    Args:
        cell_df: DataFrame of classified cells.
        group_by_gene: if True, pool all guides targeting the same gene into
            a single group (useful when per-guide cell counts are low).

    Returns list of dicts with:
        target_gene, guide_name, position, cell_barcodes, n_cells
    """
    groups = []

    target_combos = cell_df['target_combo'].unique()

    for tc in target_combos:
        if pd.isna(tc) or tc == 'ctrl_ctrl' or tc == 'unassigned':
            continue

        tc_cells = cell_df[cell_df['target_combo'] == tc]
        sample_row = tc_cells.iloc[0]

        pos1_target = str(sample_row.get('pos1_target', ''))
        pos2_target = str(sample_row.get('pos2_target', ''))

        pos1_targets = [t.strip() for t in pos1_target.split(',') if t.strip() and t.strip() != 'unassigned']
        pos2_targets = [t.strip() for t in pos2_target.split(',') if t.strip() and t.strip() != 'unassigned']

        # X_ctrl: pos1 is targeting, pos2 is ctrl
        if pos2_targets == ['ctrl'] and pos1_targets and pos1_targets != ['ctrl']:
            for guide_name in tc_cells['pos1_grna'].unique():
                if pd.isna(guide_name) or guide_name == '' or guide_name == 'unassigned':
                    continue
                guide_cells = tc_cells[tc_cells['pos1_grna'] == guide_name]
                for tg in pos1_targets:
                    groups.append({
                        'target_gene': tg,
                        'guide_name': guide_name,
                        'position': 'pos1',
                        'target_combo': tc,
                        'cell_barcodes': guide_cells.index.tolist(),
                        'n_cells': len(guide_cells),
                    })

        # ctrl_Y: pos1 is ctrl, pos2 is targeting
        elif pos1_targets == ['ctrl'] and pos2_targets and pos2_targets != ['ctrl']:
            for guide_name in tc_cells['pos2_grna'].unique():
                if pd.isna(guide_name) or guide_name == '' or guide_name == 'unassigned':
                    continue
                guide_cells = tc_cells[tc_cells['pos2_grna'] == guide_name]
                for tg in pos2_targets:
                    groups.append({
                        'target_gene': tg,
                        'guide_name': guide_name,
                        'position': 'pos2',
                        'target_combo': tc,
                        'cell_barcodes': guide_cells.index.tolist(),
                        'n_cells': len(guide_cells),
                    })

    if group_by_gene and groups:
        from collections import defaultdict
        gene_pool = defaultdict(list)
        gene_positions = {}
        for g in groups:
            gene_pool[g['target_gene']].extend(g['cell_barcodes'])
            gene_positions.setdefault(g['target_gene'], g['position'])

        groups = []
        for gene, barcodes in gene_pool.items():
            unique_bcs = list(set(barcodes))
            groups.append({
                'target_gene': gene,
                'guide_name': f'{gene}_all',
                'position': gene_positions[gene],
                'target_combo': f'{gene}_ctrl',
                'cell_barcodes': unique_bcs,
                'n_cells': len(unique_bcs),
            })

    return groups


def get_inferred_guide_cells(cell_df, group_by_gene=False):
    """Extract cells where the tested guide is at the INFERRED position.

    For inferred singlets, one position is detected and the other is inferred
    from the iBAR/pair_id. This function collects ALL cells (any partner —
    ctrl or non-ctrl) where the targeting guide was inferred, not detected.

    For each cell, the inferred position is identified by checking
    inferred_pos1_grna / inferred_pos2_grna columns:
      - If inferred_pos1_grna is populated → pos1 was inferred → test pos1 guide
      - If inferred_pos2_grna is populated → pos2 was inferred → test pos2 guide

    Barcodes for the same guide are pooled across all combos.

    Args:
        cell_df: DataFrame of classified cells (typically inferred singlets).
        group_by_gene: if True, further pool all guides targeting the same gene.

    Returns list of dicts with:
        target_gene, guide_name, position, cell_barcodes, n_cells
    """
    from collections import defaultdict

    guide_pool = defaultdict(list)  # (guide_name, target_gene, position) → [barcodes]

    has_inf_pos1 = 'inferred_pos1_grna' in cell_df.columns
    has_inf_pos2 = 'inferred_pos2_grna' in cell_df.columns

    if not has_inf_pos1 and not has_inf_pos2:
        return []

    for bc, row in cell_df.iterrows():
        tc = row.get('target_combo', '')
        if pd.isna(tc) or tc == 'ctrl_ctrl' or tc == 'unassigned':
            continue

        # Check if pos1 was inferred (and is a targeting guide, not ctrl)
        if has_inf_pos1:
            inf1 = row.get('inferred_pos1_grna', '')
            if pd.notna(inf1) and inf1 != '' and inf1 != 'unassigned':
                tg1 = row.get('inferred_pos1_target', row.get('pos1_target', ''))
                if pd.notna(tg1) and tg1 not in ('', 'unassigned', 'ctrl'):
                    guide_pool[(str(inf1), str(tg1), 'pos1')].append(bc)

        # Check if pos2 was inferred
        if has_inf_pos2:
            inf2 = row.get('inferred_pos2_grna', '')
            if pd.notna(inf2) and inf2 != '' and inf2 != 'unassigned':
                tg2 = row.get('inferred_pos2_target', row.get('pos2_target', ''))
                if pd.notna(tg2) and tg2 not in ('', 'unassigned', 'ctrl'):
                    guide_pool[(str(inf2), str(tg2), 'pos2')].append(bc)

    # Collapse into groups
    groups = []
    for (guide_name, target_gene, position), barcodes in guide_pool.items():
        unique_bcs = list(set(barcodes))
        groups.append({
            'target_gene': target_gene,
            'guide_name': guide_name,
            'position': position,
            'target_combo': f'{target_gene}_inferred',
            'cell_barcodes': unique_bcs,
            'n_cells': len(unique_bcs),
        })

    if group_by_gene and groups:
        gene_pool = defaultdict(list)
        gene_positions = {}
        for g in groups:
            gene_pool[g['target_gene']].extend(g['cell_barcodes'])
            gene_positions.setdefault(g['target_gene'], g['position'])

        groups = []
        for gene, barcodes in gene_pool.items():
            unique_bcs = list(set(barcodes))
            groups.append({
                'target_gene': gene,
                'guide_name': f'{gene}_all_inferred',
                'position': gene_positions[gene],
                'target_combo': f'{gene}_inferred',
                'cell_barcodes': unique_bcs,
                'n_cells': len(unique_bcs),
            })

    return groups


def get_grna_all_cells(cell_df, group_by_gene=False):
    """Extract ALL cells per gene, regardless of partner (ctrl or non-ctrl).

    Pools all barcodes per (guide, gene) across every combo where that gene
    is targeted. For perfect singlets, this gives maximum power per gene.

    Args:
        cell_df: DataFrame of classified cells.
        group_by_gene: if True, further pool all guides targeting the same gene.

    Returns list of dicts with:
        target_gene, guide_name, position, cell_barcodes, n_cells
    """
    from collections import defaultdict

    guide_pool = defaultdict(list)

    target_combos = cell_df['target_combo'].unique()

    for tc in target_combos:
        if pd.isna(tc) or tc == 'ctrl_ctrl' or tc == 'unassigned':
            continue

        tc_cells = cell_df[cell_df['target_combo'] == tc]
        sample_row = tc_cells.iloc[0]

        pos1_target = str(sample_row.get('pos1_target', ''))
        pos2_target = str(sample_row.get('pos2_target', ''))

        pos1_targets = [t.strip() for t in pos1_target.split(',')
                        if t.strip() and t.strip() not in ('unassigned', 'ctrl')]
        pos2_targets = [t.strip() for t in pos2_target.split(',')
                        if t.strip() and t.strip() not in ('unassigned', 'ctrl')]

        if pos1_targets:
            for guide_name in tc_cells['pos1_grna'].unique():
                if pd.isna(guide_name) or guide_name == '' or guide_name == 'unassigned':
                    continue
                guide_cells = tc_cells[tc_cells['pos1_grna'] == guide_name]
                for tg in pos1_targets:
                    guide_pool[(guide_name, tg, 'pos1')].extend(guide_cells.index.tolist())

        if pos2_targets:
            for guide_name in tc_cells['pos2_grna'].unique():
                if pd.isna(guide_name) or guide_name == '' or guide_name == 'unassigned':
                    continue
                guide_cells = tc_cells[tc_cells['pos2_grna'] == guide_name]
                for tg in pos2_targets:
                    guide_pool[(guide_name, tg, 'pos2')].extend(guide_cells.index.tolist())

    groups = []
    for (guide_name, target_gene, position), barcodes in guide_pool.items():
        unique_bcs = list(set(barcodes))
        groups.append({
            'target_gene': target_gene,
            'guide_name': guide_name,
            'position': position,
            'target_combo': f'{target_gene}_all',
            'cell_barcodes': unique_bcs,
            'n_cells': len(unique_bcs),
        })

    if group_by_gene and groups:
        gene_pool = defaultdict(list)
        gene_positions = {}
        for g in groups:
            gene_pool[g['target_gene']].extend(g['cell_barcodes'])
            gene_positions.setdefault(g['target_gene'], g['position'])

        groups = []
        for gene, barcodes in gene_pool.items():
            unique_bcs = list(set(barcodes))
            groups.append({
                'target_gene': gene,
                'guide_name': f'{gene}_all',
                'position': gene_positions[gene],
                'target_combo': f'{gene}_all',
                'cell_barcodes': unique_bcs,
                'n_cells': len(unique_bcs),
            })

    return groups


def get_detected_guide_cells(cell_df, group_by_gene=False):
    """Extract cells where the tested guide is at the DETECTED position.

    Complement of get_inferred_guide_cells: for inferred singlets, includes
    only cells where the targeting guide was directly detected (not inferred).

    For each cell:
      - If inferred_pos1_grna is populated → pos1 inferred, pos2 detected → test pos2
      - If inferred_pos2_grna is populated → pos2 inferred, pos1 detected → test pos1

    Args:
        cell_df: DataFrame of classified cells (typically inferred singlets).
        group_by_gene: if True, further pool all guides targeting the same gene.

    Returns list of dicts with:
        target_gene, guide_name, position, cell_barcodes, n_cells
    """
    from collections import defaultdict

    guide_pool = defaultdict(list)

    has_inf_pos1 = 'inferred_pos1_grna' in cell_df.columns
    has_inf_pos2 = 'inferred_pos2_grna' in cell_df.columns

    if not has_inf_pos1 and not has_inf_pos2:
        return []

    for bc, row in cell_df.iterrows():
        tc = row.get('target_combo', '')
        if pd.isna(tc) or tc == 'ctrl_ctrl' or tc == 'unassigned':
            continue

        # If pos1 was inferred → pos2 was detected → test pos2 guide
        pos1_was_inferred = False
        if has_inf_pos1:
            inf1 = row.get('inferred_pos1_grna', '')
            if pd.notna(inf1) and inf1 != '' and inf1 != 'unassigned':
                pos1_was_inferred = True
                # pos2 is the detected guide
                gn2 = row.get('pos2_grna', '')
                tg2 = row.get('pos2_target', '')
                if pd.notna(gn2) and str(gn2) not in ('', 'unassigned') and \
                   pd.notna(tg2) and str(tg2) not in ('', 'unassigned', 'ctrl'):
                    guide_pool[(str(gn2), str(tg2), 'pos2')].append(bc)

        # If pos2 was inferred → pos1 was detected → test pos1 guide
        if has_inf_pos2 and not pos1_was_inferred:
            inf2 = row.get('inferred_pos2_grna', '')
            if pd.notna(inf2) and inf2 != '' and inf2 != 'unassigned':
                gn1 = row.get('pos1_grna', '')
                tg1 = row.get('pos1_target', '')
                if pd.notna(gn1) and str(gn1) not in ('', 'unassigned') and \
                   pd.notna(tg1) and str(tg1) not in ('', 'unassigned', 'ctrl'):
                    guide_pool[(str(gn1), str(tg1), 'pos1')].append(bc)

    groups = []
    for (guide_name, target_gene, position), barcodes in guide_pool.items():
        unique_bcs = list(set(barcodes))
        groups.append({
            'target_gene': target_gene,
            'guide_name': guide_name,
            'position': position,
            'target_combo': f'{target_gene}_detected',
            'cell_barcodes': unique_bcs,
            'n_cells': len(unique_bcs),
        })

    if group_by_gene and groups:
        gene_pool = defaultdict(list)
        gene_positions = {}
        for g in groups:
            gene_pool[g['target_gene']].extend(g['cell_barcodes'])
            gene_positions.setdefault(g['target_gene'], g['position'])

        groups = []
        for gene, barcodes in gene_pool.items():
            unique_bcs = list(set(barcodes))
            groups.append({
                'target_gene': gene,
                'guide_name': f'{gene}_all_detected',
                'position': gene_positions[gene],
                'target_combo': f'{gene}_detected',
                'cell_barcodes': unique_bcs,
                'n_cells': len(unique_bcs),
            })

    return groups



def get_grna_any_cells(cell_df, group_by_gene=False):
    """Extract cells where a target gene is paired with a NON-ctrl partner.

    This is the complement of get_grna_ctrl_cells: it collects cells from
    X_Y combos where BOTH positions are targeting (excludes X_ctrl / ctrl_Y).
    Barcodes for the same guide across different non-ctrl combos are pooled
    into a single entry so each guide appears once.

    Args:
        cell_df: DataFrame of classified cells.
        group_by_gene: if True, further pool all guides targeting the same gene.

    Returns list of dicts with:
        target_gene, guide_name, position, cell_barcodes, n_cells
    """
    from collections import defaultdict

    # Accumulate per (guide_name, target_gene, position) across combos
    guide_pool = defaultdict(list)  # (guide_name, target_gene, position) → [barcodes]

    target_combos = cell_df['target_combo'].unique()

    for tc in target_combos:
        if pd.isna(tc) or tc == 'ctrl_ctrl' or tc == 'unassigned':
            continue

        tc_cells = cell_df[cell_df['target_combo'] == tc]
        sample_row = tc_cells.iloc[0]

        pos1_target = str(sample_row.get('pos1_target', ''))
        pos2_target = str(sample_row.get('pos2_target', ''))

        pos1_targets = [t.strip() for t in pos1_target.split(',')
                        if t.strip() and t.strip() not in ('unassigned', 'ctrl')]
        pos2_targets = [t.strip() for t in pos2_target.split(',')
                        if t.strip() and t.strip() not in ('unassigned', 'ctrl')]

        # Skip if either position is ctrl (those belong to ctrl_only mode)
        pos1_raw = [t.strip() for t in pos1_target.split(',') if t.strip() and t.strip() != 'unassigned']
        pos2_raw = [t.strip() for t in pos2_target.split(',') if t.strip() and t.strip() != 'unassigned']

        is_pos1_ctrl = pos1_raw == ['ctrl']
        is_pos2_ctrl = pos2_raw == ['ctrl']

        # Only include combos where NEITHER position is ctrl
        if is_pos1_ctrl or is_pos2_ctrl:
            continue

        # Pos1 targeting genes → pool barcodes per guide
        if pos1_targets:
            for guide_name in tc_cells['pos1_grna'].unique():
                if pd.isna(guide_name) or guide_name == '' or guide_name == 'unassigned':
                    continue
                guide_cells = tc_cells[tc_cells['pos1_grna'] == guide_name]
                for tg in pos1_targets:
                    guide_pool[(guide_name, tg, 'pos1')].extend(guide_cells.index.tolist())

        # Pos2 targeting genes → pool barcodes per guide
        if pos2_targets:
            for guide_name in tc_cells['pos2_grna'].unique():
                if pd.isna(guide_name) or guide_name == '' or guide_name == 'unassigned':
                    continue
                guide_cells = tc_cells[tc_cells['pos2_grna'] == guide_name]
                for tg in pos2_targets:
                    guide_pool[(guide_name, tg, 'pos2')].extend(guide_cells.index.tolist())

    # Collapse per-guide pools into unique barcode sets
    groups = []
    for (guide_name, target_gene, position), barcodes in guide_pool.items():
        unique_bcs = list(set(barcodes))
        groups.append({
            'target_gene': target_gene,
            'guide_name': guide_name,
            'position': position,
            'target_combo': f'{target_gene}_nonctrl',
            'cell_barcodes': unique_bcs,
            'n_cells': len(unique_bcs),
        })

    if group_by_gene and groups:
        gene_pool = defaultdict(list)
        gene_positions = {}
        for g in groups:
            gene_pool[g['target_gene']].extend(g['cell_barcodes'])
            gene_positions.setdefault(g['target_gene'], g['position'])

        groups = []
        for gene, barcodes in gene_pool.items():
            unique_bcs = list(set(barcodes))
            groups.append({
                'target_gene': gene,
                'guide_name': f'{gene}_all',
                'position': gene_positions[gene],
                'target_combo': f'{gene}_nonctrl',
                'cell_barcodes': unique_bcs,
                'n_cells': len(unique_bcs),
            })

    return groups


def compute_kd_metrics(adata, guide_groups, ctrl_barcodes, sample_label, min_cells=10):
    """Compute KD metrics for each guide group vs ctrl_ctrl.

    Args:
        adata: AnnData (normalized)
        guide_groups: list of dicts from get_grna_ctrl_cells
        ctrl_barcodes: list of ctrl_ctrl cell barcodes
        sample_label: 'Control' or 'Treatment'
        min_cells: minimum number of KD cells required (default 10)

    Returns:
        DataFrame with one row per guide
    """
    results = []

    # Get ctrl expression (all genes, but we'll extract per target)
    ctrl_cells = adata[adata.obs.index.isin(ctrl_barcodes)]

    if len(ctrl_cells) < 3:
        print(f"    WARNING: Only {len(ctrl_cells)} ctrl_ctrl cells for {sample_label}")
        return pd.DataFrame()

    for group in guide_groups:
        target_gene = group['target_gene']
        guide_name = group['guide_name']
        position = group['position']
        cell_bcs = group['cell_barcodes']

        # Check if target gene exists in var
        if target_gene not in adata.var_names:
            # Try appending _hg38
            if f"{target_gene}_hg38" in adata.var_names:
                target_gene = f"{target_gene}_hg38"
            else:
                # Try case-insensitive or base name match
                matches = [g for g in adata.var_names if g.split('_')[0].upper() == target_gene.split('_')[0].upper()]
                if matches:
                    target_gene = matches[0]
                else:
                    results.append({
                        'guide_name': guide_name,
                        'target_gene': group['target_gene'],
                        'position': position,
                    'target_combo': group['target_combo'],
                    'sample': sample_label,
                    'n_kd_cells': len(cell_bcs),
                    'n_ctrl_cells': len(ctrl_cells),
                    'cpm_kd': np.nan,
                    'cpm_ctrl': np.nan,
                    'pct_nonzero_kd': np.nan,
                    'pct_nonzero_ctrl': np.nan,
                    'log2FC': np.nan,
                    'pvalue': np.nan,
                    'cohens_d': np.nan,
                    'note': f'target gene {group["target_gene"]} not found in var_names',
                })
                continue

        # Extract expression
        kd_cells = adata[adata.obs.index.isin(cell_bcs)]

        if len(kd_cells) < min_cells:
            results.append({
                'guide_name': guide_name,
                'target_gene': target_gene,
                'position': position,
                'target_combo': group['target_combo'],
                'sample': sample_label,
                'n_kd_cells': len(kd_cells),
                'n_ctrl_cells': len(ctrl_cells),
                'cpm_kd': np.nan,
                'cpm_ctrl': np.nan,
                'pct_nonzero_kd': np.nan,
                'pct_nonzero_ctrl': np.nan,
                'log2FC': np.nan,
                'pvalue': np.nan,
                'cohens_d': np.nan,
                'note': f'insufficient cells (<{min_cells})',
            })
            continue

        kd_expr = kd_cells[:, target_gene].X
        ctrl_expr = ctrl_cells[:, target_gene].X

        # Convert sparse to dense
        if hasattr(kd_expr, 'toarray'):
            kd_expr = kd_expr.toarray().flatten()
        else:
            kd_expr = np.asarray(kd_expr).flatten()
        if hasattr(ctrl_expr, 'toarray'):
            ctrl_expr = ctrl_expr.toarray().flatten()
        else:
            ctrl_expr = np.asarray(ctrl_expr).flatten()

        # Revert log1p and CP10k to calculate linear Pseudobulk CPM for the output table
        # log1p(CP10k) -> expm1 -> CP10k -> *100 -> CPM
        kd_cpm = np.expm1(kd_expr) * 100
        ctrl_cpm = np.expm1(ctrl_expr) * 100
        
        cpm_kd = float(np.mean(kd_cpm))
        cpm_ctrl = float(np.mean(ctrl_cpm))
        
        pct_nonzero_kd = float(np.mean(kd_expr > 0)) * 100
        pct_nonzero_ctrl = float(np.mean(ctrl_expr > 0)) * 100

        # log2FC calculated on normalized log-scale expression directly
        mean_kd_log = float(np.mean(kd_expr))
        mean_ctrl_log = float(np.mean(ctrl_expr))
        pseudo = 0.01
        log2fc = float(np.log2((mean_kd_log + pseudo) / (mean_ctrl_log + pseudo)))

        # Wilcoxon rank-sum test
        if len(kd_expr) >= 3 and len(ctrl_expr) >= 3:
            try:
                stat, pval = stats.mannwhitneyu(kd_expr, ctrl_expr, alternative='two-sided')
            except ValueError:
                pval = np.nan
        else:
            pval = np.nan

        # Cohen's d
        d = cohens_d(kd_expr, ctrl_expr)

        results.append({
            'guide_name': guide_name,
            'target_gene': target_gene,
            'position': position,
            'target_combo': group['target_combo'],
            'sample': sample_label,
            'n_kd_cells': len(kd_cells),
            'n_ctrl_cells': len(ctrl_cells),
            'cpm_kd': cpm_kd,
            'cpm_ctrl': cpm_ctrl,
            'pct_nonzero_kd': pct_nonzero_kd,
            'pct_nonzero_ctrl': pct_nonzero_ctrl,
            'log2FC': log2fc,
            'pvalue': pval,
            'cohens_d': d,
            'note': '',
        })

    return pd.DataFrame(results)


def plot_heatmap(result_df, save_path, title_suffix='', show_pval=False):
    """Plot log2FC heatmap per guide × target gene."""
    if result_df.empty:
        return

    # Only include guides with valid log2FC
    valid = result_df.dropna(subset=['log2FC'])
    if valid.empty:
        return

    for sample in valid['sample'].unique():
        sample_df = valid[valid['sample'] == sample]

        sample_df = sample_df.copy()

        # Pivot: rows=guide_name, columns=target_gene
        pivot_fc = sample_df.pivot_table(
            index='guide_name', columns='target_gene',
            values='log2FC', aggfunc='first'
        )
        pivot_n = sample_df.pivot_table(
            index='guide_name', columns='target_gene', values='n_kd_cells', aggfunc='first'
        )
        
        pivot_cpm_kd = sample_df.pivot_table(
            index='guide_name', columns='target_gene', values='cpm_kd', aggfunc='first'
        )

        pivot_pval = None
        if show_pval and 'pvalue' in sample_df.columns:
            pivot_pval = sample_df.pivot_table(
                index='guide_name', columns='target_gene', values='pvalue', aggfunc='first'
            )
        
        # Calculate mean ctrl CPM per target gene for column labels
        avg_cpm_ctrl = sample_df.groupby('target_gene')['cpm_ctrl'].mean()

        if pivot_fc.empty or pivot_fc.shape[0] < 1:
            continue

        # Skip if all values are NaN
        if np.all(np.isnan(pivot_fc.values)):
            continue

        # Increase row height when showing pval to avoid crowding
        row_height = 1.8 if show_pval else 1.5
        fig, ax = plt.subplots(figsize=(max(8, len(pivot_fc.columns) * 0.8),
                                         max(4, len(pivot_fc.index) * row_height)))

        # Sort by target gene for readability
        pivot_fc = pivot_fc.reindex(sorted(pivot_fc.columns), axis=1)
        pivot_n = pivot_n.reindex(sorted(pivot_n.columns), axis=1)
        pivot_cpm_kd = pivot_cpm_kd.reindex(sorted(pivot_cpm_kd.columns), axis=1)
        if pivot_pval is not None:
            pivot_pval = pivot_pval.reindex(pivot_fc.index, axis=0).reindex(pivot_fc.columns, axis=1)

        # Sort guides (gRNA_1, gRNA_2, gRNA_3)
        pivot_fc = pivot_fc.reindex(sorted(pivot_fc.index), axis=0)
        pivot_n = pivot_n.reindex(sorted(pivot_n.index), axis=0)
        pivot_cpm_kd = pivot_cpm_kd.reindex(sorted(pivot_cpm_kd.index), axis=0)

        im = ax.imshow(pivot_fc.values, aspect='auto', cmap='RdBu_r',
                       vmin=-3, vmax=3)
                       
        # Build x-tick labels with target gene and baseline ctrl CPM
        x_labels = []
        for col in pivot_fc.columns:
            ctrl_val = avg_cpm_ctrl.get(col, np.nan)
            label = f"{col}\n(Ctrl CPM={ctrl_val:.1f})" if not np.isnan(ctrl_val) else col
            x_labels.append(label)
                       
        ax.set_xticks(range(len(pivot_fc.columns)))
        ax.set_xticklabels(x_labels, rotation=45, ha='right', rotation_mode='anchor', fontsize=9)
        ax.set_yticks(range(len(pivot_fc.index)))
        ax.set_yticklabels(pivot_fc.index, fontsize=8)
        ax.set_title(f'Guide KD Efficacy — {sample}{title_suffix}', fontsize=12)
        ax.set_xlabel('Target Gene')
        ax.set_ylabel('Guide')

        cbar = plt.colorbar(im, ax=ax, shrink=0.6)
        cbar.set_label('log₂FC vs ctrl_ctrl')

        # Add text annotations
        max_abs_fc = np.nanmax(np.abs(pivot_fc.values))
        font_size = 6 if show_pval else 7
        for i in range(len(pivot_fc.index)):
            for j in range(len(pivot_fc.columns)):
                value = pivot_fc.values[i, j]
                n_cells = pivot_n.values[i, j]
                cpm_k = pivot_cpm_kd.values[i, j]

                if not np.isnan(value) and not np.isnan(n_cells) and not np.isnan(cpm_k):
                    # Make text white if background is dark
                    text_color = "white" if abs(value) > max_abs_fc * 0.5 else "black"
                    
                    annot_text = f"lfc: {value:.2f}\ncpm: {cpm_k:.1f}\nn: {int(n_cells)}"
                    if show_pval and pivot_pval is not None:
                        pv = pivot_pval.values[i, j]
                        if not np.isnan(pv):
                            annot_text += f"\np: {pv:.1e}"
                    ax.text(j, i, annot_text, 
                            ha='center', va='center', color=text_color, fontsize=font_size)

        plt.tight_layout()
        suffix = sample.lower().replace(' ', '_')
        fig.savefig(save_path.replace('.png', f'_{suffix}.png'), dpi=150,
                    bbox_inches='tight')
        plt.close(fig)
        print(f"    Saved heatmap: {save_path.replace('.png', f'_{suffix}.png')}")


def plot_q1b_heatmap(result_df, save_path, title_suffix='', show_pval=True):
    """Combined Q1b heatmap: 3 fixed rows × gene columns, one figure per sample.

    Rows (in order):
        1. perfect_all         — perfect singlets, all guide combos  (positive control)
        2. detected_in_inferred — guide directly detected in inferred singlets
        3. inferred_in_inferred — guide inferred from partner iBAR in inferred singlets

    Columns: target genes, annotated with ctrl_ctrl CPM in the column label.
    Each cell shows: lfc, cpm_kd, n, and optionally p-value.

    One file per sample saved as:
        save_path → save_path.replace('.png', '_control.png') etc.
    """
    ROW_ORDER = [
        ('perfect_all',          'Perfect singlet (ctrl)'),
        ('detected_in_inferred', 'Detected in inferred'),
        ('inferred_in_inferred', 'Inferred in inferred'),
    ]

    if result_df.empty:
        return

    valid = result_df.dropna(subset=['log2FC'])
    if valid.empty:
        return

    for sample in sorted(valid['sample'].unique()):
        sample_df = valid[valid['sample'] == sample].copy()

        # Collect all target genes across all rows so columns are consistent
        all_genes = sorted(sample_df['target_gene'].unique())
        if not all_genes:
            continue

        # Build per-row pivot data — keyed by (mode_key, target_gene)
        row_data = {}  # mode_key -> {gene: {log2FC, cpm_kd, n_kd_cells, pvalue}}
        for mode_key, _label in ROW_ORDER:
            mode_df = sample_df[sample_df['test_mode'] == mode_key]
            if mode_df.empty:
                row_data[mode_key] = {}
                continue
            # Use gene-level aggregate (one row per target_gene) — take the
            # gene-pooled (guide_name == target_gene) row when present, else any.
            agg = {}
            for gene in mode_df['target_gene'].unique():
                g = mode_df[mode_df['target_gene'] == gene]
                # prefer the pooled (guide_name matches gene base) entry
                pooled = g[g['guide_name'].str.startswith(gene.split('_')[0])]
                best = pooled.iloc[0] if not pooled.empty else g.iloc[0]
                agg[gene] = {
                    'log2FC':     best['log2FC'],
                    'cpm_kd':     best.get('cpm_kd', np.nan),
                    'n_kd_cells': best.get('n_kd_cells', np.nan),
                    'pvalue':     best.get('pvalue', np.nan),
                }
            row_data[mode_key] = agg

        # Derive per-gene ctrl CPM from perfect_all row (most reliable)
        ctrl_cpm_by_gene = {}
        perf_df = sample_df[sample_df['test_mode'] == 'perfect_all']
        for gene in all_genes:
            g = perf_df[perf_df['target_gene'] == gene]
            if not g.empty and 'cpm_ctrl' in g.columns:
                ctrl_cpm_by_gene[gene] = g['cpm_ctrl'].mean()
            else:
                # fallback: any row
                g2 = sample_df[sample_df['target_gene'] == gene]
                ctrl_cpm_by_gene[gene] = g2['cpm_ctrl'].mean() if 'cpm_ctrl' in g2.columns else np.nan

        n_rows = len(ROW_ORDER)
        n_cols = len(all_genes)

        # Build value matrices
        fc_mat   = np.full((n_rows, n_cols), np.nan)
        cpm_mat  = np.full((n_rows, n_cols), np.nan)
        n_mat    = np.full((n_rows, n_cols), np.nan)
        pval_mat = np.full((n_rows, n_cols), np.nan)

        for ri, (mode_key, _) in enumerate(ROW_ORDER):
            agg = row_data.get(mode_key, {})
            for ci, gene in enumerate(all_genes):
                if gene in agg:
                    fc_mat[ri, ci]   = agg[gene]['log2FC']
                    cpm_mat[ri, ci]  = agg[gene]['cpm_kd']
                    n_mat[ri, ci]    = agg[gene]['n_kd_cells']
                    pval_mat[ri, ci] = agg[gene]['pvalue']

        # Skip if all values are NaN (nothing to show)
        if np.all(np.isnan(fc_mat)):
            continue

        row_height = 2.0 if show_pval else 1.6
        fig, ax = plt.subplots(
            figsize=(max(8, n_cols * 0.9), max(4, n_rows * row_height))
        )

        im = ax.imshow(fc_mat, aspect='auto', cmap='RdBu_r', vmin=-3, vmax=3)

        # Column labels: target gene + ctrl CPM
        x_labels = []
        for gene in all_genes:
            ctrl_val = ctrl_cpm_by_gene.get(gene, np.nan)
            if not np.isnan(ctrl_val):
                x_labels.append(f"{gene}\n(Ctrl={ctrl_val:.1f} CPM)")
            else:
                x_labels.append(gene)

        ax.set_xticks(range(n_cols))
        ax.set_xticklabels(x_labels, rotation=45, ha='right',
                            rotation_mode='anchor', fontsize=8)

        row_labels = [label for _, label in ROW_ORDER]
        ax.set_yticks(range(n_rows))
        ax.set_yticklabels(row_labels, fontsize=9)

        ax.set_title(f'Q1b KD Efficacy Across Detection Modes — {sample}{title_suffix}',
                     fontsize=12)
        ax.set_xlabel('Target Gene')
        ax.set_ylabel('Guide Detection Mode')

        cbar = plt.colorbar(im, ax=ax, shrink=0.6)
        cbar.set_label('log₂FC vs ctrl_ctrl')

        # Annotate cells
        max_abs_fc = np.nanmax(np.abs(fc_mat)) if not np.all(np.isnan(fc_mat)) else 1.0
        font_size = 6 if show_pval else 7
        for ri in range(n_rows):
            for ci in range(n_cols):
                val = fc_mat[ri, ci]
                if np.isnan(val):
                    ax.text(ci, ri, 'n/a', ha='center', va='center',
                            color='grey', fontsize=font_size)
                    continue
                text_color = 'white' if abs(val) > max_abs_fc * 0.5 else 'black'
                cpm_k = cpm_mat[ri, ci]
                n_k   = n_mat[ri, ci]
                annot = f"lfc: {val:.2f}"
                if not np.isnan(cpm_k):
                    annot += f"\ncpm: {cpm_k:.1f}"
                if not np.isnan(n_k):
                    annot += f"\nn: {int(n_k)}"
                if show_pval and not np.isnan(pval_mat[ri, ci]):
                    annot += f"\np: {pval_mat[ri, ci]:.1e}"
                ax.text(ci, ri, annot,
                        ha='center', va='center', color=text_color, fontsize=font_size)

        plt.tight_layout()
        suffix = sample.lower().replace(' ', '_')
        out = save_path.replace('.png', f'_{suffix}.png')
        fig.savefig(out, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"    Saved heatmap: {out}")


def plot_volcano(result_df, save_path, title_suffix='', max_neg_log10p=20):
    """Volcano plot: log2FC vs -log10(p-value).

    Args:
        max_neg_log10p: cap -log10(p) at this value (default 20) to avoid
            extreme values compressing the plot.
    """
    if result_df.empty:
        return

    valid = result_df.dropna(subset=['log2FC', 'pvalue'])
    valid = valid[valid['pvalue'] > 0]
    if valid.empty:
        return

    for sample in valid['sample'].unique():
        sample_df = valid[valid['sample'] == sample].copy()
        if sample_df.empty:
            continue

        sample_df['neg_log10p'] = -np.log10(sample_df['pvalue'].clip(lower=1e-300))
        sample_df['neg_log10p'] = sample_df['neg_log10p'].clip(upper=max_neg_log10p)

        fig, ax = plt.subplots(figsize=(10, 7))

        # Color by significance
        sig_mask = (sample_df['pvalue'] < 0.05) & (abs(sample_df['log2FC']) > 0.5)
        nonsig = sample_df[~sig_mask]
        sig = sample_df[sig_mask]

        ax.scatter(nonsig['log2FC'], nonsig['neg_log10p'], alpha=0.4,
                   s=30, c='grey', label='NS')
        ax.scatter(sig['log2FC'], sig['neg_log10p'], alpha=0.7,
                   s=50, c='#e74c3c', label='Significant')

        # Label significant guides
        for _, row in sig.iterrows():
            ax.annotate(row['guide_name'], (row['log2FC'], row['neg_log10p']),
                        fontsize=6, alpha=0.8, ha='center', va='bottom')

        ax.axhline(-np.log10(0.05), ls='--', color='grey', alpha=0.5, label='p=0.05')
        ax.axvline(0, ls='--', color='grey', alpha=0.3)

        ax.set_xlabel('log₂FC vs ctrl_ctrl', fontsize=12)
        ax.set_ylabel('-log₁₀(p-value)', fontsize=12)
        ax.set_title(f'Guide KD Efficacy Volcano — {sample}{title_suffix}', fontsize=13)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.2)

        plt.tight_layout()
        suffix = sample.lower().replace(' ', '_')
        fig.savefig(save_path.replace('.png', f'_{suffix}.png'), dpi=150,
                    bbox_inches='tight')
        plt.close(fig)
        print(f"    Saved volcano: {save_path.replace('.png', f'_{suffix}.png')}")


def plot_volcano_comparison(combined_df, save_path, title_suffix='', max_neg_log10p=20):
    """Multi-group volcano: overlay 3 test groups in different colors.

    Expects combined_df with a 'test_mode' column containing up to 3 groups:
      - 'perfect_all':           perfect singlets, all combos (positive control)
      - 'detected_in_inferred':  guide directly detected in inferred singlets
      - 'inferred_in_inferred':  guide inferred from partner iBAR in inferred singlets

    Args:
        combined_df: DataFrame with columns log2FC, pvalue, test_mode, target_gene, sample
        save_path: file to save the figure
        title_suffix: optional title suffix
        max_neg_log10p: cap -log10(p) at this value (default 20)
    """
    if combined_df.empty:
        return

    valid = combined_df.dropna(subset=['log2FC', 'pvalue'])
    valid = valid[valid['pvalue'] > 0]
    if valid.empty:
        return

    group_styles = {
        'perfect_all':           {'color': '#2196F3', 'label': 'Perfect singlet (all)',    'marker': 'o'},
        'detected_in_inferred':  {'color': '#4CAF50', 'label': 'Detected (inferred set)',  'marker': 's'},
        'inferred_in_inferred':  {'color': '#FF9800', 'label': 'Inferred (inferred set)',  'marker': '^'},
    }

    for sample in valid['sample'].unique():
        sample_df = valid[valid['sample'] == sample].copy()
        if sample_df.empty:
            continue

        sample_df['neg_log10p'] = -np.log10(sample_df['pvalue'].clip(lower=1e-300))
        sample_df['neg_log10p'] = sample_df['neg_log10p'].clip(upper=max_neg_log10p)

        fig, ax = plt.subplots(figsize=(12, 8))

        # Plot each group
        for mode, style in group_styles.items():
            gdf = sample_df[sample_df['test_mode'] == mode]
            if gdf.empty:
                continue

            ax.scatter(gdf['log2FC'], gdf['neg_log10p'],
                       alpha=0.7, s=50, c=style['color'],
                       marker=style['marker'], label=style['label'],
                       edgecolors='white', linewidths=0.5, zorder=3)

            # Label points
            for _, row in gdf.iterrows():
                label = row.get('target_gene', row.get('guide_name', ''))
                ax.annotate(label, (row['log2FC'], row['neg_log10p']),
                            fontsize=6, alpha=0.7, ha='center', va='bottom',
                            color=style['color'])

        # Connect the same gene across groups with faint lines
        genes_in_data = sample_df['target_gene'].unique()
        for gene in genes_in_data:
            gene_rows = sample_df[sample_df['target_gene'] == gene]
            if len(gene_rows) > 1:
                xs = gene_rows['log2FC'].values
                ys = gene_rows['neg_log10p'].values
                ax.plot(xs, ys, color='grey', alpha=0.2, linewidth=0.8, zorder=1)

        ax.axhline(-np.log10(0.05), ls='--', color='grey', alpha=0.5, label='p=0.05')
        ax.axvline(0, ls='--', color='grey', alpha=0.3)

        ax.set_xlabel('log₂FC vs ctrl_ctrl', fontsize=12)
        ax.set_ylabel('-log₁₀(p-value)', fontsize=12)
        ax.set_title(f'KD Efficacy Comparison — {sample}{title_suffix}', fontsize=13)
        ax.legend(fontsize=9, loc='upper right')
        ax.grid(True, alpha=0.2)

        plt.tight_layout()
        suffix = sample.lower().replace(' ', '_')
        out = save_path.replace('.png', f'_{suffix}.png')
        fig.savefig(out, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"    Saved comparison volcano: {out}")


def plot_detected_vs_inferred(combined_results, save_path, show=False):
    """Scatter plot comparing detected (Set 1) vs inferred (Set 2) log2FC per guide.

    Args:
        combined_results: DataFrame with 'test_set' column containing both
                          'perfect_singlet' and 'inferred_singlet' results.
        save_path: base path for saving figures.
        show: if True, call plt.show() instead of plt.close().
    """
    if combined_results.empty:
        return

    det = combined_results[combined_results['test_set'] == 'perfect_singlet'].copy()
    inf = combined_results[combined_results['test_set'] == 'inferred_singlet'].copy()

    if det.empty or inf.empty:
        print("    Cannot generate detected vs inferred plot: missing data for one set.")
        return

    # Merge on guide_name + target_gene + sample
    merge_keys = ['guide_name', 'target_gene', 'sample']
    merged = pd.merge(
        det[merge_keys + ['log2FC', 'n_kd_cells', 'pvalue']],
        inf[merge_keys + ['log2FC', 'n_kd_cells', 'pvalue']],
        on=merge_keys, suffixes=('_detected', '_inferred'), how='inner'
    )

    if merged.empty:
        print("    No overlapping guides between detected and inferred sets.")
        return

    for sample in merged['sample'].unique():
        sdf = merged[merged['sample'] == sample].dropna(
            subset=['log2FC_detected', 'log2FC_inferred']
        )
        if len(sdf) < 2:
            continue

        fig, ax = plt.subplots(figsize=(8, 8))

        ax.scatter(sdf['log2FC_detected'], sdf['log2FC_inferred'],
                   alpha=0.6, s=40, c='#2c3e50', edgecolors='white', linewidth=0.5)

        # Label each point with guide_name
        for _, row in sdf.iterrows():
            ax.annotate(row['guide_name'], (row['log2FC_detected'], row['log2FC_inferred']),
                        fontsize=6, alpha=0.7, ha='left', va='bottom')

        # Identity line
        lim_min = min(ax.get_xlim()[0], ax.get_ylim()[0])
        lim_max = max(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([lim_min, lim_max], [lim_min, lim_max], 'k--', alpha=0.3, label='y=x')

        # Correlation
        from scipy.stats import pearsonr
        r, p = pearsonr(sdf['log2FC_detected'], sdf['log2FC_inferred'])
        ax.set_title(f'Detected vs Inferred KD — {sample}\n(r={r:.3f}, p={p:.2e}, n={len(sdf)} guides)',
                     fontsize=12)

        ax.set_xlabel('log₂FC (Detected / Set 1)', fontsize=11)
        ax.set_ylabel('log₂FC (Inferred / Set 2)', fontsize=11)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.2)
        ax.set_aspect('equal', adjustable='datalim')

        plt.tight_layout()
        suffix = sample.lower().replace(' ', '_')
        out_path = save_path.replace('.png', f'_{suffix}.png')
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
        print(f"    Saved detected vs inferred scatter: {out_path}")
        if show:
            plt.show()
        else:
            plt.close(fig)



def _resolve_gene_name(gene, var_names):
    """Resolve a gene name to its match in adata.var_names (handles _hg38 suffix)."""
    if gene in var_names:
        return gene
    if f"{gene}_hg38" in var_names:
        return f"{gene}_hg38"
    matches = [g for g in var_names if g.split('_')[0].upper() == gene.upper()]
    return matches[0] if matches else None


def _extract_expr(adata, barcodes, gene_name):
    """Extract expression array for a given gene from a set of barcodes."""
    cells = adata[adata.obs.index.isin(barcodes)]
    if len(cells) == 0:
        return np.array([])
    expr = cells[:, gene_name].X
    if hasattr(expr, 'toarray'):
        expr = expr.toarray().flatten()
    else:
        expr = np.asarray(expr).flatten()
    return expr



# ============================================================================
# Generalized Expression Comparison Utility
# ============================================================================

def compare_expression(
    adata, test_barcodes, ctrl_barcodes, readout_genes,
    test_label="KD", ctrl_label="ctrl_ctrl",
    plot_type="both",        # "heatmap", "violin", "both", or None
    show_pval=True,          # annotate p-values on plots
    save_path=None,          # optional save path (without extension for 'both')
    title=None,              # optional title
    show=True,               # show or close figure
    figsize=None,            # optional (w, h) tuple
):
    """Compare gene expression between test and control cell groups.

    A general-purpose utility for visualizing expression differences.
    Computes log2FC, Mann-Whitney U p-value, and Cohen's d for each gene.

    Args:
        adata: AnnData (normalized, log1p CP10k)
        test_barcodes: list of cell barcodes for the test group
        ctrl_barcodes: list of cell barcodes for the control group
        readout_genes: list of gene names to measure
        test_label: label for the test group (used in plots)
        ctrl_label: label for the control group
        plot_type: 'heatmap', 'violin', 'both', or None (stats only)
        show_pval: if True, annotate p-values on plots
        save_path: file path to save figure. For 'both', two files are
            saved with '_heatmap' and '_violin' suffixes.
        title: plot title (auto-generated if None)
        show: if True, call plt.show(); otherwise plt.close()
        figsize: optional figure size override

    Returns:
        DataFrame with columns: readout_gene, log2FC, pvalue, cohens_d,
        n_test, n_ctrl, cpm_test, cpm_ctrl, test_label, ctrl_label
    """
    results = []

    for gene in readout_genes:
        resolved = _resolve_gene_name(gene, adata.var_names)
        if resolved is None:
            print(f"    WARNING: gene '{gene}' not found in var_names, skipping.")
            continue

        test_expr = _extract_expr(adata, test_barcodes, resolved)
        ctrl_expr = _extract_expr(adata, ctrl_barcodes, resolved)

        if len(test_expr) < 3 or len(ctrl_expr) < 3:
            results.append({
                'readout_gene': gene, 'log2FC': np.nan, 'pvalue': np.nan,
                'cohens_d': np.nan, 'n_test': len(test_expr),
                'n_ctrl': len(ctrl_expr), 'cpm_test': np.nan,
                'cpm_ctrl': np.nan, 'test_label': test_label,
                'ctrl_label': ctrl_label,
            })
            continue

        mean_test = float(np.mean(test_expr))
        mean_ctrl = float(np.mean(ctrl_expr))
        pseudo = 0.01
        lfc = float(np.log2((mean_test + pseudo) / (mean_ctrl + pseudo)))
        cpm_test = float(np.mean(np.expm1(test_expr) * 100))
        cpm_ctrl = float(np.mean(np.expm1(ctrl_expr) * 100))

        try:
            _, pval = stats.mannwhitneyu(test_expr, ctrl_expr, alternative='two-sided')
        except ValueError:
            pval = np.nan

        d = cohens_d(test_expr, ctrl_expr)

        results.append({
            'readout_gene': gene, 'log2FC': lfc, 'pvalue': pval,
            'cohens_d': d, 'n_test': len(test_expr),
            'n_ctrl': len(ctrl_expr), 'cpm_test': cpm_test,
            'cpm_ctrl': cpm_ctrl, 'test_label': test_label,
            'ctrl_label': ctrl_label,
        })

    result_df = pd.DataFrame(results)

    if plot_type is None or result_df.empty:
        return result_df

    auto_title = title or f"{test_label} vs {ctrl_label}"

    if plot_type in ('heatmap', 'both'):
        _path = save_path
        if save_path and plot_type == 'both':
            base, ext = os.path.splitext(save_path)
            _path = f"{base}_heatmap{ext or '.png'}"
        _plot_compare_heatmap(result_df, auto_title, show_pval, _path, show, figsize)

    if plot_type in ('violin', 'both'):
        _path = save_path
        if save_path and plot_type == 'both':
            base, ext = os.path.splitext(save_path)
            _path = f"{base}_violin{ext or '.png'}"
        _plot_compare_violin(
            adata, test_barcodes, ctrl_barcodes, readout_genes,
            result_df, test_label, ctrl_label, auto_title,
            show_pval, _path, show, figsize,
        )

    return result_df


def _plot_compare_heatmap(result_df, title, show_pval, save_path, show, figsize):
    """Internal: single-column heatmap for compare_expression."""
    valid = result_df.dropna(subset=['log2FC'])
    if valid.empty:
        return

    genes = valid['readout_gene'].tolist()
    lfcs = valid['log2FC'].values
    n_genes = len(genes)

    fig_w, fig_h = figsize or (4, max(2, n_genes * 0.8))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(lfcs.reshape(-1, 1), aspect='auto', cmap='RdBu_r',
                   vmin=-3, vmax=3)

    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(genes, fontsize=10)
    ax.set_xticks([0])
    ax.set_xticklabels([valid['test_label'].iloc[0]], fontsize=10)
    ax.set_title(title, fontsize=12)

    cbar = plt.colorbar(im, ax=ax, shrink=0.6)
    cbar.set_label('log₂FC')

    max_abs = max(np.nanmax(np.abs(lfcs)), 0.01)
    for i, row in valid.iterrows():
        idx = genes.index(row['readout_gene'])
        val = row['log2FC']
        color = 'white' if abs(val) > max_abs * 0.5 else 'black'
        txt = f"lfc: {val:.2f}\nn: {int(row['n_test'])}"
        if show_pval and not np.isnan(row['pvalue']):
            txt += f"\np: {row['pvalue']:.1e}"
        ax.text(0, idx, txt, ha='center', va='center', color=color, fontsize=8)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"    Saved heatmap: {save_path}")
    if show:
        plt.show()
    else:
        plt.close(fig)


def _format_pval(p):
    """Format p-value for plot annotation."""
    if p < 0.001:
        return f"p={p:.1e}"
    elif p < 0.01:
        return f"p={p:.3f}"
    elif p < 0.05:
        return f"p={p:.2f}"
    else:
        return f"p={p:.2f} (ns)"


def _plot_compare_violin(
    adata, test_barcodes, ctrl_barcodes, readout_genes,
    result_df, test_label, ctrl_label, title,
    show_pval, save_path, show, figsize,
):
    """Internal: violin plots for compare_expression."""
    import seaborn as sns

    n_genes = len(readout_genes)
    fig_w, fig_h = figsize or (max(4, n_genes * 2.5), 4)
    fig, axes = plt.subplots(1, n_genes, figsize=(fig_w, fig_h), squeeze=False)
    axes = axes.flatten()

    for i, gene in enumerate(readout_genes):
        ax = axes[i]
        resolved = _resolve_gene_name(gene, adata.var_names)
        if resolved is None:
            ax.set_title(f"{gene}\n(not found)", fontsize=9)
            ax.set_visible(False)
            continue

        test_expr = _extract_expr(adata, test_barcodes, resolved)
        ctrl_expr = _extract_expr(adata, ctrl_barcodes, resolved)

        plot_data = pd.DataFrame({
            'Expression': np.concatenate([ctrl_expr, test_expr]),
            'Group': [ctrl_label] * len(ctrl_expr) + [test_label] * len(test_expr),
        })

        sns.violinplot(
            data=plot_data, x='Group', y='Expression', ax=ax,
            palette=['#b0bec5', '#ef5350'], inner='box', cut=0,
            linewidth=0.8,
        )
        ax.set_title(gene, fontsize=11, fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel('log(CP10K+1)' if i == 0 else '')
        ax.tick_params(axis='x', rotation=30)

        # P-value bracket
        if show_pval:
            row = result_df[result_df['readout_gene'] == gene]
            if len(row) > 0 and not np.isnan(row.iloc[0]['pvalue']):
                pval = row.iloc[0]['pvalue']
                y_max = plot_data['Expression'].max()
                bracket_h = y_max * 0.05
                y_bar = y_max + bracket_h * 2

                ax.plot([0, 0, 1, 1], [y_bar - bracket_h, y_bar, y_bar, y_bar - bracket_h],
                        lw=1.0, c='black')
                ax.text(0.5, y_bar + bracket_h * 0.5, _format_pval(pval),
                        ha='center', va='bottom', fontsize=8)
                ax.set_ylim(top=y_bar + bracket_h * 4)

    fig.suptitle(title, fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"    Saved violin: {save_path}")
    if show:
        plt.show()
    else:
        plt.close(fig)


def compute_custom_combo_de(adata, cell_df, target_combos, readout_genes,
                            ctrl_barcodes, sample_label):
    """Compute DE stats for specified target combos across multiple readout genes.

    Args:
        adata: AnnData (normalized, log1p CP10k)
        cell_df: DataFrame of classified cells (e.g. perfect singlets)
        target_combos: list of target_combo strings (e.g. ['CCM2_TLNRD1', 'CCM2_MAPK7'])
        readout_genes: list of gene names to measure (e.g. ['CCM2', 'TLNRD1', 'MAPK7'])
        ctrl_barcodes: list of ctrl_ctrl cell barcodes
        sample_label: 'Control' or 'Treatment'

    Returns:
        gene_level_df: DataFrame with gene-combo level results
        guide_level_df: DataFrame with individual guide-pair level results
    """
    gene_results = []
    guide_results = []

    # Resolve gene names
    resolved = {}
    for g in readout_genes:
        rg = _resolve_gene_name(g, adata.var_names)
        if rg is None:
            print(f"    WARNING: readout gene '{g}' not found in var_names, skipping.")
        else:
            resolved[g] = rg

    if not resolved:
        return pd.DataFrame(), pd.DataFrame()

    # Get ctrl expression for each readout gene
    ctrl_expr = {}
    ctrl_cpm = {}
    for g_label, g_resolved in resolved.items():
        expr = _extract_expr(adata, ctrl_barcodes, g_resolved)
        ctrl_expr[g_label] = expr
        ctrl_cpm[g_label] = float(np.mean(np.expm1(expr) * 100)) if len(expr) > 0 else np.nan

    for tc in target_combos:
        tc_cells = cell_df[cell_df['target_combo'] == tc]
        if len(tc_cells) == 0:
            print(f"    WARNING: No cells found for combo '{tc}' in {sample_label}")
            continue

        # ---- Gene-combo level (all cells with this target_combo) ----
        combo_bcs = tc_cells.index.tolist()
        for g_label, g_resolved in resolved.items():
            kd_expr = _extract_expr(adata, combo_bcs, g_resolved)
            c_expr = ctrl_expr[g_label]

            if len(kd_expr) < 3 or len(c_expr) < 3:
                gene_results.append({
                    'target_combo': tc, 'readout_gene': g_label,
                    'sample': sample_label, 'level': 'gene_combo',
                    'guide_pair': 'all', 'n_kd': len(kd_expr), 'n_ctrl': len(c_expr),
                    'cpm_kd': np.nan, 'cpm_ctrl': ctrl_cpm[g_label],
                    'log2FC': np.nan, 'pvalue': np.nan, 'cohens_d': np.nan,
                })
                continue

            mean_kd_log = float(np.mean(kd_expr))
            mean_ctrl_log = float(np.mean(c_expr))
            pseudo = 0.01
            lfc = float(np.log2((mean_kd_log + pseudo) / (mean_ctrl_log + pseudo)))
            cpm_kd_val = float(np.mean(np.expm1(kd_expr) * 100))

            try:
                _, pval = stats.mannwhitneyu(kd_expr, c_expr, alternative='two-sided')
            except ValueError:
                pval = np.nan

            d = cohens_d(kd_expr, c_expr)

            gene_results.append({
                'target_combo': tc, 'readout_gene': g_label,
                'sample': sample_label, 'level': 'gene_combo',
                'guide_pair': 'all', 'n_kd': len(kd_expr), 'n_ctrl': len(c_expr),
                'cpm_kd': cpm_kd_val, 'cpm_ctrl': ctrl_cpm[g_label],
                'log2FC': lfc, 'pvalue': pval, 'cohens_d': d,
            })

        # ---- Guide-pair level ----
        # Build guide-pair identifiers from pos1_grna + pos2_grna
        if 'pos1_grna' in tc_cells.columns and 'pos2_grna' in tc_cells.columns:
            tc_cells = tc_cells.copy()
            tc_cells['guide_pair'] = tc_cells['pos1_grna'].astype(str).fillna('') + ' + ' + tc_cells['pos2_grna'].astype(str).fillna('')
            for gp in tc_cells['guide_pair'].unique():
                gp_cells = tc_cells[tc_cells['guide_pair'] == gp]
                gp_bcs = gp_cells.index.tolist()

                for g_label, g_resolved in resolved.items():
                    kd_expr = _extract_expr(adata, gp_bcs, g_resolved)
                    c_expr = ctrl_expr[g_label]

                    if len(kd_expr) < 10 or len(c_expr) < 3:
                        continue  # Skip guide pairs with too few cells

                    mean_kd_log = float(np.mean(kd_expr))
                    mean_ctrl_log = float(np.mean(c_expr))
                    pseudo = 0.01
                    lfc = float(np.log2((mean_kd_log + pseudo) / (mean_ctrl_log + pseudo)))
                    cpm_kd_val = float(np.mean(np.expm1(kd_expr) * 100))

                    try:
                        _, pval = stats.mannwhitneyu(kd_expr, c_expr, alternative='two-sided')
                    except ValueError:
                        pval = np.nan

                    d = cohens_d(kd_expr, c_expr)

                    guide_results.append({
                        'target_combo': tc, 'readout_gene': g_label,
                        'sample': sample_label, 'level': 'guide_pair',
                        'guide_pair': gp, 'n_kd': len(kd_expr), 'n_ctrl': len(c_expr),
                        'cpm_kd': cpm_kd_val, 'cpm_ctrl': ctrl_cpm[g_label],
                        'log2FC': lfc, 'pvalue': pval, 'cohens_d': d,
                    })

    return pd.DataFrame(gene_results), pd.DataFrame(guide_results)


def plot_custom_combo_heatmap(result_df, save_path, title='Custom Combo KD',
                              show=False):
    """Plot heatmap: rows=readout genes, columns=target combos.

    Args:
        result_df: output from compute_custom_combo_de (gene- or guide-level)
        save_path: path to save figure
        title: plot title
        show: if True, call plt.show() instead of plt.close()
    """
    if result_df.empty:
        return

    for sample in result_df['sample'].unique():
        sdf = result_df[result_df['sample'] == sample]

        pivot_fc = sdf.pivot_table(
            index='readout_gene', columns='target_combo',
            values='log2FC', aggfunc='first'
        )
        pivot_n = sdf.pivot_table(
            index='readout_gene', columns='target_combo',
            values='n_kd', aggfunc='first'
        )
        pivot_cpm = sdf.pivot_table(
            index='readout_gene', columns='target_combo',
            values='cpm_kd', aggfunc='first'
        )
        # Get ctrl CPM per readout gene (same across combos)
        ctrl_cpm_map = sdf.groupby('readout_gene')['cpm_ctrl'].first()

        if pivot_fc.empty:
            continue

        # Sort axes
        pivot_fc = pivot_fc.reindex(sorted(pivot_fc.index), axis=0)
        pivot_fc = pivot_fc.reindex(sorted(pivot_fc.columns), axis=1)
        pivot_n = pivot_n.reindex(pivot_fc.index, axis=0).reindex(pivot_fc.columns, axis=1)
        pivot_cpm = pivot_cpm.reindex(pivot_fc.index, axis=0).reindex(pivot_fc.columns, axis=1)

        n_rows, n_cols = pivot_fc.shape
        fig, ax = plt.subplots(figsize=(max(5, n_cols * 2.0), max(3, n_rows * 1.2)))

        im = ax.imshow(pivot_fc.values, aspect='auto', cmap='RdBu_r',
                       vmin=-3, vmax=3)

        # Y-axis: readout gene with ctrl CPM
        y_labels = []
        for gene in pivot_fc.index:
            c_cpm = ctrl_cpm_map.get(gene, np.nan)
            lbl = f"{gene}\n(Ctrl CPM={c_cpm:.1f})" if not np.isnan(c_cpm) else gene
            y_labels.append(lbl)

        ax.set_yticks(range(n_rows))
        ax.set_yticklabels(y_labels, fontsize=10)
        ax.set_xticks(range(n_cols))
        ax.set_xticklabels(pivot_fc.columns, rotation=30, ha='right',
                           rotation_mode='anchor', fontsize=10)
        ax.set_title(f'{title} — {sample}', fontsize=13)
        ax.set_xlabel('Target Combo')
        ax.set_ylabel('Readout Gene')

        cbar = plt.colorbar(im, ax=ax, shrink=0.7)
        cbar.set_label('log₂FC vs ctrl_ctrl')

        # Annotate
        max_abs = np.nanmax(np.abs(pivot_fc.values)) if not np.all(np.isnan(pivot_fc.values)) else 1
        for i in range(n_rows):
            for j in range(n_cols):
                val = pivot_fc.values[i, j]
                n_c = pivot_n.values[i, j]
                cpm_v = pivot_cpm.values[i, j]
                if np.isnan(val):
                    continue
                color = 'white' if abs(val) > max_abs * 0.5 else 'black'
                txt = f"lfc: {val:.2f}\ncpm: {cpm_v:.1f}\nn: {int(n_c)}"
                ax.text(j, i, txt, ha='center', va='center', color=color, fontsize=8)

        plt.tight_layout()
        suffix = sample.lower().replace(' ', '_')
        out_path = save_path.replace('.png', f'_{suffix}.png')
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
        print(f"    Saved: {out_path}")
        if show:
            plt.show()
        else:
            plt.close(fig)

# NOTE: estimate_background_noise has been moved to qc_q2_recombination.py
# Import it from there if needed:
#   from scripts.qc_q2_recombination import estimate_background_noise

def classify_inferred_sources(adata, assignment_df, umi_df, method="percentile", percentile=50, metric="total_umi", cell_category="perfect_singlet"):
    """Classify inferred singlets into Source 1 (Dropout) vs Source 2 (Recombinant).

    Uses noise threshold from Q2's estimate_background_noise to mask ambient UMIs,
    then checks whether the missing position has any remaining signal.

    Args:
        adata: AnnData object
        assignment_df: Cell assignment DataFrame
        umi_df: UMI counts DataFrame
        method: "mean", "median", or "percentile" for setting the noise threshold.
        percentile: Numeric value (e.g., 50) if method="percentile".
        metric: Metric to use for background noise estimation
        cell_category: Cell category for noise estimation

    Returns:
        DataFrame of inferred singlets with 'source' column ('dropout' or 'recombinant').
        Returns None if no inferred singlets found.
    """
    from scripts.qc_q2_recombination import estimate_background_noise

    print(f"\n  --- Q3c: Inferred Singlet Source Classification ---")

    # 1. Estimate global noise and compute threshold
    per_cell_noise = estimate_background_noise(adata, assignment_df, umi_df, metric=metric, cell_category=cell_category)

    if len(per_cell_noise) == 0:
        print("    Cannot compute noise threshold. Skipping Q3c.")
        return None

    if method == "mean":
        noise_threshold = per_cell_noise.mean()
    elif method == "median":
        noise_threshold = per_cell_noise.median()
    elif method == "percentile":
        noise_threshold = np.percentile(per_cell_noise, percentile)
    else:
        print(f"    WARNING: Unknown method '{method}', defaulting to median.")
        noise_threshold = per_cell_noise.median()

    print(f"\n    Selected Noise Threshold ({method}{' '+str(percentile) if method=='percentile' else ''}): {noise_threshold:.2f} UMIs")

    # 2. Filter valid inferred singlets
    valid_bcs = set(adata.obs.index)
    ass_df = assignment_df[assignment_df['bc_wells'].isin(valid_bcs)]
    inf_singlets = ass_df[(ass_df['detection_status'] == 'inferred') &
                          (ass_df['combo_type'] == 'singlet')].copy()

    if inf_singlets.empty:
        print("    No inferred singlets found in valid GEX cells.")
        return None

    print(f"    Analyzing {len(inf_singlets):,} valid inferred singlets.\n")

    # 3. Identify missing position (pos1 missing -> 1, pos2 missing -> 2)
    inf_singlets['missing_pos'] = np.where(inf_singlets['pos1_umi'] == 0, 1,
                                    np.where(inf_singlets['pos2_umi'] == 0, 2, np.nan))

    # 4. Map back to umi_df for inferred singlets
    inf_bcs = set(inf_singlets['bc_wells'])
    sub_umi = umi_df[umi_df['bc_wells'].isin(inf_bcs)].copy()

    # Heuristically define expected position of each observed read
    sub_umi['grna_pos'] = np.where(sub_umi['grna_name'].str.contains('_A_', na=False), 1,
                             np.where(sub_umi['grna_name'].str.contains('_B_', na=False), 2, np.nan))

    cell_to_sample = adata.obs['sample'].to_dict()
    inf_singlets['sample'] = inf_singlets['bc_wells'].map(cell_to_sample)

    # 5. Classify each cell: post-mask sum == 0 → dropout, > 0 → recombinant
    source_labels = pd.Series('dropout', index=inf_singlets.index)

    for target_pos in [1, 2]:
        mask = inf_singlets['missing_pos'] == target_pos
        if not mask.any():
            continue

        missing_bcs = inf_singlets.loc[mask, 'bc_wells']
        pos_reads = sub_umi[sub_umi['grna_pos'] == target_pos]

        post_reads = pos_reads[pos_reads['umi_count'] > noise_threshold]
        post_sum = post_reads.groupby('bc_wells')['umi_count'].sum()
        post_sum = post_sum.reindex(missing_bcs, fill_value=0)

        # Cells with residual signal → recombinant
        recomb_bcs = set(post_sum[post_sum > 0].index)
        source_labels.loc[mask & inf_singlets['bc_wells'].isin(recomb_bcs)] = 'recombinant'

    inf_singlets['source'] = source_labels.values

    # 6. Print summary
    for sample_name in ['Control', 'Treatment']:
        samp = inf_singlets[inf_singlets['sample'] == sample_name]
        if samp.empty:
            continue
        n = len(samp)
        n_drop = (samp['source'] == 'dropout').sum()
        n_recomb = (samp['source'] == 'recombinant').sum()
        print(f"    {sample_name} (n={n:,}):")
        print(f"      Source 1 (Dropout):     {n_drop:>6,} ({n_drop/n*100:5.1f}%)")
        print(f"      Source 2 (Recombinant): {n_recomb:>6,} ({n_recomb/n*100:5.1f}%)")

    return inf_singlets


def main():
    parser = argparse.ArgumentParser(
        description="QC Q3: KD Efficacy + Inferred Source Classification"
    )
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    output_dir = cfg['paths']['output_dir']

    print(f"\n{'='*60}")
    print(f"  QC Q3: KD Efficacy + Inferred Source Classification")
    print(f"{'='*60}")

    # --- 1. Load AnnData ---
    print("\n  1. Loading combined AnnData...")
    adata_path = os.path.join(output_dir, 'combined_adata.h5ad')
    adata = sc.read_h5ad(adata_path)
    print(f"     Shape: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    print(f"     obs columns: {list(adata.obs.columns)}")

    # --- 2. Normalize ---
    print("\n  2. Normalizing GEX...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    print("     Done (total-count 10K + log1p)")

    # --- 3. Classify cells ---
    print("\n  3. Classifying cells into test sets...")
    cell_sets = classify_cells(adata)

    # --- 4. Compute KD metrics per test set × sample ---
    print("\n  4. Computing KD metrics...")

    lib_to_sample = {}
    for lib in cfg['libraries']:
        lib_to_sample[lib['sublibrary_name']] = lib['sample_name']

    out_dir = os.path.join(output_dir, 'qc_analysis')
    os.makedirs(out_dir, exist_ok=True)

    all_results = []

    # Extract ctrl_ctrl reference cells ONCE per sample from perfect singlets
    perfect_singlets = cell_sets.get('perfect_singlet', (pd.DataFrame(), 10, False))[0]
    ctrl_barcodes_per_sample = {}
    for sample_name in ['Control', 'Treatment']:
        if 'sample' in perfect_singlets.columns:
            ps_sample = perfect_singlets[perfect_singlets['sample'] == sample_name]
        elif 'sublibrary' in perfect_singlets.columns:
            ps_sample = perfect_singlets[perfect_singlets['sublibrary'].map(lib_to_sample) == sample_name]
        else:
            ps_sample = pd.DataFrame()
        ctrl = get_ctrl_ctrl_cells(ps_sample)
        ctrl_barcodes_per_sample[sample_name] = ctrl.index.tolist()
        print(f"  ctrl_ctrl reference ({sample_name}): {len(ctrl):,} cells (from perfect singlets)")

    for set_name, (cell_df, mc, gbg) in cell_sets.items():
        print(f"\n  --- Test Set: {set_name} (min_cells={mc}, group_by_gene={gbg}) ---")

        if len(cell_df) == 0:
            print(f"    No cells in this set, skipping.")
            continue

        for sample_name in ['Control', 'Treatment']:
            sample_mask = cell_df['sample'] == sample_name
            if 'sample' not in cell_df.columns:
                sample_mask = cell_df['sublibrary'].map(lib_to_sample) == sample_name

            sample_cells = cell_df[sample_mask]

            if len(sample_cells) == 0:
                print(f"    No {sample_name} cells in {set_name}, skipping.")
                continue

            print(f"\n    Sample: {sample_name} ({len(sample_cells):,} cells)")

            ctrl_barcodes = ctrl_barcodes_per_sample[sample_name]
            print(f"    ctrl_ctrl reference cells: {len(ctrl_barcodes):,} (shared)")

            if len(ctrl_barcodes) < 10:
                print(f"    WARNING: Very few ctrl_ctrl cells ({len(ctrl_barcodes)}). "
                      f"Results may be unreliable.")

            guide_groups = get_grna_ctrl_cells(sample_cells, group_by_gene=gbg)
            print(f"    Guide groups (X_ctrl / ctrl_Y): {len(guide_groups)}"
                  f" {'(gene-level pooled)' if gbg else '(per-guide)'}")

            if not guide_groups:
                print(f"    No gRNA_ctrl groups found, skipping.")
                continue

            for g in guide_groups:
                if g['n_cells'] < mc:
                    print(f"      WARNING: {g['guide_name']} → {g['target_gene']} "
                          f"has only {g['n_cells']} cells (threshold={mc})")

            metrics_df = compute_kd_metrics(
                adata, guide_groups, ctrl_barcodes, sample_name, min_cells=mc
            )

            if not metrics_df.empty:
                metrics_df['test_set'] = set_name
                all_results.append(metrics_df)
                print(f"    Computed metrics for {len(metrics_df)} guide-target pairs")

    # --- 5. Combine and save ---
    print(f"\n  5. Saving results...")
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)

        csv_path = os.path.join(out_dir, 'q3_guide_efficacy.csv')
        combined.to_csv(csv_path, index=False)
        print(f"    Saved: {csv_path}")
        print(f"    Total rows: {len(combined)}")

        print(f"\n  Summary by test set:")
        for ts in combined['test_set'].unique():
            ts_df = combined[combined['test_set'] == ts]
            sig = ts_df[(ts_df['pvalue'] < 0.05) & (ts_df['log2FC'] < -0.5)]
            print(f"    {ts}:")
            print(f"      Total guide-target pairs: {len(ts_df)}")
            print(f"      Significant KD (p<0.05, log2FC<-0.5): {len(sig)}")
            if not ts_df.empty:
                print(f"      Median log2FC: {ts_df['log2FC'].median():.3f}")
                print(f"      Median Cohen's d: {ts_df['cohens_d'].median():.3f}")

        # --- 6. Plots ---
        print(f"\n  6. Generating plots...")
        for set_name in combined['test_set'].unique():
            set_df = combined[combined['test_set'] == set_name]
            plot_heatmap(
                set_df,
                os.path.join(out_dir, f'q3_heatmap_{set_name}.png'),
                title_suffix=f' ({set_name})'
            )
            plot_volcano(
                set_df,
                os.path.join(out_dir, f'q3_volcano_{set_name}.png'),
                title_suffix=f' ({set_name})'
            )
    else:
        print("    WARNING: No results computed. Check cell classification.")

    # --- 7. Q3c: Inferred Source Classification + Source 1 DE ---
    print(f"\n  7. Q3c: Inferred Source Classification...")
    assign_path = os.path.join(output_dir, 'grna_assignment', 'grna_cell_assignment.csv')
    umi_path = os.path.join(output_dir, 'grna_assignment', 'grna_umi_counts.csv')
    if os.path.exists(assign_path) and os.path.exists(umi_path):
        assignment_df = pd.read_csv(assign_path)
        umi_df = pd.read_csv(umi_path)
        qc_cfg = cfg.get('qc', {})
        metric = qc_cfg.get('q1c_noise_metric', 'total_umi')
        method = qc_cfg.get('q1c_noise_method', 'percentile')
        pct = qc_cfg.get('q1c_noise_percentile', 50)
        category = qc_cfg.get('q1c_noise_category', 'perfect_singlet')

        inf_classified = classify_inferred_sources(
            adata, assignment_df, umi_df,
            method=method, percentile=pct, metric=metric, cell_category=category
        )

        # --- 7b. Repeat DE on Source 1 (dropout) inferred singlets ---
        if inf_classified is not None and 'source' in inf_classified.columns:
            source1 = inf_classified[inf_classified['source'] == 'dropout'].copy()
            if not source1.empty:
                print(f"\n  7b. Source 1 DE validation ({len(source1):,} dropout inferred singlets)")
                source1_results = []

                for sample_name in ['Control', 'Treatment']:
                    s1_sample = source1[source1['sample'] == sample_name]
                    if s1_sample.empty:
                        continue

                    ctrl_barcodes = ctrl_barcodes_per_sample.get(sample_name, [])
                    if len(ctrl_barcodes) < 10:
                        continue

                    guide_groups = get_grna_ctrl_cells(s1_sample, group_by_gene=True)
                    if not guide_groups:
                        continue

                    metrics_df = compute_kd_metrics(
                        adata, guide_groups, ctrl_barcodes, sample_name, min_cells=5
                    )
                    if not metrics_df.empty:
                        metrics_df['test_set'] = 'source1_inferred'
                        source1_results.append(metrics_df)

                if source1_results:
                    s1_combined = pd.concat(source1_results, ignore_index=True)
                    s1_path = os.path.join(out_dir, 'q3_source1_inferred_efficacy.csv')
                    s1_combined.to_csv(s1_path, index=False)
                    print(f"    Saved Source 1 DE results: {s1_path}")
                    print(f"    Total Source 1 guide-target pairs: {len(s1_combined)}")

                    plot_heatmap(
                        s1_combined,
                        os.path.join(out_dir, 'q3_heatmap_source1_inferred.png'),
                        title_suffix=' (Source 1 Inferred)'
                    )
            else:
                print("    No Source 1 (dropout) inferred singlets found for DE validation.")
    else:
        print(f"    WARNING: Missing assignment or UMI counts CSV. Skipping Q3c.")

    print(f"\n{'='*60}")
    print(f"  QC Q3 Complete")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()

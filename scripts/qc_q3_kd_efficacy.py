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

    # --- Sets 2 & 3: Inferred singlets (pos1_missing and pos2_missing) ---
    inferred_mask = (
        (assigned['combo_type'] == 'singlet') &
        (assigned['detection_status'] == 'inferred')
    )
    inferred = assigned[inferred_mask].copy()

    if 'pos1_umi' in inferred.columns and 'pos2_umi' in inferred.columns:
        cell_sets['inferred_singlet_pos1_missing'] = inferred[inferred['pos1_umi'] == 0].copy()
        cell_sets['inferred_singlet_pos2_missing'] = inferred[inferred['pos2_umi'] == 0].copy()
    else:
        has_inf_p1 = 'inferred_pos1_grna' in inferred.columns
        has_inf_p2 = 'inferred_pos2_grna' in inferred.columns
        if has_inf_p1:
            p1_mask = inferred['inferred_pos1_grna'].fillna('') != ''
            cell_sets['inferred_singlet_pos1_missing'] = inferred[p1_mask].copy()
        else:
            cell_sets['inferred_singlet_pos1_missing'] = inferred.iloc[0:0].copy()
        if has_inf_p2:
            p2_mask = inferred['inferred_pos2_grna'].fillna('') != ''
            cell_sets['inferred_singlet_pos2_missing'] = inferred[p2_mask].copy()
        else:
            cell_sets['inferred_singlet_pos2_missing'] = inferred.iloc[0:0].copy()

    # --- Sets 4, 5, 6: Multiplets ---
    multi_mask = assigned['combo_type'] == 'multiplet'
    multi = assigned[multi_mask].copy()

    # Check for presence of 'complete' and 'inferred' in the potentially joined status strings
    has_complete = multi['detection_status'].str.contains('complete', na=False)
    has_inferred = multi['detection_status'].str.contains('inferred', na=False)

    mc_mask = has_complete & ~has_inferred
    mi_mask = has_inferred & ~has_complete
    mm_mask = has_complete & has_inferred

    cell_sets['multiplet_all_complete'] = multi[mc_mask].copy()
    cell_sets['multiplet_all_inferred'] = multi[mi_mask].copy()
    cell_sets['multiplet_mixed'] = multi[mm_mask].copy()

    # Catch remaining multiplets (e.g. unexpected missing status)
    classified_multi = mc_mask | mi_mask | mm_mask
    remaining_multi = ~classified_multi
    if remaining_multi.sum() > 0:
        cell_sets['multiplet_other'] = multi[remaining_multi].copy()

    # --- Set 7: Unassigned ---
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
    
    # Reconstruct combined inferred_singlet for legacy compatibility in Q3
    inferred_singlet = pd.concat([full['inferred_singlet_pos1_missing'], full['inferred_singlet_pos2_missing']])
    full['inferred_singlet'] = inferred_singlet

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

def filter_and_compute_inferred_gene_kd(adata, umi_df, ctrl_bcs_per_sample, backend='custom', exclude_combos=None, group_by_gene=True, min_cells=5):
    """
    Groups inferred singlets by gene, strict-filters cells that might have received 
    another guide for the 'missing' gene (UMI > 0), and computes KD metrics for both 
    the detected target and the inferred (missing) target against control cells.

    Args:
        adata: AnnData object.
        umi_df: DataFrame containing all UMI counts.
        ctrl_bcs_per_sample: Dict mapping sample_name -> list of control cell barcodes.
        backend: Which calculation backend to use ('custom' for fast scipy or 'scanpy').
        exclude_combos: List of gene combos (e.g. ['MAPK7_CCM2']) to exclude.
        group_by_gene: If False, computes on a per-guide level and falls back to genepool if < min_cells.
        min_cells: Minimum number of cells required per perturbation subset for DE testing.

    Returns:
        A combined DataFrame of KD results mapping the target combos to the DE outcome.
    """
    import pandas as pd
    import numpy as np

    # Build the inferred singlet dataframe directly from properly-filtered adata.obs
    obs_df = adata.obs
    inf = obs_df[(obs_df['detection_status'] == 'inferred') & 
                 (obs_df['combo_type'] == 'singlet')].copy()

    if inf.empty:
        return pd.DataFrame()

    inf['bc_wells'] = inf.index

    # 1. Determine the missing and detected genes and guides per cell
    p1_inf = inf.get('inferred_pos1_target', inf.get('pos1_target'))
    p2_inf = inf.get('inferred_pos2_target', inf.get('pos2_target'))
    p1_det = inf.get('pos1_target')
    p2_det = inf.get('pos2_target')

    p1_guide_inf = inf.get('inferred_pos1_grna', inf.get('pos1_grna'))
    p2_guide_inf = inf.get('inferred_pos2_grna', inf.get('pos2_grna'))
    p1_guide_det = inf.get('pos1_grna')
    p2_guide_det = inf.get('pos2_grna')

    inf['missing_gene'] = np.where(inf['pos1_umi'] == 0, p1_inf,
                                   np.where(inf['pos2_umi'] == 0, p2_inf, np.nan))
    inf['detected_gene'] = np.where(inf['pos1_umi'] == 0, p2_det,
                                    np.where(inf['pos2_umi'] == 0, p1_det, np.nan))
    
    inf['missing_guide'] = np.where(inf['pos1_umi'] == 0, p1_guide_inf,
                                    np.where(inf['pos2_umi'] == 0, p2_guide_inf, np.nan))
    inf['detected_guide'] = np.where(inf['pos1_umi'] == 0, p2_guide_det,
                                     np.where(inf['pos2_umi'] == 0, p1_guide_det, np.nan))
                                     
    # Preserve blank strings or pseudo-null 'unassigned' markers so viable partners aren't wiped
    inf['missing_gene'] = inf['missing_gene'].replace({'': 'unassigned'})
    inf['detected_gene'] = inf['detected_gene'].replace({'': 'unassigned'})
    inf['missing_guide'] = inf['missing_guide'].replace({'': 'unassigned'})
    inf['detected_guide'] = inf['detected_guide'].replace({'': 'unassigned'})
    
    inf = inf.dropna(subset=['missing_gene', 'detected_gene'])

    # Exclude confounding combos if provided
    if exclude_combos:
        exclude_set = set(exclude_combos)
        combo_str1 = inf['detected_gene'].astype(str) + '_' + inf['missing_gene'].astype(str)
        combo_str2 = inf['missing_gene'].astype(str) + '_' + inf['detected_gene'].astype(str)
        mask = (~combo_str1.isin(exclude_set)) & (~combo_str2.isin(exclude_set))
        inf = inf[mask]
        
    if inf.empty:
        return pd.DataFrame()

    # 2. Filter cells by checking umi_df for ANY expression of the 'missing_gene'
    cell_bcs = set(inf['bc_wells'])
    sub_umi = umi_df[umi_df['bc_wells'].isin(cell_bcs)].copy()
    
    missing_gene_map = inf.set_index('bc_wells')['missing_gene'].to_dict()
    sub_umi['expected_missing_gene'] = sub_umi['bc_wells'].map(missing_gene_map)
    
    sub_umi['is_missing_gene_guide'] = sub_umi.apply(
        lambda r: str(r['grna_name']).startswith(str(r['expected_missing_gene']) + '_'), axis=1
    )
    
    residual_expr = sub_umi[(sub_umi['is_missing_gene_guide']) & (sub_umi['umi_count'] > 0)]
    contaminated_bcs = set(residual_expr['bc_wells'])
    
    pure_inf = inf[~inf['bc_wells'].isin(contaminated_bcs)]
    
    print(f"Filter Inferred Singlets (Gene-level): {len(inf):,} total -> {len(pure_inf):,} pure subsets (0 UMIs for all contaminant guides)")
    
    # 3. Compute KD Metrics utilizing compute_kd_metrics (fast) or compute_custom_de (scanpy)
    results = []

    def compute_cpm(gene, group_bcs):
        query_g = gene if gene in adata.var_names else f"{gene}_hg38"
        if query_g not in adata.var_names: return np.nan
        expr = adata[adata.obs.index.isin(group_bcs), query_g].X
        if hasattr(expr, 'toarray'): expr = expr.toarray().flatten()
        else: expr = np.asarray(expr).flatten()
        return float(np.mean(np.expm1(expr) * 100))
        
    def is_control_target(g):
        if pd.isna(g): return True
        gs = str(g).lower()
        if gs in ['unassigned', '', 'none', 'nan']: return True
        if 'ctrl' in gs or 'safe' in gs or 'ntc' in gs: return True
        return False
        
    def generate_guide_groups(df, gene_col, guide_col, mode_name, pos_name):
        grp_list = []
        if group_by_gene:
            for gene, d in df.groupby(gene_col):
                if is_control_target(gene): continue
                bcs = list(set(d['bc_wells']).intersection(adata.obs.index))
                if len(bcs) < min_cells: continue
                grp_list.append({
                    'target_gene': gene,
                    'guide_name': mode_name,
                    'position': pos_name,
                    'target_combo': f"{gene}_pooled",
                    'cell_barcodes': bcs,
                    'actual_guide': f"{gene}_pooled"
                })
        else:
            pooled_bcs = {g: list(set(d['bc_wells']).intersection(adata.obs.index)) 
                          for g, d in df.groupby(gene_col) if not is_control_target(g)}
            for (gene, guide), d in df.groupby([gene_col, guide_col]):
                if is_control_target(gene): continue
                bcs = list(set(d['bc_wells']).intersection(adata.obs.index))
                if len(bcs) < min_cells:
                    fb_bcs = pooled_bcs.get(gene, [])
                    if len(fb_bcs) < min_cells: continue
                    if not any(x['actual_guide'] == f"{gene}_all (fallback)" for x in grp_list):
                        grp_list.append({
                            'target_gene': gene,
                            'guide_name': mode_name,
                            'position': pos_name,
                            'target_combo': f"{gene}_all (fallback)",
                            'cell_barcodes': fb_bcs,
                            'actual_guide': f"{gene}_all (fallback)"
                        })
                else:
                    grp_list.append({
                        'target_gene': gene,
                        'guide_name': mode_name,
                        'position': pos_name,
                        'target_combo': guide,
                        'cell_barcodes': bcs,
                        'actual_guide': guide
                    })
        return grp_list
    
    for sample, s_df in pure_inf.groupby('sample'):
        ctrl_bcs = ctrl_bcs_per_sample.get(sample, [])
        if not ctrl_bcs or len(ctrl_bcs) < 3:
            continue
            
        if backend == 'custom':
            guide_groups = []
            guide_groups.extend(generate_guide_groups(s_df, 'missing_gene', 'missing_guide', 'inferred_target_kd', 'missing'))
            guide_groups.extend(generate_guide_groups(s_df, 'detected_gene', 'detected_guide', 'detected_target_kd', 'detected'))
            
            if guide_groups:
                res_df = compute_kd_metrics(adata, guide_groups, ctrl_bcs, sample, min_cells=min_cells)
                if not res_df.empty:
                    # Rename metric columns to align natively with output schema expectations
                    res_df = res_df.rename(columns={
                        'guide_name': 'test_mode',
                        'target_gene': 'readout_gene',
                        'n_kd_cells': 'n_kd',
                        'n_ctrl_cells': 'n_ctrl'
                    })
                    if not group_by_gene:
                        ag_map = {g['target_combo']: g['actual_guide'] for g in guide_groups}
                        res_df['guide_name'] = res_df['target_combo'].map(ag_map)
                    
                    res_df['p_adj'] = np.nan
                    results.append(res_df)
                    
        else:
            all_grps = []
            all_grps.extend(generate_guide_groups(s_df, 'missing_gene', 'missing_guide', 'inferred_target_kd', 'missing'))
            all_grps.extend(generate_guide_groups(s_df, 'detected_gene', 'detected_guide', 'detected_target_kd', 'detected'))
                
            for g in all_grps:
                de_out = compute_custom_de(adata, g['cell_barcodes'], ctrl_bcs, readout_gene=g['target_gene'])
                if de_out.empty: continue
                
                row = {
                    'sample': sample,
                    'target_combo': g['target_combo'],
                    'test_mode': g['guide_name'],
                    'readout_gene': g['target_gene'],
                    'n_kd': len(g['cell_barcodes']),
                    'n_ctrl': len(ctrl_bcs),
                    'cpm_kd': compute_cpm(g['target_gene'], g['cell_barcodes']),
                    'cpm_ctrl': compute_cpm(g['target_gene'], ctrl_bcs),
                    'log2FC': de_out['logfoldchanges'].iloc[0],
                    'p_adj': de_out['pvals_adj'].iloc[0] if 'pvals_adj' in de_out.columns else np.nan,
                    'pvalue': de_out['pvals'].iloc[0],
                }
                if not group_by_gene:
                    row['guide_name'] = g['actual_guide']
                results.append(pd.DataFrame([row]))
                
    if not results:
        return pd.DataFrame()
        
    return pd.concat(results, ignore_index=True)


def compute_custom_de(adata, target_idx, ctrl_idx, readout_gene=None):
    """
    Perform differential expression (Wilcoxon) testing using Scanpy between an arbitrary 
    list of target cells and control cells.

    Args:
        adata: AnnData object containing the expression matrix.
        target_idx: List or index of cell barcodes for the target/treatment group.
        ctrl_idx: List or index of cell barcodes for the control group.
        readout_gene: (Optional) String or list of strings. If None, returns all genes tested.

    Returns:
        A pandas DataFrame of DE results with columns: names, scores, pvals, pvals_adj, logfoldchanges.
    """
    import scanpy as sc
    import pandas as pd

    # Filter to only cells that exist in the object
    target_cells = list(set(target_idx).intersection(adata.obs.index))
    ctrl_cells = list(set(ctrl_idx).intersection(adata.obs.index))

    if len(target_cells) < 3 or len(ctrl_cells) < 3:
        print(f"WARNING: Insufficient cells for DE test (Target: {len(target_cells)}, Ctrl: {len(ctrl_cells)}). Minimum 3 per group required.")
        return pd.DataFrame()

    # Subset the object for testing
    bdata = adata[target_cells + ctrl_cells, :].copy()
    
    # Assign groups
    bdata.obs['de_group'] = 'ctrl'
    bdata.obs.loc[target_cells, 'de_group'] = 'target'
    
    # Run scanpy DE
    sc.tl.rank_genes_groups(bdata, groupby='de_group', groups=['target'], reference='ctrl', method='wilcoxon')
    
    # Extract results to dataframe
    result_df = sc.get.rank_genes_groups_df(bdata, group='target')
    
    # Filter by readout gene if specified
    if readout_gene is not None:
        if isinstance(readout_gene, str):
            readout_gene = [readout_gene]
            
        valid_genes = []
        for g in readout_gene:
            if g in result_df['names'].values:
                valid_genes.append(g)
            elif f"{g}_hg38" in result_df['names'].values:
                valid_genes.append(f"{g}_hg38")
            else:
                matches = result_df[result_df['names'].str.startswith(f"{g}_")]['names'].tolist()
                if matches:
                    valid_genes.extend(matches)
        
        if len(valid_genes) > 0:
            result_df = result_df[result_df['names'].isin(valid_genes)]
        else:
            print(f"WARNING: None of the readout genes {readout_gene} were found in the dataset DE results.")
            result_df = result_df.iloc[0:0]

    return result_df


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




def plot_de_heatmap(result_df, save_path, title_suffix='', 
                    row_key="guide_name", col_key="test_mode", color_key="log2FC", 
                    split_by="sample", annotation_key=["cpm_kd", "pvalue", "n_kd_cells"]):
    """
    Generic DE heatmap plotting function.
    
    Args:
        result_df: DataFrame containing the data to plot (from compute_kd_metrics).
        save_path: Base path to save the heatmap(s).
        title_suffix: Optional string to append to the title.
        row_key: DataFrame column name to use for the y-axis.
        col_key: DataFrame column name to use for the x-axis.
        color_key: DataFrame column name for the heatmap color intensity.
        split_by: Column name to split the data into multiple plots. If None, plots all data.
        annotation_key: List of column names to annotate each cell with.
    """
    if result_df.empty:
        return

    valid = result_df.dropna(subset=[color_key])
    if valid.empty:
        return

    splits = valid[split_by].dropna().unique() if split_by else [None]

    for split_val in sorted(splits):
        if split_val is not None:
            sample_df = valid[valid[split_by] == split_val].copy()
            title = f'DE Heatmap — {split_val}{title_suffix}'
            suffix = f"_{str(split_val).lower().replace(' ', '_')}"
        else:
            sample_df = valid.copy()
            title = f'DE Heatmap{title_suffix}'
            suffix = ''

        # 3. Throw a warning if row_key and col_key do not generate unique entries
        counts = sample_df.groupby([row_key, col_key]).size()
        if (counts > 1).any():
            print(f"    ⚠ WARNING: The combination of row_key '{row_key}' and col_key '{col_key}' "
                  f"(split_by '{split_by}') does not generate unique entries. Averaging duplicates.")

        # Pivot matrix for colors
        pivot_color = sample_df.pivot_table(
            index=row_key, columns=col_key, values=color_key, aggfunc='mean'
        )
        
        # Pivot matrices for annotations
        annots = {}
        if annotation_key:
            for k in annotation_key:
                if k in sample_df.columns:
                    annots[k] = sample_df.pivot_table(
                        index=row_key, columns=col_key, values=k, aggfunc='mean'
                    )
        
        if pivot_color.empty or np.all(np.isnan(pivot_color.values)):
            continue

        # Sort indices: preserve Categorical ordering if present, else alphabetical
        if isinstance(sample_df[row_key].dtype, pd.CategoricalDtype):
            row_order = [x for x in sample_df[row_key].cat.categories if x in pivot_color.index]
        else:
            row_order = sorted(pivot_color.index)
            
        if isinstance(sample_df[col_key].dtype, pd.CategoricalDtype):
            col_order = [x for x in sample_df[col_key].cat.categories if x in pivot_color.columns]
        else:
            col_order = sorted(pivot_color.columns)

        pivot_color = pivot_color.reindex(index=row_order, columns=col_order)
        for k in annots:
            annots[k] = annots[k].reindex(index=row_order, columns=col_order)

        n_rows, n_cols = pivot_color.shape
        row_height = 1.0 + 0.35 * len(annotation_key) if annotation_key else 1.2
        fig, ax = plt.subplots(figsize=(max(8, n_cols * 1.5), max(4, n_rows * row_height)))

        im = ax.imshow(pivot_color.values, aspect='auto', cmap='RdBu_r', vmin=-3, vmax=3)

        ax.set_xticks(range(n_cols))
        ax.set_xticklabels(pivot_color.columns, rotation=45, ha='right', rotation_mode='anchor', fontsize=9)
        ax.set_yticks(range(n_rows))
        ax.set_yticklabels(pivot_color.index, fontsize=9)
        ax.set_title(title, fontsize=12)
        ax.set_xlabel(col_key)
        ax.set_ylabel(row_key)

        cbar = plt.colorbar(im, ax=ax, shrink=0.6)
        cbar.set_label(color_key)

        max_abs_val = np.nanmax(np.abs(pivot_color.values)) if not np.all(np.isnan(pivot_color.values)) else 1.0
        font_size = 8 if len(annotation_key) <= 1 else 6

        for i in range(n_rows):
            for j in range(n_cols):
                val = pivot_color.values[i, j]
                if np.isnan(val):
                    ax.text(j, i, 'n/a', ha='center', va='center', color='grey', fontsize=font_size)
                    continue

                text_color = 'white' if abs(val) > max_abs_val * 0.5 else 'black'

                # Format annotation
                annot_text = f"{color_key}: {val:.2f}"
                for k in annotation_key:
                    if k in annots:
                        v = annots[k].values[i, j]
                        if not np.isnan(v):
                            if k == 'pvalue':
                                annot_text += f"\n{k}: {v:.1e}"
                            elif 'n_' in k or (isinstance(v, (int, float, np.number)) and float(v).is_integer()):
                                annot_text += f"\n{k}: {int(v)}"
                            else:
                                annot_text += f"\n{k}: {v:.1f}"

                ax.text(j, i, annot_text, ha='center', va='center', color=text_color, fontsize=font_size)

        plt.tight_layout()
        out_path = save_path.replace('.png', f'{suffix}.png')
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"    Saved generic heatmap: {out_path}")



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
        combined_df: DataFrame with columns log2FC, pvalue, test_mode, guide_name, sample
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
                label = row.get('guide_name', row.get('target_gene', ''))
                ax.annotate(label, (row['log2FC'], row['neg_log10p']),
                            fontsize=6, alpha=0.7, ha='center', va='bottom',
                            color=style['color'])

        # Connect the same guide across groups with faint lines
        guides_in_data = sample_df['guide_name'].unique()
        for guide in guides_in_data:
            guide_rows = sample_df[sample_df['guide_name'] == guide]
            if len(guide_rows) > 1:
                xs = guide_rows['log2FC'].values
                ys = guide_rows['neg_log10p'].values
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

    gene_df = pd.DataFrame(gene_results)
    guide_df = pd.DataFrame(guide_results)
    
    # Convert 'target_combo' and 'readout_gene' to ordered Categoricals to preserve input order
    for df in [gene_df, guide_df]:
        if not df.empty:
            df['target_combo'] = pd.Categorical(df['target_combo'], categories=target_combos, ordered=True)
            df['readout_gene'] = pd.Categorical(df['readout_gene'], categories=list(resolved.keys()), ordered=True)
            
    return gene_df, guide_df




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
    from scripts.guide_design import legacy_model_eligible
    if 'guide_unit_id' in umi_df.columns and not legacy_model_eligible(umi_df):
        warnings.warn('Q3c unavailable: design-aware missing-partner calls do not identify dropout versus recombination', UserWarning)
        return None
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

def summarize_singlet_counts(adata_obs, split_by='sample', valid_pairs=None, valid_combos=None, include_recomb_multiplet=False, return_full_df=False, ref_df=None):
    """
    Summarize the number of cells per library construct (pair_id) and per combination,
    directly parsing an adata.obs dataframe containing cell-level doublet classifications.
    
    Args:
        adata_obs: The full metadata DataFrame.
        split_by: Column to split the analysis by (e.g. 'sample').
        valid_pairs: Optional explicit pair_ids to capture 0-count drops.
        valid_combos: Optional explicit gene-level combos to capture 0-count drops.
        include_recomb_multiplet: If True, evaluates scenario-1 physically recombined mismatched targets 
                                  as viable combos, ignoring their structurally broken pair_ids.
        return_full_df: If True, returns a tuple (summary_df, full_unpivoted_counts_df).
        ref_df: Optional reference design dataframe. Automatically resolves expected pair_ids and target_combos.
    """
    import pandas as pd
    import numpy as np
    
    # Check what the combo column is named naturally
    combo_col = 'target_combo' if 'target_combo' in adata_obs.columns else 'target_gene'
    
    pair_to_harmonized = {}
    if ref_df is not None:
        if 'pair_id' in ref_df.columns and 'position' in ref_df.columns and 'grna_name' in ref_df.columns:
            pos1_ref = ref_df[ref_df['position'] == 1].set_index('pair_id')['grna_name']
            pos2_ref = ref_df[ref_df['position'] == 2].set_index('pair_id')['grna_name']
            harmonized_ref = pos1_ref.astype(str) + '_' + pos2_ref.astype(str)
            pair_to_harmonized = harmonized_ref.dropna().to_dict()
            if valid_pairs is None:
                valid_pairs = harmonized_ref.dropna().unique()
        elif valid_pairs is None and 'pair_id' in ref_df.columns:
            valid_pairs = ref_df['pair_id'].dropna().unique()

        if valid_combos is None and 'target_combo' in ref_df.columns:
            valid_combos = ref_df['target_combo'].dropna().unique()
            
    # Guarantee iterable structures are isolated string sets preventing Pandas Index lookup collision
    if valid_pairs is not None:
        valid_pairs_set = set(str(x) for x in valid_pairs)
        valid_pairs = list(valid_pairs_set)
    if valid_combos is not None:
        valid_combos_set = set(str(x) for x in valid_combos)
        valid_combos = list(valid_combos_set)
    
    # Find natural subsets
    perfect_df = adata_obs[(adata_obs['combo_type'] == 'singlet') & (adata_obs['detection_status'] == 'complete')].copy()
    inferred_df = adata_obs[(adata_obs['combo_type'] == 'singlet') & (adata_obs['detection_status'] == 'inferred')].copy()
    
    perfect_df['combo_eval'] = perfect_df[combo_col]
    inferred_df['combo_eval'] = inferred_df[combo_col]
    perfect_df['pair_eval'] = perfect_df['pair_id'].map(pair_to_harmonized).fillna(perfect_df['pair_id'])
    inferred_df['pair_eval'] = inferred_df['pair_id'].map(pair_to_harmonized).fillna(inferred_df['pair_id'])
    
    if include_recomb_multiplet:
        m_mask = (adata_obs['combo_type'] == 'multiplet') & (adata_obs['detection_status'] == 'inferred')
        m_mask &= (adata_obs['pos1_grna'] != 'unassigned') & (~adata_obs['pos1_grna'].astype(str).str.contains(',', na=False))
        m_mask &= (adata_obs['pos2_grna'] != 'unassigned') & (~adata_obs['pos2_grna'].astype(str).str.contains(',', na=False))
        recomb_df = adata_obs[m_mask].copy()
        
        # Synthesize dynamically recombined combo metrics crossing pos1 and pos2 independently
        recomb_df['combo_eval'] = recomb_df['pos1_target'].astype(str) + '_' + recomb_df['pos2_target'].astype(str)
        recomb_df['pair_eval'] = recomb_df['pos1_grna'].astype(str) + '_' + recomb_df['pos2_grna'].astype(str)
        recomb_df['pair_id'] = recomb_df['pair_eval']
    else:
        recomb_df = pd.DataFrame()
        
    splits = [None] if not split_by else adata_obs[split_by].dropna().unique()
    
    res_list = []
    full_df_records = []
    
    print(f"{'='*60}")
    print(f"       Singlets Recovery Summary")
    print(f"{'='*60}")
    
    for sp in splits:
        sp_label = "All Samples" if sp is None else f"Sample {sp}"
        print(f"\n[{sp_label}]")
        
        def compute_stats(df, label):
            sub_df = df if sp is None else df[df[split_by] == sp]
            
            p_counts = sub_df['pair_eval'].value_counts()
            p_counts = p_counts[p_counts > 0]
            novel_pairs = []
            if valid_pairs is not None:
                novel_pairs = [p for p in p_counts.index if str(p) not in valid_pairs_set]
                all_valid_p = list(valid_pairs_set.union(set(str(x) for x in p_counts.index)))
                p_counts = p_counts.reindex(all_valid_p, fill_value=0)
            
            if valid_pairs is not None:
                p_stats = p_counts.loc[list(valid_pairs_set)]
                p_desc = p_stats.describe()[['min', 'max', '50%', 'mean']] if len(p_stats) > 0 else pd.Series({'min':0, 'max':0, '50%':0, 'mean':0})
                p_detected = (p_stats > 0).sum()
                p_expected = len(valid_pairs)
            else:
                p_desc = p_counts.describe()[['min', 'max', '50%', 'mean']] if len(p_counts) > 0 else pd.Series({'min':0, 'max':0, '50%':0, 'mean':0})
                p_detected = (p_counts > 0).sum()
                p_expected = len(p_counts)
            
            c_counts = sub_df['combo_eval'].value_counts()
            c_counts = c_counts[c_counts > 0]
            novel_combos = []
            if valid_combos is not None:
                novel_combos = [c for c in c_counts.index if str(c) not in valid_combos_set]
                all_valid_c = list(valid_combos_set.union(set(str(x) for x in c_counts.index)))
                c_counts = c_counts.reindex(all_valid_c, fill_value=0)
                
            if valid_combos is not None:
                c_stats = c_counts.loc[list(valid_combos_set)]
                c_desc = c_stats.describe()[['min', 'max', '50%', 'mean']] if len(c_stats) > 0 else pd.Series({'min':0, 'max':0, '50%':0, 'mean':0})
                c_detected = (c_stats > 0).sum()
                c_expected = len(valid_combos)
            else:
                c_desc = c_counts.describe()[['min', 'max', '50%', 'mean']] if len(c_counts) > 0 else pd.Series({'min':0, 'max':0, '50%':0, 'mean':0})
                c_detected = (c_counts > 0).sum()
                c_expected = len(c_counts)
            
            print(f"  {label} (Total: {len(sub_df):,})")
            print(f"    - Per Construct ({p_detected:,}/{p_expected:,} pairs) -> Mean: {p_desc['mean']:.1f} | Median: {p_desc['50%']:.1f} | Min: {p_desc['min']:.0f} | Max: {p_desc['max']:.0f}")
            print(f"    - Per Combo     ({c_detected:,}/{c_expected:,} genes) -> Mean: {c_desc['mean']:.1f} | Median: {c_desc['50%']:.1f} | Min: {c_desc['min']:.0f} | Max: {c_desc['max']:.0f}")
            
            missing_pairs = p_counts[p_counts == 0].index.tolist()
            if valid_pairs is not None:
                missing_pairs = [p for p in missing_pairs if str(p) in valid_pairs_set]
            if len(missing_pairs) > 0:
                print(f"      [!] Missing Pairs ({len(missing_pairs)}): {', '.join(map(str, missing_pairs[:10]))}{' ...' if len(missing_pairs) > 10 else ''}")
            if len(novel_pairs) > 0:
                print(f"      [+] Novel Pairs ({len(novel_pairs)}): {', '.join(map(str, novel_pairs[:10]))}{' ...' if len(novel_pairs) > 10 else ''}")
                
            missing_combos = c_counts[c_counts == 0].index.tolist()
            if valid_combos is not None:
                missing_combos = [c for c in missing_combos if str(c) in valid_combos_set]
            if len(missing_combos) > 0:
                print(f"      [!] Missing Combos ({len(missing_combos)}): {', '.join(map(str, missing_combos[:10]))}{' ...' if len(missing_combos) > 10 else ''}")
            if len(novel_combos) > 0:
                print(f"      [+] Novel Combos ({len(novel_combos)}): {', '.join(map(str, novel_combos[:10]))}{' ...' if len(novel_combos) > 10 else ''}")
            
            res_list.append({
                'split': sp,
                'subset': label,
                'total_cells': len(sub_df),
                'pairs_detected': p_detected, 'pairs_expected': p_expected,
                'pair_mean': p_desc['mean'], 'pair_median': p_desc['50%'], 
                'pair_min': p_desc['min'], 'pair_max': p_desc['max'],
                'combos_detected': c_detected, 'combos_expected': c_expected,
                'combo_mean': c_desc['mean'], 'combo_median': c_desc['50%'], 
                'combo_min': c_desc['min'], 'combo_max': c_desc['max'],
            })
            
            if return_full_df:
                for k, v in p_counts.items():
                    full_df_records.append({'split': sp, 'subset': label, 'category': 'pair_id', 'name': k, 'n_cells': v})
                for k, v in c_counts.items():
                    full_df_records.append({'split': sp, 'subset': label, 'category': 'target_combo', 'name': k, 'n_cells': v})
        
        compute_stats(perfect_df, "Perfect Singlets")
        
        combined_inf_df = pd.concat([perfect_df, inferred_df], ignore_index=True) if not inferred_df.empty else perfect_df
        compute_stats(combined_inf_df, "Perfect + Inferred")
            
        if include_recomb_multiplet and not recomb_df.empty:
            combined_recomb_df = pd.concat([combined_inf_df, recomb_df], ignore_index=True)
            compute_stats(combined_recomb_df, "Perfect + Inferred + Recomb")
            
    print(f"{'-'*60}\n")
    
    out_df = pd.DataFrame(res_list)
    if split_by is None and 'split' in out_df.columns:
        out_df = out_df.drop(columns=['split'])
        
    if return_full_df:
        full_df = pd.DataFrame(full_df_records)
        if split_by is None and 'split' in full_df.columns:
            full_df = full_df.drop(columns=['split'])
        return out_df, full_df
        
    return out_df

def evaluate_recombinant_multiplets_kd(adata, sample_name, exclude_combos=None, min_cells=5):
    """
    Computes KD metrics explicitly for Scenario-1 cells (1 guide pos1 + 1 guide pos2 mismatched),
    evaluating DE logic for the 2 detected active guides vs the 2 inferred missing structural guides.
    """
    import pandas as pd
    
    if exclude_combos is None:
        exclude_combos = []
        
    obs = adata.obs
    
    # 1. Extract Perfect Singlets for Control
    perfect_mask = (obs['combo_type'] == 'singlet') & (obs['detection_status'] == 'complete') & (obs['sample'] == sample_name)
    perfect_df = obs[perfect_mask]
    ctrl_bcs = list(perfect_df[perfect_df['target_combo'] == 'ctrl_ctrl'].index)
    
    if not ctrl_bcs:
        print(f"  [{sample_name}] Warning: No ctrl_ctrl perfect singlets found.")
        return pd.DataFrame()
        
    # 2. Identify Scenario-1 Cells
    m_mask = (obs['combo_type'] == 'multiplet') & (obs['detection_status'] == 'inferred') & (obs['sample'] == sample_name)
    m_mask &= (obs['pos1_grna'] != 'unassigned') & (~obs['pos1_grna'].str.contains(',', na=False))
    m_mask &= (obs['pos2_grna'] != 'unassigned') & (~obs['pos2_grna'].str.contains(',', na=False))
    recomb_df = obs[m_mask].copy()
    
    if exclude_combos:
        to_drop = []
        for idx, row in recomb_df.iterrows():
            combos_involved = [
                f"{row.get('pos1_target', '')}_{row.get('pos2_target', '')}", 
                f"{row.get('pos1_target', '')}_{row.get('inferred_pos2_target', '')}", 
                f"{row.get('inferred_pos1_target', '')}_{row.get('pos2_target', '')}" 
            ]
            if any(c in exclude_combos for c in combos_involved):
                to_drop.append(idx)
        recomb_df = recomb_df.drop(index=to_drop)
        
    recomb_df['bc_wells'] = recomb_df.index
    print(f"  [{sample_name}] Scenario-1 Recombined Cells: {len(recomb_df):,}")
    
    if len(recomb_df) == 0:
        return pd.DataFrame()
        
    def get_groups(df, gene_col, guide_col, mode_desc, pos_label):
        grp_list = []
        def is_control_target(g):
            if pd.isna(g): return True
            gs = str(g).lower()
            return gs in ['unassigned', '', 'none', 'nan'] or 'ctrl' in gs or 'safe' in gs or 'ntc' in gs
            
        pooled_bcs = {g: list(set(d['bc_wells']).intersection(adata.obs.index)) 
                      for g, d in df.groupby(gene_col) if not is_control_target(g)}
                      
        for (gene, guide), d in df.groupby([gene_col, guide_col]):
            if is_control_target(gene): continue
            bcs = list(set(d['bc_wells']).intersection(adata.obs.index))
            if len(bcs) < min_cells:
                fb_bcs = pooled_bcs.get(gene, [])
                if len(fb_bcs) < min_cells: continue
                if not any(x['actual_guide'] == f"{gene}_all (fallback)" for x in grp_list):
                    grp_list.append({
                        'target_gene': gene,
                        'guide_name': mode_desc,
                        'position': pos_label,
                        'target_combo': f"{gene}_all (fallback)",
                        'cell_barcodes': fb_bcs,
                        'actual_guide': f"{gene}_all (fallback)"
                    })
            else:
                grp_list.append({
                    'target_gene': gene,
                    'guide_name': mode_desc,
                    'position': pos_label,
                    'target_combo': guide,
                    'cell_barcodes': bcs,
                    'actual_guide': guide
                })
        return grp_list

    groups = []
    groups.extend(get_groups(recomb_df, 'pos1_target', 'pos1_grna', 'Detected (Recomb)', 'pos1'))
    groups.extend(get_groups(recomb_df, 'pos2_target', 'pos2_grna', 'Detected (Recomb)', 'pos2'))
    groups.extend(get_groups(recomb_df, 'inferred_pos1_target', 'inferred_pos1_grna', 'Inferred (Missing)', 'inferred_pos1'))
    groups.extend(get_groups(recomb_df, 'inferred_pos2_target', 'inferred_pos2_grna', 'Inferred (Missing)', 'inferred_pos2'))
    
    if not groups:
        return pd.DataFrame()
        
    metrics_df = compute_kd_metrics(adata, groups, ctrl_bcs, sample_name, min_cells=min_cells)
    if not metrics_df.empty:
        metrics_df['actual_guide'] = metrics_df['target_combo']
        
    return metrics_df


if __name__ == "__main__":
    main()

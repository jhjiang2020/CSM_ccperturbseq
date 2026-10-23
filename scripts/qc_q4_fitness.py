#!/usr/bin/env python3
"""
QC Q4: Fitness Analysis

Compares guide-combo frequencies between pre-lentiviral bulk library (from Q1)
and single-cell data to estimate fitness effects.

Functions:
  - analyze_fitness: cell-level frequency comparison (from old Q2).
  - identify_novel_combos: finds combos present in sc but absent from pre-lenti.
  - estimate_fitness: per-combo log2FC(sc_freq / prelenti_freq).

Usage:
  python scripts/qc_q4_fitness.py --config config.yaml \\
      --combo_freq path/to/q1_combo_frequencies.csv
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

_project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)

from scripts.qc_config_loader import load_qc_config
from scripts.guide_design import eligible_assignments


# ============================================================================
# Fitness: cell-level frequency comparison (migrated from old Q2)
# ============================================================================

def analyze_fitness(assign_path, combo_freq_path, out_dir):
    """Compare combo frequencies in cells vs pre-lentiviral library.

    Args:
        assign_path: Path to grna_cell_assignment.csv
        combo_freq_path: CSV with columns pos1, pos2, count, frequency
        out_dir: Output directory

    Returns:
        DataFrame with merged frequency and log2FC, or None if skipped.
    """
    print("\n  --- Q4 Part A: Cell-Level Fitness Analysis ---")

    if combo_freq_path is None or not os.path.exists(str(combo_freq_path)):
        print("    No pre-lentiviral combo frequency table provided. Skipping.")
        print("    Provide via --combo_freq to enable this analysis.")
        return None

    assign_df = eligible_assignments(pd.read_csv(assign_path))
    combo_freq = pd.read_csv(combo_freq_path)

    print(f"    Pre-lentiviral combos: {len(combo_freq)}")

    # Post-selection frequencies from cell assignments
    post_counts = assign_df['target_combo'].value_counts().reset_index()
    post_counts.columns = ['target_combo', 'cell_count']
    post_total = post_counts['cell_count'].sum()
    post_counts['cell_frequency'] = post_counts['cell_count'] / post_total

    # Construct target_combo if needed
    if 'target_combo' not in combo_freq.columns:
        if 'pos1' in combo_freq.columns and 'pos2' in combo_freq.columns:
            combo_freq['target_combo'] = combo_freq['pos1'] + '_' + combo_freq['pos2']
        else:
            print("    WARNING: combo_freq table doesn't have expected columns. Skipping.")
            return None

    # Normalize pre-lentiviral frequencies
    if 'frequency' not in combo_freq.columns:
        if 'count' in combo_freq.columns:
            combo_freq['frequency'] = combo_freq['count'] / combo_freq['count'].sum()
        else:
            print("    WARNING: No frequency or count column in combo_freq. Skipping.")
            return None

    merged = post_counts.merge(
        combo_freq[['target_combo', 'frequency']].rename(
            columns={'frequency': 'pre_frequency'}),
        on='target_combo', how='outer'
    )
    merged['cell_frequency'] = merged['cell_frequency'].fillna(0)
    merged['pre_frequency'] = merged['pre_frequency'].fillna(0)

    # log2FC of frequency
    pseudo = 1e-6
    merged['log2FC_fitness'] = np.log2(
        (merged['cell_frequency'] + pseudo) / (merged['pre_frequency'] + pseudo)
    )

    merged = merged.sort_values('log2FC_fitness', ascending=True)

    # Save
    fitness_path = os.path.join(out_dir, 'q4_fitness_analysis.csv')
    merged.to_csv(fitness_path, index=False)
    print(f"    Saved: {fitness_path}")

    # Plot
    fig, ax = plt.subplots(figsize=(12, max(6, len(merged) * 0.3)))
    colors = ['#e74c3c' if x < -0.5 else '#2ecc71' if x > 0.5 else '#95a5a6'
              for x in merged['log2FC_fitness']]
    ax.barh(range(len(merged)), merged['log2FC_fitness'], color=colors, alpha=0.7)
    ax.set_yticks(range(len(merged)))
    ax.set_yticklabels(merged['target_combo'], fontsize=7)
    ax.set_xlabel('log₂FC (cell frequency / library frequency)')
    ax.set_title('Combo Fitness: Post-selection vs Pre-lentiviral')
    ax.axvline(0, color='black', lw=0.5)
    ax.grid(True, alpha=0.2, axis='x')
    plt.tight_layout()
    plot_path = os.path.join(out_dir, 'q4_fitness_barplot.png')
    fig.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved: {plot_path}")

    # Print top enriched/depleted
    print(f"\n    Top 5 enriched combos:")
    for _, row in merged.tail(5).iloc[::-1].iterrows():
        print(f"      {row['target_combo']:>25s}: log2FC={row['log2FC_fitness']:+.3f}")
    print(f"    Top 5 depleted combos:")
    for _, row in merged.head(5).iterrows():
        print(f"      {row['target_combo']:>25s}: log2FC={row['log2FC_fitness']:+.3f}")

    return merged


# ============================================================================
# Novel combo detection (migrated from old Q3-fitness)
# ============================================================================

def identify_novel_combos(cell_assign_df, gseq_freq_df, adata_obs, out_dir):
    """Identify guide combinations in scRNA-seq absent from pre-lenti library.

    Args:
        cell_assign_df: Cell assignment DataFrame
        gseq_freq_df: Pre-lenti frequency DataFrame (from Q1 get_combo_frequencies)
        adata_obs: adata.obs DataFrame (for sample mapping)
        out_dir: Output directory

    Returns:
        DataFrame with merged sc/pre-lenti frequencies and novel flag.
    """
    print("\n  --- Q4 Part B: Novel Combo Detection ---")

    assigned = cell_assign_df[
        (cell_assign_df['combo_type'] == 'singlet') &
        (cell_assign_df['detection_status'] == 'complete')
    ].copy()
    if assigned.empty:
        print("    No perfect matched singlet cells found in assignment.")
        return pd.DataFrame()

    assigned['pos1_grna'] = assigned['pos1_grna'].fillna('')
    assigned['pos2_grna'] = assigned['pos2_grna'].fillna('')

    if 'sample' not in assigned.columns:
        if 'bc_wells' in assigned.columns:
            assigned['sample'] = assigned['bc_wells'].map(adata_obs['sample'])
        else:
            print("    WARNING: 'bc_wells' not found, unable to map sample.")
            assigned['sample'] = 'Unknown'

    # scRNA-seq frequency per sample
    scrna_counts = assigned.groupby(
        ['pos1_grna', 'pos2_grna', 'target_combo', 'sample']
    ).size().reset_index(name='scrna_count')
    sample_totals = scrna_counts.groupby('sample')['scrna_count'].transform('sum')
    scrna_counts['scrna_frequency'] = scrna_counts['scrna_count'] / sample_totals

    # Merge with pre-lenti
    gseq_freq_df = gseq_freq_df.copy()
    gseq_freq_df['in_prelenti'] = True
    merged = pd.merge(
        scrna_counts,
        gseq_freq_df[['pos1_grna', 'pos2_grna', 'count', 'frequency', 'in_prelenti']],
        on=['pos1_grna', 'pos2_grna'],
        how='left'
    )

    merged['in_prelenti'] = merged['in_prelenti'].fillna(False)
    merged.rename(columns={'count': 'prelenti_count', 'frequency': 'prelenti_frequency'}, inplace=True)
    merged['prelenti_count'] = merged['prelenti_count'].fillna(0)
    merged['prelenti_frequency'] = merged['prelenti_frequency'].fillna(0.0)

    novel = merged[~merged['in_prelenti']].copy()
    print(f"    Found {len(novel):,} novel guide combinations in scRNA-seq (across samples).")

    novel_path = os.path.join(out_dir, 'q4_novel_combos.csv')
    novel.sort_values('scrna_count', ascending=False).to_csv(novel_path, index=False)
    print(f"    Saved novel combos to {novel_path}")

    return merged


# ============================================================================
# Fitness estimation: per-combo log2FC (migrated from old Q3-fitness)
# ============================================================================

def estimate_fitness(merged_df, out_dir):
    """Estimate fitness effect: log2FC(scRNA-seq freq / pre-lenti freq).

    Args:
        merged_df: Output from identify_novel_combos
        out_dir: Output directory
    """
    print("\n  --- Q4 Part C: Fitness Effect Estimation ---")

    valid = merged_df[
        (merged_df['prelenti_count'] >= 10) &
        (merged_df['scrna_count'] >= 5) &
        (merged_df['target_combo'].notna()) &
        (merged_df['target_combo'] != '')
    ].copy()

    if valid.empty:
        print("    Not enough valid overlapping combos for fitness estimation.")
        return

    print(f"    Analyzing {len(valid)} valid overlapping combinations across samples.")

    valid['fitness_log2fc'] = np.log2(valid['scrna_frequency'] / valid['prelenti_frequency'])

    fitness_path = os.path.join(out_dir, 'q4_combos_fitness.csv')
    valid.sort_values(['sample', 'fitness_log2fc'], ascending=[True, False]).to_csv(fitness_path, index=False)
    print(f"    Saved combo fitness to {fitness_path}")

    # Aggregate by target gene combo + sample
    gene_fitness = valid.groupby(['target_combo', 'sample']).agg(
        n_guides=('pos1_grna', 'count'),
        total_cells=('scrna_count', 'sum'),
        mean_guide_fitness=('fitness_log2fc', 'mean'),
        median_guide_fitness=('fitness_log2fc', 'median'),
        scrna_freq_sum=('scrna_frequency', 'sum'),
        prelenti_freq_sum=('prelenti_frequency', 'sum')
    ).reset_index()

    gene_fitness['bulk_fitness_log2fc'] = np.log2(gene_fitness['scrna_freq_sum'] / gene_fitness['prelenti_freq_sum'])

    g_fitness_path = os.path.join(out_dir, 'q4_gene_target_fitness.csv')
    gene_fitness.sort_values(['sample', 'mean_guide_fitness'], ascending=[True, False]).to_csv(g_fitness_path, index=False)
    print(f"    Saved gene-level fitness to {g_fitness_path}")

    # Plot for each sample
    for sample in gene_fitness['sample'].unique():
        sample_df = gene_fitness[gene_fitness['sample'] == sample]
        if len(sample_df) >= 10:
            _plot_fitness(sample_df, out_dir, sample)


def _plot_fitness(gene_fitness, out_dir, sample_name):
    """Plot top/bottom fitness effects."""
    gf = gene_fitness.sort_values('mean_guide_fitness').copy()

    depleted = gf.head(15)
    enriched = gf.tail(15)
    ctrl = gf[gf['target_combo'] == 'ctrl_ctrl']

    to_plot = pd.concat([depleted, ctrl, enriched]).drop_duplicates()
    to_plot = to_plot.sort_values('mean_guide_fitness')

    fig, ax = plt.subplots(figsize=(10, max(6, len(to_plot) * 0.3)))

    colors = ['#e74c3c' if val < 0 else '#2ecc71' for val in to_plot['mean_guide_fitness']]
    color_mask = [c if tc != 'ctrl_ctrl' else 'grey' for c, tc in zip(colors, to_plot['target_combo'])]

    y_pos = np.arange(len(to_plot))
    ax.barh(y_pos, to_plot['mean_guide_fitness'], color=color_mask, alpha=0.8)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(to_plot['target_combo'], fontsize=8)
    ax.axvline(0, color='k', linestyle='-', linewidth=1)

    if not ctrl.empty:
        ctrl_val = ctrl['mean_guide_fitness'].iloc[0]
        ax.axvline(ctrl_val, color='grey', linestyle='--', alpha=0.5, label='ctrl_ctrl origin')
        ax.legend(fontsize=9)

    ax.set_xlabel('Fitness Proxy: log2(scRNA-seq freq / pre-lenti freq)', fontsize=11)
    ax.set_title(f'Target Combination Fitness Effects ({sample_name})\n(Mean across all guide pairs)', fontsize=13)
    ax.grid(True, alpha=0.2, axis='x')

    plt.tight_layout()
    plot_path = os.path.join(out_dir, f'q4_fitness_effects_{sample_name}.png')
    fig.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved fitness plot to {plot_path}")


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="QC Q4: Fitness Analysis"
    )
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument("--combo_freq", type=str, default=None,
                        help="Path to pre-lentiviral combo frequency CSV "
                             "(from Q1 get_combo_frequencies)")
    parser.add_argument("--guide_freq", type=str, default=None,
                        help="Path to pre-lentiviral guide frequency CSV "
                             "(from Q1 get_guide_frequencies)")
    args = parser.parse_args()

    cfg = load_qc_config(args.config)
    output_dir = cfg['paths']['output_dir']

    print(f"\n{'='*60}")
    print(f"  QC Q4: Fitness Analysis")
    print(f"{'='*60}")

    out_dir = os.path.join(output_dir, 'qc_analysis')
    os.makedirs(out_dir, exist_ok=True)

    assign_path = os.path.join(output_dir, 'grna_assignment', 'grna_cell_assignment.csv')

    # --- Part A: Cell-level fitness (old Q2 function) ---
    analyze_fitness(assign_path, args.combo_freq, out_dir)

    # --- Parts B + C: Novel combos + fitness estimation ---
    if args.combo_freq and os.path.exists(args.combo_freq):
        import scanpy as sc
        adata_path = os.path.join(output_dir, 'combined_adata.h5ad')
        adata = sc.read_h5ad(adata_path)

        assign_df = eligible_assignments(pd.read_csv(assign_path))
        combo_freq_df = pd.read_csv(args.combo_freq)

        # Standardize column names for identify_novel_combos
        if 'pos1_grna' not in combo_freq_df.columns:
            if 'pos1' in combo_freq_df.columns:
                combo_freq_df = combo_freq_df.rename(columns={'pos1': 'pos1_grna', 'pos2': 'pos2_grna'})

        merged = identify_novel_combos(assign_df, combo_freq_df, adata.obs, out_dir)
        if not merged.empty:
            estimate_fitness(merged, out_dir)

    print(f"\n{'='*60}")
    print(f"  QC Q4 Complete")
    print(f"  All results saved to: {out_dir}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()

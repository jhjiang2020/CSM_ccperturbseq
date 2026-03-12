#!/usr/bin/env python3
"""
Step 5: Build combined AnnData with GEX + gRNA information.

Merges separate GEX AnnData files (per sublibrary) with gRNA assignment data
from Step 4 into a single AnnData object.

Standardizes cell barcodes to the format "{sublib_index}_{barcode_seq}" 
to ensure global uniqueness and alignment with Step 4 output.

Usage:
  python scripts/05_build_anndata.py --config config.yaml
"""

import argparse
import os
import sys
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from scipy import io as sio
from scipy import sparse


def load_config(config_path):
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_sample_gex(output_dir, sample_name, sublibs_info):
    """
    Load combined GEX AnnData for a sample and resolve sublibrary indices.
    
    Args:
        output_dir: Path to pipeline output directory
        sample_name: Name of the sample (e.g., 'Control')
        sublibs_info: List of tuples (global_index, lib_entry) for this sample.
                      The list MUST be ordered by the sublibrary order in config.yaml,
                      which presumably matches split-pipe's __s1, __s2 assignment order.
    """
    # 1. Load Combined AnnData
    # Path: output_dir/gex_combined/{sample_name}/DGE_filtered/anndata.h5ad
    gex_path = os.path.join(output_dir, "gex_combined", sample_name, 'DGE_filtered', 'anndata.h5ad')
    if not os.path.exists(gex_path):
        print(f"  WARNING: Could not find GEX AnnData for sample '{sample_name}' at {gex_path}")
        return None

    print(f"    Loading sample {sample_name}: {gex_path}")
    adata = sc.read_h5ad(gex_path)
    
    # 2. Load Metadata CSV
    # Path: output_dir/gex_combined/{sample_name}/DGE_filtered/cell_metadata.csv
    meta_path = os.path.join(output_dir, "gex_combined", sample_name, 'DGE_filtered', 'cell_metadata.csv')
    if os.path.exists(meta_path):
        print(f"    Loading metadata: {meta_path}")
        meta_df = pd.read_csv(meta_path)
        # Assuming 'bc_wells' is the barcode column matching adata.obs.index
        if 'bc_wells' in meta_df.columns:
            meta_df.set_index('bc_wells', inplace=True)
            # Merge with existing obs, keeping left index
            cols_to_add = [c for c in meta_df.columns if c not in adata.obs.columns]
            if cols_to_add:
                adata.obs = adata.obs.join(meta_df[cols_to_add], how='left')
    
    # 3. Resolve Sublibrary Indices
    # Map suffix number (1, 2, ...) -> global_index
    # The sublibs_info list is ordered [ (idx_C1, lib_C1), (idx_C2, lib_C2) ]
    suffix_to_global = {}
    suffix_to_name = {}
    
    for i, (global_idx, lib_entry) in enumerate(sublibs_info):
        suffix_num = str(i + 1) # '1', '2'
        suffix_to_global[suffix_num] = global_idx
        suffix_to_name[suffix_num] = lib_entry['sublibrary_name']
        
    print(f"    Mapping suffixes for {sample_name}: {suffix_to_name}")

    # Create new index
    # Current index: "01_01_06__s1" -> Target: "{global_idx}_01_01_06"
    
    # Use pandas vector operations
    index_series = pd.Series(adata.obs.index)
    extracted = index_series.str.extract(r'^(.*)__s(\d+)$')
    
    # Extracted[0] = base barcode, Extracted[1] = suffix digit
    valid_mask = extracted[1].notna()
    
    if (~valid_mask).any():
        print(f"    WARNING: {(~valid_mask).sum()} cells do not match suffix pattern (e.g. __sN). Dropping them.")
        # Filter adata to valid cells
        valid_indices = index_series[valid_mask].index # numeric indices in the series
        adata = adata[valid_indices].copy()
        extracted = extracted.loc[valid_indices]
    
    # Check for unknown suffixes
    found_suffixes = extracted[1].unique()
    unknown = [s for s in found_suffixes if s not in suffix_to_global]
    if unknown:
        print(f"    WARNING: Found suffixes {unknown} not in expected {list(suffix_to_global.keys())}. Dropping these cells.")
        valid_suffix_mask = extracted[1].isin(suffix_to_global.keys())
        adata = adata[valid_suffix_mask].copy()
        extracted = extracted[valid_suffix_mask]
        
    # Construct new barcodes: global_idx + "_" + base_bc
    global_idxs = extracted[1].map(suffix_to_global).astype(str)
    new_barcodes = global_idxs + "_" + extracted[0]
    
    adata.obs.index = new_barcodes
    adata.obs['sublib_index'] = extracted[1].map(suffix_to_global).values
    adata.obs['sublibrary'] = extracted[1].map(suffix_to_name).values
    adata.obs['sample'] = sample_name
    
    # Fix sample column if needed
    if 'sample' in adata.obs.columns and not adata.obs['sample'].equals(pd.Series([sample_name]*len(adata), index=adata.obs.index)):
         adata.obs['sample_original'] = adata.obs['sample']
    adata.obs['sample'] = sample_name

    return adata


def main():
    parser = argparse.ArgumentParser(
        description="Build combined AnnData with GEX + gRNA data"
    )
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    output_dir = cfg['paths']['output_dir']
    assign_dir = os.path.join(output_dir, 'grna_assignment')

    print(f"\n{'='*60}")
    print(f"  Step 5: Build Combined AnnData")
    print(f"{'='*60}")

    # --- 1. Load and Combine GEX AnnData ---
    print("\n  1. Loading combined GEX AnnData files per sample...")
    
    libraries = cfg.get('libraries', [])
    if not libraries:
        raise ValueError("No libraries defined in config.yaml")

    # Group libraries by sample to handle combined outputs
    sample_groups = {}
    for i, lib in enumerate(libraries):
        s_name = lib['sample_name']
        if s_name not in sample_groups:
            sample_groups[s_name] = []
        sample_groups[s_name].append((i, lib))

    adatas = []
    for sample_name, sublibs in sample_groups.items():
        ad = load_sample_gex(output_dir, sample_name, sublibs)
        if ad is not None:
            adatas.append(ad)

    if not adatas:
        raise FileNotFoundError("No GEX AnnData files found for any sample.")

    print(f"\n  Concatenating {len(adatas)} samples...")
    # outer join to keep all genes if they differ (though they shouldn't for same ref)
    adata = sc.concat(adatas, join='outer', index_unique=None)
    
    # Ensure unique indices (Step 4 ensures uniqueness via prefix, but safety check)
    if not adata.obs.index.is_unique:
        print("  WARNING: Duplicate cell barcodes found after concatenation!")
        adata.obs_names_make_unique()

    print(f"     Combined Shape: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    print(f"     obs columns: {list(adata.obs.columns)}")

    # --- 2. Load gRNA cell assignments ---
    print("\n  2. Loading gRNA assignments...")
    assign_path = os.path.join(assign_dir, 'grna_cell_assignment.csv')
    if not os.path.exists(assign_path):
        raise FileNotFoundError(f"Assignment file not found: {assign_path}\nRun Step 4 first.")

    assign_df = pd.read_csv(assign_path)
    print(f"     Loaded {len(assign_df):,} assignments for "
          f"{assign_df['bc_wells'].nunique():,} cells")

    # --- 3. Load gRNA UMI matrix ---
    print("\n  3. Loading gRNA UMI matrix...")
    mtx_path = os.path.join(assign_dir, 'grna_umi_matrix.mtx')
    barcodes_path = os.path.join(assign_dir, 'grna_matrix_barcodes.txt')
    guides_path = os.path.join(assign_dir, 'grna_matrix_guides.txt')

    guide_matrix = sio.mmread(mtx_path).tocsr()
    guide_barcodes = pd.read_csv(barcodes_path, header=None)[0].tolist()
    guide_names = pd.read_csv(guides_path, header=None)[0].tolist()
    print(f"     Matrix: {guide_matrix.shape[0]:,} cells × {guide_matrix.shape[1]:,} guides")

    # --- 4. Match cell barcodes between GEX and gRNA ---
    print("\n  4. Matching cell barcodes...")
    # GEX barcodes are now "{i}_Seq"
    # gRNA barcodes from Step 4 are also "{i}_Seq"
    gex_barcodes = set(adata.obs.index)
    grna_barcodes = set(guide_barcodes)
    overlap = gex_barcodes & grna_barcodes

    print(f"     GEX cells:     {len(gex_barcodes):,}")
    print(f"     gRNA cells:    {len(grna_barcodes):,}")
    print(f"     Overlap:       {len(overlap):,} "
          f"({len(overlap)/len(gex_barcodes)*100:.1f}% of GEX)")

    if len(overlap) == 0:
        print("\n  WARNING: No overlapping barcodes found!")
        print("  This may indicate a barcode format mismatch.")
        print("  GEX barcode examples:", list(gex_barcodes)[:5])
        print("  gRNA barcode examples:", list(grna_barcodes)[:5])
        print("  Proceeding without gRNA data — saving GEX-only AnnData.")
    
    # --- 6. Add gRNA information to AnnData ---
    print("\n  6. Adding gRNA information to AnnData...")

    # Create per-cell summary from assignment data
    is_dual = 'combo_type' in assign_df.columns

    if is_dual:
        agg_dict = {
            'combo_type': ('combo_type', 'first'),  # same for all rows of a cell
            'detection_status': ('detection_status', lambda x: ','.join(sorted(x.unique()))),
            'pair_id': ('pair_id', lambda x: ','.join(sorted(x.dropna().astype(str).unique()))),
            'target_combo': ('target_combo', lambda x: ','.join(sorted(x.dropna().astype(str).unique()))),
            'gi_pair': ('gi_pair', lambda x: ','.join(sorted(x.dropna().astype(str).unique()))),
            
            'pos1_grna': ('pos1_grna', lambda x: ','.join(sorted(x.replace('', pd.NA).dropna().astype(str).unique()))),
            'pos1_ibar': ('pos1_ibar', lambda x: ','.join(sorted(x.replace('', pd.NA).dropna().astype(str).unique()))),
            'pos1_target': ('pos1_target', lambda x: ','.join(sorted(x.replace('', pd.NA).dropna().astype(str).unique()))),
            'pos1_umi': ('pos1_umi', 'sum'),
            
            'pos2_grna': ('pos2_grna', lambda x: ','.join(sorted(x.replace('', pd.NA).dropna().astype(str).unique()))),
            'pos2_ibar': ('pos2_ibar', lambda x: ','.join(sorted(x.replace('', pd.NA).dropna().astype(str).unique()))),
            'pos2_target': ('pos2_target', lambda x: ','.join(sorted(x.replace('', pd.NA).dropna().astype(str).unique()))),
            'pos2_umi': ('pos2_umi', 'sum'),
            
            'inferred_pos1_grna': ('inferred_pos1_grna', lambda x: ','.join(sorted(x.replace('', pd.NA).dropna().astype(str).unique()))),
            'inferred_pos1_target': ('inferred_pos1_target', lambda x: ','.join(sorted(x.replace('', pd.NA).dropna().astype(str).unique()))),
            
            'inferred_pos2_grna': ('inferred_pos2_grna', lambda x: ','.join(sorted(x.replace('', pd.NA).dropna().astype(str).unique()))),
            'inferred_pos2_target': ('inferred_pos2_target', lambda x: ','.join(sorted(x.replace('', pd.NA).dropna().astype(str).unique()))),
            
            'n_combos': ('n_combos', 'max'),  # same for all rows of a cell
        }
        
        cell_grna_info = assign_df.groupby('bc_wells').agg(**agg_dict).reset_index().set_index('bc_wells')
        
        str_cols = [c for c in cell_grna_info.columns if c not in ['pos1_umi', 'pos2_umi', 'n_combos']]
        for col in str_cols:
            adata.obs[col] = cell_grna_info.reindex(adata.obs.index)[col]
            adata.obs.loc[adata.obs[col] == '', col] = 'unassigned'
            adata.obs[col] = adata.obs[col].fillna('unassigned')
            
        for col in ['pos1_umi', 'pos2_umi', 'n_combos']:
            adata.obs[col] = cell_grna_info.reindex(adata.obs.index)[col]
            adata.obs[col] = adata.obs[col].fillna(0).astype(int)
            
        adata.obs['target_gene'] = adata.obs['target_combo'].apply(
            lambda x: x.split(',')[0] if x != 'unassigned' else 'unassigned'
        )
        adata.obs['has_grna'] = adata.obs['target_gene'] != 'unassigned'
        
    else:
        # Legacy single-guide aggregation
        agg_dict = {
            'grna_names': ('grna_name', lambda x: ','.join(sorted(x.unique()))),
            'target_genes': ('target_gene', lambda x: ','.join(sorted(x.unique()))),
            'total_grna_umi': ('umi_count', 'sum'),
            'n_guides': ('grna_name', 'nunique'),
            'assignment_type': ('assignment_type', 'first'),
        }
        if 'ibar_seq' in assign_df.columns:
            agg_dict['ibar_seqs'] = ('ibar_seq', lambda x: ','.join(sorted(x.unique())))
        elif 'ibar_seq_corrected' in assign_df.columns:
            agg_dict['ibar_seqs'] = ('ibar_seq_corrected', lambda x: ','.join(sorted(x.unique())))

        pair_cols = ['pair_id', 'target_combo', 'target_version', 'gi_pair', 'gi_version', 'position']
        for col in pair_cols:
            if col in assign_df.columns:
                agg_dict[f'{col}s'] = (col, lambda x, c=col: ','.join(sorted(str(v) for v in x.dropna().unique())))

        cell_grna_info = assign_df.groupby('bc_wells').agg(**agg_dict).reset_index().set_index('bc_wells')

        str_cols = [c for c in cell_grna_info.columns if c not in ['total_grna_umi', 'n_guides']]
        for col in str_cols:
            adata.obs[col] = cell_grna_info.reindex(adata.obs.index)[col]
            adata.obs[col] = adata.obs[col].fillna('unassigned')

        for col in ['total_grna_umi', 'n_guides']:
            adata.obs[col] = cell_grna_info.reindex(adata.obs.index)[col]
            adata.obs[col] = adata.obs[col].fillna(0).astype(int)

        adata.obs['target_gene'] = adata.obs['target_genes'].apply(
            lambda x: x.split(',')[0] if x != 'unassigned' else 'unassigned'
        )
        adata.obs['has_grna'] = adata.obs['target_gene'] != 'unassigned'

    print(f"     Cells with gRNA: {adata.obs['has_grna'].sum():,} "
          f"({adata.obs['has_grna'].mean()*100:.1f}%)")
    print(f"     Target gene/combo distribution:")
    for gene, count in adata.obs['target_gene'].value_counts().head(20).items():
        pct = count / adata.n_obs * 100
        print(f"       {gene:>30s}: {count:>8,} ({pct:.1f}%)")

    # --- 6. Add gRNA UMI count matrix to obsm ---
    print("\n  6. Adding gRNA UMI matrix to obsm...")

    # Build a matrix aligned to adata.obs.index
    guide_bc_to_idx = {bc: i for i, bc in enumerate(guide_barcodes)}
    n_cells = adata.n_obs
    n_guides = len(guide_names)

    aligned_rows = []
    aligned_data = []
    aligned_cols = []

    for i, bc in enumerate(adata.obs.index):
        if bc in guide_bc_to_idx:
            j = guide_bc_to_idx[bc]
            row = guide_matrix[j]
            if sparse.issparse(row):
                row = row.toarray().flatten()
            else:
                row = np.asarray(row).flatten()
            nz = np.nonzero(row)[0]
            for col_idx in nz:
                aligned_rows.append(i)
                aligned_cols.append(col_idx)
                aligned_data.append(row[col_idx])

    aligned_matrix = sparse.csr_matrix(
        (aligned_data, (aligned_rows, aligned_cols)),
        shape=(n_cells, n_guides)
    )
    adata.obsm['gRNA_counts'] = aligned_matrix

    # Store guide names as unstructured data
    adata.uns['gRNA_names'] = guide_names

    # Binary assignment (threshold-based)
    min_umi = cfg['guide_assignment']['min_umi_threshold']
    adata.obsm['gRNA_assignment'] = (aligned_matrix >= min_umi).astype(np.int8)

    print(f"     gRNA_counts shape: {aligned_matrix.shape}")
    print(f"     Non-zero entries: {aligned_matrix.nnz:,}")

    # --- 7. Add metadata ---
    print("\n  7. Adding pipeline metadata...")
    adata.uns['ccperturbseq_pipeline'] = {
        'analysis_name': cfg['analysis_name'],
        'kit': cfg['splitpipe']['kit'],
        'chemistry': cfg['splitpipe']['chemistry'],
        'guide_assignment_threshold': min_umi,
        'max_lev_dist_grna': cfg['grna_matching']['max_lev_dist'],
        'max_lev_dist_ibar': cfg['grna_matching']['max_lev_dist_ibar'],
        'grna_lengths': cfg['grna_structure']['grna_lengths'],
        'ibar_length': cfg['grna_structure']['ibar_length'],
        'n_libraries': len(cfg['libraries']),
        'library_names': [lib['sublibrary_name'] for lib in cfg['libraries']],
    }

    # Library metadata (must convert int keys to strings for h5ad serialization)
    lib_info = pd.DataFrame(cfg['libraries'])
    lib_dict = lib_info.to_dict()
    adata.uns['library_info'] = {k: {str(ki): vi for ki, vi in v.items()} 
                                 for k, v in lib_dict.items()}

    # Guide annotation table (maps guide_name → design-level info)
    guide_annot_path = os.path.join(assign_dir, 'guide_annotations.csv')
    if os.path.exists(guide_annot_path):
        guide_annot = pd.read_csv(guide_annot_path)
        adata.uns['guide_annotations'] = guide_annot.to_dict(orient='list')
        print(f"     Stored guide_annotations ({len(guide_annot)} entries) in adata.uns")

    # --- 8. Save ---
    print("\n  8. Saving combined AnnData...")
    final_path = os.path.join(output_dir, 'combined_adata.h5ad')
    adata.write_h5ad(final_path)
    print(f"     Saved: {final_path}")
    print(f"     Shape: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    print(f"     obsm keys: {list(adata.obsm.keys())}")
    print(f"     obs columns: {list(adata.obs.columns)}")

    # --- Print final summary ---
    print(f"\n{'='*60}")
    print(f"  Combined AnnData Summary")
    print(f"{'='*60}")
    print(f"  GEX: {adata.n_vars:,} genes")
    print(f"  gRNA: {n_guides} guides")
    print(f"  Cells: {adata.n_obs:,}")
    print(f"  Cells with gRNA: {adata.obs['has_grna'].sum():,} "
          f"({adata.obs['has_grna'].mean()*100:.1f}%)")
    
    count_col = 'n_guides'
    if 'n_guides' not in adata.obs.columns and 'n_combos' in adata.obs.columns:
        count_col = 'n_combos'
    
    if count_col in adata.obs.columns:
        print(f"  Mean {count_col.replace('n_', '')}/cell (assigned): "
              f"{adata.obs.loc[adata.obs['has_grna'], count_col].mean():.2f}")
    
    print(f"\n  Output: {final_path}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()

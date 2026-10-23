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
import json
from pathlib import Path

import numpy as np
import pandas as pd
try:
    import scanpy as sc
except ImportError:  # Keep sparse-matrix helpers testable outside the HPC env.
    sc = None
import yaml
from scipy import io as sio
from scipy import sparse


def load_config(config_path):
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_design_metadata(output_dir, libraries):
    metadata = {}
    for library in libraries:
        name = library['sublibrary_name']
        design_path = Path(output_dir) / 'grna_matches' / f'{name}_design.yaml'
        if design_path.is_file():
            with design_path.open() as handle:
                contract = yaml.safe_load(handle)
            metadata[name] = json.dumps({key: contract.get(key) for key in
                                         ('schema_version', 'fingerprint', 'input_fingerprint', 'profiles', 'capabilities', 'resolved_orientation')}, sort_keys=True)
    return metadata


def _guide_product_paths(output_dir, source, gex_namespace='gex_grna_identity'):
    """Return the matrix/assignment paths for one Step 4 evidence product."""
    assign_dir = Path(output_dir) / 'grna_assignment'
    if source == 'dedicated':
        return {
            'assignment': assign_dir / 'grna_cell_assignment.csv',
            'matrix': assign_dir / 'grna_umi_matrix.mtx',
            'barcodes': assign_dir / 'grna_matrix_barcodes.txt',
            'guides': assign_dir / 'grna_matrix_guides.txt',
        }
    if source == 'gex':
        gex_dir = Path(output_dir) / gex_namespace / 'grna_assignment'
        return {
            'assignment': gex_dir / 'gex_grna_assignments.csv.gz',
            'matrix': gex_dir / 'gex_grna_umi_matrix.mtx',
            'barcodes': gex_dir / 'gex_grna_matrix_barcodes.txt',
            'guides': gex_dir / 'gex_grna_matrix_guides.txt',
        }
    if source == 'combined':
        return {
            'assignment': assign_dir / 'grna_cell_assignment_combined.csv',
            'matrix': assign_dir / 'grna_umi_matrix_combined.mtx',
            'barcodes': assign_dir / 'grna_matrix_combined_barcodes.txt',
            'guides': assign_dir / 'grna_matrix_combined_guides.txt',
        }
    raise ValueError(f"Unknown guide evidence source: {source}")


def _read_guide_product(paths, required=True):
    """Load a Step 4 cells-by-guides sparse matrix and its labels."""
    required_paths = [paths['matrix'], paths['barcodes'], paths['guides']]
    missing = [str(path) for path in required_paths if not Path(path).exists()]
    if missing:
        if not required:
            return None
        raise FileNotFoundError(
            "Incomplete Step 4 guide matrix product:\n  " + "\n  ".join(missing)
        )
    matrix = sio.mmread(paths['matrix']).tocsr()
    barcodes = Path(paths['barcodes']).read_text().splitlines()
    guides = Path(paths['guides']).read_text().splitlines()
    if matrix.shape != (len(barcodes), len(guides)):
        raise ValueError(
            f"Matrix/label shape mismatch for {paths['matrix']}: "
            f"{matrix.shape} versus {len(barcodes)} x {len(guides)}"
        )
    return matrix, barcodes, guides


def _align_guide_matrix(matrix, source_barcodes, target_barcodes):
    """Sparse row reindexing without densifying individual guide rows."""
    target_index = {barcode: idx for idx, barcode in enumerate(target_barcodes)}
    source_to_target = np.fromiter(
        (target_index.get(barcode, -1) for barcode in source_barcodes),
        dtype=np.int64,
        count=len(source_barcodes),
    )
    coo = matrix.tocoo()
    mapped_rows = source_to_target[coo.row]
    keep = mapped_rows >= 0
    aligned = sparse.coo_matrix(
        (coo.data[keep], (mapped_rows[keep], coo.col[keep])),
        shape=(len(target_barcodes), matrix.shape[1]),
    ).tocsr()
    aligned.sum_duplicates()
    return aligned


def _attach_guide_matrix(adata, matrix, barcodes, guides, source, min_umi):
    """Attach one source matrix under stable source-specific AnnData keys."""
    aligned = _align_guide_matrix(matrix, barcodes, list(adata.obs.index))
    suffix = '' if source == 'final' else f'_{source}'
    adata.obsm[f'gRNA_counts{suffix}'] = aligned
    adata.obsm[f'gRNA_assignment{suffix}'] = (
        aligned >= min_umi
    ).astype(np.int8)
    adata.uns[f'gRNA_names{suffix}'] = list(guides)
    return aligned


def _assignment_comparison_table(assign_df):
    """Compact per-cell metadata used to compare dedicated and union calls."""
    if len(assign_df) == 0 or 'bc_wells' not in assign_df:
        return pd.DataFrame()

    dual_mode = 'detection_status' in assign_df.columns
    aggregations = {}
    if dual_mode:
        first_cols = ('combo_type',)
        string_cols = (
            'detection_status', 'pair_id', 'target_combo', 'gi_pair',
            'pos1_grna', 'pos2_grna',
            'inferred_pos1_grna', 'inferred_pos2_grna',
        )
        sum_cols = ('pos1_umi', 'pos2_umi')
        max_cols = ('n_combos',)
    else:
        first_cols = ('combo_type',)
        string_cols = ('grna_name', 'target_gene', 'pair_id', 'target_combo')
        sum_cols = ('umi_count',)
        max_cols = ()

    for col in first_cols:
        if col in assign_df:
            aggregations[col] = 'first'
    for col in string_cols:
        if col in assign_df:
            aggregations[col] = lambda values: ','.join(
                sorted(values.dropna().astype(str).unique())
            )
    for col in sum_cols:
        if col in assign_df:
            aggregations[col] = 'sum'
    for col in max_cols:
        if col in assign_df:
            aggregations[col] = 'max'
    return assign_df.groupby('bc_wells').agg(aggregations)


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

    if sc is None:
        raise ImportError(
            "scanpy is required to build combined_adata.h5ad. "
            "Load the pipeline Python environment before running Step 5."
        )

    cfg = load_config(args.config)
    alignment_cfg = cfg.get('gex_alignment', {}) or {}
    try:
        from scripts.gex_alignment import backend_options, uses_shared_counter
    except ModuleNotFoundError:
        from gex_alignment import backend_options, uses_shared_counter
    alignment_options = backend_options(cfg)
    if uses_shared_counter(alignment_options):
        try:
            from scripts.gex_alignment import validate_combined
        except ModuleNotFoundError:
            from gex_alignment import validate_combined
        validate_combined(cfg, args.config)
    else:
        try:
            from scripts.gex_alignment import validate_star_combined_policy
        except ModuleNotFoundError:
            from gex_alignment import validate_star_combined_policy
        validate_star_combined_policy(cfg, args.config)
    output_dir = cfg['paths']['output_dir']
    assign_dir = os.path.join(output_dir, 'grna_assignment')
    decontam = cfg.get('gex_decontamination', {})
    step1c_cfg = decontam.get('step1c', {})
    combine_cfg = decontam.get('combine_with_grna_arm', {})
    gex_namespace = step1c_cfg.get(
        'output_namespace', 'gex_grna_identity'
    )
    combined_calls = bool(combine_cfg.get('enabled', False))
    final_guide_source = 'combined' if combined_calls else 'dedicated'
    final_product_paths = _guide_product_paths(
        output_dir, final_guide_source, gex_namespace
    )

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
    assign_path = final_product_paths['assignment']
    if not assign_path.exists():
        raise FileNotFoundError(
            f"Assignment file not found: {assign_path}\n"
            f"Run Step 4 {final_guide_source} assignment first."
        )

    assign_df = pd.read_csv(assign_path)
    print(f"     Final guide-call source: {final_guide_source}")
    print(f"     Loaded {len(assign_df):,} assignments for "
          f"{assign_df['bc_wells'].nunique():,} cells")

    # --- 3. Load gRNA UMI matrix ---
    print("\n  3. Loading gRNA UMI matrix...")
    guide_matrix, guide_barcodes, guide_names = _read_guide_product(
        final_product_paths, required=True
    )
    n_guides = len(guide_names)
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
    is_dual = 'detection_status' in assign_df.columns

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
        for column in ('pair_resolution_basis', 'pos1_ibar_status', 'pos2_ibar_status', 'unpositioned_grna'):
            if column in assign_df.columns:
                agg_dict[column] = (column, lambda values: ','.join(sorted(set(values.fillna('').astype(str)))))
        if 'pair_candidates' in assign_df.columns:
            agg_dict['pair_candidates'] = ('pair_candidates', lambda values: json.dumps(sorted({candidate for value in values.fillna('[]') for candidate in json.loads(value)})))
        if 'unpositioned_umi' in assign_df.columns:
            agg_dict['unpositioned_umi'] = ('unpositioned_umi', 'sum')
        
        cell_grna_info = assign_df.groupby('bc_wells').agg(**agg_dict).reset_index().set_index('bc_wells')
        
        numeric_cols = [column for column in ('pos1_umi', 'pos2_umi', 'n_combos', 'unpositioned_umi') if column in cell_grna_info]
        str_cols = [c for c in cell_grna_info.columns if c not in numeric_cols]
        for col in str_cols:
            adata.obs[col] = cell_grna_info.reindex(adata.obs.index)[col]
            adata.obs.loc[adata.obs[col] == '', col] = 'unassigned'
            adata.obs[col] = adata.obs[col].fillna('unassigned')
            
        for col in numeric_cols:
            adata.obs[col] = cell_grna_info.reindex(adata.obs.index)[col]
            adata.obs[col] = adata.obs[col].fillna(0).astype(int)
            
        adata.obs['target_gene'] = adata.obs['target_combo'].apply(
            lambda x: x.split(',')[0] if x != 'unassigned' else 'unassigned'
        )
        adata.obs['has_grna'] = adata.obs['target_gene'] != 'unassigned'
        if 'pair_resolution_basis' in assign_df.columns:
            adata.obs['has_resolved_pair'] = adata.obs['pair_id'].ne('unassigned')
            adata.obs['has_grna'] = (adata.obs['pos1_umi'] + adata.obs['pos2_umi'] + adata.obs.get('unpositioned_umi', 0)) > 0
        
    else:
        # Single-guide metadata aggregation
        agg_dict = {
            'grna_names': ('grna_name', lambda x: ','.join(sorted(x.unique()))),
            'target_genes': ('target_gene', lambda x: ','.join(sorted(x.unique()))),
            'total_grna_umi': ('umi_count', 'sum'),
            'n_guides': ('grna_name', 'nunique'),
            'combo_type': ('combo_type', 'first'),
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

    if combined_calls:
        dedicated_assignment_path = _guide_product_paths(
            output_dir, 'dedicated', gex_namespace
        )['assignment']
        if not dedicated_assignment_path.exists():
            raise FileNotFoundError(
                "Combined calls are enabled, but the dedicated comparison "
                f"assignment is missing: {dedicated_assignment_path}"
            )
        dedicated_assign_df = pd.read_csv(dedicated_assignment_path)
        dedicated_comparison = _assignment_comparison_table(
            dedicated_assign_df
        ).add_prefix('dedicated_')
        adata.obs = adata.obs.join(dedicated_comparison, how='left')
        print(
            "     Preserved dedicated-only call metadata in "
            "adata.obs['dedicated_*'] columns"
        )

    print(f"     Cells with gRNA: {adata.obs['has_grna'].sum():,} "
          f"({adata.obs['has_grna'].mean()*100:.1f}%)")
    print(f"     Target gene/combo distribution:")
    for gene, count in adata.obs['target_gene'].value_counts().head(20).items():
        pct = count / adata.n_obs * 100
        print(f"       {gene:>30s}: {count:>8,} ({pct:.1f}%)")

    # --- 6. Add gRNA UMI count matrix to obsm ---
    print("\n  6. Adding gRNA UMI matrix to obsm...")

    min_umi = cfg['guide_assignment']['min_umi_threshold']
    aligned_matrix = _attach_guide_matrix(
        adata,
        guide_matrix,
        guide_barcodes,
        guide_names,
        source='final',
        min_umi=min_umi,
    )

    # Preserve inspectable source matrices. Combined calls remain the primary
    # unsuffixed keys only when explicitly enabled in config.
    optional_sources = []
    if combined_calls:
        optional_sources.extend(('dedicated', 'gex', 'combined'))
    elif step1c_cfg.get('enabled', False):
        optional_sources.append('gex')
    for source in optional_sources:
        product = _read_guide_product(
            _guide_product_paths(output_dir, source, gex_namespace),
            required=(source != 'gex' or combined_calls),
        )
        if product is None:
            print(f"     WARNING: Optional {source} guide matrix not found")
            continue
        source_matrix, source_barcodes, source_guides = product
        _attach_guide_matrix(
            adata,
            source_matrix,
            source_barcodes,
            source_guides,
            source=source,
            min_umi=min_umi,
        )

    print(f"     gRNA_counts shape: {aligned_matrix.shape}")
    print(f"     Non-zero entries: {aligned_matrix.nnz:,}")

    # --- 7. Add metadata ---
    print("\n  7. Adding pipeline metadata...")
    adata.uns['ccperturbseq_pipeline'] = {
        'analysis_name': cfg['analysis_name'],
        'gex_alignment_backend': alignment_cfg.get('backend', 'star'),
        'gex_umi_correction': alignment_cfg.get('umi_correction', 'parse'),
        'gex_gene_assignment': alignment_options['gene_assignment'],
        'gex_output_level': alignment_options['output_level'],
        'kit': cfg['splitpipe']['kit'],
        'chemistry': cfg['splitpipe']['chemistry'],
        'guide_assignment_threshold': min_umi,
        'max_lev_dist_grna': cfg['grna_matching']['max_lev_dist'],
        'max_lev_dist_ibar': cfg['grna_matching']['max_lev_dist_ibar'],
        'grna_lengths': cfg['grna_structure']['grna_lengths'],
        'ibar_length': cfg['grna_structure'].get('ibar_length', -1),
        'n_libraries': len(cfg['libraries']),
        'library_names': [lib['sublibrary_name'] for lib in cfg['libraries']],
        'guide_call_source': final_guide_source,
        'combined_source_used_for_final_calls': combined_calls,
        'gex_identity_namespace': gex_namespace,
    }

    # Library metadata (must convert int keys to strings for h5ad serialization)
    lib_info = pd.DataFrame(cfg['libraries'])
    lib_dict = lib_info.to_dict()
    adata.uns['library_info'] = {k: {str(ki): vi for ki, vi in v.items()} 
                                 for k, v in lib_dict.items()}
    try:
        from scripts.guide_design import enabled
    except ModuleNotFoundError:
        from guide_design import enabled
    design_metadata = load_design_metadata(output_dir, cfg['libraries']) if enabled(cfg) else {}
    if design_metadata:
        adata.uns['guide_design'] = design_metadata

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

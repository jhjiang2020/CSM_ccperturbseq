#!/usr/bin/env python3
"""Quick impact test: how many cells are affected by v0.3 stem ambiguity change?

Strategy:
  1. Load Step 2 barcodes (raw R1 sequences) for a subsample of reads
  2. Re-extract gRNA/iBAR with v0.3 code (which flags AMBIGUOUS stems)
  3. Resolve reads with BOTH v0.2 logic (tie-break to pos1) and v0.3 logic (skip fallback)
  4. Compare: which reads lose assignment? Aggregate by cell barcode.

Usage:
  python test_v03_stem_impact.py --config /path/to/config.yaml --sample C1 [--n-reads 500000]
"""
import os, sys, gzip, argparse, time
import numpy as np
from collections import Counter, defaultdict

# ============================================================================
# Import Step 3 matching functions
# ============================================================================
def _import_step3(scripts_dir):
    import importlib.util
    path = os.path.join(scripts_dir, '03_match_grna_ibar.py')
    spec = importlib.util.spec_from_file_location("match_grna_ibar_03", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

# ============================================================================
# Build lookups from reference
# ============================================================================
def build_lookups(ref_df, step3, max_lev=3, max_lev_ibar=3):
    import pandas as pd

    # Combined exact lookup: (grna_seq, ibar_seq) -> (grna_name, target_gene, pair_id, target_combo)
    combined_lookup = {}
    for _, row in ref_df.iterrows():
        key = (row['grna_seq'].upper(), row['ibar_seq'].upper())
        combined_lookup[key] = (row['grna_name'], row['target_gene'],
                                 row['pair_id'], row['target_combo'])

    # Shared gRNA sequences (appear at both pos1 and pos2)
    seq_pos = ref_df.groupby('grna_seq')['position'].apply(set)
    shared_seqs = set(seq_pos[seq_pos.apply(len) > 1].index)

    # grna_name -> list of iBARs
    grna_to_ibars = {}
    for _, row in ref_df.iterrows():
        grna_to_ibars.setdefault(row['grna_name'], []).append(row['ibar_seq'].upper())

    # grna_name -> position
    grna_name_to_position = dict(zip(ref_df['grna_name'], ref_df['position'].astype(int)))

    # grna_seq -> {position: grna_name} for shared
    grna_seq_to_pos_name = {}
    for seq in shared_seqs:
        sub = ref_df[ref_df['grna_seq'] == seq].drop_duplicates('grna_name')
        grna_seq_to_pos_name[seq] = dict(
            zip(sub['position'].astype(int), sub['grna_name'])
        )

    # gRNA matcher
    grna_matcher = step3.NeighborhoodMatcher(
        ref_df['grna_seq'].tolist(), ref_df['grna_name'].tolist(), max_dist=max_lev
    )

    # (grna_name, ibar) -> pair_id
    grna_ibar_to_pairid = {}
    for _, row in ref_df.iterrows():
        grna_ibar_to_pairid[(row['grna_name'], row['ibar_seq'].upper())] = row['pair_id']

    return {
        'combined_lookup': combined_lookup,
        'shared_seqs': shared_seqs,
        'grna_to_ibars': grna_to_ibars,
        'grna_name_to_position': grna_name_to_position,
        'grna_seq_to_pos_name': grna_seq_to_pos_name,
        'grna_matcher': grna_matcher,
        'grna_ibar_to_pairid': grna_ibar_to_pairid,
    }


def resolve_read(grna_result, grna_unique, grna_seq, ibar_seq, variant_int,
                 allow_variant_fallback, lookups, step3, max_lev_ibar):
    """Resolve a read to a pair_id using iBAR + optional variant fallback."""
    candidates = grna_result if isinstance(grna_result, list) else [grna_result]

    resolved = None
    if grna_unique:
        resolved = candidates[0]
    else:
        # iBAR disambiguation
        ibar_matches = []
        for cand in candidates:
            guide_ibars = lookups['grna_to_ibars'].get(cand, [])
            im, idist = step3._match_ibar_within_guide(ibar_seq, guide_ibars, max_lev_ibar)
            if im is not None:
                ibar_matches.append((cand, im, idist))

        if len(ibar_matches) == 1:
            resolved = ibar_matches[0][0]
        elif allow_variant_fallback and variant_int in (1, 2):
            pos_map = lookups['grna_seq_to_pos_name'].get(grna_seq, {})
            if variant_int in pos_map:
                resolved = pos_map[variant_int]

    if resolved is None:
        return None

    # iBAR verification
    guide_ibars = lookups['grna_to_ibars'].get(resolved, [])
    im, idist = step3._match_ibar_within_guide(ibar_seq, guide_ibars, max_lev_ibar)
    if im is None:
        return None
    return lookups['grna_ibar_to_pairid'].get((resolved, im))


def main():
    parser = argparse.ArgumentParser(description='v0.3 stem ambiguity impact test')
    parser.add_argument('--config', required=True, help='Path to config.yaml')
    parser.add_argument('--sample', required=True, help='Sample name (C1, C2, T1, T2)')
    parser.add_argument('--n-reads', type=int, default=0,
                        help='Number of reads to process (0 = all)')
    args = parser.parse_args()

    import yaml, pandas as pd

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    project_dir = os.path.dirname(os.path.abspath(args.config))
    scripts_dir = os.path.join(project_dir, 'scripts')
    step3 = _import_step3(scripts_dir)

    print(f"=== v0.3 Stem Ambiguity Impact Test ===")
    print(f"Sample: {args.sample}")
    print(f"Max reads: {'all' if args.n_reads == 0 else f'{args.n_reads:,}'}")
    print()

    # Load reference
    ref_path = cfg['paths']['grna_reference']
    ref_df = pd.read_csv(ref_path)
    max_lev_ibar = cfg.get('grna_matching', {}).get('max_lev_dist_ibar', 3)
    lookups = build_lookups(ref_df, step3, max_lev_ibar=max_lev_ibar)
    print(f"Shared gRNA sequences: {len(lookups['shared_seqs'])}")

    # Extraction parameters from grna_structure
    gs = cfg.get('grna_structure', {})
    tso_seq   = gs.get('fixed_seq_1', '')  # Split_TSO
    stem1     = gs.get('fixed_seq_2', '')  # CSM_stem_1
    stem2     = gs.get('fixed_seq_3', '')  # CSM_stem_2
    tracr1    = gs.get('fixed_seq_4', '')  # CSM_partial_tracr_1
    tracr2    = gs.get('fixed_seq_5', '')  # CSM_partial_tracr_2
    grna_lens = gs.get('grna_lengths', [20, 21])
    ibar_len  = gs.get('ibar_length', 7)
    tso_max   = gs.get('tso_max_edit_dist', 3)
    anch_max  = gs.get('fixed_seq_max_mismatches', 2)

    # Read Step 2 barcodes
    step2_path = os.path.join(project_dir, 'output', 'grna_barcodes',
                               f'{args.sample}_barcodes.tsv.gz')
    print(f"Reading: {step2_path}")
    t0 = time.time()

    counts = Counter()
    cell_impacts = defaultdict(lambda: {'lost': 0, 'kept': 0, 'total': 0})

    with gzip.open(step2_path, 'rt') as f:
        header = f.readline().strip().split('\t')
        for i, line in enumerate(f):
            if args.n_reads > 0 and i >= args.n_reads:
                break
            if i > 0 and i % 500000 == 0:
                elapsed = time.time() - t0
                print(f"  processed {i:,} reads ({elapsed:.0f}s)")

            fields = line.strip().split('\t')
            rec = dict(zip(header, fields))
            r1_seq = rec['r1_seq']
            cell_bc = f"{rec.get('bc_wells', '')}_{rec.get('bc_bci', '')}"
            counts['total'] += 1

            # Extract with v0.3 code (sets variant='AMBIGUOUS' for ties)
            result = step3.extract_grna_ibar(
                r1_seq, tso_seq, stem1, stem2,
                grna_lengths=grna_lens,
                ibar_length=ibar_len,
                tracr1=tracr1, tracr2=tracr2,
                tso_max_dist=tso_max,
                anchor_max_mm=anch_max,
            )
            if result is None:
                counts['extraction_failed'] += 1
                continue

            grna_seq = result['grna_seq'].upper()
            ibar_seq = result['ibar_seq'].upper()
            variant  = result.get('variant', 1)
            is_ambiguous = (variant == 'AMBIGUOUS')
            try:
                variant_int = int(variant)
            except (ValueError, TypeError):
                variant_int = None

            # Exact combined lookup — same for both versions
            if (grna_seq, ibar_seq) in lookups['combined_lookup']:
                counts['exact_combined'] += 1
                cell_impacts[cell_bc]['total'] += 1
                cell_impacts[cell_bc]['kept'] += 1
                continue

            # Fuzzy gRNA match
            grna_result, grna_dist, grna_unique = lookups['grna_matcher'].match(grna_seq)
            if grna_result is None:
                counts['grna_no_match'] += 1
                continue

            if not isinstance(grna_result, list):
                counts['unique_grna'] += 1
                cell_impacts[cell_bc]['total'] += 1
                cell_impacts[cell_bc]['kept'] += 1
                continue

            counts['shared_grna'] += 1
            if is_ambiguous:
                counts['shared_ambig_stem'] += 1

            # v0.2: always allow variant fallback (tie-break to 1)
            v02_variant = variant_int if variant_int else 1  # v0.2 set ties to 1
            pid_v02 = resolve_read(grna_result, grna_unique, grna_seq, ibar_seq,
                                   v02_variant, allow_variant_fallback=True,
                                   lookups=lookups, step3=step3, max_lev_ibar=max_lev_ibar)

            # v0.3: skip variant fallback if ambiguous
            pid_v03 = resolve_read(grna_result, grna_unique, grna_seq, ibar_seq,
                                   variant_int, allow_variant_fallback=not is_ambiguous,
                                   lookups=lookups, step3=step3, max_lev_ibar=max_lev_ibar)

            cell_impacts[cell_bc]['total'] += 1
            if pid_v02 == pid_v03:
                counts['same'] += 1
                if pid_v02 is not None:
                    cell_impacts[cell_bc]['kept'] += 1
            else:
                counts['different'] += 1
                if pid_v02 is not None and pid_v03 is None:
                    counts['lost'] += 1
                    cell_impacts[cell_bc]['lost'] += 1
                    if is_ambiguous:
                        counts['lost_ambig_stem'] += 1

    elapsed = time.time() - t0
    n_total = counts['total']
    n_shared = counts['shared_grna']

    print(f"\nDone in {elapsed:.0f}s")
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total reads processed:  {n_total:,}")
    print(f"Extraction failed:      {counts['extraction_failed']:,}")
    print(f"Exact combined:         {counts['exact_combined']:,} (unaffected)")
    print(f"Unique gRNA:            {counts['unique_grna']:,} (unaffected)")
    print(f"gRNA no match:          {counts['grna_no_match']:,}")
    print(f"Shared gRNA:            {n_shared:,}")
    if n_shared > 0:
        print(f"  Ambiguous stem:       {counts['shared_ambig_stem']:,} ({counts['shared_ambig_stem']/n_shared*100:.1f}%)")
        print(f"  Same v0.2/v0.3:       {counts['same']:,}")
        print(f"  Different:            {counts['different']:,}")
        print(f"    Lost assignment:    {counts['lost']:,}")
        print(f"      Due to ambig:     {counts['lost_ambig_stem']:,}")

    # Cell-level impact
    cells_with_loss = {bc: v for bc, v in cell_impacts.items() if v['lost'] > 0}
    print(f"\nTotal cell barcodes:    {len(cell_impacts):,}")
    print(f"Cells losing reads:     {len(cells_with_loss):,}")
    if cells_with_loss:
        losses = [v['lost'] for v in cells_with_loss.values()]
        totals = [v['total'] for v in cells_with_loss.values()]
        fracs  = [v['lost']/max(v['total'],1) for v in cells_with_loss.values()]
        print(f"  Reads lost per cell:  median={np.median(losses):.0f}  max={max(losses)}")
        print(f"  Fraction lost:        median={np.median(fracs):.3f}  max={max(fracs):.3f}")
        print(f"  Cells losing >50%:    {sum(1 for f in fracs if f > 0.5):,}")

    # Save per-cell detail
    out_path = os.path.join(project_dir, 'output', 'qc_analysis',
                             f'v03_stem_impact_{args.sample}.csv')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    rows = []
    for bc, v in cell_impacts.items():
        rows.append({'cell_bc': bc, 'total': v['total'], 'kept': v['kept'], 'lost': v['lost']})
    pd.DataFrame(rows).to_csv(out_path, index=False)
    print(f"\nPer-cell detail saved: {out_path}")


if __name__ == '__main__':
    main()

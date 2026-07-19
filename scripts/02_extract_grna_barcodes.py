#!/usr/bin/env python3
"""
Step 2: Extract cell barcodes from gRNA library R2 reads.

Uses the Parse Biosciences split-pipe barcode correction logic to extract
cell barcodes from gRNA library R2 reads. This ensures barcode identifiers
match between GEX and gRNA libraries.

Output: TSV file per gRNA library with columns:
  read_id | bc_wells | bc_bci | umi | r1_seq | r1_qual

Usage:
  python scripts/02_extract_grna_barcodes.py --config config.yaml
  python scripts/02_extract_grna_barcodes.py --config config.yaml --max_reads 10000
"""

import argparse
import gzip
import os
import sys
import time
from collections import Counter
from pathlib import Path

import yaml


def load_config(config_path):
    """Load and return the YAML configuration."""
    with open(config_path) as f:
        return yaml.safe_load(f)


def resolve_output_root(output_dir, output_namespace=None):
    """Resolve an optional output namespace without changing default paths.

    Relative namespaces live below ``paths.output_dir``; absolute namespaces
    are used as-is.  With no namespace this returns the historical output root.
    """
    if not output_namespace:
        return Path(output_dir)
    namespace = Path(output_namespace).expanduser()
    return namespace if namespace.is_absolute() else Path(output_dir) / namespace


def setup_splitpipe_imports(splitpipe_dir):
    """Add split-pipe to path and import necessary modules."""
    sys.path.insert(0, splitpipe_dir)
    from splitpipe import bcutils, pb_const, kits
    from splitpipe.bcinfo import BCInfo

    try:
        from splitpipe import cy_utils
    except ImportError:
        cy_utils = None
        print("WARNING: cy_utils (Cython) not available; falling back to pure-Python extraction.")

    return bcutils, pb_const, BCInfo, cy_utils, kits


def extract_and_check_pure_python(seq, bc_starts, bc_ends):
    """Pure Python fallback for cy_utils.extract_and_check.

    Extracts polyN (UMI) and three barcode subsequences from R2 sequence.
    Returns: (polyN, rbc1, rbc2, rbc3, no_ambigs)
    """
    # Index 0 = polyN (UMI), 1 = bc1, 2 = bc2, 3 = bc3
    parts = []
    for i in range(len(bc_starts)):
        start = bc_starts[i]
        end = bc_ends[i]
        if start >= 0 and end >= 0 and end <= len(seq):
            parts.append(seq[start:end])
        else:
            parts.append("")

    # Check for ambiguous bases (N) in barcode regions (not polyN/UMI)
    no_ambigs = True
    for i in range(1, len(parts)):
        if 'N' in parts[i]:
            no_ambigs = False
            break

    polyN = parts[0] if len(parts) > 0 else ""
    rbc1 = parts[1] if len(parts) > 1 else ""
    rbc2 = parts[2] if len(parts) > 2 else ""
    rbc3 = parts[3] if len(parts) > 3 else ""

    return polyN, rbc1, rbc2, rbc3, no_ambigs


def find_fastq_files(fastq_dir, prefix):
    """Find R1 and R2 FASTQ files matching a prefix."""
    fastq_dir = Path(fastq_dir)
    r1_files = sorted(fastq_dir.glob(f"*{prefix}*_R1_*.fastq.gz"))
    r2_files = sorted(fastq_dir.glob(f"*{prefix}*_R2_*.fastq.gz"))

    if not r1_files:
        raise FileNotFoundError(f"No R1 FASTQ files found matching *{prefix}*_R1_*.fastq.gz in {fastq_dir}")
    if not r2_files:
        raise FileNotFoundError(f"No R2 FASTQ files found matching *{prefix}*_R2_*.fastq.gz in {fastq_dir}")

    return r1_files, r2_files


def read_fastq_pairs(fq1_path, fq2_path, max_reads=0):
    """Generator yielding (header, seq1, qual1, seq2, qual2) from paired FASTQ files."""
    open_fn = gzip.open if str(fq1_path).endswith('.gz') else open

    with open_fn(fq1_path, 'rt') as f1, open_fn(fq2_path, 'rt') as f2:
        n = 0
        while True:
            # Read 4 lines from each file
            h1 = f1.readline().strip()
            s1 = f1.readline().strip()
            p1 = f1.readline().strip()
            q1 = f1.readline().strip()

            h2 = f2.readline().strip()
            s2 = f2.readline().strip()
            p2 = f2.readline().strip()
            q2 = f2.readline().strip()

            if not h1 or not h2:
                break

            n += 1
            yield h1, s1, q1, s2, q2

            if max_reads > 0 and n >= max_reads:
                break


def sample_r2_reads(fq2_path, n_sample=200_000):
    """Sample R2 reads into a DataFrame with 'seq' and 'qual' columns.

    This DataFrame is required by BCInfo.get_perfect_bc_counts() to compute
    barcode correction thresholds.
    """
    import pandas as pd
    open_fn = gzip.open if str(fq2_path).endswith('.gz') else open
    seqs = []
    quals = []
    with open_fn(fq2_path, 'rt') as f:
        n = 0
        while True:
            h = f.readline().strip()
            s = f.readline().strip()
            p = f.readline().strip()
            q = f.readline().strip()
            if not h:
                break
            seqs.append(s)
            quals.append(q)
            n += 1
            if n >= n_sample:
                break
    df = pd.DataFrame({'seq': seqs, 'qual': quals})
    print(f"  Sampled {len(df)} R2 reads for barcode threshold estimation")
    return df


def process_grna_library(
    fq1_path, fq2_path, fout,
    bcinfo, bcutils_mod, cy_utils_mod,
    bc_max_edist=1, min_r1_len=20, min_r2_len=50,
    max_reads=0,
):
    """Process a single gRNA library: extract cell barcodes from R2, keep R1.

    Args:
        fq1_path: Path to gRNA R1 FASTQ (contains gRNA sequence)
        fq2_path: Path to gRNA R2 FASTQ (contains cell barcodes)
        fout: Open file handle to write output TSV
        bcinfo: BCInfo object with barcode data
        bcutils_mod: bcutils module
        cy_utils_mod: cy_utils module (or None)
        bc_max_edist: Max edit distance for barcode correction
        min_r1_len: Minimum R1 length
        min_r2_len: Minimum R2 length
        max_reads: Max reads to process (0 = all)
    Returns:
        dict with statistics
    """
    # Sample R2 reads to compute barcode correction thresholds
    fq_df = sample_r2_reads(fq2_path)
    # Get barcode correction data
    bc_dicts = bcinfo.get_bc_dict(0)  # 0 = all rounds
    counts, ct_thresh = bcinfo.get_perfect_bc_counts(fq_df)

    # Get barcode coordinate offsets from amplicon sequence
    bc_starts, bc_ends = bcinfo.get_amp_starts_ends()

    # Get mapping dictionaries
    bc1_seq_to_windzp = bcinfo.get_seq_to_windzp(1)
    bc2_seq_to_windzp = bcinfo.get_seq_to_windzp(2)
    bc3_seq_to_windzp = bcinfo.get_seq_to_windzp(3)
    bc1_seq_to_bci = bcinfo.get_seq_to_bci(1)
    bc2_seq_to_bci = bcinfo.get_seq_to_bci(2)
    bc3_seq_to_bci = bcinfo.get_seq_to_bci(3)

    # Set up BCOffsets for cy_utils
    bc_offs = None
    if cy_utils_mod is not None:
        try:
            bc_offs = cy_utils_mod.BCOffsets(bc_starts, bc_ends)
        except Exception:
            bc_offs = None

    stats = Counter()
    t_start = time.time()

    for header1, seq1, qual1, seq2, qual2 in read_fastq_pairs(fq1_path, fq2_path, max_reads):
        stats['total_reads'] += 1

        # Check minimum lengths
        if len(seq1) < min_r1_len or len(seq2) < min_r2_len:
            stats['reads_too_short'] += 1
            continue

        # Extract barcodes from R2
        if bc_offs is not None:
            polyN, rbc1, rbc2, rbc3, no_ambigs = cy_utils_mod.extract_and_check(seq2, bc_offs)
        else:
            polyN, rbc1, rbc2, rbc3, no_ambigs = extract_and_check_pure_python(
                seq2, bc_starts, bc_ends
            )

        if not no_ambigs:
            stats['reads_ambig_bc'] += 1
            continue

        # Correct barcodes
        bc1, bc2, bc3, ed_max, ed_sum = bcutils_mod.correct_barcodes(
            rbc1, rbc2, rbc3, counts, ct_thresh, bc_dicts, bc_max_edist
        )

        if ed_max == 0:
            stats['bc_edit_dist_0'] += 1
        elif ed_max == 1:
            stats['bc_edit_dist_1'] += 1
        elif ed_max == 2:
            stats['bc_edit_dist_2'] += 1
        else:
            stats['bc_edit_dist_NA'] += 1

        # Need all three barcodes
        if not (bc1 and bc2 and bc3):
            stats['reads_no_bc'] += 1
            continue

        stats['reads_valid_barcode'] += 1

        # Build cell identifier (same format as split-pipe)
        w1 = bc1_seq_to_windzp[bc1]
        w2 = bc2_seq_to_windzp[bc2]
        w3 = bc3_seq_to_windzp[bc3]
        bc_wells = f"{w1}_{w2}_{w3}"

        i1 = bc1_seq_to_bci[bc1]
        i2 = bc2_seq_to_bci[bc2]
        i3 = bc3_seq_to_bci[bc3]
        bc_bci = f"{i1}_{i2}_{i3}"

        # Extract read_id from header
        read_id = header1.split()[0].lstrip('@')

        # Write output: read_id, bc_wells, bc_bci, umi (polyN), R1 seq, R1 qual
        fout.write(f"{read_id}\t{bc_wells}\t{bc_bci}\t{polyN}\t{seq1}\t{qual1}\n")

        if stats['total_reads'] % 1_000_000 == 0:
            elapsed = time.time() - t_start
            rate = stats['total_reads'] / elapsed
            valid_pct = stats['reads_valid_barcode'] / stats['total_reads'] * 100
            print(f"  Processed {stats['total_reads']:,} reads "
                  f"({rate:,.0f} reads/s, {valid_pct:.1f}% valid barcodes)")

    return dict(stats)


def main():
    parser = argparse.ArgumentParser(description="Extract cell barcodes from gRNA library R2 reads")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument("--max_reads", type=int, default=0,
                        help="Max reads to process per library (0=all). Overrides config.")
    parser.add_argument("--library", type=str, default=None,
                        help="Process only this sublibrary name (e.g., 'C1'). Default: all.")
    parser.add_argument(
        "--output-namespace", default=None,
        help=("Output root for this source. Relative paths are below paths.output_dir; "
              "default keeps the historical output layout."),
    )
    args = parser.parse_args()

    cfg = load_config(args.config)

    # Override max_reads if specified on command line
    max_reads = args.max_reads if args.max_reads > 0 else cfg['barcode_extraction']['max_reads']

    # Set up split-pipe imports
    splitpipe_dir = cfg['paths']['splitpipe_dir']
    bcutils_mod, pb_const, BCInfo, cy_utils_mod, kits_mod = setup_splitpipe_imports(splitpipe_dir)

    # Get barcode round list for this kit+chemistry via kits API
    chemistry = cfg['splitpipe']['chemistry']
    kit = cfg['splitpipe']['kit']

    # Parse chemistry to integer (kits API expects int, e.g. 3 not "v3")
    chem_int = kits_mod.parse_chemistry(chemistry)
    if not chem_int:
        print(f"ERROR: Invalid chemistry '{chemistry}'")
        sys.exit(1)

    # Get barcode round list: e.g. [['bc1','n141_R1_v3_8'], ['bc2','v1'], ['bc3','R3_v3']]
    bc_rlist = kits_mod.kit_bc_set_list(kit, chem_int, as_tup=True)
    if not bc_rlist:
        print(f"ERROR: Could not resolve barcode sets for kit='{kit}', chemistry='{chemistry}'")
        sys.exit(1)

    # Get amplicon sequence for this chemistry
    amp_seq = kits_mod.get_chem_amp_seq(chem_int)
    print(f"Using chemistry {chemistry}, amplicon sequence: {amp_seq}")
    print(f"Barcode rounds: {bc_rlist}")

    # Initialize BCInfo directly
    bc_path = os.path.join(splitpipe_dir, 'splitpipe', 'barcodes')
    bcinfo = BCInfo(bc_rlist, path=bc_path, amp_seq=amp_seq)
    print(f"BCInfo initialized: {bcinfo}")

    # Need to calculate perfect BC counts for correction
    # This requires reading the R2 FASTQs first to get barcode frequency
    # We'll use a two-pass approach: first count, then extract

    libraries = cfg['libraries']
    selected_library = args.library
    if selected_library:
        libraries = [lib for lib in libraries if lib['sublibrary_name'] == selected_library]
        if not libraries:
            print(f"ERROR: Library '{selected_library}' not found in config")
            sys.exit(1)

    fastq_dir = cfg['paths']['fastq_dir']
    output_dir = str(resolve_output_root(
        cfg['paths']['output_dir'], args.output_namespace
    ))
    bc_max_edist = cfg['barcode_extraction']['bc_max_edit_dist']
    min_r1_len = cfg['barcode_extraction']['min_r1_length']
    min_r2_len = cfg['barcode_extraction']['min_r2_length']

    print(f"\n{'='*60}")
    print(f"  Step 2: gRNA Barcode Extraction")
    print(f"{'='*60}")
    print(f"  Libraries: {[lib['sublibrary_name'] for lib in libraries]}")
    print(f"  Max reads: {'all' if max_reads == 0 else max_reads}")
    print(f"  BC max edit dist: {bc_max_edist}")
    print(f"")

    all_stats = {}

    for lib in libraries:
        lib_name = lib['sublibrary_name']
        grna_prefix = lib['grna_prefix']

        print(f"\n--- Processing gRNA library: {lib_name} (prefix: {grna_prefix}) ---")

        r1_files, r2_files = find_fastq_files(fastq_dir, grna_prefix)
        print(f"  R1 files: {[f.name for f in r1_files]}")
        print(f"  R2 files: {[f.name for f in r2_files]}")

        output_path = os.path.join(
            output_dir, 'grna_barcodes', f"{lib_name}_barcodes.tsv.gz"
        )
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        lib_stats = Counter()
        reads_left = max_reads

        with gzip.open(output_path, 'wt') as fout:
            # Write header
            fout.write("read_id\tbc_wells\tbc_bci\tumi\tr1_seq\tr1_qual\n")

            # Process each FASTQ pair
            for fq1, fq2 in zip(r1_files, r2_files):
                print(f"\n  Processing: {fq1.name}")
                if max_reads > 0 and reads_left <= 0:
                    print("  Reached max_reads limit. Skipping remaining files.")
                    break

                stats = process_grna_library(
                    str(fq1), str(fq2), fout,
                    bcinfo, bcutils_mod, cy_utils_mod,
                    bc_max_edist=bc_max_edist,
                    min_r1_len=min_r1_len,
                    min_r2_len=min_r2_len,
                    max_reads=reads_left,
                )
                lib_stats.update(stats)
                if max_reads > 0:
                    reads_left -= stats.get('total_reads', 0)

        all_stats[lib_name] = dict(lib_stats)

        # Print summary
        total = lib_stats.get('total_reads', 0)
        valid = lib_stats.get('reads_valid_barcode', 0)
        pct = valid / total * 100 if total > 0 else 0
        print(f"\n  Summary for {lib_name} (Combined {len(r1_files)} files):")
        print(f"    Total reads:     {total:>12,}")
        print(f"    Valid barcodes:  {valid:>12,} ({pct:.1f}%)")
        print(f"    Too short:       {lib_stats.get('reads_too_short', 0):>12,}")
        print(f"    Ambiguous BC:    {lib_stats.get('reads_ambig_bc', 0):>12,}")
        print(f"    No BC match:     {lib_stats.get('reads_no_bc', 0):>12,}")
        print(f"    BC edit dist 0:  {lib_stats.get('bc_edit_dist_0', 0):>12,}")
        print(f"    BC edit dist 1:  {lib_stats.get('bc_edit_dist_1', 0):>12,}")
        print(f"    Output: {output_path}")

    # Save overall stats
    lib_label = selected_library or "all"
    stats_path = os.path.join(output_dir, 'grna_barcodes', f'{lib_label}_extraction_stats.yaml')
    os.makedirs(os.path.dirname(stats_path), exist_ok=True)
    with open(stats_path, 'w') as f:
        yaml.dump(all_stats, f, default_flow_style=False)
    print(f"\n  Stats saved to {stats_path}")

    print(f"\n{'='*60}")
    print(f"  Step 2 complete!")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()

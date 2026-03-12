import os
import gzip
import re
import csv
from collections import defaultdict, Counter
import pandas as pd


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))


def parse_fastq(r1_path, r2_path):
    """
    Generator that yields (header, r1_seq, r2_seq) pairs.
    Handles gzip if extension is .gz.
    """
    opener = gzip.open if str(r1_path).endswith('.gz') else open
    mode = 'rt'
    
    with opener(r1_path, mode) as f1, opener(r2_path, mode) as f2:
        while True:
            # FASTQ format: 4 lines per read
            try:
                h1 = f1.readline().strip()
                s1 = f1.readline().strip()
                _  = f1.readline()
                _  = f1.readline()
                
                h2 = f2.readline().strip()
                s2 = f2.readline().strip()
                _  = f2.readline()
                _  = f2.readline()
            except ValueError:
                break # End of file usually
            
            if not h1: break
            
            yield h1, s1, s2


def compile_search_patterns(ref_df, add_flank_ibars=True):
    """
    Compiles regex patterns for Spacers and iBARs.
    Returns:
        patterns: dict of re objects
        maps: dict of sequence -> set of (Gene, Version/ID)
    """
    s1_gene_map = defaultdict(set)
    b1_gene_map = defaultdict(set)
    s2_gene_map = defaultdict(set)
    b2_gene_map = defaultdict(set)

    for idx, row in ref_df.iterrows():
        pos = row['position']
        g_name = row['grna_name']
        s_seq = str(row['grna_seq']).strip().upper()
        b_seq = str(row['ibar_seq']).strip().upper() # already rc of design

        if pos == 1:
            s1 = s_seq
            b1 = b_seq
            if add_flank_ibars:
                b1 = 'GGA' + b1 + 'GAC'
            g1 = g_name
            if s1: s1_gene_map[s1].add(g1)
            if b1: b1_gene_map[b1].add(g1)
        elif pos == 2:
            s2_orig = s_seq
            b2_orig = b_seq
            if add_flank_ibars:
                b2_orig = 'GGA' + b2_orig + 'AAC'
            g2 = g_name
            
            s2 = reverse_complement(s2_orig)
            b2 = reverse_complement(b2_orig)
            
            if s2: s2_gene_map[s2].add(g2)
            if b2: b2_gene_map[b2].add(g2)

    def make_regex(seqs):
        if not seqs: return None
        sorted_seqs = sorted(list(seqs), key=len, reverse=True)
        pattern_str = "|".join(re.escape(s) for s in sorted_seqs)
        return re.compile(pattern_str)

    patterns = {
        'S1': make_regex(s1_gene_map.keys()),
        'B1': make_regex(b1_gene_map.keys()),
        'S2': make_regex(s2_gene_map.keys()),
        'B2': make_regex(b2_gene_map.keys())
    }
    
    maps = {
        'S1_gene': s1_gene_map, 'B1_gene': b1_gene_map,
        'S2_gene': s2_gene_map, 'B2_gene': b2_gene_map,
    }
    
    return patterns, maps


def get_best_match(seq, pattern, map_dict, start_idx=0, end_idx=None):
    """
    Finds match in sequence within specific range.
    Returns: (matched_seq, set_of_genes, start_pos) or (None, None, None)
    """
    if not pattern:
        return None, None, None
    
    search_region = seq[start_idx:end_idx]
    
    m = pattern.search(search_region)
    if m:
        matched_seq = m.group()
        genes = map_dict[matched_seq]
        abs_start = start_idx + m.start()
        return matched_seq, genes, abs_start
    return None, None, None


def resolve_gene(genes, tag):
    """
    Resolves a set of genes to a single string.
    """
    if not genes:
        return None
    if len(genes) == 1:
        return list(genes)[0]
    return f"Ambiguous({len(genes)})"


def load_design_and_build_patterns(ref_file, add_flank_ibars=True):
    """
    Loads formatted gRNA reference CSV, validates it, and compiles regex patterns.
    """
    print(f"Loading gRNA reference from {ref_file}...")
    ref_df = pd.read_csv(ref_file)
    
    required_columns = ['grna_name', 'grna_seq', 'ibar_seq', 'position', 'pair_id']
    missing_cols = [c for c in required_columns if c not in ref_df.columns]
    if missing_cols:
        raise ValueError(f"Reference file is missing required columns: {missing_cols}")
        
    patterns, maps = compile_search_patterns(ref_df, add_flank_ibars=add_flank_ibars)
    return ref_df, patterns, maps


def analyze_sample(r1, r2, patterns, maps, ref_df, output_prefix):
    valid_pairs = set()
    for pair_id, group in ref_df.groupby('pair_id'):
        p1_rows = group[group['position'] == 1]
        p2_rows = group[group['position'] == 2]
        if not p1_rows.empty and not p2_rows.empty:
            p1 = p1_rows.iloc[0]['grna_name']
            p2 = p2_rows.iloc[0]['grna_name']
            valid_pairs.add((p1, p2))
        
    stats = Counter()
    
    out_file = f"{output_prefix}_detailed.csv.gz"
    print(f"Writing detailed report to {out_file}...")
    
    R1_S_START, R1_S_END = 18, 43
    R1_B_START, R1_B_END = 54, 75
    R2_S_START, R2_S_END = 49, 74
    R2_B_START, R2_B_END = 19, 40
    
    with gzip.open(out_file, 'wt', newline='') as f:
        writer = csv.writer(f)
        header = [
            "ReadID", 
            "S1_Gene", "S1_Seq", "S1_Pos",
            "B1_Gene", "B1_Seq", "B1_Pos",
            "S2_Gene", "S2_Seq", "S2_Pos",
            "B2_Gene", "B2_Seq", "B2_Pos",
            "Annotation"
        ]
        writer.writerow(header)
        
        count = 0
        for head, seq1, seq2 in parse_fastq(r1, r2):
            count += 1
            if count % 200000 == 0:
                print(f"  Processed {count} reads...", end='\r')
            
            stats['Total Pairs'] += 1
            
            s1_seq, s1_genes, s1_pos = get_best_match(seq1, patterns['S1'], maps['S1_gene'], R1_S_START, R1_S_END)
            b1_seq, b1_genes, b1_pos = get_best_match(seq1, patterns['B1'], maps['B1_gene'], R1_B_START, R1_B_END)
            s1_res = resolve_gene(s1_genes, "S1")
            b1_res = resolve_gene(b1_genes, "B1")
            
            s2_seq, s2_genes, s2_pos = get_best_match(seq2, patterns['S2'], maps['S2_gene'], R2_S_START, R2_S_END)
            b2_seq, b2_genes, b2_pos = get_best_match(seq2, patterns['B2'], maps['B2_gene'], R2_B_START, R2_B_END)
            s2_res = resolve_gene(s2_genes, "S2")
            b2_res = resolve_gene(b2_genes, "B2")
            
            annotations = []
            
            # Position 1 checking
            final_g1 = None
            if s1_res and b1_res:
                if s1_res == b1_res:
                    final_g1 = s1_res
                else:
                    annotations.append("spacer_iBAR_mismatch_pos1")
                    final_g1 = s1_res
            elif s1_res:
                final_g1 = s1_res
            elif b1_res:
                final_g1 = b1_res
                
            # Position 2 checking
            final_g2 = None
            if s2_res and b2_res:
                if s2_res == b2_res:
                    final_g2 = s2_res
                else:
                    annotations.append("spacer_iBAR_mismatch_pos2")
                    final_g2 = s2_res
            elif s2_res:
                final_g2 = s2_res
            elif b2_res:
                final_g2 = b2_res
                
            # Pair checking
            if final_g1 and final_g2:
                if "Ambiguous" in final_g1 or "Ambiguous" in final_g2:
                    annotations.append("Ambiguous_Gene_Mapping")
                elif (final_g1, final_g2) in valid_pairs:
                    pass
                else:
                    annotations.append("gene_pair_mismatch")
            else:
                if not final_g1: annotations.append("Unmapped_Pos1")
                if not final_g2: annotations.append("Unmapped_Pos2")
            
            annot_str = ";".join(annotations) if annotations else "Match"
            
            if annot_str == "Match":
                stats['Perfect Match'] += 1
            else:
                for a in annotations:
                    stats[a] += 1
            
            if s1_res or b1_res or s2_res or b2_res:
                writer.writerow([
                    head,
                    s1_res, s1_seq, s1_pos,
                    b1_res, b1_seq, b1_pos,
                    s2_res, s2_seq, s2_pos,
                    b2_res, b2_seq, b2_pos,
                    annot_str
                ])
            
    print(f"\nFinished analyzing.")
    return stats


def analyze_recombination(r1, r2, patterns, maps, ref_df, sample_name):
    print(f"--- Analyzing Recombination (iBAR linkage) for {sample_name} ---")
    
    valid_pairs = set()
    for pair_id, group in ref_df.groupby('pair_id'):
        p1_rows = group[group['position'] == 1]
        p2_rows = group[group['position'] == 2]
        if not p1_rows.empty and not p2_rows.empty:
            p1 = p1_rows.iloc[0]['grna_name']
            p2 = p2_rows.iloc[0]['grna_name']
            valid_pairs.add((p1, p2))
        
    R1_B_START, R1_B_END = 54, 75
    R2_B_START, R2_B_END = 19, 40
    
    pair_counts = Counter()
    total_reads = 0
    mapped_reads = 0
    
    for head, seq1, seq2 in parse_fastq(r1, r2):
        total_reads += 1
        if total_reads % 200000 == 0:
            print(f"  Processed {total_reads} reads...", end='\r')
            
        _, b1_genes, _ = get_best_match(seq1, patterns['B1'], maps['B1_gene'], R1_B_START, R1_B_END)
        _, b2_genes, _ = get_best_match(seq2, patterns['B2'], maps['B2_gene'], R2_B_START, R2_B_END)
        
        b1_res = resolve_gene(b1_genes, "B1")
        b2_res = resolve_gene(b2_genes, "B2")
        
        if b1_res and b2_res and "Ambiguous" not in b1_res and "Ambiguous" not in b2_res:
            pair_counts[(b1_res, b2_res)] += 1
            mapped_reads += 1
            
    print(f"\n  Finished processing. Total reads: {total_reads}, Mapped pairs: {mapped_reads}")
    
    data = []
    for (b1, b2), count in pair_counts.items():
        is_valid = (b1, b2) in valid_pairs
        data.append({
            'iBAR1_Gene': b1,
            'iBAR2_Gene': b2,
            'Count': count,
            'Is_Designed_Pair': is_valid
        })
        
    df = pd.DataFrame(data)
    if df.empty:
        print("No valid iBAR pairs found.")
        return df, {}
        
    df['Frequency'] = df['Count'] / df['Count'].sum()
    df = df.sort_values('Count', ascending=False).reset_index(drop=True)
    
    total_mapped = df['Count'].sum()
    valid_count = df[df['Is_Designed_Pair']]['Count'].sum()
    recomb_count = df[~df['Is_Designed_Pair']]['Count'].sum()
    recomb_rate = recomb_count / total_mapped if total_mapped > 0 else 0
    
    summary = {
        'Sample': sample_name,
        'Total_Mapped_Pairs': total_mapped,
        'Valid_Pairs_Count': valid_count,
        'Recombinant_Pairs_Count': recomb_count,
        'Recombination_Rate': recomb_rate
    }
    
    return df, summary

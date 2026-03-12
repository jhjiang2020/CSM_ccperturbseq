#!/usr/bin/env python3
"""
Post-Step 3 QC: Diagnose extraction-failed reads.

This script samples reads from Step 2 output (barcode TSV) and classifies
why Step 3 failed to extract gRNA/iBAR from them.

Failure modes diagnosed:
  1. TSO not found     — no TSO match within edit distance threshold
  2. Read too short    — R1 shorter than minimum structure length
  3. Stem mismatch     — TSO found but neither stem variant matches
  4. Quality issues    — high N-content, poly-G runs, etc.

For TSO-not-found reads, the script examines:
  - Best possible TSO alignment (at higher edit distances)
  - Sequence composition (GC%, N-content)
  - Whether reads contain partial gRNA sequences (off-target vs degraded)
  - Kmer signature analysis

Usage:
  python scripts/qc_step3_failures.py --config config.yaml --library C1
  python scripts/qc_step3_failures.py --config config.yaml --library C1 --n_sample 50000
"""

import argparse
import gzip
import os
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

# Try to import edlib
try:
    import edlib
    HAS_EDLIB = True
except ImportError:
    HAS_EDLIB = False


def find_tso_extended(r1_seq, tso_seq, max_dists=[3, 5, 8, 15]):
    """Find TSO with progressively looser thresholds to classify failures."""
    if not HAS_EDLIB:
        return None, -1, -1

    target = r1_seq[:80]  # Wider search window
    for md in max_dists:
        result = edlib.align(tso_seq, target, mode="HW", k=md, task="locations")
        if result['editDistance'] >= 0 and result['locations']:
            loc = result['locations'][0]
            return md, result['editDistance'], loc[1] + 1
    return None, -1, -1


def hamming_distance(s1, s2):
    """Hamming distance with length mismatch handling."""
    d = 0
    for c1, c2 in zip(s1, s2):
        if c1 != c2:
            d += 1
    return d + abs(len(s1) - len(s2))


def check_sequence_quality(seq):
    """Assess sequence quality metrics."""
    n_count = seq.count('N')
    gc_count = seq.count('G') + seq.count('C')
    gc_frac = gc_count / len(seq) if seq else 0

    # Check for poly-runs
    max_poly = 0
    curr_run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            curr_run += 1
            max_poly = max(max_poly, curr_run)
        else:
            curr_run = 1

    # Poly-G: find the START position of the first >=10bp poly-G run
    # NovaSeq poly-G typically appears at the 3' end (tail) of reads.
    # Only relevant to extraction failure if it starts within the gRNA
    # cassette region (~first 85bp).
    poly_g_start = -1
    g_run = 0
    for i in range(len(seq)):
        if seq[i] == 'G':
            g_run += 1
            if g_run >= 10 and poly_g_start < 0:
                poly_g_start = i - g_run + 1
        else:
            g_run = 0

    return {
        'length': len(seq),
        'n_count': n_count,
        'n_frac': n_count / len(seq) if seq else 0,
        'gc_frac': gc_frac,
        'max_poly_run': max_poly,
        'poly_g_start': poly_g_start,  # -1 = no poly-G >=10bp
    }


def search_grna_anywhere(r1_seq, grna_seqs, max_dist=2):
    """Check if any gRNA sequence exists anywhere in R1 (not just at expected position)."""
    if not HAS_EDLIB:
        # Simple substring search for exact matches
        for gseq in grna_seqs:
            if gseq in r1_seq:
                return True, 0
        return False, -1

    for gseq in grna_seqs:
        result = edlib.align(gseq, r1_seq, mode="HW", k=max_dist, task="locations")
        if result['editDistance'] >= 0:
            return True, result['editDistance']
    return False, -1


def classify_failed_read(r1_seq, tso_seq, stem1, stem2, grna_lengths, ibar_length,
                          grna_seqs, anchor_max_mm=2, tso_threshold=3):
    """Classify why a read failed extraction.

    Classification priority:
      1. Read too short
      2. High N content
      3. TSO analysis (not found / degraded)
      4. Stem analysis (mismatch / truncation)
      5. Poly-G only flagged as primary if within cassette region

    Returns a classification dict.
    """
    # Cassette region: offset (~7) + TSO (30) + gRNA (max 21) + stem (19) + iBAR (7) = ~84bp
    cassette_end = len(tso_seq) + max(grna_lengths) + len(stem1) + ibar_length + 10
    min_structure_len = len(tso_seq) + min(grna_lengths) + len(stem1) + ibar_length

    result = {
        'failure_class': 'unknown',
        'failure_detail': '',
        'r1_length': len(r1_seq),
    }

    # Quality check (metadata — not used for primary classification yet)
    qc = check_sequence_quality(r1_seq)
    result.update(qc)

    # 1. Read too short for any structure
    if len(r1_seq) < min_structure_len:
        result['failure_class'] = 'too_short'
        result['failure_detail'] = f'len={len(r1_seq)}, need>={min_structure_len}'
        return result

    # 2. High N content
    if qc['n_frac'] > 0.1:
        result['failure_class'] = 'high_N_content'
        result['failure_detail'] = f"N_frac={qc['n_frac']:.2f}"
        return result

    # 3. Poly-G WITHIN the cassette region (first ~85bp) — a real structural problem
    if qc['poly_g_start'] >= 0 and qc['poly_g_start'] < cassette_end:
        result['failure_class'] = 'poly_G_in_cassette'
        result['failure_detail'] = f"poly-G starts at bp {qc['poly_g_start']} (cassette ends ~{cassette_end}bp)"
        return result

    # 4. TSO analysis — try extended search
    threshold_found, best_ed, tso_end = find_tso_extended(r1_seq, tso_seq)

    if threshold_found is None or best_ed < 0:
        # No TSO at all, even at very high edit distance
        # Check if gRNA exists somewhere in the read (off-target carrier)
        grna_found, grna_ed = search_grna_anywhere(r1_seq, grna_seqs)
        if grna_found:
            result['failure_class'] = 'tso_absent_grna_present'
            result['failure_detail'] = f'gRNA found at ed={grna_ed} but no TSO'
        else:
            result['failure_class'] = 'off_target'
            result['failure_detail'] = 'no TSO or gRNA signature'
        return result

    if best_ed > tso_threshold:
        # TSO found but beyond the threshold
        result['failure_class'] = 'tso_degraded'
        result['failure_detail'] = f'best_TSO_ed={best_ed}, threshold={tso_threshold}'
        result['tso_edit_dist'] = best_ed
        result['tso_end'] = tso_end
        return result

    # 5. TSO found within threshold — stem must have failed
    # Try all gRNA lengths
    stem_dists = []
    for grna_len in grna_lengths:
        grna_end = tso_end + grna_len
        stem_end = grna_end + len(stem1)
        if stem_end <= len(r1_seq):
            stem_obs = r1_seq[grna_end:stem_end]
            d1 = hamming_distance(stem_obs, stem1)
            d2 = hamming_distance(stem_obs, stem2)
            stem_dists.append((grna_len, min(d1, d2), d1, d2, stem_obs))

    if stem_dists:
        best = min(stem_dists, key=lambda x: x[1])
        if best[1] <= anchor_max_mm:
            # Should have been extracted — maybe read truncated before iBAR
            ibar_end_pos = tso_end + best[0] + len(stem1) + ibar_length
            if ibar_end_pos > len(r1_seq):
                result['failure_class'] = 'truncated_before_ibar'
                result['failure_detail'] = f'need {ibar_end_pos}bp, have {len(r1_seq)}bp'
            else:
                result['failure_class'] = 'stem_ok_but_failed'
                result['failure_detail'] = f'unexpected: stem_d={best[1]}'
        else:
            result['failure_class'] = 'stem_mismatch'
            result['failure_detail'] = (
                f'best_stem_d={best[1]}, threshold={anchor_max_mm}, '
                f'stem_obs={best[4]}'
            )
            # Annotate if poly-G is in the tail (not the cause)
            if qc['poly_g_start'] >= cassette_end:
                result['failure_detail'] += f' (poly-G tail at bp {qc["poly_g_start"]}, NOT causal)'
    else:
        result['failure_class'] = 'too_short_after_tso'
        result['failure_detail'] = f'r1_len={len(r1_seq)}, tso_end={tso_end}'

    return result


def main():
    parser = argparse.ArgumentParser(
        description="QC diagnostic for Step 3 extraction failures"
    )
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument("--library", required=True, help="Sublibrary name")
    parser.add_argument("--n_sample", type=int, default=100_000,
                        help="Number of reads to sample for QC (default: 100000)")
    parser.add_argument("--output_dir", type=str, default=None,
                        help="Output directory (default: <output_dir>/qc/)")
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    gs = cfg['grna_structure']
    tso_seq = gs['fixed_seq_1']
    stem1 = gs['fixed_seq_2']
    stem2 = gs['fixed_seq_3']
    grna_lengths = gs['grna_lengths']
    ibar_length = gs['ibar_length']
    tso_threshold = gs.get('tso_max_edit_dist', 3)
    anchor_max_mm = gs.get('fixed_seq_max_mismatches', 2)

    # Load reference gRNA sequences for off-target detection
    ref_df = pd.read_csv(cfg['paths']['grna_reference'], comment='#')
    grna_seqs = ref_df['grna_seq'].unique().tolist()

    # Load barcode TSV from Step 2
    output_dir = cfg['paths']['output_dir']
    barcode_tsv = os.path.join(output_dir, 'grna_barcodes', f"{args.library}_barcodes.tsv")
    if not os.path.exists(barcode_tsv):
        print(f"ERROR: {barcode_tsv} not found — run Step 2 first!")
        sys.exit(1)

    print(f"{'='*60}")
    print(f"  Step 3 QC: Diagnosing extraction failures")
    print(f"{'='*60}")
    print(f"  Library:    {args.library}")
    print(f"  Barcode TSV: {barcode_tsv}")
    print(f"  Sample size: {args.n_sample:,}")
    print(f"  Edlib:       {'available' if HAS_EDLIB else 'NOT available'}")
    print()

    # Read Step 2 output
    df = pd.read_csv(barcode_tsv, sep='\t', nrows=args.n_sample)
    print(f"  Loaded {len(df):,} reads from Step 2 output")

    # Replay extraction to identify failed reads
    print(f"  Replaying extraction to identify failures...")

    failed_reads = []
    passed_reads = 0
    total = 0

    t_start = time.time()

    for _, row in df.iterrows():
        total += 1
        r1_seq = row['r1_seq']

        # Quick extraction check (reproduce Step 3 logic)
        tso_end, tso_dist = (-1, -1)
        if HAS_EDLIB:
            target = r1_seq[:50]
            result = edlib.align(tso_seq, target, mode="HW", k=tso_threshold, task="locations")
            if result['editDistance'] >= 0 and result['locations']:
                tso_end = result['locations'][0][1] + 1
                tso_dist = result['editDistance']

        extracted = False
        if tso_end >= 0:
            for grna_len in sorted(grna_lengths, reverse=True):
                grna_end = tso_end + grna_len
                stem_end = grna_end + len(stem1)
                ibar_end = stem_end + ibar_length

                if ibar_end > len(r1_seq):
                    continue

                stem_obs = r1_seq[grna_end:stem_end]
                mm1 = hamming_distance(stem_obs, stem1)
                mm2 = hamming_distance(stem_obs, stem2)

                if mm1 <= anchor_max_mm or mm2 <= anchor_max_mm:
                    extracted = True
                    break

        if extracted:
            passed_reads += 1
        else:
            failed_reads.append(r1_seq)

    elapsed = time.time() - t_start
    fail_rate = len(failed_reads) / total * 100 if total else 0

    print(f"  Replayed {total:,} reads in {elapsed:.1f}s")
    print(f"  Passed: {passed_reads:,} ({passed_reads/total*100:.1f}%)")
    print(f"  Failed: {len(failed_reads):,} ({fail_rate:.1f}%)")
    print()

    # Classify all failed reads
    print(f"  Classifying {len(failed_reads):,} failed reads...")

    classifications = []
    class_counts = Counter()

    for r1_seq in failed_reads:
        cls = classify_failed_read(
            r1_seq, tso_seq, stem1, stem2, grna_lengths, ibar_length,
            grna_seqs, anchor_max_mm, tso_threshold
        )
        classifications.append(cls)
        class_counts[cls['failure_class']] += 1

    # === Report ===
    print()
    print(f"{'='*60}")
    print(f"  FAILURE CLASSIFICATION RESULTS")
    print(f"{'='*60}")
    print(f"  Total sampled:  {total:,}")
    print(f"  Total failed:   {len(failed_reads):,} ({fail_rate:.1f}%)")
    print()

    print(f"  {'Class':<30s} {'Count':>10s} {'% of failed':>12s} {'% of total':>12s}")
    print(f"  {'-'*64}")

    for cls, count in class_counts.most_common():
        pct_fail = count / len(failed_reads) * 100 if failed_reads else 0
        pct_total = count / total * 100 if total else 0
        label = cls.replace('_', ' ')
        print(f"  {label:<30s} {count:>10,} {pct_fail:>11.1f}% {pct_total:>11.1f}%")

    # === Sequence examples for each class ===
    print()
    print(f"  EXAMPLE SEQUENCES (first 3 per class)")
    print(f"  {'-'*60}")

    class_examples = defaultdict(list)
    for i, cls in enumerate(classifications):
        if len(class_examples[cls['failure_class']]) < 3:
            class_examples[cls['failure_class']].append((i, cls, failed_reads[i]))

    for class_name in class_counts.keys():
        print(f"\n  [{class_name}] ({class_counts[class_name]:,} reads)")
        for idx, cls, seq in class_examples[class_name]:
            print(f"    Detail: {cls['failure_detail']}")
            print(f"    Seq:    {seq[:80]}{'...' if len(seq) > 80 else ''}")
            print(f"    Length: {len(seq)}, GC: {cls['gc_frac']:.2f}, N: {cls['n_count']}")
            print()

    # === Read length distribution for failed reads ===
    print(f"  READ LENGTH DISTRIBUTION (failed only)")
    print(f"  {'-'*40}")
    lengths = [len(s) for s in failed_reads]
    if lengths:
        len_hist = Counter()
        for l in lengths:
            bucket = (l // 10) * 10
            len_hist[bucket] += 1
        for bucket in sorted(len_hist.keys()):
            bar = '█' * max(1, int(len_hist[bucket] / max(len_hist.values()) * 40))
            print(f"  {bucket:>3d}-{bucket+9:<3d}bp: {len_hist[bucket]:>8,}  {bar}")

    # === TSO edit distance distribution for degraded TSOs ===
    tso_degraded = [c for c in classifications if c['failure_class'] == 'tso_degraded']
    if tso_degraded:
        print(f"\n  TSO EDIT DISTANCE (for 'tso_degraded' reads)")
        print(f"  {'-'*40}")
        ed_counts = Counter(c.get('tso_edit_dist', -1) for c in tso_degraded)
        for ed in sorted(ed_counts.keys()):
            bar = '█' * max(1, int(ed_counts[ed] / max(ed_counts.values()) * 30))
            print(f"  ed={ed}: {ed_counts[ed]:>8,}  {bar}")

    # === Save detailed results ===
    qc_dir = args.output_dir or os.path.join(output_dir, 'qc')
    os.makedirs(qc_dir, exist_ok=True)

    # Save classification summary
    summary_path = os.path.join(qc_dir, f"{args.library}_step3_qc_summary.yaml")
    summary = {
        'library': args.library,
        'total_sampled': int(total),
        'total_failed': len(failed_reads),
        'fail_rate_pct': round(fail_rate, 2),
        'classifications': {k: int(v) for k, v in class_counts.most_common()},
    }
    with open(summary_path, 'w') as f:
        yaml.dump(summary, f, default_flow_style=False)
    print(f"\n  Summary saved to: {summary_path}")

    # Save detailed per-read classifications
    detail_path = os.path.join(qc_dir, f"{args.library}_step3_qc_details.tsv")
    detail_df = pd.DataFrame(classifications)
    detail_df['r1_seq_prefix'] = [s[:60] for s in failed_reads]
    detail_df.to_csv(detail_path, sep='\t', index=False)
    print(f"  Details saved to: {detail_path}")

    # === Actionable recommendations ===
    print(f"\n{'='*60}")
    print(f"  RECOMMENDATIONS")
    print(f"{'='*60}")

    total_fail = len(failed_reads) or 1

    if class_counts.get('off_target', 0) / total_fail > 0.3:
        print("  ⚠ >30% off-target: significant non-gRNA contamination.")
        print("    → Check library prep QC, consider stricter size selection.")

    if class_counts.get('tso_degraded', 0) / total_fail > 0.3:
        print("  ⚠ >30% TSO degraded: many reads have TSO but >3 edits.")
        print("    → Consider increasing tso_max_edit_dist to 4 or 5 in config.")

    if class_counts.get('stem_mismatch', 0) / total_fail > 0.3:
        print("  ⚠ >30% stem mismatch: TSO found but stem doesn't match.")
        print("    → Check stem sequences in config, consider increasing fixed_seq_max_mismatches.")

    if class_counts.get('too_short', 0) / total_fail > 0.2:
        print("  ⚠ >20% too short: many reads truncated before structure.")
        print("    → Check sequencing length settings, quality trimming, adapter contamination.")

    if class_counts.get('poly_G_in_cassette', 0) / total_fail > 0.1:
        print("  ⚠ >10% poly-G within cassette: NovaSeq artifact overlaps gRNA region.")
        print("    → Apply poly-G trimming (fastp/cutadapt) before pipeline.")

    if class_counts.get('tso_absent_grna_present', 0) / total_fail > 0.1:
        print("  ⚠ >10% gRNA present without TSO: possible chimeric reads.")
        print("    → Check library structure, may indicate PCR chimera formation.")

    if fail_rate < 15:
        print("  ✓ Extraction failure rate is acceptable (<15%).")
    elif fail_rate < 25:
        print("  ⚠ Extraction failure rate is moderate (15-25%). Review dominant failure mode.")
    else:
        print("  ✗ Extraction failure rate is high (>25%). Action recommended.")

    # === Threshold sensitivity analysis ===
    print(f"\n{'='*60}")
    print(f"  THRESHOLD SENSITIVITY ANALYSIS")
    print(f"{'='*60}")
    print(f"  Sweeping tso_max_edit_dist × fixed_seq_max_mismatches\n"
          f"  over {total:,} sampled reads to find optimal thresholds.\n")

    tso_dists_to_test = list(range(1, 7))      # 1..6
    stem_mm_to_test = list(range(1, 6))         # 1..5

    # Precompute per-read: best TSO edit distance and best stem Hamming distance
    # so we can quickly tabulate extraction at each threshold combo.
    print("  Pre-computing TSO/stem distances for all reads...")
    t_sweep = time.time()

    read_tso_ed = []    # best TSO edit distance per read (-1 if not found)
    read_stem_hd = []   # best stem Hamming distance per read (-1 if TSO not found or read too short)

    for _, row in df.iterrows():
        r1_seq = row['r1_seq']

        # Find TSO at the most relaxed threshold we test
        tso_end_best = -1
        tso_ed_best = -1
        if HAS_EDLIB:
            target = r1_seq[:50]
            result = edlib.align(tso_seq, target, mode="HW", k=max(tso_dists_to_test), task="locations")
            if result['editDistance'] >= 0 and result['locations']:
                tso_end_best = result['locations'][0][1] + 1
                tso_ed_best = result['editDistance']

        read_tso_ed.append(tso_ed_best)

        # Find best stem distance (across all gRNA lengths)
        best_stem_d = -1
        if tso_end_best >= 0:
            for grna_len in grna_lengths:
                grna_end = tso_end_best + grna_len
                stem_end_pos = grna_end + len(stem1)
                ibar_end_pos = stem_end_pos + ibar_length
                if ibar_end_pos <= len(r1_seq):
                    stem_obs = r1_seq[grna_end:stem_end_pos]
                    d1 = hamming_distance(stem_obs, stem1)
                    d2 = hamming_distance(stem_obs, stem2)
                    d_min = min(d1, d2)
                    if best_stem_d < 0 or d_min < best_stem_d:
                        best_stem_d = d_min
        read_stem_hd.append(best_stem_d)

    print(f"  Done in {time.time() - t_sweep:.1f}s\n")

    # Build extraction rate table
    # Header
    header = f"  {'tso\\stem':<10s}"
    for sm in stem_mm_to_test:
        header += f" {'mm≤'+str(sm):>8s}"
    print(header)
    print(f"  {'-' * (10 + 9 * len(stem_mm_to_test))}")

    # Current config values for highlighting
    curr_tso = tso_threshold
    curr_stem = anchor_max_mm

    best_rate = 0
    best_params = (curr_tso, curr_stem)

    for td in tso_dists_to_test:
        row_str = f"  {'ed≤'+str(td):<10s}"
        for sm in stem_mm_to_test:
            # Count reads that would be extracted at (td, sm)
            extracted = sum(
                1 for te, sh in zip(read_tso_ed, read_stem_hd)
                if 0 <= te <= td and 0 <= sh <= sm
            )
            rate = extracted / total * 100 if total else 0
            marker = " *" if td == curr_tso and sm == curr_stem else "  "
            row_str += f" {rate:>6.1f}%{marker[1]}"
            if rate > best_rate:
                best_rate = rate
                best_params = (td, sm)
        print(row_str)

    curr_extracted = sum(
        1 for te, sh in zip(read_tso_ed, read_stem_hd)
        if 0 <= te <= curr_tso and 0 <= sh <= curr_stem
    )
    curr_rate = curr_extracted / total * 100 if total else 0

    print(f"\n  * = current config (tso_ed≤{curr_tso}, stem_mm≤{curr_stem})")
    print(f"    Current extraction rate: {curr_rate:.1f}%")
    print(f"    Best achievable rate:    {best_rate:.1f}% "
          f"(tso_ed≤{best_params[0]}, stem_mm≤{best_params[1]})")
    gain = best_rate - curr_rate
    if gain > 1.0:
        print(f"    Potential gain:          +{gain:.1f}% "
              f"(~{int(gain / 100 * total):,} more reads)")
        print(f"\n  → Consider: tso_max_edit_dist: {best_params[0]}, "
              f"fixed_seq_max_mismatches: {best_params[1]}")
    else:
        print(f"    ✓ Current thresholds are near-optimal (≤1% potential gain).")

    # Also show marginal gains: what does each parameter increase buy?
    print(f"\n  MARGINAL GAINS (reads recovered by relaxing one parameter):")
    print(f"  {'-'*55}")
    for td in tso_dists_to_test:
        if td <= curr_tso:
            continue
        gained = sum(
            1 for te, sh in zip(read_tso_ed, read_stem_hd)
            if curr_tso < te <= td and 0 <= sh <= curr_stem
        )
        pct = gained / total * 100 if total else 0
        print(f"    tso_max_edit_dist {curr_tso}→{td}: "
              f"+{gained:,} reads (+{pct:.1f}%)")

    for sm in stem_mm_to_test:
        if sm <= curr_stem:
            continue
        gained = sum(
            1 for te, sh in zip(read_tso_ed, read_stem_hd)
            if 0 <= te <= curr_tso and curr_stem < sh <= sm
        )
        pct = gained / total * 100 if total else 0
        print(f"    fixed_seq_max_mismatches {curr_stem}→{sm}: "
              f"+{gained:,} reads (+{pct:.1f}%)")

    # Save sweep results
    sweep_rows = []
    for td in tso_dists_to_test:
        for sm in stem_mm_to_test:
            ext = sum(
                1 for te, sh in zip(read_tso_ed, read_stem_hd)
                if 0 <= te <= td and 0 <= sh <= sm
            )
            sweep_rows.append({
                'tso_max_edit_dist': td,
                'stem_max_mismatches': sm,
                'extracted': ext,
                'total': total,
                'extraction_rate_pct': round(ext / total * 100, 2) if total else 0,
                'is_current': td == curr_tso and sm == curr_stem,
            })
    sweep_path = os.path.join(qc_dir, f"{args.library}_step3_threshold_sweep.tsv")
    pd.DataFrame(sweep_rows).to_csv(sweep_path, sep='\t', index=False)
    print(f"\n  Sweep results saved to: {sweep_path}")

    print(f"\n{'='*60}")
    print(f"  QC complete!")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()

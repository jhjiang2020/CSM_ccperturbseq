#!/usr/bin/env python3
"""Compare a shared-policy pilot benchmark with archived v2.1 evidence."""
import argparse
import gzip
import csv
import json
from pathlib import Path
import pandas as pd
from benchmark_grna_recovery import load


def same_table(left, right, ignore=()):
    a, b = (pd.read_csv(path, keep_default_na=False) for path in (left, right))
    columns = sorted(set(a) - set(ignore))
    if set(columns) != set(b) - set(ignore):
        return False
    try:
        pd.testing.assert_frame_equal(a[columns].sort_values(columns).reset_index(drop=True),
                                      b[columns].sort_values(columns).reset_index(drop=True), check_dtype=False)
    except AssertionError:
        return False
    return True


def run(base, previous):
    root = Path(__file__).resolve().parents[1]
    step4 = load(root/'scripts/04_assign_grna_to_cells.py', 'shared_count_validation')
    import yaml
    cfg = yaml.safe_load((base/'config.yaml').read_text())
    totals, comparisons, transitions = [], {}, []
    for index, lib in enumerate(cfg['libraries']):
        name = lib['sublibrary_name']
        current, old = base/'technical'/name, previous/'technical'/name
        assert (current/'_SUCCESS').exists(), name
        # Every input consumed by the biological auto-arm comparisons must match
        # before an archived effect estimate can be claimed for the revised run.
        comparable = ['baseline_assignments.csv.gz', 'baseline_umi_counts.csv.gz',
                      'auto_assignments.csv.gz', 'auto_umi_counts.csv.gz',
                      'additional_normal_exclusions.csv.gz']
        comparison = {file: same_table(old/file, current/file, ['recovery_policy_fingerprint'])
                      for file in comparable}
        comparisons[name] = comparison
        read = json.loads((current/'read_summary.json').read_text())
        assert read['baseline_parity']['identical'], name
        summary = pd.read_csv(current/'technical_summary.csv')
        before = pd.read_csv(current/'baseline_assignments.csv.gz', keep_default_na=False).drop_duplicates('bc_wells').set_index('bc_wells').combo_type
        after = pd.read_csv(current/'auto_assignments.csv.gz', keep_default_na=False).drop_duplicates('bc_wells').set_index('bc_wells').combo_type
        cells = sorted(set(before.index)|set(after.index))
        changes = pd.DataFrame(dict(before=before.reindex(cells).fillna('unassigned'),
                                    after=after.reindex(cells).fillna('unassigned')))
        for (a, b), n in changes.value_counts().items():
            transitions.append(dict(library=name, before=a, after=b, cells=n))
        for row in summary.to_dict('records'):
            row['singlet_cells'] = int((before if row['arm']=='baseline' else after).eq('singlet').sum())
            row['input_reads'] = read['baseline']['total_reads']
            row['all_matched_reads'] = sum(read['baseline'].get(k, 0) for k in ('exact_match','fuzzy_match','grna_only_match'))
            if row['arm']=='auto':
                row['all_matched_reads'] += read['recovery'].get('auto|recovered', 0)
            totals.append(row)
        # Recount all retained-cell molecules using the compact production
        # provenance column and compare with the scientific audit schema.
        compact = current/'compact_delta.tsv.gz'
        with gzip.open(current/'auto_delta.tsv.gz', 'rt') as source, gzip.open(compact, 'wt') as destination:
            reader = csv.DictReader(source, delimiter='\t')
            fields = [c for c in reader.fieldnames if not c.startswith('recovery_')]+['grna_seq_corrected','structure_correction']
            writer = csv.DictWriter(destination, fieldnames=fields, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            for row in reader:
                row['grna_seq_corrected'] = row['recovery_designed_spacer']
                row['structure_correction'] = json.dumps(dict(length=row['recovery_length_state'], policy=row['recovery_policy_fingerprint']))
                row['grna_dist'] = len(row['grna_seq_corrected'])-len(row['grna_seq'])
                row['match_type'] = 'relaxed_structure'
                writer.writerow(row)
        groups = step4._determine_group_cols(str(current/'baseline_matches.tsv.gz'), True)
        actual, _ = step4.stream_library_umi_dedup_sources(
            [('baseline',current/'baseline_matches.tsv.gz'),('rescue',compact)], name,
            lib['sample_name'], index, groups, cfg['guide_assignment'].get('max_umi_edits',1),
            cfg['guide_assignment'].get('min_reads_per_umi',1))
        expected = pd.read_csv(current/'auto_umi_counts.csv.gz', keep_default_na=False)
        detail = pd.DataFrame([json.loads(v) for v in actual.pop('structure_correction')])
        actual = pd.concat([actual, detail],axis=1)
        columns = groups+[c for c in actual if c not in groups and ('count' in c or c.startswith('recovery_umi_'))]
        pd.testing.assert_frame_equal(actual[columns].fillna('').sort_values(groups).reset_index(drop=True),
                                      expected[columns].fillna('').sort_values(groups).reset_index(drop=True), check_dtype=False)
        comparison['compact_full_count_parity'] = True
        print(name, comparison, flush=True)
    table = pd.DataFrame(totals)
    table.to_csv(base/'technical_summary.csv', index=False)
    metrics = ['input_reads','all_matched_reads','reads','umis','assigned_cells','pairs','complete_pairs','new_cells','additional_pairs','lost_pairs','recovered_partners','singlet_cells']
    combined = table.groupby('arm')[metrics].sum()
    combined.to_csv(base/'technical_combined.csv')
    pd.DataFrame(transitions).to_csv(base/'cell_transitions.csv',index=False)
    result = dict(libraries=comparisons, biological_auto_inputs_identical=all(all(v.values()) for v in comparisons.values()),
                  note='Biological reuse additionally requires identical frozen normal exclusions and strict training/control inputs; no expression counts changed.')
    (base/'comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    print(combined.to_string(), flush=True)
    (base/'_VALIDATION_SUCCESS').write_text('complete\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--previous', type=Path, required=True)
    args = parser.parse_args()
    run(args.output, args.previous)

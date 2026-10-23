"""Observable design compatibility QC without fully resolving iBAR assumptions."""

import csv
import gzip
import json
from collections import Counter
from pathlib import Path

import pandas as pd
import yaml

try:
    from scripts.guide_design import DesignMatcher, compile_design, extract_unit, linked_pair_status, reverse_complement
except ModuleNotFoundError:
    from guide_design import DesignMatcher, compile_design, extract_unit, linked_pair_status, reverse_complement


def cell_design_qc(assignments, molecules, out_dir, valid_bcs=None):
    directory = Path(out_dir)
    directory.mkdir(parents=True, exist_ok=True)
    if valid_bcs is not None:
        assignments = assignments[assignments['bc_wells'].isin(valid_bcs)]
        molecules = molecules[molecules['bc_wells'].isin(valid_bcs)]
    summary = {
        'status': 'observed_compatibility_only',
        'cells_with_guide_evidence': int(molecules['bc_wells'].nunique()),
        'assignment_status_counts': {str(key): int(value) for key, value in assignments['detection_status'].value_counts().items()},
        'umi_by_ibar_status': {str(key): int(value) for key, value in molecules.groupby('ibar_status')['umi_count'].sum().items()},
        'legacy_recombination_rate': None,
        'reason': 'Cell co-detection is not physical linkage; the legacy likelihood assumes fully resolving iBAR categories',
    }
    with (directory / 'q2_design_capabilities.yaml').open('w') as handle:
        yaml.safe_dump(summary, handle, sort_keys=False)
    return summary


def _extract_arm(sequence, left_anchor, contract, find_anchor, max_mismatches, read_position, right_anchor=None):
    if not left_anchor:
        return None
    _, end, _ = find_anchor(sequence, left_anchor, max_mismatches)
    if end < 0:
        return None
    observation = extract_unit(sequence, end, contract, max_mismatches, right_anchor=right_anchor)
    if observation and observation['position'] not in (0, read_position):
        observation['structure_status'] = 'conflict'
    return observation


def run_guideseq(fastq_pairs, reference, settings, out_dir, find_anchor, read_pairs,
                max_lev_dist=3, max_lev_dist_ibar=3, max_reads=0, matcher_factory=None):
    """Stream linked R1/R2 observations; report only identifiable metrics."""
    config = {
        'guide_design': settings.get('guide_design', {'enabled': True}),
        'grna_structure': {
            'fixed_seq_2': settings['stem_1'], 'fixed_seq_3': settings['stem_2'],
            'fixed_seq_4': settings['tracr_1'][:settings.get('max_anchor_length', 20) or None],
            'fixed_seq_5': settings['tracr_2'][:settings.get('max_anchor_length', 20) or None],
            'grna_lengths': settings.get('grna_lengths', [20, 21]),
            'ibar_lengths_by_position': settings.get('ibar_lengths_by_position', 'auto'),
        },
        'grna_matching': {'max_lev_dist_ibar': max_lev_dist_ibar},
    }
    contract = compile_design(reference, config)
    frozen_path = settings.get('guide_design', {}).get('resolved_contract')
    if frozen_path:
        with Path(frozen_path).expanduser().open() as handle:
            frozen = yaml.safe_load(handle)
        if frozen.get('input_reference_fingerprint', frozen.get('reference_fingerprint')) != contract['reference_fingerprint']:
            raise ValueError('Q1 frozen contract does not match the supplied reference')
        effective = compile_design(frozen['rows'], config)
        if effective['reference_fingerprint'] != frozen.get('reference_fingerprint'):
            raise ValueError('Q1 frozen contract content does not match its fingerprint')
        contract = effective
    matcher = DesignMatcher(contract)
    spacers = list(matcher.by_spacer)
    spacer_matcher = matcher_factory(spacers, spacers, max_lev_dist) if matcher_factory else None
    directory = Path(out_dir)
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / 'q1_design_capabilities.yaml').open('w') as handle:
        yaml.safe_dump(contract, handle, sort_keys=False)
    per_library = []
    extraction_files = []
    all_pairs = Counter()
    all_counts = Counter()
    max_mismatches = settings.get('stem_max_mm', 2)
    for name, first_path, second_path in fastq_pairs:
        library_dir = directory / name
        library_dir.mkdir(parents=True, exist_ok=True)
        extraction_files.append((name, {'linked_tsv_path': str(library_dir / 'design_linked_observations.tsv.gz')}))
        counts = Counter()
        pair_counts = Counter()
        with gzip.open(library_dir / 'design_linked_observations.tsv.gz', 'wt') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['read_id', 'status', 'pair_candidates', 'pos1', 'pos2'])
            for index, (read_id, first_sequence, second_sequence) in enumerate(read_pairs(first_path, second_path)):
                if max_reads and index >= max_reads:
                    break
                counts['total'] += 1
                first = _extract_arm(first_sequence, settings.get('R1_primer', ''), contract,
                                     find_anchor, max_mismatches, 1)
                second_forward = reverse_complement(second_sequence.upper())
                second = None
                anchors = [settings['R2_left_anchor']] if settings.get('R2_left_anchor') else []
                for trna in settings.get('tRNA_seqs', {}).values():
                    anchors.extend([settings.get('tRNA_leader', '') + trna, trna[-settings.get('max_anchor_length', 20):]])
                second_candidates = []
                for anchor in anchors:
                    candidate = _extract_arm(second_forward, anchor, contract, find_anchor, max_mismatches, 2,
                                             reverse_complement(settings.get('R2_primer', '').upper()) or None)
                    if candidate:
                        second_candidates.append(candidate)
                signatures = {(item['grna_seq'], item['ibar_seq'], item['position']) for item in second_candidates}
                if len(signatures) == 1:
                    second = min(second_candidates, key=lambda item: item['anchor_score'])
                resolved = []
                for observation in (first, second):
                    result = matcher.resolve(observation) if observation else None
                    if result is None and observation and spacer_matcher:
                        nearest, _, unique = spacer_matcher.match(observation['grna_seq'])
                        if unique:
                            result = matcher.resolve(observation, nearest)
                    resolved.append(result)
                status = linked_pair_status(*resolved)
                for position, result in enumerate(resolved, 1):
                    if result and result.get('ibar_decoupling_eligible'):
                        counts[f'pos{position}_decoupling_eligible'] += 1
                        counts[f'pos{position}_decoupled'] += bool(result['ibar_decoupled'])
                counts[status['status']] += 1
                if status['status'] in ('designed_unique', 'designed_ambiguous', 'recombination_compatible'):
                    counts['eligible'] += 1
                if status['status'] == 'designed_unique':
                    pair_counts[status['pair_candidates'][0]] += 1
                writer.writerow([read_id, status['status'], json.dumps(status['pair_candidates']),
                                 *(json.dumps(result, sort_keys=True) for result in resolved)])
        rates = _rates(counts)
        pd.DataFrame([{'library': name, **rates}]).to_csv(library_dir / 'q1_rates_design.csv', index=False)
        per_library.append((name, rates))
        all_pairs.update(pair_counts)
        all_counts.update(counts)
    pair_metadata = {row['pair_id']: row for row in contract['rows']}
    pair_guides = {}
    for row in contract['rows']:
        pair_guides.setdefault(row['pair_id'], {})[row['position']] = row['grna_name']
    frequencies = pd.DataFrame([
        {'pair_id': pair_id, 'target_combo': pair_metadata[pair_id]['target_combo'],
         'pos1_grna': pair_guides[pair_id].get(1, ''), 'pos2_grna': pair_guides[pair_id].get(2, ''),
         'pos1_guide': pair_guides[pair_id].get(1, ''), 'pos2_guide': pair_guides[pair_id].get(2, ''),
         'count': count, 'frequency': count / sum(all_pairs.values()),
         'eligible_unique_constructs': sum(all_pairs.values())}
        for pair_id, count in sorted(all_pairs.items())
    ], columns=['pair_id', 'target_combo', 'pos1_grna', 'pos2_grna', 'pos1_guide', 'pos2_guide', 'count', 'frequency', 'eligible_unique_constructs'])
    frequencies.to_csv(directory / 'q1_combo_frequencies.csv', index=False)
    guide_counts = Counter()
    seen = set()
    for row in contract['rows']:
        signature = (row['pair_id'], row['position'], row['grna_name'])
        if signature in seen:
            continue
        seen.add(signature)
        guide_counts[(row['grna_name'], row['target_gene'], row['position'])] += all_pairs[row['pair_id']]
    total_guides = sum(guide_counts.values())
    guide_frequencies = pd.DataFrame([
        {'guide_name': name, 'target_gene': target, 'position': position, 'count': count,
         'frequency': count / total_guides if total_guides else 0}
        for (name, target, position), count in sorted(guide_counts.items()) if count
    ], columns=['guide_name', 'target_gene', 'position', 'count', 'frequency'])
    guide_frequencies.to_csv(directory / 'q1_guide_frequencies.csv', index=False)
    pd.DataFrame([{'library': name, **rates} for name, rates in per_library]).to_csv(directory / 'q1_guideseq_summary.csv', index=False)
    return {'per_lib_rates': per_library, 'merged_rates': _rates(all_counts),
            'per_lib_extraction': extraction_files,
            'guide_frequencies': guide_frequencies, 'combo_frequencies': frequencies,
            'capabilities': contract['capabilities']}


def _rates(counts):
    eligible = counts['eligible']
    compatible = counts['designed_unique'] + counts['designed_ambiguous']
    result = {'N_total': counts['total'], 'N_eligible': eligible,
            'intact_rate': compatible / eligible if eligible else None,
            'design_compatible_rate': compatible / eligible if eligible else None,
            'uniquely_resolved_pair_rate': counts['designed_unique'] / eligible if eligible else None,
            'recombination_compatible_rate': counts['recombination_compatible'] / eligible if eligible else None,
            'interpretation': 'Observable linked design incompatibility, including sequence errors; events producing another valid signature are invisible',
            'legacy_swap_rates': None, 'legacy_swap_reason': 'The design-aware path reports observed compatibility; the historical event decomposition remains available in legacy symmetric Q1 mode',
            **{f'n_{status}': counts[status] for status in
               ('designed_unique', 'designed_ambiguous', 'recombination_compatible', 'unobserved', 'unresolved', 'structure_conflict')}}
    for position in (1, 2):
        denominator = counts[f'pos{position}_decoupling_eligible']
        result[f'decoupling_eligible_pos{position}'] = denominator
        result[f'decoupled_pos{position}'] = counts[f'pos{position}_decoupled'] / denominator if denominator else None
        result[f'decoupling_reason_pos{position}'] = 'Known iBARs in spacer groups with distinguishable forbidden combinations' if denominator else 'No informative spacer/iBAR observations; absent or indistinguishable iBARs cannot measure decoupling'
    return result

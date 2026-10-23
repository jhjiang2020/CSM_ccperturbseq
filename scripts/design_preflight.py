"""Bounded, independently anchored pilot for optional guide designs."""

import copy
import math
from collections import Counter, defaultdict
from contextlib import closing
from itertools import islice
from pathlib import Path

try:
    from scripts.guide_design import compile_design, extract_unit, reverse_complement
    from scripts.grna_preflight import _barcode_sequences, _fastq_sequences
except ModuleNotFoundError:
    from guide_design import compile_design, extract_unit, reverse_complement
    from grna_preflight import _barcode_sequences, _fastq_sequences


def check_design(barcode_path, contract, config, find_tso, sample_reads=10_000, fastq_paths=None):
    """Discover boundaries without assuming the configured iBAR lengths."""
    options = config.get('guide_design', {}) or {}
    maximum = options.get('discovery_max_ibar_length', 32)
    if isinstance(maximum, bool) or not isinstance(maximum, int) or not 0 <= maximum <= 100:
        raise ValueError('discovery_max_ibar_length must be an integer from 0 to 100')
    orientation = options.get('ibar_orientation', 'as_reference')
    if orientation not in ('auto', 'as_reference'):
        raise ValueError('ibar_orientation must be auto or as_reference')
    guide_lengths = options.get('discovery_grna_lengths', [18, 19, 20, 21, 22, 23])
    if not isinstance(guide_lengths, list) or not guide_lengths or any(isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= 100 for value in guide_lengths):
        raise ValueError('discovery_grna_lengths must be a nonempty list of lengths from 1 to 100')
    discovery_contract = dict(contract)
    discovery_contract['grna_lengths'] = sorted(set(guide_lengths) | set(contract['grna_lengths']))
    phases = {'discovery': defaultdict(Counter), 'validation': defaultdict(Counter)}
    length_counts = defaultdict(Counter)
    spacer_lengths = Counter()
    guide_reference = defaultdict(list)
    for row in contract['rows']:
        guide_reference[row['grna_seq']].append(row)
    structure = config['grna_structure']
    total, extracted, conflicts, guide_conflicts, exact_guides = 0, 0, 0, 0, 0
    paths = list(fastq_paths) if fastq_paths is not None else [None]
    if not paths:
        raise ValueError('No pilot input files')
    per_file = math.ceil(sample_reads / len(paths))
    file_evidence = {}
    for path in paths:
        file_counts = defaultdict(Counter)
        file_name = Path(path).name if path is not None else 'barcode_tsv'
        sequences = _fastq_sequences([path]) if path is not None else _barcode_sequences(barcode_path)
        with closing(sequences):
            for local_index, sequence in enumerate(islice(sequences, min(per_file, sample_reads - total))):
                sequence = sequence.upper()
                total += 1
                start, _ = find_tso(sequence, structure['fixed_seq_1'], structure.get('tso_max_edit_dist', 3))
                observation = extract_unit(sequence, start, discovery_contract, structure.get('fixed_seq_max_mismatches', 2), range(maximum + 1)) if start >= 0 else None
                if observation is None:
                    continue
                extracted += 1
                spacer_lengths[observation['grna_length']] += 1
                guides = guide_reference.get(observation['grna_seq'], [])
                exact_guides += bool(guides)
                spacer_positions = {row['position'] for row in guides}
                position = observation['position_from_scaffold'] or (next(iter(spacer_positions)) if len(spacer_positions) == 1 else 0)
                conflicts += observation['structure_status'] != 'concordant'
                if not position:
                    continue
                length_counts[position][observation['ibar_length']] += 1
                observed_positions = {value for value in (observation['position_from_length'], observation['position_from_scaffold']) if value}
                if guides and observed_positions - spacer_positions:
                    guide_conflicts += 1
                ibars = {row['ibar_seq'] for row in guides if row['position'] == position and row['ibar_seq']}
                observed = observation['ibar_seq']
                if not ibars or not observed or observation['structure_status'] != 'concordant' or (observation['position_from_length'] and observation['position_from_length'] != position):
                    continue
                phase = phases['discovery' if local_index % 2 == 0 else 'validation'][position]
                phase['informative'] += 1
                stored = observed in ibars
                reversed_match = observed in {reverse_complement(ibar) for ibar in ibars}
                phase['stored_only'] += stored and not reversed_match
                phase['rc_only'] += reversed_match and not stored
                phase['both'] += stored and reversed_match
                file_counts[position]['stored_only'] += stored and not reversed_match
                file_counts[position]['rc_only'] += reversed_match and not stored
        file_evidence[file_name] = {position: dict(counts) for position, counts in file_counts.items()}
    messages = []

    def warn(code, message):
        messages.append({'code': code, 'message': message})

    resolved_orientation = {}
    changed = set()
    for position, profile in contract['profiles'].items():
        observed_lengths = length_counts[position]
        unexpected = sum(count for length, count in observed_lengths.items() if length not in profile['ibar_lengths'])
        if sum(observed_lengths.values()) >= 100 and unexpected / sum(observed_lengths.values()) >= 0.10:
            warn('observed_ibar_length_mismatch', f'Position {position}: observed lengths {dict(observed_lengths)} disagree with design lengths {profile["ibar_lengths"]}')
        votes = []
        for phase in phases.values():
            counts = phase[position]
            stored, reversed_count = counts['stored_only'], counts['rc_only']
            winner = 'reverse_complement' if reversed_count > stored else 'as_reference'
            high, low = max(stored, reversed_count), min(stored, reversed_count)
            confident = counts['informative'] >= 50 and high >= 20 and high >= 5 * max(low, 1) and (high - low) / counts['informative'] >= 0.10
            votes.append(winner if confident else 'inconclusive')
        informative = sum(phase[position]['informative'] for phase in phases.values())
        stored_support = sum(phase[position]['stored_only'] + phase[position]['both'] for phase in phases.values())
        reverse_support = sum(phase[position]['rc_only'] + phase[position]['both'] for phase in phases.values())
        if informative >= 100 and max(stored_support, reverse_support) / informative < 0.10:
            warn('low_ibar_support', f'Position {position}: fewer than 10% of exact-guide observations match either iBAR orientation; check reference/construct identity')
        conclusion = votes[0] if votes[0] == votes[1] else 'inconclusive'
        file_votes = set()
        for counts_by_position in file_evidence.values():
            counts = counts_by_position.get(position, {})
            stored, reversed_count = counts.get('stored_only', 0), counts.get('rc_only', 0)
            if max(stored, reversed_count) >= 20 and max(stored, reversed_count) >= 5 * max(min(stored, reversed_count), 1):
                file_votes.add('reverse_complement' if reversed_count > stored else 'as_reference')
        if len(file_votes) > 1:
            conclusion = 'inconclusive'
            warn('mixed_file_orientation', f'Position {position}: input files support conflicting orientations; no automatic change applied')
        if profile['ibar_lengths'] == [0]:
            conclusion = 'not_applicable'
        resolved_orientation[position] = conclusion
        if conclusion == 'reverse_complement':
            warn('possible_reverse_complement', f'Position {position}: discovery and held-out validation support reverse-complementing the reference iBARs')
            if orientation == 'auto':
                changed.add(position)
        if conclusion == 'inconclusive':
            warn('orientation_inconclusive', f'Position {position}: insufficient or inconsistent strand evidence; reference orientation retained')
    if extracted < 100:
        warn('low_extraction_support', 'Fewer than 100 independently anchored guide units; construct validation is inconclusive')
    if extracted >= 100 and exact_guides / extracted < 0.10:
        warn('low_guide_support', 'Fewer than 10% of independently extracted spacers exactly match the reference; check library/reference identity and sequencing quality')
    unexpected_spacers = sum(count for length, count in spacer_lengths.items() if length not in contract['grna_lengths'])
    if extracted >= 100 and unexpected_spacers / extracted >= 0.10:
        warn('observed_spacer_length_mismatch', f'Observed spacer lengths {dict(spacer_lengths)} disagree with reference lengths {contract["grna_lengths"]}')
    if guide_conflicts:
        warn('guide_position_conflict', f'{guide_conflicts} exact spacer observations disagree with scaffold position; check guide/construct annotations')
    if conflicts:
        warn('scaffold_length_conflict', f'{conflicts} observations have unexpected length/scaffold combinations; these remain QC evidence')
    effective = copy.deepcopy(contract)
    effective['input_fingerprint'] = contract['fingerprint']
    if changed:
        rows = copy.deepcopy(contract['rows'])
        for row in rows:
            if row['position'] in changed:
                row['ibar_seq'] = reverse_complement(row['ibar_seq'])
        effective = compile_design(rows, config)
        effective['input_fingerprint'] = contract['fingerprint']
    effective['resolved_orientation'] = resolved_orientation
    effective['input_reference_fingerprint'] = contract['reference_fingerprint']
    return {'status': 'warning' if messages else 'ok', 'sample_limit': sample_reads,
            'sampled_reads': total, 'extracted_reads': extracted,
            'exact_guide_reads': exact_guides, 'observed_spacer_length_counts': dict(spacer_lengths),
            'observed_ibar_length_counts': {position: dict(counts) for position, counts in length_counts.items()},
            'orientation': resolved_orientation,
            'file_orientation_evidence': file_evidence,
            'orientation_evidence': {name: {position: dict(counts) for position, counts in values.items()} for name, values in phases.items()},
            'warnings': messages, 'fingerprint': effective['fingerprint'],
            'input_type': 'raw_fastq' if fastq_paths is not None else 'barcode_tsv'}, effective

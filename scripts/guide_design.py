"""Versioned guide-design capabilities shared by processing and QC."""

import hashlib
import json
import warnings
from collections import defaultdict
from itertools import combinations


SCHEMA_VERSION = 2
DESIGN_COLUMNS = [
    'guide_unit_id', 'position', 'ibar_length', 'ibar_status',
    'pair_candidates', 'pair_candidate_count', 'resolution_status',
    'position_from_length', 'position_from_scaffold', 'scaffold_variant',
    'tracr_variant', 'structure_status', 'pair_resolution_basis', 'guide_design_fingerprint',
    'legacy_model_eligible',
]


def enabled(config):
    options = config.get('guide_design', {}) or {}
    if not isinstance(options, dict):
        raise ValueError('guide_design must be a mapping')
    value = options.get('enabled', False)
    if not isinstance(value, bool):
        raise ValueError('guide_design.enabled must be true or false')
    return value or 'ibar_lengths_by_position' in config.get('grna_structure', {})


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def edit_distance(first, second):
    previous = list(range(len(second) + 1))
    for row_index, first_base in enumerate(first, 1):
        current = [row_index]
        for column_index, second_base in enumerate(second, 1):
            current.append(min(current[-1] + 1, previous[column_index] + 1,
                               previous[column_index - 1] + (first_base != second_base)))
        previous = current
    return previous[-1]


def unit_id(position, spacer, ibar):
    payload = json.dumps([int(position), spacer, ibar], separators=(',', ':'))
    return 'unit_' + hashlib.sha256(payload.encode()).hexdigest()[:20]


def _dna(value, label, allow_empty=False):
    sequence = str(value).strip().upper()
    if (not sequence and not allow_empty) or set(sequence) - set('ACGT'):
        raise ValueError(f'{label} must contain only A/C/G/T' + (' or be empty' if allow_empty else ''))
    return sequence


def _position(value):
    text = str(value).removeprefix('pos')
    if text not in ('1', '2'):
        raise ValueError(f'Invalid guide position: {value!r}')
    return int(text)


def compile_design(reference, config):
    """Resolve lengths, collision sets and conservative correction policies."""
    structure = config.get('grna_structure', {}) or {}
    options = config.get('guide_design', {}) or {}
    if not isinstance(options, dict):
        raise ValueError('guide_design must be a mapping')
    strict = options.get('strict', False)
    if not isinstance(strict, bool):
        raise ValueError('guide_design.strict must be true or false')
    messages = []

    def warn(code, message):
        if strict:
            raise ValueError(f'{code}: {message}')
        messages.append({'code': code, 'message': message})

    records = reference.fillna('').to_dict('records') if hasattr(reference, 'to_dict') else reference
    rows = []
    for raw in records:
        row = dict(raw)
        row['position'] = _position(row.get('position', 1))
        row['grna_seq'] = _dna(row['grna_seq'], 'grna_seq')
        row['ibar_seq'] = _dna(row.get('ibar_seq', ''), 'ibar_seq', allow_empty=True)
        for field in ('grna_name', 'target_gene', 'pair_id', 'target_combo'):
            row[field] = str(row.get(field, ''))
        row['guide_unit_id'] = unit_id(row['position'], row['grna_seq'], row['ibar_seq'])
        rows.append(row)
    if not rows:
        raise ValueError('Guide reference is empty')
    positions = sorted({row['position'] for row in rows})
    if config.get('guide_assignment', {}).get('dual_guide_counting_mode', False):
        if positions != [1, 2] or any(not row['pair_id'] for row in rows):
            raise ValueError('dual_guide_counting_mode requires both positions and nonempty pair_id annotations')
    if structure.get('single_guide_mode', False) and positions != [1]:
        raise ValueError('single_guide_mode conflicts with dual-position reference')
    lengths = {position: sorted({len(row['ibar_seq']) for row in rows if row['position'] == position})
               for position in positions}
    configured = structure.get('ibar_lengths_by_position', 'auto')
    if configured != 'auto':
        if not isinstance(configured, dict):
            raise ValueError('ibar_lengths_by_position must be auto or a position mapping')
        for key, values in configured.items():
            position = _position(key)
            if values == 'auto':
                continue
            values = values if isinstance(values, list) else [values]
            if any(isinstance(value, bool) or not isinstance(value, int) or value < 0 for value in values):
                raise ValueError('iBAR lengths must be nonnegative integers')
            if sorted(set(values)) != lengths.get(position, []):
                warn('ibar_length_mismatch', f'Position {position}: requested {values}, reference {lengths.get(position, [])}; using reference lengths')
    elif 'ibar_length' in structure:
        if any(values != [structure['ibar_length']] for values in lengths.values()):
            warn('legacy_ibar_length', 'Scalar ibar_length disagrees with the design; using per-position reference lengths')
    guide_lengths = sorted({len(row['grna_seq']) for row in rows})
    if set(guide_lengths) - set(structure.get('grna_lengths', guide_lengths)):
        warn('spacer_length_mismatch', f'Using reference spacer lengths {guide_lengths}')
    scaffold_mapping = options.get('scaffold_by_position', {'1': 'A', '2': 'B'})
    if not isinstance(scaffold_mapping, dict):
        raise ValueError('scaffold_by_position must be a mapping')
    scaffold_mapping = {_position(key): value for key, value in scaffold_mapping.items()}
    profiles = {}
    for position in positions:
        variant = scaffold_mapping.get(position, 'A' if position == 1 else 'B')
        if variant not in ('A', 'B'):
            raise ValueError('scaffold_by_position values must be A or B')
        stem_key, tracr_key = ('fixed_seq_2', 'fixed_seq_4') if variant == 'A' else ('fixed_seq_3', 'fixed_seq_5')
        profiles[position] = {
            'ibar_lengths': lengths[position],
            'grna_lengths': sorted({len(row['grna_seq']) for row in rows if row['position'] == position}),
            'scaffold_variant': variant,
            'stem': _dna(structure.get(stem_key, ''), stem_key),
            'tracr': _dna(structure.get(tracr_key, ''), tracr_key),
        }
    requested_distance = config.get('grna_matching', {}).get('max_lev_dist_ibar', 3)
    if isinstance(requested_distance, bool) or not isinstance(requested_distance, int) or requested_distance < 0:
        raise ValueError('max_lev_dist_ibar must be a nonnegative integer')
    unit_pairs = defaultdict(set)
    ibar_pairs = defaultdict(set)
    codebooks = defaultdict(set)
    pair_rows = defaultdict(dict)
    for row in rows:
        if row['pair_id']:
            unit_pairs[row['guide_unit_id']].add(row['pair_id'])
            ibar_pairs[(row['position'], row['ibar_seq'])].add(row['pair_id'])
            previous = pair_rows[row['pair_id']].get(row['position'])
            if previous and previous['guide_unit_id'] != row['guide_unit_id']:
                raise ValueError(f"pair_id {row['pair_id']} has multiple units at position {row['position']}")
            pair_rows[row['pair_id']][row['position']] = row
        codebooks[(row['position'], row['grna_seq'])].add(row['ibar_seq'])
    if positions == [1, 2] and any(set(pair) != {1, 2} for pair in pair_rows.values()):
        raise ValueError('Each dual pair_id must have exactly one unit at each position')
    correction = {}
    for key, barcodes in codebooks.items():
        for barcode in barcodes:
            safe_distance = 0 if len(barcode) <= 1 else requested_distance
            if len(barcodes) > 1024:
                safe_distance = 0
            for other in barcodes if safe_distance else ():
                if other != barcode:
                    safe_distance = min(safe_distance, max(0, (edit_distance(barcode, other) - 1) // 2))
            correction[unit_id(key[0], key[1], barcode)] = safe_distance
    if any(value < requested_distance for value in correction.values()):
        warn('ibar_correction_reduced', 'Reduced iBAR correction for short or insufficiently separated codewords; UMI correction is independent')
    inference = config.get('guide_assignment', {}).get('infer_missing_partner', 'auto')
    if inference not in ('auto', True, False):
        raise ValueError('infer_missing_partner must be auto, true or false')
    if inference is True and any(len(values) != 1 for values in unit_pairs.values()):
        warn('partner_inference_partial', 'Missing-partner inference is available only for units identifying one pair')
    units = sorted({row['guide_unit_id'] for row in rows})
    capabilities = {
        'units_total': len(units),
        'units_identifying_pair': sum(len(unit_pairs[key]) == 1 for key in units),
        'position_from_length': all(not set(lengths[first]) & set(lengths[second]) for first, second in combinations(positions, 2)),
        'umi_correction': True,
        'umi_grouping': 'observable_guide_unit',
        'infer_missing_partner': inference is not False,
        'legacy_recombination_model': (
            positions == [1, 2]
            and all(len(unit_pairs[key]) == 1 for key in units)
            and all(barcode and len(barcode) >= 2 and len(values) == 1 for (_, barcode), values in ibar_pairs.items())
            and len({tuple(values) for values in lengths.values()}) == 1
            and all(len(values) == 1 for values in lengths.values())
            and len({(profile['stem'], profile['tracr']) for profile in profiles.values()}) == 2
        ),
    }
    capabilities['legacy_recombination_model_reason'] = (
        'Symmetric, resolving iBAR identities retain the legacy observation categories'
        if capabilities['legacy_recombination_model'] else
        'The legacy likelihood requires symmetric, fully resolving iBAR categories; use observed design compatibility QC'
    )
    if config.get('guide_assignment', {}).get('umi_dedup_ibar') is False and any(row['ibar_seq'] for row in rows):
        warn('umi_grouping_override', 'Design-aware molecule grouping preserves iBAR identity even when legacy umi_dedup_ibar is false')
    signatures = {'ibar_tuple': defaultdict(set), 'spacer_tuple': defaultdict(set), 'full_tuple': defaultdict(set)}
    for pair_id, pair in pair_rows.items():
        ordered = [pair[position] for position in sorted(pair)]
        signatures['ibar_tuple'][tuple(row['ibar_seq'] for row in ordered)].add(pair_id)
        signatures['spacer_tuple'][tuple(row['grna_seq'] for row in ordered)].add(pair_id)
        signatures['full_tuple'][tuple(row['guide_unit_id'] for row in ordered)].add(pair_id)
    for label, lookup in signatures.items():
        capabilities[label] = {'signatures': len(lookup), 'unique': sum(len(values) == 1 for values in lookup.values())}
    capabilities['per_position'] = {
        position: {
            'ibar_lengths': lengths[position],
            'units_total': len({row['guide_unit_id'] for row in rows if row['position'] == position}),
            'units_identifying_pair': len({row['guide_unit_id'] for row in rows if row['position'] == position and len(unit_pairs[row['guide_unit_id']]) == 1}),
            'ibar_signatures_identifying_pair': sum(len(values) == 1 for (observed_position, barcode), values in ibar_pairs.items() if observed_position == position and barcode),
        } for position in positions
    }
    for position in positions:
        pool = {barcode for (observed_position, barcode) in ibar_pairs if observed_position == position and barcode}
        capabilities['per_position'][position]['decoupling_guide_groups_eligible'] = sum(
            bool(barcodes - {''}) and barcodes - {''} < pool
            for (observed_position, _), barcodes in codebooks.items() if observed_position == position
        )
    digest_rows = sorted((row['position'], row['grna_seq'], row['ibar_seq'], row['pair_id'], row['grna_name'], row['target_gene'],
                          *(str(row.get(field, '')) for field in ('target_combo', 'target_version', 'gi_pair', 'gi_version'))) for row in rows)
    matching_policy = {'max_lev_dist': config.get('grna_matching', {}).get('max_lev_dist', 3),
                       'tso': structure.get('fixed_seq_1', ''),
                       'tso_max_edit_dist': structure.get('tso_max_edit_dist', 3),
                       'fixed_seq_max_mismatches': structure.get('fixed_seq_max_mismatches', 2)}
    identity_capabilities = {key: value for key, value in capabilities.items() if key != 'infer_missing_partner'}
    digest = hashlib.sha256(json.dumps([SCHEMA_VERSION, digest_rows, profiles, correction, identity_capabilities, matching_policy], sort_keys=True).encode()).hexdigest()
    reference_digest = hashlib.sha256(json.dumps(digest_rows).encode()).hexdigest()
    return {'schema_version': SCHEMA_VERSION, 'fingerprint': digest, 'rows': rows,
            'reference_fingerprint': reference_digest,
            'profiles': profiles, 'grna_lengths': guide_lengths,
            'ibar_correction': correction, 'capabilities': capabilities,
            'requested_settings': {
                'ibar_length': structure.get('ibar_length'),
                'ibar_lengths_by_position': configured,
                'grna_lengths': structure.get('grna_lengths'),
                'max_lev_dist_ibar': requested_distance,
                'infer_missing_partner': inference,
                'umi_dedup_ibar': config.get('guide_assignment', {}).get('umi_dedup_ibar'),
            },
            'warnings': messages}


def emit_warnings(contract):
    for message in contract['warnings']:
        warnings.warn(f"[{message['code']}] {message['message']}", UserWarning, stacklevel=2)


def extract_unit(sequence, start, contract, max_mismatches=2, discover_lengths=None, right_anchor=None):
    """Measure a guide unit using independent stem and downstream boundaries."""
    sequence = sequence.upper()
    profiles = contract['profiles']
    lengths = sorted({length for profile in profiles.values() for length in profile['ibar_lengths']})
    if discover_lengths is not None:
        lengths = sorted(set(discover_lengths))
    scaffolds = {profile['scaffold_variant']: (profile['stem'], profile['tracr']) for profile in profiles.values()}
    right_scaffolds = {'unobserved': right_anchor.upper()} if right_anchor else {variant: sequences[1] for variant, sequences in scaffolds.items()}
    candidates = []
    for guide_length in contract['grna_lengths']:
        stem_start = start + guide_length
        for stem_variant, (stem, _) in scaffolds.items():
            stem_observed = sequence[stem_start:stem_start + len(stem)]
            if len(stem_observed) != len(stem):
                continue
            stem_distance = sum(first != second for first, second in zip(stem_observed, stem))
            if stem_distance > max_mismatches:
                continue
            ibar_start = stem_start + len(stem)
            for ibar_length in lengths:
                ibar_end = ibar_start + ibar_length
                for tracr_variant, tracr in right_scaffolds.items():
                    observed = sequence[ibar_end:ibar_end + len(tracr)]
                    if len(observed) != len(tracr):
                        continue
                    tracr_distance = sum(first != second for first, second in zip(observed, tracr))
                    score = stem_distance + tracr_distance
                    if score > max_mismatches:
                        continue
                    candidates.append((score, guide_length, ibar_length, stem_variant, tracr_variant,
                                       sequence[start:stem_start], sequence[ibar_start:ibar_end]))
    if not candidates:
        return None
    best_score = min(candidate[0] for candidate in candidates)
    best = sorted(set(candidate for candidate in candidates if candidate[0] == best_score))
    boundaries = {(candidate[1], candidate[2], candidate[5], candidate[6]) for candidate in best}
    if len(boundaries) != 1:
        return None
    score, guide_length, ibar_length, stem_variant, tracr_variant, spacer, ibar = best[0]
    stem_variants = {candidate[3] for candidate in best}
    tracr_variants = {candidate[4] for candidate in best}
    stem_variant = next(iter(stem_variants)) if len(stem_variants) == 1 else 'AMBIGUOUS'
    tracr_variant = next(iter(tracr_variants)) if len(tracr_variants) == 1 else 'AMBIGUOUS'
    length_positions = [position for position, profile in profiles.items() if ibar_length in profile['ibar_lengths']]
    scaffold_positions = [position for position, profile in profiles.items() if any(profile['scaffold_variant'] == candidate[3] and candidate[4] in (candidate[3], 'unobserved') for candidate in best)]
    length_position = length_positions[0] if len(length_positions) == 1 else 0
    scaffold_position = scaffold_positions[0] if len(scaffold_positions) == 1 else 0
    position = length_position or scaffold_position
    conflict = not scaffold_positions or bool(length_position and scaffold_position and length_position != scaffold_position)
    if not length_positions:
        conflict = True
    return {'grna_seq': spacer, 'ibar_seq': ibar, 'grna_length': guide_length,
            'ibar_length': ibar_length, 'position': position,
            'position_from_length': length_position, 'position_from_scaffold': scaffold_position,
            'scaffold_variant': stem_variant, 'tracr_variant': tracr_variant,
            'structure_status': 'conflict' if conflict else 'concordant',
            'variant': position or 'AMBIGUOUS', 'anchor_score': score}


class DesignMatcher:
    """Set-valued design lookup; one observation never expands into read copies."""

    def __init__(self, contract):
        self.contract = contract
        self.by_spacer = defaultdict(list)
        self.global_ibars = defaultdict(set)
        self.allowed_ibars = defaultdict(set)
        for row in contract['rows']:
            self.by_spacer[row['grna_seq']].append(row)
            if row['ibar_seq']:
                self.global_ibars[row['position']].add(row['ibar_seq'])
                self.allowed_ibars[(row['position'], row['grna_seq'])].add(row['ibar_seq'])

    def resolve(self, observation, spacer=None):
        spacer = spacer or observation['grna_seq']
        position = observation['position']
        candidates = [row for row in self.by_spacer.get(spacer, []) if not position or row['position'] == position]
        if not candidates:
            return None
        ibar = observation['ibar_seq']
        matches = []
        for row in candidates:
            distance = edit_distance(ibar, row['ibar_seq'])
            if distance <= self.contract['ibar_correction'][row['guide_unit_id']]:
                matches.append((distance, row))
        corrected = ''
        ibar_distance = -1
        status = 'present_unresolved'
        if matches:
            minimum = min(distance for distance, _ in matches)
            nearest = [row for distance, row in matches if distance == minimum]
            barcodes = {row['ibar_seq'] for row in nearest}
            if len(barcodes) == 1:
                candidates = nearest
                positions = {row['position'] for row in nearest}
                if len(positions) == 1:
                    position = next(iter(positions))
                corrected = next(iter(barcodes))
                ibar_distance = minimum
                status = 'present_resolved' if corrected else 'absent_by_design'
        names = {row['grna_name'] for row in candidates}
        targets = {row['target_gene'] for row in candidates}
        identity = unit_id(position, spacer, corrected if status != 'present_unresolved' else ibar)
        decoupled = bool(ibar and ibar in self.global_ibars[position] and (status == 'present_unresolved' or corrected != ibar))
        pair_ids = sorted({row['pair_id'] for row in candidates if row['pair_id']}) if status != 'present_unresolved' and observation['structure_status'] == 'concordant' and not decoupled else []
        result = dict(observation)
        result.update(guide_unit_id=identity, position=position,
                      ibar_decoupling_eligible=bool(ibar and ibar in self.global_ibars[position]
                                                    and self.allowed_ibars[(position, spacer)] < self.global_ibars[position]
                                                    and observation['structure_status'] == 'concordant'),
                      legacy_model_eligible=self.contract['capabilities']['legacy_recombination_model'],
                      guide_design_fingerprint=self.contract['fingerprint'],
                      grna_name=next(iter(names)) if len(names) == 1 else identity,
                      target_gene=next(iter(targets)) if len(targets) == 1 else '',
                      ibar_status=status, ibar_seq_corrected=corrected, ibar_dist=ibar_distance,
                      ibar_decoupled=decoupled,
                      pair_candidates=json.dumps(pair_ids, separators=(',', ':')),
                      pair_candidate_count=len(pair_ids),
                      pair_id=pair_ids[0] if len(pair_ids) == 1 else '',
                      target_combo=candidates[0]['target_combo'] if len(pair_ids) == 1 else '',
                      pair_resolution_basis='spacer_ibar_position' if pair_ids else '',
                      resolution_status=('structure_conflict' if observation['structure_status'] != 'concordant' else
                                         'unique_pair' if len(pair_ids) == 1 else 'ambiguous_pair' if pair_ids else 'guide_only'))
        return result


def linked_pair_status(first, second):
    """Classify a physically linked pair without requiring resolving iBARs."""
    if first is None or second is None:
        return {'status': 'unobserved', 'pair_candidates': []}
    if first['structure_status'] != 'concordant' or second['structure_status'] != 'concordant':
        return {'status': 'structure_conflict', 'pair_candidates': []}
    if first.get('ibar_decoupled') or second.get('ibar_decoupled'):
        return {'status': 'recombination_compatible', 'pair_candidates': []}
    first_pairs = set(json.loads(first['pair_candidates']))
    second_pairs = set(json.loads(second['pair_candidates']))
    overlap = sorted(first_pairs & second_pairs)
    if overlap:
        return {'status': 'designed_unique' if len(overlap) == 1 else 'designed_ambiguous', 'pair_candidates': overlap}
    eligible = first['ibar_status'] != 'present_unresolved' and second['ibar_status'] != 'present_unresolved' and bool(first_pairs and second_pairs)
    return {'status': 'recombination_compatible' if eligible else 'unresolved', 'pair_candidates': []}


def assign_cells(frame, reference, min_umi_threshold=5, infer_missing_partner=True):
    """Assign compatible pairs while allocating each observed unit at most once."""
    import pandas as pd

    pairs = defaultdict(dict)
    for row in reference.fillna('').to_dict('records'):
        pairs[str(row.get('pair_id', ''))][int(row['position'])] = row
    records = []

    def record(cell, rows, pair_id='', detection='ambiguous'):
        result = {'bc_wells': cell, 'pair_id': pair_id, 'target_combo': '', 'gi_pair': '',
                  'legacy_model_eligible': bool(rows) and all(str(row.get('legacy_model_eligible', False)).lower() in ('true', '1') for row in rows),
                  'detection_status': detection, 'inference_method': '',
                  'pair_resolution_basis': 'cell_compatibility' if detection == 'complete' else '',
                  'pair_candidate_count': 0, 'pair_candidates': '[]',
                  'unpositioned_umi': sum(row['umi_count'] for row in rows if int(row['position']) == 0),
                  'unpositioned_grna': ','.join(sorted({row['grna_name'] for row in rows if int(row['position']) == 0})),
                  'pair_umi_total': sum(row['umi_count'] for row in rows)}
        candidates = set()
        for row in rows:
            candidates.update(json.loads(row.get('pair_candidates', '[]')))
        if pair_id:
            candidates = {pair_id}
        result['pair_candidates'] = json.dumps(sorted(candidates), separators=(',', ':'))
        result['pair_candidate_count'] = len(candidates)
        for position in (1, 2):
            observed = [row for row in rows if int(row['position']) == position]
            result.update({f'pos{position}_grna': ','.join(sorted({row['grna_name'] for row in observed})),
                           f'pos{position}_ibar_status': ','.join(sorted({row.get('ibar_status', '') for row in observed})) if observed else 'unobserved',
                           f'pos{position}_target': ','.join(sorted({row['target_gene'] for row in observed})),
                           f'pos{position}_ibar': ','.join(sorted({row.get('ibar_seq_corrected', '') for row in observed})),
                           f'pos{position}_umi': sum(row['umi_count'] for row in observed),
                           f'pos{position}_umi_raw': sum(row.get('umi_count_raw', row['umi_count']) for row in observed),
                           f'inferred_pos{position}_grna': '', f'inferred_pos{position}_target': ''})
            if pair_id:
                expected = pairs[pair_id].get(position, {})
                result['target_combo'] = expected.get('target_combo', '')
                result['gi_pair'] = expected.get('gi_pair', '')
                if not observed and detection == 'inferred':
                    result[f'inferred_pos{position}_grna'] = expected.get('grna_name', '')
                    result[f'inferred_pos{position}_target'] = expected.get('target_gene', '')
                    result['inference_method'] = 'unique_observed_unit'
                    result['pair_resolution_basis'] = 'unique_observed_unit'
        return result

    for cell, group in frame.groupby('bc_wells', sort=True):
        rows = group.to_dict('records')
        candidate_sets = [set(json.loads(row.get('pair_candidates', '[]'))) if row.get('structure_status') == 'concordant' else set() for row in rows]
        all_pairs = set().union(*candidate_sets)
        complete = {pair_id for pair_id in all_pairs if {int(row['position']) for row, candidates in zip(rows, candidate_sets) if pair_id in candidates} == {1, 2}}
        allocated = set()
        cell_records = []
        for pair_id in sorted(complete):
            indices = [index for index, candidates in enumerate(candidate_sets) if pair_id in candidates and len(candidates & complete) == 1]
            observed = [rows[index] for index in indices]
            if {int(row['position']) for row in observed} == {1, 2} and sum(row['umi_count'] for row in observed) >= min_umi_threshold:
                cell_records.append(record(cell, observed, pair_id, 'complete'))
                allocated.update(indices)
        if infer_missing_partner:
            for pair_id in sorted(all_pairs):
                indices = [index for index, candidates in enumerate(candidate_sets) if index not in allocated and candidates == {pair_id} and int(rows[index]['position']) in (1, 2)]
                observed = [rows[index] for index in indices]
                if observed and sum(row['umi_count'] for row in observed) >= min_umi_threshold:
                    cell_records.append(record(cell, observed, pair_id, 'inferred'))
                    allocated.update(indices)
        remaining = [row for index, row in enumerate(rows) if index not in allocated and row['umi_count'] >= min_umi_threshold]
        for row in remaining:
            detection = 'discordant' if row.get('structure_status') != 'concordant' else 'ambiguous' if len(json.loads(row.get('pair_candidates', '[]'))) > 1 else 'guide_only'
            cell_records.append(record(cell, [row], detection=detection))
        resolved = sum(bool(row['pair_id']) for row in cell_records)
        uncertain = any(not row['pair_id'] for row in cell_records)
        for row in cell_records:
            row['n_combos'] = resolved
            row['combo_type'] = 'multiplet' if resolved > 1 else 'unresolved' if uncertain or not resolved else 'singlet'
            records.append(row)
    columns = list(record('', []).keys()) + ['n_combos', 'combo_type']
    return pd.DataFrame(records, columns=columns)


def eligible_assignments(frame):
    """Keep experimentally resolved construct calls for legacy QC consumers."""
    if 'pair_resolution_basis' not in frame.columns:
        return frame
    return frame[frame['pair_id'].fillna('').ne('') & frame['detection_status'].isin(['complete', 'inferred'])].copy()


def legacy_model_eligible(frame):
    if 'pair_resolution_basis' not in frame.columns and 'guide_unit_id' not in frame.columns:
        return True
    return bool(len(frame)) and 'legacy_model_eligible' in frame.columns and bool(frame['legacy_model_eligible'].astype(str).str.lower().isin(['true', '1']).all())

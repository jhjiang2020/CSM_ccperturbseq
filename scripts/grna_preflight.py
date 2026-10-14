"""Bounded, exact-lookup diagnostics before production guide matching."""

import csv
import gzip
from collections import Counter, defaultdict
from contextlib import closing
from itertools import islice


def preflight_options(matching_config):
    """Old configs enable a 10,000-read check without adding any keys."""
    options = matching_config.get('preflight', {})
    if not isinstance(options, dict):
        raise ValueError('grna_matching.preflight must be a mapping')
    enabled = options.get('enabled', True)
    sample_reads = options.get('sample_reads', 10_000)
    if not isinstance(enabled, bool):
        raise ValueError('grna_matching.preflight.enabled must be true or false')
    if isinstance(sample_reads, bool) or not isinstance(sample_reads, int) or not 1 <= sample_reads <= 100_000:
        raise ValueError('grna_matching.preflight.sample_reads must be an integer from 1 to 100000')
    return enabled, sample_reads


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def _barcode_sequences(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not reader.fieldnames or 'r1_seq' not in reader.fieldnames:
            raise ValueError(f'Preflight input lacks r1_seq: {path}')
        for row in reader:
            yield row['r1_seq']


def _fastq_sequences(paths):
    for path in paths:
        opener = gzip.open if str(path).endswith('.gz') else open
        with opener(path, 'rt') as handle:
            while True:
                header = handle.readline()
                if not header:
                    break
                sequence = handle.readline().rstrip('\r\n')
                plus = handle.readline()
                quality = handle.readline().rstrip('\r\n')
                if not header.startswith('@') or not plus.startswith('+') or len(sequence) != len(quality):
                    raise ValueError(f'Malformed FASTQ record during preflight: {path}')
                yield sequence


def check_ibar_reference(barcode_path, reference_pairs, config, extractor,
                         sample_reads=10_000, min_informative_reads=100,
                         fastq_paths=None):
    """Compare both iBAR strands globally and conditional on an exact guide.

    No fuzzy reference search, reference correction, or production matching is
    performed. Multiple iBARs and shared guides are represented as sets. Counts
    are read-weighted; this is a warning heuristic, not a coupling estimator.
    """
    result = {'status': 'inconclusive', 'sample_limit': sample_reads,
              'min_informative_reads': min_informative_reads, 'warnings': [],
              'input_type': 'raw_fastq' if fastq_paths is not None else 'barcode_tsv'}
    if config['ibar_length'] == 0:
        result.update(status='skipped', reason='ibar_length=0')
        return result

    def warn(code, message):
        result['warnings'].append({'code': code, 'message': message})

    global_ibars = set()
    guide_ibars = defaultdict(set)
    lengths = Counter()
    invalid = 0
    for guide, ibar in reference_pairs:
        if not isinstance(ibar, str) or not ibar or any(base not in 'ACGT' for base in ibar):
            invalid += 1
            continue
        lengths[len(ibar)] += 1
        global_ibars.add(ibar)
        if isinstance(guide, str):
            guide_ibars[guide].add(ibar)
    result['reference_ibar_length_counts'] = dict(lengths)
    result['unique_reference_ibars'] = len(global_ibars)
    if invalid:
        warn('invalid_reference_ibar', f'{invalid:,} reference rows have empty, non-string, or non-ACGT iBARs; check the formatted CSV.')
    if not global_ibars:
        warn('empty_reference_ibars', 'No valid reference iBARs are available for the configured iBAR-containing library.')
    if lengths and set(lengths) != {config['ibar_length']}:
        warn('ibar_length_mismatch',
             f"Configured ibar_length={config['ibar_length']} differs from reference lengths {sorted(lengths)}; check grna_structure and the formatted CSV.")

    global_rc = {reverse_complement(seq) for seq in global_ibars}
    guide_rc = {guide: {reverse_complement(seq) for seq in ibars}
                for guide, ibars in guide_ibars.items()}
    counts = Counter({key: 0 for key in (
        'sampled_reads', 'extracted_reads', 'exact_guide_reads',
        'global_stored_exact', 'global_rc_exact', 'global_stored_only', 'global_rc_only',
        'same_guide_stored_exact', 'same_guide_rc_exact',
        'same_guide_stored_only', 'same_guide_rc_only')})
    sequences = (_fastq_sequences(fastq_paths) if fastq_paths is not None
                 else _barcode_sequences(barcode_path))
    with closing(sequences):
        for sequence in islice(sequences, sample_reads):
            counts['sampled_reads'] += 1
            extracted = extractor(
                sequence, config['tso_seq'], config['stem1'], config['stem2'],
                config['grna_lengths'], config['ibar_length'],
                config['tracr1'], config['tracr2'], config['anchor_max_mm'],
                config['tso_max_dist'],
            )
            if extracted is None:
                continue
            counts['extracted_reads'] += 1
            observed = extracted['ibar_seq']
            stored, rc = observed in global_ibars, observed in global_rc
            counts['global_stored_exact'] += stored
            counts['global_rc_exact'] += rc
            counts['global_stored_only'] += stored and not rc
            counts['global_rc_only'] += rc and not stored
            guide = extracted['grna_seq']
            if guide not in guide_ibars:
                continue
            counts['exact_guide_reads'] += 1
            stored, rc = observed in guide_ibars[guide], observed in guide_rc[guide]
            counts['same_guide_stored_exact'] += stored
            counts['same_guide_rc_exact'] += rc
            counts['same_guide_stored_only'] += stored and not rc
            counts['same_guide_rc_only'] += rc and not stored
    result['counts'] = dict(counts)

    # Exclusive hits avoid treating palindromes or an RC-closed global pool as
    # evidence of reversal. Same-guide evidence still works for an RC-closed pool.
    evidence = []
    for prefix, denominator in (('same_guide', counts['exact_guide_reads']),
                                ('global', counts['extracted_reads'])):
        stored = counts[f'{prefix}_stored_only']
        rc = counts[f'{prefix}_rc_only']
        if (denominator >= min_informative_reads and rc >= 20
                and rc >= 5 * max(stored, 1) and (rc - stored) / denominator >= 0.10):
            evidence.append(prefix)
    if evidence:
        result['orientation_evidence'] = evidence
        warn('possible_reverse_complement',
             'Possible iBAR reference strand misspecification: reverse-complement support greatly exceeds stored-strand support '
             f"({', '.join(evidence)}). Check ibar_seq in paths.grna_reference and design.rev_comp_ibar. "
             'Step 3 reads the existing CSV; changing the design flag alone does not update it. No automatic correction was made.')

    n = counts['extracted_reads']
    if n >= min_informative_reads and global_ibars:
        if max(counts['global_stored_exact'], counts['global_rc_exact']) / n < 0.10:
            warn('low_ibar_support',
                 'Fewer than 10% of extracted iBARs exactly match either reference orientation; check library/reference identity, iBAR length/boundary, fixed anchors, and read quality.')
    if (counts['sampled_reads'] >= min_informative_reads
            and counts['extracted_reads'] / counts['sampled_reads'] < 0.10):
        warn('low_extraction_support',
             'Fewer than 10% of sampled reads passed structure extraction; iBAR orientation evidence is limited. Check fixed sequences and available R1 length.')
    if result['warnings']:
        result['status'] = 'warning'
    elif n >= min_informative_reads:
        # Exact strand evidence may be uninformative even with many structures.
        if counts['global_stored_only'] + counts['global_rc_only'] >= 20 or counts['same_guide_stored_only'] + counts['same_guide_rc_only'] >= 20:
            result['status'] = 'ok'
    if result['status'] == 'inconclusive':
        result['reason'] = 'Insufficient informative exact iBAR evidence; no strand conclusion.'
    return result


def print_preflight(result, library):
    print(f"  iBAR preflight [{library}]: {result['status']}", flush=True)
    if 'counts' in result:
        c = result['counts']
        print(f"    Sampled {c['sampled_reads']:,}; extracted {c['extracted_reads']:,}; exact guide {c['exact_guide_reads']:,}", flush=True)
        print(f"    Global exact iBAR: stored={c['global_stored_exact']:,}, reverse-complement={c['global_rc_exact']:,}", flush=True)
        print(f"    Same-guide exact iBAR: stored={c['same_guide_stored_exact']:,}, reverse-complement={c['same_guide_rc_exact']:,}", flush=True)
    for warning in result['warnings']:
        print(f"    WARNING [{warning['code']}]: {warning['message']}", flush=True)
    if 'reason' in result:
        print(f"    {result['reason']}", flush=True)

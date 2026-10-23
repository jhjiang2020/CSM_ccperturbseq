#!/usr/bin/env python3
"""Build fixed-offset, pilot-derived 7/4-nt fixtures and validate the gRNA arm."""

import argparse
from collections import Counter, defaultdict
import copy
import csv
import gzip
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import socket
import subprocess
import sys

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]


def digest(value):
    return hashlib.sha256(value.encode()).hexdigest()


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def fastq_records(path):
    with gzip.open(path, 'rt') as handle:
        ordinal = 0
        while True:
            header = handle.readline().rstrip('\r\n')
            if not header:
                break
            sequence = handle.readline().rstrip('\r\n')
            plus = handle.readline().rstrip('\r\n')
            quality = handle.readline().rstrip('\r\n')
            ordinal += 1
            if not header.startswith('@') or not plus.startswith('+') or len(sequence) != len(quality):
                raise ValueError(f'Malformed FASTQ: {path}, record {ordinal}')
            yield ordinal, header, sequence, quality


def fixed_observation(sequence, structure, spacer_lengths, tso_offsets=(0, 7)):
    """Use exact, prespecified offsets only; never call the production extractor."""
    observations = []
    tso = structure['fixed_seq_1']
    for offset in tso_offsets:
        if sequence[offset:offset + len(tso)] != tso:
            continue
        start = offset + len(tso)
        for length in spacer_lengths:
            for position in (1, 2):
                stem = structure[f'fixed_seq_{position + 1}']
                tracr = structure[f'fixed_seq_{position + 3}']
                stem_start = start + length
                ibar_start = stem_start + len(stem)
                ibar_end = ibar_start + 7
                if (sequence[stem_start:ibar_start] == stem
                        and sequence[ibar_end:ibar_end + len(tracr)] == tracr):
                    observations.append(dict(position=position, grna_seq=sequence[start:stem_start],
                                             ibar_seq=sequence[ibar_start:ibar_end],
                                             ibar_start=ibar_start, tso_offset=offset))
    return observations[0] if len(observations) == 1 else None


def shorten_ibar(sequence, quality, ibar_start, position):
    """Delete the last three expected iBAR2 bases, including their qualities."""
    if len(sequence) != len(quality):
        raise ValueError('Sequence and quality lengths differ')
    if position != 2:
        return sequence, quality
    start, end = ibar_start + 4, ibar_start + 7
    if start < 0 or end > len(sequence):
        raise ValueError('iBAR deletion lies outside the read')
    return sequence[:start] + sequence[end:], quality[:start] + quality[end:]


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def build_fixture(config_path, destination, pairs_per_library=20, max_reads=100_000):
    if destination.exists():
        raise FileExistsError(f'Refusing to replace fixture: {destination}')
    config = yaml.safe_load(config_path.read_text())
    reference_path = Path(config['paths']['grna_reference'])
    reference = pd.read_csv(reference_path, keep_default_na=False)
    lookup = defaultdict(list)
    for row in reference.to_dict('records'):
        for orientation, barcode in [('as_reference', row['ibar_seq']),
                                     ('reverse_complement', reverse_complement(row['ibar_seq']))]:
            lookup[(int(row['position']), row['grna_seq'], barcode)].append((row, orientation))
    selected, source_summary = [], []
    for library in config['libraries']:
        name = library['sublibrary_name']
        paths = sorted(Path(config['paths']['fastq_dir']).glob(f"*{library['grna_prefix']}*_R1_*.fastq.gz"))
        if not paths:
            raise FileNotFoundError(f'No gRNA R1 FASTQ for {name}')
        pools = defaultdict(list)
        counters = Counter()
        consumed = []
        for path in paths:
            scanned = 0
            for ordinal, header, sequence, quality in fastq_records(path):
                if counters['scanned_reads'] >= max_reads:
                    break
                scanned += 1
                counters['scanned_reads'] += 1
                observation = fixed_observation(sequence, config['grna_structure'], [20, 21])
                if observation is None:
                    continue
                counters['exact_fixed_layout'] += 1
                matches = lookup[(observation['position'], observation['grna_seq'], observation['ibar_seq'])]
                unique = {entry[0]['pair_id'] for entry in matches}
                if len(unique) != 1:
                    continue
                row, orientation = matches[0]
                counters[orientation] += 1
                key = (row['pair_id'], observation['position'], orientation)
                if len(pools[key]) < 3:
                    pools[key].append(dict(library=name, source_file=str(path), source_ordinal=ordinal,
                                           source_id_sha256=digest(header), source_record_sha256=digest(sequence + '\n' + quality),
                                           sequence=sequence, quality=quality, reference=row,
                                           source_reference_orientation=orientation, **observation))
            consumed.append(dict(path=str(path), records_scanned=scanned, size_bytes=path.stat().st_size))
            if counters['scanned_reads'] >= max_reads:
                break
        orientation = 'as_reference' if counters['as_reference'] >= counters['reverse_complement'] else 'reverse_complement'
        pools = {(pair, position): records for (pair, position, strand), records in pools.items() if strand == orientation}
        eligible = sorted(pair for pair, position in pools if position == 1
                          and len(pools.get((pair, 1), [])) == 3 and len(pools.get((pair, 2), [])) == 3)
        chosen = []
        for length in (20, 21):
            candidates = [pair for pair in eligible if any(len(pools[(pair, position)][0]['grna_seq']) == length for position in (1, 2))]
            if not candidates:
                raise ValueError(f'{name}: no complete construct covering spacer length {length}')
            if candidates[0] not in chosen:
                chosen.append(candidates[0])
        chosen.extend(pair for pair in eligible if pair not in chosen)
        chosen = chosen[:pairs_per_library]
        if len(chosen) < pairs_per_library:
            raise ValueError(f'{name}: only {len(chosen)} eligible pairs; increase --max-source-reads')
        for pair in chosen:
            for position in (1, 2):
                selected.extend(pools[(pair, position)])
        source_summary.append(dict(library=name, counters=dict(counters), files=consumed, selected_orientation=orientation,
                                   eligible_pairs=len(eligible), chosen_pairs=chosen))
        print(f'{name}: scanned {counters["scanned_reads"]:,}; selected {len(chosen)} pairs', flush=True)
    destination.mkdir(parents=True)
    source_rows = {}
    for entry in selected:
        row = dict(entry['reference'])
        row['ibar_seq'] = entry['ibar_seq']
        key = (row['pair_id'], int(row['position']))
        if key in source_rows and source_rows[key]['ibar_seq'] != row['ibar_seq']:
            raise ValueError('Conflicting observed iBAR orientations across libraries')
        source_rows[key] = row
    symmetric = pd.DataFrame([source_rows[key] for key in sorted(source_rows)])
    asymmetric = symmetric.copy()
    asymmetric.loc[asymmetric.position == 2, 'ibar_seq'] = asymmetric.loc[asymmetric.position == 2, 'ibar_seq'].str[:4]
    symmetric.to_csv(destination / 'reference_7_7.csv', index=False)
    asymmetric.to_csv(destination / 'reference_7_4.csv', index=False)
    observations, cells, provenance = [], [], []
    handles = {name: gzip.open(destination / f'{"control" if name == "symmetric" else name}_R1_001.fastq.gz', 'wt')
               for name in ('symmetric', 'asymmetric')}
    writers = {}
    evidence_handles = []
    try:
        for name in handles:
            handle = gzip.open(destination / f'{name}_barcodes.tsv.gz', 'wt')
            evidence_handles.append(handle)
            writers[name] = csv.writer(handle, delimiter='\t')
            writers[name].writerow(['read_id', 'bc_wells', 'bc_bci', 'umi', 'r1_seq', 'r1_qual'])
        grouped = defaultdict(list)
        for entry in selected:
            grouped[(entry['library'], entry['reference']['pair_id'], entry['position'])].append(entry)
        constructs = sorted({key[:2] for key in grouped})
        for index, (library, pair) in enumerate(constructs, 1):
            for mode, positions in [('complete', (1, 2)), ('pos1_only', (1,)), ('pos2_only', (2,))]:
                cell = f'{index:03d}_{len(positions):02d}_{positions[0]:02d}'
                cells.append(dict(bc_wells=cell, source_library=library, pair_id=pair, mode=mode,
                                  expected_umi=3 * len(positions)))
                for position in positions:
                    for molecule, entry in enumerate(grouped[(library, pair, position)]):
                        umi = ['A' * 12, 'C' * 12, 'G' * 12][molecule]
                        for copy_number in (1, 2):
                            read_id = f'pilot_{index:03d}_{mode}_p{position}_u{molecule}_r{copy_number}'
                            seq, qual = shorten_ibar(entry['sequence'], entry['quality'], entry['ibar_start'], position)
                            for name, sequence, quality in [('symmetric', entry['sequence'], entry['quality']), ('asymmetric', seq, qual)]:
                                handles[name].write(f'@{read_id}\n{sequence}\n+\n{quality}\n')
                                writers[name].writerow([read_id, cell, 'synthetic', umi, sequence, quality])
                            observations.append(dict(read_id=read_id, bc_wells=cell, umi=umi, position=position,
                                                     grna_seq=entry['grna_seq'], original_ibar=entry['ibar_seq'],
                                                     ibar_seq=entry['ibar_seq'][:4] if position == 2 else entry['ibar_seq'],
                                                     source_library=library, pair_id=pair))
                            provenance.append(dict(read_id=read_id, **{key: entry[key] for key in (
                                'source_file', 'source_ordinal', 'source_id_sha256', 'source_record_sha256',
                                'source_reference_orientation', 'ibar_start', 'tso_offset')},
                                deleted_start_0based=entry['ibar_start'] + 4 if position == 2 else None,
                                deleted_end_0based=entry['ibar_start'] + 7 if position == 2 else None))
    finally:
        for handle in list(handles.values()) + evidence_handles:
            handle.close()
    pd.DataFrame(observations).to_csv(destination / 'read_truth.tsv.gz', sep='\t', index=False)
    pd.DataFrame(provenance).to_csv(destination / 'provenance.tsv.gz', sep='\t', index=False)
    pd.DataFrame(cells).to_csv(destination / 'cell_truth.csv', index=False)
    fixture_config = {key: copy.deepcopy(config[key]) for key in ('grna_structure', 'grna_matching', 'guide_assignment')}
    fixture_config['guide_design'] = {'enabled': True, 'ibar_orientation': 'auto'}
    fixture_config['grna_structure']['ibar_lengths_by_position'] = 'auto'
    fixture_config['grna_matching']['save_failed_reads'] = True
    fixture_config['grna_matching']['preflight'] = {'enabled': True, 'sample_reads': 10_000}
    fixture_config['guide_assignment'].update(min_umi_threshold=1, max_umi_edits=1, compute_saturation_curve=False)
    fixture_config['barcode_extraction'] = {'n_processes': 1}
    fixture_config['libraries'] = [{'sublibrary_name': 'POOL', 'grna_prefix': 'asymmetric', 'sample_name': 'Synthetic'}]
    (destination / 'config_template.yaml').write_text(yaml.safe_dump(fixture_config, sort_keys=False))
    manifest = dict(schema_version=1, source_config_sha256=digest(config_path.read_text()),
                    source_reference=str(reference_path), source_reference_sha256=digest(reference_path.read_text()),
                    sources=source_summary, read_count=len(observations), synthetic_cells=len(cells),
                    unique_source_records=len(selected), pairs_in_reference=len(symmetric) // 2,
                    source_spacer_lengths=sorted(symmetric.grna_seq.str.len().unique().tolist()),
                    fixed_tso_offsets=[0, 7], original_ibar_length=7, derived_ibar_lengths={'pos1': 7, 'pos2': 4},
                    manipulation='Delete [ibar_start+4, ibar_start+7) from pos2 sequence and quality; retain all other bases',
                    caveat='Selected exact-layout reads; artificial cells, 3 UMIs per arm and 2 PCR copies per UMI. Not a biological recombination or sensitivity estimate.')
    manifest['files_sha256'] = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(destination.iterdir()) if path.is_file()}
    write_json(destination / 'manifest.json', manifest)
    return manifest


def validate_fixture(fixture, output, workers):
    if output.exists():
        raise FileExistsError(f'Refusing to replace validation: {output}')
    output.mkdir(parents=True)
    manifest = json.loads((fixture / 'manifest.json').read_text())
    for name, expected in manifest['files_sha256'].items():
        if hashlib.sha256((fixture / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f'Fixture checksum mismatch: {name}')
    truth = pd.read_csv(fixture / 'read_truth.tsv.gz', sep='\t', keep_default_na=False).set_index('read_id')
    cells = pd.read_csv(fixture / 'cell_truth.csv', keep_default_na=False).set_index('bc_wells')
    results = {'job_id': os.environ.get('JOB_ID'), 'host': socket.gethostname(),
               'python': sys.version, 'executable': sys.executable,
               'conda_prefix': os.environ.get('CONDA_PREFIX'), 'platform': platform.platform(),
               'workers': workers, 'fixture_sha256': digest((fixture / 'manifest.json').read_text()),
               'packages': {}, 'cases': {}, 'commands': []}
    for package in ('numpy', 'pandas', 'scipy', 'edlib', 'PyYAML', 'pytest', 'scanpy', 'anndata'):
        try:
            results['packages'][package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            results['packages'][package] = None
    results['source_sha256'] = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                for path in sorted((ROOT / 'scripts').glob('*.py'))}

    def run(label, arguments):
        command = [sys.executable, '-B', *map(str, arguments)]
        with (output / f'{label}.log').open('w') as handle:
            completed = subprocess.run(command, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT, timeout=1800)
        results['commands'].append({'label': label, 'argv': command, 'exit_code': completed.returncode})
        write_json(output / 'results.json', results)
        if completed.returncode:
            raise RuntimeError(f'{label} failed; see {output / (label + ".log")}')

    matched_outputs = {}
    for case in ('legacy_7_7', 'design_7_7', 'design_7_4'):
        asymmetric = case == 'design_7_4'
        prefix = 'asymmetric' if asymmetric else 'symmetric'
        config = yaml.safe_load((fixture / 'config_template.yaml').read_text())
        reference_path = fixture / ('reference_7_4.csv' if asymmetric else 'reference_7_7.csv')
        config['paths'] = {'grna_reference': str(reference_path), 'output_dir': str(output / case), 'fastq_dir': str(fixture)}
        config['libraries'][0]['grna_prefix'] = 'asymmetric' if asymmetric else 'control'
        if case == 'legacy_7_7':
            config.pop('guide_design')
            config['grna_structure'].pop('ibar_lengths_by_position')
        config_path = output / f'{case}.yaml'
        config_path.write_text(yaml.safe_dump(config, sort_keys=False))
        run(f'{case}_raw_pilot', [ROOT / 'scripts/03_match_grna_ibar.py', '--config', config_path,
                                '--preflight-only', '--preflight-fastq'])
        raw_pilot_path = output / case / 'grna_matches/POOL_preflight.yaml'
        raw_pilot = yaml.safe_load(raw_pilot_path.read_text())
        (output / f'{case}_raw_pilot.yaml').write_text(yaml.safe_dump(raw_pilot, sort_keys=False))
        run(f'{case}_step3', [ROOT / 'scripts/03_match_grna_ibar.py', '--config', config_path,
                            '--barcode-input', fixture / f'{prefix}_barcodes.tsv.gz',
                            '--n_processes', workers, '--chunk_size', 100])
        run(f'{case}_step4', [ROOT / 'scripts/04_assign_grna_to_cells.py', '--config', config_path, '--ncpus', workers])
        directory = output / case
        matched = pd.read_csv(directory / 'grna_matches/POOL_grna_matches.tsv.gz', sep='\t', keep_default_na=False).set_index('read_id').sort_index()
        if not matched.index.is_unique or set(matched.index) != set(truth.index):
            raise AssertionError(f'{case}: lost or duplicated reads')
        ordered = truth.loc[matched.index]
        expected_ibar = ordered['ibar_seq'] if asymmetric else ordered['original_ibar']
        for column, expected in [('grna_seq', ordered['grna_seq']), ('ibar_seq', expected_ibar)]:
            if matched[column].tolist() != expected.tolist():
                raise AssertionError(f'{case}: incorrect {column} boundaries')
        molecules = pd.read_csv(directory / 'grna_assignment/grna_umi_counts.csv.gz', keep_default_na=False)
        calls = pd.read_csv(directory / 'grna_assignment/grna_cell_assignment.csv', keep_default_na=False)
        expected_cells = {f'0_{key}': int(row['expected_umi']) for key, row in cells.iterrows()}
        if molecules.groupby('bc_wells')['umi_count'].sum().to_dict() != expected_cells:
            raise AssertionError(f'{case}: molecule counts differ from synthetic truth')
        if calls.groupby('bc_wells')['pair_umi_total'].sum().to_dict() != expected_cells:
            raise AssertionError(f'{case}: assignment loses or inflates molecules')
        reference = pd.read_csv(reference_path, keep_default_na=False)
        lookup = defaultdict(set)
        for row in reference.to_dict('records'):
            lookup[(int(row['position']), row['grna_seq'], row['ibar_seq'])].add(row['pair_id'])
        candidate_sets = {read_id: lookup[(int(row['position']), row['grna_seq'], expected_ibar.loc[read_id])]
                          for read_id, row in ordered.iterrows()}
        if case != 'legacy_7_7':
            for read_id, row in matched.iterrows():
                if set(json.loads(row['pair_candidates'])) != candidate_sets[read_id] or int(row['position']) != int(ordered.loc[read_id, 'position']):
                    raise AssertionError(f'{case}: candidate identity or position differs from independent truth')
            contract = yaml.safe_load((directory / 'grna_matches/POOL_design.yaml').read_text())
            if asymmetric and contract['capabilities']['legacy_recombination_model']:
                raise AssertionError('Asymmetric design incorrectly enables the legacy recombination model')
        for cell, group in ordered.groupby('bc_wells'):
            compatible = set.intersection(*(candidate_sets[read_id] for read_id in group.index))
            resolved = set(calls.loc[calls.bc_wells == f'0_{cell}', 'pair_id']) - {''}
            if not resolved <= compatible:
                raise AssertionError(f'{case}: unsupported construct assignment')
        pilot = raw_pilot
        if asymmetric and {int(position): set(counts) for position, counts in pilot['observed_ibar_length_counts'].items()} != {1: {7}, 2: {4}}:
            raise AssertionError('Raw FASTQ pilot did not independently recover the 7/4-nt design')
        results['cases'][case] = dict(reads=len(matched), cells=len(expected_cells), umi_count=int(molecules.umi_count.sum()),
                                    extraction_errors=0, count_errors=0, candidate_errors=0,
                                    detection_status_counts=calls.detection_status.value_counts().to_dict(),
                                    pilot_status=pilot['status'],
                                    observed_ibar_lengths=pilot.get('observed_ibar_length_counts'),
                                    orientation=pilot.get('orientation'))
        matched_outputs[case] = matched
        write_json(output / 'results.json', results)
    legacy = matched_outputs['legacy_7_7']
    pd.testing.assert_frame_equal(legacy, matched_outputs['design_7_7'][legacy.columns], check_dtype=False)
    results['symmetric_legacy_equivalence'] = True
    tests = ['test_guide_design.py', 'test_step3_matching.py', 'test_grna_preflight.py',
             'test_step1c_adapters.py', 'test_step1c_evidence_merge.py', 'test_gex_grna_filter.py',
             'test_gex_grna_preprocess_hook.py', 'test_gex_grna_bam_filter.py', 'test_asymmetric_pilot.py']
    results['test_sha256'] = {name: hashlib.sha256((ROOT / 'tests' / name).read_bytes()).hexdigest() for name in tests}
    regression_code = (
        'import json, sys, unittest; '
        f'names = {tests!r}; '
        'loader = unittest.TestLoader(); '
        'suite = unittest.TestSuite(loader.discover("tests", pattern=name) for name in names); '
        'result = unittest.TextTestRunner(verbosity=2).run(suite); '
        f'json.dump(dict(tests=result.testsRun, failures=len(result.failures), errors=len(result.errors), '
        f'skipped=len(result.skipped)), open({str(output / "regression_summary.json")!r}, "w"), indent=2); '
        'sys.exit(not result.wasSuccessful())'
    )
    run('regression_tests', ['-c', regression_code])
    results['regression_tests'] = json.loads((output / 'regression_summary.json').read_text())
    results['passed'] = True
    write_json(output / 'results.json', results)
    print(json.dumps(results['cases'], indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest='command', required=True)
    build = subparsers.add_parser('build')
    build.add_argument('--config', type=Path, required=True)
    build.add_argument('--destination', type=Path, required=True)
    build.add_argument('--pairs-per-library', type=int, default=20)
    build.add_argument('--max-source-reads', type=int, default=100_000)
    validate = subparsers.add_parser('validate')
    validate.add_argument('--fixture', type=Path, required=True)
    validate.add_argument('--output', type=Path, required=True)
    validate.add_argument('--workers', type=int, default=1)
    args = parser.parse_args()
    if args.command == 'build':
        if args.pairs_per_library < 1 or args.max_source_reads < 1:
            parser.error('Read and pair limits must be positive')
        build_fixture(args.config, args.destination, args.pairs_per_library, args.max_source_reads)
    elif args.command == 'validate':
        if args.workers < 1:
            parser.error('--workers must be positive')
        validate_fixture(args.fixture.resolve(), args.output.resolve(), args.workers)


if __name__ == '__main__':
    main()

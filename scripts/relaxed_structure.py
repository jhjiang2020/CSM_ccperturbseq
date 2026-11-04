"""One frozen relaxed-structure policy for every sublibrary and evidence source."""
import copy
import csv
import gzip
import hashlib
import json
import os
from collections import Counter
from pathlib import Path

try:
    from scripts import grna_recovery as recovery
    from scripts.design_preflight import check_design
    from scripts.guide_design import compile_design, reverse_complement
except ModuleNotFoundError:
    import grna_recovery as recovery
    from design_preflight import check_design
    from guide_design import compile_design, reverse_complement

COLUMNS = ['grna_seq_corrected', 'structure_correction']


def active(config):
    return config.get('grna_structure', {}).get('relaxed_structure', False) is not False


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f'.{path.name}.{os.getpid()}.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    os.replace(temporary, path)


def policy_path(config):
    return Path(config['paths']['output_dir']) / 'grna_matches/relaxed_structure.json'


def barcode_path(config, library):
    base = Path(config['paths']['output_dir']) / 'grna_barcodes'
    for suffix in ('.tsv.gz', '.tsv'):
        path = base / f'{library}_barcodes{suffix}'
        if path.exists():
            return path
    raise FileNotFoundError(f'Missing dedicated barcode evidence for {library}; run Step 2 first')


def apply_worker(worker, design):
    """Keep the existing parser, with one shared reference orientation."""
    rows = design['rows']
    if worker.get('design_contract'):
        worker['design_contract'] = design
    else:
        worker['ibar_seqs'] = list(dict.fromkeys(r['ibar_seq'] for r in rows))
        worker['ibar_names'] = worker['ibar_seqs']
        worker['combined_lookup'] = {(r['grna_seq'], r['ibar_seq']):
            (r['grna_name'], r['target_gene'], r.get('pair_id', ''), r.get('target_combo', '')) for r in rows}
        by_name = {}
        for row in rows:
            by_name.setdefault(row['grna_name'], []).append(row['ibar_seq'])
        worker['grna_to_ibars'] = by_name
    return worker


def choose_shared(reports, options):
    """Intersect validated policies; maximize pooled recovery within that set."""
    tables = [{recovery.fingerprint(s['policy']): s for s in r['settings']} for r in reports]
    keys = set(tables[0]) if tables else set()
    for table in tables:
        keys &= {key for key, value in table.items() if value['passed']}
    passing = [(sum(t[key]['discovery_identified_families'] for t in tables), tables[0][key]['policy']) for key in keys]
    if not passing:
        return dict(max_loss=0, min_retained=18, max_tso_shift=0,
                    endpoint_policy='joint', ibar_corroboration=options['ibar_corroboration'], enabled=False)
    best = max(n for n, _ in passing)
    close = [p for n, p in passing if n >= best * options['near_optimal_fraction']]
    selected = min(close, key=lambda p: (p['max_loss'], p.get('max_tso_shift', 0), -p['min_retained']))
    return dict(selected, enabled=True)


def build(config, reference, worker, initialize, baseline_processor, find_tso):
    """Run once on compute, after all Step-2 libraries are available."""
    options = recovery.settings(config)
    if not options.get('shared') or options['mode'] == 'off':
        return None
    destination = policy_path(config)
    original = compile_design(reference, config)
    source_paths = {lib['sublibrary_name']: barcode_path(config, lib['sublibrary_name']) for lib in config['libraries']}
    inputs = {name: recovery.source_stamp(path) for name, path in source_paths.items()}
    required = [p for p in original['profiles'].values() if not p['stem'] or not p['tracr']]
    samples, sampling, orientations = {}, {}, {}
    if options['mode'] == 'auto' and not required:
        for name, path in source_paths.items():
            cache = path.parent / f'{name}_recovery_sample.tsv.gz'
            rows, stamp = recovery.sample_rows(path, options['sample_reads'], options['seed'], cache)
            samples[name], sampling[name] = rows, stamp
            # Orientation checks use the same pilot, with one representative
            # per family, so PCR copies do not supply independent votes.
            families = {}
            for row in rows:
                key = row['bc_wells'], row['umi']
                if key not in families or row['read_id'] < families[key]['read_id']:
                    families[key] = row
            orient_path = destination.parent / f'{name}_structure_pilot.tsv.gz'
            orient_path.parent.mkdir(parents=True, exist_ok=True)
            with gzip.open(orient_path, 'wt') as handle:
                fields = ['read_id', 'bc_wells', 'bc_bci', 'umi', 'r1_seq', 'r1_qual']
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
                writer.writeheader()
                writer.writerows(sorted(families.values(), key=lambda row: recovery.fingerprint([row['bc_wells'], row['umi']])))
            diagnostic_config = copy.deepcopy(config)
            diagnostic_config.setdefault('guide_design', {})['ibar_orientation'] = 'auto'
            report, _ = check_design(orient_path, original, diagnostic_config, find_tso,
                                     sample_reads=options['sample_reads'])
            orientations[name] = report
            print(f'Shared pilot [{name}]: {len(rows):,} reads, {len(families):,} families; orientation {report["orientation"]}', flush=True)
    design = copy.deepcopy(original)
    chosen_orientation = {}
    for position in original['profiles']:
        conclusions = {report['orientation'][position] for report in orientations.values()
                       if report['orientation'][position] in ('as_reference', 'reverse_complement')}
        if len(conclusions) > 1:
            raise ValueError(f'Sublibraries disagree on iBAR orientation at position {position}; correct the reference/input configuration before sharing a policy')
        orientation = next(iter(conclusions), 'as_reference')
        chosen_orientation[position] = orientation
        if orientation == 'reverse_complement':
            for row in design['rows']:
                if row['position'] == position:
                    row['ibar_seq'] = reverse_complement(row['ibar_seq'])
    design = compile_design(design['rows'], config)
    design['input_fingerprint'] = original['fingerprint']
    design['resolved_orientation'] = chosen_orientation
    contract = recovery.recovery_contract(reference, config, design)
    reports = {}
    if required:
        effective = dict(max_loss=0, min_retained=18, max_tso_shift=0, endpoint_policy='joint',
                         ibar_corroboration=options['ibar_corroboration'], enabled=False)
        status = 'unavailable_missing_downstream_anchor'
    elif options['mode'] == 'auto':
        shared_worker = apply_worker(copy.deepcopy(worker), design)
        initialize(shared_worker)
        for name, rows in samples.items():
            baseline_ids = {r[0] for r in baseline_processor([
                tuple(row.get(k, '') for k in ('read_id', 'bc_wells', 'bc_bci', 'umi', 'r1_seq')) for row in rows])[0]}
            reports[name] = recovery.preflight(rows, contract, config, options, sampling[name], baseline_ids)['preflight']
            print(f'Shared pilot [{name}]: validated {sum(s["passed"] for s in reports[name]["settings"])} candidate settings', flush=True)
        effective = choose_shared(list(reports.values()), options)
        status = 'validated' if effective['enabled'] else 'inconclusive'
    else:
        effective = recovery.make_policy(contract, options)['effective']
        effective['enabled'] = True
        status = 'specified'
    policy = recovery.make_policy(contract, options, effective)
    policy.update(schema_version=2, design_contract=design, input_design_fingerprint=original['fingerprint'],
                  libraries=list(source_paths), dedicated_inputs=inputs)
    policy['fingerprint'] = recovery.fingerprint({k: v for k, v in policy.items() if k != 'fingerprint'})
    policy['preflight'] = dict(status=status, libraries=reports, orientation=orientations,
                               note='One shared policy; every sublibrary must pass. Expression is never used for selection.')
    atomic_json(destination, policy)
    print(f'Shared relaxed structure: {status}; {effective}', flush=True)
    return policy


def load(config, reference):
    path = policy_path(config)
    if not path.is_file():
        raise ValueError('Shared relaxed-structure policy is missing; run Step 3 --preflight-only across all sublibraries first')
    policy = json.loads(path.read_text())
    body = {k: v for k, v in policy.items() if k not in ('fingerprint', 'preflight')}
    if policy.get('schema_version') != 2 or recovery.fingerprint(body) != policy.get('fingerprint'):
        raise ValueError('Invalid shared relaxed-structure policy fingerprint')
    if compile_design(reference, config)['fingerprint'] != policy['input_design_fingerprint']:
        raise ValueError('Shared relaxed-structure reference/settings changed; rerun shared preflight')
    options = recovery.settings(config)
    if {k: v for k, v in options.items() if k != 'resolved_contract'} != policy['requested']:
        raise ValueError('Relaxed-structure parameters changed; rerun shared preflight')
    names = [lib['sublibrary_name'] for lib in config['libraries']]
    if names != policy['libraries']:
        raise ValueError('Sublibrary set/order changed; rerun shared preflight')
    for name in names:
        if recovery.source_stamp(barcode_path(config, name)) != policy['dedicated_inputs'][name]:
            raise ValueError(f'Dedicated input changed for {name}; rerun shared preflight')
    # JSON converts profile keys to strings; recreate the contract explicitly.
    design = compile_design(policy['design_contract']['rows'], config)
    design['input_fingerprint'] = policy['input_design_fingerprint']
    design['resolved_orientation'] = {int(k): v for k, v in policy['design_contract']['resolved_orientation'].items()}
    policy['design_contract'] = design
    return policy


def attach(worker, policy, config, reference, output, library):
    apply_worker(worker, policy['design_contract'])
    if policy['preflight']['status'].startswith('unavailable_'):
        worker.pop('recovery_contract', None)
    else:
        worker.update(recovery_contract=recovery.recovery_contract(reference, config, policy['design_contract']),
                      recovery_options=recovery.settings(config), recovery_policy=policy)
    atomic_json(Path(output) / f'{library}_recovery.json', policy)


def detail(result, original_endpoint, policy):
    return json.dumps(dict(length=result['recovery_length_state'],
                           losses=json.loads(result['recovery_possible_losses']),
                           tso_endpoints=json.loads(result['recovery_tso_endpoints']),
                           tso_shifts=[end-original_endpoint for end in json.loads(result['recovery_tso_endpoints'])],
                           ibar=result['recovery_ibar_state'], retained_edits=0,
                           policy=policy['fingerprint']), separators=(',', ':'))


def write_completion(match_path, library, policy):
    path = Path(match_path)
    atomic_json(path.parent / f'{library}_relaxed.complete.json',
                dict(complete=True, policy=policy['fingerprint'], match=recovery.source_stamp(path)))


def check_completion(match_path, library, policy):
    path = Path(match_path)
    marker = path.parent / f'{library}_relaxed.complete.json'
    if not marker.is_file():
        raise ValueError(f'Missing relaxed match completion for {library}; rerun Step 3/1c')
    saved = json.loads(marker.read_text())
    if not saved.get('complete') or saved.get('policy') != policy['fingerprint'] or saved.get('match') != recovery.source_stamp(path):
        raise ValueError(f'Stale or incomplete guide matches for {library}; rerun Step 3/1c before assignment')


def artifact_stamp(path):
    stamp = recovery.source_stamp(path)
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return dict(stamp, sha256=digest.hexdigest())


def assignment_identity(config, policy=None):
    keys = ('grna_structure', 'grna_matching', 'guide_design', 'guide_assignment', 'libraries')
    return recovery.fingerprint(dict(settings={key: config.get(key) for key in keys},
                                     reference=artifact_stamp(config['paths']['grna_reference']),
                                     policy=policy['fingerprint'] if policy else None))


def assignment_record(config, policy, tasks):
    return dict(identity=assignment_identity(config, policy), relaxed=active(config),
                inputs=[recovery.source_stamp(path) for paths, _, _, _ in tasks for _, path in paths])


def write_assignment_completion(directory, mode, record, products):
    if record is None:
        return
    for stamp in record['inputs']:
        if recovery.source_stamp(stamp['path']) != stamp:
            raise ValueError('Guide evidence changed during Step 4; rerun assignment')
    payload = dict(record, complete=True, outputs=[artifact_stamp(path) for path in products])
    atomic_json(Path(directory)/f'grna_{mode}.complete.json', payload)


def validate_assignment(config, directory, mode, policy=None):
    marker = Path(directory)/f'grna_{mode}.complete.json'
    if not marker.is_file():
        if active(config):
            raise ValueError(f'Missing relaxed {mode} assignment completion; rerun Step 4 before Step 5')
        return
    record = json.loads(marker.read_text())
    if not record.get('complete') or record.get('identity') != assignment_identity(config, policy):
        raise ValueError(f'Stale {mode} guide assignment settings/policy; rerun Step 4 before Step 5')
    for stamp in record['inputs'] + record['outputs']:
        current = artifact_stamp(stamp['path']) if 'sha256' in stamp else recovery.source_stamp(stamp['path'])
        if current != stamp:
            raise ValueError(f'Changed {mode} guide evidence/products; rerun Step 4 before Step 5')

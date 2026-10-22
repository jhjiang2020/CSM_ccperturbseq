#!/usr/bin/env python3
"""UGER configuration, dry runs and submission bookkeeping (no scientific imports).

Login-safe entry point: reads small YAML/scripts only; never opens assay data.
Run submit wrappers with ``timeout 30s`` on a login node. All numerical work is
performed by the submitted worker. Legacy configs may omit ``uger`` entirely.
"""
from __future__ import annotations

import argparse
import datetime as dt
import fcntl
import json
import math
import os
import re
import shlex
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import yaml


@dataclass(frozen=True)
class Resources:
    slots: int = 16
    memory_per_slot_gb: float = 8
    runtime_hours: float = 48
    max_concurrent: int | None = None
    queue: str | None = 'broad'
    parallel_environment: str = 'smp'
    binding: str = 'linear'
    os: str | None = None
    project: str | None = None
    extra_resources: dict | None = None

    @property
    def h_rt(self):
        seconds = round(self.runtime_hours * 3600)
        return f'{seconds // 3600:02d}:{seconds // 60 % 60:02d}:{seconds % 60:02d}'

    @property
    def h_vmem(self):
        return f'{self.memory_per_slot_gb:g}G'

    def qsub_args(self, *, array=False):
        args = ['-pe', self.parallel_environment, str(self.slots), '-binding',
                f'linear:{self.slots}' if self.binding == 'linear' else 'no',
                '-l', f'h_vmem={self.h_vmem},h_rt={self.h_rt}']
        for flag, value in (('-q', self.queue), ('-P', self.project)):
            if value:
                args += [flag, value]
        if self.os:
            args += ['-l', f'os={self.os}']
        for key, value in sorted((self.extra_resources or {}).items()):
            args += ['-l', f'{key}={value}']
        if array and self.max_concurrent is not None:
            args += ['-tc', str(self.max_concurrent)]
        return args


# Compatibility names retained for existing imports.
Step1aResources = Resources
@dataclass(frozen=True)
class Step3Resources(Resources):
    slots: int = 8
    memory_per_slot_gb: float = 4
    runtime_hours: float = 12


# The sole maintained resource-default table. Values preserve earlier requests;
# worker #$ lines remain standalone-qsub compatibility fallbacks.
# profile: (slots, GiB per slot, hours, array concurrency, required OS)
_DEFAULTS = {
    'step1a': (16, 8, 48, None, None),
    'step1b': (4, 16, 4, None, None),
    'step1c': (1, 4, 12, None, None),
    'step1c_aggregate': (4, 32, 8, None, None),
    'step2': (1, 4, 12, None, None),
    'step3': (8, 4, 12, None, None),
    'step4': (4, 32, 4, None, None),
    'step5': (1, 32, 4, None, None),
    'kallisto_reference': (8, 8, 12, None, 'RedHat8'),
    'kallisto_benchmark': (16, 8, 12, None, 'RedHat8'),
    'shared_kallisto_environment': (4, 8, 3, None, 'RedHat8'),
    'q1': (8, 4, 12, None, None),
    'gex_contamination': (8, 8, 48, 2, None),
    'gex_contamination_alignment': (8, 16, 12, 2, None),
    'gex_contamination_aggregate': (1, 4, 1, None, None),
    'c2_case_study': (4, 4, 4, None, None),
    'gex_demux_probe': (4, 4, 0.5, None, None),
    'gex_demux_run': (8, 4, 4, None, None),
    'gex_demux_qc': (2, 8, 1, None, None),
    'gex_umi_count': (2, 16, 12, 6, None),
    'gex_umi_control': (2, 16, 12, None, None),
    'gex_umi_compare': (2, 8, 2, None, None),
    'gex_carryover_scan': (2, 8, 12, 6, None),
    'gex_carryover_compare': (2, 8, 2, None, None),
    'gex_rrna_compare': (2, 8, 2, None, None),
    'rrna_reference': (4, 8, 4, None, None),
    'rrna_full': (16, 8, 168, 4, None),
    'rrna_aggregate': (2, 12, 8, None, None),
    'rrna_report': (1, 8, 1, None, None),
    'rrna_external_validate': (1, 8, 1, None, None),
    'rrna_external_preflight': (4, 4, 48, 4, None),
    'rrna_external_preaggregate': (1, 4, 2, None, None),
    'rrna_external_full': (16, 8, 168, None, None),
    'rrna_external_finalize': (1, 8, 4, None, None),
    'read_recovery_prepare': (1, 16, 4, None, None),
    'read_recovery_depth': (4, 8, 48, 4, None),
    'read_recovery_rawguide': (1, 8, 24, 4, None),
    'read_recovery_finalize': (1, 16, 4, None, None),
}
RESOURCE_FIELDS = set(Resources.__dataclass_fields__)
SCHEDULER_FIELDS = {'queue', 'parallel_environment', 'binding', 'os', 'project', 'extra_resources'}


def mapping(value, label):
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise ValueError(f'{label} must be a mapping')
    return value


def positive(value, label, maximum, *, integer=False):
    if isinstance(value, bool):
        raise ValueError(f'{label} must be a positive number')
    try:
        number = float(value)
    except (ValueError, TypeError) as exc:
        raise ValueError(f'{label} must be a positive number') from exc
    if not math.isfinite(number) or not 0 < number <= maximum or (integer and number != int(number)):
        raise ValueError(f'{label} must be positive and <= {maximum}' + (' (integer)' if integer else ''))
    return int(number) if integer else number


def token(value, label, *, nullable=True, commas=False):
    if value is None and nullable:
        return None
    pattern = r'[A-Za-z0-9_.*@/!:+,=-]+' if commas else r'[A-Za-z0-9_.*@/!:+=-]+'
    if not isinstance(value, str) or not re.fullmatch(pattern, value) or value.startswith('-'):
        raise ValueError(f'{label} must be a scheduler token' + (' or null' if nullable else ''))
    return value


def read_config(path):
    with Path(path).open() as handle:
        return mapping(yaml.safe_load(handle), 'config')


def resource_config(path):
    """QC files inherit uger from their base/sc config; local entries override it."""
    def visit(current, seen):
        current = Path(current).resolve()
        if current in seen:
            raise ValueError(f'Circular base/sc config reference: {current}')
        seen = seen | {current}
        cfg = read_config(current)
        base = cfg.get('base_config', cfg.get('sc_config'))
        result = {}
        if base:
            base_path = Path(base)
            if not base_path.is_absolute():
                base_path = current.parent / base_path
            # Q1's sc_config is an optional cross-reference; its standalone
            # worker historically tolerates an unavailable single-cell config.
            if 'base_config' in cfg or base_path.is_file():
                result = visit(base_path, seen)
        local = mapping(cfg.get('uger'), 'uger')
        unknown = set(local) - SCHEDULER_FIELDS - set(_DEFAULTS) - {'defaults'}
        if unknown:
            raise ValueError(f'Unknown uger keys: {sorted(unknown)}')
        for key, value in local.items():
            if key in _DEFAULTS or key in {'defaults', 'extra_resources'}:
                result[key] = {**mapping(result.get(key), f'uger.{key}'),
                               **mapping(value, f'uger.{key}')}
            else:
                result[key] = value
        return result
    return visit(path, set())


def load_resources(config_path, profile='step1a', *, fallback=None, uger=None):
    if profile not in _DEFAULTS:
        raise ValueError(f'Unknown resource profile: {profile}')
    slots, memory, hours, concurrent, os_name = _DEFAULTS[profile]
    defaults = Resources(slots, memory, hours, concurrent, os=os_name)
    raw_uger = resource_config(config_path) if uger is None else uger
    raw_defaults = mapping(raw_uger.get('defaults'), 'uger.defaults')
    raw_stage = mapping(raw_uger.get(profile), f'uger.{profile}')
    for label, values in (('uger.defaults', raw_defaults), (f'uger.{profile}', raw_stage)):
        if set(values) - RESOURCE_FIELDS:
            raise ValueError(f'Unknown {label} keys: {sorted(set(values) - RESOURCE_FIELDS)}')
    values = {**asdict(defaults), **(fallback or {}),
              **{k: v for k, v in raw_uger.items() if k in SCHEDULER_FIELDS},
              **raw_defaults, **raw_stage}
    label = f'uger.{profile}'
    values['slots'] = positive(values['slots'], f'{label}.slots', 256, integer=True)
    for key, maximum in (('memory_per_slot_gb', 4096), ('runtime_hours', 744)):
        values[key] = positive(values[key], f'{label}.{key}', maximum)
    if values['runtime_hours'] * 3600 < 1:
        raise ValueError(f'{label}.runtime_hours must represent at least one second')
    if values['max_concurrent'] is not None:
        values['max_concurrent'] = positive(values['max_concurrent'], f'{label}.max_concurrent', 100000, integer=True)
    for key in ('queue', 'os', 'project', 'parallel_environment'):
        values[key] = token(values[key], f'{label}.{key}', nullable=key != 'parallel_environment', commas=key == 'queue')
    if values['binding'] not in ('linear', 'none'):
        raise ValueError(f'{label}.binding must be linear or none')
    extra = mapping(values['extra_resources'], f'{label}.extra_resources')
    for key, value in extra.items():
        if not isinstance(key, str) or not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', key) or key in {'h_vmem', 'h_rt', 'os'}:
            raise ValueError(f'{label}.extra_resources: invalid/reserved resource {key!r}')
        token(str(value), f'{label}.extra_resources.{key}', nullable=False)
    values['extra_resources'] = {key: str(value) for key, value in extra.items()}
    return Resources(**values)


def load_step1a_resources(config_path):
    return load_resources(config_path, 'step1a')


def load_step3_resources(config_path):
    return Step3Resources(**asdict(load_resources(config_path, 'step3')))


def validate_all(config_path):
    uger = resource_config(config_path)
    return {key: load_resources(config_path, key, uger=uger) for key in _DEFAULTS}


def payload(resources):
    return {**asdict(resources), 'h_rt': resources.h_rt, 'h_vmem': resources.h_vmem,
            'total_memory_gb': resources.slots * resources.memory_per_slot_gb}


def split_options(args):
    """Only examine scheduler options before the worker, never its CLI arguments."""
    arity = {'-pe': 2, '-binding': 1, '-l': 1, '-q': 1, '-P': 1, '-tc': 1,
             '-N': 1, '-v': 1, '-o': 1, '-e': 1, '-j': 1, '-t': 1,
             '-hold_jid': 1, '-hold_jid_ad': 1, '-wd': 1, '-S': 1, '-V': 0,
             '-terse': 0, '-cwd': 0, '-b': 1, '-w': 1, '-m': 1, '-M': 1,
             '-p': 1, '-A': 1, '-a': 1, '-ac': 1, '-dc': 1, '-sc': 1,
             '-notify': 0, '-r': 1, '-shell': 1, '-sync': 1, '-now': 1,
             '-hard': 0, '-soft': 0, '-R': 1}
    i, options = 0, []
    while i < len(args) and args[i].startswith('-'):
        flag = args[i]
        if flag not in arity:
            raise ValueError(f'Unsupported qsub option {flag!r}; use documented options')
        n = arity[flag]
        if i + n >= len(args):
            raise ValueError(f'Missing value for {flag}')
        options.append((flag, args[i + 1:i + 1 + n]))
        i += 1 + n
    if i == len(args):
        raise ValueError('A worker script is required')
    return options, args[i:]


def resolve_command(config_path, profile, args):
    options, worker = split_options(args)
    fallback, kept = {}, []
    for flag, values in options:
        if flag == '-pe':
            fallback.update(parallel_environment=values[0], slots=int(values[1]))
        elif flag == '-binding':
            fallback['binding'] = 'none' if values[0] in ('none', 'no') else 'linear'
        elif flag in ('-q', '-P', '-tc'):
            fallback[{'-q': 'queue', '-P': 'project', '-tc': 'max_concurrent'}[flag]] = int(values[0]) if flag == '-tc' else values[0]
        elif flag == '-l':
            for expression in values[0].split(','):
                key, sep, value = expression.partition('=')
                if not sep:
                    raise ValueError(f'Resource needs a value: {expression}')
                if key == 'h_vmem':
                    match = re.fullmatch(r'([0-9]+(?:\.[0-9]+)?)([GMK])', value)
                    if not match:
                        raise ValueError('Legacy h_vmem must use G, M or K units')
                    fallback['memory_per_slot_gb'] = float(match[1]) / {'G': 1, 'M': 1024, 'K': 1024 ** 2}[match[2]]
                elif key == 'h_rt':
                    parts = [int(p) for p in value.split(':')]
                    if len(parts) != 3 or not 0 <= parts[1] < 60 or not 0 <= parts[2] < 60:
                        raise ValueError('Legacy h_rt must be HH:MM:SS')
                    fallback['runtime_hours'] = parts[0] + parts[1] / 60 + parts[2] / 3600
                elif key == 'os':
                    fallback['os'] = value
                else:
                    fallback.setdefault('extra_resources', {})[key] = value
        else:
            kept += [flag, *values]
    resources = load_resources(config_path, profile, fallback=fallback)
    return ['qsub', *resources.qsub_args(array=any(flag == '-t' for flag, _ in options)), *kept, *worker], resources


def submit(config_path, profile, args, *, dry_run=False, record=None):
    command, resources = resolve_command(config_path, profile, args)
    if dry_run:
        print('  [DRY RUN] ' + shlex.join(command))
        return None
    # Create/verify the record destination before dispatching any job.
    record = Path(record) if record else Path(read_config(config_path).get('paths', {}).get('log_dir', Path(config_path).parent / 'logs')) / 'uger_submissions.jsonl'
    record.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(record, os.O_APPEND | os.O_CREAT | os.O_WRONLY, 0o640)
    print(f'UGER {profile}: {resources.slots} slots, {resources.h_vmem}/slot, '
          f'{resources.slots * resources.memory_per_slot_gb:g}G total, {resources.h_rt}; '
          f'record: {record}', file=sys.stderr)
    try:
        actual = command if '-terse' in command else [command[0], '-terse', *command[1:]]
        output = subprocess.check_output(actual, text=True, timeout=25).strip()
        match = re.fullmatch(r'(\d+)(?:\.[\d,:-]+)?', output)
        if not match:
            raise RuntimeError(f'Cannot parse qsub job ID; scheduler output: {output!r}')
        entry = {'submitted_at': dt.datetime.now(dt.timezone.utc).isoformat(),
                 'config': str(Path(config_path).resolve()), 'profile': profile,
                 'job_id': match[1], 'qsub_output': output, 'command': actual,
                 'resources': payload(resources)}
        fcntl.flock(fd, fcntl.LOCK_EX)
        try:
            os.write(fd, (json.dumps(entry, sort_keys=True) + '\n').encode())
        finally:
            fcntl.flock(fd, fcntl.LOCK_UN)
    finally:
        os.close(fd)
    print(output)
    return match[1]


def pipeline(config_path, steps, dry_run):
    cfg = read_config(config_path)
    resources = validate_all(config_path)  # Fail before submitting a partial graph.
    selected = {'1', '2', '3', '4', '5'} if steps == 'all' else set(steps.split(','))
    if not selected or not selected <= {'1', '2', '3', '4', '5'}:
        raise ValueError('--step must be all or a comma-separated subset of 1,2,3,4,5')
    project = Path(cfg['paths']['project_dir'])
    logs = Path(cfg['paths'].get('log_dir', project / 'logs'))
    count = len(cfg.get('libraries', []))
    if count < 1:
        raise ValueError('At least one library is required')
    config_path = str(Path(config_path).resolve())
    # qsub -v uses commas as separators, even for correctly shell-quoted values.
    for value in (config_path, str(project)):
        if any(ch in value for ch in ",\n\r'"):
            raise ValueError('Config/project paths cannot contain commas/newlines (qsub -v) or apostrophes (legacy workers)')
    enabled = bool(cfg.get('gex_decontamination', {}).get('step1c', {}).get('enabled', False))
    env = f'CSM_CONFIG={config_path},CSM_PROJECT_DIR={project}'
    jobs = {}
    stages = [
        ('1', 'step1a', 'csm_step1a_gex', 'submit_step1a_gex_array.sh', True, []),
        ('1', 'step1b', 'csm_step1b_combine', 'submit_step1b_gex_combine.sh', False, ['step1a']),
        ('1c', 'step1c', 'csm_step1c_grna', 'submit_step1c_gex_grna_array.sh', True, ['step1a']),
        ('1c', 'step1c_aggregate', 'csm_step1c_aggregate', 'submit_step1c_gex_grna_aggregate.sh', False, ['step1c']),
        ('2', 'step2', 'csm_step2_bc', 'submit_step2_barcodes.sh', True, ['step1a']),
        ('3', 'step3', 'csm_step3_match', 'submit_step3_match.sh', True, ['step2']),
        ('4', 'step4', 'csm_step4_assign', 'submit_step4_assign.sh', False, ['step3', 'step1c_aggregate']),
        ('5', 'step5', 'csm_step5_anndata', 'submit_step5_anndata.sh', False, ['step1b', 'step4']),
    ]
    print(f'Config: {config_path}\nProject: {project}\nLibraries: {count}\nLogs: {logs}')
    for group, profile, name, script, array, dependencies in stages:
        if group == '1c':
            if '1' not in selected or not enabled:
                continue
        elif group not in selected:
            continue
        r = resources[profile]
        print(f'{profile}: {r.slots} slots x {r.h_vmem} = {r.slots * r.memory_per_slot_gb:g}G; {r.h_rt}', flush=True)
        extra_env = ',CSM_STEP1C_GEX_READY=true' if profile == 'step4' and 'step1c_aggregate' in jobs else ''
        options = ['-N', name, '-j', 'y', '-v', env + extra_env,
                   '-o', str(logs / (profile + '.$JOB_ID' + ('.$TASK_ID' if array else '') + '.log'))]
        if array:
            options += ['-t', f'1-{count}']
        holds = [jobs[key] for key in dependencies if key in jobs]
        if holds:
            options += ['-hold_jid', ','.join(holds)]
        job = submit(config_path, profile, [*options, str(project / 'submit' / script)],
                     dry_run=dry_run, record=logs / 'uger_submissions.jsonl')
        jobs[profile] = name if dry_run else job


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    parser.add_argument('--step', choices=tuple(_DEFAULTS), default='step1a')
    parser.add_argument('--format', choices=('json', 'tsv', 'args'), default='json')
    parser.add_argument('--pipeline', action='store_true')
    parser.add_argument('--steps', default='all')
    parser.add_argument('--submit-profile', choices=tuple(_DEFAULTS))
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--record')
    parser.add_argument('--validate', action='store_true')
    parser.add_argument('qsub_args', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.pipeline:
        pipeline(args.config, args.steps, args.dry_run)
        return
    if args.validate:
        print(json.dumps({key: payload(value) for key, value in validate_all(args.config).items()}, sort_keys=True, indent=2))
        return
    if args.submit_profile:
        validate_all(args.config)
        options = args.qsub_args[1:] if args.qsub_args[:1] == ['--'] else args.qsub_args
        submit(args.config, args.submit_profile, options, dry_run=args.dry_run, record=args.record)
        return
    r = load_resources(args.config, args.step)
    if args.format == 'tsv':
        print('\t'.join(map(str, (r.slots, f'{r.memory_per_slot_gb:g}', f'{r.runtime_hours:g}', r.h_rt))))
    elif args.format == 'args':
        print('\n'.join(r.qsub_args(array=True)))
    else:
        print(json.dumps(payload(r), sort_keys=True))


if __name__ == '__main__':
    main()

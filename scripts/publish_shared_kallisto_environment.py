#!/usr/bin/env python3
"""Inventory and seal a fresh shared kallisto tools prefix after runtime checks.

The production Parse prefix is inspected, never changed. The companion prefix
is versioned; reruns cannot overwrite a completed publication.
"""
from __future__ import annotations

import argparse
import datetime
import grp
import hashlib
import json
import os
import stat
import subprocess
import sys
import uuid
from pathlib import Path


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def command(*parts):
    try:
        return subprocess.check_output(list(map(str, parts)), text=True, stderr=subprocess.STDOUT).strip()
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f'{parts[0]} failed: {exc.output}') from exc


def write_json(path, payload):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n')
    os.chown(temporary, -1, path.parent.stat().st_gid)
    temporary.chmod(0o640)
    temporary.replace(path)


def publish_stable_alias(prefix, alias, production_prefix, record_dir):
    """Atomically point a stable public name at a completed immutable build."""
    prefix = prefix.resolve()
    alias = alias.parent.resolve() / alias.name
    if alias in (prefix, production_prefix.absolute()):
        raise ValueError('The stable alias must be separate from environment directories')
    manifest = json.loads((prefix / 'csm_shared_kallisto.complete.json').read_text())
    if manifest.get('complete') is not True or Path(manifest['prefix']).resolve() != prefix:
        raise ValueError('Only a completed tools publication may become the stable target')
    if alias.exists() and not alias.is_symlink():
        raise FileExistsError(f'Refusing to replace a directory or regular file: {alias}')
    temporary = alias.with_name(f'.{alias.name}.{uuid.uuid4().hex}.tmp')
    try:
        temporary.symlink_to(prefix.name if prefix.parent == alias.parent else prefix)
        os.lchown(temporary, -1, prefix.stat().st_gid)
        temporary.replace(alias)
    finally:
        temporary.unlink(missing_ok=True)
    write_json(record_dir / 'shared-kallisto-alias.complete.json', dict(
        complete=True, alias=str(alias), target=str(prefix), job_id=os.environ.get('JOB_ID'),
        published_at=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    print(f'Stable tools name: {alias} -> {prefix}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix', type=Path, required=True)
    parser.add_argument('--production-prefix', type=Path, required=True)
    parser.add_argument('--record-dir', type=Path, required=True)
    parser.add_argument('--project-dir', type=Path, required=True)
    parser.add_argument('--group', required=True)
    parser.add_argument('--verify-and-seal', action='store_true')
    parser.add_argument('--stable-alias', type=Path)
    parser.add_argument('--alias-only', action='store_true')
    args = parser.parse_args()
    prefix = args.prefix.resolve()
    if prefix == args.production_prefix.resolve():
        raise ValueError('Tools and production prefixes must be distinct')
    if args.alias_only:
        if args.stable_alias is None:
            parser.error('--alias-only requires --stable-alias')
        publish_stable_alias(prefix, args.stable_alias, args.production_prefix, args.record_dir)
        return
    complete = prefix / 'csm_shared_kallisto.complete.json'
    if complete.exists():
        raise FileExistsError('Refusing to overwrite a sealed environment')
    pending = prefix / 'csm_shared_kallisto.pending.json'
    if not args.verify_and_seal:
        import importlib.metadata as metadata
        for module in ('edlib', 'pandas', 'pysam', 'scipy', 'anndata', 'splitpipe', 'yaml'):
            __import__(module)
        kallisto = command(prefix / 'bin/kallisto', 'version')
        bustools = command(prefix / 'bin/bustools', 'version')
        # kb-python 0.29.5 intentionally exits 1 for `info` and bare --help.
        # Query metadata with its own interpreter and exercise the real CLI.
        kb = command(prefix / 'bin/python', '-c',
                     "from importlib.metadata import version; print(version('kb-python'))")
        command(prefix / 'bin/kb', 'ref', '--help')
        if '0.52.0' not in kallisto or '0.45.1' not in bustools or '0.29.5' not in kb:
            raise ValueError(f'Unexpected tools: {kallisto}; {bustools}; {kb}')
        artifacts = {name: {'sha256': digest(prefix / name), 'bytes': (prefix / name).stat().st_size}
                     for name in ('conda-explicit.txt', 'pip-freeze.txt', 'environment-source.yml')}
        payload = dict(schema_version=1, complete=False, publication='CSM ccPerturbSeq v2.0 tools',
            prefix=str(prefix), production_prefix=str(args.production_prefix.resolve()),
            production_manifest_sha256=digest(args.production_prefix / 'csm_shared_environment.complete.json'),
            owner_uid=os.getuid(), group=args.group, job_id=os.environ.get('JOB_ID'),
            created_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
            python=sys.version, kallisto=kallisto, bustools=bustools, kb=kb,
            production_packages={name: metadata.version(name) for name in ('pandas','scipy','anndata','pysam','edlib','pyyaml')},
            artifacts=artifacts,
            binaries={name: digest(prefix / 'bin' / name) for name in ('kallisto','bustools','kb')},
            tests={'log': str(args.record_dir / 'runtime-tests.log'),
                   'sha256': digest(args.record_dir / 'runtime-tests.log'), 'passed': True},
            runtime_source_sha256={name: digest(args.project_dir / 'scripts' / name) for name in (
                'run_kallisto_gex.py','kallisto_umi.py','kallisto_count_parallel.py','gex_alignment.py')})
        write_json(pending, payload)
        print('Inventory created; awaiting permission audit')
        return
    payload = json.loads(pending.read_text())
    group_id = grp.getgrnam(args.group).gr_gid
    files = 0
    for directory, names, filenames in os.walk(prefix, followlinks=False):
        for path in [Path(directory), *(Path(directory) / name for name in filenames),
                     *(Path(directory) / name for name in names if (Path(directory) / name).is_symlink())]:
            info = path.lstat()
            if info.st_uid != os.getuid() or info.st_gid != group_id:
                raise ValueError(f'Unexpected owner/group: {path}')
            if not stat.S_ISLNK(info.st_mode):
                if info.st_mode & 0o027:
                    raise ValueError(f'Unexpected write/other permission: {path}')
                if not info.st_mode & stat.S_IRGRP:
                    raise ValueError(f'Group cannot read: {path}')
                if stat.S_ISDIR(info.st_mode) and not info.st_mode & stat.S_IXGRP:
                    raise ValueError(f'Group cannot traverse: {path}')
            files += 1
    for name, record in payload['artifacts'].items():
        if digest(prefix / name) != record['sha256']:
            raise ValueError(f'Inventory changed: {name}')
    if digest(args.production_prefix / 'csm_shared_environment.complete.json') != payload['production_manifest_sha256']:
        raise ValueError('Production publication changed')
    payload.update(complete=True, permission_audit_entries=files,
                   published_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
    write_json(complete, payload)
    pending.unlink()
    write_json(args.record_dir / 'shared-v2.complete.json', payload)
    if args.stable_alias is not None:
        publish_stable_alias(prefix, args.stable_alias, args.production_prefix, args.record_dir)
    print(f'Published validated shared tools: {prefix}')


if __name__ == '__main__':
    main()

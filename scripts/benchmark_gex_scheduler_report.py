#!/usr/bin/env python3
"""Login-safe, stdlib-only summary of compact qacct/submission metadata.

Run with ``timeout 30s python scripts/benchmark_gex_scheduler_report.py RUN``.
Use --require-complete after jobs end; known unstarted cancellations are read
from RUN/monitor_state.json and listed separately from executed tasks.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('run', type=Path)
parser.add_argument('--require-complete', action='store_true')
args = parser.parse_args()
run = args.run
state_path = run / 'monitor_state.json'
state = json.loads(state_path.read_text()) if state_path.exists() else {}
canceled = {job: {str(t) for t in tasks} for job, tasks in state.get('canceled_before_start', {}).items()}
records = {}
for line in (run / 'logs/uger_submissions.jsonl').read_text().splitlines():
    entry = json.loads(line)
    records[entry['job_id']] = entry
for path in (run / 'logs/uger_submissions.records').glob('*.json'):
    entry = json.loads(path.read_text())
    records[entry['job_id']] = entry
rows, coverage = [], []
for job, entry in sorted(records.items()):
    command = entry['command']
    name = command[command.index('-N') + 1] if '-N' in command else ''
    tasks = command[command.index('-t') + 1] if '-t' in command else 'undefined'
    if tasks == 'undefined':
        expected = {'undefined'}
    else:
        match = re.fullmatch(r'(\d+)(?:-(\d+))?(?::(\d+))?', tasks)
        if not match:
            raise ValueError(tasks)
        first, last, step = match.groups()
        expected = {str(i) for i in range(int(first), int(last or first) + 1, int(step or 1))}
    source = run / f'logs/accounting/{job}.qacct.txt'
    blocks, current = [], {}
    if source.exists():
        for line in source.read_text().splitlines():
            if line.startswith('====='):
                if current:
                    blocks.append(current)
                current = {}
            elif len(line.split(None, 1)) == 2:
                key, value = line.split(None, 1)
                current[key] = value.strip()
        if current:
            blocks.append(current)
    observed = set()
    for block in blocks:
        if block.get('jobnumber') != job:
            raise ValueError(f'Wrong job accounting: {source}')
        task = block['taskid']
        if task in observed:
            raise ValueError(f'Duplicate attempt requires explicit audit: {job}/{task}')
        observed.add(task)
        slots = int(block['slots'])
        wall = float(block['wallclock'])
        cpu = float(block['cpu'])
        failed = int(block['failed'].split()[0])
        exit_status = int(block['exit_status'])
        rows.append(dict(job_id=job, task_id=task, job_name=name,
            status='success' if failed == 0 and exit_status == 0 else 'failed',
            failed=failed, exit_status=exit_status, hostname=block['hostname'],
            slots=slots, requested_memory_gb=entry['resources']['total_memory_gb'],
            qacct_wall_seconds=wall, qacct_ru_wallclock_seconds=float(block['ru_wallclock']),
            cpu_hours=cpu/3600, reserved_core_hours=wall*slots/3600,
            maxvmem=block['maxvmem'], raw_ru_maxrss=block['ru_maxrss'],
            start_time=block['start_time'], end_time=block['end_time'],
            accounting_source=source.relative_to(run).as_posix()))
    if observed - expected:
        raise ValueError(f'Unexpected tasks {observed - expected}: {job}')
    missing = expected - observed
    known_cancellations = missing & canceled.get(job, set())
    unresolved = missing - known_cancellations
    coverage.append(dict(job_id=job,job_name=name,submitted_tasks=tasks,
        expected_tasks=len(expected),accounted_tasks=len(observed),
        canceled_before_start=','.join(sorted(known_cancellations,key=lambda x: int(x) if x.isdigit() else -1)),
        missing_task_ids=','.join(sorted(unresolved,key=lambda x: int(x) if x.isdigit() else -1))))
output = run / 'comparison'
output.mkdir(exist_ok=True)
for name, values in [('scheduler_accounting',rows),('scheduler_coverage',coverage)]:
    with (output / f'{name}.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(values[0]));writer.writeheader();writer.writerows(values)
summary=dict(recorded_jobs=len(records),accounted_tasks=len(rows),success_tasks=sum(x['status']=='success' for x in rows),failed_tasks=sum(x['status']=='failed' for x in rows),scheduler_task_hours=sum(x['qacct_wall_seconds'] for x in rows)/3600,cpu_hours=sum(x['cpu_hours'] for x in rows),reserved_core_hours=sum(x['reserved_core_hours'] for x in rows),missing_accounting=[x for x in coverage if x['missing_task_ids']],exporter_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(output/'scheduler_accounting_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
if args.require_complete and summary['missing_accounting']:
    raise SystemExit('Missing accounting for executed or unresolved tasks; collect again after completion')

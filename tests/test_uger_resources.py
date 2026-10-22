import tempfile
import subprocess
import unittest
from pathlib import Path

import yaml

from scripts.uger_resources import load_step1a_resources, load_step3_resources


ROOT = Path(__file__).resolve().parents[1]


class UgerResourceConfigTest(unittest.TestCase):
    def write_config(self, root: Path, payload: dict) -> Path:
        path = root / "config.yaml"
        path.write_text(yaml.safe_dump(payload))
        return path

    def test_defaults_match_direct_step1a_wrapper(self):
        with tempfile.TemporaryDirectory() as temporary:
            resources = load_step1a_resources(
                self.write_config(Path(temporary), {})
            )
        self.assertEqual(resources.slots, 16)
        self.assertEqual(resources.memory_per_slot_gb, 8)
        self.assertEqual(resources.runtime_hours, 48)
        self.assertEqual(resources.h_rt, "48:00:00")
        self.assertEqual(resources.h_vmem, "8G")

    def test_custom_values_are_resolved(self):
        with tempfile.TemporaryDirectory() as temporary:
            resources = load_step1a_resources(
                self.write_config(
                    Path(temporary),
                    {
                        "uger": {
                            "step1a": {
                                "slots": 24,
                                "memory_per_slot_gb": 10,
                                "runtime_hours": 72,
                            }
                        }
                    },
                )
            )
        self.assertEqual(
            (resources.slots, resources.h_vmem, resources.h_rt),
            (24, "10G", "72:00:00"),
        )

    def test_step3_defaults_and_custom_values(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            defaults = load_step3_resources(self.write_config(root, {}))
            custom = load_step3_resources(self.write_config(root, {
                'uger': {'step3': {'slots': 16, 'memory_per_slot_gb': 3, 'runtime_hours': 48}}
            }))
        self.assertEqual((defaults.slots, defaults.h_vmem, defaults.h_rt), (8, '4G', '12:00:00'))
        self.assertEqual((custom.slots, custom.h_vmem, custom.h_rt), (16, '3G', '48:00:00'))

    def test_invalid_step3_resources_are_rejected(self):
        for raw in ({'slots': 0}, {'slots': '8;touch'}, {'runtime_hours': False}, []):
            with self.subTest(raw=raw), tempfile.TemporaryDirectory() as temporary:
                path = self.write_config(Path(temporary), {'uger': {'step3': raw}})
                with self.assertRaises(ValueError):
                    load_step3_resources(path)

    def test_step3_dry_run_preserves_dependency_and_sets_all_resources(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = self.write_config(root, {
                'paths': {'project_dir': str(ROOT), 'log_dir': str(root/'logs')},
                'libraries': [{'sublibrary_name': 'L1'}],
                'uger': {'step3': {'slots': 16, 'memory_per_slot_gb': 3, 'runtime_hours': 48}},
            })
            result = subprocess.run(['bash',str(ROOT/'submit/submit_pipeline.sh'),
                '--config',str(path),'--step','2,3','--dry-run'],
                check=True,capture_output=True,text=True)
        command = next(line for line in result.stdout.splitlines() if '[DRY RUN]' in line and 'submit_step3_match.sh' in line)
        self.assertIn('-pe smp 16 -binding linear:16', command)
        self.assertIn('-l h_vmem=3G,h_rt=48:00:00', command)
        self.assertIn('-hold_jid csm_step2_bc', command)

    def test_invalid_values_are_rejected_before_submission(self):
        invalid = (
            {"slots": 0},
            {"slots": "16;touch"},
            {"memory_per_slot_gb": True},
            {"runtime_hours": 0},
        )
        for step1a in invalid:
            with self.subTest(step1a=step1a), tempfile.TemporaryDirectory() as temporary:
                path = self.write_config(
                    Path(temporary), {"uger": {"step1a": step1a}}
                )
                with self.assertRaises(ValueError):
                    load_step1a_resources(path)

    def test_pipeline_dry_run_emits_matching_custom_binding_and_limits(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = self.write_config(
                root,
                {
                    "paths": {
                        "project_dir": str(ROOT),
                        "log_dir": str(root / "logs"),
                    },
                    "libraries": [{"sublibrary_name": "L1"}],
                    "gex_decontamination": {
                        "step1c": {"enabled": False},
                        "combine_with_grna_arm": {"enabled": False},
                    },
                    "uger": {
                        "step1a": {
                            "slots": 24,
                            "memory_per_slot_gb": 10,
                            "runtime_hours": 72,
                        }
                    },
                },
            )
            completed = subprocess.run(
                [
                    "bash",
                    str(ROOT / "submit" / "submit_pipeline.sh"),
                    "--config",
                    str(config),
                    "--step",
                    "1",
                    "--dry-run",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
        self.assertIn("-pe smp 24 -binding linear:24", completed.stdout)
        self.assertIn("-l h_vmem=10G,h_rt=72:00:00", completed.stdout)



class UnifiedUgerResourcesTest(unittest.TestCase):
    def config(self, root, payload):
        path = root / 'run.yaml'
        path.write_text(yaml.safe_dump(payload))
        return path

    def test_old_configs_keep_every_production_stage_default(self):
        from scripts.uger_resources import validate_all
        with tempfile.TemporaryDirectory() as tmp:
            resolved = validate_all(self.config(Path(tmp), {}))
        expected = {'step1a': (16, 8, 48), 'step1b': (4, 16, 4),
                    'step1c': (1, 4, 12), 'step1c_aggregate': (4, 32, 8),
                    'step2': (1, 4, 12), 'step3': (8, 4, 12),
                    'step4': (4, 32, 4), 'step5': (1, 32, 4)}
        for name, values in expected.items():
            r = resolved[name]
            self.assertEqual((r.slots, r.memory_per_slot_gb, r.runtime_hours), values)

    def test_scheduler_defaults_stage_precedence_and_array_limit(self):
        from scripts.uger_resources import load_resources, resolve_command
        with tempfile.TemporaryDirectory() as tmp:
            path = self.config(Path(tmp), {'uger': {'queue': 'other.q', 'os': 'RedHat8',
                'parallel_environment': 'threaded', 'project': 'lab',
                'defaults': {'memory_per_slot_gb': 2.5, 'max_concurrent': 3},
                'step1b': {'slots': 6, 'memory_per_slot_gb': 7, 'runtime_hours': 1.5}}})
            cmd, r = resolve_command(path, 'step1b', ['-t', '1-9', 'worker.sh'])
            single, _ = resolve_command(path, 'step1b', ['worker.sh'])
        self.assertEqual((r.slots, r.h_vmem, r.h_rt), (6, '7G', '01:30:00'))
        self.assertIn('linear:6', cmd)
        self.assertIn('os=RedHat8', cmd)
        self.assertEqual(cmd[cmd.index('-q') + 1], 'other.q')
        self.assertEqual(cmd[cmd.index('-pe') + 1], 'threaded')
        self.assertIn('-P', cmd)
        self.assertIn('-tc', cmd)
        self.assertNotIn('-tc', single)

    def test_qc_inherits_resources_and_can_override_one_field(self):
        from scripts.uger_resources import load_resources
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.config(root, {'uger': {'queue': 'inherited', 'gex_contamination': {'slots': 12}}})
            qc = root / 'qc.yaml'
            qc.write_text(yaml.safe_dump({'base_config': 'run.yaml',
                'uger': {'gex_contamination': {'memory_per_slot_gb': 6}}}))
            r = load_resources(qc, 'gex_contamination')
        self.assertEqual((r.slots, r.memory_per_slot_gb, r.queue), (12, 6, 'inherited'))

    def test_invalid_scheduler_and_resource_fields_fail(self):
        from scripts.uger_resources import validate_all
        bad = [{'step5': {'slots': True}}, {'defaults': {'slots': 1.2}},
               {'defaults': {'runtime_hours': float('nan')}}, {'queue': 'x;touch'},
               {'defaults': {'memory_gb': 64}}, {'step1b': {'max_concurrent': 0}},
               {'defaults': {'binding': 'linear:256'}}, {'step1z': {}},
               {'extra_resources': {'h_vmem': '128G'}}, {'os': ['RedHat8']}]
        for value in bad:
            with self.subTest(value=value), tempfile.TemporaryDirectory() as tmp:
                path = self.config(Path(tmp), {'uger': value})
                with self.assertRaises(ValueError):
                    validate_all(path)

    def test_standalone_q1_tolerates_unavailable_optional_sc_config(self):
        from scripts.uger_resources import load_resources
        with tempfile.TemporaryDirectory() as tmp:
            path = self.config(Path(tmp), {'sc_config': 'unavailable.yaml'})
            self.assertEqual(load_resources(path, 'q1').slots, 8)

    def test_circular_base_config_rejected(self):
        from scripts.uger_resources import validate_all
        with tempfile.TemporaryDirectory() as tmp:
            path = self.config(Path(tmp), {'base_config': 'run.yaml'})
            with self.assertRaisesRegex(ValueError, 'Circular'):
                validate_all(path)

    def test_legacy_qc_cli_values_remain_fallbacks_and_worker_flags_untouched(self):
        from scripts.uger_resources import resolve_command
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            args = ['-pe', 'smp', '3', '-l', 'h_vmem=512M,h_rt=00:30:00',
                    '-t', '1-8', '-tc', '2', 'worker.sh', '--threads', '999', '-l', 'worker-option']
            cmd, r = resolve_command(self.config(root, {}), 'gex_contamination', args)
            self.assertEqual((r.slots, r.h_vmem, r.h_rt), (3, '0.5G', '00:30:00'))
            cmd, r = resolve_command(self.config(root, {'uger': {'gex_contamination': {'slots': 7}}}), 'gex_contamination', args)
        self.assertEqual(r.slots, 7)
        self.assertIn('linear:7', cmd)
        self.assertEqual(cmd[-4:], ['--threads', '999', '-l', 'worker-option'])

    def test_actual_submissions_hold_exact_ids_and_record_resolved_resources(self):
        import json
        from unittest.mock import patch
        from scripts.uger_resources import pipeline
        with tempfile.TemporaryDirectory(prefix='csm quoted ') as tmp:
            root = Path(tmp)
            path = self.config(root, {'paths': {'project_dir': str(root), 'log_dir': str(root / 'logs')},
                'libraries': [{'sublibrary_name': 'L'}],
                'gex_decontamination': {'step1c': {'enabled': True}},
                'uger': {'step1b': {'slots': 6}}})
            with patch('scripts.uger_resources.subprocess.check_output', side_effect=[f'{100+i}.1-2:1\n' for i in range(8)]) as qsub:
                pipeline(path, 'all', False)
            commands = [call.args[0] for call in qsub.call_args_list]
            rows = [json.loads(row) for row in (root / 'logs/uger_submissions.jsonl').read_text().splitlines()]
        self.assertEqual(len(rows), 8)
        self.assertEqual(rows[1]['resources']['slots'], 6)
        self.assertIn('100', commands[1])
        self.assertIn('105,103', commands[6])
        self.assertIn('101,106', commands[7])
        self.assertIn('CSM_STEP1C_GEX_READY=true', ','.join(commands[6]))
        self.assertEqual(rows[1]['resources']['total_memory_gb'], 96)

    def test_invalid_later_stage_prevents_partial_submission(self):
        from unittest.mock import patch
        from scripts.uger_resources import pipeline
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = self.config(root, {'paths': {'project_dir': str(root)},
                'libraries': [{'sublibrary_name': 'L'}], 'uger': {'step5': {'slots': 0}}})
            with patch('scripts.uger_resources.subprocess.check_output') as qsub:
                with self.assertRaises(ValueError):
                    pipeline(path, 'all', False)
                qsub.assert_not_called()

    def test_resource_defaults_match_standalone_worker_directives(self):
        import re
        from scripts.uger_resources import load_resources
        workers = {'step1a': 'step1a_gex_array', 'step1b': 'step1b_gex_combine',
            'step1c': 'step1c_gex_grna_array', 'step1c_aggregate': 'step1c_gex_grna_aggregate',
            'step2': 'step2_barcodes', 'step3': 'step3_match', 'step4': 'step4_assign',
            'step5': 'step5_anndata', 'kallisto_reference': 'kallisto_reference',
            'kallisto_benchmark': 'kallisto_benchmark', 'q1': 'q1_guideseq',
            'shared_kallisto_environment': 'shared_kallisto_environment'}
        with tempfile.TemporaryDirectory() as tmp:
            path = self.config(Path(tmp), {})
            for profile, worker in workers.items():
                with self.subTest(profile=profile):
                    r = load_resources(path, profile)
                    text = (ROOT / f'submit/submit_{worker}.sh').read_text()
                    self.assertEqual(int(re.search(r'^#\$ -pe smp (\d+)', text, re.M)[1]), r.slots)
                    self.assertEqual(float(re.search(r'^#\$ -l h_vmem=(\d+)G', text, re.M)[1]), r.memory_per_slot_gb)
                    hours, minutes, seconds = map(int, re.search(r'^#\$ -l h_rt=([\d:]+)', text, re.M)[1].split(':'))
                    self.assertEqual(hours + minutes/60 + seconds/3600, r.runtime_hours)


if __name__ == "__main__":
    unittest.main()

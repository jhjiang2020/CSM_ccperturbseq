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


if __name__ == "__main__":
    unittest.main()

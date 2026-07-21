import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import yaml

from scripts.write_fastq_input_identity import build_payload


ROOT = Path(__file__).resolve().parents[1]


class ReleasePackagingTest(unittest.TestCase):
    def test_development_and_audit_tools_are_not_packaged(self):
        tracked = subprocess.check_output(
            ["git", "ls-files"], cwd=ROOT, text=True
        ).splitlines()
        forbidden = re.compile(
            r"(^|/)(v14_audit_|submit_v14_|validate_v14_|.*_candidate\.py$|"
            r"benchmark_gex_grna_filter\.py$|debug_env\.sh$)"
        )
        exposed = [path for path in tracked if forbidden.search(path)]
        self.assertEqual(exposed, [])

    def test_every_literal_wrapper_runtime_asset_is_git_tracked(self):
        tracked = set(subprocess.check_output(
            ["git", "ls-files"], cwd=ROOT, text=True
        ).splitlines())
        patterns = (
            (re.compile(r'\$\{SCRIPTS_DIR\}/([A-Za-z0-9_.-]+)'), "scripts"),
            (re.compile(r'\$\{PROJECT_DIR\}/scripts/([A-Za-z0-9_.-]+)'), "scripts"),
            (re.compile(r'\$\{PROJECT_DIR\}/patches/([A-Za-z0-9_.-]+)'), "patches"),
        )
        for wrapper in sorted((ROOT / "submit").glob("*.sh")):
            text = wrapper.read_text()
            for pattern, directory in patterns:
                for name in pattern.findall(text):
                    relative = f"{directory}/{name}"
                    self.assertTrue((ROOT / relative).is_file(), (wrapper.name, relative))
                    self.assertIn(relative, tracked, (wrapper.name, relative))

    def test_step1a_runtime_assets_are_git_tracked(self):
        assets = (
            "patches/parsebio_splitpipe_1.6.4_gex_decontamination.patch",
            "scripts/splitpipe_gex_wrapper.py",
            "scripts/gex_grna_preprocess_hook.py",
            "scripts/gex_grna_filter.py",
            "scripts/gex_grna_bam_filter.py",
            "scripts/write_fastq_input_identity.py",
        )
        tracked = set(subprocess.check_output(
            ["git", "ls-files"], cwd=ROOT, text=True
        ).splitlines())
        for relative in assets:
            self.assertTrue((ROOT / relative).is_file(), relative)
            self.assertIn(relative, tracked, relative)

    def test_existing_alignment_step1a_supplies_mol_statistics(self):
        text = (ROOT / "submit/submit_step1a_gex_array.sh").read_text()
        self.assertIn('SOURCE_PRE_ALIGN_STATS="${SOURCE_SUBLIB}/process/pre_align_stats.csv"', text)
        self.assertIn('--source-pre-align-stats "${SOURCE_PRE_ALIGN_STATS}"', text)
        self.assertLess(
            text.index('--source-pre-align-stats "${SOURCE_PRE_ALIGN_STATS}"'),
            text.index("--mode mol"),
        )

    def test_c2_case_study_release_assets_are_tracked(self):
        tracked = set(subprocess.check_output(
            ["git", "ls-files"], cwd=ROOT, text=True
        ).splitlines())
        assets = (
            "docs/C2_GEX_GRNA_MISALIGNMENT_CASE_STUDY.md",
            "scripts/c2_gex_grna_case_study.py",
            "submit/submit_c2_gex_grna_case_study.sh",
            "docs/assets/c2_gex_grna_misalignment/read_transformation.png",
            "docs/assets/c2_gex_grna_misalignment/read_transformation.svg",
            "docs/assets/c2_gex_grna_misalignment/klf2_locus.png",
            "docs/assets/c2_gex_grna_misalignment/klf2_locus.svg",
            "docs/assets/c2_gex_grna_misalignment/nfkb1_locus.png",
            "docs/assets/c2_gex_grna_misalignment/nfkb1_locus.svg",
        )
        for relative in assets:
            self.assertTrue((ROOT / relative).is_file(), relative)
            self.assertIn(relative, tracked, relative)
        evidence_root = ROOT / "docs/assets/c2_gex_grna_misalignment/data"
        evidence_manifest = yaml.safe_load(
            (evidence_root / "extraction.complete.json").read_text()
        )
        evidence_files = {
            "extraction.complete.json", *evidence_manifest["output_files"]
        }
        for name in evidence_files:
            relative = f"docs/assets/c2_gex_grna_misalignment/data/{name}"
            self.assertTrue((ROOT / relative).is_file(), relative)
            self.assertIn(relative, tracked, relative)

    def test_c2_wrapper_accepts_absolute_production_environment(self):
        text = (ROOT / "submit/submit_c2_gex_grna_case_study.sh").read_text()
        self.assertNotIn('CONDA_ENV}" != "splitpipe', text)
        self.assertIn('conda activate "${CONDA_ENV}"', text)
        self.assertIn("import pandas, pysam, splitpipe, yaml", text)

    def test_quickstart_publishes_cross_user_environment(self):
        text = (ROOT / "README.md").read_text()
        shared = "/medpop/rgupta/Jiahao/conda/envs/splitpipe-share"
        self.assertEqual(text.count(f'conda_prefix: "{shared}"'), 1)
        self.assertEqual(text.count(f'conda_env: "{shared}"'), 1)
        self.assertIn("installation-free shared environment", text)
        self.assertNotIn("/home/unix/jjiahao/", text)

    def test_every_multislot_production_wrapper_has_matching_binding(self):
        for path in sorted((ROOT / "submit").glob("*.sh")):
            text = path.read_text()
            match = re.search(r"^#\$ -pe smp (\d+)\s*$", text, re.MULTILINE)
            if not match:
                continue
            if int(match.group(1)) <= 1:
                continue
            binding = re.search(r"^#\$ -binding linear:(\d+)\s*$", text, re.MULTILINE)
            self.assertIsNotNone(binding, path.name)
            self.assertEqual(match.group(1), binding.group(1), path.name)
            self.assertIn("NSLOTS", text, path.name)

    def test_production_wrappers_have_no_personal_conda_fallback(self):
        for path in sorted((ROOT / "submit").glob("*.sh")):
            self.assertNotIn("$HOME/my-conda", path.read_text(), path.name)

    def test_input_identity_preserves_numeric_wrapper_order(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = {}
            for read in (1, 2):
                for lane in (2, 10):
                    path = root / f"sample_L{lane}_R{read}_001.fastq.gz"
                    path.touch()
                    paths[(read, lane)] = path
            payload = build_payload(
                "L1",
                [paths[(1, 2)], paths[(1, 10)]],
                [paths[(2, 2)], paths[(2, 10)]],
            )
            self.assertEqual(payload["pair_count"], 2)
            self.assertTrue(payload["ordered_identity_sha256"])
            self.assertIn("L2", payload["pairs"][0]["r1"]["path"])
            self.assertIn("L10", payload["pairs"][1]["r1"]["path"])

    def test_relocated_project_step1_dry_run(self):
        with tempfile.TemporaryDirectory() as temporary:
            relocated = Path(temporary) / "relocated_ccperturbseq"
            (relocated / "submit").mkdir(parents=True)
            (relocated / "scripts").mkdir()
            shutil.copy(ROOT / "submit/submit_pipeline.sh", relocated / "submit")
            shutil.copy(ROOT / "scripts/uger_resources.py", relocated / "scripts")
            config = {
                "environment": {
                    "conda_prefix": "/configured/conda",
                    "conda_env": "splitpipe",
                },
                "uger": {"step1a": {
                    "slots": 16, "memory_per_slot_gb": 8, "runtime_hours": 48,
                }},
                "paths": {
                    "project_dir": str(relocated),
                    "log_dir": str(relocated / "logs"),
                },
                "libraries": [{
                    "sublibrary_name": "L1", "gex_prefix": "L1",
                    "grna_prefix": "G1", "sample_name": "sample",
                }],
                "gex_decontamination": {
                    "step1c": {"enabled": False},
                    "combine_with_grna_arm": {"enabled": False},
                },
            }
            config_path = relocated / "config.yaml"
            config_path.write_text(yaml.safe_dump(config))
            result = subprocess.run(
                [
                    "bash", str(relocated / "submit/submit_pipeline.sh"),
                    "--config", str(config_path), "--step", "1", "--dry-run",
                ],
                text=True,
                capture_output=True,
                check=True,
            )
            self.assertIn(str(relocated), result.stdout)
            self.assertNotIn("/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq", result.stdout)
            self.assertIn("-pe smp 16 -binding linear:16", result.stdout)


if __name__ == "__main__":
    unittest.main()

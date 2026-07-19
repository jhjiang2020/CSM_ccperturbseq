import importlib.util
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]


def load_script(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


step2 = load_script("step2_adapters", "02_extract_grna_barcodes.py")
step3 = load_script("step3_adapters", "03_match_grna_ibar.py")
step1c = load_script("step1c_runner", "01c_process_gex_grna.py")
step1c_validator = load_script(
    "step1c_validator", "validate_step1c_manifests.py"
)


def write_step1a_evidence_manifest(root, library="C1", rows=1):
    evidence = root / "valid_positive_barcode_evidence.tsv.gz"
    with gzip.open(evidence, "wt") as handle:
        handle.write("read_id\tbc_wells\tbc_bci\tumi\tr1_seq\tr1_qual\n")
        for index in range(rows):
            handle.write(
                f"read{index}\tA1_A1_A1\t1_1_1\tAAAAAAAAAA\tACGT\tIIII\n"
            )
    manifest = root / "rejected_pairs.complete.json"
    manifest.write_text(json.dumps({
        "schema_version": 4,
        "status": "complete",
        "complete": True,
        "library": library,
        "barcode_evidence": str(evidence),
        "barcode_evidence_sha256": hashlib.sha256(evidence.read_bytes()).hexdigest(),
        "barcode_evidence_rows": rows,
        "statistics": {
            "csm_gex_filter_positive_valid_barcode_reads": rows,
        },
    }))
    return manifest, evidence


class Step1cAdapterTests(unittest.TestCase):
    def test_default_output_layout_is_unchanged(self):
        output = "/project/output"
        self.assertEqual(step2.resolve_output_root(output), Path(output))
        self.assertEqual(step3.resolve_output_root(output), Path(output))
        self.assertEqual(
            step2.resolve_output_root(output, "gex_grna_identity"),
            Path(output) / "gex_grna_identity",
        )
        self.assertEqual(
            step3.resolve_output_root(output, "/scratch/identity"),
            Path("/scratch/identity"),
        )

    def test_runner_runs_step3_directly_on_step1a_barcode_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config_path = root / "config.yaml"
            config_path.write_text(yaml.safe_dump({
                "paths": {"output_dir": str(root / "output")},
            }))
            manifest, evidence = write_step1a_evidence_manifest(root)
            commands = step1c.build_commands(
                ROOT / "scripts", "/env/python", config_path, "C1", manifest,
                "gex_grna_identity", 8,
            )
            step3_cmd, namespace, barcode_input, effective_workers = commands
            self.assertIn("--barcode-input", step3_cmd)
            barcode_arg = Path(step3_cmd[step3_cmd.index("--barcode-input") + 1])
            self.assertEqual(barcode_arg, evidence.resolve())
            self.assertEqual(barcode_input, evidence.resolve())
            self.assertEqual(namespace, root / "output/gex_grna_identity")
            self.assertEqual(effective_workers, 1)

    def test_runner_dry_run_needs_no_splitpipe_import(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config_path = root / "config.yaml"
            config_path.write_text(yaml.safe_dump({
                "paths": {"output_dir": str(root / "output")},
                "libraries": [{"sublibrary_name": "C1"}],
                "barcode_extraction": {"n_processes": 3},
            }))
            manifest, _ = write_step1a_evidence_manifest(root)
            result = subprocess.run(
                [
                    sys.executable, str(ROOT / "scripts/01c_process_gex_grna.py"),
                    "--config", str(config_path), "--library", "C1",
                    "--input-manifest", str(manifest), "--dry-run",
                ],
                check=True, text=True, capture_output=True,
            )
            self.assertIn("03_match_grna_ibar.py", result.stdout)
            self.assertNotIn("02_extract_grna_barcodes.py", result.stdout)
            self.assertIn("--n_processes 1", result.stdout)

    def test_step1c_rejects_stale_or_recalibrated_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest, evidence = write_step1a_evidence_manifest(root)
            payload = json.loads(manifest.read_text())
            payload["schema_version"] = 3
            manifest.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError, "schema 4"):
                step1c.load_barcode_evidence(manifest, "C1")
            payload["schema_version"] = 4
            manifest.write_text(json.dumps(payload))
            with gzip.open(evidence, "at") as handle:
                handle.write("extra\tA1_A1_A1\t1_1_1\tAAAAAAAAAA\tACGT\tIIII\n")
            with self.assertRaisesRegex(ValueError, "SHA-256"):
                step1c.load_barcode_evidence(manifest, "C1")

    def test_strict_step1c_product_preserves_raw_step3_matches(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            match = root / "C2_grna_matches.tsv.gz"
            fields = [
                "read_id", "bc_wells", "bc_bci", "umi", "grna_name",
                "grna_seq", "ibar_seq", "target_gene", "grna_dist",
                "ibar_dist", "ibar_seq_corrected", "ibar_decoupled",
            ]
            rows = [
                ["keep", "A1_A1_A1", "1_1_1", "AAAAAAAAAA", "g1", "ACGT", "AAAAAAA", "G", "1", "1", "AAAAAAA", "False"],
                ["guide2", "A1_A1_A1", "1_1_1", "CCCCCCCCCC", "g1", "ACGT", "AAAAAAA", "G", "2", "0", "AAAAAAA", "False"],
                ["noibar", "A1_A1_A1", "1_1_1", "GGGGGGGGGG", "g1", "ACGT", "", "G", "0", "-1", "", "False"],
                ["decoupled", "A1_A1_A1", "1_1_1", "TTTTTTTTTT", "g1", "ACGT", "CCCCCCC", "G", "0", "1", "CCCCCCC", "True"],
            ]
            with gzip.open(match, "wt") as handle:
                handle.write("\t".join(fields) + "\n")
                for row in rows:
                    handle.write("\t".join(row) + "\n")
            config = {
                "guide_assignment": {"dual_guide_counting_mode": True},
                "gex_decontamination": {"step1c": {"strict_matching": {
                    "enabled": True,
                    "max_guide_edit_distance": 1,
                    "max_ibar_edit_distance": 1,
                    "require_resolved_ibar_in_dual_mode": True,
                }}},
            }
            strict, raw, stats = step1c.install_strict_match_product(match, config, "C2")
            with gzip.open(strict, "rt") as handle:
                strict_lines = handle.readlines()
            with gzip.open(raw, "rt") as handle:
                raw_lines = handle.readlines()
            self.assertEqual(len(strict_lines), 2)
            self.assertEqual(len(raw_lines), 5)
            self.assertEqual(stats["retained_matches"], 1)
            self.assertEqual(stats["dropped_guide_distance"], 1)
            self.assertEqual(stats["dropped_unresolved_ibar"], 1)
            self.assertEqual(stats["dropped_ibar_decoupled"], 1)

    def test_dual_short_r1_emits_ibar_warning(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            evidence = root / "valid_positive_barcode_evidence.tsv.gz"
            with gzip.open(evidence, "wt") as handle:
                handle.write("read_id\tbc_wells\tbc_bci\tumi\tr1_seq\tr1_qual\n")
                handle.write(
                    "r1\tA1_A1_A1\t1_1_1\tAAAAAAAAAA\t"
                    + "A" * 54 + "\t" + "I" * 54 + "\n"
                )
            manifest = root / "rejected_pairs.complete.json"
            manifest.write_text(json.dumps({"barcode_evidence": str(evidence)}))
            config = {
                "grna_structure": {
                    "fixed_seq_1": "A" * 30,
                    "fixed_seq_2": "C" * 19,
                    "fixed_seq_3": "G" * 19,
                    "grna_lengths": [20, 21],
                    "ibar_length": 7,
                },
                "guide_assignment": {"dual_guide_counting_mode": True},
            }
            message = step1c.warn_if_ibar_unobservable(config, manifest)
            self.assertIn("54 < 76 bp", message)

    def test_step1c_validator_rejects_stale_config(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output = root / "output"
            source = root / "step1a.json"
            source.write_text("{}")
            barcode = root / "barcodes.tsv.gz"
            match = root / "matches.tsv.gz"
            barcode.touch()
            match.touch()
            config_path = root / "config.yaml"
            config = {
                "paths": {"output_dir": str(output)},
                "libraries": [{"sublibrary_name": "C1"}],
                "gex_decontamination": {
                    "step1c": {
                        "enabled": True,
                        "output_namespace": "gex_grna_identity",
                        "strict_matching": {"enabled": False},
                    }
                },
            }
            config_path.write_text(yaml.safe_dump(config))
            config_hash = hashlib.sha256(config_path.read_bytes()).hexdigest()
            source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
            completion = output / "gex_grna_identity/manifests/C1.step1c.complete.json"
            completion.parent.mkdir(parents=True)
            completion.write_text(json.dumps({
                "schema_version": 2,
                "complete": True,
                "status": "complete",
                "source_library": "C1",
                "config_sha256": config_hash,
                "step1a_manifest": str(source),
                "step1a_manifest_sha256": source_hash,
                "barcode_output": str(barcode),
                "barcode_output_sha256": hashlib.sha256(barcode.read_bytes()).hexdigest(),
                "match_output": str(match),
                "match_output_sha256": hashlib.sha256(match.read_bytes()).hexdigest(),
            }))
            manifests = step1c_validator.validate(config_path)
            self.assertEqual(len(manifests), 1)
            summary_path = step1c_validator.write_audit_summary(
                config_path, manifests
            )
            summary = yaml.safe_load(summary_path.read_text())
            self.assertIn("C1", summary["libraries"])
            config["analysis_name"] = "changed"
            config_path.write_text(yaml.safe_dump(config))
            with self.assertRaisesRegex(ValueError, "stale config"):
                step1c_validator.validate(config_path)


if __name__ == "__main__":
    unittest.main()

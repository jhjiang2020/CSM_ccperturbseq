import tempfile
import unittest
import subprocess
from pathlib import Path

from scripts.qc_fastq_resolver import resolve_all_library_pairs
from scripts.qc_gex_grna_contamination import load_merged_config
from tests.release_test_helpers import (
    library_aggregate_case,
    resolver_generality,
)


ROOT = Path(__file__).resolve().parents[1]


def touch_pair(root, stem):
    r1 = root / f"{stem}_R1_001.fastq.gz"
    r2 = root / f"{stem}_R2_001.fastq.gz"
    r1.touch()
    r2.touch()
    return r1, r2


class QcFastqResolverTest(unittest.TestCase):
    def test_qc_launcher_can_gate_alignment_and_aggregate_without_rescan(self):
        launcher = ROOT / "submit/submit_qc_gex_grna_contamination_pipeline.sh"
        alignment = subprocess.run(
            ["bash", str(launcher), "--mode", "alignment", "--dry-run"],
            text=True,
            capture_output=True,
            check=True,
        ).stdout
        self.assertIn("Alignment command:", alignment)
        self.assertNotIn("Array command:", alignment)
        aggregate = subprocess.run(
            ["bash", str(launcher), "--mode", "aggregate", "--dry-run"],
            text=True,
            capture_output=True,
            check=True,
        ).stdout
        self.assertIn("Aggregate: qsub submit_", aggregate)
        self.assertNotIn("Array command:", aggregate)

    def test_production_qc_aggregates_one_four_and_twelve_libraries(self):
        base = load_merged_config(str(ROOT / "qc_config_gex_grna_contamination.yaml"))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cases = {
                count: library_aggregate_case(
                    root=root,
                    base_cfg=base,
                    library_count=count,
                )
                for count in (1, 4, 12)
            }
        for count, result in cases.items():
            self.assertEqual(result["library_count"], count)
            self.assertEqual(result["aggregate_summary_rows"], count)
            self.assertEqual(result["complete_aggregate_rows"], count)

    def test_stage5_executes_all_one_and_four_library_resolver_gates(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = resolver_generality(Path(temporary))
        self.assertEqual(result["one_library"]["pair_count"], 3)
        self.assertEqual(result["four_libraries"]["library_count"], 4)
        self.assertTrue(result["four_libraries"]["overlapping_prefixes_isolated"])
        self.assertTrue(result["four_libraries"]["explicit_pairs_supported"])
        self.assertTrue(result["duplicate_names_rejected"])
        self.assertTrue(result["shared_fastq_pairs_rejected"])

    def test_qc_submitters_have_no_user_specific_environment_fallback(self):
        submitters = (
            "submit_q1_guideseq.sh",
            "submit_qc_gex_grna_contamination.sh",
            "submit_qc_gex_grna_contamination_aggregate.sh",
            "submit_qc_gex_grna_alignment_array.sh",
            "submit_qc_gex_grna_contamination_pipeline.sh",
        )
        forbidden = ("/home/unix/jjiahao", "$HOME/my-conda", "#$ -o /medpop")
        for name in submitters:
            text = (ROOT / "submit" / name).read_text()
            for value in forbidden:
                self.assertNotIn(value, text, (name, value))

    def test_overlapping_prefixes_are_token_bounded(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            c1 = touch_pair(root, "lane_C1")
            c10 = touch_pair(root, "lane_C10")
            cfg = {
                "paths": {"fastq_dir": str(root)},
                "libraries": [
                    {"sublibrary_name": "C1", "gex_prefix": "C1"},
                    {"sublibrary_name": "C10", "gex_prefix": "C10"},
                ],
            }
            resolved = resolve_all_library_pairs(cfg, str(root / "config.yaml"))
            self.assertEqual(
                resolved["C1"],
                [tuple(str(path.resolve()) for path in c1)],
            )
            self.assertEqual(
                resolved["C10"],
                [tuple(str(path.resolve()) for path in c10)],
            )

    def test_explicit_pairs_and_multiple_lanes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            lane1 = touch_pair(root, "unrelated_lane1")
            lane2 = touch_pair(root, "unrelated_lane2")
            cfg = {
                "paths": {"fastq_dir": str(root)},
                "libraries": [
                    {
                        "sublibrary_name": "LIB",
                        "gex_fastq_pairs": [
                            {"r1": str(lane1[0]), "r2": str(lane1[1])},
                            [str(lane2[0]), str(lane2[1])],
                        ],
                    }
                ],
            }
            resolved = resolve_all_library_pairs(cfg, str(root / "config.yaml"))
            self.assertEqual(len(resolved["LIB"]), 2)

    def test_duplicate_names_and_shared_pairs_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            pair = touch_pair(root, "lane_SHARED")
            duplicate_names = {
                "paths": {"fastq_dir": str(root)},
                "libraries": [
                    {"sublibrary_name": "X", "gex_prefix": "SHARED"},
                    {"sublibrary_name": "X", "gex_prefix": "SHARED"},
                ],
            }
            with self.assertRaisesRegex(ValueError, "duplicate sublibrary_name"):
                resolve_all_library_pairs(
                    duplicate_names,
                    str(root / "config.yaml"),
                )
            shared = {
                "paths": {"fastq_dir": str(root)},
                "libraries": [
                    {
                        "sublibrary_name": "X",
                        "gex_fastq_pairs": [[str(pair[0]), str(pair[1])]],
                    },
                    {
                        "sublibrary_name": "Y",
                        "gex_fastq_pairs": [[str(pair[0]), str(pair[1])]],
                    },
                ],
            }
            with self.assertRaisesRegex(ValueError, "assigned to both"):
                resolve_all_library_pairs(shared, str(root / "config.yaml"))


if __name__ == "__main__":
    unittest.main()

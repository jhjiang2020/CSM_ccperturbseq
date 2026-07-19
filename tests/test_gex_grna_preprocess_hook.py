import gzip
import importlib.util
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import yaml


ROOT = Path(__file__).resolve().parents[1]
HOOK_PATH = ROOT / "scripts" / "gex_grna_preprocess_hook.py"
TSO = "AAGCAGTGGTATCAACGCAGAGTGAATGGG"
GUIDE = "ACGTACGTACGTACGTACGT"


def load_hook(name):
    spec = importlib.util.spec_from_file_location(name, HOOK_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fastq_count(path):
    with gzip.open(path, "rt") as handle:
        lines = sum(1 for _ in handle)
    if lines % 4:
        raise AssertionError("malformed FASTQ")
    return lines // 4


def fastq_ids(path):
    result = []
    with gzip.open(path, "rt") as handle:
        while True:
            header = handle.readline().strip()
            if not header:
                break
            result.append(header.split()[0].lstrip("@"))
            handle.readline(); handle.readline(); handle.readline()
    return result


class PreprocessHookTest(unittest.TestCase):
    def make_config(self, root, step1c=True, existing_alignment=False):
        reference = root / "reference.csv"
        reference.write_text("grna_seq,grna_name\n" f"{GUIDE},anonymous\n")
        config = {
            "paths": {"grna_reference": str(reference)},
            "grna_structure": {
                "fixed_seq_1": TSO,
                "grna_lengths": [20],
                "tso_max_edit_dist": 3,
            },
            "gex_decontamination": {
                "enabled": True,
                "tso_max_edit_dist": 1,
                "max_tso_start": 2,
                "tso_matching": {"backend": "myers"},
                "filter_matching": {
                    "backend": "packed_inline",
                    "guide_distance": "hamming",
                    "max_guide_mismatches": 1,
                    "allow_guide_indels": False,
                    "full_guide_required": True,
                    "remove_reference_supported_ambiguous_hits": True,
                },
                "io": {"write_rejected_raw_pairs": "auto"},
                "step1c": {"enabled": step1c},
                "existing_alignment": {"enabled": existing_alignment},
            },
        }
        config_path = root / "config.yaml"
        config_path.write_text(yaml.safe_dump(config))
        return config_path

    def environment(self, config, output, run_id):
        return {
            "CSM_GEX_DECONTAM_CONFIG": str(config),
            "CSM_GEX_DECONTAM_LIBRARY": "C1",
            "CSM_GEX_DECONTAM_OUTPUT_DIR": str(output),
            "CSM_GEX_DECONTAM_RUN_ID": run_id,
        }

    def test_positive_pairs_and_statistics_are_finalized_atomically(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = self.make_config(root, step1c=True)
            output = root / "filter" / "C1"
            with patch.dict(os.environ, self.environment(config, output, "run1"), clear=False):
                hook = load_hook("csm_hook_test_positive")
                hook.preload_from_environment()
                hook.start_worker(fastq_index=0, chunk=0)
                positive = hook.classify_and_record(
                    header1="@cluster1 1:N:0:IDX",
                    seq1=TSO + GUIDE + "GTTT",
                    quality1="I" * (len(TSO) + len(GUIDE) + 4),
                    header2="@cluster1 2:N:0:IDX",
                    seq2="A" * 58,
                    quality2="I" * 58,
                    barcode_valid=True,
                )
                hook.record_barcode_evidence(
                    header1="@cluster1 1:N:0:IDX",
                    bc_wells="A1_A1_A1",
                    bc_bci="1_1_1",
                    umi="AAAAAAAAAA",
                    seq1=TSO + GUIDE + "GTTT",
                    quality1="I" * (len(TSO) + len(GUIDE) + 4),
                )
                negative_decision = hook.classify_and_record_decision(
                    header1="@cluster2 1:N:0:IDX",
                    seq1=TSO + "T" * 24,
                    quality1="I" * (len(TSO) + 24),
                    header2="@cluster2 2:N:0:IDX",
                    seq2="C" * 58,
                    quality2="I" * 58,
                    barcode_valid=True,
                )
                trim_end, trim_route = hook.resolve_survivor_trim(
                    decision=negative_decision,
                    seq1=TSO + "T" * 24,
                    parse_tso_sequence=TSO[13:],
                    parse_search_end=32,
                )
                self.assertEqual(trim_end, len(TSO))
                hook.record_survivor_trim(trim_route)
                negative = bool(negative_decision.filter_positive)
                worker_stats = hook.close_worker()
                manifest_path = Path(hook.finalize_from_environment())

            self.assertTrue(positive)
            self.assertFalse(negative)
            self.assertEqual(worker_stats["csm_gex_filter_positive_reads"], 1)
            manifest = json.loads(manifest_path.read_text())
            self.assertTrue(manifest["complete"])
            self.assertEqual(manifest["library"], "C1")
            self.assertEqual(manifest["rejected_pair_count"], 1)
            self.assertEqual(manifest["schema_version"], 4)
            self.assertEqual(manifest["barcode_evidence_rows"], 1)
            self.assertTrue(manifest["barcode_evidence_sha256"])
            self.assertEqual(fastq_count(Path(manifest["rejected_r1"])), 1)
            self.assertEqual(fastq_count(Path(manifest["rejected_r2"])), 1)
            self.assertEqual(
                manifest["statistics"]["csm_gex_filter_positive_valid_barcode_reads"],
                1,
            )
            self.assertEqual(
                manifest["rates"]["filter_positive_fraction_of_examined_reads"],
                0.5,
            )
            self.assertEqual(
                manifest["rates"]["removed_from_star_fraction_of_valid_barcode_reads"],
                0.5,
            )
            self.assertFalse(list(output.glob(".run1.*")))

    def test_auto_rejected_output_is_off_when_step1c_is_disabled(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = self.make_config(root, step1c=False)
            output = root / "filter" / "C1"
            with patch.dict(os.environ, self.environment(config, output, "run2"), clear=False):
                hook = load_hook("csm_hook_test_no_pairs")
                hook.preload_from_environment()
                hook.start_worker(fastq_index=0, chunk=None)
                self.assertTrue(hook.classify_and_record(
                    header1="@cluster1 1:N:0:IDX",
                    seq1=TSO + GUIDE,
                    quality1="I" * (len(TSO) + len(GUIDE)),
                    header2="@cluster1 2:N:0:IDX",
                    seq2="A" * 58,
                    quality2="I" * 58,
                    barcode_valid=False,
                ))
                hook.close_worker()
                manifest = json.loads(Path(hook.finalize_from_environment()).read_text())
            self.assertFalse(manifest["write_rejected_raw_pairs"])
            self.assertIsNone(manifest["rejected_r1"])
            self.assertEqual(
                manifest["statistics"]["csm_gex_filter_positive_invalid_barcode_reads"],
                1,
            )

    def test_step1c_evidence_does_not_require_raw_rejected_fastqs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = self.make_config(root, step1c=True)
            payload = yaml.safe_load(config.read_text())
            payload["gex_decontamination"]["io"]["write_rejected_raw_pairs"] = False
            config.write_text(yaml.safe_dump(payload))
            output = root / "filter" / "C1"
            with patch.dict(
                os.environ,
                self.environment(config, output, "evidence-only"),
                clear=False,
            ):
                hook = load_hook("csm_hook_test_evidence_only")
                hook.preload_from_environment()
                hook.start_worker(fastq_index=0, chunk=0)
                self.assertTrue(hook.classify_and_record(
                    header1="@cluster1 1:N:0:IDX",
                    seq1=TSO + GUIDE,
                    quality1="I" * (len(TSO) + len(GUIDE)),
                    header2="@cluster1 2:N:0:IDX",
                    seq2="A" * 58,
                    quality2="I" * 58,
                    barcode_valid=True,
                ))
                hook.record_barcode_evidence(
                    header1="@cluster1 1:N:0:IDX",
                    bc_wells="A1_A1_A1",
                    bc_bci="1_1_1",
                    umi="AAAAAAAAAA",
                    seq1=TSO + GUIDE,
                    quality1="I" * (len(TSO) + len(GUIDE)),
                )
                hook.close_worker()
                manifest = json.loads(
                    Path(hook.finalize_from_environment()).read_text()
                )
            self.assertFalse(manifest["write_rejected_raw_pairs"])
            self.assertIsNone(manifest["rejected_r1"])
            self.assertEqual(manifest["barcode_evidence_rows"], 1)
            self.assertTrue(Path(manifest["barcode_evidence"]).is_file())

    def test_existing_alignment_auto_retains_positive_ids(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = self.make_config(
                root, step1c=False, existing_alignment=True
            )
            output = root / "filter" / "C1"
            with patch.dict(os.environ, self.environment(config, output, "legacy"), clear=False):
                hook = load_hook("csm_hook_test_legacy_pairs")
                hook.preload_from_environment()
                hook.start_worker(fastq_index=0, chunk=None)
                hook.classify_and_record(
                    header1="@legacy1 1:N:0:IDX",
                    seq1=TSO + GUIDE,
                    quality1="I" * (len(TSO) + len(GUIDE)),
                    header2="@legacy1 2:N:0:IDX",
                    seq2="A" * 58,
                    quality2="I" * 58,
                    barcode_valid=True,
                )
                hook.record_barcode_evidence(
                    header1="@legacy1 1:N:0:IDX",
                    bc_wells="A1_A1_A1",
                    bc_bci="1_1_1",
                    umi="AAAAAAAAAA",
                    seq1=TSO + GUIDE,
                    quality1="I" * (len(TSO) + len(GUIDE)),
                )
                hook.close_worker()
                manifest = json.loads(
                    Path(hook.finalize_from_environment()).read_text()
                )
            self.assertTrue(manifest["write_rejected_raw_pairs"])
            self.assertEqual(fastq_count(Path(manifest["rejected_r1"])), 1)

    def test_survivor_trim_routes_full_partial_and_no_tso(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = self.make_config(root, step1c=False)
            output = root / "filter" / "C1"
            with patch.dict(
                os.environ,
                self.environment(config, output, "trim-routes"),
                clear=False,
            ):
                hook = load_hook("csm_hook_test_trim_routes")
                hook.preload_from_environment()
                hook.start_worker(fastq_index=0, chunk=0)

                # An error in Parse's terminal 17-mer: the fuzzy full-TSO
                # detector succeeds while Parse's exact suffix search fails.
                replacement = "A" if TSO[20] != "A" else "C"
                fuzzy_tso = TSO[:20] + replacement + TSO[21:]
                fuzzy_seq = fuzzy_tso + "T" * 24
                fuzzy_decision = hook.classify_and_record_decision(
                    header1="@fuzzy-full 1:N:0:IDX",
                    seq1=fuzzy_seq,
                    quality1="I" * len(fuzzy_seq),
                    header2="@fuzzy-full 2:N:0:IDX",
                    seq2="A" * 58,
                    quality2="I" * 58,
                    barcode_valid=True,
                )
                self.assertTrue(fuzzy_decision.tso_found)
                self.assertFalse(fuzzy_decision.filter_positive)
                self.assertNotIn(TSO[13:], fuzzy_seq[:32])
                trim_end, route = hook.resolve_survivor_trim(
                    decision=fuzzy_decision,
                    seq1=fuzzy_seq,
                    parse_tso_sequence=TSO[13:],
                    parse_search_end=32,
                )
                self.assertEqual((trim_end, route), (len(TSO), "full_tso"))
                hook.record_survivor_trim(route)

                # A partial TSO retaining only Parse's exact 17-mer remains a
                # normal survivor and uses the original Parse fallback.
                partial_seq = TSO[13:] + "C" * 30
                partial_decision = hook.classify_and_record_decision(
                    header1="@partial 1:N:0:IDX",
                    seq1=partial_seq,
                    quality1="I" * len(partial_seq),
                    header2="@partial 2:N:0:IDX",
                    seq2="C" * 58,
                    quality2="I" * 58,
                    barcode_valid=True,
                )
                self.assertFalse(partial_decision.tso_found)
                trim_end, route = hook.resolve_survivor_trim(
                    decision=partial_decision,
                    seq1=partial_seq,
                    parse_tso_sequence=TSO[13:],
                    parse_search_end=32,
                )
                self.assertEqual((trim_end, route), (17, "parse_exact_17mer"))
                hook.record_survivor_trim(route)

                no_tso_seq = "C" * 54
                no_tso_decision = hook.classify_and_record_decision(
                    header1="@no-tso 1:N:0:IDX",
                    seq1=no_tso_seq,
                    quality1="I" * len(no_tso_seq),
                    header2="@no-tso 2:N:0:IDX",
                    seq2="G" * 58,
                    quality2="I" * 58,
                    barcode_valid=True,
                )
                trim_end, route = hook.resolve_survivor_trim(
                    decision=no_tso_decision,
                    seq1=no_tso_seq,
                    parse_tso_sequence=TSO[13:],
                    parse_search_end=32,
                )
                self.assertEqual((trim_end, route), (0, "none"))
                hook.record_survivor_trim(route)

                worker_stats = hook.close_worker()
                manifest = json.loads(
                    Path(hook.finalize_from_environment()).read_text()
                )

            self.assertEqual(
                worker_stats["csm_gex_survivor_full_tso_trimmed_reads"], 1
            )
            self.assertEqual(
                worker_stats["csm_gex_survivor_parse_17mer_trimmed_reads"], 1
            )
            self.assertEqual(worker_stats["csm_gex_survivor_untrimmed_reads"], 1)
            self.assertEqual(manifest["schema_version"], 4)
            self.assertEqual(
                manifest["survivor_tso_trim_policy"],
                "full_tso_endpoint_then_parse_exact_17mer",
            )
            self.assertAlmostEqual(
                manifest["rates"]["full_tso_trim_fraction_of_star_survivors"],
                1 / 3,
            )
            self.assertAlmostEqual(
                manifest["rates"]["parse_17mer_trim_fraction_of_star_survivors"],
                1 / 3,
            )

    def test_fused_compact_route_preserves_trim_and_raw_artifacts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = self.make_config(root, step1c=True)
            output = root / "filter" / "C1"
            with patch.dict(
                os.environ,
                self.environment(config, output, "fused-route"),
                clear=False,
            ):
                hook = load_hook("csm_hook_test_fused_route")
                hook.preload_from_environment()
                hook.start_worker(
                    fastq_index=0,
                    chunk=0,
                    parse_tso_sequence=TSO[13:],
                    parse_search_end=32,
                )
                positive_seq = TSO + GUIDE + "GTT"
                decision, trim_end = hook.process_compact(
                    seq1=positive_seq, barcode_valid=True
                )
                self.assertTrue(decision.filter_positive)
                self.assertEqual(trim_end, 0)
                hook.record_rejected_pair(
                    header1=b"@positive 1:N:0:IDX\n",
                    seq1=positive_seq.encode(),
                    quality1=("I" * len(positive_seq)).encode(),
                    header2=b"@positive 2:N:0:IDX\n",
                    seq2=b"A" * 58,
                    quality2=b"I" * 58,
                )
                hook.record_barcode_evidence(
                    header1=b"@positive 1:N:0:IDX\n",
                    bc_wells="A1_A1_A1",
                    bc_bci="1_1_1",
                    umi="AAAAAAAAAA",
                    seq1=positive_seq.encode(),
                    quality1=("I" * len(positive_seq)).encode(),
                )

                partial_seq = TSO[13:] + "C" * 30
                partial, trim_end = hook.process_compact(
                    seq1=partial_seq, barcode_valid=True
                )
                self.assertFalse(partial.tso_found)
                self.assertEqual(trim_end, 17)
                worker_stats = hook.close_worker()
                manifest = json.loads(
                    Path(hook.finalize_from_environment()).read_text()
                )

            self.assertEqual(worker_stats["csm_gex_filter_positive_reads"], 1)
            self.assertEqual(
                worker_stats["csm_gex_survivor_parse_17mer_trimmed_reads"], 1
            )
            self.assertEqual(fastq_count(Path(manifest["rejected_r1"])), 1)

    def test_manifest_reuse_detects_removal_policy_changes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = self.make_config(root, step1c=True)
            output = root / "filter" / "C1"
            environment = self.environment(config, output, "fingerprint")
            with patch.dict(os.environ, environment, clear=False):
                hook = load_hook("csm_hook_test_fingerprint_original")
                hook.preload_from_environment()
                hook.start_worker(fastq_index=0, chunk=None)
                hook.close_worker()
                hook.finalize_from_environment()
                self.assertTrue(hook.existing_filter_manifest_is_compatible())

                manifest_path = output / "rejected_pairs.complete.json"
                old_schema = json.loads(manifest_path.read_text())
                old_schema["schema_version"] = 2
                manifest_path.write_text(json.dumps(old_schema))
                stale = load_hook("csm_hook_test_fingerprint_old_schema")
                self.assertFalse(stale.existing_filter_manifest_is_compatible())

                # Restore the current manifest for the evidence and filter-
                # setting compatibility checks below.
                old_schema["schema_version"] = 3
                manifest_path.write_text(json.dumps(old_schema))

            manifest_payload = json.loads(
                (output / "rejected_pairs.complete.json").read_text()
            )
            Path(manifest_payload["rejected_r1"]).unlink()
            with patch.dict(os.environ, environment, clear=False):
                missing_evidence = load_hook(
                    "csm_hook_test_fingerprint_missing_evidence"
                )
                self.assertFalse(
                    missing_evidence.existing_filter_manifest_is_compatible()
                )

            # Restore the evidence so the policy-fingerprint assertion below
            # tests the intended mismatch rather than the missing file.
            with gzip.open(manifest_payload["rejected_r1"], "wt"):
                pass

            payload = yaml.safe_load(config.read_text())
            payload["gex_decontamination"]["max_tso_start"] = 1
            config.write_text(yaml.safe_dump(payload))
            with patch.dict(os.environ, environment, clear=False):
                changed = load_hook("csm_hook_test_fingerprint_changed")
                self.assertFalse(changed.existing_filter_manifest_is_compatible())

    def test_worker_artifacts_use_numeric_fastq_chunk_order(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = self.make_config(root, step1c=True)
            output = root / "filter" / "C1"
            with patch.dict(
                os.environ, self.environment(config, output, "numeric"), clear=False
            ):
                hook = load_hook("csm_hook_test_numeric_order")
                hook.preload_from_environment()
                for fastq_index, chunk, read_id in ((10, 11, "late"), (2, 3, "early")):
                    sequence = TSO + GUIDE
                    hook.start_worker(fastq_index=fastq_index, chunk=chunk)
                    hook.classify_and_record(
                        header1=f"@{read_id} 1:N:0:IDX",
                        seq1=sequence,
                        quality1="I" * len(sequence),
                        header2=f"@{read_id} 2:N:0:IDX",
                        seq2="A" * 58,
                        quality2="I" * 58,
                        barcode_valid=True,
                    )
                    hook.record_barcode_evidence(
                        header1=f"@{read_id} 1:N:0:IDX",
                        bc_wells="A1_A1_A1",
                        bc_bci="1_1_1",
                        umi="AAAAAAAAAA",
                        seq1=sequence,
                        quality1="I" * len(sequence),
                    )
                    hook.close_worker()
                manifest = json.loads(Path(hook.finalize_from_environment()).read_text())
            self.assertEqual(fastq_ids(Path(manifest["rejected_r1"])), ["early", "late"])
            self.assertIn(".f2.c3.", manifest["worker_parts"][0])
            self.assertIn(".f10.c11.", manifest["worker_parts"][1])


if __name__ == "__main__":
    unittest.main()

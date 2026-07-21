import importlib.util
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]


def load_script(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


validator = load_script(
    "gex_filter_manifest_validator", "validate_gex_filter_manifests.py"
)
filter_module = load_script("gex_filter_for_validation", "gex_grna_filter.py")


def write_stats(path, values):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "statistic,value\n"
        + "".join(f"{key},{value}\n" for key, value in values.items())
    )


class FilterManifestValidationTest(unittest.TestCase):
    def test_historical_bam_contract_rehashes_source_and_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output = root / "output_v1_4_bamtest"
            source = root / "production/gex_C1/process/source.bam"
            source.parent.mkdir(parents=True)
            source.write_bytes(b"production-source")
            reference = root / "reference.csv"
            reference.write_text(
                "grna_seq,grna_name\nACGTACGTACGTACGTACGT,g1\n"
            )
            config = {
                "paths": {
                    "output_dir": str(output),
                    "grna_reference": str(reference),
                },
                "libraries": [{"sublibrary_name": "C1"}],
                "grna_structure": {
                    "fixed_seq_1": "AAGCAGTGGTATCAACGCAGAGTGAATGGG",
                    "grna_lengths": [20],
                },
                "gex_decontamination": {
                    "enabled": True,
                    "tso_matching": {"backend": "myers"},
                    "filter_matching": {
                        "backend": "packed_inline",
                        "guide_distance": "hamming",
                        "max_guide_mismatches": 1,
                        "allow_guide_indels": False,
                        "full_guide_required": True,
                    },
                    "existing_alignment": {
                        "enabled": True,
                        "sandbox_output_dir_template": "{output_dir}/gex_{library}",
                    },
                },
            }
            config_path = root / "config.yaml"
            config_path.write_text(yaml.safe_dump(config))
            classifier = filter_module.load_filter_from_config(config_path)
            reference_hash = classifier.index.stats.reference_sha256
            settings_hash = filter_module.filter_settings_sha256(
                classifier.settings
            )
            filter_dir = output / "gex_grna_filter/C1"
            filter_dir.mkdir(parents=True)
            (filter_dir / "rejected_pairs.complete.json").write_text(
                json.dumps(
                    {
                        "schema_version": 4,
                        "complete": True,
                        "status": "complete",
                        "library": "C1",
                        "filter_backend": "packed_inline",
                        "survivor_tso_trim_policy": "full_tso_endpoint_then_parse_exact_17mer",
                        "reference_sha256": reference_hash,
                        "filter_settings_sha256": settings_hash,
                        "write_rejected_raw_pairs": False,
                        "statistics": {
                            "csm_gex_filter_positive_reads": 2,
                            "csm_gex_filter_positive_valid_barcode_reads": 1
                        },
                    }
                )
            )
            process = output / "gex_C1/process"
            process.mkdir(parents=True)
            bam = process / "barcode_headAligned_sorted.bam"
            bam.write_bytes(b"filtered-output")
            source_stats = source.parent / "pre_align_stats.csv"
            write_stats(source_stats, {
                "number_of_reads": 4,
                "reads_valid_barcode": 3,
                "reads_align_input": 3,
                "reads_align_unique": 2,
                "reads_align_multimap": 0,
            })
            sandbox_stats = process / "pre_align_stats.csv"
            write_stats(sandbox_stats, {
                "number_of_reads": 4,
                "reads_valid_barcode": 2,
                "reads_align_input": 3,
                "reads_align_unique": 1,
                "reads_align_multimap": 0,
                "reads_map_transcriptome": 1,
            })
            digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
            harmonized = process / "csm_harmonized_bam.complete.json"
            harmonized.write_text(
                json.dumps(
                    {
                        "schema_version": 2,
                        "complete": True,
                        "library": "C1",
                        "alignment_reused": True,
                        "source_bam": str(source),
                        "source_bam_sha256": digest(source),
                        "output_bam": str(bam),
                        "output_bam_sha256": digest(bam),
                        "filter_reference_sha256": reference_hash,
                        "filter_settings_sha256": settings_hash,
                        "survivor_tso_trim_policy": "not_applied_existing_alignment",
                        "name_sort_validated": True,
                        "records_before": 3,
                        "records_removed": 2,
                        "records_after": 1,
                        "qnames_before": 2,
                        "qnames_removed": 1,
                        "qnames_after": 1,
                        "qnames_with_original_header": 2,
                        "positive_id_records": 2,
                        "positive_unique_ids": 2,
                        "positive_valid_barcode_records": 1,
                        "positive_qnames_seen": 1,
                        "positive_ids_missing_from_bam": 1,
                        "source_unique_qnames": 2,
                        "source_multimap_qnames": 0,
                        "removed_unique_qnames": 1,
                        "removed_multimap_qnames": 0,
                        "alignment_stats_contract": {
                            "schema_version": 1,
                            "source_pre_align_stats": str(source_stats),
                            "source_pre_align_stats_sha256": digest(source_stats),
                            "sandbox_pre_align_stats": str(sandbox_stats),
                            "sandbox_pre_align_stats_sha256_at_handoff": digest(sandbox_stats),
                            "number_of_reads": 4,
                            "source_reads_valid_barcode": 3,
                            "sandbox_reads_valid_barcode": 2,
                            "valid_positive_records": 1,
                            "source_reads_align_input": 3,
                            "source_reads_align_unique": 2,
                            "source_reads_align_multimap": 0,
                            "removed_unique_qnames": 1,
                            "removed_multimap_qnames": 0,
                            "filtered_reads_align_unique": 1,
                            "filtered_reads_align_multimap": 0,
                        },
                    }
                )
            )
            self.assertEqual(
                len(validator.validate(config_path, only_library="C1")), 1
            )
            payload = json.loads(harmonized.read_text())
            payload["positive_ids_missing_from_bam"] = 0
            harmonized.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError, "missing-positive QNAME arithmetic"):
                validator.validate(config_path, only_library="C1")
            payload["positive_ids_missing_from_bam"] = 1
            payload["positive_qnames_seen"] = 0
            harmonized.write_text(json.dumps(payload))
            with self.assertRaisesRegex(
                ValueError, "removed QNAMEs do not equal observed positive QNAMEs"
            ):
                validator.validate(config_path, only_library="C1")
            payload["positive_qnames_seen"] = 1
            harmonized.write_text(json.dumps(payload))
            source.write_bytes(b"modified-production")
            with self.assertRaisesRegex(ValueError, "source BAM SHA-256 changed"):
                validator.validate(config_path, only_library="C1")

    def test_filter_and_harmonized_bam_contracts_validate_together(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output = root / "output"
            reference = root / "reference.csv"
            reference.write_text(
                "grna_seq,grna_name\nACGTACGTACGTACGTACGT,g1\n"
            )
            config = {
                "paths": {
                    "output_dir": str(output),
                    "grna_reference": str(reference),
                },
                "libraries": [{"sublibrary_name": "C1"}],
                "grna_structure": {
                    "fixed_seq_1": "AAGCAGTGGTATCAACGCAGAGTGAATGGG",
                    "grna_lengths": [20],
                },
                "gex_decontamination": {
                    "enabled": True,
                    "tso_matching": {"backend": "myers"},
                    "filter_matching": {
                        "backend": "packed_inline",
                        "guide_distance": "hamming",
                        "max_guide_mismatches": 1,
                        "allow_guide_indels": False,
                        "full_guide_required": True,
                    },
                    "existing_alignment": {"enabled": False},
                },
            }
            config_path = root / "config.yaml"
            config_path.write_text(yaml.safe_dump(config))
            classifier = filter_module.load_filter_from_config(config_path)
            reference_hash = classifier.index.stats.reference_sha256
            settings_hash = filter_module.filter_settings_sha256(
                classifier.settings
            )

            filter_dir = output / "gex_grna_filter/C1"
            filter_dir.mkdir(parents=True)
            (filter_dir / "rejected_pairs.complete.json").write_text(json.dumps({
                "schema_version": 4,
                "complete": True,
                "status": "complete",
                "library": "C1",
                "filter_backend": "packed_inline",
                "survivor_tso_trim_policy": "full_tso_endpoint_then_parse_exact_17mer",
                "reference_sha256": reference_hash,
                "filter_settings_sha256": settings_hash,
                "write_rejected_raw_pairs": False,
            }))
            process = output / "gex_C1/process"
            process.mkdir(parents=True)
            bam = process / "barcode_headAligned_sorted.bam"
            bam.write_bytes(b"fixture")
            harmonized = process / "csm_harmonized_bam.complete.json"
            harmonized.write_text(json.dumps({
                "complete": True,
                "library": "C1",
                "filter_reference_sha256": reference_hash,
                "filter_settings_sha256": settings_hash,
                "survivor_tso_trim_policy": "full_tso_endpoint_then_parse_exact_17mer",
                "output_bam": str(bam),
            }))
            self.assertEqual(len(validator.validate(config_path)), 1)

            filter_manifest = filter_dir / "rejected_pairs.complete.json"
            filter_payload = json.loads(filter_manifest.read_text())
            filter_payload["schema_version"] = 2
            filter_manifest.write_text(json.dumps(filter_payload))
            with self.assertRaisesRegex(ValueError, "schema is stale"):
                validator.validate(config_path)
            filter_payload["schema_version"] = 4
            filter_manifest.write_text(json.dumps(filter_payload))

            payload = json.loads(harmonized.read_text())
            payload["filter_settings_sha256"] = "stale"
            harmonized.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError, "removal policy is stale"):
                validator.validate(config_path)


if __name__ == "__main__":
    unittest.main()

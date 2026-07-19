import gzip
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from io import StringIO
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/gex_grna_bam_filter.py"
SPEC = importlib.util.spec_from_file_location("gex_grna_bam_filter", SCRIPT)
bam_filter = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = bam_filter
SPEC.loader.exec_module(bam_filter)


def sam_record(qname, flag=0, pos=1, tags="GX:Z:GENE1"):
    return (
        f"{qname}\t{flag}\tchr1\t{pos}\t255\t10M\t*\t0\t0\t"
        f"ACGTACGTAC\tFFFFFFFFFF\t{tags}\n"
    )


class BamFilterUnitTests(unittest.TestCase):
    def test_register_fresh_harmonized_bam_is_idempotent_and_checksummed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sandbox = root / "gex_C1"
            source = sandbox / "process/barcode_headAligned_sorted.bam"
            source.parent.mkdir(parents=True)
            source.write_bytes(b"fresh-bam-fixture")
            filter_manifest = root / "rejected_pairs.complete.json"
            filter_manifest.write_text(json.dumps({
                "schema_version": 4,
                "complete": True,
                "library": "C1",
                "run_id": "run1",
                "reference_sha256": "reference",
                "filter_settings_sha256": "settings",
                "survivor_tso_trim_policy": "full_tso_endpoint_then_parse_exact_17mer",
                "statistics": {
                    "csm_gex_filter_positive_reads": 4,
                    "csm_gex_filter_positive_valid_barcode_reads": 3,
                },
            }))
            payload, completion = bam_filter.register_fresh_harmonized_bam(
                source, sandbox, filter_manifest, library="C1"
            )
            self.assertFalse(payload["alignment_reused"])
            self.assertEqual(payload["valid_barcode_reads_removed_before_star"], 3)
            self.assertEqual(
                payload["survivor_tso_trim_policy"],
                "full_tso_endpoint_then_parse_exact_17mer",
            )
            self.assertTrue(completion.is_file())
            same, _ = bam_filter.register_fresh_harmonized_bam(
                source, sandbox, filter_manifest, library="C1"
            )
            self.assertEqual(same["output_bam_sha256"], payload["output_bam_sha256"])

    def test_raw_header_and_splitpipe_qname_normalize_to_same_id(self):
        raw = "@A00488:1:FC:1:1101:1000:2000 1:N:0:INDEX"
        qname = "01_02_03__T__1_2_3__AAA_CCC_GGG__UMI__STAMP__OH_@A00488:1:FC:1:1101:1000:2000"
        self.assertEqual(
            bam_filter.normalize_raw_cluster_id(raw),
            bam_filter.raw_cluster_from_splitpipe_qname(qname),
        )

    def test_stream_removes_every_record_and_preserves_survivors(self):
        qa = "01_01_01__T__1_1_1__A_C_G__U1__STAMP__OH_@readA"
        qb = "02_01_01__T__2_1_1__A_C_G__U2__STAMP__OH_@readB"
        lines = [
            "@HD\tVN:1.6\tSO:queryname\n",
            "@SQ\tSN:chr1\tLN:1000\n",
            sam_record(qa, 0, 10),
            sam_record(qa, 256, 20),
            sam_record(qa, 2048, 30),
            sam_record(qb, 0, 40, "GX:Z:GENE2\tCB:Z:CELL"),
        ]
        output = StringIO()
        stats, header = bam_filter.filter_sam_lines(
            lines, output, frozenset({"readA"})
        )
        self.assertEqual(stats.records_before, 4)
        self.assertEqual(stats.records_removed, 3)
        self.assertEqual(stats.records_after, 1)
        self.assertEqual(stats.qnames_removed, 1)
        self.assertEqual(stats.positive_qnames_seen, 1)
        self.assertEqual(header, "".join(lines[:2]))
        self.assertEqual(output.getvalue(), "".join(lines[:2] + [lines[-1]]))

    def test_unsorted_source_is_rejected(self):
        output = StringIO()
        with self.assertRaisesRegex(ValueError, "not query-name sorted"):
            bam_filter.filter_sam_lines(
                ["@HD\tSO:queryname\n", sam_record("read10"), sam_record("read2")],
                output, frozenset(),
            )

    def test_source_without_original_header_qnames_is_rejected(self):
        output = StringIO()
        with self.assertRaisesRegex(ValueError, "do not retain.*__OH_"):
            bam_filter.filter_sam_lines(
                ["@HD\tSO:queryname\n", sam_record("encoded_without_original")],
                output, frozenset({"readA"}),
            )

    def test_manifest_generates_ids_from_rejected_fastq(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            r1 = root / "rejected_R1.fastq.gz"
            with gzip.open(r1, "wt") as handle:
                handle.write("@readA 1:N:0:I\nACGT\n+\nFFFF\n")
            manifest = root / "rejected_pairs.complete.json"
            manifest.write_text(json.dumps({
                "status": "complete", "complete": True,
                "rejected_pair_count": 1, "rejected_r1": str(r1),
            }))
            evidence = bam_filter.ids_from_filter_manifest(manifest)
            self.assertEqual(evidence.ids, frozenset({"readA"}))
            self.assertEqual(evidence.record_count, 1)

    def test_manifest_errors_when_positive_ids_are_unavailable(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "rejected_pairs.complete.json"
            manifest.write_text(json.dumps({
                "status": "complete", "complete": True,
                "rejected_pair_count": 3,
            }))
            with self.assertRaisesRegex(ValueError, "does not provide positive IDs"):
                bam_filter.ids_from_filter_manifest(manifest)

            # rejected_pair_count may be zero solely because rejected FASTQ
            # output was disabled; positive filter statistics take priority.
            manifest.write_text(json.dumps({
                "status": "complete", "complete": True,
                "rejected_pair_count": 0,
                "statistics": {"csm_gex_filter_positive_reads": 3},
            }))
            with self.assertRaisesRegex(ValueError, "does not provide positive IDs"):
                bam_filter.ids_from_filter_manifest(manifest)

    def test_completed_zero_positive_manifest_is_valid(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "rejected_pairs.complete.json"
            manifest.write_text(json.dumps({
                "status": "complete", "complete": True,
                "library": "C1", "rejected_pair_count": 0,
            }))
            evidence = bam_filter.ids_from_filter_manifest(manifest)
            self.assertEqual(evidence.ids, frozenset())
            self.assertEqual(evidence.library, "C1")


class SamtoolsFallbackIntegrationTest(unittest.TestCase):
    def test_cli_writes_new_bam_and_completion_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fake_samtools = root / "samtools"
            fake_samtools.write_text("""#!/usr/bin/env python3
import pathlib, sys
args = sys.argv[1:]
if '-o' in args:
    out = pathlib.Path(args[args.index('-o') + 1])
    out.write_text(sys.stdin.read())
elif '-H' in args:
    for line in pathlib.Path(args[-1]).read_text().splitlines(True):
        if line.startswith('@'):
            sys.stdout.write(line)
elif '-c' in args:
    print(sum(not line.startswith('@') and bool(line.strip()) for line in pathlib.Path(args[-1]).read_text().splitlines()))
elif '-h' in args:
    sys.stdout.write(pathlib.Path(args[-1]).read_text())
else:
    raise SystemExit(2)
""")
            fake_samtools.chmod(0o755)

            qa = "01_01_01__T__1_1_1__A_C_G__U1__STAMP__OH_@readA"
            qb = "02_01_01__T__2_1_1__A_C_G__U2__STAMP__OH_@readB"
            source = root / "barcode_headAligned_anno.bam"
            source_text = (
                "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:1000\n"
                + sam_record(qa, 0, 10)
                + sam_record(qa, 256, 20)
                + sam_record(qb, 0, 30, "GX:Z:GENE2\tGN:Z:G2")
            )
            source.write_text(source_text)
            source_hash = bam_filter._sha256(source)

            r1 = root / "rejected_R1.fastq.gz"
            with gzip.open(r1, "wt") as handle:
                handle.write("@readA 1:N:0:I\nACGT\n+\nFFFF\n")
            manifest = root / "rejected_pairs.complete.json"
            manifest.write_text(json.dumps({
                "status": "complete", "complete": True,
                "rejected_pair_count": 1, "rejected_r1": str(r1),
            }))
            sandbox = root / "sandbox"

            subprocess.run([
                sys.executable, str(SCRIPT),
                "--source-bam", str(source),
                "--output-dir", str(sandbox),
                "--filter-manifest", str(manifest),
                "--backend", "samtools", "--samtools", str(fake_samtools),
            ], check=True, text=True, capture_output=True)

            output = sandbox / "process/barcode_headAligned_sorted.bam"
            expected = (
                "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:1000\n"
                + sam_record(qb, 0, 30, "GX:Z:GENE2\tGN:Z:G2")
            )
            self.assertEqual(output.read_text(), expected)
            self.assertEqual(bam_filter._sha256(source), source_hash)
            completion = json.loads(
                (sandbox / "process/csm_harmonized_bam.complete.json").read_text()
            )
            self.assertTrue(completion["complete"])
            self.assertTrue(completion["alignment_reused"])
            self.assertEqual(completion["source_bam_kind"], "annotated")
            self.assertEqual(completion["records_before"], 3)
            self.assertEqual(completion["records_removed"], 2)
            self.assertEqual(completion["records_after"], 1)

            # Safe reruns never overwrite the installed output.
            rerun = subprocess.run([
                sys.executable, str(SCRIPT),
                "--source-bam", str(source), "--output-dir", str(sandbox),
                "--filter-manifest", str(manifest),
                "--backend", "samtools", "--samtools", str(fake_samtools),
            ], text=True, capture_output=True)
            self.assertNotEqual(rerun.returncode, 0)
            self.assertIn("choose a new output directory", rerun.stderr)


if __name__ == "__main__":
    unittest.main()

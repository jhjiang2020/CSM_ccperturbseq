#!/usr/bin/env python3
"""Focused tests for the C2 gRNA-to-GEX case-study helpers."""

import gzip
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.c2_gex_grna_case_study import (
    IBAR_LENGTH,
    PUBLISHED_CONTROLLED_COLUMNS,
    PUBLISHED_GEX_COLUMNS,
    PUBLISHED_REPRESENTATIVE_COLUMNS,
    R2_LINKER_2,
    R2_LINKER_2_RC,
    SPLITPIPE_TSO,
    STEM_VARIANTS,
    TARGETS,
    TSO,
    TRACR_VARIANTS,
    aligned_base_components,
    as_bool,
    deserialize_blocks,
    intervals_overlap,
    merge_intervals,
    normalize_gene_name,
    parse_cigar_blocks,
    load_bundle,
    read_fastq_records,
    raw_component_segments,
    select_representative,
    serialize_blocks,
    splitpipe_trim,
    summarize_positions,
    write_deterministic_dataframe,
)


ROOT = Path(__file__).resolve().parents[1]
PUBLISHED_DATA = (
    ROOT / "docs/assets/c2_gex_grna_misalignment/data"
)


class CigarBlockTest(unittest.TestCase):
    def test_soft_clips_do_not_consume_reference(self):
        self.assertEqual(parse_cigar_blocks(100, "3S21M5S"), [(100, 120)])

    def test_splice_creates_two_reference_blocks(self):
        self.assertEqual(
            parse_cigar_blocks(100, "13M260539N13M"),
            [(100, 112), (260652, 260664)],
        )

    def test_deletion_excludes_unaligned_reference_bases(self):
        self.assertEqual(parse_cigar_blocks(10, "5M2D4M"), [(10, 14), (17, 20)])

    def test_serialization_and_overlap(self):
        blocks = [(100, 120), (200, 220)]
        self.assertEqual(deserialize_blocks(serialize_blocks(blocks)), blocks)
        self.assertTrue(intervals_overlap(blocks, 119, 150))
        self.assertFalse(intervals_overlap(blocks, 121, 199))


class RepresentativeSelectionTest(unittest.TestCase):
    def setUp(self):
        self.guide = "G" * 20
        self.stem = "TTTTCCCCAAAAGGGGTTT"

    def row(self, truth_id, cigar="20M4S", guide=None):
        guide = guide or self.guide
        return {
            "truth_id": truth_id,
            "original_read_id": f"raw-{truth_id}",
            "guide_query": guide,
            "cigar": cigar,
            "mapq": 255,
            "nm": 0,
            "target_concordant": True,
            "production_filter_removed": True,
            "tso_offset": 0,
        }

    def raw(self, truth_id, exact=True, guide=None):
        guide = guide or self.guide
        prefix = TSO if exact else "C" + TSO[1:]
        sequence = (prefix + guide + self.stem + "ACGT" * 30)[:150]
        return {f"raw-{truth_id}": {"sequence": sequence, "quality": "I" * 150}}

    def test_dominant_pattern_then_lowest_truth_id(self):
        rows = [self.row("GQ3"), self.row("GQ2"), self.row("GQ1", cigar="19M5S")]
        raw = {}
        for truth_id in ("GQ1", "GQ2", "GQ3"):
            raw.update(self.raw(truth_id))
        selected = select_representative(rows, raw, [self.stem])
        self.assertEqual(selected["truth_id"], "GQ2")

    def test_falls_back_when_dominant_pattern_lacks_exact_tso(self):
        rows = [self.row("GQ1"), self.row("GQ2"), self.row("GQ3", cigar="19M5S")]
        raw = {}
        raw.update(self.raw("GQ1", exact=False))
        raw.update(self.raw("GQ2", exact=False))
        raw.update(self.raw("GQ3", exact=True))
        selected = select_representative(rows, raw, [self.stem])
        self.assertEqual(selected["truth_id"], "GQ3")


class SmallHelperTest(unittest.TestCase):
    def test_splitpipe_exact_tso_trim(self):
        sequence = TSO + "G" * 20 + "TTTT"
        trimmed, index = splitpipe_trim(sequence)
        self.assertEqual(index, 13)
        self.assertEqual(trimmed, "G" * 20 + "TTTT")
        self.assertIn(SPLITPIPE_TSO, TSO)

    def test_quote_normalization(self):
        self.assertEqual(normalize_gene_name('""KLF2""'), "KLF2")
        self.assertEqual(normalize_gene_name('"NFKB1"'), "NFKB1")

    def test_serialized_booleans(self):
        self.assertTrue(as_bool(True))
        self.assertTrue(as_bool("true"))
        self.assertFalse(as_bool(False))
        self.assertFalse(as_bool("False"))

    def test_truth_fastq_is_keyed_before_metadata(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "truth.fastq"
            path.write_text("@GQ000000027|C=cell|T=KLF2\nACGT\n+\nIIII\n")
            records = read_fastq_records(path)
        self.assertEqual(records["GQ000000027"]["sequence"], "ACGT")

    def test_raw_components_include_opposite_end_barcode_readthrough(self):
        guide = "G" * 20
        sequence = (
            TSO
            + guide
            + STEM_VARIANTS[1]
            + "ACGTACG"
            + TRACR_VARIANTS[1]
            + "A" * 28
            + "C" * 8
            + R2_LINKER_2_RC
            + "G" * 8
            + "CT"
        )
        self.assertEqual(len(sequence), 150)
        self.assertEqual(
            raw_component_segments(sequence, guide),
            [
                (0, 30, "TSO"),
                (30, 50, "guide"),
                (50, 69, "stem"),
                (69, 76, "ibar"),
                (76, 92, "tracr"),
                (92, 120, "bridge"),
                (120, 128, "r2_barcode_rc"),
                (128, 140, "r2_linker_rc"),
                (140, 148, "r2_barcode_rc"),
                (148, 150, "r2_linker_rc"),
            ],
        )
        self.assertEqual(IBAR_LENGTH, 7)
        complement = str.maketrans("ACGT", "TGCA")
        self.assertEqual(R2_LINKER_2.translate(complement)[::-1], R2_LINKER_2_RC)

    def test_aligned_bases_inherit_guide_then_stem_annotation(self):
        self.assertEqual(aligned_base_components(20, 20), ["guide"] * 20)
        self.assertEqual(
            aligned_base_components(21, 20),
            ["guide"] * 20 + ["stem"],
        )

    def test_position_summary_is_sorted_and_uses_dominant_guide(self):
        frame = pd.DataFrame(
            [
                {"assigned_gene_name": "KLF2", "target_concordant": "True", "start_1": 20, "guide_query": "B"},
                {"assigned_gene_name": "KLF2", "target_concordant": True, "start_1": 10, "guide_query": "A"},
                {"assigned_gene_name": "KLF2", "target_concordant": True, "start_1": 20, "guide_query": "B"},
                {"assigned_gene_name": "KLF2", "target_concordant": False, "start_1": 30, "guide_query": "C"},
            ]
        )
        positions, guides = summarize_positions(frame, TARGETS[0])
        self.assertIn("hg38_chr19:10", positions)
        self.assertLess(positions.index(":10"), positions.index(":20"))
        self.assertEqual(guides, "`A`<br>`B`")

    def test_merge_intervals(self):
        self.assertEqual(
            merge_intervals([(10, 20), (18, 30), (31, 40), (100, 110)]),
            [(10, 40), (100, 110)],
        )

    def test_deterministic_gzip_dataframe(self):
        frame = pd.DataFrame([{"name": "example", "count": 3}])
        with tempfile.TemporaryDirectory() as temporary:
            first = Path(temporary) / "first.tsv.gz"
            second = Path(temporary) / "second.tsv.gz"
            write_deterministic_dataframe(frame, first)
            write_deterministic_dataframe(frame, second)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            with gzip.open(first, "rt") as handle:
                self.assertEqual(handle.read(), "name\tcount\nexample\t3\n")

    def test_published_bundle_is_minimal_and_sanitized(self):
        manifest = json.loads(
            (PUBLISHED_DATA / "extraction.complete.json").read_text()
        )
        self.assertTrue(manifest["complete"])
        self.assertTrue(manifest["sanitized"])
        self.assertEqual(
            set(manifest["output_files"]),
            {path.name for path in PUBLISHED_DATA.iterdir()
             if path.name != "extraction.complete.json"},
        )
        for name, metadata in manifest["output_files"].items():
            digest = hashlib.sha256((PUBLISHED_DATA / name).read_bytes()).hexdigest()
            self.assertEqual(digest, metadata["sha256"], name)
        self.assertEqual(
            set(
                pd.read_csv(
                    PUBLISHED_DATA / "controlled_short54_alignments.tsv.gz",
                    sep="\t",
                ).columns
            ),
            set(PUBLISHED_CONTROLLED_COLUMNS),
        )
        self.assertEqual(
            set(
                pd.read_csv(
                    PUBLISHED_DATA / "representative_reads.tsv", sep="\t"
                ).columns
            ),
            set(PUBLISHED_REPRESENTATIVE_COLUMNS),
        )
        self.assertEqual(
            set(
                pd.read_csv(
                    PUBLISHED_DATA / "putative_genuine_gex_reads.tsv.gz",
                    sep="\t",
                ).columns
            ),
            set(PUBLISHED_GEX_COLUMNS),
        )
        for path in sorted(PUBLISHED_DATA.iterdir()):
            raw = gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes()
            text = raw.decode("utf-8")
            self.assertNotIn("LH00146:", text, path.name)
            self.assertNotIn("/medpop/", text, path.name)
            self.assertNotIn("/broad/", text, path.name)

        forbidden_columns = {
            "original_read_id",
            "cell",
            "barcode_indices",
            "umi",
            "query_sequence",
        }
        self.assertTrue(
            forbidden_columns.isdisjoint(PUBLISHED_CONTROLLED_COLUMNS)
        )
        self.assertTrue(
            forbidden_columns.isdisjoint(PUBLISHED_REPRESENTATIVE_COLUMNS)
        )
        self.assertTrue(forbidden_columns.isdisjoint(PUBLISHED_GEX_COLUMNS))

    def test_bundle_loader_rejects_hash_mismatch(self):
        with tempfile.TemporaryDirectory() as temporary:
            copied = Path(temporary) / "data"
            shutil.copytree(PUBLISHED_DATA, copied)
            with (copied / "global_summary.json").open("a") as handle:
                handle.write(" ")
            with self.assertRaisesRegex(RuntimeError, "SHA256 global_summary.json"):
                load_bundle(copied)

    def test_bundle_loader_requires_hash_for_every_input(self):
        with tempfile.TemporaryDirectory() as temporary:
            copied = Path(temporary) / "data"
            shutil.copytree(PUBLISHED_DATA, copied)
            manifest_path = copied / "extraction.complete.json"
            manifest = json.loads(manifest_path.read_text())
            del manifest["output_files"]["global_summary.json"]
            manifest_path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(
                RuntimeError, "lacks hashes.*global_summary.json"
            ):
                load_bundle(copied)

    def test_committed_report_and_figures_render_deterministically(self):
        with tempfile.TemporaryDirectory() as temporary:
            rendered = Path(temporary)
            subprocess.check_call(
                [
                    sys.executable,
                    str(ROOT / "scripts/c2_gex_grna_case_study.py"),
                    "render",
                    "--data-dir",
                    str(PUBLISHED_DATA),
                    "--docs-dir",
                    str(rendered),
                ]
            )
            expected = [
                Path("C2_GEX_GRNA_MISALIGNMENT_CASE_STUDY.md"),
                *(
                    Path("assets/c2_gex_grna_misalignment") / f"{stem}.{suffix}"
                    for stem in ("klf2_locus", "nfkb1_locus", "read_transformation")
                    for suffix in ("png", "svg")
                ),
            ]
            for relative in expected:
                self.assertEqual(
                    (ROOT / "docs" / relative).read_bytes(),
                    (rendered / relative).read_bytes(),
                    str(relative),
                )


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Synthetic regression tests for GEX gRNA contamination QC."""

import copy
import gzip
import json
import tempfile
import unittest
from collections import Counter
from pathlib import Path

import pandas as pd
import yaml

from scripts.qc_gex_grna_contamination import (
    BarcodeCall,
    FastqRecord,
    Confusion,
    GeneAnnotator,
    GuideModel,
    HyperLogLog,
    PilotBarcodeCorrector,
    aggregate,
    benchmark_library,
    cluster_bootstrap_ratio_interval,
    collapse_umi_payloads,
    find_tso,
    filter_library,
    iter_fastq_pairs,
    parse_alignment_truth_qname,
    paired_cluster_bootstrap_change,
    splitpipe_processed_record,
    stable_fraction,
    validate_step1a_barcode_evidence_with_qc_decoder,
    _cross_source_overlap,
    _production_assignment_sets,
)
from tests.release_test_helpers import (
    artifact_semantic_digests,
    benchmark_artifact_digests,
)


ROOT = Path(__file__).resolve().parents[1]


def mutate_base(sequence, index):
    replacement = "A" if sequence[index] != "A" else "C"
    return sequence[:index] + replacement + sequence[index + 1 :]


def write_fastq(path, records):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as handle:
        for read_id, sequence in records:
            handle.write(f"@{read_id}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def count_fastq(path):
    return sum(1 for _ in iter_fastq_pairs(str(path)))


class GuideModelTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(ROOT / "config.yaml") as handle:
            cls.base_cfg = yaml.safe_load(handle)
        with open(ROOT / "examples/qc/qc_config_gex_grna_contamination.example.yaml") as handle:
            cls.qc_cfg = yaml.safe_load(handle)["gex_grna_contamination_qc"]
        cls.reference = pd.read_csv(ROOT / "references/grna_ibar_reference.csv")
        cls.model = GuideModel(
            cls.reference,
            cls.base_cfg["grna_structure"],
            cls.qc_cfg["full_length_truth"],
        )

    def construct(self, row, suffix="ACGTACGTACGT"):
        variant = int(row["position"])
        return (
            self.model.tso
            + row["grna_seq"]
            + self.model.stems[variant]
            + row["ibar_seq"]
            + self.model.tracrs[variant]
            + suffix
        )

    def test_full_length_truth_both_variants(self):
        for variant in (1, 2):
            row = self.reference[self.reference["position"] == variant].iloc[0]
            evidence = self.model.classify_full(self.construct(row).lower())
            self.assertTrue(evidence.structure_truth)
            self.assertTrue(evidence.known_structure_truth)
            self.assertTrue(evidence.complete_reference_truth)
            self.assertEqual(evidence.variant, variant)
            self.assertEqual(evidence.tier, "complete_known_construct")

    def test_step1a_evidence_matches_qc_decoder_without_full_rescan(self):
        class FakeCorrector:
            def call(self, sequence):
                if sequence.startswith("A"):
                    return BarcodeCall(
                        True,
                        cell="A1_B2_C3",
                        barcode_indices="1_2_3",
                        umi="AAAAAAAAAA",
                    )
                return BarcodeCall(False, umi="CCCCCCCCCC")

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            r1 = root / "rejected_R1.fastq.gz"
            r2 = root / "rejected_R2.fastq.gz"
            write_fastq(r1, [("read1", "ACGT"), ("read2", "TGCA")])
            write_fastq(r2, [("read1", "ATTT"), ("read2", "CTTT")])
            evidence = root / "valid_positive_barcode_evidence.tsv.gz"
            with gzip.open(evidence, "wt") as handle:
                handle.write("read_id\tbc_wells\tbc_bci\tumi\tr1_seq\tr1_qual\n")
                handle.write(
                    "read1\tA1_B2_C3\t1_2_3\tAAAAAAAAAA\tACGT\tIIII\n"
                )

            def digest(path):
                import hashlib
                return hashlib.sha256(path.read_bytes()).hexdigest()

            manifest = root / "rejected_pairs.complete.json"
            payload = {
                "complete": True,
                "library": "L1",
                "rejected_pair_count": 2,
                "barcode_evidence_rows": 1,
                "rejected_r1": str(r1),
                "rejected_r1_sha256": digest(r1),
                "rejected_r2": str(r2),
                "rejected_r2_sha256": digest(r2),
                "barcode_evidence": str(evidence),
                "barcode_evidence_sha256": digest(evidence),
            }
            manifest.write_text(json.dumps(payload))
            result = validate_step1a_barcode_evidence_with_qc_decoder(
                {}, str(root / "unused.yaml"), "L1", manifest, payload,
                barcode_corrector=FakeCorrector(),
            )
            self.assertTrue(result["all_fields_identical"])
            self.assertEqual(result["rejected_pairs_decoded"], 2)
            self.assertEqual(result["valid_rows_compared"], 1)

    def test_all_reference_constructs_detected_at_54_offsets_zero_to_two(self):
        rows = self.reference.drop_duplicates(["grna_seq", "position"])
        short_cfg = self.qc_cfg["short_read_filters"]
        for row in rows.to_dict("records"):
            construct = self.construct(row)
            for offset in range(3):
                short = self.model.short_features("C" * offset + construct, 54, short_cfg)
                self.assertTrue(
                    short["calls"]["tso_template_boundary_d0"],
                    msg=f"missed {row['grna_name']} position={row['position']} offset={offset}",
                )
                self.assertGreaterEqual(
                    short["template_prefix_hit"].minimum_boundary_bases, 1
                )

    def test_offset_three_is_diagnostic_not_primary(self):
        row = self.reference.iloc[0]
        short = self.model.short_features(
            "CCC" + self.construct(row), 54, self.qc_cfg["short_read_filters"]
        )
        self.assertFalse(short["calls"]["tso_template_prefix_d0"])
        self.assertTrue(short["calls"]["tso_anywhere_template_prefix_d0"])

    def test_tso_and_spacer_sequencing_errors(self):
        row = self.reference.iloc[0]
        sequence = self.construct(row)
        sequence = mutate_base(sequence, 4)
        sequence = mutate_base(sequence, len(self.model.tso) + 5)
        short = self.model.short_features(
            sequence, 54, self.qc_cfg["short_read_filters"]
        )
        self.assertFalse(short["calls"]["tso_template_boundary_d0"])
        self.assertTrue(short["calls"]["tso_template_boundary_d1"])

    def test_tso_positive_cdna_is_not_called_by_guide_filter(self):
        sequence = self.model.tso + "ACGT" * 35
        full = self.model.classify_full(sequence)
        short = self.model.short_features(
            sequence, 54, self.qc_cfg["short_read_filters"]
        )
        self.assertTrue(full.tso_found)
        self.assertFalse(full.structure_truth)
        self.assertFalse(short["calls"]["tso_template_boundary_d2"])
        self.assertTrue(short["calls"]["tso_only_d3"])

    def test_21mer_ending_g_is_not_forced_to_20mer(self):
        candidates = self.reference[
            (self.reference["grna_seq"].str.len() == 21)
            & self.reference["grna_seq"].str.endswith("G")
        ]
        self.assertFalse(candidates.empty)
        row = candidates.iloc[0]
        short = self.model.short_features(
            self.construct(row), 54, self.qc_cfg["short_read_filters"]
        )
        self.assertEqual(short["guide_hit"].length, 21)
        self.assertTrue(short["calls"]["tso_template_boundary_d0"])

    def test_shared_spacer_remains_binary_positive(self):
        shared = (
            self.reference.groupby("grna_seq")["grna_name"]
            .nunique()
            .loc[lambda values: values > 1]
        )
        self.assertFalse(shared.empty)
        row = self.reference[self.reference["grna_seq"] == shared.index[0]].iloc[0]
        short = self.model.short_features(
            self.construct(row), 54, self.qc_cfg["short_read_filters"]
        )
        self.assertTrue(short["calls"]["tso_template_boundary_d0"])
        self.assertTrue(short["template_prefix_hit"].identity_ambiguous)

    def test_post_splitpipe_orientation_and_structure(self):
        row = self.reference.iloc[0]
        raw = self.construct(row)
        post = raw[len(self.model.tso) :]
        evidence = self.model.classify_post_splitpipe(post)
        self.assertTrue(evidence.structure_truth)
        self.assertTrue(evidence.complete_reference_truth)


class UtilityTest(unittest.TestCase):
    def test_dependency_free_tso_matching_supports_offsets_and_edits(self):
        tso = "AAGCAGTGGTATCAACGCAGAGTGAATGGG"
        start, end, distance = find_tso("CC" + tso + "ACGT", tso, 3, 50)
        self.assertEqual((start, end, distance), (2, 2 + len(tso), 0))

        substituted = list(tso)
        for index in (2, 13, 25):
            substituted[index] = "A" if substituted[index] != "A" else "C"
        start, end, distance = find_tso(
            "C" + "".join(substituted) + "ACGT", tso, 3, 50
        )
        self.assertEqual((start, end, distance), (1, 1 + len(tso), 3))

        inserted = tso[:12] + "C" + tso[12:]
        start, end, distance = find_tso("C" + inserted + "ACGT", tso, 3, 50)
        self.assertEqual((start, end, distance), (1, 2 + len(tso), 1))

        deleted = tso[:12] + tso[13:]
        start, end, distance = find_tso("C" + deleted + "ACGT", tso, 3, 50)
        self.assertEqual((start, end, distance), (1, len(tso), 1))

    def test_splitpipe_tso_processing_truncates_first(self):
        tso = "AAGCAGTGGTATCAACGCAGAGTGAATGGG"
        guide_and_stem = "G" * 21 + "ACGTACGTACGT"
        raw_sequence = "CC" + tso + guide_and_stem
        record = FastqRecord(
            "@truth", raw_sequence, "+", "I" * len(raw_sequence)
        )
        sequence, quality, trim_index = splitpipe_processed_record(
            record, raw_length=54
        )
        self.assertEqual(trim_index, 15)
        self.assertEqual(sequence, guide_and_stem[:22])
        self.assertEqual(len(sequence), len(quality))

        terminal_error = mutate_base(raw_sequence, 2 + 20)
        error_record = FastqRecord(
            "@error", terminal_error, "+", "I" * len(terminal_error)
        )
        sequence, quality, trim_index = splitpipe_processed_record(
            error_record, raw_length=54
        )
        self.assertEqual(trim_index, -1)
        self.assertEqual(len(sequence), 54)
        self.assertEqual(len(sequence), len(quality))

    def test_hyperloglog_cardinality_and_duplicate_stability(self):
        sketch = HyperLogLog(precision=12)
        for index in range(20_000):
            sketch.add("cell", index)
            sketch.add("cell", index)
        self.assertLess(abs(sketch.estimate() - 20_000) / 20_000, 0.08)

    def test_barcode_correction_preserves_pb_identity(self):
        rows = []
        round_sequences = {
            1: [(1, 1, "AAAAAAAA"), (97, 1, "CCCCCCCC")],
            2: [(1, 1, "GGGGGGGG"), (2, 2, "TTTTTTTT")],
            3: [(1, 1, "ACACACAC"), (2, 2, "TGTGTGTG")],
        }
        for bc_round, values in round_sequences.items():
            for bci, wind, sequence in values:
                rows.append(
                    {
                        "bci": bci,
                        "sequence": sequence,
                        "wind": wind,
                        "bc_round": bc_round,
                    }
                )
        corrector = PilotBarcodeCorrector(
            pd.DataFrame(rows),
            "NNNNNNNNNN33333333ATGAGGGGTCAG22222222TCCAACCACCTC11111111",
        )

        def r2(bc1):
            return (
                "ACGTACGTAC"
                + "ACACACAC"
                + "ATGAGGGGTCAG"
                + "GGGGGGGG"
                + "TCCAACCACCTC"
                + bc1
            )

        r_call = corrector.call(r2("AAAAAAAA"))
        t_call = corrector.call(r2("CCCCCCCC"))
        self.assertTrue(r_call.valid and t_call.valid)
        self.assertEqual(r_call.cell, t_call.cell)
        self.assertEqual(r_call.barcode_indices, "1_1_1")
        self.assertEqual(t_call.barcode_indices, "97_1_1")

        corrector.perfect_combination_counts[(
            "AAAAAAAA",
            "GGGGGGGG",
            "ACACACAC",
        )] = 10
        corrector.perfect_count_threshold = 5
        rescued = corrector.call(r2("AAAAAAAT"))
        self.assertTrue(rescued.valid)
        self.assertEqual(rescued.max_distance, 1)
        self.assertEqual(rescued.barcode_indices, "1_1_1")

    def test_barcode_correction_matches_splitpipe_indel_rescue(self):
        rows = []
        references = {
            1: (114, 18, "CTTCATCA"),
            2: (92, 92, "GAACAGGC"),
            3: (20, 20, "CTGTCTAT"),
        }
        for bc_round, (bci, wind, sequence) in references.items():
            rows.append({
                "bci": bci,
                "sequence": sequence,
                "wind": wind,
                "bc_round": bc_round,
            })
        corrector = PilotBarcodeCorrector(
            pd.DataFrame(rows),
            "NNNNNNNNNN33333333ATGAGGGGTCAG22222222TCCAACCACCTC11111111",
        )
        combination = tuple(references[index][2] for index in (1, 2, 3))
        corrector.perfect_combination_counts[combination] = 10
        corrector.perfect_count_threshold = 5
        observed = (
            "TTTTTATACC"
            + "TGTCTATA"  # round 3: terminal shift, Levenshtein 2
            + "ATGAGGGGTCAG"
            + "AACAGGCT"  # round 2: terminal shift, Levenshtein 2
            + "TCCAACCACCTC"
            + "CTTCATAC"  # round 1: two substitutions
        )
        call = corrector.call(observed)
        self.assertTrue(call.valid)
        self.assertEqual(call.cell, "18_92_20")
        self.assertEqual(call.barcode_indices, "114_92_20")
        self.assertEqual(call.max_distance, 2)

    def test_umi_collapse_matches_abundance_directed_rule(self):
        payloads = {
            "AAAAAAAAAA": Counter(total=10, guide=10, guide_target=5),
            "AAAAAAAAAT": Counter(total=2, guide=2, guide_target=1),
            "CCCCCCCCCC": Counter(total=3, guide=3),
            "CCCCCCCCCA": Counter(total=3),
        }
        clusters = collapse_umi_payloads(payloads)
        self.assertEqual(sorted(cluster["total"] for cluster in clusters), [3, 3, 12])
        merged = next(cluster for cluster in clusters if cluster["total"] == 12)
        self.assertEqual(merged["guide"], 12)
        self.assertEqual(merged["guide_target"], 6)

    def test_umi_collapse_large_neighbor_index_branch(self):
        alphabet = "ACGT"

        def block_code(value):
            digits = []
            for divisor in (16, 4, 1):
                digits.append(alphabet[(value // divisor) % 4] * 4)
            return "".join(digits)

        payloads = {
            "A" * 12: Counter(total=10, guide=10),
            "A" * 11 + "C": Counter(total=2, guide=2),
        }
        for value in range(1, 50):
            payloads[block_code(value)] = Counter(total=1)
        self.assertEqual(len(payloads), 51)
        clusters = collapse_umi_payloads(payloads)
        self.assertEqual(len(clusters), 50)
        self.assertIn(12, [cluster["total"] for cluster in clusters])
        self.assertEqual(sum(cluster["total"] for cluster in clusters), 61)

    def test_alignment_truth_qname_round_trip(self):
        truth = parse_alignment_truth_qname(
            "GQ000000001|C=01_02_03|B=1_2_3|U=ACGTACGTAC|T=TP53,MDM2|Q=guide|X=1|O=2"
        )
        self.assertEqual(truth["truth_id"], "GQ000000001")
        self.assertEqual(truth["targets"], ("TP53", "MDM2"))
        self.assertTrue(truth["tso_trimmed"])

    def test_cluster_bootstrap_and_paired_change(self):
        counts = {
            (f"cell{index}", "1_1_1"): Counter(
                input_reads=10,
                false_gene_assigned_reads=2 if index % 2 == 0 else 4,
            )
            for index in range(20)
        }
        interval = cluster_bootstrap_ratio_interval(
            counts,
            "false_gene_assigned_reads",
            "input_reads",
            200,
            "test",
        )
        self.assertLess(interval[0], 0.3)
        self.assertGreater(interval[1], 0.3)
        self.assertEqual(
            interval,
            cluster_bootstrap_ratio_interval(
                counts,
                "false_gene_assigned_reads",
                "input_reads",
                200,
                "test",
            ),
        )

        full = pd.DataFrame(
            {
                "cell": [f"cell{index}" for index in range(10)],
                "barcode_indices": ["1_1_1"] * 10,
                "input_reads": [10] * 10,
                "false_gene_assigned_reads": [1] * 10,
            }
        )
        short = full.copy()
        short["false_gene_assigned_reads"] = 3
        change = paired_cluster_bootstrap_change(
            full,
            short,
            "false_gene_assigned_reads",
            "input_reads",
            100,
            "paired-test",
        )
        self.assertAlmostEqual(
            change["absolute_delta_cluster_bootstrap_95_low"], 0.2
        )
        self.assertAlmostEqual(
            change["absolute_delta_cluster_bootstrap_95_high"], 0.2
        )

    def test_gene_annotator_strand_span_and_exon_assignment(self):
        gene_info = {
            "gene_bins": {
                "chr1:0:+": ["gene_plus"],
                "chr1:0:-": ["gene_minus"],
            },
            "gtf_dict_stepsize": 100,
            "gene_starts": {"gene_plus": 10, "gene_minus": 10},
            "gene_ends": {"gene_plus": 90, "gene_minus": 90},
            "gene_id_to_name": {
                "gene_plus": '\"PLUS\"',
                "gene_minus": "MINUS",
            },
            "gene_id_to_biotype": {
                "gene_plus": "protein_coding",
                "gene_minus": "protein_coding",
            },
            "genes_to_exons": {
                "gene_plus": {20: 30},
                "gene_minus": {40: 50},
            },
        }

        class FakeRead:
            def __init__(self, positions, reverse=False, unmapped=False):
                self._positions = positions
                self.is_reverse = reverse
                self.is_unmapped = unmapped
                self.reference_name = "chr1"

            def get_reference_positions(self, full_length=False):
                del full_length
                return list(self._positions)

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "gene_info.json"
            path.write_text(json.dumps(gene_info))
            annotator = GeneAnnotator(str(path))
            self.assertEqual(annotator.names["gene_plus"], "PLUS")
            gene, name, exonic = annotator.assign(FakeRead(list(range(20, 31))))
            self.assertEqual((gene, name, exonic), ("gene_plus", "PLUS", 1.0))
            gene, name, exonic = annotator.assign(
                FakeRead(list(range(40, 51)), reverse=True)
            )
            self.assertEqual((gene, name, exonic), ("gene_minus", "MINUS", 1.0))
            self.assertEqual(annotator.assign(FakeRead([5, 6, 7]))[0], "")
            self.assertEqual(annotator.assign(FakeRead([], unmapped=True))[0], "")

    def test_paired_fastq_id_validation_and_truncation_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            r1 = root / "r1.fastq.gz"
            r2 = root / "r2.fastq.gz"
            write_fastq(r1, [("x/1", "A" * 60), ("y/1", "C" * 60)])
            write_fastq(r2, [("x/2", "G" * 58), ("y/2", "T" * 58)])
            pairs = list(iter_fastq_pairs(str(r1), str(r2)))
            self.assertEqual(len(pairs), 2)
            self.assertEqual(pairs[0][0].quality, "I" * 60)

            write_fastq(r2, [("not_x/2", "G" * 58)])
            with self.assertRaisesRegex(ValueError, "ID mismatch"):
                list(iter_fastq_pairs(str(r1), str(r2)))

    def test_fastq_sequence_quality_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.fastq"
            path.write_text("@x\nACGT\n+\nIII\n")
            with self.assertRaisesRegex(ValueError, "sequence/quality"):
                list(iter_fastq_pairs(str(path)))

    def test_stable_fraction_is_deterministic(self):
        self.assertEqual(stable_fraction("read-1"), stable_fraction("read-1"))
        self.assertNotEqual(stable_fraction("read-1"), stable_fraction("read-2"))

    def test_hard_negative_metric(self):
        confusion = Confusion()
        confusion.update(False, True, hard_negative=True)
        confusion.update(False, False, hard_negative=True)
        metrics = confusion.as_dict()
        self.assertEqual(metrics["hard_negative_n"], 2)
        self.assertEqual(metrics["hard_negative_fpr"], 0.5)


class MiniPipelineTest(unittest.TestCase):
    def test_production_assignment_sets_preserve_dual_pair_calls(self):
        frame = pd.DataFrame([
            {"bc_wells": "cell1", "pair_id": "pairA", "grna_name": "g1"},
            {"bc_wells": "cell1", "pair_id": "pairA", "grna_name": "g2"},
            {"bc_wells": "cell2", "pair_id": "pairB", "grna_name": "g3"},
        ])
        self.assertEqual(
            _production_assignment_sets(frame, dual_mode=True),
            {"cell1": {"pairA"}, "cell2": {"pairB"}},
        )
        self.assertEqual(
            _production_assignment_sets(frame.iloc[[0, 2]], dual_mode=False),
            {"cell1": {"g1"}, "cell2": {"g3"}},
        )

    def test_cross_source_overlap_distinguishes_ibar_pb_and_guide_only(self):
        columns = ["bc_wells", "bc_bci", "umi", "grna_seq", "ibar_seq_corrected"]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gex = root / "gex.tsv.gz"
            dedicated = root / "dedicated.tsv.gz"
            pd.DataFrame([
                ["A1_A1_A1", "1_1_1", "AAAAAAAAAA", "ACGT", "AAAAAAA"],
                ["A1_A1_A1", "2_2_2", "CCCCCCCCCC", "TGCA", ""],
            ], columns=columns).to_csv(gex, sep="\t", index=False, compression="gzip")
            pd.DataFrame([
                ["A1_A1_A1", "1_1_1", "AAAAAAAAAA", "ACGT", "AAAAAAA"],
                ["A1_A1_A1", "9_9_9", "CCCCCCCCCC", "TGCA", "CCCCCCC"],
            ], columns=columns).to_csv(
                dedicated, sep="\t", index=False, compression="gzip"
            )
            rows, summary = _cross_source_overlap(gex, dedicated, "C2")
            self.assertEqual(summary["exact_cross_source_guide_ibar_umi_overlap"], 1)
            self.assertEqual(summary["supporting_overlap_also_matching_pb"], 1)
            self.assertEqual(summary["guide_only_cross_source_overlap_without_ibar"], 1)
            self.assertEqual(set(rows["overlap_type"]), {
                "guide_ibar_umi", "guide_ibar_umi_pb", "guide_only_no_ibar"
            })

    def test_benchmark_and_paired_filter_end_to_end(self):
        with open(ROOT / "config.yaml") as handle:
            base = yaml.safe_load(handle)
        with open(ROOT / "examples/qc/qc_config_gex_grna_contamination.example.yaml") as handle:
            qc = yaml.safe_load(handle)["gex_grna_contamination_qc"]
        reference = pd.read_csv(ROOT / "references/grna_ibar_reference.csv")
        row = reference.iloc[0]
        variant = int(row["position"])
        structure = base["grna_structure"]
        positive = (
            structure["fixed_seq_1"]
            + row["grna_seq"]
            + structure[f"fixed_seq_{variant + 1}"]
            + row["ibar_seq"]
            + structure[f"fixed_seq_{variant + 3}"]
            + "ACGT" * 20
        )
        tso_negative = structure["fixed_seq_1"] + "ACGT" * 35
        no_tso = "GCTA" * 45

        with tempfile.TemporaryDirectory() as tmp:
            temp = Path(tmp)
            fastq_dir = temp / "fastq"
            fastq_dir.mkdir()
            r1 = fastq_dir / "lane_TEST_R1_001.fastq.gz"
            r2 = fastq_dir / "lane_TEST_R2_001.fastq.gz"
            records = [("positive", positive), ("tso_cdna", tso_negative), ("gex", no_tso)]
            write_fastq(r1, records)
            valid_r2 = (
                "ACGTACGTAC"
                + "ACACACAC"
                + "ATGAGGGGTCAG"
                + "GGGGGGGG"
                + "TCCAACCACCTC"
                + "AAAAAAAA"
            )
            write_fastq(r2, [(read_id, valid_r2) for read_id, _ in records])

            barcode_data = temp / "barcode_data.csv"
            pd.DataFrame(
                [
                    {"bci": 1, "sequence": "AAAAAAAA", "wind": 1, "bc_round": 1},
                    {"bci": 1, "sequence": "GGGGGGGG", "wind": 1, "bc_round": 2},
                    {"bci": 1, "sequence": "ACACACAC", "wind": 1, "bc_round": 3},
                ]
            ).to_csv(barcode_data, index=False)

            run_cfg = {
                "paths": {
                    "fastq_dir": str(fastq_dir),
                    "output_dir": str(temp / "splitpipe_output"),
                    "grna_reference": str(ROOT / "references/grna_ibar_reference.csv"),
                },
                "libraries": [{"sublibrary_name": "L1", "gex_prefix": "TEST"}],
                "grna_structure": structure,
                "gex_decontamination": copy.deepcopy(base["gex_decontamination"]),
                "gex_grna_contamination_qc": copy.deepcopy(qc),
            }
            run_cfg["gex_decontamination"]["tso_matching"]["backend"] = "myers"
            run_cfg["gex_grna_contamination_qc"]["output_dir"] = str(
                temp / "qc_serial"
            )
            run_cfg["gex_grna_contamination_qc"]["raw_benchmark"].update(
                {"max_reads": 0, "sample_fraction": 1.0, "read_table_max": 10}
            )
            run_cfg["gex_grna_contamination_qc"]["synthetic_fastq"].update(
                {"sample_fraction": 1.0, "max_reads": 0, "write_r2": True}
            )
            run_cfg["gex_grna_contamination_qc"]["bam_audit"].update(
                {"enabled": True, "required_for_audit": False}
            )
            run_cfg["gex_grna_contamination_qc"]["molecule_estimation"].update(
                {
                    "enabled": True,
                    "barcode_data_template": str(barcode_data),
                    "barcode_sample_start": 0,
                    "barcode_sample_end": 3,
                }
            )
            run_cfg["gex_grna_contamination_qc"]["alignment_benchmark"].update(
                {
                    "enabled": True,
                    "sample_fraction": 1.0,
                    "max_truth_reads": 0,
                }
            )
            run_cfg["gex_grna_contamination_qc"]["plots"]["enabled"] = False

            serial_summary = benchmark_library(
                run_cfg,
                str(temp / "config.yaml"),
                "L1",
                threads=1,
                execution_backend="serial",
                batch_size=1,
            )
            run_cfg["gex_grna_contamination_qc"]["output_dir"] = str(temp / "qc")
            summary = benchmark_library(
                run_cfg,
                str(temp / "config.yaml"),
                "L1",
                threads=4,
                execution_backend="parallel",
                batch_size=1,
            )
            self.assertEqual(summary["execution"]["workers"], 4)
            for key in (
                "raw_read_prevalence",
                "full_length_surrogate_truth",
                "short_read_hazard",
                "selected_filter",
                "recommended_filter_from_benchmark",
                "synthetic_reads_written",
            ):
                self.assertEqual(summary[key], serial_summary[key], key)
            for relative in (
                "read_evidence.tsv.gz",
                "short54_filter_metrics.csv",
                "short54_metrics_by_tso_offset.csv",
                "contaminant_guide_counts.csv",
                "contaminant_target_counts.csv",
                "contaminant_tso_offsets.csv",
                "per_cycle_base_composition.csv",
                "cell_read_contamination.csv",
                "synthetic_54/L1_all_R1_54.fastq.gz",
                "synthetic_54/L1_all_R2.fastq.gz",
                "synthetic_54/L1_clean_R1_54.fastq.gz",
                "synthetic_54/L1_clean_R2.fastq.gz",
                "synthetic_54/L1_rejected_R1_54.fastq.gz",
                "synthetic_54/L1_rejected_R2.fastq.gz",
                "alignment_input/full/L1_ground_truth_R1_full.fastq.gz",
                "alignment_input/short54/L1_ground_truth_R1_54.fastq.gz",
                "alignment_input/ground_truth_manifest.tsv.gz",
            ):
                serial_path = temp / "qc_serial" / "L1" / relative
                parallel_path = temp / "qc" / "L1" / relative
                if relative.endswith(".gz"):
                    with gzip.open(serial_path, "rt") as left, gzip.open(
                        parallel_path, "rt"
                    ) as right:
                        self.assertEqual(left.read(), right.read(), relative)
                else:
                    self.assertEqual(
                        serial_path.read_text(),
                        parallel_path.read_text(),
                        relative,
                    )
            self.assertEqual(
                artifact_semantic_digests(temp / "qc_serial", "L1"),
                artifact_semantic_digests(temp / "qc", "L1"),
            )
            self.assertEqual(
                benchmark_artifact_digests(temp / "qc_serial", "L1"),
                benchmark_artifact_digests(temp / "qc", "L1"),
            )
            aggregate_path = (
                temp / "qc_serial" / "L1" / "cell_read_contamination.csv"
            )
            aggregate_lines = aggregate_path.read_text().splitlines()
            aggregate_path.write_text(
                "\n".join([aggregate_lines[0], *reversed(aggregate_lines[1:])])
                + "\n"
            )
            self.assertEqual(
                benchmark_artifact_digests(temp / "qc_serial", "L1"),
                benchmark_artifact_digests(temp / "qc", "L1"),
            )
            self.assertEqual(summary["records_evaluated"], 3)
            self.assertEqual(
                summary["full_length_surrogate_truth"]["structure_contaminants"], 1
            )
            self.assertEqual(summary["synthetic_reads_written"], 3)

            synthetic = temp / "qc" / "L1" / "synthetic_54"
            self.assertEqual(count_fastq(synthetic / "L1_all_R1_54.fastq.gz"), 3)
            self.assertEqual(count_fastq(synthetic / "L1_rejected_R1_54.fastq.gz"), 1)
            self.assertEqual(count_fastq(synthetic / "L1_clean_R1_54.fastq.gz"), 2)

            filtered = filter_library(run_cfg, str(temp / "config.yaml"), "L1")
            self.assertEqual(filtered["total"], 3)
            self.assertEqual(filtered["rejected"], 1)
            self.assertEqual(count_fastq(filtered["clean_r1"]), 2)
            self.assertEqual(count_fastq(filtered["clean_r2"]), 2)
            self.assertEqual(count_fastq(filtered["rejected_r1"]), 1)
            self.assertEqual(count_fastq(filtered["rejected_r2"]), 1)
            for r1_record, r2_record in iter_fastq_pairs(
                filtered["clean_r1"], filtered["clean_r2"]
            ):
                self.assertEqual(len(r1_record.sequence), 54)
                self.assertEqual(len(r2_record.sequence), 58)

            truth_meta = summary["alignment_ground_truth_input"]
            truth_meta["reads_written_per_variant"] = 10
            with open(temp / "qc" / "L1" / "benchmark_summary.yaml", "w") as handle:
                yaml.safe_dump(summary, handle, sort_keys=False)
            for variant, false_reads, false_molecules in (
                ("full", 1, 1),
                ("short54", 3, 2),
            ):
                impact_dir = temp / "qc" / "L1" / "alignment_star" / variant
                impact_dir.mkdir(parents=True)
                impact = {
                    "library": "L1",
                    "variant": variant,
                    "truth_panel_token": truth_meta["truth_panel_token"],
                    "benchmark_qc_module_sha256": truth_meta[
                        "qc_module_sha256"
                    ],
                    "benchmark_config_sha256": truth_meta[
                        "benchmark_config_sha256"
                    ],
                    "alignment_audit_qc_module_sha256": truth_meta[
                        "qc_module_sha256"
                    ],
                    "input_reads": 10,
                    "mapped_primary_reads": 5,
                    "unique_reads": 4,
                    "multimapping_primary_reads": 1,
                    "false_gene_assigned_reads": false_reads,
                    "guide_target_assigned_reads": 0,
                    "false_gene_assignment_fraction_of_input_reads": false_reads / 10,
                    "false_gene_assignment_fraction_of_input_reads_wilson_95": [
                        0.0,
                        0.5,
                    ],
                    "guide_target_assignment_fraction_of_input_reads": 0.0,
                    "estimated_exact_key_molecules_hll": {
                        "input": 8.0,
                        "gene_assigned": float(false_molecules),
                    },
                    "false_gene_assignment_fraction_of_input_molecules": (
                        false_molecules / 8
                    ),
                    "splitpipe_equivalent_hamming1_molecule_impact": {
                        "enabled": True,
                        "ground_truth_input_molecules_by_guide": 8,
                        "false_gene_assigned_molecules": false_molecules,
                        "guide_target_molecules_any_read": 0,
                        "false_gene_assignment_fraction_of_input_molecules": (
                            false_molecules / 8
                        ),
                    },
                    "predicted_tso_trim_input_counts": {"True": 8, "False": 2},
                    "outcomes_by_predicted_tso_trim": {
                        "True": {"false_gene_assigned_reads": false_reads},
                        "False": {"false_gene_assigned_reads": 0},
                    },
                }
                with open(impact_dir / "alignment_impact_summary.yaml", "w") as handle:
                    yaml.safe_dump(impact, handle)
            aggregate(run_cfg, str(temp / "config.yaml"))
            self.assertTrue(
                (temp / "qc" / "all_libraries_contamination_summary.csv").exists()
            )
            self.assertTrue(
                (temp / "qc" / "all_libraries_short54_filter_metrics.csv").exists()
            )
            self.assertTrue(
                (temp / "qc" / "all_libraries_alignment_impact.csv").exists()
            )
            delta_path = (
                temp
                / "qc"
                / "all_libraries_full_vs_short54_alignment_delta.csv"
            )
            self.assertTrue(delta_path.exists())
            delta = pd.read_csv(delta_path).iloc[0]
            self.assertAlmostEqual(
                delta["read_false_assignment_absolute_delta"], 0.2
            )
            self.assertTrue(
                (temp / "qc" / "GEX_GRNA_CONTAMINATION_AUDIT_REPORT.md").exists()
            )
            self.assertTrue((temp / "qc" / "audit_completeness.csv").exists())
            completeness = pd.read_csv(
                temp / "qc" / "audit_completeness.csv"
            ).iloc[0]
            self.assertTrue(bool(completeness["bam_audit_expected"]))
            self.assertFalse(bool(completeness["bam_audit_required"]))
            self.assertFalse(bool(completeness["bam_audit_complete"]))
            self.assertTrue(bool(completeness["audit_complete"]))


if __name__ == "__main__":
    unittest.main()

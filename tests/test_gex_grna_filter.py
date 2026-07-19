#!/usr/bin/env python3
"""Focused tests for the removal-only packed Step 1a classifier."""

import csv
import random
import tempfile
import unittest
from dataclasses import FrozenInstanceError
from pathlib import Path
from unittest.mock import patch

import yaml

from scripts.gex_grna_filter import (
    FilterSettings,
    GexGrnaFilter,
    PackedGuideIndex,
    TsoMatcher,
    encode_dna,
    find_tso,
    load_filter_from_config,
    packed_hamming_distance,
)


TSO = "AAGCAGTGGTATCAACGCAGAGTGAATGGG"
GUIDE20 = "ACGTACGTACGTACGTACGT"
GUIDE21 = "TGCATGCATGCATGCATGCAT"


def mutate(sequence, position, replacement=None):
    if replacement is None:
        replacement = "A" if sequence[position] != "A" else "C"
    return sequence[:position] + replacement + sequence[position + 1 :]


def settings(lengths=(20, 21), remove_ambiguous=True):
    return FilterSettings(
        tso_sequence=TSO,
        tso_max_edit_distance=3,
        max_tso_start=2,
        guide_lengths=tuple(lengths),
        max_guide_mismatches=1,
        remove_ambiguous_hits=remove_ambiguous,
        max_seed_bucket_size=10_000,
        reference_path="fixture.csv",
    )


class PackedEncodingTest(unittest.TestCase):
    def test_two_bit_hamming_counts_bases_not_bits(self):
        reference = "A" * 20
        query = "C" + "A" * 18 + "T"
        self.assertEqual(
            packed_hamming_distance(encode_dna(reference), encode_dna(query), 20),
            2,
        )

    def test_every_single_substitution_is_retrieved_for_20_and_21mers(self):
        for guide in (GUIDE20, GUIDE21):
            index = PackedGuideIndex([guide], guide_lengths=(len(guide),))
            for position in range(len(guide)):
                match = index.match(mutate(guide, position))
                self.assertTrue(match.supported, msg=f"position {position}")
                self.assertEqual(match.distance, 1)
                self.assertEqual(len(match.sequence_group_ids), 1)

    def test_two_substitutions_do_not_pass_seed_candidate_verification(self):
        index = PackedGuideIndex([GUIDE20], guide_lengths=(20,))
        query = mutate(mutate(GUIDE20, 1), 17)
        match = index.match(query)
        self.assertFalse(match.supported)
        self.assertEqual(match.reason, "no_match")

    def test_exact_lookup_precedes_seed_search_and_collapses_duplicate_rows(self):
        index = PackedGuideIndex(
            [GUIDE20, GUIDE20, GUIDE21], guide_lengths=(20, 21)
        )
        exact = index.match(GUIDE20)
        self.assertTrue(exact.supported)
        self.assertEqual(exact.distance, 0)
        self.assertEqual(exact.candidate_count, 1)
        self.assertTrue(exact.shared_reference_sequence)
        self.assertFalse(exact.sequence_ambiguous)
        self.assertEqual(exact.reference_row_count, 2)
        self.assertEqual(index.stats.reference_rows, 3)
        self.assertEqual(index.stats.unique_sequence_groups, 2)
        self.assertEqual(index.stats.duplicate_reference_rows, 1)

    def test_complementary_masks_and_stats_are_immutable(self):
        index = PackedGuideIndex([GUIDE20, GUIDE21], guide_lengths=(20, 21))
        for length_stats in index.stats.by_length:
            first, second = length_stats.seed_masks
            full_mask = (1 << (2 * length_stats.guide_length)) - 1
            self.assertEqual(first & second, 0)
            self.assertEqual(first | second, full_mask)
            self.assertEqual(sum(length_stats.seed_position_counts), length_stats.guide_length)
        with self.assertRaises(FrozenInstanceError):
            index.stats.reference_rows = 10

    def test_tied_neighbor_reports_anonymous_sequence_ambiguity(self):
        guide_a = "A" * 20
        guide_c = "C" + "A" * 19
        query = "G" + "A" * 19
        index = PackedGuideIndex([guide_a, guide_c], guide_lengths=(20,))
        match = index.match(query)
        self.assertTrue(match.supported)
        self.assertEqual(match.distance, 1)
        self.assertTrue(match.sequence_ambiguous)
        self.assertEqual(len(match.sequence_group_ids), 2)
        self.assertEqual(match.reason, "tied_sequences_hamming_1")

    def test_reference_validation_rejects_lowercase_and_wrong_lengths(self):
        with self.assertRaisesRegex(ValueError, "uppercase A/C/G/T"):
            PackedGuideIndex([GUIDE20.lower()], guide_lengths=(20,))
        with self.assertRaisesRegex(ValueError, "not configured guide_lengths"):
            PackedGuideIndex([GUIDE21], guide_lengths=(20,))


class TsoDetectionTest(unittest.TestCase):
    def test_exact_tso_is_found_at_each_configured_offset(self):
        for offset in range(3):
            result = find_tso("C" * offset + TSO + GUIDE20, TSO, 3, 2)
            self.assertTrue(result.found)
            self.assertEqual((result.start, result.end), (offset, offset + len(TSO)))
            self.assertEqual(result.distance, 0)
            self.assertEqual(result.match_type, "exact")

    def test_tso_beyond_configured_offset_is_not_found(self):
        result = find_tso("CCC" + TSO + GUIDE20, TSO, 3, 2)
        self.assertFalse(result.found)

    def test_fuzzy_full_tso_supports_substitution_insertion_and_deletion(self):
        substituted = mutate(TSO, 12)
        inserted = TSO[:12] + "C" + TSO[12:]
        deleted = TSO[:12] + TSO[13:]
        cases = (
            (substituted, len(TSO)),
            (inserted, len(TSO) + 1),
            (deleted, len(TSO) - 1),
        )
        for observed_tso, observed_length in cases:
            result = find_tso("C" + observed_tso + GUIDE20, TSO, 3, 2)
            self.assertTrue(result.found)
            self.assertEqual(result.start, 1)
            self.assertEqual(result.end, 1 + observed_length)
            self.assertEqual(result.distance, 1)
            self.assertEqual(result.match_type, "fuzzy")

    def test_precomputed_myers_matcher_equals_reference(self):
        matcher = TsoMatcher(TSO, 3, 2, backend="myers")
        rng = random.Random(1402)
        alphabet = "ACGTN"
        cases = [
            "".join(rng.choice(alphabet) for _ in range(rng.randint(25, 80)))
            for _ in range(2000)
        ]
        cases.extend("C" * offset + TSO + GUIDE20 for offset in range(7))
        for sequence in cases:
            self.assertEqual(matcher.find(sequence), find_tso(sequence, TSO, 3, 2))

    def test_required_edlib_backend_fails_fast_when_unavailable(self):
        with patch("scripts.gex_grna_filter.HAS_EDLIB", False):
            with self.assertRaisesRegex(RuntimeError, "requires edlib"):
                TsoMatcher(TSO, 3, 2, backend="edlib")


class RemovalClassifierTest(unittest.TestCase):
    def test_exact_and_one_mismatch_reads_are_positive(self):
        classifier = GexGrnaFilter(
            PackedGuideIndex([GUIDE20, GUIDE21], guide_lengths=(20, 21)),
            settings(),
        )
        exact = classifier.classify("exact", "CC" + TSO + GUIDE21 + "GTT")
        fuzzy = classifier.classify(
            "fuzzy", TSO + mutate(GUIDE20, 17) + "GTT"
        )
        self.assertTrue(exact.filter_positive)
        self.assertEqual(exact.confidence_tier, "exact")
        self.assertEqual(exact.guide_length, 21)
        self.assertTrue(fuzzy.filter_positive)
        self.assertEqual(fuzzy.guide_distance, 1)
        self.assertFalse(fuzzy.sequence_ambiguous)

    def test_no_tso_short_guide_and_tso_cdna_are_negative(self):
        classifier = GexGrnaFilter(
            PackedGuideIndex([GUIDE20], guide_lengths=(20,)), settings((20,))
        )
        cases = (
            ("no_tso", "A" * 80, "no_tso"),
            ("short", TSO + GUIDE20[:19], "insufficient_full_guide"),
            ("cdna", TSO + "T" * 20, "no_reference_match"),
        )
        for read_id, sequence, reason in cases:
            decision = classifier.classify(read_id, sequence)
            self.assertFalse(decision.filter_positive)
            self.assertEqual(decision.filter_reason, reason)

    def test_fuzzy_tso_boundary_is_used_for_guide_extraction(self):
        classifier = GexGrnaFilter(
            PackedGuideIndex([GUIDE20], guide_lengths=(20,)), settings((20,))
        )
        inserted_tso = TSO[:12] + "C" + TSO[12:]
        decision = classifier.classify("insertion", inserted_tso + GUIDE20)
        self.assertTrue(decision.filter_positive)
        self.assertEqual(decision.tso_end, len(TSO) + 1)
        self.assertEqual(decision.guide_distance, 0)

    def test_ambiguity_removal_policy_changes_boolean_not_support(self):
        guide_a = "A" * 20
        guide_c = "C" + "A" * 19
        query = "G" + "A" * 19
        index = PackedGuideIndex([guide_a, guide_c], guide_lengths=(20,))
        remove = GexGrnaFilter(index, settings((20,), True)).classify(
            "remove", TSO + query
        )
        retain = GexGrnaFilter(index, settings((20,), False)).classify(
            "retain", TSO + query
        )
        self.assertTrue(remove.reference_supported)
        self.assertTrue(remove.filter_positive)
        self.assertTrue(remove.sequence_ambiguous)
        self.assertTrue(retain.reference_supported)
        self.assertFalse(retain.filter_positive)
        self.assertIn("not_removed", retain.filter_reason)

    def test_cross_length_results_do_not_depend_on_config_order(self):
        guide20 = "A" * 20
        guide21 = guide20 + "C"
        index = PackedGuideIndex([guide20, guide21], guide_lengths=(21, 20))
        classifier = GexGrnaFilter(index, settings((20, 21), True))
        decision = classifier.classify("both", TSO + guide21)
        self.assertTrue(decision.filter_positive)
        self.assertEqual(decision.guide_lengths, (20, 21))
        self.assertTrue(decision.sequence_ambiguous)

    def test_decision_has_no_biological_identity_fields(self):
        classifier = GexGrnaFilter(
            PackedGuideIndex([GUIDE20], guide_lengths=(20,)), settings((20,))
        )
        payload = classifier.classify("read", TSO + GUIDE20).to_dict()
        forbidden = {
            "grna_name",
            "target_gene",
            "ibar_seq",
            "pair_id",
            "umi",
            "position",
            "variant",
        }
        self.assertTrue(forbidden.isdisjoint(payload))


class ConfigLoadingTest(unittest.TestCase):
    def test_loads_relative_reference_and_comments(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference = root / "reference.csv"
            with reference.open("w", newline="") as handle:
                handle.write("# generated fixture\n")
                writer = csv.DictWriter(handle, fieldnames=("grna_name", "grna_seq"))
                writer.writeheader()
                writer.writerow({"grna_name": "not_exposed", "grna_seq": GUIDE20})
            config = {
                "paths": {"grna_reference": "reference.csv"},
                "grna_structure": {
                    "fixed_seq_1": TSO,
                    "grna_lengths": [20],
                    "tso_max_edit_dist": 2,
                },
                "gex_decontamination": {
                    "tso_max_edit_dist": 1,
                    "max_tso_start": 2,
                    "tso_matching": {"backend": "myers"},
                    "filter_matching": {
                        "guide_distance": "hamming",
                        "max_guide_mismatches": 1,
                        "allow_guide_indels": False,
                        "full_guide_required": True,
                        "backend": "packed_inline",
                        "index": "two_bit_two_seed",
                        "remove_reference_supported_ambiguous_hits": True,
                    },
                    "resource_limits": {"max_seed_bucket_size": 100},
                },
            }
            config_path = root / "config.yaml"
            with config_path.open("w") as handle:
                yaml.safe_dump(config, handle)

            classifier = load_filter_from_config(config_path)
            self.assertEqual(classifier.settings.reference_path, str(reference.resolve()))
            self.assertEqual(classifier.settings.tso_max_edit_distance, 1)
            self.assertTrue(classifier.classify("read", TSO + GUIDE20).filter_positive)

    def test_rejects_incompatible_packed_policy(self):
        config = {
            "paths": {"grna_reference": "unused.csv"},
            "grna_structure": {"fixed_seq_1": TSO, "grna_lengths": [20]},
            "gex_decontamination": {
                "filter_matching": {
                    "backend": "packed_inline",
                    "guide_distance": "levenshtein",
                }
            },
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "config.yaml"
            with path.open("w") as handle:
                yaml.safe_dump(config, handle)
            with self.assertRaisesRegex(ValueError, "guide_distance"):
                load_filter_from_config(path)


if __name__ == "__main__":
    unittest.main()

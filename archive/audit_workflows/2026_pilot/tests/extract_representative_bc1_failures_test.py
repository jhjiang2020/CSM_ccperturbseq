"""Archived, non-collected tests for the representative-read extractor."""

import unittest
import random

from scripts.extract_representative_bc1_failures import (
    KNOWN_SEQUENCE_ARTIFACTS,
    reconstruct_paired_fragment,
    reverse_complement,
    screen_sequence_artifacts,
)


class RepresentativeBc1FailureTests(unittest.TestCase):
    def test_example9_contains_reverse_tso_and_read2_adapter(self):
        sequence = (
            "ATTACAGCACTCTTCACAATGAGGGGTCAGTTCACGCATCCAACCACCTC"
            "CCATTCACTCTGCGTTGATACCACTGCTTAGATCGGAAGAGCGTCGTGTA"
        )
        result = screen_sequence_artifacts(sequence)
        hits = {
            (row["artifact"], row["orientation"]): row["length"]
            for row in result["supported_artifact_hits"]
        }
        self.assertEqual(hits[("Split_TSO", "reverse_complement")], 30)
        self.assertEqual(hits[("Illumina_Read2_adapter", "forward")], 21)

    def test_example1_has_no_supported_artifact_in_first_100_bases(self):
        sequence = (
            "ACTTGGAGGACTGGCATAATGCGGGGTCAGACATTGGCTCCAACATGCCA"
            "GAGGCTGTTCACCTTGGAGACCTGCTGCGGATATGGGTACGGCCCGGCGC"
        )
        result = screen_sequence_artifacts(sequence)
        self.assertEqual(result["supported_artifact_hits"], [])
        self.assertLess(
            result["best_artifact_overlap"]["length"],
            result["minimum_supported_exact_match"],
        )

    def test_reconstructs_short_insert_after_adapter_trimming(self):
        fragment = ("ACGT" * 25)
        r1 = fragment + KNOWN_SEQUENCE_ARTIFACTS["Illumina_Read1_adapter"]
        r2 = (
            reverse_complement(fragment)
            + KNOWN_SEQUENCE_ARTIFACTS["Illumina_Read2_adapter"]
        )
        result = reconstruct_paired_fragment(r1, r2)
        self.assertEqual(result["status"], "supported_overlap")
        self.assertEqual(result["r1_adapter_start"], len(fragment))
        self.assertEqual(result["r2_adapter_start"], len(fragment))
        self.assertEqual(result["estimated_fragment_length"], len(fragment))
        self.assertEqual(
            result["reconstructed_fragment_r1_orientation"], fragment
        )

    def test_reconstructs_ordinary_paired_overlap(self):
        fragment = "".join(random.Random(7).choices("ACGT", k=240))
        r1 = fragment[:150]
        r2 = reverse_complement(fragment[-150:])
        result = reconstruct_paired_fragment(r1, r2)
        self.assertEqual(result["status"], "supported_overlap")
        self.assertEqual(result["overlap_length"], 60)
        self.assertEqual(result["estimated_fragment_length"], len(fragment))
        self.assertEqual(
            result["reconstructed_fragment_r1_orientation"], fragment
        )


if __name__ == "__main__":
    unittest.main()

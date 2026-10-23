import csv
import gzip
from pathlib import Path
import unittest

from scripts.validate_asymmetric_pilot import fastq_records, fixed_observation, shorten_ibar


STRUCTURE = {
    'fixed_seq_1': 'AAGCAGTGGTATCAACGCAGAGTGAATGGG',
    'fixed_seq_2': 'GTTTGAGAGCTAAGCAGGA',
    'fixed_seq_3': 'GTTTCAGAGCTATGCTGGA',
    'fixed_seq_4': 'GACTGCTTAGCAAGTT',
    'fixed_seq_5': 'AACAGCATAGCAAGTT',
}


class PilotFixtureTest(unittest.TestCase):
    def test_fixed_coordinates_and_synchronized_quality_deletion(self):
        for position in (1, 2):
            for length in (20, 21):
                for offset in (0, 7):
                    sequence = ('N' * offset + STRUCTURE['fixed_seq_1'] + 'A' * length
                                + STRUCTURE[f'fixed_seq_{position + 1}'] + 'ACGTGCA'
                                + STRUCTURE[f'fixed_seq_{position + 3}'] + 'TTTT')
                    quality = ''.join(chr(33 + index % 60) for index in range(len(sequence)))
                    observation = fixed_observation(sequence, STRUCTURE, [20, 21])
                    self.assertEqual(observation['position'], position)
                    self.assertEqual(observation['ibar_seq'], 'ACGTGCA')
                    self.assertEqual(observation['tso_offset'], offset)
                    start = observation['ibar_start']
                    edited, edited_quality = shorten_ibar(sequence, quality, start, position)
                    if position == 2:
                        self.assertEqual(edited, sequence[:start + 4] + sequence[start + 7:])
                        self.assertEqual(edited_quality, quality[:start + 4] + quality[start + 7:])
                    else:
                        self.assertEqual((edited, edited_quality), (sequence, quality))

    def test_bad_layout_is_not_used_as_ground_truth(self):
        self.assertIsNone(fixed_observation('ACGT' * 40, STRUCTURE, [20, 21]))
        with self.assertRaisesRegex(ValueError, 'quality'):
            shorten_ibar('ACGT', 'III', 0, 2)
        with self.assertRaisesRegex(ValueError, 'outside'):
            shorten_ibar('ACGT', 'IIII', 3, 2)

    def test_published_fixture_changes_only_expected_ibar2_bases(self):
        fixture = Path(__file__).resolve().parent / 'fixtures/asymmetric_pilot_v2'
        control = list(fastq_records(fixture / 'control_R1_001.fastq.gz'))
        asymmetric = list(fastq_records(fixture / 'asymmetric_R1_001.fastq.gz'))
        with gzip.open(fixture / 'read_truth.tsv.gz', 'rt') as handle:
            truth = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(len(control), 1920)
        self.assertEqual(len(asymmetric), len(control))
        self.assertEqual(len(truth), len(control))
        for original, edited, expected in zip(control, asymmetric, truth):
            self.assertEqual(original[1], edited[1])
            self.assertEqual(edited[1], '@' + expected['read_id'])
            observation = fixed_observation(original[2], STRUCTURE, [20, 21])
            self.assertEqual(observation['position'], int(expected['position']))
            self.assertEqual(observation['grna_seq'], expected['grna_seq'])
            self.assertEqual(observation['ibar_seq'], expected['original_ibar'])
            self.assertEqual(shorten_ibar(original[2], original[3], observation['ibar_start'], observation['position']), edited[2:])


if __name__ == '__main__':
    unittest.main()

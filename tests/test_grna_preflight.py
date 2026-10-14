import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from scripts.grna_preflight import check_ibar_reference, preflight_options


class IbarPreflightTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.config = dict(ibar_length=4, tso_seq='', stem1='', stem2=None,
                           grna_lengths=[20], tracr1='', tracr2=None,
                           anchor_max_mm=2, tso_max_dist=3)

    def check(self, observations, references=(('G1', 'AACG'),), limit=10000,
              compressed=False, extractor=None):
        path = self.root / ('input.tsv.gz' if compressed else 'input.tsv')
        opener = gzip.open if compressed else open
        with opener(path, 'wt') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['r1_seq'])
            writer.writerows([[value] for value in observations])
        def extract(value, *args):
            guide, ibar = value.split(':')
            return {'grna_seq': guide, 'ibar_seq': ibar}
        return check_ibar_reference(path, references, self.config,
                                    extractor or extract, sample_reads=limit)

    def codes(self, result):
        return {warning['code'] for warning in result['warnings']}

    def test_correct_orientation_has_no_warning(self):
        result = self.check(['G1:AACG'] * 150)
        self.assertEqual(result['status'], 'ok')
        self.assertFalse(result['warnings'])
        self.assertEqual(result['counts']['same_guide_stored_exact'], 150)

    def test_reversed_orientation_warns_without_changing_reference(self):
        reference = [('G1', 'AACG')]
        result = self.check(['G1:CGTT'] * 150, reference, compressed=True)
        self.assertIn('possible_reverse_complement', self.codes(result))
        self.assertEqual(result['counts']['same_guide_rc_exact'], 150)
        self.assertEqual(reference, [('G1', 'AACG')])

    def test_same_guide_evidence_handles_rc_closed_global_reference(self):
        result = self.check(['G1:CGTT'] * 150,
                            [('G1', 'AACG'), ('G2', 'CGTT')])
        self.assertEqual(result['counts']['global_stored_only'], 0)
        self.assertEqual(result['orientation_evidence'], ['same_guide'])
        self.assertIn('possible_reverse_complement', self.codes(result))

    def test_palindromes_and_sparse_samples_are_inconclusive(self):
        palindrome = self.check(['G1:ATAT'] * 150, [('G1', 'ATAT')])
        sparse = self.check(['G1:CGTT'] * 30)
        for result in (palindrome, sparse):
            self.assertEqual(result['status'], 'inconclusive')
            self.assertFalse(result['warnings'])

    def test_shared_guides_and_multiple_ibars_do_not_double_count(self):
        result = self.check(['G1:AACG'] * 150,
                            [('G1', 'AACG'), ('G1', 'AACG'), ('G1', 'CCCC')])
        self.assertEqual(result['counts']['same_guide_stored_exact'], 150)
        self.assertFalse(result['warnings'])

    def test_length_low_support_and_structure_diagnostics(self):
        self.config['ibar_length'] = 5
        result = self.check(['G1:GGGG'] * 150)
        self.assertIn('ibar_length_mismatch', self.codes(result))
        self.assertIn('low_ibar_support', self.codes(result))
        result = self.check(['G1:AACG'] * 150, extractor=lambda *args: None)
        self.assertIn('low_extraction_support', self.codes(result))
        self.assertNotIn('possible_reverse_complement', self.codes(result))

    def test_sample_is_bounded_and_no_ibar_skips_input(self):
        result = self.check(['G1:CGTT'] * 150, limit=25)
        self.assertEqual(result['counts']['sampled_reads'], 25)
        self.config['ibar_length'] = 0
        result = check_ibar_reference(self.root/'missing.tsv', [], self.config,
                                      lambda *args: self.fail('must not extract'))
        self.assertEqual(result['status'], 'skipped')

    def test_invalid_reference_is_reported(self):
        result = self.check(['G1:AACG'] * 150,
                            [('G1', None), ('G1', 'acgt'), ('G1', 'AACG')])
        self.assertIn('invalid_reference_ibar', self.codes(result))

    def test_options_preserve_old_configs_and_validate_new_values(self):
        self.assertEqual(preflight_options({}), (True, 10000))
        self.assertEqual(preflight_options({'preflight': {'enabled': False}}), (False, 10000))
        for invalid in (None, True, {'enabled': 'false'}, {'sample_reads': True},
                        {'sample_reads': 0}, {'sample_reads': 100001}):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                preflight_options({'preflight': invalid})


if __name__ == '__main__':
    unittest.main()

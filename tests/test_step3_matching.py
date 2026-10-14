import csv
import gzip
import importlib.util
import itertools
import random
import subprocess
import sys
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts/03_match_grna_ibar.py'
spec = importlib.util.spec_from_file_location('step3_under_test', SCRIPT)
step3 = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = step3
spec.loader.exec_module(step3)


class MatchingScalabilityTest(unittest.TestCase):
    def test_fuzzy_matching_agrees_with_exhaustive_distance_search(self):
        rng = random.Random(41)
        sequences = [''.join(x) for x in itertools.product('ACGT', repeat=4)]
        references = rng.sample(sequences, 24)
        references += [references[0], references[1][:-1], references[2]+'A']
        names = [f'g{i}' for i in range(len(references))]
        queries = sequences + [seq[:3] for seq in sequences[::4]] + [seq+'A' for seq in sequences[::4]]
        for max_dist in range(4):
            matcher = step3.NeighborhoodMatcher(references, names, max_dist)
            for query in queries:
                distances = [step3.levenshtein_distance(query, seq) for seq in references]
                best = min(distances)
                indices = [i for i,distance in enumerate(distances) if distance == best]
                actual = matcher.match(query)
                if best > max_dist:
                    expected = (None, None, False)
                elif len({references[i] for i in indices}) > 1:
                    expected = (None, best, False)
                elif len(indices) == 1:
                    expected = (names[indices[0]], best, True)
                else:
                    expected = (sorted(names[i] for i in indices), best, False)
                    actual = (sorted(actual[0]), actual[1], actual[2])
                self.assertEqual(actual, expected, (query, max_dist))

    def test_work_and_completed_results_stay_bounded(self):
        consumed = []
        def items():
            for item in range(20):
                consumed.append(item)
                yield item
        with ThreadPoolExecutor(2) as executor:
            output = []
            for value in step3._bounded_results(executor, abs, items(), max_pending=3):
                output.append(value)
                self.assertLessEqual(len(consumed) - len(output), 3)
        self.assertEqual(sorted(output), list(range(20)))

    def test_worker_errors_are_not_swallowed(self):
        with ThreadPoolExecutor(1) as executor, self.assertRaises(ValueError):
            list(step3._bounded_results(executor, int, ['bad'], 2))


class Step3CommandTest(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.guide = 'ACGTGCACTGATCGTACGTA'
        tso = 'AAGCAGTGGTATCAACGCAGAGTGAATGGG'
        stem = 'GTTTGAGAGCTAGAAATAG'
        tracr = 'GACTGCTT'
        self.ref = self.root / 'reference.csv'
        with self.ref.open('w') as handle:
            writer = csv.writer(handle)
            writer.writerow(['grna_name','grna_seq','ibar_seq','target_gene','position'])
            writer.writerow(['guide1', self.guide, 'AACG', 'GENE1', 1])
        self.barcode = self.root / 'barcodes.tsv.gz'
        with gzip.open(self.barcode, 'wt') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['read_id','bc_wells','bc_bci','umi','r1_seq','r1_qual'])
            for i in range(120):
                r1 = 'GATACGT'+tso+self.guide+stem+'CGTT'+tracr
                writer.writerow([f'r{i}', 'A1_B1_C1', '001', 'ACGTACGT', r1, 'F'*len(r1)])
        with gzip.open(self.root/'run_guide1_S1_R1_001.fastq.gz', 'wt') as handle:
            for i in range(120):
                handle.write(f'@r{i}\n{r1}\n+\n'+ 'F'*len(r1)+'\n')
        self.config = dict(paths={'grna_reference':str(self.ref)},
            grna_structure=dict(single_guide_mode=True, fixed_seq_1=tso,
                fixed_seq_2=stem, fixed_seq_4=tracr, grna_lengths=[20],
                ibar_length=4, fixed_seq_max_mismatches=0, tso_max_edit_dist=0),
            grna_matching=dict(max_lev_dist=1, max_lev_dist_ibar=1, save_failed_reads=True),
            barcode_extraction={'n_processes':1}, libraries=[{'sublibrary_name':'L1'}])
        self.config['paths']['fastq_dir'] = str(self.root)
        self.config['libraries'][0]['grna_prefix'] = '_guide1_'

    def run_command(self, name, workers=1, extra=()):
        output = self.root/name
        self.config['paths']['output_dir'] = str(output)
        path = self.root/f'{name}.yaml'
        path.write_text(yaml.safe_dump(self.config))
        barcode_args = [] if '--preflight-fastq' in extra else ['--barcode-input', str(self.barcode)]
        result = subprocess.run([sys.executable, '-B', str(SCRIPT), '--config',str(path),
            '--library','L1',*barcode_args,'--n_processes',str(workers),
            '--chunk_size','17',*extra], capture_output=True, text=True, check=True)
        return output, result.stdout

    def test_preflight_only_warns_and_never_opens_production_outputs(self):
        output, log = self.run_command('preflight', extra=['--preflight-only'])
        self.assertIn('WARNING [possible_reverse_complement]', log)
        self.assertTrue((output/'grna_matches/L1_preflight.yaml').exists())
        self.assertFalse((output/'grna_matches/L1_grna_matches.tsv.gz').exists())
        self.assertFalse((output/'grna_matches/L1_matching_stats.yaml').exists())

    def test_automatic_warning_precedes_matching_and_parallel_calls_agree(self):
        output1, log1 = self.run_command('single', workers=1)
        output2, log2 = self.run_command('parallel', workers=2)
        self.assertLess(log2.index('WARNING [possible_reverse_complement]'), log2.index('Matching: 2 worker'))
        self.assertIn('Chunks 1:', log2)
        def read(output):
            with gzip.open(output/'grna_matches/L1_grna_matches.tsv.gz','rt') as handle:
                return sorted(csv.DictReader(handle,delimiter='\t'),key=lambda r:r['read_id'])
        rows1, rows2 = read(output1), read(output2)
        self.assertEqual(rows1, rows2)
        self.assertEqual(len(rows1), 120)
        self.assertEqual(rows1[0]['bc_bci'], '001')
        self.assertEqual(self.ref.read_text().count('AACG'), 1)

    def test_explicit_preflight_overrides_disabled_automatic_check(self):
        self.config['grna_matching']['preflight'] = {'enabled':False}
        output, log = self.run_command('disabled')
        self.assertIn('disabled by config', log)
        self.assertFalse((output/'grna_matches/L1_preflight.yaml').exists())
        _, log = self.run_command('forced', extra=['--preflight-only'])
        self.assertIn('WARNING [possible_reverse_complement]', log)

    def test_raw_fastq_preflight_works_before_barcode_extraction(self):
        self.barcode.unlink()
        output, log = self.run_command('raw', extra=['--preflight-only','--preflight-fastq'])
        report = yaml.safe_load((output/'grna_matches/L1_preflight.yaml').read_text())
        self.assertEqual(report['input_type'], 'raw_fastq')
        self.assertEqual(report['counts']['sampled_reads'], 120)
        self.assertIn('WARNING [possible_reverse_complement]', log)
        self.assertFalse((output/'grna_matches/L1_grna_matches.tsv.gz').exists())


if __name__ == '__main__':
    unittest.main()

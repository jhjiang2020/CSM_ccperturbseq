"""Sparse recovery annotations must not change parsing or molecule identity."""
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

PATH = Path(__file__).resolve().parents[1] / 'scripts/04_assign_grna_to_cells.py'
spec = importlib.util.spec_from_file_location('step4_stream', PATH)
step4 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(step4)


class MatchStreamTest(unittest.TestCase):
    def test_sparse_recovery_text_across_internal_and_outer_chunks(self):
        # The annotation is beyond the selected-column count, as in Step 3.
        # Crossing the C parser's internal buffer triggered pandas IndexError.
        groups = ['bc_wells', 'grna_name', 'target_gene', 'guide_unit_id', 'position']
        columns = groups + ['umi', 'read_id'] + [f'unused{i}' for i in range(40)] + ['structure_correction']
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'matches.tsv'
            with path.open('w') as handle:
                writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t')
                writer.writeheader()
                for i in range(20000):
                    writer.writerow(dict(bc_wells='cell', grna_name='guide', target_gene='GENE',
                        guide_unit_id='unit', position=1, umi='ACGTACGT', read_id=f'r{i}',
                        structure_correction='' if i < 12000 else json.dumps({'length':'shortened_only','policy':'frozen'})))
            for size in (500000, 10000):
                for multi in (False, True):
                    with self.subTest(chunk_size=size, multi=multi):
                        if multi:
                            frame, stats = step4.stream_library_umi_dedup_sources(
                                [('dedicated',path)], 'L1', 'sample', 0, groups, chunk_size=size)
                        else:
                            frame, stats = step4.stream_library_umi_dedup(
                                path, 'L1', 'sample', 0, groups, chunk_size=size)
                        self.assertEqual(len(frame), 1)
                        self.assertEqual(int(frame.iloc[0].umi_count), 1)
                        self.assertEqual(int(frame.iloc[0].read_count), 20000)
                        detail = json.loads(frame.iloc[0].structure_correction)
                        self.assertEqual(detail['recovery_umi_mixed'], 1)
                        self.assertEqual(detail['recovery_umi_rescue_only'], 0)


if __name__ == '__main__':
    unittest.main()

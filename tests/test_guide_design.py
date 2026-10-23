import copy
import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pandas as pd
import yaml

from scripts.design_preflight import check_design
from scripts.guide_design import (
    DesignMatcher, assign_cells, compile_design, extract_unit,
    legacy_model_eligible, linked_pair_status, reverse_complement,
)
from scripts.guide_design_qc import cell_design_qc, run_guideseq
from scripts.gex_grna_filter import settings_from_config
from scripts.qc_config_loader import load_guideseq_config


ROOT = Path(__file__).resolve().parents[1]


def load_script(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


STEP3 = load_script('design_step3', '03_match_grna_ibar.py')
STEP4 = load_script('design_step4', '04_assign_grna_to_cells.py')
STEP1C = load_script('design_step1c', '01c_process_gex_grna.py')
STEP5 = load_script('design_step5', '05_build_anndata.py')
Q1 = load_script('design_q1', 'qc_q1_guideseq.py')


def configuration():
    return {
        'guide_design': {'enabled': True},
        'grna_structure': {
            'fixed_seq_1': 'AAGCAGTGGTATCAACGCAGAGTGAATGGG',
            'fixed_seq_2': 'GTTTGAGAGCTAAGCAGGA',
            'fixed_seq_3': 'GTTTCAGAGCTATGCTGGA',
            'fixed_seq_4': 'GACTGCTTAGCAAGTT',
            'fixed_seq_5': 'AACAGCATAGCAAGTT',
            'single_guide_mode': False, 'grna_lengths': [19, 20],
            'ibar_lengths_by_position': 'auto',
            'fixed_seq_max_mismatches': 2, 'tso_max_edit_dist': 0,
        },
        'grna_matching': {'max_lev_dist': 1, 'max_lev_dist_ibar': 1, 'save_failed_reads': True},
        'guide_assignment': {'dual_guide_counting_mode': True, 'min_umi_threshold': 1,
                             'umi_dedup_ibar': True, 'max_umi_edits': 1,
                             'compute_saturation_curve': False},
        'barcode_extraction': {'n_processes': 1},
        'libraries': [{'sublibrary_name': 'L1', 'sample_name': 'sample'}],
    }


def reference(second_ibar='C', first_ibar='ACCTGTA'):
    return pd.DataFrame([
        dict(grna_name='g1', grna_seq='ACGTACGTACGTACGTACG', ibar_seq=first_ibar,
             target_gene='GENE1', position=1, pair_id='pair1', target_combo='GENE1_GENE2', stem_variant='A', target_version='v1'),
        dict(grna_name='g2', grna_seq='TGCATGCATGCATGCATGCA', ibar_seq=second_ibar,
             target_gene='GENE2', position=2, pair_id='pair1', target_combo='GENE1_GENE2', stem_variant='B', target_version='v1'),
    ])


def sequence(row, config, ibar=None, scaffold=None):
    variant = scaffold or row['position']
    structure = config['grna_structure']
    return (structure['fixed_seq_1'] + row['grna_seq'] + structure[f'fixed_seq_{variant + 1}']
            + (row['ibar_seq'] if ibar is None else ibar) + structure[f'fixed_seq_{variant + 3}'])


def observed(row, contract, config, **kwargs):
    raw = sequence(row, config, **kwargs)
    return extract_unit(raw, len(config['grna_structure']['fixed_seq_1']), contract)


class DesignContractTest(unittest.TestCase):
    def test_model_capability_depends_on_design_identifiability(self):
        config = configuration()
        self.assertTrue(compile_design(reference('TGACTAG'), config)['capabilities']['legacy_recombination_model'])
        self.assertFalse(compile_design(reference('C'), config)['capabilities']['legacy_recombination_model'])
        self.assertFalse(compile_design(reference(''), config)['capabilities']['legacy_recombination_model'])
        self.assertFalse(legacy_model_eligible(pd.DataFrame([{'guide_unit_id': 'unit', 'legacy_model_eligible': False}])))

    def test_length_identifies_position_with_identical_scaffolds(self):
        config = configuration()
        config['grna_structure']['fixed_seq_3'] = config['grna_structure']['fixed_seq_2']
        config['grna_structure']['fixed_seq_5'] = config['grna_structure']['fixed_seq_4']
        table = reference('')
        contract = compile_design(table, config)
        extraction = observed(table.iloc[1], contract, config)
        self.assertEqual(extraction['position'], 2)
        self.assertEqual(extraction['position_from_scaffold'], 0)
        self.assertEqual(extraction['structure_status'], 'concordant')

    def test_fingerprint_covers_annotations_and_matching_policy(self):
        config = configuration()
        table = reference()
        expected = compile_design(table, config)['fingerprint']
        self.assertEqual(compile_design(table.iloc[::-1], config)['fingerprint'], expected)
        config['guide_assignment']['infer_missing_partner'] = False
        self.assertEqual(compile_design(table, config)['fingerprint'], expected)
        config['grna_matching']['max_lev_dist'] = 2
        self.assertNotEqual(compile_design(table, config)['fingerprint'], expected)
        config = configuration()
        table['target_combo'] = 'new_label'
        self.assertNotEqual(compile_design(table, config)['fingerprint'], expected)

    def test_close_codewords_cannot_be_fuzzy_corrected(self):
        table = reference()
        other = table.copy()
        other['pair_id'] = 'pair2'
        other.loc[other.position == 1, 'ibar_seq'] = 'ACCTGTT'
        contract = compile_design(pd.concat([table, other], ignore_index=True), configuration())
        self.assertTrue(all(distance == 0 for distance in contract['ibar_correction'].values()))

    def test_partner_inference_can_be_disabled_without_losing_guide_counts(self):
        config = configuration()
        table = reference('')
        contract = compile_design(table, config)
        match = DesignMatcher(contract).resolve(observed(table.iloc[0], contract, config))
        match.update(bc_wells='cell', umi_count=4)
        calls = assign_cells(pd.DataFrame([match]), table, 1, infer_missing_partner=False)
        self.assertEqual(calls['pair_umi_total'].sum(), 4)
        self.assertEqual(calls.iloc[0]['inferred_pos2_grna'], '')

    def test_invalid_bases_and_incomplete_pairs_are_rejected(self):
        table = reference()
        table.loc[0, 'ibar_seq'] = 'N'
        with self.assertRaisesRegex(ValueError, 'A/C/G/T'):
            compile_design(table, configuration())
        with self.assertRaisesRegex(ValueError, 'both positions'):
            compile_design(reference().iloc[:1], configuration())

    def test_supported_layouts_and_absence(self):
        config = configuration()
        for first, second in [('ACCTGTA', 'TGACTAG'), ('ACCTGTA', 'C'), ('ACCTGTA', ''), ('', 'ACCTGTA'), ('', '')]:
            with self.subTest(lengths=(len(first), len(second))):
                table = reference(second, first)
                contract = compile_design(table, config)
                matcher = DesignMatcher(contract)
                for row in table.to_dict('records'):
                    extraction = observed(row, contract, config)
                    self.assertEqual(extraction['grna_seq'], row['grna_seq'])
                    self.assertEqual(extraction['ibar_seq'], row['ibar_seq'])
                    match = matcher.resolve(extraction)
                    self.assertEqual(match['pair_id'], 'pair1')
                    self.assertEqual(match['ibar_status'], 'present_resolved' if row['ibar_seq'] else 'absent_by_design')
                    self.assertEqual(match['position'], row['position'])

    def test_short_ibar_correction_and_explicit_conflicts(self):
        config = configuration()
        config['grna_structure']['ibar_lengths_by_position'] = {'pos1': [7], 'pos2': [7]}
        contract = compile_design(reference(), config)
        self.assertTrue(any(warning['code'] == 'ibar_length_mismatch' for warning in contract['warnings']))
        second_id = contract['rows'][1]['guide_unit_id']
        self.assertEqual(contract['ibar_correction'][second_id], 0)
        config['guide_design']['strict'] = True
        with self.assertRaisesRegex(ValueError, 'ibar_length_mismatch'):
            compile_design(reference(), config)

    def test_ambiguous_units_and_joint_pair_resolution(self):
        config = configuration()
        table = reference('')
        repeated = table.copy()
        repeated['pair_id'] = 'pair2'
        repeated.loc[repeated.position == 1, 'ibar_seq'] = 'TGAACCT'
        table = pd.concat([table, repeated], ignore_index=True)
        contract = compile_design(table, config)
        matcher = DesignMatcher(contract)
        first = matcher.resolve(observed(table.iloc[0], contract, config))
        second = matcher.resolve(observed(table.iloc[1], contract, config))
        self.assertEqual(first['pair_candidate_count'], 1)
        self.assertEqual(second['pair_candidate_count'], 2)
        self.assertEqual(second['pair_id'], '')
        self.assertEqual(linked_pair_status(first, second)['pair_candidates'], ['pair1'])
        second.update(bc_wells='cell', umi_count=3)
        calls = assign_cells(pd.DataFrame([second]), table, 1)
        self.assertEqual(calls.iloc[0]['detection_status'], 'ambiguous')
        first.update(bc_wells='cell', umi_count=2)
        calls = assign_cells(pd.DataFrame([first, second]), table, 1)
        self.assertEqual(calls.iloc[0]['detection_status'], 'complete')
        self.assertEqual(calls.iloc[0]['pair_candidate_count'], 1)
        self.assertEqual(calls['pair_umi_total'].sum(), 5)

    def test_conflict_is_retained_and_truncation_is_not_absence(self):
        config = configuration()
        table = reference('')
        contract = compile_design(table, config)
        row = table.iloc[0]
        extraction = observed(row, contract, config, scaffold=2)
        self.assertEqual(extraction['position_from_length'], 1)
        self.assertEqual(extraction['position_from_scaffold'], 2)
        self.assertEqual(extraction['structure_status'], 'conflict')
        match = DesignMatcher(contract).resolve(extraction)
        self.assertEqual(match['pair_id'], '')
        raw = sequence(table.iloc[1], config)[:-10]
        self.assertIsNone(extract_unit(raw, len(config['grna_structure']['fixed_seq_1']), contract))

    def test_ibar_indel_is_not_reassigned_to_the_other_position(self):
        config = configuration()
        table = reference('TGACTA')
        contract = compile_design(table, config)
        extraction = observed(table.iloc[0], contract, config, ibar='ACCTGT')
        self.assertEqual(extraction['position_from_length'], 2)
        self.assertEqual(extraction['position_from_scaffold'], 1)
        self.assertEqual(extraction['structure_status'], 'conflict')
        match = DesignMatcher(contract).resolve(extraction)
        self.assertTrue(match is None or not match['pair_id'])

    def test_linked_recombination_without_ibar_pair_uniqueness(self):
        config = configuration()
        first_pair = reference('', '')
        second_pair = reference('', '')
        second_pair['pair_id'] = 'pair2'
        second_pair['grna_seq'] = ['TTTTACGTACGTACGTACG', 'CCCCCGCATGCATGCATGCA']
        table = pd.concat([first_pair, second_pair], ignore_index=True)
        contract = compile_design(table, config)
        self.assertEqual(contract['capabilities']['ibar_tuple']['unique'], 0)
        self.assertEqual(contract['capabilities']['spacer_tuple']['unique'], 2)
        matcher = DesignMatcher(contract)
        first = matcher.resolve(observed(table.iloc[0], contract, config))
        second = matcher.resolve(observed(table.iloc[3], contract, config))
        self.assertEqual(linked_pair_status(first, second)['status'], 'recombination_compatible')

    def test_foreign_known_ibar_remains_decoupling_evidence(self):
        config = configuration()
        table = reference()
        other = table.copy()
        other['pair_id'] = 'pair2'
        other['grna_name'] = ['g3', 'g4']
        other.loc[other.position == 1, 'grna_seq'] = 'TTTTACGTACGTACGTACG'
        other.loc[other.position == 1, 'ibar_seq'] = 'TGAACCT'
        table = pd.concat([table, other], ignore_index=True)
        contract = compile_design(table, config)
        matcher = DesignMatcher(contract)
        first = matcher.resolve(observed(table.iloc[0], contract, config, ibar='TGAACCT'))
        second = matcher.resolve(observed(table.iloc[1], contract, config))
        self.assertTrue(first['ibar_decoupled'])
        self.assertTrue(first['ibar_decoupling_eligible'])
        self.assertEqual(first['pair_id'], '')
        self.assertEqual(linked_pair_status(first, second)['status'], 'recombination_compatible')
        self.assertEqual(contract['capabilities']['per_position'][2]['decoupling_guide_groups_eligible'], 0)

    def test_shared_evidence_never_inflates_counts(self):
        config = configuration()
        table = reference('')
        other = table.copy()
        other['pair_id'] = 'pair2'
        other.loc[other.position == 1, 'ibar_seq'] = 'TGAACCT'
        table = pd.concat([table, other], ignore_index=True)
        contract = compile_design(table, config)
        matcher = DesignMatcher(contract)
        observations = []
        for index in (0, 1, 2):
            match = matcher.resolve(observed(table.iloc[index], contract, config))
            match.update(bc_wells='cell', umi_count=3)
            observations.append(match)
        calls = assign_cells(pd.DataFrame(observations), table, 1)
        self.assertEqual(calls['pair_umi_total'].sum(), 9)
        self.assertNotIn('complete', set(calls['detection_status']))

    def test_duplicate_full_signatures_stay_ambiguous(self):
        table = reference()
        other = table.copy()
        other['pair_id'] = 'pair2'
        table = pd.concat([table, other], ignore_index=True)
        config = configuration()
        contract = compile_design(table, config)
        matcher = DesignMatcher(contract)
        first, second = [matcher.resolve(observed(row, contract, config)) for row in table.iloc[:2].to_dict('records')]
        self.assertEqual(linked_pair_status(first, second)['status'], 'designed_ambiguous')

    def test_position_ambiguous_guides_retain_molecules(self):
        config = configuration()
        config['grna_structure']['fixed_seq_3'] = config['grna_structure']['fixed_seq_2']
        config['grna_structure']['fixed_seq_5'] = config['grna_structure']['fixed_seq_4']
        table = reference('', '')
        table.loc[1, 'grna_seq'] = table.iloc[0]['grna_seq']
        contract = compile_design(table, config)
        match = DesignMatcher(contract).resolve(observed(table.iloc[0], contract, config))
        self.assertEqual(match['position'], 0)
        self.assertEqual(match['ibar_status'], 'absent_by_design')
        match.update(bc_wells='cell', umi_count=3)
        calls = assign_cells(pd.DataFrame([match]), table, 1)
        self.assertEqual(calls['unpositioned_umi'].sum(), 3)
        self.assertEqual(calls.iloc[0]['pair_id'], '')


class DesignWorkflowTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.config = configuration()

    def write_reads(self, table, reverse=False, per_unit=120):
        path = self.root / 'barcodes.tsv.gz'
        with gzip.open(path, 'wt') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['read_id', 'bc_wells', 'bc_bci', 'umi', 'r1_seq'])
            for row in table.to_dict('records'):
                for index in range(per_unit):
                    ibar = reverse_complement(row['ibar_seq']) if reverse else row['ibar_seq']
                    writer.writerow([f"r{row['position']}_{index}", 'A1_B1_C1', '1_1_1', 'ACGTACGT', sequence(row, self.config, ibar=ibar)])
        return path

    def test_pilot_discovers_lengths_orientation_and_zero_ibar(self):
        table = reference('')
        self.config['guide_design']['ibar_orientation'] = 'auto'
        contract = compile_design(table, self.config)
        path = self.write_reads(table, reverse=True)
        result, effective = check_design(path, contract, self.config, STEP3.find_tso)
        self.assertEqual(result['orientation'][1], 'reverse_complement')
        self.assertEqual(result['orientation'][2], 'not_applicable')
        self.assertEqual(effective['rows'][0]['ibar_seq'], reverse_complement(table.iloc[0]['ibar_seq']))
        self.assertEqual(effective['input_fingerprint'], contract['fingerprint'])
        self.assertEqual(result['observed_ibar_length_counts'][2], {0: 120})

    def test_cli_matching_and_counting_for_asymmetric_and_no_ibar(self):
        for first, second in [('ACCTGTA', 'C'), ('ACCTGTA', ''), ('', '')]:
            with self.subTest(lengths=(len(first), len(second))):
                table = reference(second, first)
                reference_path = self.root / 'reference.csv'
                table.to_csv(reference_path, index=False)
                reads = self.write_reads(table, per_unit=3)
                output = self.root / f'output_{len(first)}_{len(second)}'
                self.config['paths'] = {'grna_reference': str(reference_path), 'output_dir': str(output)}
                config_path = self.root / 'config.yaml'
                config_path.write_text(yaml.safe_dump(self.config))
                subprocess.run([sys.executable, '-B', str(ROOT / 'scripts/03_match_grna_ibar.py'),
                                '--config', str(config_path), '--barcode-input', str(reads),
                                '--n_processes', '2', '--chunk_size', '2'], check=True, capture_output=True, text=True)
                matched = output / 'grna_matches/L1_grna_matches.tsv.gz'
                rows = pd.read_csv(matched, sep='\t', keep_default_na=False)
                self.assertEqual(len(rows), 6)
                columns = STEP4._determine_group_cols(matched, True)
                counted, _ = STEP4.stream_library_umi_dedup(matched, 'L1', 'sample', None, columns)
                self.assertEqual(counted['umi_count'].sum(), 2)
                annotated, calls = STEP4._annotate_and_assign(counted, table, [], columns, self.config['guide_assignment'])
                self.assertEqual(len(annotated), 2)
                self.assertEqual(calls.iloc[0]['detection_status'], 'complete')
                self.assertEqual(calls.iloc[0]['pair_umi_total'], 2)
                merged, _ = STEP4.stream_library_umi_dedup_sources([('dedicated', matched), ('gex', matched)], 'L1', 'sample', None,
                                                                   columns, require_resolved_ibar_sources={'gex'})
                self.assertEqual(merged['umi_count'].sum(), 2)
                self.assertEqual(merged['read_count_gex'].sum(), 6)

    def test_symmetric_matching_preserves_legacy_calls(self):
        table = reference('TGACTAG')
        ref_path = self.root / 'reference.csv'
        table.to_csv(ref_path, index=False)
        reads = self.write_reads(table, per_unit=2)
        outputs = []
        for design_aware in (False, True):
            config = copy.deepcopy(self.config)
            config['guide_design']['enabled'] = design_aware
            config['grna_structure']['ibar_length'] = 7
            config['grna_structure'].pop('ibar_lengths_by_position')
            output = self.root / str(design_aware)
            config['paths'] = {'grna_reference': str(ref_path), 'output_dir': str(output)}
            path = self.root / f'config_{design_aware}.yaml'
            path.write_text(yaml.safe_dump(config))
            subprocess.run([sys.executable, '-B', str(ROOT / 'scripts/03_match_grna_ibar.py'), '--config', str(path),
                            '--barcode-input', str(reads)], check=True, capture_output=True, text=True)
            outputs.append(pd.read_csv(output / 'grna_matches/L1_grna_matches.tsv.gz', sep='\t', keep_default_na=False))
        pd.testing.assert_frame_equal(outputs[0], outputs[1][STEP3.OUTPUT_COLS])

    def test_frozen_orientation_is_used_when_production_pilot_is_disabled(self):
        table = reference('')
        self.config['guide_design']['ibar_orientation'] = 'auto'
        reads = self.write_reads(table, reverse=True)
        contract = compile_design(table, self.config)
        _, effective = check_design(reads, contract, self.config, STEP3.find_tso)
        frozen_path = self.root / 'frozen.yaml'
        frozen_path.write_text(yaml.safe_dump(effective))
        self.config['guide_design']['resolved_contract'] = str(frozen_path)
        self.config['grna_matching']['preflight'] = {'enabled': False}
        ref_path = self.root / 'reference.csv'
        table.to_csv(ref_path, index=False)
        self.config['paths'] = {'grna_reference': str(ref_path), 'output_dir': str(self.root / 'output')}
        path = self.root / 'config.yaml'
        path.write_text(yaml.safe_dump(self.config))
        result = subprocess.run([sys.executable, '-B', str(ROOT / 'scripts/03_match_grna_ibar.py'), '--config', str(path),
                                 '--barcode-input', str(reads)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        matches = pd.read_csv(self.root / 'output/grna_matches/L1_grna_matches.tsv.gz', sep='\t')
        self.assertEqual(len(matches), 240)
        self.assertTrue(matches['pair_id'].eq('pair1').all())

    def test_step0_writes_capabilities_and_preserves_empty_ibar(self):
        design = pd.DataFrame([{
            'spacer_1': reference().iloc[0]['grna_seq'], 'spacer_2': reference().iloc[1]['grna_seq'],
            'iBAR_1': 'ACCTGTA', 'iBAR_2': '', 'target_pos1': 'A', 'target_pos2': 'B',
            'spacer_version_pos1': 'v1', 'spacer_version_pos2': 'v1',
            'spacer_1_target': 'GENE1', 'spacer_2_target': 'GENE2',
            'target': 'GENE1_GENE2', 'target_version': 'v1',
        }])
        design_path = self.root / 'design.csv'
        design.to_csv(design_path, index=False)
        output = self.root / 'reference.csv'
        self.config['design'] = {'design_file': str(design_path), 'rev_comp_ibar': False}
        self.config['paths'] = {'grna_reference': str(output)}
        path = self.root / 'config.yaml'
        path.write_text(yaml.safe_dump(self.config))
        subprocess.run([sys.executable, '-B', str(ROOT / 'scripts/00_format_reference.py'), '--config', str(path)],
                       check=True, capture_output=True, text=True)
        formatted = pd.read_csv(output, keep_default_na=False)
        self.assertEqual(list(formatted['ibar_length']), [7, 0])
        self.assertEqual(formatted.iloc[1]['ibar_seq'], '')
        capabilities = yaml.safe_load(output.with_suffix('.design.yaml').read_text())
        self.assertEqual(capabilities['capabilities']['units_identifying_pair'], 2)

    def test_pilot_reports_wrong_lengths_without_rewriting_the_design(self):
        table = reference('')
        reads = self.write_reads(table)
        expected = compile_design(reference('C'), self.config)
        result, effective = check_design(reads, expected, self.config, STEP3.find_tso)
        self.assertIn('observed_ibar_length_mismatch', {message['code'] for message in result['warnings']})
        self.assertEqual(effective['fingerprint'], expected['fingerprint'])
        expected_table = reference('C')
        expected_table.loc[0, 'grna_seq'] += 'A'
        expected = compile_design(expected_table, self.config)
        result, _ = check_design(reads, expected, self.config, STEP3.find_tso)
        self.assertIn('observed_spacer_length_mismatch', {message['code'] for message in result['warnings']})

    def test_mixed_file_orientation_does_not_auto_flip(self):
        table = reference('')
        self.config['guide_design']['ibar_orientation'] = 'auto'
        paths = []
        for name, reverse, count in [('forward', False, 600), ('reverse', True, 120)]:
            path = self.root / f'{name}.fastq.gz'
            with gzip.open(path, 'wt') as handle:
                for index in range(count):
                    raw = sequence(table.iloc[0], self.config, ibar=reverse_complement(table.iloc[0]['ibar_seq']) if reverse else None)
                    handle.write(f'@r{index}\n{raw}\n+\n' + 'F' * len(raw) + '\n')
            paths.append(path)
        expected = compile_design(table, self.config)
        result, effective = check_design(None, expected, self.config, STEP3.find_tso, fastq_paths=paths)
        self.assertEqual(result['orientation'][1], 'inconclusive')
        self.assertEqual(effective['fingerprint'], expected['fingerprint'])

    def test_full_step4_cli_and_stale_design_rejection(self):
        table = reference('')
        ref_path = self.root / 'reference.csv'
        table.to_csv(ref_path, index=False)
        reads = self.write_reads(table, per_unit=3)
        output = self.root / 'output'
        self.config['paths'] = {'grna_reference': str(ref_path), 'output_dir': str(output)}
        path = self.root / 'config.yaml'
        path.write_text(yaml.safe_dump(self.config))
        environment = {**os.environ, 'MPLCONFIGDIR': '/tmp/csm-design-mpl', 'XDG_CACHE_HOME': '/tmp/csm-design-cache'}
        subprocess.run([sys.executable, '-B', str(ROOT / 'scripts/03_match_grna_ibar.py'), '--config', str(path), '--barcode-input', str(reads)], check=True, capture_output=True, text=True, env=environment)
        result = subprocess.run([sys.executable, '-B', str(ROOT / 'scripts/04_assign_grna_to_cells.py'), '--config', str(path)], capture_output=True, text=True, env=environment)
        self.assertEqual(result.returncode, 0, result.stderr)
        calls = pd.read_csv(output / 'grna_assignment/grna_cell_assignment.csv')
        self.assertEqual(calls.iloc[0]['detection_status'], 'complete')
        molecules = pd.read_csv(output / 'grna_assignment/grna_umi_counts.csv.gz')
        self.assertEqual(molecules['umi_count'].sum(), 2)
        self.config['grna_matching']['max_lev_dist'] = 2
        path.write_text(yaml.safe_dump(self.config))
        result = subprocess.run([sys.executable, '-B', str(ROOT / 'scripts/04_assign_grna_to_cells.py'), '--config', str(path)], capture_output=True, text=True, env=environment)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('changed since Step 3', result.stderr)

    def test_linked_qc_uses_spacers_with_no_ibars(self):
        table = reference('', '')
        second_pair = table.copy()
        second_pair['pair_id'] = 'pair2'
        second_pair['grna_seq'] = ['TTTTACGTACGTACGTACG', 'CCCCCGCATGCATGCATGCA']
        table = pd.concat([table, second_pair], ignore_index=True)
        structure = self.config['grna_structure']
        left = 'GATTACAACGTAGCTACGTA'
        settings = {'guide_design': {'enabled': True}, 'stem_1': structure['fixed_seq_2'],
                    'stem_2': structure['fixed_seq_3'], 'tracr_1': structure['fixed_seq_4'],
                    'tracr_2': structure['fixed_seq_5'], 'grna_lengths': [19, 20],
                    'R1_primer': structure['fixed_seq_1'], 'R2_left_anchor': left,
                    'R2_primer': reverse_complement(structure['fixed_seq_5']), 'stem_max_mm': 2}
        pairs = []
        for index, second_index in enumerate((1, 3)):
            first = sequence(table.iloc[0], self.config)
            second_forward = left + sequence(table.iloc[second_index], self.config)[len(structure['fixed_seq_1']):]
            pairs.append((f'r{index}', first, reverse_complement(second_forward)))
        def read_pairs(first, second):
            return iter(pairs)
        def find_anchor(raw, anchor, maximum):
            index = raw.upper().find(anchor.upper())
            return (index, index + len(anchor), 0) if index >= 0 else (-1, -1, -1)
        result = run_guideseq([('library', 'unused1', 'unused2')], table, settings, self.root,
                              find_anchor, read_pairs, matcher_factory=STEP3.NeighborhoodMatcher)
        self.assertEqual(result['merged_rates']['N_eligible'], 2)
        self.assertEqual(result['merged_rates']['recombination_compatible_rate'], 0.5)
        self.assertEqual(result['combo_frequencies']['count'].sum(), 1)
        first_path, second_path = self.root / 'R1.fastq.gz', self.root / 'R2.fastq.gz'
        for path, position in ((first_path, 1), (second_path, 2)):
            with gzip.open(path, 'wt') as handle:
                for read in pairs:
                    raw = read[position]
                    handle.write(f'@{read[0]}\n{raw}\n+\n' + 'F' * len(raw) + '\n')
        result = Q1.run_guideseq_qc([('library', str(first_path), str(second_path))], table, settings, self.root / 'wrapper')
        self.assertEqual(result['merged_rates']['recombination_compatible_rate'], 0.5)
        duplicate = table.iloc[:2].copy()
        duplicate['pair_id'] = 'pair3'
        repeated = pd.concat([table, duplicate], ignore_index=True)
        result = run_guideseq([('library', 'unused1', 'unused2')], repeated, settings, self.root / 'ambiguous',
                              find_anchor, read_pairs, matcher_factory=STEP3.NeighborhoodMatcher)
        self.assertEqual(result['merged_rates']['design_compatible_rate'], 0.5)
        self.assertEqual(result['merged_rates']['uniquely_resolved_pair_rate'], 0)
        self.assertIsNone(result['merged_rates']['decoupled_pos1'])

    def test_unavailable_qc_is_null_with_reason(self):
        calls = pd.DataFrame([{'bc_wells': 'cell', 'detection_status': 'ambiguous'}])
        molecules = pd.DataFrame([{'bc_wells': 'cell', 'ibar_status': 'absent_by_design', 'umi_count': 3}])
        result = cell_design_qc(calls, molecules, self.root)
        self.assertIsNone(result['legacy_recombination_rate'])
        self.assertTrue(result['reason'])

    def test_gex_filter_uses_reference_lengths_without_ibar_requirements(self):
        path = self.root / 'reference.csv'
        reference('').to_csv(path, index=False)
        config = copy.deepcopy(self.config)
        config['paths'] = {'grna_reference': str(path)}
        config['grna_structure']['grna_lengths'] = [20, 21]
        settings = settings_from_config(config, self.root / 'config.yaml')
        self.assertEqual(settings.guide_lengths, (19, 20))

    def test_step1c_accepts_design_absence_and_rejects_scaffold_conflict(self):
        path = self.root / 'L1_grna_matches.tsv.gz'
        rows = [dict(grna_dist=0, ibar_dist=0, ibar_seq_corrected='', ibar_decoupled=False,
                     ibar_status='absent_by_design', structure_status=status) for status in ('concordant', 'conflict')]
        pd.DataFrame(rows).to_csv(path, sep='\t', index=False)
        _, _, counts = STEP1C.install_strict_match_product(path, self.config, 'L1')
        self.assertEqual(counts['retained_matches'], 1)
        self.assertEqual(counts['dropped_structure_conflict'], 1)
        self.assertEqual(counts['dropped_unresolved_ibar'], 0)

    def test_q1_loader_inherits_design_and_resolves_contract_path(self):
        base = self.root / 'base.yaml'
        config = copy.deepcopy(self.config)
        config['guide_design']['resolved_contract'] = 'frozen.yaml'
        base.write_text(yaml.safe_dump(config))
        qc = self.root / 'qc.yaml'
        qc.write_text(yaml.safe_dump({'sc_config': 'base.yaml', 'guideseq': {'stem_1': 'AAAA'}}))
        merged = load_guideseq_config(qc)
        self.assertTrue(merged['guideseq']['guide_design']['enabled'])
        self.assertEqual(merged['guideseq']['guide_design']['resolved_contract'], str(self.root / 'frozen.yaml'))
        self.assertEqual(merged['guideseq']['ibar_lengths_by_position'], 'auto')

    def test_anndata_contract_summary_has_serializable_values(self):
        directory = self.root / 'grna_matches'
        directory.mkdir()
        contract = compile_design(reference(''), self.config)
        (directory / 'L1_design.yaml').write_text(yaml.safe_dump(contract))
        metadata = STEP5.load_design_metadata(self.root, self.config['libraries'])
        self.assertIsInstance(metadata['L1'], str)
        parsed = json.loads(metadata['L1'])
        self.assertEqual(parsed['fingerprint'], contract['fingerprint'])
        self.assertEqual(parsed['profiles']['2']['ibar_lengths'], [0])


if __name__ == '__main__':
    unittest.main()

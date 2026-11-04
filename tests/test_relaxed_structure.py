"""Shared recovery policy, compact provenance, and safe guide-only restarts."""
import copy
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import yaml

from scripts import relaxed_structure as shared
from scripts.grna_recovery import settings, fingerprint
from scripts.guide_design import reverse_complement
from scripts.validate_step1c_manifests import validate as validate_step1c
from scripts.gex_alignment import settings_identity, guide_restart_compatible
from scripts.grna_restart import write_anndata
from scripts.uger_resources import pipeline
from test_guide_design import configuration, reference, sequence, STEP3, STEP4, STEP1C, STEP5


class SharedRecoveryTests(unittest.TestCase):
    def test_config_and_conservative_shared_selection(self):
        cfg = configuration()
        cfg['grna_structure']['relaxed_structure'] = 'auto'
        opts = settings(cfg)
        self.assertTrue(opts['shared'])
        self.assertEqual(opts['max_tso_shift'], 3)
        p = dict(max_loss=3, min_retained=17, max_tso_shift=2)
        q = dict(max_loss=3, min_retained=17, max_tso_shift=3)
        reports = [{'settings': [dict(policy=p, passed=True, discovery_identified_families=n),
                                  dict(policy=q, passed=i == 0, discovery_identified_families=n+10)]}
                   for i, n in enumerate((100, 200))]
        self.assertEqual(shared.choose_shared(reports, opts), dict(p, enabled=True))
        reports[1]['settings'][0]['passed'] = False
        self.assertFalse(shared.choose_shared(reports, opts)['enabled'])
        cfg['grna_recovery'] = {'mode': 'auto'}
        with self.assertRaises(ValueError):
            settings(cfg)
        cfg.pop('grna_recovery')
        for value in ('manual', 1, None):
            cfg['grna_structure']['relaxed_structure'] = value
            with self.assertRaises(ValueError):
                settings(cfg)

    def test_auto_shares_validated_orientation_and_rejects_conflicts(self):
        cfg = configuration()
        cfg['grna_structure']['relaxed_structure'] = 'auto'
        cfg['libraries'].append(dict(sublibrary_name='L2', sample_name='sample'))
        ref = reference('TGACTAG')
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cfg['paths'] = dict(output_dir=tmp)
            (root/'grna_barcodes').mkdir()
            raw = sequence(ref.iloc[0], cfg)
            for name in ('L1', 'L2'):
                pd.DataFrame([dict(read_id='r',bc_wells='cell',bc_bci='bci',umi='ACGT',
                                   r1_seq=raw,r1_qual='I'*len(raw))]).to_csv(
                    root/f'grna_barcodes/{name}_barcodes.tsv',sep='\t',index=False)
            orientation = dict(orientation={1:'reverse_complement',2:'as_reference'})
            report = dict(preflight=dict(settings=[dict(policy=dict(max_loss=2,min_retained=18,max_tso_shift=2),
                         passed=True,discovery_identified_families=10)]))
            with patch.object(shared, 'check_design', return_value=(orientation, {})), \
                 patch.object(shared.recovery, 'preflight', return_value=report):
                policy = shared.build(cfg,ref,{},lambda worker:None,lambda rows:([],[],{}),STEP3.find_tso)
            self.assertEqual(policy['effective']['max_loss'],2)
            self.assertEqual(policy['libraries'],['L1','L2'])
            oriented = {row['grna_name']:row['ibar_seq'] for row in policy['design_contract']['rows']}
            self.assertEqual(oriented['g1'],reverse_complement(ref.iloc[0].ibar_seq))
            self.assertEqual(shared.load(cfg,ref)['fingerprint'],policy['fingerprint'])
            # Apply the validated orientation to real reads, then ensure Step 4
            # resolves against that same reference rather than the input strand.
            cfg['paths']['grna_reference'] = str(root/'reference.csv')
            ref.to_csv(root/'reference.csv',index=False)
            for name in ('L1','L2'):
                reads = []
                for i,row in enumerate(ref.to_dict('records')):
                    raw = sequence(row,cfg,ibar=reverse_complement(row['ibar_seq']) if row['position']==1 else row['ibar_seq'])
                    reads.append(dict(read_id=f'r{i}',bc_wells='cell',bc_bci='bci',umi='ACGT',r1_seq=raw,r1_qual='I'*len(raw)))
                pd.DataFrame(reads).to_csv(root/f'grna_barcodes/{name}_barcodes.tsv',sep='\t',index=False)
            with patch.object(shared,'check_design',return_value=(orientation,{})), \
                 patch.object(shared.recovery,'preflight',return_value=report):
                policy = shared.build(cfg,ref,{},lambda worker:None,lambda rows:([],[],{}),STEP3.find_tso)
            config_path = root/'config.yaml'
            config_path.write_text(yaml.safe_dump(cfg))
            for name in ('L1','L2'):
                result = subprocess.run([sys.executable,STEP3.__file__,'--config',str(config_path),
                    '--library',name,'--n_processes','1'],capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            product = STEP4._stream_assignment_product(cfg,str(root),'dedicated','gex_grna_identity',
                True,1,1,1,ref,[])
            self.assertEqual(set(product['assignment'].pair_id),{'pair1'})
            self.assertEqual(set(product['assignment'].detection_status),{'complete'})
            conflict = dict(orientation={1:'as_reference',2:'as_reference'})
            with patch.object(shared,'check_design',side_effect=[(orientation,{}),(conflict,{})]):
                with self.assertRaisesRegex(ValueError,'Sublibraries disagree'):
                    shared.build(cfg,ref,{},lambda worker:None,lambda rows:([],[],{}),STEP3.find_tso)

    def test_cli_shared_policy_compact_correction_and_stale_inputs(self):
        cfg = configuration()
        cfg['gex_decontamination'] = dict(step1c=dict(enabled=True, output_namespace='captured_guides'))
        cfg['grna_structure'].update(relaxed_structure=True, relaxed_min_retained=16,
                                      tso_max_edit_dist=3)
        cfg['libraries'].append(dict(sublibrary_name='L2', sample_name='sample'))
        end = len(cfg['grna_structure']['fixed_seq_1'])
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table = reference('')  # Asymmetric iBAR must remain supported.
            table.to_csv(root/'ref.csv', index=False)
            cfg['paths'] = dict(grna_reference=str(root/'ref.csv'), output_dir=str(root))
            reads = []
            for i, row in enumerate(table.to_dict('records')):
                raw = sequence(row, cfg)
                for kind, seq in [('normal', raw), ('short', raw[:end]+raw[end+3:])]:
                    reads.append(dict(read_id=f'{i}_{kind}', bc_wells='01_01_01', bc_bci='1_1_1',
                                      umi='ACGTACGT', r1_seq=seq, r1_qual='I'*len(seq)))
            (root/'grna_barcodes').mkdir()
            for lib in ('L1', 'L2'):
                pd.DataFrame(reads).to_csv(root/f'grna_barcodes/{lib}_barcodes.tsv', sep='\t', index=False)
            config_path = root/'config.yaml'
            config_path.write_text(yaml.safe_dump(cfg))
            command = [sys.executable, STEP3.__file__, '--config', str(config_path)]
            result = subprocess.run(command+['--preflight-only'], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            frozen = shared.load(cfg, table)
            for lib in ('L1', 'L2'):
                result = subprocess.run(command+['--library', lib, '--n_processes', '1'], capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
                path = root/f'grna_matches/{lib}_grna_matches.tsv.gz'
                shared.check_completion(path, lib, frozen)
                calls = pd.read_csv(path, sep='\t', keep_default_na=False)
                self.assertEqual(set(calls.read_id), {r['read_id'] for r in reads})
                self.assertNotIn('recovery_class', calls)
                short = calls[calls.read_id.str.endswith('_short')]
                self.assertEqual(set(short.match_type), {'relaxed_structure'})
                for row in short.itertuples():
                    self.assertTrue(row.grna_seq_corrected.endswith(row.grna_seq))
                    self.assertEqual(row.grna_dist, len(row.grna_seq_corrected)-len(row.grna_seq))
                    self.assertEqual(json.loads(row.structure_correction)['policy'], frozen['fingerprint'])
                STEP1C.install_strict_match_product(path, cfg, lib)
                strict = pd.read_csv(path, sep='\t')
                self.assertEqual(set(strict.read_id), set(calls.read_id))
                with self.assertRaises(ValueError):
                    shared.check_completion(path, lib, frozen)
            # Full Step-1c reruns must rebind a changed config hash, retain prior
            # evidence, and reuse the same policy in the configured namespace.
            commands = []
            for lib in ('L1','L2'):
                evidence = root/f'{lib}_evidence.tsv.gz'
                pd.DataFrame(reads).to_csv(evidence,sep='\t',index=False,compression='gzip')
                manifest = root/f'{lib}_source.json'
                manifest.write_text(json.dumps(dict(complete=True,status='complete',schema_version=4,
                    library=lib,barcode_evidence=str(evidence),barcode_evidence_sha256=STEP1C.sha256_file(evidence),
                    barcode_evidence_rows=len(reads),statistics=dict(csm_gex_filter_positive_valid_barcode_reads=len(reads)))))
                commands.append([sys.executable, STEP1C.__file__, '--config',str(config_path),
                                 '--library',lib,'--input-manifest',str(manifest),'--n_processes','1'])
            for command in commands:
                result = subprocess.run(command,capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            validate_step1c(config_path)
            config_path.write_text(config_path.read_text()+'\n# guide restart config revision\n')
            with self.assertRaisesRegex(ValueError, 'stale config'):
                validate_step1c(config_path)
            for command in commands:
                result = subprocess.run(command,capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                self.assertIn('previous guide evidence is preserved',result.stdout)
            validate_step1c(config_path)
            for lib in ('L1','L2'):
                shared.check_completion(root/f'captured_guides/grna_matches/{lib}_grna_matches.tsv.gz',lib,frozen)
            self.assertEqual(len(list((root/'captured_guides/restarts').glob('*/relocation.json'))),2)
            for lib in ('L1','L2'):
                result = subprocess.run([sys.executable, STEP3.__file__, '--config',str(config_path),
                                         '--library',lib,'--n_processes','1'],capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            product = STEP4._stream_assignment_product(cfg,str(root),'combined','captured_guides',
                                                       True,1,1,1,table,[])
            directory = STEP4._write_assignment_product(product,str(root),'captured_guides')
            shared.validate_assignment(cfg,directory,'combined',frozen)
            self.assertEqual(product['deduped'].read_count.sum(),8)
            self.assertEqual(product['deduped'].umi_count.sum(),4)
            changed = copy.deepcopy(cfg)
            changed['grna_structure']['relaxed_max_tso_shift'] = 1
            with self.assertRaisesRegex(ValueError, 'parameters changed'):
                shared.load(changed, table)
            with (root/'grna_barcodes/L1_barcodes.tsv').open('a') as handle:
                handle.write('\n')
            with self.assertRaisesRegex(ValueError, 'input changed'):
                shared.load(cfg, table)

    def test_compact_provenance_uses_joint_umi_clusters(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            normal = dict(read_id='r1', bc_wells='cell', grna_name='g', target_gene='GENE',
                          umi='ACGTACGT', structure_correction='')
            short = dict(normal, read_id='r2', structure_correction=json.dumps(dict(length='shortened_only', policy='fp')))
            paths = []
            for source, rows in [('dedicated', [normal]), ('gex', [short])]:
                path = root/f'{source}.tsv'
                pd.DataFrame(rows).to_csv(path, sep='\t', index=False)
                paths.append((source, path))
            counts, _ = STEP4.stream_library_umi_dedup_sources(paths, 'L', 'sample', None,
                              ['bc_wells', 'grna_name', 'target_gene'], min_reads_per_umi=2)
            self.assertEqual(counts.umi_count.sum(), 1)
            self.assertEqual(counts.read_count.sum(), 2)
            detail = json.loads(counts.iloc[0].structure_correction)
            self.assertEqual(detail['recovery_umi_mixed'], 1)

    def test_missing_library_cannot_create_partial_relaxed_matrix(self):
        cfg = dict(libraries=[dict(sublibrary_name='L1', sample_name='sample')],
                   grna_structure=dict(relaxed_structure='auto'))
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(FileNotFoundError, 'Incomplete dedicated'):
                STEP4._source_mode_tasks(cfg, tmp, 'dedicated', 'gex_grna_identity')

    def test_hash_reuse_is_limited_to_guide_switches(self):
        cfg = dict(libraries=[dict(sublibrary_name='L1')], gex_alignment=dict(backend='star'),
                   gex_decontamination=dict(enabled=False, step1c=dict(enabled=False),
                                            combine_with_grna_arm=dict(enabled=False)))
        old = settings_identity(cfg)
        cfg['gex_decontamination']['step1c']['enabled'] = True
        cfg['gex_decontamination']['combine_with_grna_arm']['enabled'] = True
        cfg['grna_structure'] = dict(relaxed_structure='auto')
        self.assertTrue(guide_restart_compatible(old, cfg))
        for key, value in [('enabled', True), ('min_anchor', 8)]:
            modified = copy.deepcopy(cfg)
            modified['gex_decontamination'][key] = value
            self.assertFalse(guide_restart_compatible(old, modified))
        cfg['libraries'].append(dict(sublibrary_name='L2'))
        self.assertFalse(guide_restart_compatible(old, cfg))

    def test_restart_keeps_old_annotations_and_compact_metadata(self):
        import anndata as ad
        import numpy as np
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            old = ad.AnnData(np.array([[1, 2]], dtype='int32'),
                             obs=pd.DataFrame(index=pd.Index(['cell'], dtype=object)),
                             var=pd.DataFrame(index=pd.Index(['g1', 'g2'], dtype=object)))
            old.obs['user_label'] = pd.Series(['test'], index=old.obs_names, dtype=object)
            path = root/'combined_adata.h5ad'
            old.write_h5ad(path)
            fresh = ad.AnnData(old.X.copy(), obs=old.obs.drop(columns='user_label').copy(), var=old.var.copy())
            (root/'grna_matches').mkdir()
            (root/'grna_matches/relaxed_structure.json').write_text('{"fingerprint":"fp"}')
            (root/'grna_assignment').mkdir()
            pd.DataFrame([dict(bc_wells='cell', recovery_umi_mixed=1, recovery_umi_shortened_only=2)]).to_csv(
                root/'grna_assignment/grna_recovery_cells_dedicated.csv.gz',index=False)
            STEP5.attach_recovery_metadata(fresh, root, [], compact=True)
            write_anndata(fresh, path)
            saved = ad.read_h5ad(path)
            np.testing.assert_array_equal(saved.X, old.X)
            self.assertNotIn('user_label', saved.obs)
            self.assertIn('grna_structure', saved.uns)
            self.assertEqual(saved.uns['grna_structure']['molecule_counts']['dedicated'].loc['cell','recovery_umi_mixed'], 1)
            backup = list((root/'grna_restarts').glob('*/combined_adata.h5ad'))
            self.assertEqual(len(backup), 1)
            self.assertEqual(ad.read_h5ad(backup[0]).obs.user_label.tolist(), ['test'])

    def test_step5_rejects_stale_count_policy_settings_and_products(self):
        cfg = configuration()
        cfg['grna_structure']['relaxed_structure'] = True
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ref, evidence = root/'reference.csv',root/'matches.tsv'
            reference().to_csv(ref,index=False)
            evidence.write_text('reads')
            cfg['paths'] = dict(grna_reference=str(ref))
            policy = dict(fingerprint='shared-policy')
            record = shared.assignment_record(cfg,policy,[([('dedicated',evidence)],'L1','sample',0)])
            product = root/'counts.csv'
            product.write_text('counted')
            shared.write_assignment_completion(root,'dedicated',record,[product])
            shared.validate_assignment(cfg,root,'dedicated',policy)
            with self.assertRaisesRegex(ValueError,'settings/policy'):
                shared.validate_assignment(cfg,root,'dedicated',dict(fingerprint='different'))
            cfg['guide_assignment']['min_umi_threshold'] = 4
            with self.assertRaisesRegex(ValueError,'settings/policy'):
                shared.validate_assignment(cfg,root,'dedicated',policy)
            cfg['guide_assignment']['min_umi_threshold'] = 1
            product.write_text('changed')
            with self.assertRaisesRegex(ValueError,'products'):
                shared.validate_assignment(cfg,root,'dedicated',policy)

    def test_shared_scheduling_and_guide_only_restart(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cfg = dict(paths=dict(project_dir=tmp, output_dir=str(root/'output')),
                       libraries=[dict(sublibrary_name='L1'), dict(sublibrary_name='L2')],
                       grna_structure=dict(relaxed_structure='auto'),
                       gex_decontamination=dict(step1c=dict(enabled=True)))
            path = root/'config.yaml'
            path.write_text(yaml.safe_dump(cfg))
            for steps, number in [('all', 9), ('3,4,5', 6)]:
                with patch('scripts.uger_resources.submit', side_effect=[str(i) for i in range(number)]) as submit:
                    pipeline(path, steps, False)
                names = [c.args[1] for c in submit.call_args_list]
                self.assertLess(names.index('step3_pilot'), names.index('step3'))
                self.assertLess(names.index('step3'), names.index('step1c'))
                pilot = submit.call_args_list[names.index('step3_pilot')].args[2]
                self.assertIn('CSM_GRNA_PREFLIGHT=true', ','.join(pilot))
                self.assertNotIn('-t', pilot)
                if steps != 'all':
                    self.assertNotIn('step1a', names)
                    self.assertNotIn('step1b', names)


if __name__ == '__main__':
    unittest.main()

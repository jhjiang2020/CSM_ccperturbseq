"""Assignment invariants that distinguish observed guides from construct inference."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np
import pandas as pd

PATH=Path(__file__).resolve().parents[1]/'scripts/benchmark_guide_methods.py'
spec=importlib.util.spec_from_file_location('bench_methods',PATH)
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)


class AssignmentContract(unittest.TestCase):
    def test_guide_identity_uses_spacer_and_position(self):
        self.assertNotEqual(m.canonical_id('AAAA',1,'ctrl'),m.canonical_id('CCCC',1,'ctrl'))
        self.assertNotEqual(m.canonical_id('AAAA',1,'ctrl'),m.canonical_id('AAAA',2,'ctrl'))

    def test_matrix_retains_zeros_and_sums_distinct_ibar_units(self):
        f=pd.DataFrame({'bc_wells':['c1','c1','other'],'grna_name':['A','A','B'],'umi_count':[2,3,9]})
        np.testing.assert_array_equal(m.counts_matrix(f,['c1','c2'],['A','B']).toarray(),[[5,0],[0,0]])

    def test_unknown_guide_is_error(self):
        with self.assertRaises(ValueError):
            m.counts_matrix(pd.DataFrame({'bc_wells':['c'],'grna_name':['bad'],'umi_count':[1]}),['c'],['A'])

    def test_unpositioned_units_are_audited_without_imputing_arm(self):
        frame=pd.DataFrame({'bc_wells':['c','c'],'grna_name':['A','unit_x'],
                            'guide_unit_id':['unit_a','unit_x'],'position':[1,0],'umi_count':[5,1]})
        named,excluded=m.identifiable_counts(frame,['A','B'])
        self.assertEqual(excluded.grna_name.tolist(),['unit_x'])
        np.testing.assert_array_equal(m.counts_matrix(named,['c'],['A','B']).toarray(),[[5,0]])
        frame.loc[1,'position']=1
        with self.assertRaises(ValueError): m.identifiable_counts(frame,['A','B'])

    def test_pair_sum_and_imputation_are_different_from_per_guide_threshold(self):
        ref=pd.DataFrame({'pair_id':['p','p'],'position':[1,2],'grna_name':['A','B']})
        f=pd.DataFrame({'bc_wells':['c1','c2'],'pair_id':['p','p'],'pos1_umi':[3,5],'pos2_umi':[3,0]})
        observed=m.pair_calls(f,['c1','c2'],['A','B'],ref)
        inferred=m.pair_calls(f,['c1','c2'],['A','B'],ref,True)
        np.testing.assert_array_equal(observed,[[True,True],[True,False]])
        np.testing.assert_array_equal(inferred,[[True,True],[True,True]])
        self.assertFalse((np.array([3,3])>=5).any())

    def test_positive_jaccard_ignores_shared_absences(self):
        a=np.zeros((10,10),bool); b=a.copy(); a[0,0]=1; b[0,1]=1
        self.assertEqual(m.agreement(a,b)['positive_jaccard'],0)
        self.assertEqual(m.agreement(a,b)['union'],2)

    def test_crispat_native_ties_and_threshold(self):
        import tempfile
        import crispat
        import anndata as ad
        from scipy import sparse
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); inp=root/'input.h5ad'
            ad.AnnData(sparse.csr_matrix([[0,0],[3,3],[5,1]]),
                       obs=pd.DataFrame(index=['zero','tie','unique']),
                       var=pd.DataFrame(index=['A','B'])).write_h5ad(inp)
            dest=root/'max'; crispat.ga_max(str(inp),str(dest)+'/',UMI_threshold=1)
            x=m.read_crispat_assignments(dest,['zero','tie','unique'],['A','B'])
            np.testing.assert_array_equal(x,[[False,False],[False,False],[True,False]])
            dest=root/'umi'; crispat.ga_umi(str(inp),[5],str(dest)+'/')
            x=m.read_crispat_assignments(dest,['zero','tie','unique'],['A','B'])
            np.testing.assert_array_equal(x,[[False,False],[False,False],[True,False]])

    def test_export_and_concordance_keep_observed_and_imputed_separate(self):
        import tempfile
        import json
        from scipy import sparse
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); (root/'_COUNTS_SUCCESS').touch()
            cells=pd.DataFrame({'cell_id':['C1a','C2a','T1a','T2a'],'sublibrary':m.LIBS,
                                'sample':['Control','Control','Treatment','Treatment']})
            cells.to_csv(root/'cells.csv.gz',index=False)
            ref=pd.DataFrame({'grna_name':['A','B'],'original_grna_name':['old_A','old_B'],
                              'grna_seq':['AAAA','CCCC'],'position':[1,2],
                              'target_gene':['GENE','ctrl'],'pair_id':['p','p']})
            ref.to_csv(root/'reference.csv',index=False)
            units=pd.DataFrame({'bc_wells':['C1a','C1a','C2a'], 'grna_name':['A','B','A'],
                                'umi_count':[5,1,5], 'pair_id':['p','p','p'], 'position':[1,2,1],
                                'structure_status':'concordant','ibar_status':'present_resolved'})
            pairs=pd.DataFrame({'bc_wells':['C1a','C2a'],'pair_id':['p','p'],
                                'pos1_umi':[5,5],'pos2_umi':[1,0]})
            for arm in ['standard','recovered']:
                dest=root/arm/'grna_assignment'; dest.mkdir(parents=True)
                units.to_csv(dest/'grna_umi_counts.csv.gz',index=False)
                pairs.to_csv(dest/'grna_cell_assignment.csv',index=False)
            units.to_csv(root/'recovered/grna_assignment/grna_umi_counts_by_source.csv.gz',index=False)
            pairs.to_csv(root/'recovered/grna_assignment/grna_cell_assignment_combined.csv',index=False)
            m.export(root,include_historical=False)
            observed=sparse.load_npz(root/'assignments/combined_pair5_observed/calls.npz').toarray()
            imputed=sparse.load_npz(root/'assignments/combined_pair5_imputed/calls.npz').toarray()
            rescue=sparse.load_npz(root/'assignments/ibar_observed_rescue1/calls.npz').toarray()
            np.testing.assert_array_equal(observed,[[1,1],[1,0],[0,0],[0,0]])
            np.testing.assert_array_equal(imputed,[[1,1],[1,1],[0,0],[0,0]])
            np.testing.assert_array_equal(rescue,observed)
            for family in m.FAMILIES:
                dest=root/'fits'/family; dest.mkdir(parents=True)
                (dest/'_FINISHED').touch(); (dest/'status.json').write_text(json.dumps({'failures':[]}))
            m.concordance(root)
            completeness=pd.read_csv(root/'concordance/pair_completeness.csv')
            row=completeness[(completeness.method=='combined_pair5_imputed') & (completeness.library=='all')].iloc[0]
            self.assertEqual(row.both_observed_called,1)
            self.assertEqual(row.one_observed_called,1)
            self.assertEqual(row.called_without_same_pair_observation,1)
            self.assertTrue((root/'concordance/_SUCCESS').is_file())

    def test_all_mixture_interfaces(self):
        import tempfile
        import crispat
        import anndata as ad
        from scipy import sparse
        rng=np.random.default_rng(812)
        counts=rng.poisson(.2,(400,4))
        for i in range(400): counts[i,i%4]+=rng.poisson(25)
        names=[f'c{i}' for i in range(400)]; guides=['A','B','C','D']
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); inp=root/'input.h5ad'
            ad.AnnData(sparse.csr_matrix(counts),obs=pd.DataFrame({'batch':1},index=names),
                       var=pd.DataFrame(index=guides)).write_h5ad(inp)
            for family in m.FAMILIES:
                if family in ['umi','max']: continue
                with self.subTest(family=family):
                    options={'UMI_threshold':1}
                    if family in ['ratio','quantiles']: options['thresholds']=[.5]
                    else: options['n_iter']=25
                    if family in ['gauss','poisson_gauss']: options.update(n_jobs=1,make_plots=False)
                    if family in ['poisson','negative_binomial','binomial']: options.update(parallelize=False,subsample_size=400)
                    dest=root/family
                    getattr(crispat,'ga_'+family)(str(inp),output_dir=str(dest)+'/',**options)
                    x=m.read_crispat_assignments(dest,names,guides)
                    self.assertFalse((x & (counts==0)).any())


if __name__=='__main__': unittest.main()

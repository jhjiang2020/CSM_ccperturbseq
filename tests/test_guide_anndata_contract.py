import tempfile
import unittest
import importlib.util
from pathlib import Path

import anndata as ad
import pandas as pd
from scipy import sparse

from scripts.guide_anndata import assignment_matrix, attach_contract, guide_table


class GuideAnnDataContract(unittest.TestCase):
    def test_pair_sum_and_missing_arm_survive_roundtrip(self):
        counts = sparse.csr_matrix([[3,3],[5,0],[9,0],[0,0]])
        cells = ['c1','c2','c3','c4']
        guides = ['a','b']
        call_rows = pd.DataFrame(dict(bc_wells=cells[:3], pair_id=['ab','ab',''],
            pos1_grna=['a','a','a'],pos2_grna=['b','',''],pos1_umi=[3,5,9],pos2_umi=[3,0,0],
            inferred_pos2_grna=['','b','']))
        primary = assignment_matrix(call_rows,counts,cells,guides,True)
        self.assertEqual(primary.toarray().tolist(),[[1,1],[1,0],[0,0],[0,0]])
        ref = pd.DataFrame(dict(grna_name=guides,grna_seq=['AAAA','CCCC'],target_gene=['A','B'],
                                position=[1,2],ibar_seq=['AG',''],pair_id=['ab','ab']))
        evidence = pd.DataFrame(dict(bc_wells=['c1','c1','c2','c3','c4'],grna_name=['a','b','a','a','unit_x'],
                                     umi_count=[3,3,5,9,1]))
        obj = ad.AnnData(sparse.csr_matrix([[1],[2],[3],[4]]),obs=pd.DataFrame(index=cells))
        obj.uns['gRNA_names'] = guides
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'counts.csv.gz'
            evidence.to_csv(path,index=False)
            attach_contract(obj,counts,primary,evidence,call_rows,ref,path,
                            {'dual_guide_counting_mode':True,'min_umi_threshold':5},'combined',[],path.name)
            output = Path(temp)/'test.h5ad'
            if hasattr(getattr(ad,'settings',None),'allow_write_nullable_strings'):
                with ad.settings.override(allow_write_nullable_strings=True):
                    obj.write_h5ad(output)
            else:
                obj.write_h5ad(output)
            check = ad.read_h5ad(output)
        self.assertEqual(check.obs.grna_primary_n_guides.tolist(),[2,1,0,0])
        self.assertEqual(check.obs.grna_observed_umis.tolist(),[6,5,9,1])
        self.assertEqual(check.obs.grna_unindexed_umis.tolist(),[0,0,0,1])
        self.assertEqual(check.uns['grna_inferred_partners'].guide_id.tolist(),['b'])
        self.assertEqual(check.uns['grna_features'].ibar_design.tolist(),['present','none'])

    def test_ambiguous_named_guide_is_marked(self):
        reference = pd.DataFrame(dict(grna_name=['a','a'],grna_seq=['AAAA','CCCC'],target_gene=['A','A']))
        self.assertEqual(guide_table(['a'],reference,reference).identity_status.tolist(),['ambiguous'])

    def test_single_honors_accepted_rows(self):
        counts = sparse.csr_matrix([[8,6],[1,5],[0,0]])
        frame = pd.DataFrame(dict(bc_wells=['c1','c2'],grna_name=['a','b'],umi_count=[8,5]))
        self.assertEqual(assignment_matrix(frame,counts,['c1','c2','c3'],['a','b'],False).toarray().tolist(),
                         [[1,0],[0,1],[0,0]])

    def test_fixed_guide_cutoff_sums_units_before_top_or_discard(self):
        script = Path(__file__).resolve().parents[1]/'scripts/04_assign_grna_to_cells.py'
        spec = importlib.util.spec_from_file_location('step4',script)
        step4 = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(step4)
        evidence = pd.DataFrame(dict(bc_wells=['c1','c1','c1','c2','c2'],grna_name=['a','a','b','a','a'],
                                     umi_count=[3,3,5,4,2],ibar_seq=['AG','AT','CA','AG','AT']))
        top = step4.assign_single_guides(evidence,5,'top')
        self.assertEqual(top.groupby('bc_wells').umi_count.sum().to_dict(),{'c1':6,'c2':6})
        self.assertEqual(set(top.grna_name),{'a'})
        discarded = step4.assign_single_guides(evidence,5,'discard')
        self.assertEqual(set(discarded.bc_wells),{'c2'})


if __name__=='__main__':
    unittest.main()

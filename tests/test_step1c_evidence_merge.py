import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]


def load_script(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


step4 = load_script("step4_step1c_merge", "04_assign_grna_to_cells.py")
step5 = load_script("step5_step1c_merge", "05_build_anndata.py")


class Step1cEvidenceMergeTests(unittest.TestCase):
    def test_gex_product_allows_complete_zero_match_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output = root / "output"
            match_dir = output / "gex_grna_identity/grna_matches"
            match_dir.mkdir(parents=True)
            pd.DataFrame(columns=[
                "bc_wells", "grna_name", "target_gene", "umi"
            ]).to_csv(
                match_dir / "L1_grna_matches.tsv.gz",
                sep="\t", index=False, compression="gzip",
            )
            cfg = {
                "paths": {"output_dir": str(output)},
                "libraries": [{
                    "sublibrary_name": "L1", "sample_name": "S1"
                }],
                "guide_assignment": {
                    "dual_guide_counting_mode": False,
                    "min_umi_threshold": 1,
                    "multiplet_strategy": "all",
                    "min_dominance_fraction": 0.0,
                },
            }
            product = step4._stream_assignment_product(
                cfg, str(output), "gex", "gex_grna_identity",
                False, 1, 1, 1, None, [],
            )
            self.assertEqual(product["matrix"].shape, (0, 0))
            self.assertEqual(len(product["assignment"]), 0)
            assign_dir = step4._write_assignment_product(
                product, str(output), "gex_grna_identity"
            )
            self.assertTrue(
                Path(assign_dir, "step1c_summary.yaml").is_file()
            )
            loaded = step5._read_guide_product(
                step5._guide_product_paths(
                    output, "gex", "gex_grna_identity"
                )
            )
            self.assertEqual(loaded[0].shape, (0, 0))

    @staticmethod
    def _write_matches(path, rows):
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False)

    def test_sources_combine_before_read_filter_and_count_overlap_once(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            dedicated = root / "dedicated.tsv"
            gex = root / "gex.tsv"
            shared = {
                "bc_wells": "01_01_01",
                "grna_name": "guide_A",
                "target_gene": "GENE_A",
                "umi": "AAAAAAAA",
            }
            self._write_matches(dedicated, [shared])
            self._write_matches(gex, [shared, {**shared, "umi": "CCCCCCCC"}])

            result, stats = step4.stream_library_umi_dedup_sources(
                [("dedicated", dedicated), ("gex", gex)],
                "L1",
                "sample",
                0,
                ["bc_wells", "grna_name", "target_gene"],
                max_umi_edits=0,
                min_reads_per_umi=2,
                chunk_size=1,
            )

            self.assertEqual(len(result), 1)
            row = result.iloc[0]
            self.assertEqual(row["bc_wells"], "0_01_01_01")
            self.assertEqual(row["umi_count"], 1)
            self.assertEqual(row["exact_umi_count_overlap"], 1)
            self.assertEqual(row["umi_count_both_sources"], 1)
            self.assertEqual(row["source"], "both")
            self.assertEqual(stats["umis_filtered_low_reads"], 1)

    def test_dual_gex_evidence_requires_resolved_ibar(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "gex.tsv"
            base = {
                "bc_wells": "01_01_01",
                "grna_name": "guide_A",
                "target_gene": "GENE_A",
                "umi": "AAAAAAAA",
                "ibar_seq": "RAWIBAR",
            }
            self._write_matches(path, [
                {**base, "ibar_seq_corrected": ""},
                {
                    **base,
                    "umi": "CCCCCCCC",
                    "ibar_seq_corrected": "IBAR001",
                },
            ])

            result, stats = step4.stream_library_umi_dedup_sources(
                [("gex", path)],
                "L1",
                "sample",
                0,
                [
                    "bc_wells",
                    "grna_name",
                    "target_gene",
                    "ibar_seq_corrected",
                    "ibar_seq_unresolved",
                ],
                require_resolved_ibar_sources={"gex"},
            )

            self.assertEqual(len(result), 1)
            self.assertEqual(result.iloc[0]["ibar_seq_corrected"], "IBAR001")
            self.assertEqual(stats["source_stats"]["gex"]["input_reads"], 2)
            self.assertEqual(stats["source_stats"]["gex"]["eligible_reads"], 1)

    def test_gex_product_uses_separate_namespace_and_documented_names(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            product = {
                "mode": "gex",
                "deduped": pd.DataFrame({
                    "bc_wells": ["0_bc"],
                    "grna_name": ["guide_A"],
                    "umi_count": [1],
                }),
                "assignment": pd.DataFrame({
                    "bc_wells": ["0_bc"],
                    "grna_name": ["guide_A"],
                }),
                "matrix": sparse.csr_matrix([[1]]),
                "barcodes": ["0_bc"],
                "guides": ["guide_A"],
                "lib_stats": {"L1": {"total_reads": 1}},
            }
            step4._write_assignment_product(
                product, str(root), "gex_grna_identity"
            )
            assign_dir = root / "gex_grna_identity/grna_assignment"
            for name in (
                "gex_grna_assignments.csv.gz",
                "gex_grna_umi_counts.csv.gz",
                "gex_grna_umi_matrix.mtx",
                "gex_grna_matrix_barcodes.txt",
                "gex_grna_matrix_guides.txt",
                "step1c_summary.yaml",
            ):
                self.assertTrue((assign_dir / name).exists(), name)

    def test_combined_product_resolves_both_namespaces_and_calls_guides(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "output"
            dedicated_dir = output / "grna_matches"
            gex_dir = output / "gex_grna_identity/grna_matches"
            dedicated_dir.mkdir(parents=True)
            gex_dir.mkdir(parents=True)
            row = {
                "bc_wells": "01_01_01",
                "grna_name": "guide_A",
                "target_gene": "GENE_A",
                "umi": "AAAAAAAA",
            }
            self._write_matches(dedicated_dir / "L1_grna_matches.tsv", [row])
            self._write_matches(gex_dir / "L1_grna_matches.tsv", [row])
            cfg = {
                "libraries": [{
                    "sublibrary_name": "L1",
                    "sample_name": "sample",
                }],
                "guide_assignment": {
                    "dual_guide_counting_mode": False,
                    "min_umi_threshold": 1,
                    "multiplet_strategy": "all",
                    "min_dominance_fraction": 0.0,
                },
            }
            product = step4._stream_assignment_product(
                cfg,
                str(output),
                "combined",
                "gex_grna_identity",
                dedup_ibar=False,
                max_umi_edits=0,
                min_reads_per_umi=1,
                ncpus=1,
                ref_df=pd.DataFrame({"grna_name": ["guide_A"]}),
                available_cols=["grna_name"],
            )
            self.assertEqual(product["matrix"].shape, (1, 1))
            self.assertEqual(product["matrix"][0, 0], 1)
            self.assertEqual(product["assignment"].iloc[0]["source"], "both")

    def test_step5_sparse_alignment_preserves_order_without_densifying(self):
        matrix = sparse.csr_matrix([[1, 0], [0, 3]])
        aligned = step5._align_guide_matrix(
            matrix,
            ["cell_B", "cell_A"],
            ["cell_A", "cell_missing", "cell_B"],
        )
        np.testing.assert_array_equal(
            aligned.toarray(),
            np.array([[0, 3], [0, 0], [1, 0]]),
        )

    def test_step5_preserves_dedicated_pair_metadata_for_comparison(self):
        assignments = pd.DataFrame({
            "bc_wells": ["cell_A", "cell_A"],
            "combo_type": ["multiplet", "multiplet"],
            "detection_status": ["complete", "inferred"],
            "pair_id": ["pair_1", "pair_2"],
            "target_combo": ["A_B", "C_D"],
            "pos1_grna": ["g1", "g3"],
            "pos2_grna": ["g2", ""],
            "pos1_umi": [3, 2],
            "pos2_umi": [3, 0],
            "n_combos": [2, 2],
        })
        summary = step5._assignment_comparison_table(assignments)
        self.assertEqual(summary.loc["cell_A", "pair_id"], "pair_1,pair_2")
        self.assertEqual(summary.loc["cell_A", "pos1_umi"], 5)
        self.assertEqual(summary.loc["cell_A", "n_combos"], 2)


if __name__ == "__main__":
    unittest.main()

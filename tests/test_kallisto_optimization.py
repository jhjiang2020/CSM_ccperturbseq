"""Exact benchmark comparisons must reject even one-bit fractional changes."""
import tempfile
import csv
import json
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as sio, sparse

from scripts.benchmark_kallisto_optimization import exact_product, serial_repaired_reference
from scripts.run_kallisto_gex import write_product


class ExactOptimizationTests(unittest.TestCase):
    def test_serial_reference_is_exact_and_limited_to_documented_restored_wells(self):
        with tempfile.TemporaryDirectory() as tmp:
            baseline = Path(tmp)
            arm = "kallisto_fractional"
            process = baseline / arm / "gex_C1/process"
            process.mkdir(parents=True)
            archive = baseline / "fractional_metadata_originals" / arm / "C1"
            archive.mkdir(parents=True)
            features = pd.DataFrame({"gene_id": [f"G{i}" for i in range(7)]})
            obs = pd.DataFrame({"tscp_count": [3, 3]}, index=pd.Index(["01_02_03", "01_02_04"], name="bc_wells"))
            obs.iloc[1:].to_csv(archive / "aw_uf_cell_metadata.csv")
            qc = {"molecules_below_read_threshold": 0, "metadata_roundoff_repair": {
                "complete": True, "existing_gene_weights_unchanged": True, "restored_awuf_wells": 1}}
            (process / "kallisto_qc.json").write_text(json.dumps(qc))
            with (process / "tscp_assignment.csv").open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["bc_wells", "cell_barcode", "gene", "assignment_weight"])
                for barcode, families in (("AA_AA_AA", 2), ("TT_AA_AA", 1)):
                    for _ in range(families):
                        for gene in features.gene_id:
                            writer.writerow(["01_02_03", barcode, gene, 1 / 7])
            saved = sparse.csr_matrix([[float("0.1428571428571428") * 3] * 7, [3.] + [0.] * 6])
            reference, details = serial_repaired_reference(saved, obs, features, baseline=baseline, arm=arm, library="C1")
            np.testing.assert_array_equal(reference.toarray(), [[3 / 7] * 7, [3.] + [0.] * 6])
            self.assertEqual(details["baseline_roundtrip_wells"], 1)
            self.assertEqual(details["baseline_roundtrip_entries"], 7)
            qc["metadata_roundoff_repair"]["complete"] = False
            (process / "kallisto_qc.json").write_text(json.dumps(qc))
            with self.assertRaisesRegex(AssertionError, "No documented"):
                serial_repaired_reference(saved, obs, features, baseline=baseline, arm=arm, library="C1")

    def test_exact_checker_validates_h5ad_and_rejects_fractional_roundoff(self):
        with tempfile.TemporaryDirectory() as tmp:
            old, new = Path(tmp) / "old", Path(tmp) / "new"
            matrix = sparse.csr_matrix([[1 / 3, 2 / 3]])
            features = pd.DataFrame(dict(gene_id=["A", "B"], gene_name=["A", "B"],
                                         genome=["hg38", "hg38"]))
            obs = pd.DataFrame({"tscp_count": [1]}, index=pd.Index(["01_02_03"], name="bc_wells"))
            for path in (old, new):
                write_product(path, matrix, features, obs, {"backend": "kallisto"})
            self.assertTrue(exact_product(old, new)["exact_counts"])
            modified = matrix.copy()
            modified.data[0] = np.nextafter(modified.data[0], 1)
            sio.mmwrite(new / "count_matrix.mtx", modified)
            with self.assertRaisesRegex(AssertionError, "matrix entries changed"):
                exact_product(old, new)
            sio.mmwrite(new / "count_matrix.mtx", matrix)
            import h5py
            with h5py.File(new / "anndata.h5ad", "r+") as handle:
                handle["layers/counts/data"][0] = modified.data[0]
            with self.assertRaisesRegex(AssertionError, "H5AD raw counts disagree"):
                exact_product(old, new)


if __name__ == "__main__":
    unittest.main()

"""The nine-combination report must keep correction choices and cohorts explicit."""

import importlib.util
import json
from pathlib import Path
import tempfile
import types
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy import io as sio, sparse
import yaml

from scripts.benchmark_gex_backends import ARMS, arm_options
from scripts.benchmark_gex_report import report


@unittest.skipUnless(importlib.util.find_spec("matplotlib"), "matplotlib is required for report rendering")
class BenchmarkReportTests(unittest.TestCase):
    def test_all_nine_combinations_and_within_backend_comparisons(self):
        import scripts.benchmark_gex_report as module
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "manifests").mkdir()
            (root / "configs").mkdir()
            libraries = [{"sublibrary_name": name, "sample_name": "Control"}
                         for name in ("C1", "C2", "T1", "T2")]
            config = root / "config.yaml"
            config.write_text(yaml.safe_dump({"libraries": libraries}))
            for arm in ARMS:
                (root / "configs" / f"{arm}.yaml").write_text("{}\n")
                base, flavor = arm_options(arm)
                for lib in libraries:
                    name = lib["sublibrary_name"]
                    source = root / arm / f"gex_{name}"
                    process = source / "process"
                    dge = source / "Control/DGE_filtered"
                    process.mkdir(parents=True)
                    dge.mkdir(parents=True)
                    matrix = sparse.csr_matrix([[1., 4., 3.], [7., 2., 1.], [2., 5., 4.]])
                    if flavor == "none":
                        matrix = matrix * 1.1
                    elif flavor == "directional":
                        matrix = matrix * 0.9
                    if base == "kallisto_fractional":
                        matrix = matrix * 1.2
                    obs = pd.DataFrame(index=pd.Index(["01_01_01", "01_01_02", "01_01_03"], name="bc_wells"))
                    genes = pd.DataFrame({"gene_id": ["A", "B", "C"]})
                    for directory, matrix_name, obs_name in ((dge, "count_matrix.mtx", "cell_metadata.csv"),
                                (process, "aw_uf_count_matrix.mtx", "aw_uf_cell_metadata.csv")):
                        sio.mmwrite(directory / matrix_name, matrix)
                        obs.to_csv(directory / obs_name)
                        genes.to_csv(directory / "all_genes.csv", index=False)
                    (process / "barcode_headLog.final.out").write_text(
                        "Started job on | Oct 22 19:00:00\nFinished on | Oct 22 19:00:20\n")
                    (process / "kallisto_qc.json").write_text(json.dumps({"pseudoaligned_reads": 10, "timing": {
                        "alignment_seconds": 30, "input_format_seconds": 5,
                        "umi_count_seconds": 60, "publish_seconds": 2}}))
                    timing = {
                        "complete": True, "threads": 16, "elapsed_seconds": 100,
                        "shared_preprocess_seconds": 20, "survivor_sha256": "identical-input"}
                    if base == "kallisto_fractional" and flavor != "none":
                        timing["metadata_roundoff_repair"] = {"complete": True, "elapsed_seconds": 2}
                    (root / "manifests" / f"{name}.{arm}.json").write_text(json.dumps(timing))
            args = types.SimpleNamespace(config=config, work_dir=root, pairs=10_000_000)
            with patch("scripts.combine_kallisto_gex.combine") as combination, patch(
                    "scripts.gex_alignment.validate_combined"), patch.object(
                    module, "__file__", str(root / "scripts/benchmark_gex_report.py")), patch(
                    "scripts.gex_alignment.validate_star_policies") as star_policies:
                report(args)
            self.assertEqual(combination.call_count, 6)
            self.assertEqual(star_policies.call_count, 3)
            summary = pd.read_csv(root / "comparison/library_summary.csv")
            agreement = pd.read_csv(root / "comparison/agreement.csv")
            runtime = pd.read_csv(root / "comparison/runtime.csv")
            self.assertEqual(len(summary), 36)
            self.assertEqual(len(agreement), 60)
            self.assertEqual(len(runtime), 36)
            within = agreement[(agreement["first"] == "star") & (agreement["second"] == "star_none")]
            np.testing.assert_allclose(within["fixed_first_cohort_median_umi_ratio"], 1.1)
            np.testing.assert_allclose(agreement["shared_umi_spearman"], 1)
            np.testing.assert_allclose(runtime["core_speedup_vs_same_flavor_star"], 1)
            self.assertEqual(runtime["metadata_republish_seconds"].sum(), 16)
            success = json.loads((root / "BENCHMARK_SUCCESS.json").read_text())
            self.assertEqual(success["arms"], list(ARMS))
            public = (root / "docs/KALLISTO_PILOT_BENCHMARK_REPORT.md").read_text()
            self.assertIn("Core runtimes retain the original", public)
            for suffix in ("", "_directional", "_none"):
                image = f"umi_agreement{suffix}.png"
                self.assertTrue((root / "docs/assets/kallisto_benchmark_10m" / image).is_file())
                self.assertIn(f"](assets/kallisto_benchmark_10m/{image})", public)


if __name__ == "__main__":
    unittest.main()

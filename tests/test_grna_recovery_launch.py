"""Release integration: opt-in scheduling and audit count preservation."""

import copy
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import yaml

from scripts.grna_recovery import settings
from scripts.uger_resources import pipeline
from test_guide_design import STEP4


class RecoveryLaunchTests(unittest.TestCase):
    def test_standard_templates_annotate_one_opt_in(self):
        project = Path(__file__).resolve().parents[1]
        for name in ('config_single.yaml', 'config_dual.yaml', 'config_rgupta.yaml',
                     'config_rgupta_single.yaml', 'config_rgupta_single_no_ibar.yaml'):
            cfg = yaml.safe_load((project / name).read_text())
            self.assertNotIn('grna_recovery', cfg)
            self.assertFalse(cfg['grna_structure']['relaxed_structure'])
            self.assertEqual(settings(cfg)['mode'], 'off')
            for key in ('step1c', 'combine_with_grna_arm'):
                self.assertFalse(cfg['gex_decontamination'][key]['enabled'])
            cfg['grna_structure']['relaxed_structure'] = 'auto'
            self.assertEqual(settings(cfg)['mode'], 'auto')

    def config(self, root, mode):
        path = root / "config.yaml"
        path.write_text(yaml.safe_dump({
            "paths": {"project_dir": str(root), "output_dir": str(root / "output")},
            "libraries": [{"sublibrary_name": "L1", "sample_name": "S1"}],
            "gex_decontamination": {"step1c": {"enabled": True}},
            "grna_recovery": {"mode": mode},
        }))
        return path

    def test_recovery_step1c_waits_for_dedicated_policy(self):
        for mode in ("audit", "auto", "manual"):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as tmp:
                path = self.config(Path(tmp), mode)
                with patch("scripts.uger_resources.submit", side_effect=[str(i) for i in range(100, 108)]) as submit:
                    pipeline(path, "all", False)
                calls = {c.args[1]: c.args[2] for c in submit.call_args_list}
                self.assertEqual(calls["step1c"][calls["step1c"].index("-hold_jid") + 1], "100,103")
                self.assertEqual(calls["step4"][calls["step4"].index("-hold_jid") + 1], "103,105")

    def test_missing_policy_stops_partial_submission(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = self.config(root, "auto")
            with patch("scripts.uger_resources.submit") as submit:
                with self.assertRaisesRegex(ValueError, "frozen dedicated policy"):
                    pipeline(path, "1", False)
                submit.assert_not_called()
            policy = root / "output/grna_matches/L1_recovery.json"
            policy.parent.mkdir(parents=True)
            policy.write_text("{}")  # Content validation belongs to the compute worker.
            with patch("scripts.uger_resources.submit", return_value="100") as submit:
                pipeline(path, "1", False)
                self.assertEqual(submit.call_count, 4)

    def test_yaml_off_preserves_standard_order(self):
        cfg = yaml.safe_load("grna_recovery:\n  mode: off\n")
        self.assertEqual(settings(cfg)["mode"], "off")
        with tempfile.TemporaryDirectory() as tmp:
            path = self.config(Path(tmp), cfg["grna_recovery"]["mode"])
            with patch("scripts.uger_resources.submit", side_effect=[str(i) for i in range(100, 108)]) as submit:
                pipeline(path, "all", False)
            calls = {c.args[1]: c.args[2] for c in submit.call_args_list}
            self.assertEqual(calls["step1c"][calls["step1c"].index("-hold_jid") + 1], "100")
            self.assertEqual(calls["step4"][calls["step4"].index("-hold_jid") + 1], "105,103")

    def test_audit_preserves_standard_counts_at_read_support_boundary(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            row = dict(read_id="same-read", bc_wells="01_01_01", grna_name="g1",
                       target_gene="GENE", umi="ACGTACGT")
            paths = [root / "grna_matches/L1_grna_matches.tsv",
                     root / "gex_grna_identity/grna_matches/L1_grna_matches.tsv"]
            cfg = dict(libraries=[dict(sublibrary_name="L1", sample_name="S1")],
                       guide_assignment=dict(dual_guide_counting_mode=False, min_umi_threshold=1,
                                             multiplet_strategy="all", min_dominance_fraction=0))
            for path in paths:
                path.parent.mkdir(parents=True)
                pd.DataFrame([row]).to_csv(path, sep="\t", index=False)
            baseline = STEP4._stream_assignment_product(
                cfg, str(root), "combined", "gex_grna_identity", False, 0, 2, 1, None, [])
            audit_cfg = copy.deepcopy(cfg)
            audit_cfg["grna_recovery"] = {"mode": "audit"}
            for path in paths:
                pd.DataFrame([{**row, "recovery_class": "baseline",
                               "recovery_length_state": "not_reassessed",
                               "recovery_policy_fingerprint": "audit-policy"}]).to_csv(path, sep="\t", index=False)
            audited = STEP4._stream_assignment_product(
                audit_cfg, str(root), "combined", "gex_grna_identity", False, 0, 2, 1, None, [])
            self.assertEqual(baseline["deduped"].read_count.sum(), 2)
            pd.testing.assert_frame_equal(
                baseline["assignment"], audited["assignment"][baseline["assignment"].columns])
            for column in ("umi_count", "umi_count_raw", "read_count"):
                pd.testing.assert_series_equal(baseline["deduped"][column], audited["deduped"][column])


if __name__ == "__main__":
    unittest.main()

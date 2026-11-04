"""Identity, boundary, layout and molecule-accounting regression tests."""

import copy
import json
import tempfile
import unittest
from pathlib import Path
import pandas as pd
from scripts.grna_recovery import (
    RecoveryEngine,
    settings,
    recovery_contract,
    make_policy,
    validate_saved_policy,
    preflight,
)
from test_guide_design import configuration, reference, sequence, STEP3, STEP4


class RecoveryTests(unittest.TestCase):
    def engine(self, table=None, cfg=None, **changes):
        cfg = cfg or configuration()
        cfg["grna_recovery"] = {
            "mode": "manual",
            "spacer": {"min_retained": 16},
            "tso": {"max_edit_distance": 3},
        }
        opts = settings(cfg)
        opts.update(changes)
        contract = recovery_contract(reference() if table is None else table, cfg)
        return RecoveryEngine(contract, cfg["grna_structure"]["fixed_seq_1"], 2, opts)

    def test_layouts_and_loss_are_relative_to_design(self):
        cfg = configuration()
        end = len(cfg["grna_structure"]["fixed_seq_1"])
        for first, second in [
            ("ACCTGTA", "TGACTAG"),
            ("ACCTGTA", "C"),
            ("ACCTGTA", ""),
            ("", "TGACTAG"),
            ("", ""),
        ]:
            table = reference(second, first)
            engine = self.engine(table)
            for row in table.to_dict("records"):
                seq = sequence(row, cfg)
                for loss in (1, 2, 3):
                    damaged = seq[:end] + seq[end + loss :]
                    result, reason = engine.match(damaged, "I" * len(damaged))
                    self.assertIsNotNone(result, (first, second, loss, reason))
                    self.assertEqual(result["recovery_designed_spacer"], row["grna_seq"])
                    self.assertEqual(
                        result["recovery_ibar_state"],
                        "agree" if row["ibar_seq"] else "absent_by_design",
                    )
                damaged = seq[:end] + seq[end + 4 :]
                self.assertIsNone(engine.match(damaged, "I" * len(damaged))[0])

    def test_quality_truncation_and_errors_fail_closed(self):
        cfg = configuration()
        seq = sequence(reference().iloc[0], cfg)
        engine = self.engine()
        self.assertIsNone(engine.match(seq, None)[0])
        self.assertIsNone(engine.match(seq, "!" * len(seq))[0])
        self.assertIsNone(engine.match(seq[:-10], "I" * (len(seq) - 10))[0])
        end = len(engine.tso)
        seq = seq[: end + 9] + "N" + seq[end + 10 :]
        self.assertIsNone(engine.match(seq, "I" * len(seq))[0])

    def test_foreign_ibar_optional_never_changes_spacer(self):
        cfg = configuration()
        table = reference("TGACTAG")
        extra = table.copy()
        extra["pair_id"] = "pair2"
        extra["grna_name"] = ["g3", "g4"]
        extra["grna_seq"] = ["CCTATCGGATCAGTCCACT", "ATGACGTCTTCGAGATCTGA"]
        extra["ibar_seq"] = ["TGAACCT", "CACATGA"]
        engine = self.engine(pd.concat([table, extra], ignore_index=True))
        seq = sequence(table.iloc[0], cfg, ibar="TGAACCT")
        result, _ = engine.match(seq, "I" * len(seq))
        self.assertEqual(result["recovery_ibar_state"], "discordant")
        self.assertEqual(result["pair_id"], "")
        self.assertEqual(result["recovery_designed_spacer"], table.iloc[0].grna_seq)
        self.assertIsNone(
            engine.match(seq, "I" * len(seq), dict(engine.options, ibar_corroboration="require"))[0]
        )

    def test_ambiguity_cannot_be_resolved_by_ibar_or_length_floor(self):
        table = reference("TGACTAG")
        extra = table.copy()
        extra["pair_id"] = "pair2"
        extra["grna_name"] = ["g3", "g4"]
        extra.loc[0, "grna_seq"] = "TT" + table.iloc[0].grna_seq[2:]
        extra["ibar_seq"] = ["TGAACCT", "CACATGA"]
        engine = self.engine(pd.concat([table, extra], ignore_index=True))
        seq = sequence(table.iloc[0], configuration())
        end = len(engine.tso)
        seq = seq[:end] + seq[end + 2 :]
        result, reason = engine.match(seq, "I" * len(seq))
        self.assertIsNone(result)
        self.assertEqual(reason, "ambiguous_spacer_or_position")

    def test_full_and_short_layouts_are_not_called_physical_damage(self):
        engine = self.engine()
        seq = sequence(reference().iloc[0], configuration())
        result, _ = engine.match(seq, "I" * len(seq))
        self.assertEqual(result["recovery_length_state"], "uncertain")
        self.assertIn(0, json.loads(result["recovery_possible_losses"]))

    def test_frozen_contract_tamper_and_stale_settings(self):
        engine = self.engine()
        policy = make_policy(engine.contract, engine.options)
        self.assertEqual(validate_saved_policy(policy, engine.contract, engine.options), policy)
        tampered = copy.deepcopy(policy)
        tampered["effective"]["max_loss"] = 0
        with self.assertRaises(ValueError):
            validate_saved_policy(tampered, engine.contract, engine.options)
        with self.assertRaises(ValueError):
            validate_saved_policy(policy, engine.contract, dict(engine.options, min_quality=20))

    def test_insufficient_pilot_does_not_enable(self):
        engine = self.engine()
        cfg = configuration()
        opts = dict(engine.options, mode="auto")
        p = preflight([], engine.contract, cfg, opts)
        self.assertFalse(p["effective"]["enabled"])

    def test_original_endpoint_negative_controls_are_actually_evaluated(self):
        cfg = configuration()
        table = reference("TGACTAG")
        table.loc[0, "grna_seq"] = "A" * 19
        engine = self.engine(table, endpoint_policy="original")
        seq = sequence(table.iloc[0], cfg)
        rows = [
            dict(
                read_id=str(i),
                bc_wells=f"cell{i}",
                umi="ACGTACGT",
                r1_seq=seq,
                r1_qual="I" * len(seq),
            )
            for i in range(120)
        ]
        policy = preflight(rows, engine.contract, cfg, dict(engine.options, mode="auto"))
        self.assertTrue(
            all(s["negative_controls"]["accepted"] > 0 for s in policy["preflight"]["settings"])
        )
        self.assertFalse(policy["effective"]["enabled"])

    def test_recovery_does_not_split_umi_clusters(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "matches.tsv"
            rows = [
                dict(
                    bc_wells="cell",
                    grna_name="g",
                    target_gene="GENE",
                    umi=u,
                    recovery_class=c,
                    recovery_length_state=l,
                    recovery_policy_fingerprint="fp",
                )
                for u, c, l in [
                    ("AAAA", "baseline", "not_reassessed"),
                    ("AAAA", "short_spacer", "shortened_only"),
                    ("AAAT", "short_spacer", "shortened_only"),
                ]
            ]
            pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
            df, _ = STEP4.stream_library_umi_dedup(
                str(path), "L", "s", None, ["bc_wells", "grna_name", "target_gene"]
            )
            self.assertEqual(df.umi_count.sum(), 1)
            self.assertEqual(df.recovery_umi_mixed.sum(), 1)
            self.assertEqual(df.recovery_umi_rescue_only.sum(), 0)


class RecoveryWorkflowTests(unittest.TestCase):
    def test_cli_off_audit_manual_and_unvalidated_auto(self):
        import subprocess
        import sys
        import yaml

        cfg = configuration()
        cfg["grna_structure"]["tso_max_edit_dist"] = 3
        table = reference("TGACTAG")
        end = len(cfg["grna_structure"]["fixed_seq_1"])
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            ref = root / "ref.csv"
            table.to_csv(ref, index=False)
            reads = []
            for i, row in enumerate(table.to_dict("records")):
                seq = sequence(row, cfg)
                for label, raw in [("intact", seq), ("short", seq[:end] + seq[end + 3 :])]:
                    reads.append(
                        dict(
                            read_id=f"{i}_{label}",
                            bc_wells="01_01_01",
                            bc_bci="1_1_1",
                            umi="ACGTACGT" + str(i),
                            r1_seq=raw,
                            r1_qual="I" * len(raw),
                        )
                    )
            source = root / "reads.tsv"
            pd.DataFrame(reads).to_csv(source, sep="\t", index=False)
            outputs = {}
            for mode in ("off", "manual", "audit", "auto"):
                output = root / mode
                cfg["paths"] = {"grna_reference": str(ref), "output_dir": str(output)}
                cfg["grna_recovery"] = {"mode": mode, "spacer": {"min_retained": 16}}
                config = root / f"{mode}.yaml"
                config.write_text(yaml.safe_dump(cfg))
                subprocess.run(
                    [
                        sys.executable,
                        str(Path(STEP3.__file__)),
                        "--config",
                        str(config),
                        "--library",
                        "L1",
                        "--barcode-input",
                        str(source),
                        "--n_processes",
                        "1",
                    ],
                    check=True,
                    capture_output=True,
                    text=True,
                )
                outputs[mode] = pd.read_csv(
                    output / "grna_matches/L1_grna_matches.tsv.gz", sep="\t", keep_default_na=False
                )
            self.assertEqual(set(outputs["manual"].read_id), {r["read_id"] for r in reads})
            self.assertEqual(set(outputs["audit"].read_id), set(outputs["off"].read_id))
            self.assertEqual(set(outputs["auto"].read_id), set(outputs["off"].read_id))
            self.assertEqual(set(outputs["manual"].recovery_class), {"baseline", "short_spacer"})
            baseline = outputs["off"].set_index("read_id")
            recovered = outputs["manual"].set_index("read_id")
            pd.testing.assert_frame_equal(
                baseline, recovered.loc[baseline.index, baseline.columns], check_dtype=False
            )

    def test_step2_reservoir_is_reusable_and_rejects_changed_source(self):
        from scripts.grna_recovery import ReservoirWriter, sample_rows

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "input.tsv"
            cache = root / "sample.tsv.gz"
            with source.open("w") as handle:
                writer = ReservoirWriter(handle, 7, 21)
                writer.write("read_id\tbc_wells\tumi\tr1_seq\tr1_qual\n")
                for i in range(100):
                    writer.write(f"r{i}\tc{i}\tu{i}\tACGT\tIIII\n")
            writer.save(cache, source)
            rows, meta = sample_rows(source, 7, 21, cache)
            self.assertEqual(len(rows), 7)
            self.assertEqual(meta["sampling"], "uniform_read_reservoir_during_step2")
            with source.open("a") as handle:
                handle.write("new\tcell\tumi\tACGT\tIIII\n")
            _, meta = sample_rows(source, 7, 21, cache)
            self.assertEqual(meta["reads_scanned"], 101)


class RecoverySourceTests(unittest.TestCase):
    def test_shared_read_is_counted_once_across_source_copies(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            paths = []
            row = dict(
                read_id="instrument:1",
                bc_wells="cell",
                grna_name="g",
                target_gene="GENE",
                umi="ACGTACGT",
                recovery_class="endpoint",
                recovery_length_state="uncertain",
                recovery_policy_fingerprint="fp",
            )
            for source in ["dedicated", "gex"]:
                path = root / f"{source}.tsv"
                pd.DataFrame([row]).to_csv(path, sep="\t", index=False)
                paths.append((source, path))
            frame, stats = STEP4.stream_library_umi_dedup_sources(
                paths, "L", "sample", None, ["bc_wells", "grna_name", "target_gene"]
            )
            self.assertEqual(frame.read_count.sum(), 1)
            self.assertEqual(frame.umi_count.sum(), 1)
            self.assertEqual(frame.umi_count_both_sources.sum(), 1)
            self.assertEqual(stats["total_reads"], 1)
            row["grna_name"] = "wrong"
            pd.DataFrame([row]).to_csv(paths[1][1], sep="\t", index=False)
            with self.assertRaises(ValueError):
                STEP4.stream_library_umi_dedup_sources(
                    paths, "L", "sample", None, ["bc_wells", "grna_name", "target_gene"]
                )


class RecoveryAnnDataTests(unittest.TestCase):
    def test_provenance_roundtrip_preserves_expression(self):
        import anndata as ad
        import numpy as np
        from test_guide_design import STEP5

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "grna_matches").mkdir()
            (root / "grna_assignment").mkdir()
            (root / "grna_matches/L1_recovery.json").write_text(
                json.dumps({"effective": {"max_loss": 3}, "fingerprint": "fp"})
            )
            pd.DataFrame(
                [
                    dict(
                        bc_wells="cell1",
                        recovery_umi_shortened_only=2,
                        recovery_umi_length_uncertain=1,
                    )
                ]
            ).to_csv(root / "grna_assignment/grna_recovery_cells_dedicated.csv.gz", index=False)
            counts = np.asarray([[3, 0], [0, 5]], dtype=np.int32)
            data = ad.AnnData(
                counts.copy(),
                obs=pd.DataFrame(index=pd.Index(["cell1", "cell2"], dtype=object)),
                var=pd.DataFrame(index=pd.Index(["g1", "g2"], dtype=object)),
            )
            STEP5.attach_recovery_metadata(data, root, [{"sublibrary_name": "L1"}])
            data.write_h5ad(root / "test.h5ad")
            saved = ad.read_h5ad(root / "test.h5ad")
            np.testing.assert_array_equal(saved.X, counts)
            self.assertEqual(saved.obs.dedicated_recovery_umi_shortened_only.tolist(), [2, 0])
            self.assertEqual(
                json.loads(saved.uns["grna_recovery_policies"]["L1"])["fingerprint"], "fp"
            )


if __name__ == "__main__":
    unittest.main()

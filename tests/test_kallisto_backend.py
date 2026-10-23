"""Molecule parity, native BUS integration, and fraction-safe DGE contracts."""

import ast
import copy
import csv
import importlib.util
import json
import os
import random
import struct
import subprocess
import tempfile
import types
import unittest
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy import io as sio, sparse

from scripts import gex_alignment as contract
from scripts import kallisto_umi as umi
from scripts.combine_kallisto_gex import combine
from scripts.run_kallisto_gex import count_native, publish_samples, seal, write_product, whole_molecule_totals


ROOT = Path(__file__).resolve().parents[1]


def installed_parse_kernel():
    source = ROOT / "ParseBiosciences-Pipeline.1.6.4/splitpipe/molinfo.py"
    if not source.is_file():
        raise unittest.SkipTest("Locally obtained Parse 1.6.4 source is unavailable")
    tree = ast.parse(source.read_text())
    function = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "collapse_polyN")
    namespace = {"defaultdict": defaultdict, "combinations": combinations,
                 "cy_utils": types.SimpleNamespace(
                     single_mut_seqs=umi.single_mutations,
                     hamdist2=lambda a, b: sum(x != y for x, y in zip(a, b)))}
    # Execute the vendor's unchanged function body, without loading its HPC-only
    # compiled dependencies or changing any source in the Parse distribution.
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(source), "exec"), namespace)
    return namespace["collapse_polyN"]


def write_bus(path, records, bc_len=6, umi_len=4):
    with open(path, "wb") as handle:
        handle.write(struct.pack("<4sIIII", b"BUS\0", 1, bc_len, umi_len, 0))
        for barcode, sequence, ec, count in records:
            handle.write(umi.BUS_RECORD.pack(umi.dna_encode(barcode), umi.dna_encode(sequence), ec, count, 0, 0))


def reference_fixture(root):
    reference = root / "reference"
    reference.mkdir()
    (reference / "index.idx").write_bytes(b"test index")
    (reference / "t2g.txt").write_text("tx1\tA\n" "tx2\tA\n" "tx3\tB\n")
    features = pd.DataFrame([
        ["tx1", "A", "Alpha", "human", "transcript"],
        ["tx2", "A", "Alpha", "human", "transcript"],
        ["tx3", "B", "Beta", "human", "transcript"],
    ], columns=["transcript_id", "gene_id", "gene_name", "genome", "feature_type"])
    features.to_csv(reference / "features.tsv", sep="\t", index=False)
    contract.atomic_json(reference / "reference.json", {
        "schema_version": 1, "complete": True, "workflow": "standard",
        "artifacts": {name: contract.sha256_file(reference / name)
                      for name in ("index.idx", "t2g.txt", "features.tsv")},
    })
    return reference


def native_fixture(root, options, min_reads=1, records=None, threads=1):
    reference = reference_fixture(root)
    bus = root / "output.bus"
    write_bus(bus, records if records is not None else [("AACCGG", "AAAA", 0, 5), ("AACCGG", "AAAT", 0, 2),
                    ("AACCGG", "CCCC", 1, 3), ("AACCGG", "GGGG", 2, 2),
                    ("TTCCGG", "AAAA", 0, 2)])
    (root / "transcripts.txt").write_text("tx1\ntx2\ntx3\n")
    (root / "matrix.ec").write_text("0\t0\n1\t0,1\n2\t0,2\n")
    lookup = {"AACCGG": ("01_02_03", "1_2_3", "T", "AA_CC_GG"),
              "TTCCGG": ("01_02_03", "4_2_3", "R", "TT_CC_GG")}
    genes = pd.DataFrame([["A", '"Alpha"', "human"], ["B", '"Beta"', "human"]],
                         columns=["gene_id", "gene_name", "genome"])
    result = count_native(bus, root / "transcripts.txt", root / "matrix.ec", reference, genes,
                          lookup, 6, 4, root / "reads.sqlite", root / "tscp_assignment.csv", options,
                          min_reads=min_reads, threads=threads)
    return result, reference


class FakeSample:
    def get_name(self):
        return "Control"

    def bc_df_for_sample(self, frame):
        return frame[frame["bc1_wind"] == "01"]


class FakeBarcodeInfo:
    def get_bc_data(self, round_number):
        return pd.DataFrame({"sequence": ["AA"], "well": ["A1"]})

    def get_seq_to_windzp(self, round_number):
        return {"AA": f"{round_number:02d}"}


class FakeParse:
    def get_par_val(self, key, **kwargs):
        return {"mol_awuf_tscp_thresh": 1, "dge_species_purity_thresh": 0.9,
                "dge_tscp_cut_umin": 1}[key]

    def get_samples(self):
        return [FakeSample()]

    def get_bc_info(self):
        return FakeBarcodeInfo()

    def get_tscp_cutoff(self, counts, sample):
        return True, {"cell_tscp_cutoff": 1, "cell_call": "parse_test_fixture"}


class UmiParityTests(unittest.TestCase):
    def test_benchmark_covers_all_backend_and_correction_combinations(self):
        from scripts.benchmark_gex_backends import ARMS, arm_options, runtime_config
        config = {"paths": {}, "gex_decontamination": {"existing_alignment": {},
                  "step1c": {}, "combine_with_grna_arm": {}}}
        self.assertEqual(len(ARMS), 9)
        self.assertEqual(len(set(ARMS)), 9)
        with tempfile.TemporaryDirectory() as tmp:
            args = types.SimpleNamespace(work_dir=Path(tmp), config="unused")
            for arm in ARMS:
                base, flavor = arm_options(arm)
                with patch("scripts.benchmark_gex_backends.load_config", return_value=config):
                    cfg, _ = runtime_config(args, arm)
                self.assertEqual(cfg["gex_alignment"].get("umi_correction", "parse"), flavor)
                self.assertEqual(cfg["gex_alignment"]["backend"], "star" if base == "star" else "kallisto")
                if base != "star":
                    self.assertEqual(cfg["gex_alignment"]["kallisto"]["gene_assignment"],
                                     "fractional" if base == "kallisto_fractional" else "unique")

    def test_three_correction_flavors_have_distinct_semantics(self):
        counts = {"AAAA": 5, "AAAT": 2, "AATT": 1}
        for flavor, expected in (("parse", 2), ("directional", 1), ("none", 3)):
            families = umi.collapse_umis(counts, flavor=flavor)
            self.assertEqual(len(families), expected)
            self.assertEqual(sum(reads for _, reads, _ in families), 8)
            self.assertEqual({member for _, _, members in families for member in members}, set(counts))
        # Strict abundance can merge 5 vs 4; the directional condition cannot.
        self.assertEqual(len(umi.collapse_umis({"AAAA": 5, "AAAT": 4}, flavor="parse")), 1)
        self.assertEqual(len(umi.collapse_umis({"AAAA": 5, "AAAT": 4}, flavor="directional")), 2)
        # Singleton ties are directional edges; larger equal counts are not.
        self.assertEqual(len(umi.collapse_umis({"AAAA": 1, "AAAT": 1}, flavor="directional")), 1)
        self.assertEqual(len(umi.collapse_umis({"AAAA": 2, "AAAT": 2}, flavor="directional")), 2)

    def test_directional_uses_fixed_counts_and_bustools_root_order(self):
        counts = {"AAAA": 4, "AAAT": 1, "AATT": 1, "ATTT": 6}
        self.assertEqual({u: n for u, n, _ in umi.collapse_umis(counts, flavor="directional")},
                         {"AAAA": 6, "ATTT": 6})
        fixed = umi.collapse_umis({"AAAA": 3, "AAAT": 2, "AATA": 3}, flavor="directional")
        self.assertEqual(len(fixed), 2)

    def test_none_keeps_exact_duplicates_and_gene_collisions_distinct(self):
        products = list(umi.count_barcode([("AAAA", 0, 3), ("AAAA", 1, 2), ("AAAT", 0, 1)],
                                         {0: (0,), 1: (1,)}, ["A", "B"], umi_correction="none"))
        self.assertEqual(len(products), 3)
        self.assertEqual(sum(reads for _, _, reads, _ in products), 6)

    def test_strict_abundance_ties_and_nontransitive_chains(self):
        kernel = installed_parse_kernel()
        for counts in ({"AAAA": 5, "AAAT": 2, "AATT": 1},
                       {"AAAA": 2, "AAAT": 2}, {"AAAA": 4, "AAAT": 3, "AATA": 3}):
            result = umi.collapse_umis(counts, kernel=kernel)
            ordered = sorted(counts, key=lambda x: -counts[x])
            self.assertEqual({u: n for u, n, _ in result}, kernel(ordered, [counts[u] for u in ordered]))
            self.assertEqual(sum(n for _, n, _ in result), sum(counts.values()))
        self.assertEqual(len(umi.collapse_umis({"AAAA": 2, "AAAT": 2})), 2)
        self.assertEqual(len(umi.collapse_umis({"AAAA": 5, "AAAT": 2, "AATT": 1})), 2)

    def test_randomized_parity_in_both_parse_neighbor_branches(self):
        rng = random.Random(2701)
        kernel = installed_parse_kernel()
        for size in (1, 7, 50, 51, 100):
            for _ in range(8):
                counts = {}
                while len(counts) < size:
                    counts["".join(rng.choices("ACGT", k=6))] = rng.randint(1, 20)
                umi.collapse_umis(counts, kernel=kernel)

    def test_gene_assignment_precedes_collapse_and_umi_collisions_survive(self):
        classes = {0: (0,), 1: (1,), 2: (0, 1)}
        rows = [("AAAA", 0, 4), ("AAAA", 1, 3), ("AAAA", 2, 2)]
        products = list(umi.count_barcode(rows, classes, ["A", "B"], fractional=True))
        self.assertEqual({genes for _, genes, _, _ in products}, {("A",), ("B",)})
        self.assertEqual(len(products), 2)

    def test_ambiguous_support_does_not_double_count_anchored_molecule(self):
        rows = [("AAAA", 0, 3), ("AAAA", 1, 2)]
        products = list(umi.count_barcode(rows, {0: (0,), 1: (0, 1)}, ["A", "B"], fractional=True))
        self.assertEqual(products, [("AAAA", ("A",), 5, (0,))])

    def test_read_threshold_is_applied_after_umi_collapse(self):
        products = list(umi.count_barcode([("AAAA", 0, 3), ("AAAT", 0, 2)], {0: (0,)}, ["A"], min_reads=5))
        self.assertEqual(len(products), 1)


class NativeBusTests(unittest.TestCase):
    def test_parallel_matches_serial_for_all_correction_and_assignment_modes(self):
        rng = random.Random(72)
        records = [(bc, umi.dna_decode(rng.randrange(256), 4), rng.randrange(3), rng.randrange(1, 5))
                   for bc in ("AACCGG", "TTCCGG") for _ in range(350)]
        for flavor in ("parse", "directional", "none"):
            for assignment in ("unique", "fractional"):
                for transcript_assignment in ("unique", "fractional"):
                    with self.subTest(flavor=flavor, assignment=assignment, tx=transcript_assignment), tempfile.TemporaryDirectory() as tmp:
                        products, texts = [], []
                        for threads in (1, 2):
                            root = Path(tmp) / str(threads)
                            root.mkdir()
                            result, _ = native_fixture(root, self.options(
                                umi_correction=flavor, gene_assignment=assignment,
                                transcript_assignment=transcript_assignment),
                                min_reads=3, records=records, threads=threads)
                            products.append(result)
                            texts.append((root / "tscp_assignment.csv").read_bytes())
                        self.assertEqual(texts[0], texts[1])
                        for index in (0, 2):
                            np.testing.assert_array_equal(products[0][index].toarray(), products[1][index].toarray())
                        for index in (1, 3, 4):
                            pd.testing.assert_frame_equal(products[0][index], products[1][index])
                        self.assertEqual(products[0][5], products[1][5])

    def test_parallel_preserves_thirds_across_many_rt_identities(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference = reference_fixture(root)
            features = pd.read_csv(reference / "features.tsv", sep="\t")
            features.loc[len(features)] = ["tx4", "C", "Gamma", "human", "transcript"]
            features.to_csv(reference / "features.tsv", sep="\t", index=False)
            (root / "transcripts.txt").write_text("tx1\ntx2\ntx3\ntx4\n")
            (root / "matrix.ec").write_text("0\t0,2,3\n")
            genes = pd.DataFrame([[g, g, "human"] for g in "ABC"], columns=["gene_id", "gene_name", "genome"])
            lookup = {umi.dna_decode(i, 6): (f"01_02_{i // 4:02d}", str(i), "R" if i % 2 else "T", str(i))
                      for i in range(32)}
            records = [(barcode, sequence, 0, n) for barcode in lookup
                       for sequence, n in (("AAAA", 1), ("AAAT", 1), ("AATT", 3), ("CCCC", 1), ("GGGG", 1))]
            write_bus(root / "output.bus", records)
            for flavor in ("parse", "directional", "none"):
                results, assignments = [], []
                for threads in (1, 4):
                    prefix = f"{flavor}-{threads}"
                    diagnostics = {}
                    result = count_native(root / "output.bus", root / "transcripts.txt", root / "matrix.ec", reference,
                        genes, lookup, 6, 4, root / (prefix + ".sqlite"), root / (prefix + ".csv"),
                        self.options(gene_assignment="fractional", umi_correction=flavor), threads=threads,
                        kernel=installed_parse_kernel(), diagnostics=diagnostics)
                    self.assertEqual(diagnostics["workers"], threads)
                    results.append(result)
                    assignments.append((root / (prefix + ".csv")).read_bytes())
                np.testing.assert_array_equal(results[0][0].data, results[1][0].data)
                np.testing.assert_array_equal(results[0][2].data, results[1][2].data)
                self.assertEqual(assignments[0], assignments[1])
                pd.testing.assert_frame_equal(results[0][4], results[1][4])

    def test_merge_never_writes_the_worker_read_spool(self):
        import sqlite3
        from scripts import run_kallisto_gex as native
        original_spool, readers = native.spool_bus, []
        def hold_reader(*args):
            connection, total = original_spool(*args)
            reader = sqlite3.connect(args[1])
            reader.execute("BEGIN")
            reader.execute("SELECT COUNT(*) FROM reads").fetchone()
            readers.append(reader)
            return connection, total
        try:
            with tempfile.TemporaryDirectory() as tmp, patch.object(native, "spool_bus", side_effect=hold_reader):
                result, _ = native_fixture(Path(tmp), self.options(), threads=2)
                self.assertEqual(result[0].sum(), 3)
        finally:
            for reader in readers:
                reader.close()

    def test_parallel_worker_error_propagates_and_removes_scratch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with patch("scripts.kallisto_count_parallel.count_barcode", side_effect=RuntimeError("worker failure")):
                with self.assertRaisesRegex(RuntimeError, "worker failure"):
                    native_fixture(root, self.options(), threads=2)
            self.assertFalse(list(root.glob("umi-shards-*")))

    def test_fractional_roundoff_cannot_drop_cells_at_whole_umi_cutoffs(self):
        matrix = sparse.csr_matrix([[np.nextafter(3., 0), 0], [1.5, 1.5]])
        totals = whole_molecule_totals(matrix)
        np.testing.assert_array_equal(totals, [3, 3])
        self.assertTrue(np.issubdtype(totals.dtype, np.integer))
        with self.assertRaisesRegex(ValueError, "conserve whole-molecule"):
            whole_molecule_totals(sparse.csr_matrix([[1., 1.1]]))

    def options(self, **overrides):
        result = {"gene_assignment": "unique", "output_level": "both",
                  "transcript_assignment": "fractional", "cell_umi_threshold": None}
        result.update(overrides)
        return result

    def test_barcode_indices_are_collapsed_after_molecule_counting(self):
        with tempfile.TemporaryDirectory() as tmp:
            result, _ = native_fixture(Path(tmp), self.options())
            matrix, _, tx, _, obs, stats = result
            np.testing.assert_array_equal(matrix.toarray(), [[3, 0]])
            np.testing.assert_allclose(tx.toarray(), [[2.5, 0.5, 0]])
            self.assertEqual(obs.index.tolist(), ["01_02_03"])
            self.assertEqual(stats["discarded_multigene_reads"], 2)
            self.assertEqual(obs.loc["01_02_03", "mread_count"], 12)

    def test_all_flavors_reach_gene_and_transcript_products(self):
        records = [("AACCGG", sequence, 0, n) for sequence, n in (("AAAA", 5), ("AAAT", 2), ("AATT", 1))]
        for flavor, expected in (("parse", 2), ("directional", 1), ("none", 3)):
            with self.subTest(flavor=flavor), tempfile.TemporaryDirectory() as tmp:
                result, _ = native_fixture(Path(tmp), self.options(umi_correction=flavor), records=records)
                self.assertEqual(result[0].sum(), expected)
                self.assertEqual(result[2].sum(), expected)

    def test_fractional_gene_molecules_conserve_total_weight(self):
        with tempfile.TemporaryDirectory() as tmp:
            result, _ = native_fixture(Path(tmp), self.options(gene_assignment="fractional"))
            np.testing.assert_allclose(result[0].toarray(), [[3.5, 0.5]])
            self.assertEqual(result[0].sum(), result[-1]["molecules"])
            self.assertEqual(result[-1]["multigene_molecules"], 1)

    def test_read_threshold_follows_collapse_but_assignment_retains_molecules(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result, _ = native_fixture(root, self.options(), min_reads=3)
            np.testing.assert_array_equal(result[0].toarray(), [[2, 0]])
            assignments = pd.read_csv(root / "tscp_assignment.csv")
            self.assertEqual(sorted(assignments["count"]), [2, 3, 7])
            self.assertEqual(result[-1]["molecules_below_read_threshold"], 1)

    def test_conflicting_transcript_support_remains_unresolved(self):
        products = list(umi.count_barcode([("AAAA", 0, 4), ("AAAT", 1, 2)],
                                          {0: (0,), 1: (1,)}, ["A", "A"]))
        self.assertEqual(products, [("AAAA", ("A",), 6, ())])

    def test_truncated_bus_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.bus"
            write_bus(path, [("AACCGG", "AAAA", 0, 1)])
            path.write_bytes(path.read_bytes()[:-1])
            with self.assertRaisesRegex(ValueError, "Truncated BUS record"):
                list(umi.bus_records(path, 6, 4))

    def test_fastq_adapter_preserves_primer_indices_and_reports_invalid_umis(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "barcode_head.fastq"
            source.write_text("@01_02_03__T__1_2_3__AA_CC_GG__AAAA__time\nACGT\n+\nIIII\n"
                              "@01_02_03__R__4_2_3__TT_CC_GG__AANA__time\nACGT\n+\nIIII\n")
            lookup, bc_len, umi_len, stats = umi.prepare_fastqs(source, root / "cdna", root / "tags")
            self.assertEqual((bc_len, umi_len), (6, 4))
            self.assertEqual(lookup["TTCCGG"][1:3], ("4_2_3", "R"))
            self.assertEqual(stats["invalid_umi_reads"], 1)
            self.assertIn("AACCGGAAAA", (root / "tags").read_text())


@unittest.skipUnless(importlib.util.find_spec("anndata"), "AnnData is required for DGE integration tests")
class DgeContractTests(unittest.TestCase):
    options = NativeBusTests.options
    def test_parse_string_inference_uses_legacy_h5ad_encoding(self):
        import anndata as ad
        with tempfile.TemporaryDirectory() as tmp, pd.option_context("future.infer_string", True):
            obs = pd.DataFrame({"sample": pd.array(["Control", None], dtype="string")},
                               index=pd.Index(["01_02_03", "01_02_04"], name="bc_wells"))
            features = pd.DataFrame({"feature_id": ["A"], "gene_id": ["A"],
                                     "gene_name": ["Alpha"], "genome": ["human"]})
            with ad.settings.override(allow_write_nullable_strings=False):
                write_product(Path(tmp), sparse.csr_matrix([[1], [2]]), features, obs, {})
            result = ad.read_h5ad(Path(tmp) / "anndata.h5ad")
            self.assertEqual(result.obs_names.tolist(), ["01_02_03", "01_02_04"])
            self.assertTrue(pd.isna(result.obs.iloc[1]["sample"]))

    def test_h5ad_and_mtx_preserve_fractional_counts(self):
        import anndata as ad
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            options = self.options(gene_assignment="fractional")
            result, _ = native_fixture(root, options)
            publish_samples(FakeParse(), root / "gex_C1", result, options, {"backend": "kallisto"})
            path = root / "gex_C1/Control/DGE_filtered"
            np.testing.assert_allclose(sio.mmread(path / "count_matrix.mtx").toarray(), [[3.5, 0.5]])
            adata = ad.read_h5ad(path / "anndata.h5ad")
            self.assertEqual(adata.var_names.tolist(), ["Alpha_human", "Beta_human"])
            np.testing.assert_allclose(adata.layers["counts"].toarray(), [[3.5, 0.5]])
            self.assertEqual(adata.obs.loc["01_02_03", "matrix_count"], 4)
            transcript = ad.read_h5ad(root / "gex_C1/Control/DGE_transcript_filtered/anndata.h5ad")
            self.assertEqual(transcript.var_names.tolist(), ["tx1", "tx2", "tx3"])

    def test_transcript_primary_product_and_gene_based_cell_calls(self):
        import anndata as ad
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            options = self.options(output_level="transcript")
            result, _ = native_fixture(root, options)
            publish_samples(FakeParse(), root / "gex_C1", result, options, {"backend": "kallisto"})
            directory = root / "gex_C1/Control/DGE_filtered"
            adata = ad.read_h5ad(directory / "anndata.h5ad")
            self.assertEqual(adata.var_names.tolist(), ["tx1", "tx2", "tx3"])
            self.assertEqual(adata.obs_names.tolist(), ["01_02_03"])
            self.assertTrue((directory / "all_transcripts.csv").is_file())
            self.assertFalse((directory / "all_genes.csv").exists())

    def test_combination_suffixes_fractional_counts_and_stale_outputs(self):
        import anndata as ad
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = root / "fixture"
            fixture.mkdir()
            options = self.options(gene_assignment="fractional")
            result, reference = native_fixture(fixture, options)
            config = {"paths": {"output_dir": str(root / "output")},
                      "gex_alignment": {"backend": "kallisto", "kallisto": dict(options, reference_dir=str(reference))},
                      "libraries": [{"sublibrary_name": name, "sample_name": "Control"} for name in ("C1", "C2")]}
            config_path = root / "config.yaml"
            for name in ("C1", "C2"):
                library = root / "output" / f"gex_{name}"
                publish_samples(FakeParse(), library, result, options, {"backend": "kallisto"})
                process = library / "process"
                (process / "barcode_data.csv").write_text("fixture\n")
                (process / "run_proc_def.json").write_text("{}\n")
                (process / "tscp_assignment.csv").write_text("fixture\n")
                (process / "kallisto").mkdir()
                (process / "kallisto/input.json").write_text("{}\n")
                contract.atomic_json(process / "kallisto.pending.json", {
                    "schema_version": 1, "library": name, "backend": "kallisto",
                    "settings_sha256": contract.settings_identity(config),
                    "reference_sha256": contract.sha256_file(reference / "reference.json"),
                    "barcode_data_sha256": contract.sha256_file(process / "barcode_data.csv"),
                    "provenance": {"backend": "kallisto"},
                })
                seal(config, config_path, name)
            destination = combine(config, config_path)
            contract.validate_combined(config, config_path)
            adata = ad.read_h5ad(destination / "Control/DGE_filtered/anndata.h5ad")
            self.assertEqual(adata.obs_names.tolist(), ["01_02_03__s1", "01_02_03__s2"])
            np.testing.assert_allclose(adata.X.toarray(), [[3.5, 0.5], [3.5, 0.5]])
            spec = importlib.util.spec_from_file_location("step5", ROOT / "scripts/05_build_anndata.py")
            step5 = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(step5)
            with patch.object(step5, "sc", ad):
                joined = step5.load_sample_gex(str(root / "output"), "Control", [
                    (0, config["libraries"][0]), (2, config["libraries"][1])])
            self.assertEqual(joined.obs_names.tolist(), ["0_01_02_03", "2_01_02_03"])
            guide_counts = sparse.csr_matrix([[2], [4]])
            step5._attach_guide_matrix(joined, guide_counts, joined.obs_names.tolist(), ["guide1"], "final", 1)
            np.testing.assert_array_equal(joined.obsm["gRNA_counts"].toarray(), [[2], [4]])
            config["paths"]["genome_dir"] = "a_different_parse_reference"
            with self.assertRaisesRegex(ValueError, "Parse genome configuration differs"):
                contract.validate_combined(config, config_path)
            del config["paths"]["genome_dir"]
            (destination / "Control/DGE_filtered/count_matrix.mtx").write_text("changed")
            with self.assertRaisesRegex(ValueError, "Changed combined"):
                contract.validate_combined(config, config_path)


class ReferenceBuilderTests(unittest.TestCase):
    def test_build_uses_private_tmp_and_pins_runtime_assets(self):
        from scripts.build_kallisto_reference import build
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, gtf = root / "genome.fa", root / "genes.gtf"
            fasta.write_text(">chr1\nACGT\n")
            gtf.write_text('chr1\tx\tgene\t1\t4\t.\t+\t.\tgene_id "A"; gene_version "3"; gene_name "Alpha";\n')
            args = types.SimpleNamespace(output=root / "reference", fasta=fasta, gtf=gtf,
                                         kb="kb", kallisto="kallisto", workflow="nac", threads=2, genome="human")
            def fake_run(command, **kwargs):
                if command[0] == "kb":
                    work = Path(command[command.index("-i") + 1]).parent
                    self.assertEqual(Path(command[command.index("--tmp") + 1]).parent, work)
                    (work / "index.idx").write_bytes(b"index")
                    (work / "t2g.txt").write_text("tx1\tA.3\nA.3\tA.3\n")
                    (work / "nascent.txt").write_text("A.3\n")
                    (work / "cdna.txt").write_text("tx1\n")
                    kwargs["stdout"].write("build completed\n")
                return types.SimpleNamespace(returncode=0)
            with patch("scripts.build_kallisto_reference.subprocess.run", side_effect=fake_run), patch(
                    "scripts.build_kallisto_reference.subprocess.check_output", return_value="kallisto, version 0.52.0"):
                output = build(args)
            contract.load_reference(output)
            features = pd.read_csv(output / "features.tsv", sep="\t")
            self.assertEqual(features["feature_type"].tolist(), ["transcript", "nascent"])
            self.assertEqual(features["gene_id"].tolist(), ["A", "A"])
            self.assertEqual((output / "t2g.txt").read_text(), "tx1\tA\nA.3\tA\n")
            self.assertFalse((output / "kb_tmp").exists())
            with self.assertRaises(FileExistsError):
                build(args)

    def test_failed_build_retains_diagnostic_without_publishing_reference(self):
        from scripts.build_kallisto_reference import build
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, gtf = root / "genome.fa", root / "genes.gtf"
            fasta.write_text("source")
            gtf.write_text("source")
            args = types.SimpleNamespace(output=root / "reference", fasta=fasta, gtf=gtf,
                                         kb="kb", kallisto="kallisto", workflow="standard", threads=1)
            def failed_run(command, **kwargs):
                kwargs["stdout"].write("specific index failure\n")
                return types.SimpleNamespace(returncode=3)
            with patch("scripts.build_kallisto_reference.subprocess.run", side_effect=failed_run):
                with self.assertRaisesRegex(RuntimeError, "exit code 3"):
                    build(args)
            self.assertFalse((root / "reference").exists())
            self.assertIn("specific index failure", (root / "reference.failed.log").read_text())


class BenchmarkMathTests(unittest.TestCase):
    def test_sparse_expression_correlation_matches_dense_including_zeros(self):
        from scripts.benchmark_gex_report import sparse_cell_correlations
        first = np.array([[0, 2, 0, 5], [1, 0, 4, 0], [0, 0, 0, 0]])
        second = np.array([[0, 1, 3, 6], [2, 0, 3, 0], [0, 0, 0, 0]])
        actual = sparse_cell_correlations(sparse.csr_matrix(first), sparse.csr_matrix(second))
        expected = [np.corrcoef(np.log1p(a), np.log1p(b))[0, 1] for a, b in zip(first[:2], second[:2])]
        np.testing.assert_allclose(actual[:2], expected)
        self.assertTrue(np.isnan(actual[2]))

    def test_star_time_includes_index_loading_and_handles_midnight(self):
        from scripts.benchmark_gex_report import star_alignment_seconds
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "Log.final.out"
            path.write_text("Started job on | Oct 22 23:59:30\nStarted mapping on | Oct 22 23:59:50\nFinished on | Oct 23 00:00:30\n")
            self.assertEqual(star_alignment_seconds(path), 60)


class ConfigurationTests(unittest.TestCase):
    def test_common_umi_option_is_valid_for_star_and_kallisto(self):
        for backend in ("star", "kallisto"):
            for flavor in ("parse", "directional", "none"):
                options = contract.backend_options({"gex_alignment": {
                    "backend": backend, "umi_correction": flavor, "kallisto": {"reference_dir": "ref"}}})
                self.assertEqual(options["umi_correction"], flavor)
        with self.assertRaisesRegex(ValueError, "umi_correction"):
            contract.backend_options({"gex_alignment": {"umi_correction": "unsupported"}})

    def test_count_threads_is_a_nonnegative_integer(self):
        for value in (0, 1, 16):
            options = contract.backend_options({"gex_alignment": {"backend": "kallisto", "kallisto": {
                "reference_dir": "ref", "count_threads": value}}})
            self.assertEqual(options["count_threads"], value)
        for value in (-1, True, 1.5, "2"):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "count_threads"):
                contract.backend_options({"gex_alignment": {"backend": "kallisto", "kallisto": {
                    "reference_dir": "ref", "count_threads": value}}})

    def test_backend_mixing_is_rejected_even_without_star_bam(self):
        with tempfile.TemporaryDirectory() as tmp:
            process = Path(tmp)
            (process / "aw_uf_count_matrix.mtx").write_text("prior counts")
            with self.assertRaisesRegex(ValueError, "prior GEX counts"):
                contract.check_output_backend({"backend": "kallisto"}, process)
            (process / "kallisto").mkdir()
            contract.check_output_backend({"backend": "kallisto"}, process)
            with self.assertRaisesRegex(ValueError, "kallisto products"):
                contract.check_output_backend({"backend": "star"}, process)

    def test_star_is_default_and_options_are_validated(self):
        defaults = contract.backend_options({})
        self.assertEqual((defaults["backend"], defaults["umi_correction"], defaults["gene_assignment"]),
                         ("star", "parse", "unique"))
        options = contract.backend_options({"gex_alignment": {"backend": "kallisto", "kallisto": {"reference_dir": "ref"}}})
        self.assertEqual((options["output_level"], options["gene_assignment"]), ("gene", "unique"))
        with self.assertRaises(ValueError):
            contract.backend_options({"gex_alignment": {"backend": "unknown"}})
        with self.assertRaisesRegex(ValueError, "existing STAR"):
            contract.backend_options({"gex_alignment": {"backend": "kallisto", "kallisto": {"reference_dir": "ref"}},
                                      "gex_decontamination": {"existing_alignment": {"enabled": True}}})

    def test_reference_checksum_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            reference = reference_fixture(Path(tmp))
            contract.load_reference(reference)
            (reference / "index.idx").write_bytes(b"different index")
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                contract.load_reference(reference)


class StarUmiPolicyTests(unittest.TestCase):
    def test_policy_changes_require_fresh_star_output_and_recombination(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config_path = root / "config.yaml"
            config = {"paths": {"output_dir": str(root)}, "gex_alignment": {"backend": "star", "umi_correction": "directional"},
                      "libraries": [{"sublibrary_name": "C1", "sample_name": "Control"}]}
            process = root / "gex_C1/process"
            process.mkdir(parents=True)
            for name in ("aw_uf_count_matrix.mtx", "tscp_assignment.csv.gz", "barcode_data.csv", "run_proc_def.json"):
                (process / name).write_text("fixture")
            with self.assertRaisesRegex(ValueError, "another UMI correction policy"):
                contract.check_output_backend(contract.backend_options(config), process)
            contract.write_star_policy(config, process.parent)
            contract.check_output_backend(contract.backend_options(config), process)
            contract.validate_star_policies(config, config_path)
            # Existing Step 1b may update sample names without changing molecules.
            (process / "run_proc_def.json").write_text("renamed sample metadata")
            contract.validate_star_policies(config, config_path)
            destination = root / "gex_combined/Control/DGE_filtered"
            destination.mkdir(parents=True)
            for name in ("anndata.h5ad", "count_matrix.mtx", "cell_metadata.csv"):
                (destination / name).write_text("combined fixture")
            contract.write_star_combined_policy(config, config_path)
            contract.validate_star_combined_policy(config, config_path)
            config["gex_alignment"]["umi_correction"] = "none"
            with self.assertRaisesRegex(ValueError, "policy changed"):
                contract.validate_star_combined_policy(config, config_path)
            config["gex_alignment"]["umi_correction"] = "directional"
            (process / "aw_uf_count_matrix.mtx").write_text("changed")
            with self.assertRaisesRegex(ValueError, "molecule output changed"):
                contract.validate_star_policies(config, config_path)
            contract.write_star_policy(config, process.parent)
            with self.assertRaisesRegex(ValueError, "inputs changed since Step 1b"):
                contract.validate_star_combined_policy(config, config_path)

    def test_legacy_parse_output_remains_usable_but_custom_policy_requires_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = {"paths": {"output_dir": tmp}, "libraries": [{"sublibrary_name": "C1"}]}
            contract.validate_star_combined_policy(config, Path(tmp) / "config.yaml")
            config["gex_alignment"] = {"backend": "star", "umi_correction": "none"}
            with self.assertRaisesRegex(ValueError, "Missing UMI correction provenance"):
                contract.validate_star_combined_policy(config, Path(tmp) / "config.yaml")


@unittest.skipUnless(os.environ.get("CSM_TEST_BUSTOOLS"), "Set CSM_TEST_BUSTOOLS for directional graph parity")
class RealBustoolsTests(unittest.TestCase):
    def test_directional_matches_bustools_for_gene_resolved_evidence(self):
        executable = os.environ["CSM_TEST_BUSTOOLS"]
        rng = random.Random(5804)
        cases = [{"AAAAAA": 5, "AAAAAT": 2, "AAAATT": 1},
                 {"AAAAAA": 4, "AAAAAT": 1, "AAAATT": 1, "AAATTT": 6},
                 {"AAAAAA": 1, "AAAAAT": 1, "AAAATT": 1}]
        for n in (5, 30, 80):
            for _ in range(6):
                sequences = {"AAAAAA"}
                while len(sequences) < n:
                    parent = rng.choice(sorted(sequences))
                    sequences.add(rng.choice(list(umi.single_mutations(parent))))
                cases.append({sequence: rng.randint(1, 15) for sequence in sorted(sequences)})
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "transcripts.txt").write_text("tx1\n")
            (root / "t2g.txt").write_text("tx1\tA\n")
            (root / "matrix.ec").write_text("0\t0\n")
            for counts in cases:
                write_bus(root / "input.bus", [("AACCGG", seq, 0, count) for seq, count in sorted(counts.items())], umi_len=6)
                subprocess.run([executable, "umicorrect", "-o", str(root / "output.bus"),
                                "-g", str(root / "t2g.txt"), "-e", str(root / "matrix.ec"),
                                "-t", str(root / "transcripts.txt"), str(root / "input.bus")], check=True, capture_output=True)
                expected = {sequence: count for _, sequence, _, count in umi.bus_records(root / "output.bus", 6, 6)}
                actual = {sequence: count for sequence, count, _ in umi.collapse_umis(counts, flavor="directional")}
                self.assertEqual(actual, expected)


@unittest.skipUnless(os.environ.get("CSM_TEST_KALLISTO"), "Set CSM_TEST_KALLISTO for real executable smoke test")
class RealKallistoTests(unittest.TestCase):
    def test_real_pseudoalignment_and_custom_counting(self):
        executable = os.environ["CSM_TEST_KALLISTO"]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference = reference_fixture(root)
            rng = random.Random(42)
            sequences = ["".join(rng.choices("ACGT", k=240)) for _ in range(3)]
            fasta = root / "targets.fa"
            fasta.write_text("".join(f">tx{i + 1}\n{seq}\n" for i, seq in enumerate(sequences)))
            subprocess.run([executable, "index", "-i", str(reference / "index.idx"), str(fasta)],
                           check=True, capture_output=True)
            source = root / "barcode_head.fastq"
            with source.open("w") as handle:
                for index, sequence_umi in enumerate(["AAAA"] * 4 + ["AAAT"]):
                    seq = sequences[0][30:110]
                    handle.write(f"@01_02_03__T__1_2_3__AA_CC_GG__{sequence_umi}__r{index}\n{seq}\n+\n{'I' * len(seq)}\n")
            lookup, bc_len, umi_len, _ = umi.prepare_fastqs(source, root / "cdna.fq", root / "tag.fq")
            bus = root / "bus"
            subprocess.run([executable, "bus", "-i", str(reference / "index.idx"), "-o", str(bus),
                            "-x", "1,0,6:1,6,10:0,0,0", "--fr-stranded", "-t", "1",
                            str(root / "cdna.fq"), str(root / "tag.fq")], check=True, capture_output=True)
            genes = pd.DataFrame([["A", "Alpha", "human"], ["B", "Beta", "human"]],
                                 columns=["gene_id", "gene_name", "genome"])
            result = count_native(bus / "output.bus", bus / "transcripts.txt", bus / "matrix.ec", reference,
                                  genes, lookup, bc_len, umi_len, root / "reads.db", root / "molecules.csv",
                                  NativeBusTests().options(), kernel=installed_parse_kernel())
            np.testing.assert_array_equal(result[0].toarray(), [[1, 0]])
            self.assertEqual(result[-1]["pseudoaligned_reads"], 5)


if __name__ == "__main__":
    unittest.main()

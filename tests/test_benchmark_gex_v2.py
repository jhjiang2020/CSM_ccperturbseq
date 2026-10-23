"""Cache safety, scientific cohorts and the twelve-policy experiment graph."""
import importlib.util
import ast
from collections import defaultdict
import json
import os
import re
import shlex
import sqlite3
from pathlib import Path
import subprocess
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy import io, sparse
import yaml

from scripts import benchmark_gex_v2 as bench
from scripts import benchmark_gex_v2_report as reports
from scripts.gex_alignment import write_star_parameters

ROOT = Path(__file__).resolve().parents[1]


class BenchmarkV2Tests(unittest.TestCase):
    def test_code_seal_follows_local_dependencies_and_ignores_other_workflows(self):
        import hashlib
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            for name in ("benchmark_gex_v2.py", "benchmark_gex_v2_report.py", "splitpipe_gex_wrapper.py",
                         "write_fastq_input_identity.py", "gex_grna_preprocess_hook.py"):
                (directory / name).write_text("# entry point\n")
            (directory / "benchmark_gex_v2.py").write_text('from scripts.counter import count\nhelper = "worker.py"\n')
            (directory / "counter.py").write_text('from rules import policy\n')
            (directory / "rules.py").write_text('policy = 1\n')
            (directory / "worker.py").write_text('import counter\n')
            unrelated = directory / "03_match_grna_ibar.py"
            unrelated.write_text('workflow = 1\n')
            def digest():
                return hashlib.sha256(b''.join(p.read_bytes() for p in bench.benchmark_sources(directory))).hexdigest()
            original = digest()
            unrelated.write_text('workflow = 2\n')
            self.assertEqual(digest(), original)
            self.assertIn(directory / "worker.py", bench.benchmark_sources(directory))
            (directory / "rules.py").write_text('policy = 2\n')
            self.assertNotEqual(digest(), original)

    def test_read_spool_staging_uses_an_independent_local_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "shared.sqlite"
            with sqlite3.connect(source) as db:
                db.execute("CREATE TABLE reads (n INTEGER)")
                db.execute("INSERT INTO reads VALUES (1)")
            local = root / "local"
            local.mkdir()
            staged, details = bench.stage_read_database(source, local)
            self.assertFalse(staged.is_symlink())
            self.assertNotEqual(source.stat().st_ino, staged.stat().st_ino)
            self.assertEqual(source.read_bytes(), staged.read_bytes())
            self.assertEqual(details["cache_storage"], "job_local_copy_v2")
            with sqlite3.connect(staged) as db:
                db.execute("INSERT INTO reads VALUES (2)")
            with sqlite3.connect(source) as db:
                self.assertEqual(db.execute("SELECT SUM(n) FROM reads").fetchone()[0], 1)

    def test_retiming_never_archives_a_partial_or_live_counting_stage(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "manifests").mkdir()
            (root / "manifests/C1.kallisto_unique_parse.molecules.started").touch()
            args = types.SimpleNamespace(work_dir=root, repair_note="local cache staging")
            with self.assertRaisesRegex(ValueError, "Partial/live"):
                bench.archive_kallisto_counting(args)
            self.assertFalse((root / "timing_baselines/shared_sqlite_v1").exists())

    def test_retiming_requires_success_of_its_archive_dependency(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = types.SimpleNamespace(work_dir=Path(tmp), library="C1", arm="kallisto_unique_parse", require_count_baseline=True)
            with self.assertRaisesRegex(ValueError, "archive completion"):
                bench.run_arm(args)

    def test_archive_scope_preserves_unrelated_local_counting_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "manifests").mkdir()
            selected = root / "kallisto_unique_parse/gex_C1"
            selected.mkdir(parents=True)
            other = root / "kallisto_unique_none/gex_C1"
            other.mkdir(parents=True)
            (root / "manifests/C1.kallisto_unique_none.molecules.started").touch()
            (root / "experiment.json").write_text(json.dumps({"code_sha256": "revision"}))
            args = types.SimpleNamespace(work_dir=root, repair_note="retime shared storage", archive_arms=["kallisto_unique_parse"])
            bench.archive_kallisto_counting(args)
            self.assertTrue(other.exists())
            self.assertFalse(selected.exists())
            archive = json.loads((root / "timing_baselines/shared_sqlite_v1/archive.json").read_text())
            self.assertTrue(archive["complete"])
            self.assertEqual(archive["relocated_libraries"], [{"arm": "kallisto_unique_parse", "library": "C1"}])

    def test_retiming_equality_rejects_a_one_bit_count_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            arm = "kallisto_fractional_parse"
            baseline = work / "timing_baselines/shared_sqlite_v1"
            old_process = baseline / arm / "gex_C1/process"
            process = work / arm / "gex_C1/process"
            old_process.mkdir(parents=True)
            process.mkdir(parents=True)
            obs = pd.DataFrame({"tscp_count": [2]}, index=["01_01_01"])
            features = pd.DataFrame({"gene_id": ["A", "B"]})
            values = (sparse.csr_matrix([[1., 1.]]), features, None, features, obs, {"molecules": 2})
            bench.dump_pickle(old_process / "matrices.pickle", values)
            changed = list(values)
            changed[0] = values[0].copy()
            changed[0].data[0] = np.nextafter(1., 2.)
            bench.dump_pickle(process / "matrices.pickle", changed)
            for directory in (baseline / "manifests", work / "manifests"):
                directory.mkdir(parents=True)
                (directory / f"C1.{arm}.molecules.json").write_text(json.dumps({"artifacts": [
                    {"path": "/process/tscp_assignment.csv", "sha256": "identical"}]}))
            args = types.SimpleNamespace(work_dir=work, library="C1", arm=arm)
            with self.assertRaisesRegex(AssertionError, "scientific counts"):
                bench.verify_retimed_counts(args, process.parent)
            bench.dump_pickle(process / "matrices.pickle", values)
            bench.verify_retimed_counts(args, process.parent)
            proof = json.loads((process / "count_staging_equality.json").read_text())
            self.assertTrue(proof["exact_matrices_and_metadata"])
            self.assertTrue(proof["identical_molecule_assignments"])

    def test_vendor_star_command_has_one_multimap_limit_and_all_secondary_hits(self):
        vendor = ROOT / "ParseBiosciences-Pipeline.1.6.4/splitpipe/align.py"
        if not vendor.exists():
            self.skipTest("Parse installation is separate")
        tree = ast.parse(vendor.read_text())
        function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "get_STAR_align_call")
        scope = {"re": re}
        exec(compile(ast.Module(body=[function], type_ignores=[]), str(vendor), "exec"), scope)
        with tempfile.TemporaryDirectory() as tmp:
            for limit in (1, 3, 7):
                parameters = write_star_parameters({"gex_alignment": {"backend": "star", "star": {"max_multimaps": limit}}}, Path(tmp) / "star.txt")
                values = dict(line.split(None, 1) for line in parameters.read_text().splitlines())
                values.update(genome_dir="/reference", nthreads=16, align_star_linker_clip3p=False, align_star_keep_no_map=False)
                class Parse:
                    def get_exec(self, key):
                        return "STAR"
                    def filepath(self, key, *args, **kwargs):
                        return "/data/cdna.fastq" if key == "PF_FASTQ_BC" else "/data/barcode_head"
                    def get_par_val(self, key, **kwargs):
                        return values[key]
                tokens = shlex.split(scope["get_STAR_align_call"](Parse()))
                self.assertEqual(tokens.count("--outFilterMultimapNmax"), 1)
                self.assertEqual(tokens[tokens.index("--outFilterMultimapNmax") + 1], str(limit))
                self.assertEqual(tokens[tokens.index("--outSAMmultNmax") + 1], "-1")

    def test_twelve_policies_and_only_two_alignments(self):
        self.assertEqual(len(set(bench.ARMS)), 12)
        for arm in bench.ARMS:
            self.assertIn(bench.arm_options(arm)[0], bench.BACKENDS)
        stages = {s for arm in bench.ARMS for s in reports.arm_stages(arm)}
        self.assertEqual({s for s in stages if s.endswith(".align")}, {"star.align", "kallisto.align"})
        self.assertEqual({s for s in stages if s.endswith(".assign")},
                         {"star.unique.assign", "star.fractional.assign", "kallisto.assign"})
        self.assertEqual(len(reports.relaxation_edges()), 18)

    def test_fixed_cohort_alignment_fills_missing_cells_and_reorders_genes(self):
        product = (sparse.csr_matrix([[3., 1.], [7., 2.]]),
                   pd.DataFrame(index=["b", "a"]), pd.DataFrame({"gene_id": ["Y", "X"]}))
        matrix = reports.align_matrix(product, pd.Index(["a", "missing", "b"]), pd.Index(["X", "Y"]))
        np.testing.assert_array_equal(matrix.toarray(), [[2, 7], [0, 0], [1, 3]])
        gains = reports.gain_metrics([10, 20, 30], [12, 17, 30])
        self.assertEqual(gains["delta_umis"], -1)
        self.assertEqual((gains["cells_gain"], gains["cells_lose"], gains["cells_equal"]), (1, 1, 1))

    def test_cache_changes_rejected_and_mutable_metadata_copied(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source"
            source.mkdir()
            bam = source / "reads.bam"
            bam.write_bytes(b"immutable")
            metadata = source / "run_proc_def.json"
            metadata.write_text("{}")
            record = {"complete": True, "artifacts": [bench.fingerprint(bam)]}
            bench.validate_artifacts(record, full=True)
            bench.clone_process(source, root / "destination")
            self.assertEqual(bam.stat().st_ino, (root / "destination/reads.bam").stat().st_ino)
            self.assertNotEqual(metadata.stat().st_ino, (root / "destination/run_proc_def.json").stat().st_ino)
            bam.write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "Changed"):
                bench.validate_artifacts(record)

    @unittest.skipUnless(hasattr(os, "wait4"), "per-child resource accounting needs wait4")
    def test_per_stage_resource_accounting_is_not_cumulative_and_failures_not_sealed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "manifests").mkdir()
            args = types.SimpleNamespace(work_dir=root, library="C1", threads=2, verify_cache=True)
            command = [sys.executable, "-c", "sum(i*i for i in range(10000000))"]
            first = bench.stage(args, "first", command, [])
            second = bench.stage(args, "second", [sys.executable, "-c", "pass"], [])
            self.assertGreater(first["metrics"]["user_seconds"], second["metrics"]["user_seconds"])
            with self.assertRaises(subprocess.CalledProcessError):
                bench.stage(args, "failed", [sys.executable, "-c", "raise SystemExit(3)"], [])
            self.assertFalse((root / "manifests/C1.failed.json").exists())
            with self.assertRaises(FileExistsError):
                bench.stage(args, "failed", command, [])

    def test_native_gene_cache_requires_annotations(self):
        class Read:
            def has_tag(self, key):
                return key in self.tags
            def get_tag(self, key):
                return self.tags[key]
        read = Read()
        read.tags = {"GX": "GENE", "RE": 1}
        self.assertEqual(bench.cached_gene(read), ("GENE", 1.))
        read.tags = {}
        with self.assertRaises(ValueError):
            bench.cached_gene(read)

    @unittest.skipUnless(importlib.util.find_spec("pysam"), "pysam needed")
    def test_cached_native_annotation_retains_vendor_gene_ranking_and_gates(self):
        import pysam
        source = ROOT / "ParseBiosciences-Pipeline.1.6.4/splitpipe/molinfo.py"
        if not source.exists():
            self.skipTest("Parse 1.6.4 sources are separately installed")
        tree = ast.parse(source.read_text())
        functions = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in ("get_gene", "check_exon_alignment")]
        scope = {"np": np}
        exec(compile(ast.Module(body=functions, type_ignores=[]), str(source), "exec"), scope)
        info = {"gene_bins": defaultdict(list, {"chr1:0:+": ["noncoding", "coding"], "chr1:0:-": ["reverse"]}),
                "gtf_dict_stepsize": 1000, "gene_starts": {g: 0 for g in ("noncoding", "coding", "reverse")},
                "gene_ends": {g: 900 for g in ("noncoding", "coding", "reverse")},
                "gene_id_to_name": {g: g for g in ("noncoding", "coding", "reverse")},
                "gene_id_to_biotype": {"noncoding": "lncRNA", "coding": "protein_coding", "reverse": "protein_coding"},
                "genes_merged_exons": {"noncoding": [(0, 900)], "coding": [(0, 120)], "reverse": [(0, 900)]}}
        with tempfile.TemporaryDirectory() as tmp:
            first, second = Path(tmp) / "input.bam", Path(tmp) / "cache.bam"
            with pysam.AlignmentFile(first, "wb", header={"HD": {"SO": "queryname"}, "SQ": [{"SN": "chr1", "LN": 1000}]}) as output:
                for i, (flag, mapq, start) in enumerate(((0, 255, 100), (0, 255, 400), (16, 255, 100), (0, 3, 100), (256, 255, 100), (4, 0, -1))):
                    read = pysam.AlignedSegment()
                    read.query_name = f"read{i}"
                    read.query_sequence = "A" * 20
                    read.flag = flag
                    read.reference_id = -1 if flag == 4 else 0
                    read.reference_start = start
                    read.mapping_quality = mapq
                    if flag != 4:
                        read.cigarstring = "20M"
                    output.write(read)
            bench.annotate_unique(first, second, info, 255, scope["get_gene"])
            with pysam.AlignmentFile(second, "rb") as cached:
                results = [(r.query_name, bench.cached_gene(r)) for r in cached]
            self.assertEqual(results, [("read0", ("coding", 1.)), ("read1", ("coding", 0.)),
                ("read2", ("reverse", 1.)), ("read3", ("", 0.)), ("read5", ("", 0.))])

    def test_launcher_schedules_4_8_8_48_and_resource_profile(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = root / "config.yaml"
            config.write_text(yaml.safe_dump({"uger": {"kallisto_benchmark": {"slots": 8, "memory_per_slot_gb": 6, "runtime_hours": 5}}}))
            shim = root / "timeout"
            shim.write_text('#!/bin/bash\nshift\nexec "$@"\n')
            shim.chmod(0o755)
            env = dict(os.environ, PATH=f"{root}:{os.environ['PATH']}", CSM_UGER_PYTHON=sys.executable)
            result = subprocess.run(["bash", str(ROOT / "submit/submit_gex_benchmark_v2.sh"), "--config", str(config),
                "--work-dir", str(root / "run"), "--reference", str(root / "reference"), "--kallisto", str(root / "kallisto"), "--dry-run"],
                text=True, capture_output=True, env=env, check=True)
            lines = result.stdout.splitlines()
            self.assertEqual(len(lines), 7)
            for line in lines[:-1]:
                self.assertIn("-pe smp 8 -binding linear:8", line)
                self.assertIn("h_vmem=6G,h_rt=05:00:00", line)
            self.assertEqual(sum("-t 1-8" in line for line in lines), 2)
            self.assertEqual(sum("-t 1-48" in line for line in lines), 1)
            self.assertIn("-hold_jid csm_v2_bench_run", lines[-2])
            self.assertIn("-hold_jid csm_v2_bench_report", lines[-1])
            # Explicit profile settings override fallback publisher sizing too.
            self.assertIn("-pe smp 8 -binding linear:8", lines[-1])
            qsub = root / "qsub"
            qsub.write_text("#!/bin/bash\nprintf '100.1-4:1\\n'\n")
            qsub.chmod(0o755)
            actual = subprocess.run(["bash", str(ROOT / "submit/submit_gex_benchmark_v2.sh"), "--config", str(config),
                "--work-dir", str(root / "run"), "--reference", str(root / "reference"), "--kallisto", str(root / "kallisto")],
                text=True, capture_output=True, env=env, check=True)
            self.assertEqual(len(actual.stdout.splitlines()), 7)
            records = [json.loads(line) for line in (root / "run/logs/uger_submissions.jsonl").read_text().splitlines()]
            self.assertEqual(len(records), 7)
            for record in records[1:]:
                command = record["command"]
                self.assertEqual(command[command.index("-hold_jid") + 1], "100")

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "matplotlib needed")
    def test_report_all_pairs_fixed_cohort_and_shared_stage_accounting(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "manifests").mkdir()
            libraries = [{"sublibrary_name": name, "sample_name": "Sample"} for name in bench.LIBRARIES]
            (root / "experiment.json").write_text(json.dumps({"config": {"libraries": libraries}}))
            features = pd.DataFrame({"gene_id": ["A", "B", "C"]})
            for library in bench.LIBRARIES:
                stages = {s for arm in bench.ARMS for s in reports.arm_stages(arm)} | {"sample"}
                for stage in stages:
                    record = {"complete": True, "dependencies": {}, "artifacts": [], "metrics": {
                        "wall_seconds": 1., "cpu_hours": 0.01, "reserved_core_hours": .02, "max_process_rss_kb": 1000}}
                    (root / "manifests" / f"{library}.{stage}.json").write_text(json.dumps(record))
                for arm in bench.ARMS:
                    base = root / arm / f"gex_{library}"
                    obs = pd.DataFrame(index=pd.Index(["01_01_01", "01_01_02", "01_01_03"], name="bc_wells"))
                    matrix = sparse.csr_matrix([[1., 4., 2.], [6., 2., 3.], [3., 5., 7.]])
                    if arm.endswith("none"):
                        matrix = matrix * 2
                    for path, mtx, cells in ((base / "process", "aw_uf_count_matrix.mtx", "aw_uf_cell_metadata.csv"),
                                             (base / "Sample/DGE_filtered", "count_matrix.mtx", "cell_metadata.csv")):
                        path.mkdir(parents=True)
                        io.mmwrite(path / mtx, matrix)
                        obs.to_csv(path / cells)
                        features.to_csv(path / "all_genes.csv", index=False)
            args = types.SimpleNamespace(work_dir=root, pairs=50_000_000, verify_cache=True)
            with patch.object(reports, "plot"):
                reports.report(args)
            self.assertEqual(len(pd.read_csv(root / "comparison/agreement.csv")), 4 * 66)
            self.assertEqual(len(pd.read_csv(root / "comparison/library_summary.csv")), 48)
            gains = pd.read_csv(root / "comparison/umi_gains.csv")
            self.assertEqual(len(gains), 4 * 18)
            np.testing.assert_array_equal(gains[(gains["factor"] == "umi_policy") & gains["second"].str.endswith("none")].gain_percent, 100)
            runtime = pd.read_csv(root / "comparison/runtime.csv")
            np.testing.assert_array_equal(runtime.reconstructed_step1a_seconds, 6)
            np.testing.assert_array_equal(runtime.cached_policy_seconds, 2)
            resources = pd.read_csv(root / "comparison/stage_resources.csv")
            self.assertEqual(sum(resources.stage == "star.align"), 4)
            success = json.loads((root / "BENCHMARK_SUCCESS.json").read_text())
            self.assertEqual(success["configuration_pairs_per_library"], 66)
            (root / "manifests/C1.star_unique_parse.publish.json").unlink()
            with self.assertRaises(FileNotFoundError):
                reports.report(args)
            self.assertFalse((root / "BENCHMARK_SUCCESS.json").exists())


if __name__ == "__main__":
    unittest.main()

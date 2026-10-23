"""BAM ambiguity, shared policies, fractional products and provenance handoffs."""
import copy
import json
import os
import random
import sqlite3
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pysam

from scripts import gex_alignment as contract
from scripts.gex_counting import count_evidence
from scripts.gex_runtime import seal, run_parse_counting
from scripts.combine_kallisto_gex import combine
from scripts.run_star_fractional import candidate_genes, run_counting, spool_bam


def gene_info():
    return dict(all_genes=pd.DataFrame({"gene_id": ["A", "B"], "gene_name": ["Alpha", "Beta"],
                                       "genome": ["human", "human"]}),
                gtf_dict_stepsize=100, gene_bins={"chr1:0:+": ["A", "B"], "chr1:100:+": ["B"]},
                gene_starts={"A": 0, "B": 40}, gene_ends={"A": 89, "B": 199})


def bam_fixture(path, records=None, sort_order="queryname"):
    """Records are (UMI, barcode, original ID, [(position, flag, NH), ...])."""
    if records is None:
        records = [(umi, "AA_CC_GG", f"r{i}_{umi}", [(10, 0, 1)])
                   for umi, count in (("AAAA", 5), ("AAAT", 4)) for i in range(count)]
        records += [("CCCC", "AA_CC_GG", f"multi{i}", [(10, 0, 2), (110, 256, 2)]) for i in range(3)]
        records += [("GGGG", "AA_CC_GG", "beta", [(110, 0, 1)]),
                    ("TTTT", "AA_CC_GG", "same_gene", [(10, 0, 2), (20, 256, 2)]),
                    ("AAAA", "TT_CC_GG", "other_primer", [(10, 0, 1)])]
    with pysam.AlignmentFile(path, "wb", header={"HD": {"VN": "1.6", "SO": sort_order},
                                                 "SQ": [{"SN": "chr1", "LN": 1000}]}) as bam:
        for umi, barcode, original, hits in sorted(records, key=lambda r: (r[1], r[0], r[2])):
            bci, primer = ("1_2_3", "T") if barcode == "AA_CC_GG" else ("4_2_3", "R")
            for position, flag, nh in hits:
                read = pysam.AlignedSegment(bam.header)
                read.query_name = f"01_02_03__{primer}__{bci}__{barcode}__{umi}__stamp__OH_@{original}"
                read.query_sequence = "A" * 10
                read.query_qualities = pysam.qualitystring_to_array("I" * 10)
                read.flag = flag
                if not flag & 4:
                    read.reference_id = 0
                    read.reference_start = position
                    read.cigarstring = "10M"
                    read.mapping_quality = 255 if nh == 1 else 3
                    read.set_tag("NH", nh)
                bam.write(read)
    return path


def count_fixture(root, flavor="parse", threads=1, min_reads=1):
    options = contract.backend_options({"gex_alignment": {"backend": "star", "gene_assignment": "fractional",
                                                          "umi_correction": flavor}})
    bam = bam_fixture(root / "reads.bam")
    database = root / "reads.sqlite"
    conn, classes, features, lookup, stats = spool_bam(bam, database, gene_info())
    matrices = count_evidence(conn, classes, features, gene_info()["all_genes"], lookup,
                              database, root / "tscp_assignment.csv", options,
                              aligned_reads=stats["gene_compatible_reads"], threads=threads, min_reads=min_reads)
    return matrices, options


class Sample:
    def get_name(self):
        return "Control"

    def bc_df_for_sample(self, frame):
        return frame[frame["bc1_wind"] == "01"]


class ParseFixture:
    def get_par_val(self, key, **kwargs):
        return {"mol_awuf_read_thresh": 1, "mol_awuf_tscp_thresh": 1,
                "dge_species_purity_thresh": 0.9, "dge_tscp_cut_umin": 1}[key]

    def get_gene_info(self):
        return gene_info()

    def get_samples(self):
        return [Sample()]

    def get_tscp_cutoff(self, counts, sample):
        return True, {"cell_tscp_cutoff": 1, "cell_call": "fixture"}

    def get_bc_info(self):
        return types.SimpleNamespace(
            get_bc_data=lambda i: pd.DataFrame({"sequence": ["AA"], "well": ["A1"]}),
            get_seq_to_windzp=lambda i: {"AA": f"{i:02d}"})


class SharedPolicyTests(unittest.TestCase):
    def test_all_twelve_gene_policy_combinations_and_legacy_keys(self):
        for backend in ("star", "kallisto"):
            for assignment in ("unique", "fractional"):
                for correction in ("parse", "directional", "none"):
                    cfg = {"gex_alignment": {"backend": backend, "gene_assignment": assignment,
                                             "umi_correction": correction, "kallisto": {"reference_dir": "ref"}}}
                    options = contract.backend_options(cfg)
                    self.assertEqual((options["backend"], options["gene_assignment"], options["umi_correction"]),
                                     (backend, assignment, correction))
        old = {"gex_alignment": {"backend": "kallisto", "kallisto": {
            "reference_dir": "ref", "gene_assignment": "fractional", "output_level": "both"}}}
        self.assertEqual(contract.backend_options(old)["output_level"], "both")
        old["gex_alignment"]["gene_assignment"] = "unique"
        with self.assertRaisesRegex(ValueError, "Conflicting"):
            contract.backend_options(old)

    def test_star_rejects_transcript_output_and_invalid_policies(self):
        for level in ("transcript", "both"):
            with self.assertRaisesRegex(ValueError, "STAR supports gene"):
                contract.backend_options({"gex_alignment": {"output_level": level}})
        for maximum in (True, 0, -1, 1.5):
            with self.assertRaisesRegex(ValueError, "max_multimaps"):
                contract.backend_options({"gex_alignment": {"star": {"max_multimaps": maximum}}})

    def test_unique_star_rejects_reusing_old_alignment_with_changed_mapping_limit(self):
        with tempfile.TemporaryDirectory() as tmp:
            process = Path(tmp)
            (process / "aw_uf_count_matrix.mtx").write_text("historical matrix")
            options = contract.backend_options({"gex_alignment": {"star": {"max_multimaps": 10}}})
            with self.assertRaisesRegex(ValueError, "mapping limit changed"):
                contract.check_output_backend(options, process)


class StarEvidenceTests(unittest.TestCase):
    def test_input_spool_closes_on_reference_validation_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            database = root / "reads.sqlite"
            conn, classes, features, lookup, stats = spool_bam(
                bam_fixture(root / "reads.bam"), database, gene_info())
            features["gene_id"] = "missing-from-reference"
            options = contract.backend_options({"gex_alignment": {"gene_assignment": "fractional"}})
            with self.assertRaisesRegex(ValueError, "absent from Parse reference"):
                count_evidence(conn, classes, features, gene_info()["all_genes"], lookup,
                               database, root / "assignments.csv", options,
                               aligned_reads=stats["gene_compatible_reads"])
            with self.assertRaises(sqlite3.ProgrammingError):
                conn.execute("SELECT COUNT(*) FROM reads")

    def test_stranded_containment_retains_overlap_and_spliced_reads(self):
        with tempfile.TemporaryDirectory() as tmp:
            bam = bam_fixture(Path(tmp) / "one.bam", [("AAAA", "AA_CC_GG", "a", [(50, 0, 1)])])
            with pysam.AlignmentFile(bam) as handle:
                read = next(handle)
                self.assertEqual(candidate_genes(read, gene_info()), {"A", "B"})
                read.flag = 16
                self.assertEqual(candidate_genes(read, gene_info()), set())
                read.flag = 0
                read.reference_start = 80
                self.assertEqual(candidate_genes(read, gene_info()), {"A", "B"})
                read.reference_start = 85
                self.assertEqual(candidate_genes(read, gene_info()), {"B"})
                read.reference_start = 50
                read.cigarstring = "5M100N5M"
                self.assertEqual(candidate_genes(read, gene_info()), {"B"})

    def test_multilocus_same_gene_overlap_and_read_accounting(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            records = [("AAAA", "AA_CC_GG", "overlap", [(50, 0, 1)]),
                       ("CCCC", "AA_CC_GG", "same", [(10, 0, 2), (20, 256, 2)]),
                       ("GGGG", "AA_CC_GG", "reverse", [(10, 16, 1)]),
                       ("TTTT", "AA_CC_GG", "unmapped", [(-1, 4, 0)]),
                       ("NNNN", "AA_CC_GG", "invalid", [(10, 0, 1)])]
            conn, classes, features, lookup, stats = spool_bam(bam_fixture(root / "reads.bam", records),
                                                             root / "reads.sqlite", gene_info())
            self.assertEqual(set(classes.values()), {(0,), (0, 1)})
            self.assertEqual(conn.execute("SELECT SUM(n) FROM reads").fetchone()[0], 2)
            self.assertEqual((stats["input_reads"], stats["mapped_reads"], stats["invalid_umi_reads"]), (5, 4, 1))
            conn.close()

    def test_incomplete_or_invalid_bams_fail_before_counting(self):
        cases = [([(10, 0, 2)], "NH mismatch"),
                 ([(10, 0, 2), (110, 0, 2)], "one primary"),
                 ([(10, 0, 1), (110, 2048, 1)], "supplementary")]
        for hits, error in cases:
            with self.subTest(hits=hits), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                bam = bam_fixture(root / "reads.bam", [("AAAA", "AA_CC_GG", "bad", hits)])
                with self.assertRaisesRegex(ValueError, error):
                    spool_bam(bam, root / "reads.sqlite", gene_info())
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bam = bam_fixture(root / "reads.bam", sort_order="coordinate")
            with self.assertRaisesRegex(ValueError, "queryname-sorted"):
                spool_bam(bam, root / "reads.sqlite", gene_info())

    def test_multimap_limit_excludes_entire_reads(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            conn, _, _, _, stats = spool_bam(bam_fixture(root / "reads.bam"),
                                            root / "reads.sqlite", gene_info(), max_multimaps=1)
            self.assertEqual(stats["over_multimap_limit_reads"], 4)
            self.assertEqual(conn.execute("SELECT SUM(n) FROM reads").fetchone()[0], 11)
            conn.close()

    def test_all_umi_policies_conserve_fractional_counts_and_parallel_equality(self):
        for flavor, expected in (("parse", [3.5, 1.5]), ("directional", [4.5, 1.5]), ("none", [4.5, 1.5])):
            outputs = []
            for threads in (1, 2):
                with self.subTest(flavor=flavor, threads=threads), tempfile.TemporaryDirectory() as tmp:
                    root = Path(tmp)
                    matrices, _ = count_fixture(root, flavor, threads)
                    np.testing.assert_array_equal(matrices[0].toarray(), [expected])
                    self.assertEqual(matrices[4]["tscp_count"].iloc[0], sum(expected))
                    outputs.append((matrices[0].toarray(), (root / "tscp_assignment.csv").read_bytes()))
            np.testing.assert_array_equal(outputs[0][0], outputs[1][0])
            self.assertEqual(outputs[0][1], outputs[1][1])

    def test_read_cutoff_follows_umi_correction(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            matrices, _ = count_fixture(root, min_reads=6)
            np.testing.assert_array_equal(matrices[0].toarray(), [[1, 0]])
            assignments = pd.read_csv(root / "tscp_assignment.csv")
            self.assertTrue((assignments["count"] < 6).any())


@unittest.skipUnless(os.environ.get("CSM_TEST_STAR"), "Set CSM_TEST_STAR for real alignment smoke test")
class RealStarTests(unittest.TestCase):
    def test_real_secondary_alignments_feed_all_three_umi_policies(self):
        executable = os.environ["CSM_TEST_STAR"]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rng = random.Random(43129)
            genome = list("".join(rng.choices("ACGT", k=2000)))
            shared = "".join(rng.choices("ACGT", k=100))
            genome[100:200] = shared
            genome[700:800] = shared
            genome = "".join(genome)
            fasta = root / "genome.fa"
            fasta.write_text(f">chr1\n{genome}\n")
            index = root / "index"
            index.mkdir()
            subprocess.run([executable, "--runMode", "genomeGenerate", "--runThreadN", "1",
                            "--genomeDir", str(index), "--genomeFastaFiles", str(fasta),
                            "--genomeSAindexNbases", "3", "--genomeChrBinNbits", "8",
                            "--outFileNamePrefix", str(root / "build-")],
                           check=True, capture_output=True, text=True)
            fastq = root / "cdna.fastq"
            with fastq.open("w") as handle:
                for umi, sequence, count in (("AAAA", shared[10:90], 5),
                                             ("AAAT", shared[10:90], 4),
                                             ("CCCC", genome[230:310], 1)):
                    for read_index in range(count):
                        name = f"01_02_03__T__1_2_3__AA_CC_GG__{umi}__stamp__OH_@{umi}_{read_index}"
                        handle.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")
            prefix = str(root / "align-")
            subprocess.run([executable, "--genomeDir", str(index), "--runThreadN", "1",
                            "--readFilesIn", str(fastq), "--outFileNamePrefix", prefix,
                            "--outSAMtype", "BAM", "Unsorted", "--outFilterMultimapNmax", "3",
                            "--outSAMmultNmax", "-1"], check=True, capture_output=True, text=True)
            bam = root / "sorted.bam"
            pysam.sort("-n", "-o", str(bam), prefix + "Aligned.out.bam")
            with pysam.AlignmentFile(bam) as handle:
                alignments = list(handle)
            self.assertEqual(len(alignments), 19)
            self.assertEqual(sum(read.is_secondary for read in alignments), 9)
            self.assertEqual({read.get_tag("NH") for read in alignments}, {1, 2})
            info = gene_info()
            info.update(gene_starts={"A": 100, "B": 700}, gene_ends={"A": 399, "B": 999},
                        gene_bins={"chr1:100:+": ["A"], "chr1:200:+": ["A"], "chr1:700:+": ["B"]})
            for flavor, expected in (("parse", [1.5, 0.5]), ("directional", [2.0, 1.0]), ("none", [2.0, 1.0])):
                with self.subTest(flavor=flavor):
                    database = root / f"{flavor}.sqlite"
                    conn, classes, features, lookup, stats = spool_bam(bam, database, info)
                    self.assertEqual(stats["gene_compatible_reads"], 10)
                    options = contract.backend_options({"gex_alignment": {
                        "gene_assignment": "fractional", "umi_correction": flavor}})
                    matrices = count_evidence(conn, classes, features, info["all_genes"], lookup,
                                              database, root / f"{flavor}.csv", options,
                                              aligned_reads=stats["gene_compatible_reads"])
                    np.testing.assert_array_equal(matrices[0].toarray(), [expected])
            # Parse's native molinfo drops secondaries. Such a BAM must fail,
            # rather than silently treating all ambiguous reads as gene A.
            incomplete = root / "primary-only.bam"
            with pysam.AlignmentFile(bam) as source, pysam.AlignmentFile(incomplete, "wb", template=source) as target:
                for read in source:
                    if not read.is_secondary:
                        target.write(read)
            with self.assertRaisesRegex(ValueError, "NH mismatch"):
                spool_bam(incomplete, root / "incomplete.sqlite", info)


class StarHandoffTests(unittest.TestCase):
    def test_callback_failure_removes_old_completion_and_restores_parse(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = {"paths": {"output_dir": str(root), "splitpipe_dir": str(root / "parse")},
                      "gex_alignment": {"backend": "star", "gene_assignment": "fractional"},
                      "libraries": [{"sublibrary_name": "C1"}]}
            process = root / "gex_C1/process"
            process.mkdir(parents=True)
            marker = process / "star.complete.json"
            contract.atomic_json(marker, {"settings_sha256": contract.settings_identity(config)})
            original = lambda spipe: True
            package = types.ModuleType("splitpipe")
            package.dge = types.SimpleNamespace(run_dge=original)

            def fail(*args):
                raise RuntimeError("counting failed")

            previous_argv = sys.argv
            with patch.dict("sys.modules", {"splitpipe": package}), \
                 patch("scripts.gex_runtime.validate_splitpipe_tree"), \
                 patch("scripts.gex_runtime.runpy.run_path", side_effect=lambda *a, **k: package.dge.run_dge(None)):
                with self.assertRaisesRegex(RuntimeError, "counting failed"):
                    run_parse_counting(config, root / "config.yaml", "C1", 1, fail)
            self.assertFalse(marker.exists())
            self.assertIs(package.dge.run_dge, original)
            self.assertIs(sys.argv, previous_argv)

    def test_callback_seal_combine_h5ad_and_stale_inputs(self):
        import anndata as ad
        from scripts.kallisto_umi import collapse_umis, single_mutations
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference = root / "reference"
            reference.mkdir()
            for name in ("gene_info.json", "all_genes.csv"):
                (reference / name).write_text("fixture")
            config = {"paths": {"genome_dir": str(reference), "output_dir": str(root / "output")},
                      "gex_alignment": {"backend": "star", "gene_assignment": "fractional"},
                      "libraries": [{"sublibrary_name": n, "sample_name": "Control"} for n in ("C1", "C2")]}
            path = root / "config.yaml"
            kernel = lambda umis, counts: {u: n for u, n, _ in collapse_umis(dict(zip(umis, counts)))}
            package = types.ModuleType("splitpipe")
            package.cy_utils = types.SimpleNamespace(single_mut_seqs=single_mutations)
            molinfo = types.ModuleType("splitpipe.molinfo")
            molinfo.collapse_polyN = kernel
            for library in ("C1", "C2"):
                process = contract.library_output_dir(config, path, library) / "process"
                process.mkdir(parents=True)
                bam_fixture(process / "barcode_headAligned_sorted.bam")
                (process / "barcode_data.csv").write_text("fixture")
                (process / "pre_align_stats.csv").write_text("statistic,value\nreads_align_input,15\n")
                contract.atomic_json(process / "run_proc_def.json", {"genome": {"genome_dir": str(reference)}})
                with patch.dict("sys.modules", {"splitpipe": package, "splitpipe.molinfo": molinfo}):
                    self.assertTrue(run_counting(ParseFixture(), config, path, library, 1))
                seal(config, path, library)
                contract.validate_completion(config, path, library)
                contract.check_output_backend(contract.backend_options(config), process)
            destination = combine(config, path)
            contract.validate_combined(config, path)
            adata = ad.read_h5ad(destination / "Control/DGE_filtered/anndata.h5ad")
            self.assertEqual(adata.obs_names.tolist(), ["01_02_03__s1", "01_02_03__s2"])
            np.testing.assert_array_equal(adata.X.toarray(), [[3.5, 1.5], [3.5, 1.5]])
            np.testing.assert_array_equal(adata.layers["counts"].toarray(), adata.X.toarray())
            self.assertEqual(adata.uns["gex_quantification"]["backend"], "star")
            changed = copy.deepcopy(config)
            changed["gex_alignment"]["gene_assignment"] = "unique"
            with self.assertRaisesRegex(ValueError, "another gene assignment"):
                contract.check_output_backend(contract.backend_options(changed), process)
            changed["gex_alignment"]["gene_assignment"] = "fractional"
            changed["gex_alignment"]["umi_correction"] = "none"
            with self.assertRaisesRegex(ValueError, "another UMI correction"):
                contract.check_output_backend(contract.backend_options(changed), process)
            changed = copy.deepcopy(config)
            changed["gex_alignment"]["star"] = {"max_multimaps": 10}
            with self.assertRaisesRegex(ValueError, "mapping limit changed"):
                contract.check_output_backend(contract.backend_options(changed), process)
            contract.check_counting_settings(config, process)
            with self.assertRaisesRegex(ValueError, "counting settings changed"):
                contract.check_counting_settings(changed, process)
            with (process / "barcode_headAligned_sorted.bam").open("ab") as handle:
                handle.write(b"changed")
            with self.assertRaisesRegex(ValueError, "Changed STAR alignment"):
                contract.validate_completion(config, path, "C2")
            with self.assertRaisesRegex(ValueError, "Changed STAR alignment"):
                combine(config, path)


if __name__ == "__main__":
    unittest.main()

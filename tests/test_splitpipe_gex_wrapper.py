import hashlib
import sys
import tempfile
import unittest
from pathlib import Path

from scripts.splitpipe_gex_wrapper import (
    EXPECTED_SOURCE_SHA256,
    MODULE_TO_SOURCE,
    OverlayError,
    SplitpipeOverlayFinder,
    apply_hunks,
    install_runtime_overlay,
    parse_unified_diff,
    validate_splitpipe_tree,
)


ROOT = Path(__file__).resolve().parents[1]
SPLITPIPE_DIR = ROOT / "ParseBiosciences-Pipeline.1.6.4"
PATCH_FILE = ROOT / "patches" / "parsebio_splitpipe_1.6.4_gex_decontamination.patch"


class SplitpipeGexRuntimeOverlayTest(unittest.TestCase):
    def test_tracked_patch_transforms_both_pristine_modules(self):
        sections = parse_unified_diff(PATCH_FILE.read_text(encoding="utf-8"))
        self.assertEqual(set(MODULE_TO_SOURCE.values()), set(sections))

        transformed = {}
        for relative_path in MODULE_TO_SOURCE.values():
            source = (SPLITPIPE_DIR / relative_path).read_text(encoding="utf-8")
            transformed[relative_path] = apply_hunks(
                source, sections[relative_path], relative_path
            )
            compile(transformed[relative_path], relative_path, "exec")

        self.assertIn("_csm_load_gex_hook", transformed["splitpipe/preprocess.py"])
        self.assertIn(
            "classify_and_record_decision", transformed["splitpipe/preprocess.py"]
        )
        self.assertIn("resolve_survivor_trim", transformed["splitpipe/preprocess.py"])
        self.assertIn("record_survivor_trim", transformed["splitpipe/preprocess.py"])
        self.assertIn("_worker_error", transformed["splitpipe/preprocess.py"])
        self.assertIn("CSM_PRESERVE_HARMONIZED_BAM", transformed["splitpipe/molinfo.py"])

    def test_installing_overlay_does_not_write_vendor_sources(self):
        before = {
            path: hashlib.sha256((SPLITPIPE_DIR / path).read_bytes()).hexdigest()
            for path in MODULE_TO_SOURCE.values()
        }
        try:
            install_runtime_overlay(SPLITPIPE_DIR, PATCH_FILE)
        finally:
            sys.meta_path[:] = [
                finder
                for finder in sys.meta_path
                if not isinstance(finder, SplitpipeOverlayFinder)
            ]
        after = {
            path: hashlib.sha256((SPLITPIPE_DIR / path).read_bytes()).hexdigest()
            for path in MODULE_TO_SOURCE.values()
        }
        self.assertEqual(before, after)
        self.assertEqual(EXPECTED_SOURCE_SHA256, after)

    def test_checksum_guard_rejects_modified_parse_source(self):
        with tempfile.TemporaryDirectory() as temporary:
            splitpipe_dir = Path(temporary)
            package_dir = splitpipe_dir / "splitpipe"
            package_dir.mkdir()
            (splitpipe_dir / "split-pipe").write_text(
                'VER_NUM = "v1.6.4"\n', encoding="utf-8"
            )
            for relative_path in MODULE_TO_SOURCE.values():
                target = splitpipe_dir / relative_path
                target.write_bytes((SPLITPIPE_DIR / relative_path).read_bytes())
            (package_dir / "preprocess.py").write_text(
                (package_dir / "preprocess.py").read_text(encoding="utf-8")
                + "\n# local modification\n",
                encoding="utf-8",
            )

            with self.assertRaises(OverlayError):
                validate_splitpipe_tree(splitpipe_dir)


if __name__ == "__main__":
    unittest.main()

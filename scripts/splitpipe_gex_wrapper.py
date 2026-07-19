#!/usr/bin/env python3
"""Run Parse split-pipe with the CSM v1.4 modules overlaid in memory.

This wrapper never writes into the Parse distribution.  It installs an import
finder for the two modules touched by the tracked v1.4 unified diff, validates
the pristine Parse 1.6.4 sources, and executes ``split-pipe`` normally.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.abc
import importlib.util
import multiprocessing
import os
from pathlib import Path
import re
import runpy
import sys
from dataclasses import dataclass


EXPECTED_VERSION = "v1.6.4"
EXPECTED_SOURCE_SHA256 = {
    "splitpipe/preprocess.py": (
        "482901fb39c7ede7e100ad9bfdce889a8a960e203af0bfd94e9c21b217e60efa"
    ),
    "splitpipe/molinfo.py": (
        "fca58f2fb7aeb0d0e2d9a8f9c239e465093463432381280ca6edd6c4ab245fa7"
    ),
}
MODULE_TO_SOURCE = {
    "splitpipe.preprocess": "splitpipe/preprocess.py",
    "splitpipe.molinfo": "splitpipe/molinfo.py",
}
HUNK_RE = re.compile(
    r"^@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@(?:.*)\n?$"
)


class OverlayError(RuntimeError):
    """Raised when the runtime overlay cannot be applied safely."""


@dataclass(frozen=True)
class Hunk:
    old_start: int
    old_count: int
    new_start: int
    new_count: int
    body: tuple[str, ...]


def _normalise_diff_path(value: str) -> str:
    path = value.split("\t", 1)[0].strip()
    if path.startswith(("a/", "b/")):
        path = path[2:]
    return path


def parse_unified_diff(patch_text: str) -> dict[str, tuple[Hunk, ...]]:
    """Parse the file/hunk subset used by the tracked CSM unified diff."""
    lines = patch_text.splitlines(keepends=True)
    parsed: dict[str, list[Hunk]] = {}
    current_path: str | None = None
    index = 0

    while index < len(lines):
        line = lines[index]
        if line.startswith("--- "):
            if index + 1 >= len(lines) or not lines[index + 1].startswith("+++ "):
                raise OverlayError("malformed unified diff: missing +++ file header")
            current_path = _normalise_diff_path(lines[index + 1][4:])
            parsed.setdefault(current_path, [])
            index += 2
            continue

        if line.startswith("@@ "):
            if current_path is None:
                raise OverlayError("malformed unified diff: hunk has no file header")
            match = HUNK_RE.match(line)
            if match is None:
                raise OverlayError(f"malformed unified diff hunk: {line.rstrip()}")
            old_start, old_count, new_start, new_count = match.groups()
            body: list[str] = []
            index += 1
            while index < len(lines):
                body_line = lines[index]
                if body_line.startswith(("@@ ", "--- ")):
                    break
                if body_line.startswith("\\ No newline at end of file"):
                    index += 1
                    continue
                if not body_line.startswith((" ", "+", "-")):
                    raise OverlayError(
                        "malformed unified diff body line: " + body_line.rstrip()
                    )
                body.append(body_line)
                index += 1
            parsed[current_path].append(
                Hunk(
                    old_start=int(old_start),
                    old_count=int(old_count or "1"),
                    new_start=int(new_start),
                    new_count=int(new_count or "1"),
                    body=tuple(body),
                )
            )
            continue
        index += 1

    return {path: tuple(hunks) for path, hunks in parsed.items()}


def apply_hunks(source_text: str, hunks: tuple[Hunk, ...], source_name: str) -> str:
    """Apply hunks to text with exact old-side context validation."""
    source_lines = source_text.splitlines(keepends=True)
    line_delta = 0
    for hunk_number, hunk in enumerate(hunks, start=1):
        old_lines = [line[1:] for line in hunk.body if line[0] in " -"]
        new_lines = [line[1:] for line in hunk.body if line[0] in " +"]
        if len(old_lines) != hunk.old_count or len(new_lines) != hunk.new_count:
            raise OverlayError(
                f"{source_name}: inconsistent line counts in hunk {hunk_number}"
            )

        start = hunk.old_start - 1 + line_delta
        end = start + hunk.old_count
        if start < 0 or source_lines[start:end] != old_lines:
            # GNU patch permits a line-number offset when the complete context
            # is still exact. Do the same, but accept only one exact location.
            locations = [
                candidate
                for candidate in range(len(source_lines) - len(old_lines) + 1)
                if source_lines[candidate : candidate + len(old_lines)] == old_lines
            ]
            if len(locations) != 1:
                raise OverlayError(
                    f"{source_name}: source context for hunk {hunk_number} has "
                    f"{len(locations)} exact matches; refusing an unsafe runtime patch"
                )
            start = locations[0]
            end = start + hunk.old_count
        source_lines[start:end] = new_lines
        line_delta += hunk.new_count - hunk.old_count
    return "".join(source_lines)


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def validate_splitpipe_tree(splitpipe_dir: Path) -> None:
    splitpipe_bin = splitpipe_dir / "split-pipe"
    if not splitpipe_bin.is_file():
        raise OverlayError(f"split-pipe executable not found: {splitpipe_bin}")
    launcher = splitpipe_bin.read_text(encoding="utf-8")
    if f'VER_NUM = "{EXPECTED_VERSION}"' not in launcher:
        raise OverlayError(
            f"runtime overlay supports Parse split-pipe {EXPECTED_VERSION} only: "
            f"{splitpipe_bin}"
        )

    for relative_path, expected_hash in EXPECTED_SOURCE_SHA256.items():
        source_path = splitpipe_dir / relative_path
        if not source_path.is_file():
            raise OverlayError(f"required Parse source not found: {source_path}")
        observed_hash = _sha256(source_path.read_bytes())
        if observed_hash != expected_hash:
            detail = ""
            if "_csm_load_gex_hook" in source_path.read_text(
                encoding="utf-8", errors="replace"
            ):
                detail = " (the source appears to have the old physical CSM patch installed)"
            raise OverlayError(
                f"unexpected Parse 1.6.4 source checksum for {source_path}{detail}; "
                f"expected {expected_hash}, observed {observed_hash}"
            )


class _OverlayLoader(importlib.abc.Loader):
    def __init__(self, source_path: Path, hunks: tuple[Hunk, ...]):
        self.source_path = source_path
        self.hunks = hunks

    def create_module(self, spec):  # noqa: ANN001
        return None

    def exec_module(self, module) -> None:  # noqa: ANN001
        original = self.source_path.read_text(encoding="utf-8")
        overlaid = apply_hunks(original, self.hunks, str(self.source_path))
        module.__file__ = str(self.source_path)
        exec(compile(overlaid, str(self.source_path), "exec"), module.__dict__)


class SplitpipeOverlayFinder(importlib.abc.MetaPathFinder):
    """Intercept only the Parse modules named in ``MODULE_TO_SOURCE``."""

    def __init__(self, splitpipe_dir: Path, patch_sections: dict[str, tuple[Hunk, ...]]):
        self.splitpipe_dir = splitpipe_dir
        self.patch_sections = patch_sections

    def find_spec(self, fullname, path=None, target=None):  # noqa: ANN001
        relative_path = MODULE_TO_SOURCE.get(fullname)
        if relative_path is None:
            return None
        source_path = self.splitpipe_dir / relative_path
        loader = _OverlayLoader(source_path, self.patch_sections[relative_path])
        return importlib.util.spec_from_loader(
            fullname, loader, origin=str(source_path)
        )


def install_runtime_overlay(splitpipe_dir: Path, patch_file: Path) -> None:
    """Validate and install the import-time overlay without touching disk sources."""
    validate_splitpipe_tree(splitpipe_dir)
    if not patch_file.is_file():
        raise OverlayError(f"runtime overlay patch not found: {patch_file}")
    sections = parse_unified_diff(patch_file.read_text(encoding="utf-8"))
    missing = sorted(set(MODULE_TO_SOURCE.values()) - set(sections))
    if missing:
        raise OverlayError("runtime overlay patch is missing: " + ", ".join(missing))

    # Validate every transformation before Parse imports either module.
    for relative_path in MODULE_TO_SOURCE.values():
        source_path = splitpipe_dir / relative_path
        overlaid = apply_hunks(
            source_path.read_text(encoding="utf-8"),
            sections[relative_path],
            str(source_path),
        )
        compile(overlaid, str(source_path), "exec")

    if os.environ.get("CSM_GEX_DECONTAM_HOOK"):
        start_method = multiprocessing.get_start_method()
        if start_method != "fork":
            raise OverlayError(
                "the inline preprocessing overlay requires multiprocessing start "
                f"method 'fork', observed {start_method!r}; UGER/Linux uses 'fork'"
            )

    sys.meta_path.insert(0, SplitpipeOverlayFinder(splitpipe_dir, sections))


def _parse_wrapper_args(argv: list[str]) -> tuple[argparse.Namespace, list[str]]:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--wrapper_splitpipe_dir", required=True)
    parser.add_argument("--wrapper_patch", required=True)
    return parser.parse_known_args(argv)


def main(argv: list[str] | None = None) -> int:
    wrapper_args, splitpipe_args = _parse_wrapper_args(
        list(sys.argv[1:] if argv is None else argv)
    )
    splitpipe_dir = Path(wrapper_args.wrapper_splitpipe_dir).resolve()
    patch_file = Path(wrapper_args.wrapper_patch).resolve()

    sys.path.insert(0, str(splitpipe_dir))
    install_runtime_overlay(splitpipe_dir, patch_file)
    splitpipe_bin = splitpipe_dir / "split-pipe"
    sys.argv = [str(splitpipe_bin), *splitpipe_args]
    print(
        "CSM v1.4: applying GEX decontamination as an in-memory Parse 1.6.4 "
        "module overlay (vendor files unchanged)",
        flush=True,
    )
    try:
        runpy.run_path(str(splitpipe_bin), run_name="__main__")
    except SystemExit as exc:
        if exc.code is None:
            return 0
        if isinstance(exc.code, int):
            return exc.code
        print(exc.code, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except OverlayError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)

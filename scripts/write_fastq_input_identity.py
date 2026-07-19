#!/usr/bin/env python3
"""Write the exact ordered FASTQ identity consumed by one Step 1a task."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
from pathlib import Path


def pair_key(path: Path) -> str:
    name = path.name
    key, substitutions = re.subn(r"([._-])R[12]([._-])", r"\1R?\2", name, count=1)
    if substitutions != 1:
        raise ValueError(f"cannot identify R1/R2 token in FASTQ name: {path}")
    return key


def identity(path: Path) -> dict:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise ValueError(f"FASTQ does not exist: {resolved}")
    stat = resolved.stat()
    return {
        "path": str(resolved),
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }


def build_payload(library: str, r1_files: list[Path], r2_files: list[Path]) -> dict:
    if not r1_files or len(r1_files) != len(r2_files):
        raise ValueError("Step 1a requires the same nonzero number of R1 and R2 files")
    pairs = []
    digest = hashlib.sha256()
    for pair_index, (r1, r2) in enumerate(zip(r1_files, r2_files)):
        if pair_key(r1) != pair_key(r2):
            raise ValueError(f"FASTQ pair-name mismatch: {r1.name} versus {r2.name}")
        row = {
            "pair_index": pair_index,
            "r1": identity(r1),
            "r2": identity(r2),
        }
        pairs.append(row)
        digest.update(json.dumps(row, sort_keys=True, separators=(",", ":")).encode())
        digest.update(b"\n")
    return {
        "schema_version": 1,
        "library": library,
        "pair_count": len(pairs),
        "ordering": "numeric_version_sort_then_pair_index",
        "ordered_identity_sha256": digest.hexdigest(),
        "pairs": pairs,
    }


def atomic_write(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--r1", type=Path, action="append", required=True)
    parser.add_argument("--r2", type=Path, action="append", required=True)
    args = parser.parse_args()
    atomic_write(args.output, build_payload(args.library, args.r1, args.r2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Project-agnostic FASTQ-pair resolution for QC libraries."""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Dict, List, Sequence, Tuple


FASTQ_R1_RE = re.compile(r"_R1_(?P<suffix>.+\.f(?:ast)?q(?:\.gz)?)$")


def resolve_config_path(value: str, config_path: str) -> Path:
    candidate = Path(os.path.expandvars(os.path.expanduser(str(value))))
    if not candidate.is_absolute():
        candidate = Path(config_path).resolve().parent / candidate
    return candidate.resolve()


def prefix_is_token(filename: str, prefix: str) -> bool:
    """Require non-alphanumeric boundaries so C1 never selects C10."""
    return re.search(
        rf"(?<![A-Za-z0-9]){re.escape(prefix)}(?![A-Za-z0-9])",
        filename,
    ) is not None


def validate_pair(r1: Path, r2: Path) -> Tuple[str, str]:
    if not r1.is_file():
        raise FileNotFoundError(f"R1 FASTQ does not exist: {r1}")
    if not r2.is_file():
        raise FileNotFoundError(f"R2 FASTQ does not exist: {r2}")
    if r1.resolve() == r2.resolve():
        raise ValueError(f"R1 and R2 resolve to the same file: {r1}")
    return str(r1), str(r2)


def explicit_pairs(entry: dict, config_path: str) -> List[Tuple[str, str]]:
    values = entry.get("gex_fastq_pairs")
    if values is None:
        return []
    if not isinstance(values, list) or not values:
        raise ValueError("gex_fastq_pairs must be a non-empty list")
    pairs = []
    for index, value in enumerate(values, start=1):
        if isinstance(value, dict):
            if set(value) != {"r1", "r2"}:
                raise ValueError(
                    f"explicit pair {index} must contain exactly r1 and r2"
                )
            r1_value, r2_value = value["r1"], value["r2"]
        elif isinstance(value, (list, tuple)) and len(value) == 2:
            r1_value, r2_value = value
        else:
            raise ValueError(
                f"explicit pair {index} must be {{r1, r2}} or a two-item list"
            )
        pairs.append(
            validate_pair(
                resolve_config_path(str(r1_value), config_path),
                resolve_config_path(str(r2_value), config_path),
            )
        )
    if len(set(pairs)) != len(pairs):
        raise ValueError("gex_fastq_pairs contains a duplicate pair")
    return pairs


def discover_token_bounded_pairs(
    fastq_dir: str,
    prefix: str,
) -> List[Tuple[str, str]]:
    if not prefix:
        raise ValueError("gex_prefix must be non-empty when explicit pairs are absent")
    root = Path(fastq_dir).resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"FASTQ directory does not exist: {root}")
    r1_files = set()
    for pattern in (
        "*_R1_*.fastq.gz",
        "*_R1_*.fq.gz",
        "*_R1_*.fastq",
        "*_R1_*.fq",
    ):
        for candidate in root.rglob(pattern):
            match = FASTQ_R1_RE.search(candidate.name)
            if match and prefix_is_token(candidate.name[: match.start()], prefix):
                r1_files.add(candidate)
    pairs = []
    for r1 in sorted(r1_files):
        r2 = r1.with_name(r1.name.replace("_R1_", "_R2_", 1))
        pairs.append(validate_pair(r1, r2))
    if not pairs:
        raise FileNotFoundError(
            f"No token-bounded GEX FASTQ pairs for prefix {prefix!r} under {root}"
        )
    return pairs


def resolve_library_pairs(
    cfg: dict,
    config_path: str,
    library: dict,
) -> List[Tuple[str, str]]:
    pairs = explicit_pairs(library, config_path)
    if pairs:
        return pairs
    fastq_dir = resolve_config_path(cfg["paths"]["fastq_dir"], config_path)
    return discover_token_bounded_pairs(str(fastq_dir), str(library.get("gex_prefix", "")))


def resolve_all_library_pairs(
    cfg: dict,
    config_path: str,
) -> Dict[str, List[Tuple[str, str]]]:
    libraries = cfg.get("libraries", [])
    if not isinstance(libraries, list) or not libraries:
        raise ValueError("config must define at least one library")
    names = [str(entry.get("sublibrary_name", "")) for entry in libraries]
    if any(not name for name in names):
        raise ValueError("every library requires a non-empty sublibrary_name")
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError("duplicate sublibrary_name values: " + ", ".join(duplicates))
    resolved = {}
    owners = {}
    for entry, name in zip(libraries, names):
        pairs = resolve_library_pairs(cfg, config_path, entry)
        for pair in pairs:
            canonical = tuple(str(Path(value).resolve()) for value in pair)
            prior = owners.get(canonical)
            if prior is not None:
                raise ValueError(
                    f"FASTQ pair is assigned to both {prior!r} and {name!r}: {pair}"
                )
            owners[canonical] = name
        resolved[name] = pairs
    return resolved

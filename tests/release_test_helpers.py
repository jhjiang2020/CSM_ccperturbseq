"""Small, release-facing test helpers with no HPC/audit dependencies."""

from __future__ import annotations

import copy
import csv
import gzip
import hashlib
import json
from pathlib import Path

import yaml

from scripts import qc_gex_grna_contamination as qc
from scripts.qc_fastq_resolver import resolve_all_library_pairs


BENCHMARK_TOP_LEVEL_ARTIFACTS = {
    "cell_read_contamination.csv",
    "contaminant_guide_counts.csv",
    "contaminant_target_counts.csv",
    "contaminant_tso_offsets.csv",
    "per_cycle_base_composition.csv",
    "read_evidence.tsv.gz",
    "short54_filter_metrics.csv",
    "short54_metrics_by_tso_offset.csv",
}
ROW_ORDER_INDEPENDENT_ARTIFACTS = {
    "cell_read_contamination.csv",
    "contaminant_guide_counts.csv",
    "contaminant_target_counts.csv",
}


def artifact_semantic_digests(output_root: Path, library: str) -> dict:
    library_dir = output_root / library
    artifacts = {}
    for path in sorted(item for item in library_dir.rglob("*") if item.is_file()):
        relative = path.relative_to(library_dir).as_posix()
        if relative == "benchmark_summary.yaml":
            continue
        digest = hashlib.sha256()
        logical_bytes = 0
        opener = gzip.open if path.name.endswith(".gz") else open
        with opener(path, "rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
                logical_bytes += len(block)
        artifacts[relative] = {
            "sha256": digest.hexdigest(),
            "logical_size_bytes": logical_bytes,
        }
    return artifacts


def _selected_artifact_paths(root: Path, library: str) -> list[Path]:
    library_dir = root / library
    paths = [library_dir / name for name in sorted(BENCHMARK_TOP_LEVEL_ARTIFACTS)]
    paths.extend(sorted((library_dir / "synthetic_54").glob("*.fastq.gz")))
    alignment = library_dir / "alignment_input"
    paths.append(alignment / "ground_truth_manifest.tsv.gz")
    paths.extend(sorted((alignment / "full").glob("*.fastq.gz")))
    paths.extend(sorted((alignment / "short54").glob("*.fastq.gz")))
    missing = [path for path in paths if not path.is_file()]
    if missing:
        raise RuntimeError(f"release artifact set is incomplete: {missing}")
    return paths


def benchmark_artifact_digests(root: Path, library: str) -> dict:
    library_dir = root / library
    result = {}
    for path in _selected_artifact_paths(root, library):
        relative = path.relative_to(library_dir).as_posix()
        if relative in ROW_ORDER_INDEPENDENT_ARTIFACTS:
            with path.open(newline="") as handle:
                reader = csv.reader(handle)
                header = next(reader)
                rows = sorted(tuple(row) for row in reader)
            canonical = json.dumps(
                {"header": header, "rows": rows}, separators=(",", ":")
            ).encode()
            result[relative] = {
                "sha256": hashlib.sha256(canonical).hexdigest(),
                "logical_size_bytes": len(canonical),
                "comparison": "row_order_independent_csv",
                "row_count": len(rows),
            }
            continue
        digest = hashlib.sha256()
        logical_size = 0
        opener = gzip.open if path.name.endswith(".gz") else open
        with opener(path, "rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
                logical_size += len(block)
        result[relative] = {
            "sha256": digest.hexdigest(),
            "logical_size_bytes": logical_size,
            "comparison": "exact_logical_bytes",
        }
    return result


def write_fastq(path: Path, records: list[tuple[str, str]]) -> None:
    with gzip.open(path, "wt") as handle:
        for read_id, sequence in records:
            handle.write(f"@{read_id}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def write_pair(root: Path, stem: str, records: list[tuple[str, str]]):
    r1 = root / f"{stem}_R1_001.fastq.gz"
    r2 = root / f"{stem}_R2_001.fastq.gz"
    write_fastq(r1, records)
    write_fastq(r2, [(read_id, "A" * 58) for read_id, _ in records])
    return r1, r2


def resolver_generality(root: Path) -> dict:
    config_path = root / "resolver_config.yaml"
    one_dir = root / "one_library"
    one_dir.mkdir()
    for lane in (1, 2, 3):
        write_pair(one_dir, f"lane{lane}_ONLY", [(f"only_{lane}", "A" * 80)])
    one = resolve_all_library_pairs(
        {
            "paths": {"fastq_dir": str(one_dir)},
            "libraries": [{"sublibrary_name": "ONLY", "gex_prefix": "ONLY"}],
        },
        str(config_path),
    )

    four_dir = root / "four_libraries"
    four_dir.mkdir()
    c1 = write_pair(four_dir, "lane_C1", [("c1", "A" * 80)])
    write_pair(four_dir, "lane_C10", [("c10", "A" * 80)])
    write_pair(four_dir, "lane_C100", [("c100", "A" * 80)])
    explicit = write_pair(four_dir, "arbitrary_filename", [("explicit", "A" * 80)])
    four_cfg = {
        "paths": {"fastq_dir": str(four_dir)},
        "libraries": [
            {"sublibrary_name": "C1", "gex_prefix": "C1"},
            {"sublibrary_name": "C10", "gex_prefix": "C10"},
            {"sublibrary_name": "C100", "gex_prefix": "C100"},
            {
                "sublibrary_name": "EXPLICIT",
                "gex_fastq_pairs": [{"r1": str(explicit[0]), "r2": str(explicit[1])}],
            },
        ],
    }
    four = resolve_all_library_pairs(four_cfg, str(config_path))
    duplicate = copy.deepcopy(four_cfg)
    duplicate["libraries"][1]["sublibrary_name"] = "C1"
    try:
        resolve_all_library_pairs(duplicate, str(config_path))
    except ValueError:
        duplicate_rejected = True
    else:
        duplicate_rejected = False
    shared = copy.deepcopy(four_cfg)
    shared["libraries"][-1]["gex_fastq_pairs"] = [
        {"r1": str(c1[0]), "r2": str(c1[1])}
    ]
    try:
        resolve_all_library_pairs(shared, str(config_path))
    except ValueError:
        shared_rejected = True
    else:
        shared_rejected = False
    return {
        "one_library": {"library_count": 1, "pair_count": len(one["ONLY"])},
        "four_libraries": {
            "library_count": len(four),
            "overlapping_prefixes_isolated": len(four["C1"]) == 1,
            "explicit_pairs_supported": len(four["EXPLICIT"]) == 1,
        },
        "duplicate_names_rejected": duplicate_rejected,
        "shared_fastq_pairs_rejected": shared_rejected,
    }


def library_aggregate_case(*, root: Path, base_cfg: dict, library_count: int) -> dict:
    case_root = root / f"aggregate_{library_count}"
    case_root.mkdir()
    fastq_dir = case_root / "fastqs"
    fastq_dir.mkdir()
    output_dir = case_root / "output"
    guide = "ACGTACGTACGTACGTACGT"
    ibar = "ACGTACG"
    structure = base_cfg["grna_structure"]
    positive = (
        structure["fixed_seq_1"] + guide + structure["fixed_seq_2"]
        + ibar + structure["fixed_seq_4"] + "ACGT" * 20
    )
    names = (["L1", "L10", "L100"] + [f"LX{i:02d}" for i in range(1, 10)])[
        :library_count
    ]
    libraries = []
    for index, name in enumerate(names, start=1):
        token = "arbitrary_explicit" if library_count > 1 and index == len(names) else name
        lane1 = write_pair(fastq_dir, f"lane1_{token}", [(f"{name}_p", positive)])
        lane2 = write_pair(fastq_dir, f"lane2_{token}", [(f"{name}_n", "GCTA" * 45)])
        if token == "arbitrary_explicit":
            libraries.append({
                "sublibrary_name": name,
                "gex_fastq_pairs": [
                    {"r1": str(lane1[0]), "r2": str(lane1[1])},
                    {"r1": str(lane2[0]), "r2": str(lane2[1])},
                ],
            })
        else:
            libraries.append({"sublibrary_name": name, "gex_prefix": name})
    reference = case_root / "reference.csv"
    reference.write_text(
        "grna_seq,grna_name,target_gene,position,ibar_seq\n"
        f"{guide},guide_1,GENE1,1,{ibar}\n"
    )
    cfg = copy.deepcopy(base_cfg)
    cfg["paths"].update({
        "fastq_dir": str(fastq_dir),
        "output_dir": str(case_root / "splitpipe_output"),
        "grna_reference": str(reference),
    })
    cfg["libraries"] = libraries
    qc_cfg = cfg["gex_grna_contamination_qc"]
    qc_cfg["output_dir"] = str(output_dir)
    qc_cfg["raw_benchmark"].update({"max_reads": 0, "sample_fraction": 1.0, "read_table_max": 2})
    qc_cfg["synthetic_fastq"]["enabled"] = False
    qc_cfg["molecule_estimation"]["enabled"] = False
    qc_cfg["alignment_benchmark"]["enabled"] = False
    qc_cfg["bam_audit"]["enabled"] = False
    qc_cfg["guide_assignment_impact"]["enabled"] = False
    qc_cfg["plots"]["enabled"] = False
    cfg["gex_decontamination"]["tso_matching"]["backend"] = "myers"
    config_path = case_root / "config.yaml"
    config_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    resolved = resolve_all_library_pairs(cfg, str(config_path))
    for library in libraries:
        qc.benchmark_library(cfg, str(config_path), library["sublibrary_name"])
    qc.aggregate(cfg, str(config_path))
    with (output_dir / "all_libraries_contamination_summary.csv").open(newline="") as handle:
        summaries = list(csv.DictReader(handle))
    with (output_dir / "audit_completeness.csv").open(newline="") as handle:
        complete = list(csv.DictReader(handle))
    return {
        "library_count": len(resolved),
        "aggregate_summary_rows": len(summaries),
        "complete_aggregate_rows": sum(row["audit_complete"] == "True" for row in complete),
    }

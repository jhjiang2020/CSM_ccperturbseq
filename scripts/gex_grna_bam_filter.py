#!/usr/bin/env python3
"""Safely remove Step-1a-positive GEX reads from an existing split-pipe BAM.

The source BAM is opened read-only and is never replaced.  A filtered,
name-sorted BAM is written atomically below a separate sandbox together with
``process/csm_harmonized_bam.complete.json``.  Positive raw cluster IDs may be
read directly, selected from Step-1a decision rows, or regenerated from the
rejected R1 FASTQ named by a completed Step-1a manifest.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator, Optional, TextIO


PIPELINE_VERSION = "1.4"
OH_MARKER = "__OH_"
FILTER_MANIFEST_SCHEMA = 4
SURVIVOR_TRIM_POLICY = "full_tso_endpoint_then_parse_exact_17mer"


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open("r")


def normalize_raw_cluster_id(value: str) -> str:
    """Normalize a raw FASTQ ID or full split-pipe QNAME to its cluster token.

    split-pipe appends ``__OH_`` plus the original R1 header to its encoded
    QNAME. STAR retains the first whitespace-delimited word.  We therefore
    compare the whitespace-free raw token and tolerate the embedded ``@``.
    """
    token = str(value).strip()
    if OH_MARKER in token:
        token = token.rsplit(OH_MARKER, 1)[1]
    token = token.split()[0] if token else ""
    token = token.lstrip("@")
    if not token:
        raise ValueError(f"empty raw cluster ID derived from {value!r}")
    return token


def raw_cluster_from_splitpipe_qname(qname: str) -> Optional[str]:
    """Extract the normalized original cluster ID from a split-pipe QNAME."""
    if OH_MARKER not in qname:
        return None
    return normalize_raw_cluster_id(qname)


def _truthy(value) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


@dataclass(frozen=True)
class PositiveEvidence:
    ids: frozenset[str]
    record_count: int
    source_type: str
    source_path: str
    expected_positive_records: Optional[int] = None
    library: Optional[str] = None
    filter_run_id: Optional[str] = None
    reference_sha256: Optional[str] = None
    filter_settings_sha256: Optional[str] = None


def ids_from_fastq(path: Path) -> PositiveEvidence:
    ids = set()
    records = 0
    with open_text(path) as handle:
        while True:
            header = handle.readline().rstrip("\r\n")
            if not header:
                break
            sequence = handle.readline().rstrip("\r\n")
            plus = handle.readline().rstrip("\r\n")
            quality = handle.readline().rstrip("\r\n")
            if not sequence or not plus.startswith("+") or len(sequence) != len(quality):
                raise ValueError(f"malformed rejected FASTQ record {records + 1} in {path}")
            ids.add(normalize_raw_cluster_id(header))
            records += 1
    return PositiveEvidence(frozenset(ids), records, "rejected_r1_fastq", str(path))


def ids_from_plain_list(path: Path) -> PositiveEvidence:
    ids = set()
    records = 0
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            value = line.strip()
            if not value or value.startswith("#"):
                continue
            # A one-column file is canonical; accepting the first TSV column
            # makes scheduler-generated lists convenient without ambiguity.
            value = value.split("\t", 1)[0]
            try:
                ids.add(normalize_raw_cluster_id(value))
            except ValueError as exc:
                raise ValueError(f"invalid positive ID at {path}:{line_number}: {exc}") from exc
            records += 1
    if records == 0:
        raise ValueError(
            f"positive ID list is empty, so positive IDs are unavailable: {path}"
        )
    return PositiveEvidence(frozenset(ids), records, "positive_id_list", str(path))


def _decision_id(row: dict, path: Path, row_number: int) -> str:
    for key in ("read_id", "raw_read_id", "qname", "query_name", "header1"):
        value = row.get(key)
        if value not in (None, ""):
            return normalize_raw_cluster_id(value)
    raise ValueError(
        f"positive decision row {row_number} in {path} has no read ID column"
    )


def ids_from_decisions(path: Path) -> PositiveEvidence:
    """Load filter-positive IDs from JSONL(.gz) or headered TSV(.gz)."""
    ids = set()
    records = 0
    rows_seen = 0
    with open_text(path) as handle:
        first = handle.readline()
        if not first:
            raise ValueError(f"decision rows are empty; positive IDs unavailable: {path}")
        handle.seek(0)
        is_jsonl = first.lstrip().startswith("{")
        if is_jsonl:
            rows: Iterable[dict] = (
                json.loads(line) for line in handle if line.strip()
            )
        else:
            rows = csv.DictReader(handle, delimiter="\t")
        for rows_seen, row in enumerate(rows, start=1):
            if not isinstance(row, dict):
                raise ValueError(f"invalid decision row {rows_seen} in {path}")
            if not _truthy(row.get("filter_positive", False)):
                continue
            ids.add(_decision_id(row, path, rows_seen))
            records += 1
    if rows_seen == 0:
        raise ValueError(f"decision rows are empty; positive IDs unavailable: {path}")
    return PositiveEvidence(frozenset(ids), records, "filter_decisions", str(path))


def _resolve_manifest_path(manifest_path: Path, raw) -> Optional[Path]:
    if raw in (None, ""):
        return None
    path = Path(raw).expanduser()
    if not path.is_absolute():
        path = manifest_path.parent / path
    return path.resolve()


def ids_from_filter_manifest(path: Path) -> PositiveEvidence:
    path = path.expanduser().resolve()
    with path.open() as handle:
        manifest = json.load(handle)
    if not isinstance(manifest, dict):
        raise ValueError(f"filter manifest must be a JSON object: {path}")
    if manifest.get("complete") is not True or str(manifest.get("status", "")).lower() != "complete":
        raise ValueError(f"Step-1a filter manifest is not complete: {path}")

    statistics = manifest.get("statistics", {}) or {}
    # Filter-positive records are the removal truth. rejected_pair_count is
    # only equivalent when rejected FASTQ output was enabled; otherwise the
    # hook legitimately reports zero rejected pairs despite positive calls.
    expected = statistics.get("csm_gex_filter_positive_reads")
    if expected is None:
        expected = manifest.get("rejected_pair_count")
    expected = int(expected) if expected is not None else None

    # A completed zero-positive manifest is authoritative even when the hook
    # also emitted empty decision/FASTQ files.
    if expected == 0:
        return PositiveEvidence(
            frozenset(), 0, "completed_zero_positive_manifest", str(path), 0,
            manifest.get("library"), manifest.get("run_id"),
            manifest.get("reference_sha256"),
            manifest.get("filter_settings_sha256"),
        )

    positive_list = _resolve_manifest_path(
        path, manifest.get("positive_qnames") or manifest.get("positive_ids")
    )
    decisions = _resolve_manifest_path(path, manifest.get("filter_decisions"))
    rejected_r1 = _resolve_manifest_path(path, manifest.get("rejected_r1"))

    if positive_list and positive_list.is_file():
        evidence = ids_from_plain_list(positive_list)
    elif decisions and decisions.is_file():
        evidence = ids_from_decisions(decisions)
    elif rejected_r1 and rejected_r1.is_file():
        evidence = ids_from_fastq(rejected_r1)
    else:
        raise ValueError(
            "completed Step-1a manifest does not provide positive IDs: expected "
            "positive_qnames, filter_decisions, or rejected_r1"
        )

    if expected is not None and evidence.record_count != expected:
        raise ValueError(
            f"positive evidence has {evidence.record_count} records but manifest expects {expected}"
        )
    if expected is None and evidence.record_count == 0:
        raise ValueError(
            "completed Step-1a manifest has empty evidence but no zero-positive count; "
            "positive IDs are unavailable"
        )
    return PositiveEvidence(
        evidence.ids, evidence.record_count, evidence.source_type,
        evidence.source_path, expected, manifest.get("library"),
        manifest.get("run_id"), manifest.get("reference_sha256"),
        manifest.get("filter_settings_sha256"),
    )


def load_positive_evidence(*, filter_manifest=None, positive_ids=None,
                           decision_rows=None) -> PositiveEvidence:
    selected = [value for value in (filter_manifest, positive_ids, decision_rows) if value]
    if len(selected) != 1:
        raise ValueError(
            "provide exactly one of --filter-manifest, --positive-ids, or --decision-rows"
        )
    if filter_manifest:
        return ids_from_filter_manifest(Path(filter_manifest))
    if positive_ids:
        return ids_from_plain_list(Path(positive_ids).expanduser().resolve())
    return ids_from_decisions(Path(decision_rows).expanduser().resolve())


def _queryname_sort_key(name: str):
    """Approximate samtools natural query-name ordering for validation."""
    return tuple(
        (0, int(part)) if part.isdigit() else (1, part)
        for part in re.split(r"(\d+)", name)
        if part != ""
    )


@dataclass
class FilterStats:
    records_before: int = 0
    records_removed: int = 0
    records_after: int = 0
    qnames_before: int = 0
    qnames_removed: int = 0
    qnames_after: int = 0
    positive_qnames_seen: int = 0
    qnames_with_original_header: int = 0


class QnameFilter:
    """Stateful QNAME filter that also validates name-sort order."""

    def __init__(self, positive_ids: frozenset[str]):
        self.positive_ids = positive_ids
        self.stats = FilterStats()
        self.seen_positive = set()
        self.previous_qname = None
        self.previous_key = None

    def keep(self, qname: str) -> bool:
        sort_key = _queryname_sort_key(qname)
        if self.previous_key is not None and sort_key < self.previous_key:
            raise ValueError(
                f"source BAM is not query-name sorted: {qname!r} follows {self.previous_qname!r}"
            )
        new_qname = qname != self.previous_qname
        if new_qname:
            self.stats.qnames_before += 1
        self.previous_qname = qname
        self.previous_key = sort_key

        self.stats.records_before += 1
        raw_id = raw_cluster_from_splitpipe_qname(qname)
        if new_qname and raw_id is not None:
            self.stats.qnames_with_original_header += 1
        remove = raw_id is not None and raw_id in self.positive_ids
        if remove:
            self.stats.records_removed += 1
            if new_qname:
                self.stats.qnames_removed += 1
            self.seen_positive.add(raw_id)
            return False
        self.stats.records_after += 1
        if new_qname:
            self.stats.qnames_after += 1
        return True

    def finalize(self) -> FilterStats:
        if self.stats.records_before != self.stats.records_removed + self.stats.records_after:
            raise RuntimeError("BAM record invariant failed: before != removed + after")
        if self.stats.qnames_before != self.stats.qnames_removed + self.stats.qnames_after:
            raise RuntimeError("BAM QNAME invariant failed: before != removed + after")
        if self.stats.records_before and not self.stats.qnames_with_original_header:
            raise ValueError(
                "source BAM QNAMEs do not retain split-pipe's __OH_ original-header token"
            )
        self.stats.positive_qnames_seen = len(self.seen_positive)
        return self.stats


def filter_sam_lines(lines: Iterable[str], output: TextIO,
                     positive_ids: frozenset[str]):
    """Filter SAM text while preserving every header/survivor line exactly."""
    qfilter = QnameFilter(positive_ids)
    header_lines = []
    saw_record = False
    for line in lines:
        if not saw_record and line.startswith("@"):
            header_lines.append(line)
            output.write(line)
            continue
        saw_record = True
        if not line.strip():
            continue
        fields = line.split("\t", 1)
        if len(fields) != 2:
            raise ValueError("malformed SAM alignment line")
        if qfilter.keep(fields[0]):
            output.write(line)
    _validate_header_sort_text("".join(header_lines))
    return qfilter.finalize(), "".join(header_lines)


def _validate_header_sort_dict(header: dict) -> None:
    sort_order = str((header.get("HD") or {}).get("SO", "")).lower()
    if sort_order and sort_order not in {"queryname", "unknown", "unsorted"}:
        raise ValueError(
            f"source BAM header sort order is {sort_order!r}, not query-name order"
        )


def _validate_header_sort_text(header_text: str) -> None:
    for line in header_text.splitlines():
        if not line.startswith("@HD\t"):
            continue
        fields = dict(
            field.split(":", 1) for field in line.split("\t")[1:] if ":" in field
        )
        _validate_header_sort_dict({"HD": fields})
        break


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _filter_with_pysam(source: Path, temporary: Path,
                       positive_ids: frozenset[str]):
    import pysam

    qfilter = QnameFilter(positive_ids)
    with pysam.AlignmentFile(str(source), "rb") as source_bam:
        header_dict = source_bam.header.to_dict()
        header_text = str(source_bam.header)
        _validate_header_sort_dict(header_dict)
        with pysam.AlignmentFile(str(temporary), "wb", template=source_bam) as output:
            for record in source_bam.fetch(until_eof=True):
                if qfilter.keep(record.query_name):
                    output.write(record)
    stats = qfilter.finalize()

    with pysam.AlignmentFile(str(temporary), "rb") as check:
        if str(check.header) != header_text:
            raise RuntimeError("filtered BAM header differs from source BAM header")
        output_count = sum(1 for _ in check.fetch(until_eof=True))
    if output_count != stats.records_after:
        raise RuntimeError(
            f"filtered BAM contains {output_count} records, expected {stats.records_after}"
        )
    return stats, header_text


def _run_checked(command, **kwargs):
    result = subprocess.run(command, text=True, capture_output=True, **kwargs)
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(command)}\n{result.stderr.strip()}"
        )
    return result


def _filter_with_samtools(source: Path, temporary: Path,
                          positive_ids: frozenset[str], samtools: str):
    view_err = tempfile.TemporaryFile(mode="w+")
    write_err = tempfile.TemporaryFile(mode="w+")
    view = subprocess.Popen(
        [samtools, "view", "--no-PG", "-h", str(source)],
        stdout=subprocess.PIPE, stderr=view_err, text=True,
    )
    writer = subprocess.Popen(
        [samtools, "view", "--no-PG", "-b", "-o", str(temporary), "-"],
        stdin=subprocess.PIPE, stderr=write_err, text=True,
    )
    try:
        if view.stdout is None or writer.stdin is None:
            raise RuntimeError("failed to open samtools streaming pipes")
        stats, header_text = filter_sam_lines(view.stdout, writer.stdin, positive_ids)
        writer.stdin.close()
        view_rc = view.wait()
        writer_rc = writer.wait()
        if view_rc != 0 or writer_rc != 0:
            view_err.seek(0)
            write_err.seek(0)
            raise RuntimeError(
                f"samtools streaming failed: reader={view_rc} writer={writer_rc}; "
                f"reader stderr={view_err.read().strip()!r}; "
                f"writer stderr={write_err.read().strip()!r}"
            )
    except Exception:
        view.terminate()
        writer.terminate()
        view.wait()
        writer.wait()
        raise
    finally:
        view_err.close()
        write_err.close()

    check_header = _run_checked(
        [samtools, "view", "--no-PG", "-H", str(temporary)]
    ).stdout
    if check_header != header_text:
        raise RuntimeError("filtered BAM header differs from source BAM header")
    output_count = int(_run_checked([samtools, "view", "-c", str(temporary)]).stdout.strip())
    if output_count != stats.records_after:
        raise RuntimeError(
            f"filtered BAM contains {output_count} records, expected {stats.records_after}"
        )
    return stats, header_text


def atomic_json(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def register_fresh_harmonized_bam(source_bam, output_dir, filter_manifest,
                                  *, library=None):
    """Register the BAM produced after pre-STAR filtering in a fresh run."""
    source = Path(source_bam).expanduser().resolve()
    sandbox = Path(output_dir).expanduser().resolve()
    filter_path = Path(filter_manifest).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"fresh harmonized BAM does not exist: {source}")
    expected = (sandbox / "process" / "barcode_headAligned_sorted.bam").resolve()
    if source != expected:
        raise ValueError(
            f"fresh BAM must use the split-pipe sorted-BAM contract: {expected}"
        )
    with filter_path.open() as handle:
        filter_payload = json.load(handle)
    if filter_payload.get("complete") is not True:
        raise ValueError(f"Step-1a filter manifest is incomplete: {filter_path}")
    if filter_payload.get("schema_version") != FILTER_MANIFEST_SCHEMA:
        raise ValueError(f"Step-1a filter manifest schema is stale: {filter_path}")
    if filter_payload.get("survivor_tso_trim_policy") != SURVIVOR_TRIM_POLICY:
        raise ValueError(f"Step-1a survivor TSO trim policy is stale: {filter_path}")
    filter_library = filter_payload.get("library")
    if library and filter_library and library != filter_library:
        raise ValueError(
            f"--library {library!r} does not match filter manifest "
            f"library {filter_library!r}"
        )

    manifest_path = sandbox / "process" / "csm_harmonized_bam.complete.json"
    bam_sha256 = _sha256(source)
    filter_sha256 = _sha256(filter_path)
    if manifest_path.exists():
        existing = json.loads(manifest_path.read_text())
        if (
            existing.get("complete") is True
            and existing.get("alignment_reused") is False
            and existing.get("output_bam_sha256") == bam_sha256
            and existing.get("filter_manifest_sha256") == filter_sha256
        ):
            return existing, manifest_path
        # A changed filter manifest makes split-pipe preprocessing set
        # reuse=false, so the BAM reaching this point was rebuilt and may
        # atomically supersede the stale completion record.

    stats = filter_payload.get("statistics", {}) or {}
    payload = {
        "schema_version": 1,
        "pipeline_version": PIPELINE_VERSION,
        "status": "complete",
        "complete": True,
        "alignment_reused": False,
        "library": library or filter_library,
        "output_bam": str(source),
        "output_bam_sha256": bam_sha256,
        "source_bam_kind": "fresh_name_sorted",
        "filter_manifest": str(filter_path),
        "filter_manifest_sha256": filter_sha256,
        "filter_run_id": filter_payload.get("run_id"),
        "filter_reference_sha256": filter_payload.get("reference_sha256"),
        "filter_settings_sha256": filter_payload.get("filter_settings_sha256"),
        "survivor_tso_trim_policy": filter_payload.get(
            "survivor_tso_trim_policy"
        ),
        "filter_positive_reads": int(
            stats.get("csm_gex_filter_positive_reads", 0)
        ),
        "valid_barcode_reads_removed_before_star": int(
            stats.get("csm_gex_filter_positive_valid_barcode_reads", 0)
        ),
        "completed_at": datetime.now(timezone.utc).isoformat(),
    }
    atomic_json(manifest_path, payload)
    return payload, manifest_path


def filter_existing_bam(source_bam, output_dir, evidence: PositiveEvidence,
                        *, backend="auto", samtools=None, library=None):
    source = Path(source_bam).expanduser().resolve()
    sandbox = Path(output_dir).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"source BAM does not exist: {source}")
    process_dir = sandbox / "process"
    output = process_dir / "barcode_headAligned_sorted.bam"
    manifest_path = process_dir / "csm_harmonized_bam.complete.json"
    if output == source:
        raise ValueError("output BAM resolves to the source BAM; in-place filtering is forbidden")
    if output.parent == source.parent:
        raise ValueError(
            "output sandbox resolves to the source run's process directory; "
            "choose a separate sandbox"
        )
    if output.exists() or manifest_path.exists():
        raise FileExistsError(
            f"sandbox already contains a harmonized BAM/manifest: {process_dir}; "
            "choose a new output directory"
        )
    process_dir.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp.{os.getpid()}")

    selected_backend = backend
    if selected_backend == "auto":
        try:
            import pysam  # noqa: F401
            selected_backend = "pysam"
        except ImportError:
            selected_backend = "samtools"

    source_sha256 = _sha256(source)
    try:
        if selected_backend == "pysam":
            try:
                stats, header_text = _filter_with_pysam(source, temporary, evidence.ids)
            except ImportError as exc:
                raise RuntimeError("pysam backend requested but pysam is unavailable") from exc
        elif selected_backend == "samtools":
            executable = samtools or shutil.which("samtools")
            if not executable:
                raise RuntimeError(
                    "neither pysam nor samtools is available; activate the splitpipe "
                    "environment or provide --samtools"
                )
            stats, header_text = _filter_with_samtools(
                source, temporary, evidence.ids, str(executable)
            )
        else:
            raise ValueError(f"unknown BAM backend: {selected_backend}")

        if stats.records_before != stats.records_removed + stats.records_after:
            raise RuntimeError("BAM record-count invariant failed before installation")
        os.replace(temporary, output)
    except Exception:
        temporary.unlink(missing_ok=True)
        # Do not leave an apparently usable empty sandbox contract.
        raise

    output_sha256 = _sha256(output)
    missing_ids = len(evidence.ids) - stats.positive_qnames_seen
    payload = {
        "schema_version": 1,
        "pipeline_version": PIPELINE_VERSION,
        "status": "complete",
        "complete": True,
        "alignment_reused": True,
        "library": library or evidence.library,
        "source_bam": str(source),
        "source_bam_sha256": source_sha256,
        "source_bam_kind": (
            "annotated" if source.name.endswith("Aligned_anno.bam") else "name_sorted"
        ),
        "output_bam": str(output),
        "output_bam_sha256": output_sha256,
        "backend": selected_backend,
        "name_sort_validated": True,
        "header_sha256": hashlib.sha256(header_text.encode()).hexdigest(),
        "positive_evidence_type": evidence.source_type,
        "positive_evidence_path": evidence.source_path,
        "filter_run_id": evidence.filter_run_id,
        "filter_reference_sha256": evidence.reference_sha256,
        "filter_settings_sha256": evidence.filter_settings_sha256,
        # Historical alignments can remove positive QNAMEs but cannot apply a
        # new sequence-trimming boundary without rerunning STAR.
        "survivor_tso_trim_policy": "not_applied_existing_alignment",
        "positive_id_records": evidence.record_count,
        "positive_unique_ids": len(evidence.ids),
        "positive_qnames_seen": stats.positive_qnames_seen,
        "positive_ids_missing_from_bam": missing_ids,
        "records_before": stats.records_before,
        "records_removed": stats.records_removed,
        "records_after": stats.records_after,
        "qnames_before": stats.qnames_before,
        "qnames_removed": stats.qnames_removed,
        "qnames_after": stats.qnames_after,
        "qnames_with_original_header": stats.qnames_with_original_header,
        "completed_at": datetime.now(timezone.utc).isoformat(),
    }
    atomic_json(manifest_path, payload)
    return payload, manifest_path


def main():
    parser = argparse.ArgumentParser(
        description="Filter Step-1a-positive QNAMEs from an existing split-pipe BAM"
    )
    parser.add_argument("--source-bam", required=True)
    parser.add_argument(
        "--output-dir", required=True,
        help="New sandbox root; output is process/barcode_headAligned_sorted.bam",
    )
    parser.add_argument(
        "--register-fresh", action="store_true",
        help="Register an already filtered fresh-run BAM; do not filter records",
    )
    evidence = parser.add_mutually_exclusive_group(required=False)
    evidence.add_argument("--filter-manifest")
    evidence.add_argument("--positive-ids")
    evidence.add_argument("--decision-rows")
    parser.add_argument("--backend", choices=("auto", "pysam", "samtools"), default="auto")
    parser.add_argument("--samtools", help="Explicit samtools executable for fallback mode")
    parser.add_argument(
        "--library", help="Sublibrary label (inferred from --filter-manifest when available)"
    )
    args = parser.parse_args()

    try:
        if args.register_fresh:
            if not args.filter_manifest:
                raise ValueError("--register-fresh requires --filter-manifest")
            payload, manifest_path = register_fresh_harmonized_bam(
                args.source_bam, args.output_dir, args.filter_manifest,
                library=args.library,
            )
            print(f"Registered fresh harmonized BAM: {manifest_path}")
            return
        positive = load_positive_evidence(
            filter_manifest=args.filter_manifest,
            positive_ids=args.positive_ids,
            decision_rows=args.decision_rows,
        )
        if args.library and positive.library and args.library != positive.library:
            raise ValueError(
                f"--library {args.library!r} does not match filter manifest library "
                f"{positive.library!r}"
            )
        payload, manifest_path = filter_existing_bam(
            args.source_bam, args.output_dir, positive,
            backend=args.backend, samtools=args.samtools, library=args.library,
        )
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")

    print(
        f"Filtered BAM: {payload['records_before']:,} records before, "
        f"{payload['records_removed']:,} removed, {payload['records_after']:,} after"
    )
    print(f"Completion manifest: {manifest_path}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Environment-gated split-pipe preprocessing hook for GEX decontamination.

The in-memory split-pipe preprocessing overlay loads this file only when
``CSM_GEX_DECONTAM_HOOK`` is set by Step 1a. A packed guide index is constructed
in the parent before split-pipe forks preprocessing workers. Each worker writes
private rejected-pair/statistic parts; the parent combines them atomically
before split-pipe advances to human-genome STAR alignment.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
import sys
import time
from collections import Counter
from dataclasses import asdict
from importlib import metadata
from pathlib import Path
from typing import Optional

import yaml


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from gex_grna_filter import (  # noqa: E402
    GexGrnaFilter,
    filter_settings_sha256,
    load_filter_from_config,
)


_CLASSIFIER: Optional[GexGrnaFilter] = None
_CONFIG: Optional[dict] = None
_CONFIG_PATH: Optional[Path] = None
_OUTPUT_DIR: Optional[Path] = None
_LIBRARY = ""
_RUN_ID = ""
_WRITE_REJECTED = False
_WRITE_DECISIONS = False
_KEEP_PARTS = False
_INPUT_MANIFEST: Optional[Path] = None
_WORKER = None

SURVIVOR_TRIM_POLICY = "full_tso_endpoint_then_parse_exact_17mer"
SURVIVOR_TRIM_STAT_KEYS = {
    "full_tso": "csm_gex_survivor_full_tso_trimmed_reads",
    "parse_exact_17mer": "csm_gex_survivor_parse_17mer_trimmed_reads",
    "none": "csm_gex_survivor_untrimmed_reads",
}

_STAT_KEYS = (
    "csm_gex_filter_reads_examined",
    "csm_gex_filter_valid_barcode_reads_examined",
    "csm_gex_filter_invalid_barcode_reads_examined",
    "csm_gex_tso_candidate_reads",
    "csm_gex_tso_exact_reads",
    "csm_gex_tso_fuzzy_reads",
    "csm_gex_filter_positive_reads",
    "csm_gex_filter_negative_reads",
    "csm_gex_filter_positive_valid_barcode_reads",
    "csm_gex_filter_positive_invalid_barcode_reads",
    "csm_gex_filter_positive_ambiguous_reads",
    "csm_gex_filter_positive_unique_sequence_reads",
    "csm_gex_filter_positive_exact_reads",
    "csm_gex_filter_positive_hamming_1_reads",
    "csm_gex_step1c_rejected_raw_pairs",
    "csm_gex_step1c_barcode_evidence_rows",
    "csm_gex_survivor_full_tso_trimmed_reads",
    "csm_gex_survivor_parse_17mer_trimmed_reads",
    "csm_gex_survivor_untrimmed_reads",
)
_STAT_INDEX = {key: index for index, key in enumerate(_STAT_KEYS)}
_SURVIVOR_STAT_NAMES = frozenset(SURVIVOR_TRIM_STAT_KEYS.values())


def _as_text(value) -> str:
    """Convert a raw split-pipe FASTQ field only when an artifact needs it."""

    if isinstance(value, bytes):
        return value.decode().strip()
    return str(value).strip()


def _atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def _enabled_rejected_output(decontamination: dict) -> bool:
    raw = (decontamination.get("io", {}) or {}).get(
        "write_rejected_raw_pairs", True
    )
    if isinstance(raw, str) and raw.lower() == "auto":
        return bool(
            (decontamination.get("step1c", {}) or {}).get("enabled", False)
            or (decontamination.get("existing_alignment", {}) or {}).get(
                "enabled", False
            )
        )
    return bool(raw)


def preload_from_environment() -> GexGrnaFilter:
    """Validate config and construct the immutable classifier once."""
    global _CLASSIFIER, _CONFIG, _CONFIG_PATH, _OUTPUT_DIR
    global _LIBRARY, _RUN_ID, _WRITE_REJECTED, _WRITE_DECISIONS, _KEEP_PARTS
    global _INPUT_MANIFEST
    if _CLASSIFIER is not None:
        return _CLASSIFIER

    config_raw = os.environ.get("CSM_GEX_DECONTAM_CONFIG", "")
    output_raw = os.environ.get("CSM_GEX_DECONTAM_OUTPUT_DIR", "")
    _LIBRARY = os.environ.get("CSM_GEX_DECONTAM_LIBRARY", "")
    if not config_raw or not output_raw or not _LIBRARY:
        raise RuntimeError(
            "CSM GEX hook requires CONFIG, OUTPUT_DIR, and LIBRARY environment values"
        )

    _CONFIG_PATH = Path(config_raw).expanduser().resolve()
    _OUTPUT_DIR = Path(output_raw).expanduser().resolve()
    with _CONFIG_PATH.open() as handle:
        _CONFIG = yaml.safe_load(handle) or {}
    decontamination = _CONFIG.get("gex_decontamination", {}) or {}
    if not bool(decontamination.get("enabled", False)):
        raise RuntimeError("CSM hook was requested but gex_decontamination.enabled is false")
    integration = decontamination.get("integration", "splitpipe_pre_router")
    if integration not in {"splitpipe_pre_router", "splitpipe_pre_candidate_router"}:
        raise RuntimeError(
            "v1.4 supports integration: splitpipe_pre_router only "
            "(splitpipe_pre_candidate_router is accepted as a legacy alias)"
        )
    if integration == "splitpipe_pre_candidate_router":
        print(
            "WARNING: integration splitpipe_pre_candidate_router is deprecated; "
            "use splitpipe_pre_router",
            file=sys.stderr,
        )
    existing_config = decontamination.get("existing_alignment", {}) or {}
    if existing_config.get("policy", "auto_sandbox") != "auto_sandbox":
        raise RuntimeError(
            "v1.4 supports existing_alignment.policy: auto_sandbox only"
        )

    policy_path = decontamination.get("qc_policy_file")
    if bool(decontamination.get("require_qc_policy", False)):
        if not policy_path:
            raise RuntimeError("require_qc_policy is true but qc_policy_file is empty")
        candidate = Path(policy_path).expanduser()
        if not candidate.is_absolute():
            candidate = _CONFIG_PATH.parent / candidate
        if not candidate.is_file():
            raise RuntimeError(f"required QC filter policy does not exist: {candidate}")
        with candidate.open() as handle:
            policy = yaml.safe_load(handle) or {}
        if not bool(policy.get("approved", False)):
            raise RuntimeError(f"QC filter policy is not approved: {candidate}")

    _OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    _RUN_ID = os.environ.get("CSM_GEX_DECONTAM_RUN_ID", "") or (
        f"{int(time.time())}.{os.getpid()}"
    )
    input_manifest_raw = os.environ.get("CSM_GEX_DECONTAM_INPUT_MANIFEST", "")
    _INPUT_MANIFEST = (
        Path(input_manifest_raw).expanduser().resolve()
        if input_manifest_raw
        else None
    )
    if _INPUT_MANIFEST is not None and not _INPUT_MANIFEST.is_file():
        raise RuntimeError(
            f"CSM input identity manifest does not exist: {_INPUT_MANIFEST}"
        )
    _WRITE_REJECTED = _enabled_rejected_output(decontamination)
    step1c_enabled = bool((decontamination.get("step1c", {}) or {}).get("enabled", False))
    combine_enabled = bool(
        (decontamination.get("combine_with_grna_arm", {}) or {}).get("enabled", False)
    )
    if combine_enabled and not step1c_enabled:
        raise RuntimeError(
            "combine_with_grna_arm.enabled requires gex_decontamination.step1c.enabled"
        )
    if combine_enabled and not bool(
        (decontamination.get("combine_with_grna_arm", {}) or {}).get(
            "retain_source_provenance", True
        )
    ):
        raise RuntimeError(
            "v1.4 combined counts require retain_source_provenance: true"
        )
    io_config = decontamination.get("io", {}) or {}
    if not bool(io_config.get("keep_rejected_raw_pairs_after_step1c", True)):
        raise RuntimeError(
            "v1.4 retains raw rejected pairs for audit/resume; "
            "keep_rejected_raw_pairs_after_step1c must be true"
        )
    stats_config = decontamination.get("filter_stats", {}) or {}
    _WRITE_DECISIONS = bool(
        io_config.get("write_read_decisions", stats_config.get("write_read_decisions", False))
    )
    _KEEP_PARTS = bool(io_config.get("keep_worker_partials", False))
    _CLASSIFIER = load_filter_from_config(_CONFIG_PATH)
    return _CLASSIFIER


def existing_filter_manifest_is_compatible() -> bool:
    """Return whether a reusable preprocessed FASTQ has matching filter state."""
    preload_from_environment()
    if _OUTPUT_DIR is None or _CLASSIFIER is None:
        return False
    path = _OUTPUT_DIR / "rejected_pairs.complete.json"
    if not path.is_file():
        return False
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return False
    if payload.get("complete") is not True or payload.get("library") != _LIBRARY:
        return False
    if payload.get("schema_version") != 4:
        return False
    if payload.get("survivor_tso_trim_policy") != SURVIVOR_TRIM_POLICY:
        return False
    if payload.get("reference_sha256") != _CLASSIFIER.index.stats.reference_sha256:
        return False
    if payload.get("filter_settings_sha256") != filter_settings_sha256(
        _CLASSIFIER.settings
    ):
        return False
    if payload.get("filter_backend") != "packed_inline":
        return False
    if payload.get("tso_match_backend") != _CLASSIFIER.settings.tso_backend:
        return False
    if _WRITE_REJECTED and not payload.get("write_rejected_raw_pairs"):
        return False
    for key in ("barcode_evidence",):
        raw = payload.get(key)
        if not raw or not Path(raw).is_file():
            return False
    if _WRITE_REJECTED:
        for key in ("rejected_r1", "rejected_r2"):
            raw = payload.get(key)
            if not raw or not Path(raw).is_file():
                return False
    if _WRITE_DECISIONS:
        raw = payload.get("filter_decisions")
        if not raw or not Path(raw).is_file():
            return False
    return True


class _WorkerSink:
    def __init__(
        self,
        fastq_index: int,
        chunk,
        parse_tso_sequence: str = "",
        parse_search_end: int = 0,
    ) -> None:
        if _OUTPUT_DIR is None:
            raise RuntimeError("hook was not preloaded")
        chunk_label = "all" if chunk is None else str(chunk)
        self.part_key = f"{_RUN_ID}.f{fastq_index}.c{chunk_label}.p{os.getpid()}"
        self.prefix = _OUTPUT_DIR / f".{self.part_key}"
        self.counts = [0] * len(_STAT_KEYS)
        self.parse_tso_sequence = str(parse_tso_sequence)
        self.parse_search_end = int(parse_search_end)
        self.fastq_index = int(fastq_index)
        self.chunk_index = None if chunk is None else int(chunk)
        self.r1_handle = None
        self.r2_handle = None
        self.decision_handle = None
        self.barcode_handle = None

    def _open_rejected(self) -> None:
        if self.r1_handle is not None:
            return
        self.r1_handle = gzip.open(f"{self.prefix}.R1.fastq.gz", "wt")
        self.r2_handle = gzip.open(f"{self.prefix}.R2.fastq.gz", "wt")

    def _open_decisions(self) -> None:
        if self.decision_handle is None:
            self.decision_handle = gzip.open(
                f"{self.prefix}.decisions.jsonl.gz", "wt"
            )

    def _open_barcode_evidence(self) -> None:
        if self.barcode_handle is None:
            self.barcode_handle = gzip.open(
                f"{self.prefix}.barcode_evidence.tsv.gz", "wt"
            )

    def _increment(self, key: str) -> None:
        self.counts[_STAT_INDEX[key]] += 1

    def _record_decision(self, decision, barcode_valid: bool) -> None:
        counts = self.counts
        counts[_STAT_INDEX["csm_gex_filter_reads_examined"]] += 1
        if barcode_valid:
            counts[_STAT_INDEX["csm_gex_filter_valid_barcode_reads_examined"]] += 1
        else:
            counts[_STAT_INDEX["csm_gex_filter_invalid_barcode_reads_examined"]] += 1
        if decision.tso_found:
            counts[_STAT_INDEX["csm_gex_tso_candidate_reads"]] += 1
            if decision.tso_match_type == "exact":
                counts[_STAT_INDEX["csm_gex_tso_exact_reads"]] += 1
            elif decision.tso_match_type == "fuzzy":
                counts[_STAT_INDEX["csm_gex_tso_fuzzy_reads"]] += 1
            else:
                raise RuntimeError(
                    f"unknown TSO match type {decision.tso_match_type!r}"
                )
        if decision.filter_positive:
            counts[_STAT_INDEX["csm_gex_filter_positive_reads"]] += 1
            if barcode_valid:
                counts[
                    _STAT_INDEX["csm_gex_filter_positive_valid_barcode_reads"]
                ] += 1
            else:
                counts[
                    _STAT_INDEX["csm_gex_filter_positive_invalid_barcode_reads"]
                ] += 1
            if decision.sequence_ambiguous:
                counts[
                    _STAT_INDEX["csm_gex_filter_positive_ambiguous_reads"]
                ] += 1
            else:
                counts[
                    _STAT_INDEX["csm_gex_filter_positive_unique_sequence_reads"]
                ] += 1
            if decision.confidence_tier == "exact":
                counts[_STAT_INDEX["csm_gex_filter_positive_exact_reads"]] += 1
            elif decision.confidence_tier == "hamming_1":
                counts[
                    _STAT_INDEX["csm_gex_filter_positive_hamming_1_reads"]
                ] += 1
            else:
                raise RuntimeError(
                    f"unknown positive confidence tier {decision.confidence_tier!r}"
                )
        else:
            counts[_STAT_INDEX["csm_gex_filter_negative_reads"]] += 1

    def _record_rejected_pair(
        self,
        *,
        header1,
        seq1,
        quality1,
        header2,
        seq2,
        quality2,
    ) -> None:
        self._open_rejected()
        self.r1_handle.write(
            f"{_as_text(header1)}\n{_as_text(seq1)}\n+\n{_as_text(quality1)}\n"
        )
        self.r2_handle.write(
            f"{_as_text(header2)}\n{_as_text(seq2)}\n+\n{_as_text(quality2)}\n"
        )
        self.counts[_STAT_INDEX["csm_gex_step1c_rejected_raw_pairs"]] += 1

    def _record_barcode_evidence(
        self,
        *,
        header1,
        bc_wells: str,
        bc_bci: str,
        umi: str,
        seq1,
        quality1,
    ) -> None:
        """Store split-pipe-corrected barcode evidence for a valid positive."""

        read_id = _as_text(header1).split()[0].lstrip("@")
        fields = (
            read_id,
            str(bc_wells),
            str(bc_bci),
            str(umi),
            _as_text(seq1),
            _as_text(quality1),
        )
        if any("\t" in value or "\n" in value for value in fields):
            raise RuntimeError("barcode-evidence fields cannot contain tabs/newlines")
        self._open_barcode_evidence()
        self.barcode_handle.write("\t".join(fields) + "\n")
        self.counts[_STAT_INDEX["csm_gex_step1c_barcode_evidence_rows"]] += 1

    def process_compact(self, *, seq1: str, barcode_valid: bool):
        """Fuse default classification, statistics, and survivor trimming."""

        if _CLASSIFIER is None:
            raise RuntimeError("hook classifier was not preloaded")
        if _WRITE_DECISIONS:
            raise RuntimeError("compact hook route cannot write full decisions")
        decision = _CLASSIFIER.route(seq1)
        self._record_decision(decision, barcode_valid)
        if decision.filter_positive or not barcode_valid:
            return decision, 0
        if decision.tso_found:
            trim_end = int(decision.tso_end)
            if trim_end <= 0 or trim_end > len(seq1):
                raise RuntimeError(
                    f"invalid full-TSO trim endpoint {trim_end} for R1 length {len(seq1)}"
                )
            self.counts[
                _STAT_INDEX["csm_gex_survivor_full_tso_trimmed_reads"]
            ] += 1
            return decision, trim_end
        index = seq1.find(
            self.parse_tso_sequence, 0, self.parse_search_end
        )
        if index >= 0:
            self.counts[
                _STAT_INDEX["csm_gex_survivor_parse_17mer_trimmed_reads"]
            ] += 1
            return decision, index + len(self.parse_tso_sequence)
        self.counts[_STAT_INDEX["csm_gex_survivor_untrimmed_reads"]] += 1
        return decision, 0

    def classify_and_record_decision(
        self,
        *,
        header1: str,
        seq1: str,
        quality1: str,
        header2: str,
        seq2: str,
        quality2: str,
        barcode_valid: bool,
    ):
        if _CLASSIFIER is None:
            raise RuntimeError("hook classifier was not preloaded")
        read_id = _as_text(header1).split()[0].lstrip("@")
        decision = (
            _CLASSIFIER.classify(read_id, seq1, self.part_key)
            if _WRITE_DECISIONS
            else _CLASSIFIER.route(seq1)
        )
        self._record_decision(decision, barcode_valid)
        if decision.filter_positive:
            if _WRITE_REJECTED:
                self._record_rejected_pair(
                    header1=header1,
                    seq1=seq1,
                    quality1=quality1,
                    header2=header2,
                    seq2=seq2,
                    quality2=quality2,
                )

        if _WRITE_DECISIONS:
            self._open_decisions()
            payload = decision.to_dict()
            payload["barcode_valid"] = bool(barcode_valid)
            self.decision_handle.write(json.dumps(payload, sort_keys=True) + "\n")
        return decision

    def record_survivor_trim(self, route: str) -> None:
        try:
            stat_key = SURVIVOR_TRIM_STAT_KEYS[route]
        except KeyError as exc:
            raise ValueError(f"unknown CSM survivor trim route: {route!r}") from exc
        self.counts[_STAT_INDEX[stat_key]] += 1

    def close(self) -> dict:
        for handle in (
            self.r1_handle,
            self.r2_handle,
            self.decision_handle,
            self.barcode_handle,
        ):
            if handle is not None:
                handle.close()
        stats = {
            key: self.counts[index]
            for index, key in enumerate(_STAT_KEYS)
            if self.counts[index] or key in _SURVIVOR_STAT_NAMES
        }
        payload = {
            "run_id": _RUN_ID,
            "part_key": self.part_key,
            "fastq_index": self.fastq_index,
            "chunk_index": self.chunk_index,
            "statistics": stats,
        }
        _atomic_json(Path(f"{self.prefix}.stats.json"), payload)
        return stats


def start_worker(
    *,
    fastq_index: int,
    chunk=None,
    parse_tso_sequence: str = "",
    parse_search_end: int = 0,
) -> None:
    global _WORKER
    preload_from_environment()
    _WORKER = _WorkerSink(
        fastq_index,
        chunk,
        parse_tso_sequence=parse_tso_sequence,
        parse_search_end=parse_search_end,
    )


def write_read_decisions_enabled() -> bool:
    preload_from_environment()
    return _WRITE_DECISIONS


def write_rejected_raw_pairs_enabled() -> bool:
    preload_from_environment()
    return _WRITE_REJECTED


def process_compact(*, seq1: str, barcode_valid: bool):
    if _WORKER is None:
        raise RuntimeError("CSM GEX hook worker was not started")
    return _WORKER.process_compact(seq1=seq1, barcode_valid=barcode_valid)


def record_rejected_pair(**kwargs) -> None:
    if _WORKER is None:
        raise RuntimeError("CSM GEX hook worker was not started")
    if not _WRITE_REJECTED:
        raise RuntimeError("rejected raw-pair output is disabled")
    _WORKER._record_rejected_pair(**kwargs)


def record_barcode_evidence(**kwargs) -> None:
    if _WORKER is None:
        raise RuntimeError("CSM GEX hook worker was not started")
    _WORKER._record_barcode_evidence(**kwargs)


def classify_and_record(**kwargs) -> bool:
    if _WORKER is None:
        raise RuntimeError("CSM GEX hook worker was not started")
    return bool(_WORKER.classify_and_record_decision(**kwargs).filter_positive)


def classify_and_record_decision(**kwargs):
    """Classify/record one raw pair and return its full routing decision."""
    if _WORKER is None:
        raise RuntimeError("CSM GEX hook worker was not started")
    return _WORKER.classify_and_record_decision(**kwargs)


def resolve_survivor_trim(
    *, decision, seq1: str, parse_tso_sequence: str, parse_search_end: int
) -> tuple[int, str]:
    """Choose trimming for a filter-negative, valid-barcode STAR survivor.

    A full TSO recognized by the v1.4 detector owns the trim boundary, including
    fuzzy matches. Otherwise this reproduces Parse's exact 17-mer search.
    """
    if decision is None:
        raise RuntimeError("CSM survivor trimming requires a filter decision")
    if decision.filter_positive:
        raise RuntimeError("filter-positive reads cannot enter survivor trimming")
    if decision.tso_found:
        trim_end = int(decision.tso_end)
        if trim_end <= 0 or trim_end > len(seq1):
            raise RuntimeError(
                f"invalid full-TSO trim endpoint {trim_end} for R1 length {len(seq1)}"
            )
        return trim_end, "full_tso"

    index = seq1.find(parse_tso_sequence, 0, parse_search_end)
    if index >= 0:
        return index + len(parse_tso_sequence), "parse_exact_17mer"
    return 0, "none"


def record_survivor_trim(route: str) -> None:
    if _WORKER is None:
        raise RuntimeError("CSM GEX hook worker was not started")
    _WORKER.record_survivor_trim(route)


def close_worker() -> dict:
    global _WORKER
    if _WORKER is None:
        return {}
    stats = _WORKER.close()
    _WORKER = None
    return stats


def _combine_binary(parts: list[Path], destination: Path) -> None:
    temporary = destination.with_name(f".{destination.name}.tmp.{os.getpid()}")
    with temporary.open("wb") as output:
        for part in parts:
            with part.open("rb") as source:
                shutil.copyfileobj(source, output, length=1024 * 1024)
        output.flush()
        os.fsync(output.fileno())
    os.replace(temporary, destination)


def _fastq_digest_and_count(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    count = 0
    with gzip.open(path, "rt") as handle:
        while True:
            header = handle.readline().rstrip("\n")
            if not header:
                break
            sequence = handle.readline().rstrip("\n")
            plus = handle.readline().rstrip("\n")
            quality = handle.readline().rstrip("\n")
            if not sequence or not plus.startswith("+") or len(sequence) != len(quality):
                raise RuntimeError(f"malformed rejected FASTQ record in {path}")
            digest.update(header.split()[0].lstrip("@").encode("utf-8"))
            digest.update(b"\n")
            count += 1
    return count, digest.hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _barcode_evidence_count(path: Path) -> int:
    count = 0
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n")
        expected = "read_id\tbc_wells\tbc_bci\tumi\tr1_seq\tr1_qual"
        if header != expected:
            raise RuntimeError(f"invalid barcode-evidence header in {path}")
        for line in handle:
            if len(line.rstrip("\n").split("\t")) != 6:
                raise RuntimeError(f"malformed barcode-evidence row in {path}")
            count += 1
    return count


def _part_order(payload: dict) -> tuple[int, int, str]:
    chunk = payload.get("chunk_index")
    return (
        int(payload["fastq_index"]),
        -1 if chunk is None else int(chunk),
        str(payload.get("part_key", "")),
    )


def finalize_from_environment() -> str:
    """Combine worker parts and atomically install the completion manifest."""
    preload_from_environment()
    if _OUTPUT_DIR is None or _CLASSIFIER is None or _CONFIG_PATH is None:
        raise RuntimeError("hook was not preloaded")

    stats_parts = sorted(_OUTPUT_DIR.glob(f".{_RUN_ID}.*.stats.json"))
    if not stats_parts:
        raise RuntimeError(f"no CSM GEX worker statistics found for run {_RUN_ID}")
    totals = Counter()
    part_keys = []
    ordered_parts = []
    for path in stats_parts:
        with path.open() as handle:
            payload = json.load(handle)
        if payload.get("run_id") != _RUN_ID:
            raise RuntimeError(f"worker part run ID mismatch: {path}")
        ordered_parts.append((payload, path))
    ordered_parts.sort(key=lambda item: _part_order(item[0]))
    for payload, path in ordered_parts:
        totals.update(payload.get("statistics", {}))
        part_keys.append(payload.get("part_key", path.name))

    def ordered_artifact_parts(suffix: str) -> list[Path]:
        parts = []
        for payload, _ in ordered_parts:
            candidate = _OUTPUT_DIR / f".{payload['part_key']}.{suffix}"
            if candidate.is_file():
                parts.append(candidate)
        return parts

    rejected_r1 = rejected_r2 = None
    rejected_count = 0
    ordered_id_sha256 = None
    if _WRITE_REJECTED:
        r1_parts = ordered_artifact_parts("R1.fastq.gz")
        r2_parts = ordered_artifact_parts("R2.fastq.gz")
        if len(r1_parts) != len(r2_parts):
            raise RuntimeError("rejected R1/R2 worker part counts differ")
        rejected_r1 = _OUTPUT_DIR / "rejected_R1.fastq.gz"
        rejected_r2 = _OUTPUT_DIR / "rejected_R2.fastq.gz"
        if r1_parts:
            _combine_binary(r1_parts, rejected_r1)
            _combine_binary(r2_parts, rejected_r2)
        else:
            with gzip.open(rejected_r1, "wt"):
                pass
            with gzip.open(rejected_r2, "wt"):
                pass
        r1_count, ordered_id_sha256 = _fastq_digest_and_count(rejected_r1)
        r2_count, r2_digest = _fastq_digest_and_count(rejected_r2)
        if r1_count != r2_count or ordered_id_sha256 != r2_digest:
            raise RuntimeError("final rejected R1/R2 names or counts differ")
        rejected_count = r1_count
        expected = totals.get("csm_gex_step1c_rejected_raw_pairs", 0)
        if rejected_count != expected:
            raise RuntimeError(
                f"rejected FASTQ count {rejected_count} != worker statistic {expected}"
            )

    barcode_evidence = _OUTPUT_DIR / "valid_positive_barcode_evidence.tsv.gz"
    barcode_parts = ordered_artifact_parts("barcode_evidence.tsv.gz")
    header_part = _OUTPUT_DIR / f".{_RUN_ID}.barcode_evidence.header.tsv.gz"
    with gzip.open(header_part, "wt") as handle:
        handle.write("read_id\tbc_wells\tbc_bci\tumi\tr1_seq\tr1_qual\n")
    _combine_binary([header_part, *barcode_parts], barcode_evidence)
    header_part.unlink()
    barcode_evidence_rows = _barcode_evidence_count(barcode_evidence)
    expected_barcode_rows = totals.get(
        "csm_gex_step1c_barcode_evidence_rows", 0
    )
    if barcode_evidence_rows != expected_barcode_rows:
        raise RuntimeError(
            "barcode-evidence row count "
            f"{barcode_evidence_rows} != worker statistic {expected_barcode_rows}"
        )
    expected_valid_positives = totals.get(
        "csm_gex_filter_positive_valid_barcode_reads", 0
    )
    if barcode_evidence_rows != expected_valid_positives:
        raise RuntimeError(
            "barcode-evidence row count "
            f"{barcode_evidence_rows} != valid positive count {expected_valid_positives}"
        )

    decisions_path = None
    if _WRITE_DECISIONS:
        decision_parts = ordered_artifact_parts("decisions.jsonl.gz")
        decisions_path = _OUTPUT_DIR / "filter_decisions.jsonl.gz"
        if decision_parts:
            _combine_binary(decision_parts, decisions_path)
        else:
            with gzip.open(decisions_path, "wt"):
                pass

    manifest = {
        "schema_version": 4,
        "pipeline_version": "1.4",
        "status": "complete",
        "complete": True,
        "library": _LIBRARY,
        "run_id": _RUN_ID,
        "config": str(_CONFIG_PATH),
        "filter_backend": "packed_inline",
        "tso_match_backend": _CLASSIFIER.settings.tso_backend,
        "tso_match_backend_version": (
            metadata.version("edlib")
            if _CLASSIFIER.settings.tso_backend == "edlib"
            else None
        ),
        "survivor_tso_trim_policy": SURVIVOR_TRIM_POLICY,
        "reference": _CLASSIFIER.settings.reference_path,
        "reference_sha256": _CLASSIFIER.index.stats.reference_sha256,
        "filter_settings_sha256": filter_settings_sha256(_CLASSIFIER.settings),
        "index_statistics": asdict(_CLASSIFIER.index.stats),
        "statistics": dict(totals),
        "write_rejected_raw_pairs": _WRITE_REJECTED,
        "rejected_pair_count": rejected_count,
        "ordered_read_id_sha256": ordered_id_sha256,
        "rejected_r1": str(rejected_r1) if rejected_r1 else None,
        "rejected_r2": str(rejected_r2) if rejected_r2 else None,
        "rejected_r1_sha256": _sha256_file(rejected_r1) if rejected_r1 else None,
        "rejected_r2_sha256": _sha256_file(rejected_r2) if rejected_r2 else None,
        "rejected_r1_size_bytes": rejected_r1.stat().st_size if rejected_r1 else 0,
        "rejected_r2_size_bytes": rejected_r2.stat().st_size if rejected_r2 else 0,
        "barcode_evidence": str(barcode_evidence),
        "barcode_evidence_rows": barcode_evidence_rows,
        "barcode_evidence_sha256": _sha256_file(barcode_evidence),
        "barcode_evidence_size_bytes": barcode_evidence.stat().st_size,
        "filter_decisions": str(decisions_path) if decisions_path else None,
        "worker_parts": part_keys,
    }
    if _INPUT_MANIFEST is not None:
        with _INPUT_MANIFEST.open() as handle:
            input_identity = json.load(handle)
        if input_identity.get("library") != _LIBRARY:
            raise RuntimeError("input identity manifest library mismatch")
        manifest["input_identity"] = input_identity
        manifest["input_identity_manifest"] = str(_INPUT_MANIFEST)
        manifest["input_identity_manifest_sha256"] = _sha256_file(_INPUT_MANIFEST)
    examined = totals.get("csm_gex_filter_reads_examined", 0)
    tso_candidates = totals.get("csm_gex_tso_candidate_reads", 0)
    valid_examined = totals.get("csm_gex_filter_valid_barcode_reads_examined", 0)
    positives = totals.get("csm_gex_filter_positive_reads", 0)
    valid_positives = totals.get(
        "csm_gex_filter_positive_valid_barcode_reads", 0
    )
    star_survivors = max(0, valid_examined - valid_positives)
    full_tso_survivor_trims = totals.get(
        "csm_gex_survivor_full_tso_trimmed_reads", 0
    )
    parse_17mer_survivor_trims = totals.get(
        "csm_gex_survivor_parse_17mer_trimmed_reads", 0
    )
    untrimmed_survivors = totals.get("csm_gex_survivor_untrimmed_reads", 0)
    routed_survivors = (
        full_tso_survivor_trims
        + parse_17mer_survivor_trims
        + untrimmed_survivors
    )
    if routed_survivors != star_survivors:
        raise RuntimeError(
            "CSM STAR-survivor routing count mismatch: "
            f"expected {star_survivors}, observed {routed_survivors}"
        )
    manifest["rates"] = {
        "filter_positive_fraction_of_examined_reads": (
            positives / examined if examined else 0.0
        ),
        "filter_positive_fraction_of_tso_candidates": (
            positives / tso_candidates if tso_candidates else 0.0
        ),
        "removed_from_star_fraction_of_valid_barcode_reads": (
            valid_positives / valid_examined if valid_examined else 0.0
        ),
        "full_tso_trim_fraction_of_star_survivors": (
            full_tso_survivor_trims / star_survivors if star_survivors else 0.0
        ),
        "parse_17mer_trim_fraction_of_star_survivors": (
            parse_17mer_survivor_trims / star_survivors if star_survivors else 0.0
        ),
    }
    manifest_path = _OUTPUT_DIR / "rejected_pairs.complete.json"
    _atomic_json(manifest_path, manifest)
    _atomic_json(_OUTPUT_DIR / "filter_summary.json", manifest)

    if not _KEEP_PARTS:
        for pattern in (
            f".{_RUN_ID}.*.stats.json",
            f".{_RUN_ID}.*.R1.fastq.gz",
            f".{_RUN_ID}.*.R2.fastq.gz",
            f".{_RUN_ID}.*.decisions.jsonl.gz",
            f".{_RUN_ID}.*.barcode_evidence.tsv.gz",
        ):
            for part in _OUTPUT_DIR.glob(pattern):
                part.unlink()
    return str(manifest_path)

#!/usr/bin/env python3
"""Bounded multiprocessing engine for production GEX-contamination QC.

The parent process
streams paired gzip FASTQs once, keeps at most ``2 * workers`` batches in
flight, and merges worker-local counters and HLL registers.  Forked workers
share the read-only guide/barcode models copy-on-write, and only one aggregate
per worker crosses the process boundary.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import multiprocessing
import resource
import time
import traceback
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Iterator, Optional, Sequence, Tuple

try:
    from . import qc_gex_grna_contamination as qc
except ImportError:  # Direct execution from scripts/ on the HPC.
    import qc_gex_grna_contamination as qc


_MODEL = None
_BARCODE_CORRECTOR = None
_READ_LENGTH = 54
_SHORT_CONFIG = None
_SELECTED_METHOD = ""
_MOLECULE_ENABLED = False
_HLL_PRECISION = 18
_BASE_COMPOSITION_LENGTH = 100
_ALIGNMENT_ENABLED = False
_ALIGNMENT_FRACTION = 1.0


def _open_fastq_binary(path: str):
    return gzip.open(path, "rb") if str(path).endswith(".gz") else open(path, "rb")


def _read_fastq_record_binary(handle, path: str, record_number: int):
    header = handle.readline()
    if not header:
        return None
    sequence = handle.readline()
    plus = handle.readline()
    quality = handle.readline()
    if not sequence or not plus or not quality:
        raise ValueError(f"Truncated FASTQ record {record_number} in {path}")
    header = header.rstrip(b"\r\n")
    sequence = sequence.rstrip(b"\r\n").upper()
    quality = quality.rstrip(b"\r\n")
    if not header.startswith(b"@"):
        raise ValueError(f"Invalid FASTQ header in {path}: {header[:80]!r}")
    if not plus.startswith(b"+"):
        raise ValueError(f"Invalid FASTQ plus line in {path}: {header[:80]!r}")
    if len(sequence) != len(quality):
        raise ValueError(
            f"FASTQ sequence/quality mismatch at record {record_number} in {path}"
        )
    read_id = header.split(maxsplit=1)[0].lstrip(b"@")
    if read_id.endswith((b"/1", b"/2")):
        read_id = read_id[:-2]
    return read_id, sequence


def iter_fastq_sequence_pairs(r1_path: str, r2_path: str):
    """Yield validated uppercase sequence pairs without allocating record objects."""
    with _open_fastq_binary(r1_path) as r1_handle, _open_fastq_binary(
        r2_path
    ) as r2_handle:
        record_number = 0
        while True:
            next_number = record_number + 1
            r1 = _read_fastq_record_binary(r1_handle, r1_path, next_number)
            r2 = _read_fastq_record_binary(r2_handle, r2_path, next_number)
            if r1 is None and r2 is None:
                return
            if r1 is None or r2 is None:
                raise ValueError(f"FASTQ length mismatch: {r1_path} versus {r2_path}")
            record_number = next_number
            if r1[0] != r2[0]:
                raise ValueError(
                    f"FASTQ ID mismatch at record {record_number}: "
                    f"{r1[0]!r} != {r2[0]!r}"
                )
            yield r1[1].decode("ascii"), r2[1].decode("ascii")


def _fit_perfect_counts_fast(corrector, r2_paths, sample_start: int, sample_end: int):
    """Fit the exact production threshold with a sequence-only binary parser."""
    if sample_start < 0 or sample_end <= sample_start:
        raise ValueError("barcode calibration interval must satisfy 0 <= start < end")
    counts = Counter()
    index = 0
    for r2_path in qc.normalize_fastq_paths(r2_paths):
        if index >= sample_end:
            break
        with _open_fastq_binary(r2_path) as handle:
            record_number = 0
            while index < sample_end:
                record_number += 1
                record = _read_fastq_record_binary(
                    handle,
                    r2_path,
                    record_number,
                )
                if record is None:
                    break
                if index >= sample_start:
                    sequence = record[1].decode("ascii")
                    umi, bc1, bc2, bc3 = corrector.raw_parts(sequence)
                    if (
                        umi
                        and bc1 in corrector.sequences_by_round[1]
                        and bc2 in corrector.sequences_by_round[2]
                        and bc3 in corrector.sequences_by_round[3]
                    ):
                        counts[(bc1, bc2, bc3)] += 1
                index += 1
    corrector.perfect_combination_counts = counts
    if counts:
        ordered = sorted(counts.values(), reverse=True)
        total = sum(ordered)
        cumulative = 0
        threshold = ordered[-1]
        best_delta = float("inf")
        for value in ordered:
            cumulative += value
            delta = abs(cumulative / total - corrector.perfect_cell_fraction)
            if delta < best_delta:
                best_delta = delta
                threshold = value
        corrector.perfect_count_threshold = max(5, int(threshold))


def make_barcode_corrector_fast(
    cfg: dict,
    config_path: str,
    library_name: str,
    r2_paths,
):
    """Construct the production corrector with an equivalent faster calibration scan."""
    molecule_cfg = cfg["gex_grna_contamination_qc"]["molecule_estimation"]
    template = molecule_cfg.get(
        "barcode_data_template",
        "{output_dir}/gex_{library}/process/barcode_data.csv",
    )
    barcode_path = qc.resolve_path(
        template.format(
            output_dir=cfg["paths"]["output_dir"],
            library=library_name,
        ),
        config_path,
    )
    if not Path(barcode_path).is_file():
        raise FileNotFoundError(
            f"Molecule estimation requires split-pipe barcode data: {barcode_path}"
        )
    barcode_df = qc.pd.read_csv(barcode_path, comment="#")
    corrector = qc.PilotBarcodeCorrector(
        barcode_df,
        str(
            molecule_cfg.get(
                "amplicon_sequence",
                "NNNNNNNNNN33333333ATGAGGGGTCAG22222222TCCAACCACCTC11111111",
            )
        ),
        max_mismatches=int(molecule_cfg.get("barcode_max_mismatches", 2)),
        perfect_cell_fraction=float(
            molecule_cfg.get("perfect_barcode_cell_fraction", 0.92)
        ),
    )
    _fit_perfect_counts_fast(
        corrector,
        r2_paths,
        int(molecule_cfg.get("barcode_sample_start", 3_000_000)),
        int(molecule_cfg.get("barcode_sample_end", 4_000_000)),
    )
    return corrector


def canonical_json_sha256(payload: object) -> str:
    serialized = json.dumps(
        payload,
        allow_nan=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(serialized).hexdigest()


class ScanAccumulator:
    """Commutative QC state that can be merged exactly across workers."""

    def __init__(self, *, molecule_enabled: bool, hll_precision: int):
        self.records_evaluated = 0
        self.truth_counts = Counter()
        self.confusion_by_truth = {
            "structure_truth": defaultdict(qc.Confusion),
            "known_structure_truth": defaultdict(qc.Confusion),
            "inclusive_truth": defaultdict(qc.Confusion),
        }
        self.molecule_estimator = (
            qc.MoleculeEstimator(hll_precision) if molecule_enabled else None
        )
        self.guide_counts = Counter()
        self.target_counts = Counter()
        self.tso_offsets = Counter()
        self.positive_stratum_totals = Counter()
        self.positive_stratum_calls = Counter()
        self.hard_negative_stratum_totals = Counter()
        self.hard_negative_stratum_calls = Counter()
        self.base_composition = qc.initialize_base_composition(
            _BASE_COMPOSITION_LENGTH
        )

    def add(
        self,
        r1_sequence: str,
        full,
        short: dict,
        barcode_call,
    ) -> None:
        self.records_evaluated += 1
        calls = short["calls"]
        splitpipe_tso_index = short["splitpipe_tso_index"]
        if _SELECTED_METHOD not in calls:
            raise KeyError(
                f"selected_method {_SELECTED_METHOD!r} is unavailable: "
                f"{sorted(calls)}"
            )
        if self.molecule_estimator is not None:
            self.molecule_estimator.add(barcode_call, r1_sequence, full)
        self.truth_counts[full.tier] += 1
        self.truth_counts["structure_truth"] += int(full.structure_truth)
        self.truth_counts["known_structure_truth"] += int(
            full.known_structure_truth
        )
        self.truth_counts["complete_reference_truth"] += int(
            full.complete_reference_truth
        )
        self.truth_counts["tso_known_guide_truth"] += int(
            full.tso_known_guide_truth
        )
        self.truth_counts["tso_found"] += int(full.tso_found)
        inclusive_truth = bool(full.structure_truth or full.tso_known_guide_truth)
        self.truth_counts["inclusive_truth"] += int(inclusive_truth)
        hard_negative = full.tso_found and not full.structure_truth
        qc.update_base_composition(self.base_composition, "all", r1_sequence)
        if full.structure_truth:
            qc.update_base_composition(
                self.base_composition, "structure_truth", r1_sequence
            )
            self.tso_offsets[full.tso_start] += 1
            if full.guide_names:
                for guide_name in full.guide_names:
                    self.guide_counts[guide_name] += 1 / len(full.guide_names)
            if full.target_genes:
                for target in full.target_genes:
                    self.target_counts[target] += 1 / len(full.target_genes)
        elif full.tso_found:
            qc.update_base_composition(
                self.base_composition, "tso_positive_gex_like", r1_sequence
            )
        for truth_name, truth in (
            ("structure_truth", full.structure_truth),
            ("known_structure_truth", full.known_structure_truth),
            ("inclusive_truth", inclusive_truth),
        ):
            for method, call in calls.items():
                self.confusion_by_truth[truth_name][method].update(
                    truth,
                    bool(call),
                    hard_negative=hard_negative and not truth,
                )
        selected_call = bool(calls[_SELECTED_METHOD])
        if full.structure_truth:
            stratum = (full.tso_start, full.guide_length, full.variant)
            self.positive_stratum_totals[stratum] += 1
            for method, call in calls.items():
                self.positive_stratum_calls[stratum + (method,)] += int(call)
            self.truth_counts["structure_selected_filter_call"] += int(
                selected_call
            )
            self.truth_counts["structure_predicted_splitpipe_tso_trim"] += int(
                splitpipe_tso_index >= 0
            )
            self.truth_counts["structure_splitpipe_trim_known_guide_d2"] += int(
                calls.get("splitpipe_trim_guide_d2", False)
            )
        elif hard_negative:
            stratum = (full.tso_start, 0, 0)
            self.hard_negative_stratum_totals[stratum] += 1
            for method, call in calls.items():
                self.hard_negative_stratum_calls[stratum + (method,)] += int(
                    call
                )
            self.truth_counts["hard_negative_selected_filter_call"] += int(
                selected_call
            )

    def merge(self, other: "ScanAccumulator") -> None:
        self.records_evaluated += other.records_evaluated
        self.truth_counts.update(other.truth_counts)
        self.guide_counts.update(other.guide_counts)
        self.target_counts.update(other.target_counts)
        self.tso_offsets.update(other.tso_offsets)
        self.positive_stratum_totals.update(other.positive_stratum_totals)
        self.positive_stratum_calls.update(other.positive_stratum_calls)
        self.hard_negative_stratum_totals.update(
            other.hard_negative_stratum_totals
        )
        self.hard_negative_stratum_calls.update(
            other.hard_negative_stratum_calls
        )
        for group, positions in other.base_composition.items():
            for index, counts in enumerate(positions):
                self.base_composition[group][index].update(counts)
        for truth_name, source_methods in other.confusion_by_truth.items():
            for method, source in source_methods.items():
                target = self.confusion_by_truth[truth_name][method]
                for field in (
                    "tp",
                    "fp",
                    "tn",
                    "fn",
                    "hard_negative_n",
                    "hard_negative_fp",
                ):
                    setattr(target, field, getattr(target, field) + getattr(source, field))
        if (self.molecule_estimator is None) != (
            other.molecule_estimator is None
        ):
            raise RuntimeError("cannot merge molecule-enabled and disabled scans")
        if self.molecule_estimator is not None:
            _merge_molecule_estimators(
                self.molecule_estimator,
                other.molecule_estimator,
            )


def _merge_molecule_estimators(target, source) -> None:
    if target.precision != source.precision or set(target.sketches) != set(
        source.sketches
    ):
        raise RuntimeError("incompatible molecule estimators")
    for key, source_sketch in source.sketches.items():
        target_registers = target.sketches[key].registers
        for index, value in enumerate(source_sketch.registers):
            if value > target_registers[index]:
                target_registers[index] = value
    target.valid_read_counts.update(source.valid_read_counts)
    target.barcode_edit_counts.update(source.barcode_edit_counts)
    for cell, counts in source.cell_read_counts.items():
        target.cell_read_counts[cell].update(counts)


def _new_accumulator() -> ScanAccumulator:
    if _MODEL is None or _SHORT_CONFIG is None or not _SELECTED_METHOD:
        raise RuntimeError("parallel QC worker was not configured before fork")
    return ScanAccumulator(
        molecule_enabled=_MOLECULE_ENABLED,
        hll_precision=_HLL_PRECISION,
    )


def _accumulate_batch(
    accumulator: ScanAccumulator,
    records: list[Tuple[str, str]],
) -> None:
    for r1_sequence, r2_sequence in records:
        full = _MODEL.classify_full(r1_sequence)
        short = _MODEL.short_features(
            r1_sequence,
            _READ_LENGTH,
            _SHORT_CONFIG,
        )
        barcode = (
            _BARCODE_CORRECTOR.call(r2_sequence)
            if _BARCODE_CORRECTOR is not None
            else None
        )
        accumulator.add(
            r1_sequence,
            full,
            short,
            barcode,
        )


def iter_sequence_batches(
    pairs: Sequence[Tuple[str, str]],
    *,
    max_reads: int,
    batch_size: int,
) -> Iterator[list[Tuple[str, str]]]:
    batch = []
    seen = 0
    for r1_path, r2_path in pairs:
        for r1_sequence, r2_sequence in iter_fastq_sequence_pairs(r1_path, r2_path):
            batch.append((r1_sequence, r2_sequence))
            seen += 1
            if len(batch) >= batch_size:
                yield batch
                batch = []
            if max_reads and seen >= max_reads:
                if batch:
                    yield batch
                return
    if batch:
        yield batch


def _aggregate_worker_loop(task_queue, result_queue) -> None:
    accumulator = None
    error = None
    while True:
        batch = task_queue.get()
        if batch is None:
            break
        if error is not None:
            continue
        try:
            if accumulator is None:
                accumulator = _new_accumulator()
            _accumulate_batch(accumulator, batch)
        except Exception as exc:  # Keep draining so the producer cannot deadlock.
            error = {
                "error_type": type(exc).__name__,
                "error": str(exc),
                "traceback": traceback.format_exc(),
            }
    if error is None:
        if accumulator is None:
            accumulator = _new_accumulator()
        result_queue.put(("ok", accumulator))
    else:
        result_queue.put(("error", error))


def aggregate_batches(
    batches: Iterable[list[Tuple[str, str]]],
    *,
    workers: int,
) -> ScanAccumulator:
    """Classify with bounded queues and one exact aggregate per worker."""
    if workers == 1:
        accumulator = _new_accumulator()
        for batch in batches:
            _accumulate_batch(accumulator, batch)
        return accumulator
    context = multiprocessing.get_context("fork")
    maximum_pending = 2 * workers
    task_queue = context.Queue(maxsize=maximum_pending)
    result_queue = context.Queue()
    processes = [
        context.Process(
            target=_aggregate_worker_loop,
            args=(task_queue, result_queue),
        )
        for _ in range(workers)
    ]
    try:
        for process in processes:
            process.start()
        for batch in batches:
            task_queue.put(batch)
        for _ in processes:
            task_queue.put(None)
        worker_results = [result_queue.get() for _ in processes]
        for process in processes:
            process.join()
        errors = [payload for status, payload in worker_results if status != "ok"]
        abnormal = [process.exitcode for process in processes if process.exitcode != 0]
        if errors or abnormal:
            raise RuntimeError(
                f"parallel QC worker failure: errors={errors}, exitcodes={abnormal}"
            )
        merged = _new_accumulator()
        for _, accumulator in worker_results:
            merged.merge(accumulator)
        return merged
    finally:
        for process in processes:
            if process.is_alive():
                process.terminate()
            process.join()
        task_queue.close()
        result_queue.close()
        task_queue.join_thread()
        result_queue.join_thread()


def configure_runtime(
    *,
    model,
    barcode_corrector,
    read_length: int,
    short_config: dict,
    selected_method: str,
    molecule_enabled: bool,
    hll_precision: int,
    base_composition_length: int,
    alignment_enabled: bool,
    alignment_fraction: float,
) -> None:
    """Configure the fork-inherited immutable production scan state."""

    global _MODEL, _BARCODE_CORRECTOR, _READ_LENGTH, _SHORT_CONFIG
    global _SELECTED_METHOD, _MOLECULE_ENABLED, _HLL_PRECISION
    global _BASE_COMPOSITION_LENGTH, _ALIGNMENT_ENABLED, _ALIGNMENT_FRACTION
    _MODEL = model
    _BARCODE_CORRECTOR = barcode_corrector
    _READ_LENGTH = int(read_length)
    _SHORT_CONFIG = short_config
    _SELECTED_METHOD = str(selected_method)
    _MOLECULE_ENABLED = bool(molecule_enabled)
    _HLL_PRECISION = int(hll_precision)
    _BASE_COMPOSITION_LENGTH = int(base_composition_length)
    _ALIGNMENT_ENABLED = bool(alignment_enabled)
    _ALIGNMENT_FRACTION = float(alignment_fraction)


def _read_evidence_payload(full, short: dict) -> dict:
    guide_hit = short["guide_hit"]
    template_hit = short["template_prefix_hit"]
    inclusive_truth = bool(full.structure_truth or full.tso_known_guide_truth)
    row = {
        "full_length": None,
        "truth_tier": full.tier,
        "structure_truth": full.structure_truth,
        "known_structure_truth": full.known_structure_truth,
        "complete_reference_truth": full.complete_reference_truth,
        "tso_known_guide_truth": full.tso_known_guide_truth,
        "inclusive_truth": inclusive_truth,
        "full_tso_offset": full.tso_start,
        "full_tso_distance": full.tso_distance,
        "full_guide_query": full.guide_query,
        "full_guide_names": ",".join(full.guide_names),
        "full_target_genes": ",".join(full.target_genes),
        "full_variant": full.variant,
        "short_tso_distance": short["tso_distance"],
        "short_tso_at_expected_start": short["tso_at_expected_start"],
        "predicted_splitpipe_tso_index": short["splitpipe_tso_index"],
        "short_guide_distance": guide_hit.distance if guide_hit else -1,
        "short_guide_length": guide_hit.length if guide_hit else 0,
        "short_template_prefix_length": template_hit.length if template_hit else 0,
        "short_template_prefix_distance": (
            template_hit.distance if template_hit else -1
        ),
        "short_template_boundary_bases": (
            template_hit.minimum_boundary_bases if template_hit else 0
        ),
        "short_template_identity_ambiguous": (
            template_hit.identity_ambiguous if template_hit else False
        ),
        "short_template_candidate_guides": (
            ",".join(sorted(template_hit.guide_names)) if template_hit else ""
        ),
        "short_partial_stem_length": short["partial_stem_length"],
        "short_partial_stem_mismatches": (
            short["partial_stem_mismatches"]
            if short["partial_stem_mismatches"] is not None
            else -1
        ),
        "short_template_hamming": (
            short["template_hamming"]
            if short["template_hamming"] is not None
            else -1
        ),
    }
    row.update({f"call_{key}": bool(value) for key, value in short["calls"].items()})
    return row


def _accumulate_artifact_batch(
    accumulator: ScanAccumulator,
    records: list[Tuple[str, str, bool, bool]],
) -> dict:
    read_rows = []
    synthetic = []
    alignment = []
    for local_index, (
        r1_sequence,
        r2_sequence,
        include_read_row,
        include_synthetic,
    ) in enumerate(records):
        full = _MODEL.classify_full(r1_sequence)
        short = _MODEL.short_features(r1_sequence, _READ_LENGTH, _SHORT_CONFIG)
        barcode = (
            _BARCODE_CORRECTOR.call(r2_sequence)
            if _BARCODE_CORRECTOR is not None
            else None
        )
        accumulator.add(r1_sequence, full, short, barcode)
        if include_read_row:
            row = _read_evidence_payload(full, short)
            row["full_length"] = len(r1_sequence)
            read_rows.append((local_index, row))
        if include_synthetic:
            synthetic.append(
                (local_index, bool(short["calls"][_SELECTED_METHOD]))
            )
        if (
            _ALIGNMENT_ENABLED
            and full.structure_truth
            and barcode is not None
            and barcode.valid
            and qc.stable_fraction(
                barcode.cell + "|" + barcode.barcode_indices + ":alignment"
            )
            < _ALIGNMENT_FRACTION
        ):
            alignment.append((local_index, full, barcode))
    return {
        "read_rows": read_rows,
        "synthetic": synthetic,
        "alignment": alignment,
    }


def _artifact_worker_loop(task_queue, result_queue) -> None:
    accumulator = _new_accumulator()
    prior_error = None
    while True:
        item = task_queue.get()
        if item is None:
            break
        batch_id, records = item
        if prior_error is not None:
            result_queue.put(("error", batch_id, prior_error))
            continue
        try:
            result = _accumulate_artifact_batch(accumulator, records)
        except Exception as exc:
            prior_error = {
                "error_type": type(exc).__name__,
                "error": str(exc),
                "traceback": traceback.format_exc(),
            }
            result_queue.put(("error", batch_id, prior_error))
        else:
            result_queue.put(("batch", batch_id, result))
    result_queue.put(
        ("final", None, accumulator if prior_error is None else prior_error)
    )


def aggregate_artifact_batches(
    batches,
    *,
    workers: int,
    consume_batch,
) -> ScanAccumulator:
    """Process numbered batches with bounded raw-record retention.

    ``batches`` yields ``(batch_id, task_records, raw_records)``.  Only compact
    sequence tuples cross to workers; the parent retains at most ``2*workers``
    raw batches and writes selected artifacts after the matching batch result
    returns.  Completed batches are consumed in input order.
    """

    if workers < 1:
        raise ValueError("workers must be positive")
    if workers == 1:
        accumulator = _new_accumulator()
        for batch_id, task_records, raw_records in batches:
            result = _accumulate_artifact_batch(accumulator, task_records)
            consume_batch(batch_id, result, raw_records)
        return accumulator

    context = multiprocessing.get_context("fork")
    maximum_pending = 2 * workers
    task_queue = context.Queue(maxsize=maximum_pending)
    result_queue = context.Queue()
    processes = [
        context.Process(target=_artifact_worker_loop, args=(task_queue, result_queue))
        for _ in range(workers)
    ]
    pending_raw = {}
    completed = {}
    next_to_consume = 0
    source = iter(batches)
    source_exhausted = False
    try:
        for process in processes:
            process.start()
        while not source_exhausted or pending_raw:
            while not source_exhausted and len(pending_raw) < maximum_pending:
                try:
                    batch_id, task_records, raw_records = next(source)
                except StopIteration:
                    source_exhausted = True
                    break
                if batch_id in pending_raw or batch_id in completed:
                    raise RuntimeError(f"duplicate artifact batch ID: {batch_id}")
                pending_raw[batch_id] = raw_records
                task_queue.put((batch_id, task_records))
            if not pending_raw:
                continue
            status, batch_id, payload = result_queue.get()
            if status != "batch":
                raise RuntimeError(
                    f"parallel QC artifact worker failure for batch {batch_id}: {payload}"
                )
            completed[batch_id] = payload
            while next_to_consume in completed:
                consume_batch(
                    next_to_consume,
                    completed.pop(next_to_consume),
                    pending_raw.pop(next_to_consume),
                )
                next_to_consume += 1

        for _ in processes:
            task_queue.put(None)
        final_results = []
        while len(final_results) < len(processes):
            status, _, payload = result_queue.get()
            if status != "final":
                raise RuntimeError(f"unexpected QC worker message: {status}")
            final_results.append(payload)
        for process in processes:
            process.join()
        abnormal = [process.exitcode for process in processes if process.exitcode != 0]
        errors = [item for item in final_results if isinstance(item, dict)]
        if abnormal or errors:
            raise RuntimeError(
                f"parallel QC worker failure: errors={errors}, exitcodes={abnormal}"
            )
        merged = _new_accumulator()
        for accumulator in final_results:
            merged.merge(accumulator)
        return merged
    finally:
        for process in processes:
            if process.is_alive():
                process.terminate()
            process.join()
        task_queue.close()
        result_queue.close()
        task_queue.join_thread()
        result_queue.join_thread()


def _filter_count_rows(confusion_by_truth, library_name: str) -> list[dict]:
    rows = []
    for truth_name, methods in confusion_by_truth.items():
        for method, confusion in methods.items():
            rows.append(
                {
                    "library": library_name,
                    "truth_definition": truth_name,
                    "method": method,
                    "tp": int(confusion.tp),
                    "fp": int(confusion.fp),
                    "tn": int(confusion.tn),
                    "fn": int(confusion.fn),
                    "hard_negative_n": int(confusion.hard_negative_n),
                    "hard_negative_fp": int(confusion.hard_negative_fp),
                }
            )
    rows.sort(key=lambda row: (row["truth_definition"], row["method"]))
    return rows


def _cell_rows(molecule_estimator) -> list[dict]:
    if molecule_estimator is None:
        return []
    rows = []
    for cell, counts in molecule_estimator.cell_read_counts.items():
        all_reads = int(counts["all"])
        tso_reads = int(counts["tso"])
        contaminants = int(counts["structure_contaminant"])
        rows.append(
            {
                "cell": cell,
                "valid_barcode_reads": all_reads,
                "tso_reads": tso_reads,
                "structure_contaminant_reads": contaminants,
                "contaminant_fraction_of_all": (
                    contaminants / all_reads if all_reads else 0.0
                ),
                "contaminant_fraction_of_tso": (
                    contaminants / tso_reads if tso_reads else 0.0
                ),
            }
        )
    rows.sort(key=lambda row: row["cell"])
    return rows


def _truth_summaries(truth_counts: Counter, records_evaluated: int) -> tuple[dict, dict]:
    contamination_count = int(truth_counts["structure_truth"])
    tso_count = int(truth_counts["tso_found"])
    ci_low, ci_high = qc.wilson_interval(contamination_count, records_evaluated)
    tso_contam_low, tso_contam_high = qc.wilson_interval(
        contamination_count,
        tso_count,
    )
    raw_prevalence = {
        "all_reads": records_evaluated,
        "tso_reads": tso_count,
        "structure_contaminant_reads": contamination_count,
        "contaminant_fraction_of_all_reads": (
            contamination_count / records_evaluated if records_evaluated else 0.0
        ),
        "contaminant_fraction_of_all_reads_wilson_95": [ci_low, ci_high],
        "contaminant_fraction_of_tso_reads": (
            contamination_count / tso_count if tso_count else 0.0
        ),
        "contaminant_fraction_of_tso_reads_wilson_95": [
            tso_contam_low,
            tso_contam_high,
        ],
    }
    confidence_tiers = {}
    for tier_name, count_key, interpretation in (
        (
            "definitive_reference_construct",
            "complete_reference_truth",
            "TSO + structure + known spacer + compatible reference iBAR",
        ),
        (
            "reference_supported_structure",
            "known_structure_truth",
            "TSO + complete structure + known spacer",
        ),
        (
            "primary_structural_surrogate",
            "structure_truth",
            "TSO + spacer-sized segment + paired stem/iBAR-sized block/tracr",
        ),
        (
            "inclusive_possible_contaminant",
            "inclusive_truth",
            "Primary structural call or TSO followed by a known spacer",
        ),
    ):
        tier_count = int(truth_counts[count_key])
        tier_low, tier_high = qc.wilson_interval(tier_count, records_evaluated)
        confidence_tiers[tier_name] = {
            "count": tier_count,
            "fraction": tier_count / records_evaluated if records_evaluated else 0.0,
            "wilson_95_low": tier_low,
            "wilson_95_high": tier_high,
            "interpretation": interpretation,
        }
    excluded = {
        "structure_truth",
        "known_structure_truth",
        "complete_reference_truth",
        "tso_known_guide_truth",
        "inclusive_truth",
        "tso_found",
        "structure_selected_filter_call",
        "structure_predicted_splitpipe_tso_trim",
        "structure_splitpipe_trim_known_guide_d2",
        "hard_negative_selected_filter_call",
    }
    full_truth = {
        "structure_contaminants": contamination_count,
        "fraction": contamination_count / records_evaluated if records_evaluated else 0.0,
        "wilson_95_low": ci_low,
        "wilson_95_high": ci_high,
        "known_structure_contaminants": int(truth_counts["known_structure_truth"]),
        "complete_reference_contaminants": int(
            truth_counts["complete_reference_truth"]
        ),
        "tso_known_guide_contaminants": int(
            truth_counts["tso_known_guide_truth"]
        ),
        "inclusive_contaminants": int(truth_counts["inclusive_truth"]),
        "tso_positive_reads": tso_count,
        "confidence_tiers": confidence_tiers,
        "tier_counts": {
            key: int(value)
            for key, value in truth_counts.items()
            if key not in excluded
        },
    }
    return raw_prevalence, full_truth


def scan_library(
    *,
    cfg: dict,
    config_path: str,
    library_name: str,
    pairs: Sequence[Tuple[str, str]],
    workers: int,
    batch_size: int,
    max_reads: int,
) -> dict:
    """Run the production scan and return a serial-comparable fingerprint."""
    global _MODEL, _BARCODE_CORRECTOR, _READ_LENGTH, _SHORT_CONFIG
    global _SELECTED_METHOD, _MOLECULE_ENABLED, _HLL_PRECISION
    global _BASE_COMPOSITION_LENGTH
    if workers < 1:
        raise ValueError("workers must be positive")
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    if not pairs:
        raise ValueError("at least one explicit FASTQ pair is required")

    started = time.perf_counter()
    qc_config = cfg["gex_grna_contamination_qc"]
    raw_config = qc_config["raw_benchmark"]
    if float(raw_config.get("sample_fraction", 1.0)) != 1.0:
        raise ValueError("parallel production scan requires sample_fraction=1.0")
    short_config = qc_config["short_read_filters"]
    molecule_config = qc_config.get("molecule_estimation", {})
    alignment_config = qc_config.get("alignment_benchmark", {})
    molecule_enabled = bool(molecule_config.get("enabled", True))
    need_barcodes = molecule_enabled or bool(alignment_config.get("enabled", True))

    _MODEL = qc.make_model(cfg, config_path)
    _BARCODE_CORRECTOR = (
        make_barcode_corrector_fast(
            cfg,
            config_path,
            library_name,
            [r2_path for _, r2_path in pairs],
        )
        if need_barcodes
        else None
    )
    _READ_LENGTH = int(qc_config.get("short_read_length", 54))
    _SHORT_CONFIG = short_config
    _SELECTED_METHOD = str(
        short_config.get("selected_method", "tso_guide_d1")
    )
    _MOLECULE_ENABLED = molecule_enabled
    _HLL_PRECISION = int(molecule_config.get("hll_precision", 18))
    _BASE_COMPOSITION_LENGTH = int(
        raw_config.get("base_composition_length", 100)
    )
    scan_started = time.perf_counter()
    batches = iter_sequence_batches(
        pairs,
        max_reads=max_reads,
        batch_size=batch_size,
    )
    accumulator = aggregate_batches(batches, workers=workers)
    records_evaluated = accumulator.records_evaluated
    truth_counts = accumulator.truth_counts
    confusion_by_truth = accumulator.confusion_by_truth
    molecule_estimator = accumulator.molecule_estimator
    if max_reads and records_evaluated != max_reads:
        raise RuntimeError(
            f"parallel scan evaluated {records_evaluated:,}, expected {max_reads:,}"
        )

    raw_prevalence, full_truth = _truth_summaries(
        truth_counts,
        records_evaluated,
    )
    molecule_summary = (
        molecule_estimator.summary() if molecule_estimator is not None else None
    )
    short_hazard = {
        "structure_truth_with_exact_splitpipe_tso_trim": int(
            truth_counts["structure_predicted_splitpipe_tso_trim"]
        ),
        "structure_truth_with_splitpipe_trim_and_known_guide_d2": int(
            truth_counts["structure_splitpipe_trim_known_guide_d2"]
        ),
        "structure_truth_called_by_selected_filter": int(
            truth_counts["structure_selected_filter_call"]
        ),
        "tso_positive_structure_negative_called_by_selected_filter": int(
            truth_counts["hard_negative_selected_filter_call"]
        ),
    }
    barcode_summary = None
    if _BARCODE_CORRECTOR is not None:
        barcode_summary = {
            "perfect_combination_sample_count": int(
                sum(_BARCODE_CORRECTOR.perfect_combination_counts.values())
            ),
            "perfect_combination_count_threshold": int(
                _BARCODE_CORRECTOR.perfect_count_threshold
            ),
            "maximum_barcode_mismatches": _BARCODE_CORRECTOR.max_mismatches,
        }
    canonical_values = {
        "records_evaluated": records_evaluated,
        "raw_read_prevalence": raw_prevalence,
        "full_length_surrogate_truth": full_truth,
        "short_read_hazard": short_hazard,
        "selected_filter": _SELECTED_METHOD,
        "filter_confusion_counts": _filter_count_rows(
            confusion_by_truth,
            library_name,
        ),
        "barcode_correction": barcode_summary,
        "molecule_summary": molecule_summary,
        "cell_read_contamination": _cell_rows(molecule_estimator),
    }
    scan_elapsed = time.perf_counter() - scan_started
    elapsed = time.perf_counter() - started
    index_stats = {
        "guide": _MODEL.guide_index.stats.__dict__,
        "template_prefix": _MODEL.template_prefix_index.stats.__dict__,
        "barcodes": {
            str(round_number): index.stats.__dict__
            for round_number, index in (
                _BARCODE_CORRECTOR.index_by_round.items()
                if _BARCODE_CORRECTOR is not None
                else []
            )
        },
    }
    return {
        "records": records_evaluated,
        "workers": workers,
        "batch_size": batch_size,
        "setup_seconds": scan_started - started,
        "scan_wall_seconds": scan_elapsed,
        "scan_read_pairs_per_second": records_evaluated / scan_elapsed,
        "wall_seconds": elapsed,
        "read_pairs_per_second": records_evaluated / elapsed,
        "parent_peak_rss_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        / 1024.0,
        "children_peak_rss_mib": resource.getrusage(
            resource.RUSAGE_CHILDREN
        ).ru_maxrss
        / 1024.0,
        "maximum_pending_batches": 1 if workers == 1 else 2 * workers,
        "maximum_pending_records": batch_size * (1 if workers == 1 else 2 * workers),
        "worker_local_aggregation": True,
        "index_stats": index_stats,
        "canonical_results": {
            "schema_version": 1,
            "values": canonical_values,
            "sha256": canonical_json_sha256(canonical_values),
        },
    }

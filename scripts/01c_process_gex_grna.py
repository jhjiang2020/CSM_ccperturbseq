#!/usr/bin/env python3
"""Run ordinary Step 3 on barcode evidence emitted by Step 1a.

Step 1a corrects cell barcodes against the full GEX calibration interval and
stores valid filter-positive rows in the ordinary Step-3 input schema.  Step 1c
must consume that exact evidence; recalibrating from the small rejected subset
is deliberately unsupported because it is not scientifically equivalent.
"""

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _manifest_fastq_path(manifest_path, payload, key):
    raw = payload.get(key)
    if not raw:
        return None
    path = Path(raw).expanduser()
    if not path.is_absolute():
        path = manifest_path.parent / path
    return path.resolve()


def warn_if_ibar_unobservable(config, manifest_path, sample_records=10_000):
    """Warn when Step-1a evidence R1s cannot contain a complete iBAR.

    This is diagnostic only. Step 1c still invokes the ordinary Step 3 code,
    which remains authoritative and will not invent a missing iBAR.
    """
    try:
        from scripts.guide_design import enabled
    except ModuleNotFoundError:
        from guide_design import enabled
    if enabled(config):
        return None
    structure = config.get("grna_structure", {}) or {}
    assignment = config.get("guide_assignment", {}) or {}
    ibar_length = int(structure.get("ibar_length", 0))
    if ibar_length <= 0 or not assignment.get("dual_guide_counting_mode", False):
        return None

    guide_lengths = [int(value) for value in structure.get("grna_lengths", [])]
    stems = [
        str(structure.get(key, ""))
        for key in ("fixed_seq_2", "fixed_seq_3")
        if structure.get(key)
    ]
    tso = str(structure.get("fixed_seq_1", ""))
    if not guide_lengths or not stems or not tso:
        return None
    full_ibar_minimum = len(tso) + min(guide_lengths) + min(map(len, stems)) + ibar_length

    payload = json.loads(Path(manifest_path).read_text())
    lengths = []
    evidence_path = _manifest_fastq_path(
        Path(manifest_path), payload, "barcode_evidence"
    )
    if evidence_path is not None and evidence_path.is_file():
        with gzip.open(evidence_path, "rt", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                sequence = str(row.get("r1_seq", ""))
                quality = str(row.get("r1_qual", ""))
                if not sequence or len(sequence) != len(quality):
                    break
                lengths.append(len(sequence))
                if len(lengths) >= sample_records:
                    break
    if lengths and max(lengths) < full_ibar_minimum:
        message = (
            "WARNING: all sampled GEX-derived gRNA R1 reads are shorter than "
            f"the approximate full-iBAR requirement ({max(lengths)} < "
            f"{full_ibar_minimum} bp). Step 1a removal remains valid, but "
            "ordinary dual-guide Step 3 will generally reject these reads; "
            "no guide-only dual identity will be invented."
        )
        print(message, file=sys.stderr)
        return message
    return None


def output_root(config, namespace):
    root = Path(config["paths"]["output_dir"])
    candidate = Path(namespace).expanduser()
    return candidate if candidate.is_absolute() else root / candidate


def load_barcode_evidence(manifest_path, library):
    manifest_path = Path(manifest_path).expanduser().resolve()
    payload = json.loads(manifest_path.read_text())
    if payload.get("complete") is not True or payload.get("status") != "complete":
        raise ValueError(f"Step 1a manifest is not complete: {manifest_path}")
    if payload.get("schema_version") != 4:
        raise ValueError(
            "Step 1c requires Step 1a manifest schema 4 with corrected "
            "barcode evidence"
        )
    if payload.get("library") != library:
        raise ValueError(
            f"Step 1a manifest library {payload.get('library')!r} does not "
            f"match requested library {library!r}"
        )
    raw = payload.get("barcode_evidence")
    if not raw:
        raise ValueError("Step 1a manifest has no barcode_evidence artifact")
    evidence = Path(raw).expanduser()
    if not evidence.is_absolute():
        evidence = manifest_path.parent / evidence
    evidence = evidence.resolve()
    if not evidence.is_file():
        raise ValueError(f"Step 1a barcode evidence does not exist: {evidence}")
    expected_hash = payload.get("barcode_evidence_sha256")
    observed_hash = sha256_file(evidence)
    if not expected_hash or observed_hash != expected_hash:
        raise ValueError("Step 1a barcode-evidence SHA-256 mismatch")
    expected_rows = int(payload.get("barcode_evidence_rows", -1))
    stat_rows = int(
        (payload.get("statistics", {}) or {}).get(
            "csm_gex_filter_positive_valid_barcode_reads", -2
        )
    )
    if expected_rows < 0 or expected_rows != stat_rows:
        raise ValueError(
            "Step 1a barcode-evidence count does not match valid-positive statistics"
        )
    with gzip.open(evidence, "rt") as handle:
        header = handle.readline().rstrip("\n")
        if header != "read_id\tbc_wells\tbc_bci\tumi\tr1_seq\tr1_qual":
            raise ValueError("Step 1a barcode evidence has an incompatible header")
        observed_rows = sum(1 for line in handle if line.rstrip("\n"))
    if observed_rows != expected_rows:
        raise ValueError(
            f"Step 1a barcode-evidence rows {observed_rows} != {expected_rows}"
        )
    return evidence, payload


def build_commands(script_dir, python, config_path, library, manifest,
                   namespace, n_processes, max_reads=0, chunk_size=100_000):
    if max_reads:
        raise ValueError("Step 1c no longer supports truncating Step 1a evidence")
    barcode_input, evidence_manifest = load_barcode_evidence(manifest, library)
    evidence_rows = int(evidence_manifest["barcode_evidence_rows"])
    # Do not start idle worker processes for the usual ~30-CPM input. Larger
    # projects still scale by one worker per Step-3 chunk, capped by NSLOTS.
    effective_processes = min(
        max(1, int(n_processes)),
        max(1, math.ceil(evidence_rows / max(1, int(chunk_size)))),
    )
    namespace_root = output_root(yaml.safe_load(config_path.read_text()), namespace)
    step3 = [
        python, str(script_dir / "03_match_grna_ibar.py"),
        "--config", str(config_path),
        "--library", library,
        "--barcode-input", str(barcode_input),
        "--output-namespace", namespace,
        "--n_processes", str(effective_processes),
        "--chunk_size", str(chunk_size),
    ]
    return step3, namespace_root, barcode_input, effective_processes


def atomic_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def install_strict_match_product(match_output, config, library):
    """Preserve raw Step-3 matches and install a conservative merge product."""
    step1c_cfg = ((config.get("gex_decontamination", {}) or {}).get("step1c", {}) or {})
    strict = step1c_cfg.get("strict_matching", {}) or {}
    if not bool(strict.get("enabled", True)):
        return Path(match_output), None, {"enabled": False}

    raw_path = Path(match_output)
    unfiltered = raw_path.with_name(f"{library}_grna_matches_unfiltered.tsv.gz")
    if unfiltered.exists():
        raise RuntimeError(f"refusing existing Step 1c unfiltered product: {unfiltered}")
    os.replace(raw_path, unfiltered)
    temporary = raw_path.with_name(f".{raw_path.name}.tmp.{os.getpid()}")
    max_guide = int(strict.get("max_guide_edit_distance", 1))
    max_ibar = int(strict.get("max_ibar_edit_distance", 1))
    dual_mode = bool((config.get("guide_assignment", {}) or {}).get("dual_guide_counting_mode", False))
    require_ibar = bool(strict.get("require_resolved_ibar_in_dual_mode", True))
    reject_decoupled = bool(strict.get("reject_ibar_decoupled", True))
    counts = {
        "enabled": True,
        "max_guide_edit_distance": max_guide,
        "max_ibar_edit_distance": max_ibar,
        "require_resolved_ibar": bool(dual_mode and require_ibar),
        "reject_ibar_decoupled": bool(dual_mode and reject_decoupled),
        "input_matches": 0,
        "retained_matches": 0,
        "dropped_guide_distance": 0,
        "dropped_ibar_distance": 0,
        "dropped_unresolved_ibar": 0,
        "dropped_ibar_decoupled": 0,
        "dropped_structure_conflict": 0,
    }
    with gzip.open(unfiltered, "rt", newline="") as source, gzip.open(
        temporary, "wt", newline=""
    ) as destination:
        reader = csv.DictReader(source, delimiter="\t")
        required = {"grna_dist", "ibar_dist", "ibar_seq_corrected"}
        if dual_mode and reject_decoupled:
            required.add("ibar_decoupled")
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise RuntimeError(f"ordinary Step 3 output lacks strict-match fields: {sorted(missing)}")
        writer = csv.DictWriter(
            destination,
            fieldnames=reader.fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in reader:
            counts["input_matches"] += 1
            if row.get('structure_status', 'concordant') != 'concordant':
                counts['dropped_structure_conflict'] += 1
                continue
            guide_distance = int(row["grna_dist"])
            ibar_distance = int(row["ibar_dist"])
            resolved_ibar = str(row.get("ibar_seq_corrected", "")).strip()
            ibar_decoupled = str(row.get("ibar_decoupled", "")).strip().lower() in {
                "1", "true", "yes",
            }
            if guide_distance < 0 or guide_distance > max_guide:
                counts["dropped_guide_distance"] += 1
                continue
            absent_by_design = row.get('ibar_status') == 'absent_by_design' and row.get('structure_status') == 'concordant'
            if dual_mode and require_ibar and not resolved_ibar and not absent_by_design:
                counts["dropped_unresolved_ibar"] += 1
                continue
            if dual_mode and reject_decoupled and ibar_decoupled:
                counts["dropped_ibar_decoupled"] += 1
                continue
            if resolved_ibar and (ibar_distance < 0 or ibar_distance > max_ibar):
                counts["dropped_ibar_distance"] += 1
                continue
            writer.writerow(row)
            counts["retained_matches"] += 1
    os.replace(temporary, raw_path)
    return raw_path, unfiltered, counts


def main():
    parser = argparse.ArgumentParser(
        description="Step 1c: conservatively identify guides in Step 1a evidence"
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--library", required=True)
    parser.add_argument("--input-manifest", required=True)
    parser.add_argument("--output-namespace", default="gex_grna_identity")
    parser.add_argument("--n_processes", type=int, default=None)
    parser.add_argument("--max_reads", type=int, default=0)
    parser.add_argument("--chunk_size", type=int, default=100_000)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    config_path = Path(args.config).expanduser().resolve()
    manifest = Path(args.input_manifest).expanduser().resolve()
    with config_path.open() as handle:
        config = yaml.safe_load(handle)
    configured_libraries = {
        row["sublibrary_name"] for row in config.get("libraries", [])
    }
    if args.library not in configured_libraries:
        parser.error(f"library {args.library!r} is not present in the config")
    if not manifest.is_file():
        parser.error(f"input manifest does not exist: {manifest}")

    n_processes = args.n_processes or config["barcode_extraction"].get(
        "n_processes", 1
    )
    script_dir = Path(__file__).resolve().parent
    commands = build_commands(
        script_dir, sys.executable, config_path, args.library, manifest,
        args.output_namespace, n_processes, args.max_reads, args.chunk_size,
    )
    step3, namespace_root, barcode_input, effective_processes = commands
    completion_path = (
        namespace_root / "manifests" / f"{args.library}.step1c.complete.json"
    )

    if args.dry_run:
        print(" ".join(step3))
        return
    if completion_path.exists():
        raise RuntimeError(
            f"refusing stale/existing Step 1c completion record: {completion_path}"
        )

    print(f"Step 1c source library: {args.library}")
    print(f"Step 1a evidence manifest: {manifest}")
    print(
        f"Step 1c workers: {effective_processes} effective "
        f"({n_processes} scheduler slots available)"
    )
    warn_if_ibar_unobservable(config, manifest)
    subprocess.run(step3, check=True)

    barcode_output = barcode_input
    match_output = (
        namespace_root / "grna_matches" / f"{args.library}_grna_matches.tsv.gz"
    )
    missing = [path for path in (barcode_output, match_output) if not path.is_file()]
    if missing:
        raise RuntimeError(
            "Step 1c child commands returned successfully but outputs are missing: "
            + ", ".join(str(path) for path in missing)
        )
    match_output, unfiltered_match_output, strict_statistics = install_strict_match_product(
        match_output, config, args.library
    )

    step3_stats_path = (
        namespace_root / "grna_matches" / f"{args.library}_matching_stats.yaml"
    )
    statistics = {}
    for stage, stats_path in (("step3", step3_stats_path),):
        if stats_path.is_file():
            with stats_path.open() as handle:
                stage_stats = yaml.safe_load(handle) or {}
            statistics[stage] = stage_stats.get(args.library, stage_stats)

    completion = {
        "schema_version": 2,
        "status": "complete",
        "complete": True,
        "source": "gex_step1a_corrected_barcode_evidence",
        "source_library": args.library,
        "step1a_manifest": str(manifest),
        "config": str(config_path),
        "config_sha256": sha256_file(config_path),
        "step1a_manifest_sha256": sha256_file(manifest),
        "output_namespace": str(namespace_root.resolve()),
        "barcode_output": str(barcode_output.resolve()),
        "barcode_output_sha256": sha256_file(barcode_output),
        "barcode_evidence_rows": int(
            json.loads(manifest.read_text())["barcode_evidence_rows"]
        ),
        "scheduler_slots": int(n_processes),
        "effective_workers": int(effective_processes),
        "match_output": str(match_output.resolve()),
        "match_output_sha256": sha256_file(match_output),
        "match_output_size_bytes": match_output.stat().st_size,
        "unfiltered_match_output": (
            str(unfiltered_match_output.resolve())
            if unfiltered_match_output is not None
            else None
        ),
        "unfiltered_match_output_sha256": (
            sha256_file(unfiltered_match_output)
            if unfiltered_match_output is not None
            else None
        ),
        "strict_matching": strict_statistics,
        "step3_stats": str(step3_stats_path.resolve()),
        "statistics": statistics,
        "completed_at": datetime.now(timezone.utc).isoformat(),
    }
    atomic_json(completion_path, completion)
    print(f"Step 1c completion manifest: {completion_path}")


if __name__ == "__main__":
    main()

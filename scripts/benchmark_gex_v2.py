#!/usr/bin/env python3
"""Cached 50M-pair four-library Step 1a benchmark. Run data stages on compute nodes.

The native STAR unique path retains Parse's molecule/DGE implementations.
Only its repeated get_gene call reads cached GX/RE tags instead of reannotating.
STAR fractional and kallisto use the production shared molecule counter.
"""
from __future__ import annotations

import argparse
import ast
import copy
import datetime
import gzip
import hashlib
import json
import os
import pickle
import runpy
import shutil
import socket
import sqlite3
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

import yaml

try:
    from scripts.gex_alignment import atomic_json, backend_options, sha256_file, load_reference, star_reference_identity, write_star_parameters
    from scripts.kallisto_umi import fastq_records
    from scripts.qc_fastq_resolver import resolve_all_library_pairs
except ModuleNotFoundError:
    from gex_alignment import atomic_json, backend_options, sha256_file, load_reference, star_reference_identity, write_star_parameters
    from kallisto_umi import fastq_records
    from qc_fastq_resolver import resolve_all_library_pairs

LIBRARIES = ("C1", "C2", "T1", "T2")
BACKENDS = ("star", "kallisto")
ASSIGNMENTS = ("unique", "fractional")
UMIS = ("parse", "directional", "none")
ARMS = tuple(f"{b}_{a}_{u}" for b in BACKENDS for a in ASSIGNMENTS for u in UMIS)
SCRIPT = Path(__file__).resolve()
SCRIPTS = SCRIPT.parent
PATCH = SCRIPTS.parent / "patches/parsebio_splitpipe_1.6.4_gex_decontamination.patch"


def arm_options(arm):
    if arm not in ARMS:
        raise ValueError(f"Unknown arm: {arm}")
    return tuple(arm.split("_"))


def read_json(path):
    return json.loads(Path(path).read_text())


def dump_pickle(path, value):
    with Path(path).open("wb") as handle:
        pickle.dump(value, handle, protocol=pickle.HIGHEST_PROTOCOL)


def load_pickle(path):
    # Only benchmark-generated local caches, never user-supplied serialized objects.
    with Path(path).open("rb") as handle:
        return pickle.load(handle)


def stage_read_database(source, directory):
    """Stage a sealed read spool onto job-local storage, as production does."""
    source, directory = Path(source), Path(directory)
    destination = directory / "reads.sqlite"
    if destination.exists():
        raise FileExistsError(destination)
    before = source.stat()
    started = time.perf_counter()
    shutil.copyfile(source, destination)
    after = source.stat()
    if ((before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns)
            or destination.stat().st_size != before.st_size):
        raise ValueError("Immutable read cache changed during local staging")
    return destination, {"cache_stage_seconds": time.perf_counter() - started,
                         "cache_stage_bytes": before.st_size, "cache_storage": "job_local_copy_v2"}


def fingerprint(path):
    path = Path(path)
    stat = path.stat()
    return {"path": str(path.resolve()), "bytes": stat.st_size,
            "mtime_ns": stat.st_mtime_ns, "sha256": sha256_file(path)}


def validate_artifacts(record, *, full=False):
    if record.get("complete") is not True:
        raise ValueError("Incomplete benchmark cache")
    for artifact in record["artifacts"]:
        path = Path(artifact["path"])
        stat = path.stat()
        if stat.st_size != artifact["bytes"] or stat.st_mtime_ns != artifact["mtime_ns"]:
            raise ValueError(f"Changed benchmark artifact: {path}")
        if full and sha256_file(path) != artifact["sha256"]:
            raise ValueError(f"Checksum mismatch: {path}")


def config_for(args, arm=None):
    cfg = copy.deepcopy(yaml.safe_load(Path(args.config).read_text()))
    if tuple(lib["sublibrary_name"] for lib in cfg["libraries"]) != LIBRARIES:
        raise ValueError("Benchmark requires the ordered four-library pilot: C1,C2,T1,T2")
    # Resolve input paths before writing derived YAMLs in the run directory.
    for key, value in cfg["paths"].items():
        if isinstance(value, str) and key.endswith(("_dir", "_reference")):
            cfg["paths"][key] = str((Path(args.config).resolve().parent / value).resolve())
    decontam = cfg.setdefault("gex_decontamination", {})
    for name in ("existing_alignment", "step1c", "combine_with_grna_arm"):
        decontam.setdefault(name, {})["enabled"] = False
    b, a, u = arm_options(arm or "star_unique_parse")
    cfg["gex_alignment"] = {"backend": b, "gene_assignment": a, "umi_correction": u,
                            "output_level": "gene", "count_threads": 0,
                            "star": {"max_multimaps": args.max_multimaps},
                            "kallisto": {"reference_dir": str(Path(args.reference).resolve()),
                                         "executable": args.kallisto, "strand": "forward"}}
    cfg["paths"]["output_dir"] = str(args.work_dir / (arm or "preprocess"))
    cfg["paths"]["fastq_dir"] = str(args.work_dir / "sampled_fastq")
    cfg["paths"]["log_dir"] = str(args.work_dir / "logs")
    return cfg


def benchmark_sources(directory=None):
    """Seal local imports and explicitly executed helpers, excluding other workflows."""
    directory = Path(directory) if directory is not None else SCRIPTS
    pending = {"benchmark_gex_v2.py", "benchmark_gex_v2_report.py", "splitpipe_gex_wrapper.py",
               "write_fastq_input_identity.py", "gex_grna_preprocess_hook.py"}
    found = set()
    while pending:
        name = pending.pop()
        if name in found:
            continue
        found.add(name)
        tree = ast.parse((directory / name).read_text(), filename=name)
        candidates = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                candidates.update(alias.name.split(".")[-1] + ".py" for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                candidates.add((node.module or "").split(".")[-1] + ".py")
                candidates.update(alias.name + ".py" for alias in node.names)
            elif isinstance(node, ast.Constant) and isinstance(node.value, str) and node.value.endswith(".py"):
                candidates.add(Path(node.value).name)
        pending.update(n for n in candidates if (directory / n).is_file() and n not in found)
    return [directory / name for name in sorted(found)]


def contract(args):
    cfg = config_for(args)
    # Scientific settings stay sealed; repairs require an explicit recorded revision.
    vendor = Path(cfg["paths"]["splitpipe_dir"])
    files = benchmark_sources() + [PATCH, vendor / "split-pipe"]
    files += [vendor / "splitpipe" / name for name in ("preprocess.py", "align.py", "molinfo.py", "dge.py")]
    digest = hashlib.sha256()
    for path in files:
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    payload = {"schema": 2, "pairs": args.pairs, "threads": args.threads,
               "config": cfg, "code_sha256": digest.hexdigest(), "arms": list(ARMS),
               "reference_manifest": fingerprint(Path(args.reference) / "reference.json"),
               "star_annotation_sha256": star_reference_identity(cfg, args.config),
               "kallisto": fingerprint(shutil.which(args.kallisto) or args.kallisto)}
    path = args.work_dir / "experiment.json"
    # Launcher creates the contract before concurrent arrays start.
    if path.exists():
        previous = read_json(path)
        if previous != payload:
            old_scientific = {k: v for k, v in previous.items() if k != "code_sha256"}
            new_scientific = {k: v for k, v in payload.items() if k != "code_sha256"}
            if args.mode != "accept-revision" or not args.repair_note or old_scientific != new_scientific:
                raise ValueError("Experiment settings/code changed; record a reviewed repair revision or choose a new directory")
            directory = args.work_dir / "code_revisions" / payload["code_sha256"][:16]
            directory.mkdir(parents=True, exist_ok=False)
            atomic_json(directory / "previous_experiment.json", previous)
            source_hashes = {}
            for source in files:
                relative = source.relative_to(SCRIPTS.parent) if source.is_relative_to(SCRIPTS.parent) else Path("vendor") / source.name
                target = directory / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
                source_hashes[str(relative)] = sha256_file(source)
            event = {"utc": datetime.datetime.now(datetime.timezone.utc).isoformat(), "reason": args.repair_note,
                     "previous_code_sha256": previous["code_sha256"], "code_sha256": payload["code_sha256"],
                     "scientific_settings_unchanged": True, "source_hashes": source_hashes}
            atomic_json(directory / "revision.json", event)
            with (args.work_dir / "code_revisions.jsonl").open("a") as handle:
                handle.write(json.dumps(event, sort_keys=True) + "\n")
            atomic_json(path, payload)
    else:
        if args.mode != "init":
            raise ValueError("Run init once before scheduling benchmark stages")
        load_reference(Path(args.reference))
        atomic_json(path, payload)
    return cfg


def cfg_path(args, arm=None):
    path = args.work_dir / "configs" / f"{args.library}.{arm or 'preprocess'}.yaml"
    path.parent.mkdir(parents=True, exist_ok=True)
    atomic = path.with_suffix(f".{os.getpid()}.tmp")
    atomic.write_text(yaml.safe_dump(config_for(args, arm), sort_keys=False))
    atomic.replace(path)
    return path


def manifest(args, name):
    return args.work_dir / "manifests" / f"{args.library}.{name}.json"


def require(args, name):
    record = read_json(manifest(args, name))
    validate_artifacts(record, full=args.verify_cache)
    return record


def stage(args, name, command, outputs, dependencies=(), env=None):
    """Per-child wait4 accounting per stage; cached work is never timed again."""
    parents = {d: sha256_file(manifest(args, d)) for d in dependencies}
    for dep in dependencies:
        require(args, dep)
    path = manifest(args, name)
    if path.exists():
        old = require(args, name)
        if old["dependencies"] != parents:
            raise ValueError(f"Stale dependency for {name}")
        return old
    # A partial stage is not resumable: prevents vendor auto-skip contaminating timings.
    started_path = path.with_suffix(".started")
    if started_path.exists():
        raise FileExistsError(f"Partial stage {name}; inspect logs and use a new work directory")
    started_path.touch()
    logs = args.work_dir / "logs"
    logs.mkdir(exist_ok=True)
    usage = logs / f"{args.library}.{name}.usage"
    log = logs / f"{args.library}.{name}.log"
    timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
    before = time.perf_counter()
    peaks = {"peak_tree_rss_gb": 0., "peak_temp_bytes": 0, "peak_working_file_bytes": 0}
    working_paths = outputs if name == "sample" else {Path(p).parent for p in outputs}
    with tempfile.TemporaryDirectory(prefix="csm-bench-stage-") as temporary, log.open("w") as handle:
        child_env = dict(os.environ if env is None else env, TMPDIR=temporary)
        process = subprocess.Popen(command, stdout=handle, stderr=subprocess.STDOUT, env=child_env)
        stop = threading.Event()
        def sample_resources():
            while not stop.is_set():
                peaks["peak_tree_rss_gb"] = max(peaks["peak_tree_rss_gb"], tree_rss(process.pid) / 1024**3)
                peaks["peak_temp_bytes"] = max(peaks["peak_temp_bytes"], tree_bytes([Path(temporary)]))
                peaks["peak_working_file_bytes"] = max(peaks["peak_working_file_bytes"], tree_bytes(working_paths))
                stop.wait(2)
        sampler = threading.Thread(target=sample_resources, daemon=True)
        sampler.start()
        # Unlike getrusage(RUSAGE_CHILDREN), wait4 returns this child's usage
        # (and descendants it waited for), not earlier stages in this worker.
        _, status, child_usage = os.wait4(process.pid, 0)
        exit_code = os.waitstatus_to_exitcode(status)
        process.returncode = exit_code
        stop.set()
        sampler.join()
        peaks["peak_temp_bytes"] = max(peaks["peak_temp_bytes"], tree_bytes([Path(temporary)]))
        peaks["peak_working_file_bytes"] = max(peaks["peak_working_file_bytes"], tree_bytes(working_paths))
    wall = time.perf_counter() - before
    metrics = {"user_seconds": child_usage.ru_utime, "system_seconds": child_usage.ru_stime,
               "max_process_rss_kb": child_usage.ru_maxrss / (1024 if sys.platform == "darwin" else 1),
               "major_faults": child_usage.ru_majflt, "minor_faults": child_usage.ru_minflt,
               "fs_inputs": child_usage.ru_inblock, "fs_outputs": child_usage.ru_oublock}
    atomic_json(usage, metrics)
    metrics.update(wall_seconds=wall, reserved_core_hours=wall * args.threads / 3600,
                   cpu_hours=(metrics.get("user_seconds", 0) + metrics.get("system_seconds", 0)) / 3600,
                   host=socket.gethostname(), slots=args.threads, exit_code=exit_code, **peaks,
                   requested_memory_per_slot_gb=os.environ.get("CSM_BENCH_MEMORY_PER_SLOT_GB", "8"),
                   requested_runtime_hours=os.environ.get("CSM_BENCH_RUNTIME_HOURS", "12"),
                   job_id=os.environ.get("JOB_ID"), task_id=os.environ.get("SGE_TASK_ID"),
                   cpu_affinity=sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else [],
                   command=command, started_utc=timestamp)
    atomic_json(logs / f"{args.library}.{name}.metrics.json", metrics)
    if exit_code:
        raise subprocess.CalledProcessError(exit_code, command)
    payload = {"complete": True, "stage": name, "library": args.library,
               "code_sha256": read_json(args.work_dir / "experiment.json")["code_sha256"] if (args.work_dir / "experiment.json").exists() else None,
               "dependencies": parents, "metrics": metrics,
               "artifacts": [fingerprint(p) for p in outputs]}
    atomic_json(path, payload)
    started_path.unlink()
    return payload


def tree_rss(pid):
    """Sample Linux descendants. Sum RSS includes shared pages more than once."""
    try:
        status = Path(f"/proc/{pid}/status").read_text()
        rss = next((int(line.split()[1]) * 1024 for line in status.splitlines() if line.startswith("VmRSS:")), 0)
        children = Path(f"/proc/{pid}/task/{pid}/children").read_text().split()
        return rss + sum(tree_rss(int(child)) for child in children)
    except (OSError, ValueError):
        return 0


def tree_bytes(paths):
    total, seen = 0, set()
    for root in paths:
        entries = [(str(Path(root).parent), [], [Path(root).name])] if Path(root).is_file() else os.walk(root)
        for directory, _, names in entries:
            for name in names:
                try:
                    stat = (Path(directory) / name).stat()
                    identity = (stat.st_dev, stat.st_ino)
                    if identity not in seen:
                        total += stat.st_size
                        seen.add(identity)
                except FileNotFoundError:
                    pass
    return total


def internal_command(args, action, arm=None):
    command = [sys.executable, str(SCRIPT), "internal", "--action", action,
               "--config", str(args.config), "--work-dir", str(args.work_dir),
               "--library", args.library, "--reference", args.reference,
               "--kallisto", args.kallisto, "--pairs", str(args.pairs),
               "--threads", str(args.threads), "--max-multimaps", str(args.max_multimaps)]
    if arm:
        command += ["--arm", arm]
    return command


def parse_command(args, mode, root, extra=()):
    cfg = config_for(args)
    return [sys.executable, str(SCRIPTS / "splitpipe_gex_wrapper.py"),
            "--wrapper_splitpipe_dir", cfg["paths"]["splitpipe_dir"],
            "--wrapper_patch", str(PATCH), "--mode", mode, "--one_step",
            "--output_dir", str(root), "--nthreads", str(args.threads), "--kit_score_skip", *extra]


def clone_process(source, destination):
    """Copy mutable Parse bookkeeping; link large immutable intermediates only."""
    destination.mkdir(parents=True, exist_ok=True)
    for path in source.iterdir():
        if not path.is_file():
            continue
        target = destination / path.name
        if target.exists():
            continue
        if path.suffix == ".bam" or path.name == "barcode_head.fastq":
            os.link(path, target)
        else:
            shutil.copy2(path, target)


def roots(args):
    lib = f"gex_{args.library}"
    return args.work_dir / "preprocess" / lib, args.work_dir / "alignment/star" / lib


def prepare(args):
    cfg = config_for(args)
    library = next(lib for lib in cfg["libraries"] if lib["sublibrary_name"] == args.library)
    sampled = args.work_dir / "sampled_fastq"
    sampled.mkdir(exist_ok=True)
    fqs = [sampled / f"{library['gex_prefix']}_R{r}_001.fastq.gz" for r in (1, 2)]
    stage(args, "sample", internal_command(args, "sample"), fqs + [sampled / f"{args.library}.json"])
    pre, _ = roots(args)
    extras = ["--kit", cfg["splitpipe"]["kit"], "--chemistry", cfg["splitpipe"]["chemistry"],
              "--genome_dir", cfg["paths"]["genome_dir"], "--fq1", str(fqs[0]), "--fq2", str(fqs[1]),
              "--sample", library["sample_name"], cfg["splitpipe"]["sample_wells"]]
    env = {k: v for k, v in os.environ.items() if not k.startswith("CSM_GEX_DECONTAM")}
    if cfg.get("gex_decontamination", {}).get("enabled", False):
        filter_dir = args.work_dir / "shared_filter" / args.library
        filter_dir.mkdir(parents=True, exist_ok=True)
        identity = filter_dir / "step1a_input_identity.json"
        subprocess.run([sys.executable, str(SCRIPTS / "write_fastq_input_identity.py"), "--library", args.library,
                        "--output", str(identity), "--r1", str(fqs[0]), "--r2", str(fqs[1])], check=True)
        env.update(CSM_GEX_DECONTAM_CONFIG=str(cfg_path(args)), CSM_GEX_DECONTAM_LIBRARY=args.library,
                   CSM_GEX_DECONTAM_OUTPUT_DIR=str(filter_dir), CSM_GEX_DECONTAM_INPUT_MANIFEST=str(identity),
                   CSM_GEX_DECONTAM_HOOK=str(SCRIPTS / "gex_grna_preprocess_hook.py"),
                   CSM_GEX_DECONTAM_RUN_ID=f"benchmark.v2.{args.library}")
    outputs = [pre / "process" / name for name in ("barcode_head.fastq", "barcode_data.csv", "run_proc_def.json")]
    stage(args, "preprocess", parse_command(args, "pre", pre, extras), outputs, ("sample",), env)


def align(args):
    pre, star = roots(args)
    require(args, "preprocess")
    if args.backend == "star":
        clone_process(pre / "process", star / "process")
        parameters = write_star_parameters(config_for(args), args.work_dir / "configs" / f"{args.library}.star.parameters.txt")
        stage(args, "star.align", parse_command(args, "align", star,
              ["--parfile", str(parameters)]),
              [star / "process/barcode_headAligned.out.bam"], ("preprocess",))
        stage(args, "star.sort", parse_command(args, "post", star),
              [star / "process/barcode_headAligned_sorted.bam"], ("star.align",))
    else:
        cache = args.work_dir / "cache" / args.library / "kallisto"
        cache.mkdir(parents=True, exist_ok=True)
        stage(args, "kallisto.format", internal_command(args, "format"),
              [cache / n for n in ("cdna.fastq", "tags.fastq", "format.pickle")], ("preprocess",))
        _, bc_len, umi_len, _ = load_pickle(cache / "format.pickle")
        command = [args.kallisto, "bus", "-i", str(Path(args.reference) / "index.idx"),
                   "-o", str(cache / "bus"), "-x", f"1,0,{bc_len}:1,{bc_len},{bc_len + umi_len}:0,0,0",
                   "-t", str(args.threads), "--fr-stranded", str(cache / "cdna.fastq"), str(cache / "tags.fastq")]
        stage(args, "kallisto.align", command,
              [cache / "bus" / n for n in ("output.bus", "transcripts.txt", "matrix.ec", "run_info.json")],
              ("kallisto.format",))


def retry_star(args):
    """Archive a failed pre-mapping attempt and restart only STAR from Parse cache."""
    if manifest(args, "star.align").exists():
        raise ValueError("STAR alignment already succeeded; never retry a valid alignment")
    started = manifest(args, "star.align").with_suffix(".started")
    metrics = args.work_dir / "logs" / f"{args.library}.star.align.metrics.json"
    if not started.exists() or not metrics.exists() or read_json(metrics)["exit_code"] == 0:
        raise ValueError("Retry requires a recorded failed STAR attempt and confirmation the job ended")
    archive = args.work_dir / "failed_attempts" / args.library / ("star.align." + str(time.time_ns()))
    archive.mkdir(parents=True)
    _, star = roots(args)
    if star.exists():
        star.rename(archive / "library_output")
    for suffix in ("log", "usage", "metrics.json"):
        path = args.work_dir / "logs" / f"{args.library}.star.align.{suffix}"
        if path.exists():
            path.rename(archive / path.name)
    started.rename(archive / started.name)
    atomic_json(archive / "repair.json", {"reason": args.repair_note, "library": args.library,
                "code_sha256": read_json(args.work_dir / "experiment.json")["code_sha256"]})
    args.backend = "star"
    align(args)


def assign(args):
    cache = args.work_dir / "cache" / args.library
    if args.backend == "star":
        (cache / "star_unique").mkdir(parents=True, exist_ok=True)
        (cache / "star_fractional").mkdir(exist_ok=True)
        stage(args, "star.unique.assign", internal_command(args, "annotate"),
              [cache / "star_unique/assigned.bam"], ("star.sort",))
        stage(args, "star.fractional.assign", internal_command(args, "evidence", "star_fractional_parse"),
              [cache / "star_fractional" / n for n in ("reads.sqlite", "evidence.pickle")], ("star.sort",))
    else:
        stage(args, "kallisto.assign", internal_command(args, "evidence", "kallisto_unique_parse"),
              [cache / "kallisto" / n for n in ("reads.sqlite", "evidence.pickle")], ("kallisto.align",))


def run_arm(args):
    b, a, _ = arm_options(args.arm)
    if getattr(args, "require_count_baseline", False):
        path = args.work_dir / "timing_baselines/shared_sqlite_v1/archive.json"
        record = read_json(path) if path.exists() else {}
        if (record.get("complete") is not True or {"arm": args.arm, "library": args.library} not in record.get("relocated_libraries", [])):
            raise ValueError("Retiming requires a successful, matching archive completion record")
    dep = f"star.{a}.assign" if b == "star" else "kallisto.assign"
    require(args, dep)
    cfg_path(args, args.arm)
    pre, star = roots(args)
    root = args.work_dir / args.arm / f"gex_{args.library}"
    clone_process((star if b == "star" else pre) / "process", root / "process")
    if b == "star" and a == "unique":
        bam = root / "process/barcode_headAligned_sorted.bam"
        cache_bam = args.work_dir / "cache" / args.library / "star_unique/assigned.bam"
        if not manifest(args, f"{args.arm}.molecules").exists():
            bam.unlink(missing_ok=True)
            os.link(cache_bam, bam)
        outputs = [root / "process" / n for n in ("tscp_assignment.csv.gz", "aw_uf_count_matrix.mtx",
                                                               "aw_uf_cell_metadata.csv", "aw_uf_bc_counts.csv")]
        stage(args, f"{args.arm}.molecules", internal_command(args, "native", args.arm), outputs, (dep,))
        command = parse_command(args, "dge", root, ["--save_anndata"])
    else:
        stage(args, f"{args.arm}.molecules", internal_command(args, "count", args.arm),
              [root / "process" / n for n in ("matrices.pickle", "tscp_assignment.csv")], (dep,))
        verify_retimed_counts(args, root)
        command = internal_command(args, "publish", args.arm)
    sample = next(lib["sample_name"] for lib in config_for(args)["libraries"] if lib["sublibrary_name"] == args.library)
    outputs = [root / sample / "DGE_filtered" / n for n in ("count_matrix.mtx", "cell_metadata.csv", "all_genes.csv")]
    outputs += [root / "process" / n for n in ("aw_uf_count_matrix.mtx", "aw_uf_cell_metadata.csv", "all_genes.csv")]
    stage(args, f"{args.arm}.publish", command, outputs, (f"{args.arm}.molecules",))


def archive_kallisto_counting(args):
    """Preserve completed shared-storage measurements before consistent retiming.

    Submit only after the prior counting/publication jobs have terminated.
    Refuse partial stages so live outputs can never be moved by this operation.
    """
    if not args.repair_note:
        raise ValueError("Archiving timings requires a recorded reason")
    baseline = args.work_dir / "timing_baselines/shared_sqlite_v1"
    if baseline.exists():
        raise FileExistsError(baseline)
    entries = []
    for arm in getattr(args, "archive_arms", None) or ARMS:
        if not arm.startswith("kallisto_"):
            continue
        for library in LIBRARIES:
            for phase in ("molecules", "publish"):
                if (args.work_dir / "manifests" / f"{library}.{arm}.{phase}.started").exists():
                    raise ValueError(f"Partial/live stage must finish before retiming: {library}/{arm}/{phase}")
            source = args.work_dir / arm / f"gex_{library}"
            if source.exists():
                entries.append((arm, library, source))
    baseline.mkdir(parents=True)
    for arm, library, source in entries:
        target = baseline / arm / source.name
        target.parent.mkdir(parents=True, exist_ok=True)
        source.rename(target)
        for phase in ("molecules", "publish"):
            stem = f"{library}.{arm}.{phase}"
            for path in [args.work_dir / "manifests" / f"{stem}.json", *list((args.work_dir / "logs").glob(f"{stem}.*"))]:
                if path.exists():
                    parent = baseline / path.parent.name
                    parent.mkdir(exist_ok=True)
                    path.rename(parent / path.name)
    atomic_json(baseline / "archive.json", {"complete": True, "reason": args.repair_note,
                "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                "archive_code_sha256": read_json(args.work_dir / "experiment.json")["code_sha256"],
                "relocated_libraries": [{"arm": arm, "library": library} for arm, library, _ in entries],
                "relocation_policy": "original manifests retained unchanged; artifact paths map below this baseline directory"})


def verify_retimed_counts(args, root):
    """Require bit-exact matrices/metadata and byte-identical molecule assignments."""
    baseline = args.work_dir / "timing_baselines/shared_sqlite_v1"
    original = baseline / args.arm / f"gex_{args.library}" / "process/matrices.pickle"
    if not original.exists():
        return
    proof = root / "process/count_staging_equality.json"
    current_manifest = manifest(args, f"{args.arm}.molecules")
    identity = sha256_file(current_manifest)
    if proof.exists() and read_json(proof).get("molecules_manifest_sha256") == identity:
        return
    from scipy import sparse
    import pandas as pd
    started = time.perf_counter()
    old, new = load_pickle(original), load_pickle(root / "process/matrices.pickle")
    if len(old) != len(new):
        raise AssertionError("Retiming changed molecule products")
    for first, second in zip(old, new):
        if sparse.issparse(first):
            equal = sparse.issparse(second) and first.shape == second.shape and first.dtype == second.dtype and (first != second).nnz == 0
        elif isinstance(first, pd.DataFrame):
            equal = isinstance(second, pd.DataFrame) and first.equals(second)
        else:
            equal = first == second
        if not equal:
            raise AssertionError(f"Retiming changed scientific counts or metadata: {args.library}/{args.arm}")
    old_record = read_json(baseline / "manifests" / current_manifest.name)
    new_record = read_json(current_manifest)
    old_hash = next(a["sha256"] for a in old_record["artifacts"] if a["path"].endswith("/tscp_assignment.csv"))
    new_hash = next(a["sha256"] for a in new_record["artifacts"] if a["path"].endswith("/tscp_assignment.csv"))
    if old_hash != new_hash:
        raise AssertionError("Retiming changed molecule assignment bytes")
    atomic_json(proof, {"complete": True, "library": args.library, "arm": args.arm,
                "exact_matrices_and_metadata": True, "identical_molecule_assignments": True,
                "molecules_manifest_sha256": identity, "baseline": str(original),
                "verification_seconds": time.perf_counter() - started})


def with_parse(args, callback):
    """Initialize vendor reference/sample objects without production completion sealing."""
    cfg = config_for(args, args.arm)
    splitpipe = Path(cfg["paths"]["splitpipe_dir"])
    sys.path.insert(0, str(splitpipe))
    from splitpipe import dge
    original = dge.run_dge
    called = False
    def invoke(spipe):
        nonlocal called
        called = True
        callback(spipe)
        return True
    dge.run_dge = invoke
    pre, _ = roots(args)
    root = args.work_dir / (args.arm or "initialize") / f"gex_{args.library}"
    clone_process(pre / "process", root / "process")
    old = sys.argv
    sys.argv = [str(splitpipe / "split-pipe"), "--mode", "dge", "--one_step", "--output_dir", str(root),
                "--nthreads", str(args.threads), "--kit_score_skip"]
    try:
        try:
            runpy.run_path(str(splitpipe / "split-pipe"), run_name="__main__")
        except SystemExit as exc:
            if exc.code not in (None, 0):
                raise
        if not called:
            raise RuntimeError("Parse did not initialize benchmark callback")
    finally:
        dge.run_dge, sys.argv = original, old


def cached_gene(read, *_):
    if not read.has_tag("GX") or not read.has_tag("RE"):
        raise ValueError("Native STAR cache is missing GX/RE annotations")
    return read.get_tag("GX"), float(read.get_tag("RE"))


def annotate_unique(source_path, target_path, info, min_mapq, get_gene):
    """Cache native Parse gene decisions once, preserving primary read order."""
    import pysam
    keys = ("gene_bins", "gtf_dict_stepsize", "gene_starts", "gene_ends", "gene_id_to_name", "gene_id_to_biotype", "genes_merged_exons")
    with pysam.AlignmentFile(source_path, "rb") as source, \
         pysam.AlignmentFile(target_path, "wb", template=source) as target:
        chrom = dict(enumerate(source.references))
        for read in source:
            if read.is_secondary:
                continue
            gene, frac = (get_gene(read, chrom, *(info[k] for k in keys))
                          if not read.is_unmapped and read.mapping_quality >= min_mapq else ("", 0))
            read.set_tag("GX", gene)
            read.set_tag("RE", int(frac > 0.5))
            target.write(read)


def internal(args):
    pre, star = roots(args)
    cache = args.work_dir / "cache" / args.library
    cfg = config_for(args, args.arm)
    if args.action == "sample":
        original = yaml.safe_load(Path(args.config).read_text())
        pairs = resolve_all_library_pairs(original, args.config)[args.library]
        quotas = [args.pairs // len(pairs) + (i < args.pairs % len(pairs)) for i in range(len(pairs))]
        library = next(lib for lib in original["libraries"] if lib["sublibrary_name"] == args.library)
        sampled = args.work_dir / "sampled_fastq"
        provenance = []
        with gzip.open(sampled / f"{library['gex_prefix']}_R1_001.fastq.gz", "wt", compresslevel=1) as out1, \
             gzip.open(sampled / f"{library['gex_prefix']}_R2_001.fastq.gz", "wt", compresslevel=1) as out2:
            for (r1, r2), quota in zip(pairs, quotas):
                readers = fastq_records(r1), fastq_records(r2)
                try:
                    for _ in range(quota):
                        first, second = next(readers[0]), next(readers[1])
                        if first[0].split()[0].removesuffix("/1") != second[0].split()[0].removesuffix("/2"):
                            raise ValueError("Sampled paired FASTQ names disagree")
                        for out, (h, s, q) in zip((out1, out2), (first, second)):
                            out.write(f"{h}\n{s}\n+\n{q}\n")
                except StopIteration as exc:
                    raise ValueError(f"Lane shorter than requested prefix ({quota} pairs): {r1}") from exc
                finally:
                    for reader in readers:
                        reader.close()
                provenance.append({"r1": str(r1), "r2": str(r2), "pairs": quota,
                                   "source_stat": [{"bytes": Path(p).stat().st_size, "mtime_ns": Path(p).stat().st_mtime_ns} for p in (r1, r2)]})
        atomic_json(sampled / f"{args.library}.json", {"pairs": args.pairs, "selection": "balanced_lane_prefix", "lanes": provenance})
    elif args.action == "format":
        try:
            from scripts.kallisto_umi import prepare_fastqs
        except ModuleNotFoundError:
            from kallisto_umi import prepare_fastqs
        folder = cache / "kallisto"
        dump_pickle(folder / "format.pickle", prepare_fastqs(pre / "process/barcode_head.fastq", folder / "cdna.fastq", folder / "tags.fastq"))
    elif args.action == "native":
        try:
            from scripts.splitpipe_gex_wrapper import install_runtime_overlay
        except ModuleNotFoundError:
            from splitpipe_gex_wrapper import install_runtime_overlay
        splitpipe = Path(cfg["paths"]["splitpipe_dir"])
        sys.path.insert(0, str(splitpipe))
        install_runtime_overlay(splitpipe, PATCH, arm_options(args.arm)[2])
        from splitpipe import molinfo
        molinfo.get_gene = cached_gene
        sys.argv = [str(splitpipe / "split-pipe"), "--mode", "mol", "--one_step", "--output_dir",
                    str(args.work_dir / args.arm / f"gex_{args.library}"), "--nthreads", str(args.threads), "--kit_score_skip"]
        runpy.run_path(str(splitpipe / "split-pipe"), run_name="__main__")
    else:
        def callback(spipe):
            import pandas as pd
            try:
                from scripts.gex_counting import count_evidence, publish_samples, read_features
                from scripts.kallisto_umi import read_equivalence_classes, spool_bus
                from scripts.run_star_fractional import spool_bam
            except ModuleNotFoundError:
                from gex_counting import count_evidence, publish_samples, read_features
                from kallisto_umi import read_equivalence_classes, spool_bus
                from run_star_fractional import spool_bam
            if args.action == "annotate":
                from splitpipe.molinfo import get_gene
                annotate_unique(star / "process/barcode_headAligned_sorted.bam", cache / "star_unique/assigned.bam",
                                spipe.get_gene_info(), spipe.get_par_val("mol_min_mapq_score", as_int=True), get_gene)
                return
            b, a, _ = arm_options(args.arm)
            folder = cache / ("kallisto" if b == "kallisto" else "star_fractional")
            if args.action == "evidence":
                if b == "star":
                    connection, classes, features, lookup, stats = spool_bam(
                        star / "process/barcode_headAligned_sorted.bam", folder / "reads.sqlite", spipe.get_gene_info(), max_multimaps=args.max_multimaps)
                    aligned = stats["gene_compatible_reads"]
                else:
                    lookup, bc_len, umi_len, stats = load_pickle(folder / "format.pickle")
                    features = read_features(Path(args.reference), folder / "bus/transcripts.txt")
                    classes = read_equivalence_classes(folder / "bus/matrix.ec", len(features))
                    connection, aligned = spool_bus(folder / "bus/output.bus", folder / "reads.sqlite", lookup, classes, bc_len, umi_len)
                connection.close()
                dump_pickle(folder / "evidence.pickle", (classes, features, lookup, aligned, stats))
            else:
                root = args.work_dir / args.arm / f"gex_{args.library}"
                if args.action == "count":
                    import tempfile
                    from splitpipe import cy_utils
                    from splitpipe.molinfo import collapse_polyN
                    classes, features, lookup, aligned, _ = load_pickle(folder / "evidence.pickle")
                    with tempfile.TemporaryDirectory(prefix="csm-bench-count-") as tmp:
                        database, execution = stage_read_database(folder / "reads.sqlite", tmp)
                        connection = sqlite3.connect(database.resolve().as_uri() + "?mode=ro", uri=True)
                        matrices = count_evidence(connection, classes, features, spipe.get_gene_info()["all_genes"], lookup,
                            database, root / "process/tscp_assignment.csv", backend_options(cfg), aligned_reads=aligned,
                            min_reads=spipe.get_par_val("mol_awuf_read_thresh", as_int=True), kernel=collapse_polyN,
                            mutation_function=cy_utils.single_mut_seqs, threads=args.threads, diagnostics=execution)
                    dump_pickle(root / "process/matrices.pickle", matrices)
                    atomic_json(root / "process/count_execution.json", execution)
                elif args.action == "publish":
                    matrices = load_pickle(root / "process/matrices.pickle")
                    publish_samples(spipe, root, matrices, backend_options(cfg), {"benchmark": "v2_cached_50m", **cfg["gex_alignment"]})
                else:
                    raise ValueError(args.action)
        with_parse(args, callback)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("init", "prepare", "align", "assign", "run", "report", "internal", "accept-revision", "retry-star", "archive-kallisto-counting"))
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--kallisto", required=True)
    parser.add_argument("--library", choices=LIBRARIES)
    parser.add_argument("--backend", choices=BACKENDS)
    parser.add_argument("--arm", choices=ARMS)
    parser.add_argument("--pairs", type=int, default=50_000_000)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--max-multimaps", type=int, default=3)
    parser.add_argument("--verify-cache", action="store_true")
    parser.add_argument("--action")
    parser.add_argument("--repair-note")
    parser.add_argument("--archive-arm", action="append", dest="archive_arms", choices=ARMS)
    parser.add_argument("--require-count-baseline", action="store_true")
    args = parser.parse_args()
    if args.pairs < 1 or args.threads < 1 or args.max_multimaps < 1:
        parser.error("pairs, threads and max-multimaps must be positive")
    args.config = args.config.resolve()
    args.work_dir = args.work_dir.resolve()
    args.reference = str(Path(args.reference).resolve())
    args.kallisto = shutil.which(args.kallisto) or str(Path(args.kallisto).resolve())
    args.work_dir.mkdir(parents=True, exist_ok=True)
    (args.work_dir / "manifests").mkdir(exist_ok=True)
    if args.mode not in ("init", "report", "accept-revision", "archive-kallisto-counting") and not args.library:
        parser.error("this mode requires --library")
    if args.mode in ("align", "assign") and not args.backend:
        parser.error("align/assign require --backend")
    if args.mode == "run" and not args.arm:
        parser.error("run requires --arm")
    contract(args)
    if args.mode in ("init", "accept-revision"):
        return
    if args.mode == "report":
        try:
            from scripts.benchmark_gex_v2_report import report
        except ModuleNotFoundError:
            from benchmark_gex_v2_report import report
        report(args)
    elif args.mode == "internal":
        internal(args)
    else:
        {"prepare": prepare, "align": align, "assign": assign, "run": run_arm, "retry-star": retry_star,
         "archive-kallisto-counting": archive_kallisto_counting}[args.mode](args)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Matched-depth pilot comparison of V1.5 STAR and native kallisto counting.

All data access belongs on compute nodes. Shared preprocessing makes the
alignment/counting comparison independent of barcode-calibration differences.
"""

from __future__ import annotations

import argparse
import copy
import datetime
import gzip
import json
import os
import resource
import shutil
import socket
import subprocess
import sys
import time
from pathlib import Path

import yaml

try:
    from scripts.gex_alignment import atomic_json, sha256_file
    from scripts.kallisto_umi import fastq_records
    from scripts.qc_fastq_resolver import resolve_all_library_pairs
except ModuleNotFoundError:
    from gex_alignment import atomic_json, sha256_file
    from kallisto_umi import fastq_records
    from qc_fastq_resolver import resolve_all_library_pairs


BASE_ARMS = ("star", "kallisto_unique", "kallisto_fractional")
UMI_FLAVORS = ("parse", "directional", "none")
ARMS = tuple(base if flavor == "parse" else f"{base}_{flavor}"
             for flavor in UMI_FLAVORS for base in BASE_ARMS)


def arm_options(arm):
    for flavor in UMI_FLAVORS[1:]:
        if arm.endswith("_" + flavor):
            return arm.removesuffix("_" + flavor), flavor
    return arm, "parse"


def timed_command(command, log_path, resource_path, env=None):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    with log_path.open("w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=env)
    elapsed = time.perf_counter() - start
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    scale = 1024**3 if sys.platform == "darwin" else 1024**2
    payload = {"elapsed_seconds": elapsed, "exit_code": result.returncode, "command": command,
               "peak_rss_gb": usage.ru_maxrss / scale,
               "host": socket.gethostname(),
               "cpu_affinity": sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else [],
               "child_user_seconds": usage.ru_utime, "child_system_seconds": usage.ru_stime}
    atomic_json(resource_path, payload)
    atomic_json(resource_path.with_suffix(".json"), payload)
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, command)
    return payload


def load_config(args):
    return yaml.safe_load(Path(args.config).read_text())


def runtime_config(args, arm):
    base, flavor = arm_options(arm)
    cfg = copy.deepcopy(load_config(args))
    cfg["paths"]["output_dir"] = str(args.work_dir / arm)
    cfg["paths"]["fastq_dir"] = str(args.work_dir / "sampled_fastq")
    cfg["paths"]["log_dir"] = str(args.work_dir / "logs")
    cfg["gex_decontamination"]["existing_alignment"]["enabled"] = False
    cfg["gex_decontamination"]["step1c"]["enabled"] = False
    cfg["gex_decontamination"]["combine_with_grna_arm"]["enabled"] = False
    if base == "star":
        cfg["gex_alignment"] = {"backend": "star"}
    else:
        cfg["gex_alignment"] = {"backend": "kallisto", "kallisto": {
            "reference_dir": str(args.work_dir / "reference"),
            "executable": str(args.work_dir / "tools/kallisto/kallisto"),
            "gene_assignment": "fractional" if base.endswith("fractional") else "unique",
            "output_level": "gene", "strand": "forward",
        }}
    # Preserve the identity of historical Parse-default benchmark outputs.
    if flavor != "parse":
        cfg["gex_alignment"]["umi_correction"] = flavor
    path = args.work_dir / "configs" / f"{arm}.yaml"
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    temporary.write_text(yaml.safe_dump(cfg, sort_keys=False))
    temporary.replace(path)
    return cfg, path


def preprocess(args):
    config = load_config(args)
    library = next(lib for lib in config["libraries"] if lib["sublibrary_name"] == args.library)
    pairs = resolve_all_library_pairs(config, args.config)[args.library]
    sampled = args.work_dir / "sampled_fastq"
    sampled.mkdir(parents=True, exist_ok=True)
    prefix = library["gex_prefix"]
    r1_out = sampled / f"{prefix}_R1_001.fastq.gz"
    r2_out = sampled / f"{prefix}_R2_001.fastq.gz"
    quotas = [args.pairs // len(pairs) + (index < args.pairs % len(pairs)) for index in range(len(pairs))]
    counts = []
    start_sampling = time.perf_counter()
    sampling_manifest = args.work_dir / "manifests" / f"{args.library}.sampling.json"
    reuse = False
    if sampling_manifest.is_file():
        previous = json.loads(sampling_manifest.read_text())
        reuse = previous.get("paired_reads") == args.pairs and all(path.is_file() for path in (r1_out, r2_out))
        if reuse:
            counts = previous["sampling"]
    if not reuse:
        with gzip.open(r1_out, "wt", compresslevel=1) as out1, gzip.open(r2_out, "wt", compresslevel=1) as out2:
            for (r1, r2), quota in zip(pairs, quotas):
                n = 0
                first, second = fastq_records(r1), fastq_records(r2)
                try:
                    while n < quota:
                        h1, s1, q1 = next(first)
                        h2, s2, q2 = next(second)
                        if h1.split()[0].removesuffix("/1") != h2.split()[0].removesuffix("/2"):
                            raise ValueError("R1/R2 read names disagree during matched sampling")
                        out1.write(f"{h1}\n{s1}\n+\n{q1}\n")
                        out2.write(f"{h2}\n{s2}\n+\n{q2}\n")
                        n += 1
                except StopIteration as exc:
                    raise ValueError(f"Source pair has fewer than its {quota} requested reads: {r1}") from exc
                finally:
                    first.close()
                    second.close()
                counts.append({"r1": r1, "r2": r2, "read_pairs": n, "selection": "deterministic_prefix_per_lane"})
    sampling_seconds = time.perf_counter() - start_sampling
    if reuse:
        sampling_seconds = previous["sampling_seconds"]
    if not reuse:
        atomic_json(sampling_manifest, {"paired_reads": args.pairs, "sampling": counts,
                                       "sampling_seconds": sampling_seconds})
    pre_root = args.work_dir / "preprocess" / f"gex_{args.library}"
    config_path = args.work_dir / "configs" / "preprocess.yaml"
    config_path.parent.mkdir(exist_ok=True)
    temporary = config_path.with_name(config_path.name + f".{os.getpid()}.tmp")
    temporary.write_text(yaml.safe_dump(config, sort_keys=False))
    temporary.replace(config_path)
    scripts = Path(__file__).resolve().parent
    splitpipe = Path(config["paths"]["splitpipe_dir"])
    command = [sys.executable, str(scripts / "splitpipe_gex_wrapper.py"),
               "--wrapper_splitpipe_dir", str(splitpipe),
               "--wrapper_patch", str(scripts.parent / "patches/parsebio_splitpipe_1.6.4_gex_decontamination.patch"),
               "--mode", "pre", "--one_step", "--kit", config["splitpipe"]["kit"],
               "--chemistry", config["splitpipe"]["chemistry"],
               "--genome_dir", config["paths"]["genome_dir"],
               "--fq1", str(r1_out), "--fq2", str(r2_out), "--output_dir", str(pre_root),
               "--sample", library["sample_name"], config["splitpipe"]["sample_wells"],
               "--nthreads", str(args.threads), "--kit_score_skip"]
    env = dict(os.environ)
    filter_dir = args.work_dir / "shared_filter" / args.library
    filter_dir.mkdir(parents=True, exist_ok=True)
    identity_path = filter_dir / "step1a_input_identity.json"
    subprocess.run([sys.executable, str(scripts / "write_fastq_input_identity.py"),
                    "--library", args.library, "--output", str(identity_path),
                    "--r1", str(r1_out), "--r2", str(r2_out)], check=True)
    env.update(CSM_GEX_DECONTAM_CONFIG=str(config_path), CSM_GEX_DECONTAM_LIBRARY=args.library,
               CSM_GEX_DECONTAM_OUTPUT_DIR=str(filter_dir),
               CSM_GEX_DECONTAM_HOOK=str(scripts / "gex_grna_preprocess_hook.py"),
               CSM_GEX_DECONTAM_RUN_ID=f"benchmark.{args.library}.{os.getpid()}",
               CSM_GEX_DECONTAM_INPUT_MANIFEST=str(identity_path), PYTHONDONTWRITEBYTECODE="1")
    timing = timed_command(command, args.work_dir / "logs" / f"{args.library}.preprocess.log",
                           args.work_dir / "logs" / f"{args.library}.preprocess.time", env)
    survivor = pre_root / "process/barcode_head.fastq"
    if not survivor.is_file():
        raise FileNotFoundError(survivor)
    for arm in ARMS:
        cfg, path = runtime_config(args, arm)
        destination = args.work_dir / arm / f"gex_{args.library}"
        if destination.exists():
            raise FileExistsError(f"Benchmark arm already exists: {destination}")
        def copy_pre(source, target):
            if Path(source).name == "barcode_head.fastq":
                os.link(source, target)
            else:
                shutil.copy2(source, target)
            return target
        shutil.copytree(pre_root, destination, copy_function=copy_pre)
        # Dedicated copies allow each arm's validation to retain source paths.
        arm_filter = args.work_dir / arm / "gex_grna_filter" / args.library
        shutil.copytree(filter_dir, arm_filter)
    manifest = {"complete": True, "library": args.library, "paired_reads": args.pairs,
                "sampling_seconds": sampling_seconds, "sampling": counts, "preprocessing": timing,
                "survivor_sha256": sha256_file(survivor),
                "sampled_fastq_sha256": [sha256_file(r1_out), sha256_file(r2_out)]}
    atomic_json(args.work_dir / "manifests" / f"{args.library}.preprocess.json", manifest)


def run_arm(args):
    # A failed rerun must not leave a successful timing record from an older run.
    completion = args.work_dir / "manifests" / f"{args.library}.{args.arm}.json"
    completion.unlink(missing_ok=True)
    config, config_path = runtime_config(args, args.arm)
    base, flavor = arm_options(args.arm)
    scripts = Path(__file__).resolve().parent
    root = args.work_dir / args.arm / f"gex_{args.library}"
    if not root.exists():
        # Extend a completed preparation without resampling or re-filtering.
        pre_root = args.work_dir / "preprocess" / f"gex_{args.library}"
        def copy_pre(source, target):
            if Path(source).name == "barcode_head.fastq":
                os.link(source, target)
            else:
                shutil.copy2(source, target)
            return target
        shutil.copytree(pre_root, root, copy_function=copy_pre)
        shutil.copytree(args.work_dir / "shared_filter" / args.library,
                        args.work_dir / args.arm / "gex_grna_filter" / args.library)
    source = root / "process/barcode_head.fastq"
    preparation = json.loads((args.work_dir / "manifests" / f"{args.library}.preprocess.json").read_text())
    if sha256_file(source) != preparation["survivor_sha256"]:
        raise ValueError("Benchmark arms do not share exactly the same preprocessed reads")
    if base == "star":
        command = [sys.executable, str(scripts / "splitpipe_gex_wrapper.py"),
                   "--wrapper_splitpipe_dir", config["paths"]["splitpipe_dir"],
                   "--wrapper_patch", str(scripts.parent / "patches/parsebio_splitpipe_1.6.4_gex_decontamination.patch"),
                   "--wrapper_umi_correction", flavor,
                   "--mode", "align", "--until_step", "dge", "--output_dir", str(root),
                   "--nthreads", str(args.threads), "--kit_score_skip"]
    else:
        command = [sys.executable, str(scripts / "run_kallisto_gex.py"), "--config", str(config_path),
                   "--library", args.library, "--threads", str(args.threads)]
    # Preprocessing is already complete. Its hook must not reopen old evidence.
    env = {key: value for key, value in os.environ.items() if not key.startswith("CSM_GEX_DECONTAM")}
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    log = args.work_dir / "logs" / f"{args.library}.{args.arm}.log"
    resource = args.work_dir / "logs" / f"{args.library}.{args.arm}.time"
    timing = timed_command(command, log, resource, env)
    if base == "star":
        try:
            from scripts.gex_alignment import write_star_policy
        except ModuleNotFoundError:
            from gex_alignment import write_star_policy
        write_star_policy(config, root)
    payload = dict(timing, complete=True, library=args.library, arm=args.arm,
                   backend=base, umi_correction=flavor,
                   shared_preprocess_seconds=preparation["preprocessing"]["elapsed_seconds"],
                   reads=args.pairs, survivor_sha256=preparation["survivor_sha256"], threads=args.threads)
    atomic_json(completion, payload)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("prepare", "run", "report"))
    parser.add_argument("--config", required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--library")
    parser.add_argument("--arm", choices=ARMS)
    parser.add_argument("--pairs", type=int, default=10_000_000)
    parser.add_argument("--threads", type=int, default=16)
    args = parser.parse_args()
    args.work_dir = args.work_dir.resolve()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    if args.mode == "prepare":
        preprocess(args)
    elif args.mode == "run":
        if args.arm is None:
            parser.error("run requires --arm")
        run_arm(args)
    else:
        try:
            from scripts.benchmark_gex_report import report
        except ModuleNotFoundError:
            from benchmark_gex_report import report
        report(args)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Full-library v2.1 ablations with shared, unchanged baseline evidence.

All read denominators use the complete Step-2 input. Cell comparisons use the
frozen GEX-good cohort. Each arm is counted by the production Step-4 code.
"""

import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd
import yaml
from grna_recovery import (
    RecoveryEngine,
    recovery_contract,
    settings,
    RECOVERY_COLUMNS,
    make_policy,
)


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def worker_config(cfg, ref):
    s = cfg["grna_structure"]
    m = cfg["grna_matching"]
    seq_names = defaultdict(list)
    perguide = defaultdict(list)
    for row in ref.drop_duplicates("grna_name").to_dict("records"):
        seq_names[row["grna_seq"]].append(row["grna_name"])
    unique = list(seq_names)
    ibars = list(ref.ibar_seq.unique())
    rows = ref.to_dict("records")
    for r in rows:
        perguide[r["grna_name"]].append(r["ibar_seq"])
    return dict(
        tso_seq=s["fixed_seq_1"],
        stem1=s["fixed_seq_2"],
        stem2=s["fixed_seq_3"],
        tracr1=s["fixed_seq_4"],
        tracr2=s["fixed_seq_5"],
        grna_lengths=s["grna_lengths"],
        ibar_length=s["ibar_length"],
        anchor_max_mm=s["fixed_seq_max_mismatches"],
        tso_max_dist=s["tso_max_edit_dist"],
        max_lev_grna=m["max_lev_dist"],
        max_lev_ibar=m["max_lev_dist_ibar"],
        grna_seqs=unique,
        grna_names=list(range(len(unique))),
        seq_id_to_names={i: seq_names[g] for i, g in enumerate(unique)},
        ibar_seqs=ibars,
        ibar_names=ibars,
        combined_lookup={
            (r["grna_seq"], r["ibar_seq"]): (
                r["grna_name"],
                r["target_gene"],
                r["pair_id"],
                r["target_combo"],
            )
            for r in rows
        },
        ref_grna_to_target={r["grna_name"]: r["target_gene"] for r in rows},
        grna_to_ibars=dict(perguide),
        grna_name_to_position={r["grna_name"]: int(r["position"]) for r in rows},
        single_guide_mode=False,
        save_failed_reads=False,
    )


STATE = {}


def init(root, cfg, ref, good, policies, options, contract):
    step3 = load(Path(root) / "scripts/03_match_grna_ibar.py", "bench_step3")
    wc = worker_config(cfg, ref)
    step3._init_worker(wc)
    engine = RecoveryEngine(contract, wc["tso_seq"], wc["anchor_max_mm"], options)
    STATE.update(
        step3=step3,
        engine=engine,
        good=good,
        policies=policies,
        wc=wc,
        contract=contract,
        options=options,
    )


def process(rows):
    st = STATE
    mapper = st["step3"]
    engine = st["engine"]
    good = st["good"]
    stats = Counter()
    profiles = Counter()
    baseline, _, base_stats = mapper._process_rows_baseline([r[:5] for r in rows])
    matched = {r[0] for r in baseline}
    baseline_good = [r for r in baseline if r[1] in good]
    deltas = []
    exclusions = set()
    # Production wrapper parity is tested in every chunk on a fixed 100-row
    # subset, without running all ablations through the baseline matcher again.
    for row in rows:
        if row[0] in matched:
            continue
        candidates = engine.candidates(row[4], row[5])
        end, _ = mapper.find_tso(row[4], engine.tso, st["wc"]["tso_max_dist"])
        for arm, policy in st["policies"].items():
            if not policy.get("enabled", True):
                continue
            result, reason = engine.resolve(candidates, policy, end)
            stats[arm + "|" + reason] += 1
            if result is None:
                continue
            profiles[
                arm,
                result["recovery_length_state"],
                result["recovery_ibar_state"],
                result["grna_name"],
            ] += 1
            if row[1] in good:
                deltas.append((arm, row[:4], result))
        if row[1] in good:
            for c in candidates:
                if c["loss"] == 0:
                    for r in engine.by_guide[c["position"], c["designed_spacer"]]:
                        exclusions.add((row[1], r["target_gene"]))
    sample = rows[:100]
    wc = dict(
        st["wc"],
        recovery_contract=st["contract"],
        recovery_options=st["options"],
        recovery_policy=make_policy(st["contract"], st["options"], st["policies"]["auto"]),
    )
    mapper._worker_state.update(
        recovery=engine, recovery_options=st["options"], recovery_policy=wc["recovery_policy"]
    )
    actual, _, _ = mapper._process_rows(sample)
    expected = {r[0] for r in baseline if r[0] in {r[0] for r in sample}}
    expected.update(
        r[0] for arm, r, _ in deltas if arm == "auto" and r[0] in {r[0] for r in sample}
    )
    # Restrict parity to retained cells: deltas intentionally omit other cells.
    actual_ids = {r[0] for r in actual if r[1] in good}
    expected = {r for r in expected if r in {x[0] for x in sample if x[1] in good}}
    assert actual_ids == expected, "Production wrapper and ablation engine disagree"
    mapper._worker_state["recovery"] = None
    return baseline_good, deltas, base_stats, stats, profiles, exclusions


def run(args):
    root = Path(__file__).resolve().parents[1]
    base = Path(args.output)
    out = base / "technical" / args.library
    out.mkdir(parents=True, exist_ok=False)
    cfg = yaml.safe_load((base / "configs/auto.yaml").read_text())
    ref = pd.read_csv(cfg["paths"]["grna_reference"]).fillna("")
    options = settings(cfg)
    contract = recovery_contract(ref, cfg)
    frozen = json.loads((base / f"pilot/grna_matches/{args.library}_recovery.json").read_text())
    common = dict(max_loss=3, min_retained=17, endpoint_policy="joint", ibar_corroboration="report")
    policies = {
        "endpoint": dict(common, max_loss=0),
        "spacer": dict(common, endpoint_policy="original"),
        "combined18": dict(common, min_retained=18),
        "combined17": common,
        "auto": frozen["effective"],
        "ibar_required": dict(common, ibar_corroboration="require"),
    }
    meta = pd.read_csv(
        root / f"output/grna_gex_invalid_union_20261023_v3/{args.library}/cell_metadata.csv.gz",
        keep_default_na=False,
    )
    good = set(meta.base_barcode)
    index = [l["sublibrary_name"] for l in cfg["libraries"]].index(args.library)
    sample = cfg["libraries"][index]["sample_name"]
    step3 = load(root / "scripts/03_match_grna_ibar.py", "bench_main_step3")
    step4 = load(root / "scripts/04_assign_grna_to_cells.py", "bench_step4")
    columns = step3.OUTPUT_COLS + RECOVERY_COLUMNS
    files = {arm: gzip.open(out / f"{arm}_delta.tsv.gz", "wt") for arm in policies}
    writers = {
        arm: csv.writer(handle, delimiter="\t", lineterminator="\n")
        for arm, handle in files.items()
    }
    for w in writers.values():
        w.writerow(columns)
    fbase = gzip.open(out / "baseline_matches.tsv.gz", "wt")
    wb = csv.writer(fbase, delimiter="\t", lineterminator="\n")
    wb.writerow(step3.OUTPUT_COLS)
    stats = Counter()
    base_stats = Counter()
    profiles = Counter()
    exclusions = set()
    started = time.monotonic()
    digest = 0
    baseline_rows = 0
    try:
        with ProcessPoolExecutor(
            max_workers=args.workers,
            initializer=init,
            initargs=(str(root), cfg, ref, good, policies, options, contract),
        ) as pool:
            chunks = step3._chunk_generator(
                str(root / f"output/grna_barcodes/{args.library}_barcodes.tsv.gz"), 100000, True
            )
            for bi, (baseline, deltas, bs, ss, pp, xx) in enumerate(
                step3._bounded_results(pool, process, chunks, 2 * args.workers), 1
            ):
                wb.writerows(baseline)
                baseline_rows += len(baseline)
                for r in baseline:
                    digest ^= int.from_bytes(
                        hashlib.blake2b(("\t".join(map(str, r))).encode(), digest_size=16).digest(),
                        "big",
                    )
                for arm, prefix, r in deltas:
                    r["recovery_policy_fingerprint"] = (
                        frozen["fingerprint"]
                        if arm == "auto"
                        else make_policy(contract, options, policies[arm])["fingerprint"]
                    )
                    writers[arm].writerow(
                        (
                            *prefix,
                            r["grna_name"],
                            r["grna_seq"],
                            r["ibar_seq"],
                            r["target_gene"],
                            0,
                            r["ibar_dist"],
                            "recovery_" + r["recovery_class"],
                            r["variant"],
                            r["grna_length"],
                            r["pair_id"],
                            r["target_combo"],
                            r["ibar_seq_corrected"],
                            r["ibar_decoupled"],
                            *(r[c] for c in RECOVERY_COLUMNS),
                        )
                    )
                base_stats.update(bs)
                stats.update(ss)
                profiles.update(pp)
                exclusions.update(xx)
                if bi % 10 == 0:
                    print(
                        json.dumps(
                            dict(
                                chunks=bi,
                                reads=base_stats["total_reads"],
                                seconds=round(time.monotonic() - started),
                            )
                        ),
                        flush=True,
                    )
    finally:
        fbase.close()
        for f in files.values():
            f.close()
    scan_seconds = time.monotonic() - started
    # Freeze baseline equality against the existing production read calls.
    old_digest = 0
    old_rows = 0
    with gzip.open(root / f"output/grna_matches/{args.library}_grna_matches.tsv.gz", "rt") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["bc_wells"] not in good:
                continue
            old_rows += 1
            old_digest ^= int.from_bytes(
                hashlib.blake2b(
                    ("\t".join(row[c] for c in step3.OUTPUT_COLS)).encode(), digest_size=16
                ).digest(),
                "big",
            )
    parity = dict(
        reconstructed_rows=baseline_rows,
        production_rows=old_rows,
        reconstructed_digest=str(digest),
        production_digest=str(old_digest),
        identical=baseline_rows == old_rows and digest == old_digest,
    )
    (out / "baseline_parity.json").write_text(json.dumps(parity, indent=2) + "\n")
    pd.DataFrame(
        [dict(bc_wells=f"{index}_{cell}", target_gene=gene) for cell, gene in sorted(exclusions)]
    ).to_csv(out / "additional_normal_exclusions.csv.gz", index=False)
    pd.DataFrame(
        [
            dict(arm=a, length_state=l, ibar_state=i, guide=g, reads=n)
            for (a, l, i, g), n in profiles.items()
        ]
    ).to_csv(out / "read_profiles.csv", index=False)
    (out / "read_summary.json").write_text(
        json.dumps(
            dict(
                library=args.library,
                baseline=dict(base_stats),
                recovery=dict(stats),
                policies=policies,
                scan_seconds=scan_seconds,
                retained_cells=len(good),
                baseline_parity=parity,
            ),
            indent=2,
        )
        + "\n"
    )
    group_cols = step4._determine_group_cols(str(out / "baseline_matches.tsv.gz"), True)
    if isinstance(group_cols, tuple):
        group_cols = group_cols[0]
    ga = cfg["guide_assignment"]
    available = [
        c
        for c in [
            "grna_name",
            "pair_id",
            "position",
            "target_combo",
            "target_version",
            "gi_pair",
            "gi_version",
        ]
        if c in ref
    ]
    assignments = {}
    summaries = []
    for arm in ["baseline"] + list(policies):
        paths = [("baseline", out / "baseline_matches.tsv.gz")]
        if arm != "baseline":
            paths.append(("rescue", out / f"{arm}_delta.tsv.gz"))
        counts, libstats = step4.stream_library_umi_dedup_sources(
            paths,
            args.library,
            sample,
            index,
            group_cols,
            ga.get("max_umi_edits", 1),
            ga.get("min_reads_per_umi", 1),
        )
        counts, calls = step4._annotate_and_assign(counts, ref, available, group_cols, ga)
        counts.to_csv(out / f"{arm}_umi_counts.csv.gz", index=False)
        calls.to_csv(out / f"{arm}_assignments.csv.gz", index=False)
        assignments[arm] = calls
        pairs = {
            (r.bc_wells, str(r.pair_id))
            for r in calls.itertuples()
            if pd.notna(r.pair_id) and str(r.pair_id)
        }
        cells = {c for c, _ in pairs}
        complete = {
            (r.bc_wells, str(r.pair_id))
            for r in calls.itertuples()
            if r.detection_status == "complete"
        }
        if arm == "baseline":
            base_pairs = pairs
            base_cells = cells
            base_complete = complete
        summaries.append(
            dict(
                arm=arm,
                library=args.library,
                retained_cells=len(good),
                reads=libstats["total_reads"],
                umis=int(counts.umi_count.sum()),
                assigned_cells=len(cells),
                pairs=len(pairs),
                complete_pairs=len(complete),
                new_cells=len(cells - base_cells),
                additional_pairs=len(pairs - base_pairs),
                lost_pairs=len(base_pairs - pairs),
                recovered_partners=len((complete - base_complete) & base_pairs),
                seconds=libstats["elapsed_s"],
            )
        )
    pd.DataFrame(summaries).to_csv(out / "technical_summary.csv", index=False)
    (out / "_SUCCESS").write_text("complete\n")
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--library", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--workers", type=int, default=4)
    run(p.parse_args())

#!/usr/bin/env python3
"""Test frozen policies on archived GEX guide evidence; preserve expression.

This archive was selected by the earlier structure scan. It is not a new scan
of all GEX reads, so yields describe the archive only.
"""

import argparse
import csv
import gzip
import json
import re
from collections import Counter
from pathlib import Path
import pandas as pd
import yaml
from grna_recovery import recovery_contract, settings, RecoveryEngine, RECOVERY_COLUMNS
from benchmark_grna_recovery import load, worker_config


def run(library):
    root = Path(__file__).resolve().parents[1]
    base = root / "output/grna_spacer_scan_20261023/fiveprime_20261104/v2_1"
    out = base / "gex_validation" / library
    out.mkdir(parents=True, exist_ok=True)
    cfg = yaml.safe_load((base / "configs/auto.yaml").read_text())
    ref = pd.read_csv(cfg["paths"]["grna_reference"]).fillna("")
    contract = recovery_contract(ref, cfg)
    options = settings(cfg)
    policy = json.loads((base / f"pilot/grna_matches/{library}_recovery.json").read_text())
    step3 = load(root / "scripts/03_match_grna_ibar.py", "gex_validation_step3")
    wc = worker_config(cfg, ref)
    wc.update(recovery_contract=contract, recovery_options=options, recovery_policy=policy)
    step3._init_worker(wc)
    source = root / f"output/gex_invalid_guides_20261023/full/{library}/gex_guide_evidence.tsv.gz"
    index = [l["sublibrary_name"] for l in cfg["libraries"]].index(library)
    good = set(
        pd.read_csv(
            root / f"output/grna_gex_invalid_union_20261023_v3/{library}/cell_metadata.csv.gz"
        ).base_barcode
    )
    stems = re.compile("|".join(re.escape(p["stem"]) for p in contract["profiles"].values()))
    genes = {r["grna_seq"]: r["target_gene"] for r in contract["rows"]}
    lengths = sorted({len(g) for g in genes})
    normal = set()
    stats = Counter()
    rows = []
    normal_engine = RecoveryEngine(
        contract, cfg["grna_structure"]["fixed_seq_1"], 2, dict(options, min_quality=0)
    )
    with gzip.open(source, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            stats["archived_reads"] += 1
            rows.append(
                tuple(row[k] for k in ["read_id", "bc_wells", "bc_bci", "umi", "r1_seq", "r1_qual"])
            )
            if row["bc_wells"] not in good:
                continue
            cell = f"{index}_{row['bc_wells']}"
            seq = row["r1_seq"]
            for match in stems.finditer(seq):
                for length in lengths:
                    tag = seq[max(0, match.start() - length) : match.start()]
                    if tag in genes:
                        normal.add((cell, genes[tag]))
            for candidate in normal_engine.candidates(seq, row["r1_qual"]):
                if candidate["loss"] == 0:
                    normal.add((cell, genes[candidate["designed_spacer"]]))
    matched, _, metrics = step3._process_rows(rows)
    columns = step3.OUTPUT_COLS + RECOVERY_COLUMNS
    with gzip.open(out / "matches.tsv.gz", "wt") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(matched)
    stats.update(metrics)
    stats["retained_cell_recovered_reads"] = sum(
        r[1] in good and r[17] != "baseline" for r in matched
    )
    pd.DataFrame(sorted(normal), columns=["bc_wells", "target_gene"]).to_csv(
        out / "normal_target_exclusions.csv.gz", index=False
    )
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                library=library,
                counts=stats,
                policy_fingerprint=policy["fingerprint"],
                scope="Previously archived barcode-valid, structurally eligible GEX guide reads; not all raw GEX",
                expression_changed=False,
            ),
            indent=2,
        )
        + "\n"
    )
    (out / "_SUCCESS").write_text("complete\n")
    print(json.dumps(dict(library=library, counts=stats)), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--library", required=True)
    run(p.parse_args().library)

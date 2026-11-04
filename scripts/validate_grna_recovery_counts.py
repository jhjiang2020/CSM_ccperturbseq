#!/usr/bin/env python3
"""Verify final Step-4 implementation against frozen full-pilot auto counts."""

import argparse
import json
from pathlib import Path
import pandas as pd
import yaml
from benchmark_grna_recovery import load


def run(library):
    root = Path(__file__).resolve().parents[1]
    base = root / "output/grna_spacer_scan_20261023/fiveprime_20261104/v2_1"
    out = base / "technical" / library
    assert (out / "_SUCCESS").exists()
    cfg = yaml.safe_load((base / "configs/auto.yaml").read_text())
    step4 = load(root / "scripts/04_assign_grna_to_cells.py", "final_step4")
    index = [lib["sublibrary_name"] for lib in cfg["libraries"]].index(library)
    groups = step4._determine_group_cols(str(out / "baseline_matches.tsv.gz"), True)
    ga = cfg["guide_assignment"]
    current, stats = step4.stream_library_umi_dedup_sources(
        [("baseline", out / "baseline_matches.tsv.gz"), ("rescue", out / "auto_delta.tsv.gz")],
        library,
        cfg["libraries"][index]["sample_name"],
        index,
        groups,
        ga.get("max_umi_edits", 1),
        ga.get("min_reads_per_umi", 1),
    )
    previous = pd.read_csv(out / "auto_umi_counts.csv.gz", keep_default_na=False)
    columns = groups + [
        c for c in current if c not in groups and ("count" in c or c.startswith("recovery_umi_"))
    ]
    a = current[columns].fillna("").sort_values(groups).reset_index(drop=True)
    b = previous[columns].fillna("").sort_values(groups).reset_index(drop=True)
    pd.testing.assert_frame_equal(a, b, check_dtype=False)
    result = dict(
        library=library,
        identical=True,
        groups=len(a),
        reads=int(a.read_count.sum()),
        umis=int(a.umi_count.sum()),
        shared_duplicates=stats["source_stats"]["rescue"].get("duplicate_shared_reads", 0),
    )
    (base / f"validation/{library}_final_count_parity.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(result, flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--library", required=True)
    run(p.parse_args().library)

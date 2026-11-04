#!/usr/bin/env python3
"""Compact technical transition tables from completed recovery benchmarks."""

import json
from pathlib import Path
import pandas as pd


def run():
    root = Path(__file__).resolve().parents[1]
    base = root / "output/grna_spacer_scan_20261023/fiveprime_20261104/v2_1"
    libraries = ["C1", "C2", "T1", "T2"]
    arms = ["endpoint", "spacer", "combined18", "combined17", "auto", "ibar_required"]
    summaries = []
    transitions = []
    thresholds = []
    validation = []
    readflow = []
    for lib in libraries:
        out = base / "technical" / lib
        assert (out / "_SUCCESS").exists(), lib
        library_summary = pd.read_csv(out / "technical_summary.csv").set_index("arm")
        summaries.append(library_summary)
        previous = pd.read_csv(out / "baseline_assignments.csv.gz", keep_default_na=False)
        before = previous.drop_duplicates("bc_wells").set_index("bc_wells").combo_type
        base_pairs = set(zip(previous.bc_wells, previous.pair_id))
        baseline_singlets = int(before.eq("singlet").sum())
        library_summary.loc[
            "baseline",
            ["singlet_cells", "new_singlets", "net_singlet_gain", "singlet_to_multiplet"],
        ] = [baseline_singlets, 0, 0, 0]
        read = json.loads((out / "read_summary.json").read_text())
        total = read["baseline"]["total_reads"]
        matched = sum(
            read["baseline"].get(k, 0) for k in ["exact_match", "fuzzy_match", "grna_only_match"]
        )
        for arm in arms:
            add = read["recovery"].get(arm + "|recovered", 0)
            readflow.append(
                dict(
                    library=lib,
                    arm=arm,
                    input_reads=total,
                    baseline_matched_reads=matched,
                    added_reads=add,
                    additional_percentage_points=100 * add / total,
                )
            )
            current = pd.read_csv(out / f"{arm}_assignments.csv.gz", keep_default_na=False)
            after = current.drop_duplicates("bc_wells").set_index("bc_wells").combo_type
            cells = sorted(set(before.index) | set(after.index))
            frame = pd.DataFrame(
                {
                    "before": before.reindex(cells).fillna("unassigned"),
                    "after": after.reindex(cells).fillna("unassigned"),
                }
            )

            singlets = int(after.eq("singlet").sum())
            library_summary.loc[
                arm, ["singlet_cells", "new_singlets", "net_singlet_gain", "singlet_to_multiplet"]
            ] = [
                singlets,
                int((frame.before.eq("unassigned") & frame.after.eq("singlet")).sum()),
                singlets - baseline_singlets,
                int((frame.before.eq("singlet") & frame.after.eq("multiplet")).sum()),
            ]
            for (a, b), n in frame.value_counts().items():
                transitions.append(dict(library=lib, arm=arm, before=a, after=b, cells=n))
            counts = pd.read_csv(out / f"{arm}_umi_counts.csv.gz", keep_default_na=False)
            pair_support = (
                counts[counts.pair_id.ne("")]
                .groupby(["bc_wells", "pair_id"])
                .umi_count_rescue_only.sum()
            )
            new_pairs = set(zip(current.bc_wells, current.pair_id)) - base_pairs
            for cutoff in [1, 2, 3]:
                supported = {key for key, n in pair_support.items() if n >= cutoff} & new_pairs
                thresholds.append(
                    dict(
                        library=lib,
                        arm=arm,
                        min_incremental_rescue_umis=cutoff,
                        new_pairs=len(supported),
                        new_cells=len({cell for cell, pair in supported} - set(before.index)),
                    )
                )
        p = json.loads((base / f"pilot/grna_matches/{lib}_recovery.json").read_text())
        chosen = {k: v for k, v in p["effective"].items() if k != "enabled"}
        for setting in p["preflight"]["settings"]:
            if setting["policy"] != chosen:
                continue
            for label, metric in setting["validation"].items():
                validation.append(dict(library=lib, challenge=label, **metric))
    keys = ["bc_wells", "grna_name", "target_gene", "ibar_seq_corrected", "ibar_seq_unresolved"]
    metrics = ["umi_count", "umi_count_raw", "read_count"]
    reconstructed = pd.concat(
        [
            pd.read_csv(
                base / f"technical/{lib}/baseline_umi_counts.csv.gz",
                usecols=keys + metrics,
                keep_default_na=False,
            )
            for lib in libraries
        ]
    )
    production = pd.read_csv(
        root / "output/grna_assignment/grna_umi_counts.csv.gz",
        usecols=keys + metrics,
        keep_default_na=False,
    )
    production = production[production.bc_wells.isin(set(reconstructed.bc_wells))]
    compared = reconstructed.merge(
        production, on=keys, how="outer", suffixes=("_reconstructed", "_production")
    ).fillna(0)
    mismatches = {
        column: int(compared[column + "_reconstructed"].ne(compared[column + "_production"]).sum())
        for column in metrics
    }
    (base / "baseline_count_parity.json").write_text(
        json.dumps(
            dict(
                groups=len(compared), mismatches=mismatches, identical=not any(mismatches.values())
            ),
            indent=2,
        )
        + "\n"
    )
    summary = pd.concat([s.reset_index() for s in summaries], ignore_index=True)
    summary.to_csv(base / "technical_summary.csv", index=False)
    numeric = [
        "retained_cells",
        "reads",
        "umis",
        "assigned_cells",
        "pairs",
        "complete_pairs",
        "new_cells",
        "additional_pairs",
        "lost_pairs",
        "recovered_partners",
        "singlet_cells",
        "new_singlets",
        "net_singlet_gain",
        "singlet_to_multiplet",
    ]
    summary.groupby("arm")[numeric].sum().reset_index().to_csv(
        base / "technical_combined.csv", index=False
    )
    pd.DataFrame(transitions).to_csv(base / "cell_transitions.csv", index=False)
    pd.DataFrame(thresholds).to_csv(base / "minimum_rescue_support.csv", index=False)
    pd.DataFrame(validation).to_csv(base / "automatic_validation.csv", index=False)
    pd.DataFrame(readflow).to_csv(base / "read_recovery.csv", index=False)
    (base / "_TECHNICAL_SUMMARY_SUCCESS").write_text("complete\n")
    print(summary.groupby("arm")[numeric].sum().to_string(), flush=True)


if __name__ == "__main__":
    run()

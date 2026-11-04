#!/usr/bin/env python3
"""Reference-wide error challenges independent of expression outcomes."""

import copy
import json
from collections import Counter
from pathlib import Path
import pandas as pd
import yaml
from grna_recovery import RecoveryEngine, recovery_contract, settings
from benchmark_grna_recovery import load, worker_config


def run():
    root = Path(__file__).resolve().parents[1]
    base = root / "output/grna_spacer_scan_20261023/fiveprime_20261104/v2_1"
    cfg = yaml.safe_load((base / "configs/auto.yaml").read_text())
    opts = settings(cfg)
    opts["min_retained"] = 17
    results = []
    for dataset, refpath in [
        ("pilot", root / "references/grna_ibar_reference.csv"),
        ("gavin", root / "output/gavin_vs_pilot_20261020/gavin.reference.snapshot.csv"),
    ]:
        ref = pd.read_csv(refpath).fillna("")
        design_cfg = copy.deepcopy(cfg)
        if dataset == "gavin":
            design_cfg["guide_assignment"]["dual_guide_counting_mode"] = False
            design_cfg["grna_structure"].update(single_guide_mode=True, ibar_length=10)
        contract = recovery_contract(ref, design_cfg)
        engine = RecoveryEngine(contract, cfg["grna_structure"]["fixed_seq_1"], 2, opts)
        tally = Counter()
        examples = []
        for row in ref.drop_duplicates(["position", "grna_seq", "ibar_seq"]).to_dict("records"):
            pos = int(row["position"])
            g = row["grna_seq"]
            profile = contract["profiles"][pos]
            prefix = "AAAAAAA" + engine.tso
            tail = profile["stem"] + row["ibar_seq"] + profile["tracr"]
            cases = []
            for loss in range(4):
                tag = g[loss:]
                cases.append((f"loss{loss}", prefix + tag + tail, True))
                if loss == 3:
                    for n in (1, 2, 3):
                        mutated = list(prefix)
                        for ix in (8, 16, 24)[:n]:
                            mutated[ix] = next(b for b in "ACGT" if b != mutated[ix])
                        cases.append((f"loss3_TSO_{n}sub", "".join(mutated) + tag + tail, True))
                    cases.extend(
                        [
                            (
                                "loss3_TSO_insertion",
                                prefix[:17] + "A" + prefix[17:] + tag + tail,
                                True,
                            ),
                            ("loss3_TSO_deletion", prefix[:17] + prefix[18:] + tag + tail, True),
                            (
                                "retained_spacer_error",
                                prefix + tag[:8] + "N" + tag[9:] + tail,
                                False,
                            ),
                            ("partial_TSO", prefix[-12:] + tag + tail, False),
                            ("truncated_downstream", prefix + tag + tail[:-10], False),
                        ]
                    )
            other = contract["profiles"].get(3 - pos)
            if other and other["tracr"] != profile["tracr"]:
                cases.append(
                    (
                        "cross_scaffold",
                        prefix + g[3:] + profile["stem"] + row["ibar_seq"] + other["tracr"],
                        False,
                    )
                )
            for label, seq, identity_expected in cases:
                pred, reason = engine.match(seq, "I" * len(seq))
                tally[label, "tested"] += 1
                if pred:
                    tally[label, "accepted"] += 1
                    correct = (pred["position"], pred["recovery_designed_spacer"]) == (pos, g)
                    tally[label, "wrong_identity"] += not correct
                    tally[label, "negative_accepted"] += not identity_expected
                    if not correct and len(examples) < 20:
                        examples.append(
                            dict(label=label, truth=g, predicted=pred["recovery_designed_spacer"])
                        )
                else:
                    tally[label, "rejected_" + reason] += 1
        for label in sorted({k[0] for k in tally}):
            results.append(
                dict(
                    dataset=dataset,
                    challenge=label,
                    **{k[1]: n for k, n in tally.items() if k[0] == label},
                )
            )
        (base / f"validation/{dataset}_ambiguity_examples.json").write_text(
            json.dumps(
                dict(suffix_collision_groups=engine.reference_safety, incorrect_examples=examples),
                indent=2,
            )
            + "\n"
        )
    pd.DataFrame(results).fillna(0).to_csv(
        base / "validation/reference_challenges.csv", index=False
    )
    # Compare the off implementation with the exact effective v2.0 snapshot on
    # all four independent whole-file pilot samples (400k reads).
    ref = pd.read_csv(cfg["paths"]["grna_reference"]).fillna("")
    wc = worker_config(cfg, ref)
    old = load(base / "baseline/scripts/03_match_grna_ibar.py", "snapshot_step3")
    new = load(root / "scripts/03_match_grna_ibar.py", "current_step3")
    old._init_worker(wc)
    new._init_worker(wc)
    parity = []
    for lib in ["C1", "C2", "T1", "T2"]:
        path = base / f"pilot/grna_matches/{lib}_recovery_sample.tsv.gz"
        rows = []
        for chunk in new._chunk_generator(str(path), 100000, True):
            rows.extend(chunk)
        previous = old._process_rows([r[:5] for r in rows])
        current = new._process_rows([r[:5] for r in rows])
        same = previous == current
        parity.append(
            dict(library=lib, reads=len(rows), identical=same, matched_reads=len(current[0]))
        )
        assert same, lib
    (base / "validation/off_snapshot_parity.json").write_text(json.dumps(parity, indent=2) + "\n")
    (base / "validation/_CHALLENGES_SUCCESS").write_text("complete\n")
    print(pd.DataFrame(results).fillna(0).to_string(index=False), flush=True)


if __name__ == "__main__":
    run()

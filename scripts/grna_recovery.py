"""Conservative, reference-aware 5-prime guide recovery (v2.1).

Identity and physical length are separate. The full reference is the competitor
set; iBARs never select the spacer. Existing successful calls are not rewritten.
"""

import copy
import csv
import gzip
import hashlib
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

try:
    import edlib
except ImportError:
    edlib = None

try:
    from scripts.guide_design import DesignMatcher, compile_design
except ModuleNotFoundError:
    from guide_design import DesignMatcher, compile_design

SCHEMA_VERSION = 1
RECOVERY_COLUMNS = [
    "recovery_class",
    "recovery_length_state",
    "recovery_designed_spacer",
    "recovery_possible_losses",
    "recovery_tso_endpoints",
    "recovery_ibar_state",
    "recovery_policy_fingerprint",
]
BASELINE_METADATA = ("baseline", "not_reassessed", "", "[]", "[]", "not_reassessed")


def fingerprint(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def settings(config):
    structure = config.get("grna_structure", {}) or {}
    shared = "relaxed_structure" in structure
    raw = config.get("grna_recovery", {}) or {}
    if shared:
        value = structure["relaxed_structure"]
        if value not in (False, True, "auto") or isinstance(value, (int, float)) and not isinstance(value, bool):
            raise ValueError("grna_structure.relaxed_structure must be false, true or auto")
        if raw:
            raise ValueError("Use grna_structure.relaxed_structure; do not combine it with grna_recovery")
        raw = {
            "mode": "auto" if value == "auto" else "manual" if value else "off",
            "spacer": {"max_5prime_loss": structure.get("relaxed_max_5prime_loss", 3),
                       "min_retained": structure.get("relaxed_min_retained", 17),
                       "min_quality": structure.get("relaxed_min_quality", 30)},
            "ibar_corroboration": "require" if structure.get("relaxed_require_ibar", False) else "report",
        }
        if value == "auto" and "relaxed_min_retained" not in structure:
            raw["spacer"]["min_retained"] = "auto"
    if not isinstance(raw, dict):
        raise ValueError("grna_recovery must be a mapping")
    # PyYAML reads an unquoted YAML `off` as False.
    mode = raw.get("mode", "off")
    if mode is False:
        mode = "off"
    if mode not in ("off", "audit", "auto", "manual"):
        raise ValueError("Invalid grna_recovery.mode")
    tso = raw.get("tso", {}) or {}
    spacer = raw.get("spacer", {}) or {}
    pilot = raw.get("pilot", {}) or {}
    result = dict(
        mode=mode,
        endpoint_policy=tso.get("endpoint_policy", "joint"),
        tso_max_edits=tso.get(
            "max_edit_distance", config.get("grna_structure", {}).get("tso_max_edit_dist", 3)
        ),
        search_end=tso.get("search_end", 50),
        max_loss=spacer.get("max_5prime_loss", 3),
        min_retained=spacer.get("min_retained", "auto"),
        min_quality=spacer.get("min_quality", 30),
        ibar_corroboration=raw.get("ibar_corroboration", "report"),
        sample_reads=pilot.get("sample_reads", 100000),
        seed=pilot.get("seed", 21),
        min_validation_families=pilot.get("min_validation_families", 1000),
        min_accuracy=pilot.get("min_accuracy", 0.995),
        near_optimal_fraction=pilot.get("near_optimal_fraction", 0.99),
        resolved_contract=raw.get("resolved_contract"),
    )
    if result["endpoint_policy"] not in ("joint", "original"):
        raise ValueError("Invalid recovery endpoint_policy")
    if result["ibar_corroboration"] not in ("report", "require"):
        raise ValueError("Invalid ibar_corroboration")
    for key, lo, hi in [
        ("tso_max_edits", 0, 5),
        ("search_end", 20, 100),
        ("max_loss", 0, 3),
        ("min_quality", 0, 40),
        ("sample_reads", 1, 100000),
        ("min_validation_families", 1, 100000),
    ]:
        v = result[key]
        if isinstance(v, bool) or not isinstance(v, int) or not lo <= v <= hi:
            raise ValueError(f"Invalid recovery {key}")
    if result["min_retained"] != "auto" and (
        isinstance(result["min_retained"], bool)
        or not isinstance(result["min_retained"], int)
        or not 15 <= result["min_retained"] <= 30
    ):
        raise ValueError("Invalid recovery min_retained")
    if not 0.9 <= result["min_accuracy"] < 1 or not 0 < result["near_optimal_fraction"] <= 1:
        raise ValueError("Invalid recovery validation threshold")
    if tso.get("partial_tso", False):
        raise ValueError("Partial TSO rescue is not validated or implemented")
    if spacer.get("retained_edit_distance", 0) != 0:
        raise ValueError("v2.1 recovery requires exact retained spacer sequence")
    if pilot.get("split_by", "cell_umi") != "cell_umi":
        raise ValueError("Recovery validation must split by cell_umi")
    if shared:
        result["shared"] = True
        result["max_tso_shift"] = structure.get("relaxed_max_tso_shift", 3)
        if isinstance(result["max_tso_shift"], bool) or not isinstance(result["max_tso_shift"], int) or not 0 <= result["max_tso_shift"] <= 10:
            raise ValueError("relaxed_max_tso_shift must be an integer from 0 to 10")
        if not isinstance(structure.get("relaxed_require_ibar", False), bool):
            raise ValueError("relaxed_require_ibar must be true or false")
    return result


def recovery_contract(reference, config, effective_design=None):
    # Recovery needs a complete downstream anchor; old short-read parsing remains
    # available when this new route is unavailable.
    contract = copy.deepcopy(effective_design or compile_design(reference, config))
    contract["ibar_correction"] = {key: 0 for key in contract["ibar_correction"]}
    return contract


def source_stamp(path):
    p = Path(path)
    s = p.stat()
    return dict(path=str(p.resolve()), size=s.st_size, mtime_ns=s.st_mtime_ns)


def sample_rows(path, limit=100000, seed=21, cache_path=None):
    """Uniform bounded read reservoir over the complete existing Step-2 file.

    Sampling is read weighted. Validation subsequently uses one representative
    per cell/UMI family and never splits a family between discovery/validation.
    """
    stamp = source_stamp(path)
    if cache_path:
        cache_path = Path(cache_path)
        meta = cache_path.with_suffix(".json")
        if cache_path.exists() and meta.exists():
            saved = json.loads(meta.read_text())
            if (
                saved.get("source") == stamp
                and saved.get("limit") == limit
                and saved.get("seed") == seed
            ):
                with gzip.open(cache_path, "rt") as f:
                    return list(csv.DictReader(f, delimiter="\t")), saved
    rng = random.Random(seed)
    sample = []
    total = 0
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not {"read_id", "bc_wells", "umi", "r1_seq", "r1_qual"} <= set(reader.fieldnames or []):
            return [], dict(source=stamp, status="unavailable_missing_sequence_quality_or_UMI")
        for total, row in enumerate(reader, 1):
            if total <= limit:
                sample.append(row)
            else:
                j = rng.randrange(total)
                if j < limit:
                    sample[j] = row
    info = dict(
        source=stamp,
        limit=limit,
        seed=seed,
        reads_scanned=total,
        sampled_reads=len(sample),
        sampling="uniform_read_reservoir_whole_file",
    )
    if cache_path:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(cache_path, "wt") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=list(sample[0])
                if sample
                else ["read_id", "bc_wells", "bc_bci", "umi", "r1_seq", "r1_qual"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(sample)
        meta.write_text(json.dumps(info, indent=2) + "\n")
    return sample, info


class RecoveryEngine:
    def __init__(self, contract, tso, anchor_max_mismatches=2, options=None):
        if edlib is None:
            raise RuntimeError("gRNA recovery requires edlib")
        self.contract = contract
        self.tso = tso
        self.mm = anchor_max_mismatches
        self.options = options or dict(
            max_loss=3,
            min_retained=17,
            tso_max_edits=3,
            search_end=50,
            min_quality=30,
            endpoint_policy="joint",
            ibar_corroboration="report",
        )
        self.index = defaultdict(set)
        self.by_guide = defaultdict(list)
        self.ibar = defaultdict(set)
        for row in contract["rows"]:
            p, g = row["position"], row["grna_seq"]
            self.by_guide[p, g].append(row)
            if row["ibar_seq"]:
                self.ibar[p, row["ibar_seq"]].add(g)
            # Shadow competitors are retained beyond assignment's loss cap.
            for loss in range(min(6, len(g) - 1) + 1):
                if loss == 0 or len(g) - loss >= 15:
                    self.index[p, g[loss:]].add((g, loss))
        self.lengths = sorted({len(tag) for _, tag in self.index})
        self.resolver = DesignMatcher(contract)
        self.reference_safety = []
        for (p, tag), hits in self.index.items():
            if len({g for g, _ in hits}) > 1:
                self.reference_safety.append(
                    dict(position=p, observed_spacer=tag, candidates=sorted(g for g, _ in hits))
                )

    @staticmethod
    def _starts(sequence, pattern, lo, hi, mismatches):
        """Pigeonhole seeds enumerate all Hamming-distance-bounded locations."""
        if len(pattern) <= mismatches:
            return range(max(0, lo), min(hi, len(sequence) - len(pattern)) + 1)
        starts = set()
        parts = mismatches + 1
        for i in range(parts):
            left = i * len(pattern) // parts
            right = (i + 1) * len(pattern) // parts
            seed = pattern[left:right]
            at = sequence.find(seed, max(0, lo + left), min(len(sequence), hi + right))
            while at >= 0:
                start = at - left
                if lo <= start <= hi and start >= 0 and start + len(pattern) <= len(sequence):
                    starts.add(start)
                at = sequence.find(seed, at + 1, min(len(sequence), hi + right))
        return sorted(starts)

    def candidates(self, sequence, quality=None):
        sequence = sequence.upper()
        o = self.options
        if quality is None or len(quality) != len(sequence):
            return []
        lo = max(0, len(self.tso) - o["tso_max_edits"]) + min(self.lengths)
        hi = o["search_end"] + max(self.lengths)
        distances = {}
        found = []
        for p, profile in self.contract["profiles"].items():
            stem, tracr = profile["stem"], profile["tracr"]
            if not stem or not tracr:
                continue
            for ss in self._starts(sequence, stem, lo, hi, self.mm):
                observed = sequence[ss : ss + len(stem)]
                sm = sum(a != b for a, b in zip(stem, observed))
                if len(observed) != len(stem) or sm > self.mm:
                    continue
                se = ss + len(stem)
                for ilen in profile["ibar_lengths"]:
                    ts = se + ilen
                    obs = sequence[ts : ts + len(tracr)]
                    if len(obs) != len(tracr):
                        continue
                    tm = sum(a != b for a, b in zip(tracr, obs))
                    if sm + tm > self.mm:
                        continue
                    bar = sequence[se:ts]
                    for length in self.lengths:
                        end = ss - length
                        if end < 1 or end > o["search_end"]:
                            continue
                        tag = sequence[end:ss]
                        hits = self.index.get((p, tag))
                        if (
                            not hits
                            or min(map(ord, quality[end:ss]), default=0) - 33 < o["min_quality"]
                        ):
                            continue
                        if end not in distances:
                            distances[end] = edlib.align(
                                self.tso[::-1],
                                sequence[:end][::-1],
                                mode="SHW",
                                k=o["tso_max_edits"],
                                task="distance",
                            )["editDistance"]
                        td = distances[end]
                        if td < 0:
                            continue
                        for g, loss in hits:
                            state = (
                                "absent_by_design"
                                if ilen == 0
                                else "agree"
                                if self.ibar.get((p, bar)) == {g}
                                else "discordant"
                                if (p, bar) in self.ibar
                                else "unrecognized"
                            )
                            found.append(
                                dict(
                                    position=p,
                                    designed_spacer=g,
                                    loss=loss,
                                    retained_length=length,
                                    spacer=tag,
                                    ibar=bar,
                                    ibar_length=ilen,
                                    ibar_state=state,
                                    tso_end=end,
                                    tso_distance=td,
                                    stem_start=ss,
                                    anchor_score=sm + tm,
                                )
                            )
        return found

    def resolve(self, candidates, policy=None, original_endpoint=None):
        p = policy or self.options
        credible = [
            c
            for c in candidates
            if p.get("endpoint_policy", "joint") != "original" or c["tso_end"] == original_endpoint
        ]
        if not credible:
            return None, "no_supported_layout"
        identities = {(c["position"], c["designed_spacer"]) for c in credible}
        # Include shadow candidates in ambiguity detection before length policy.
        if len(identities) != 1:
            return None, "ambiguous_spacer_or_position"
        eligible = [
            c
            for c in credible
            if c["loss"] <= p["max_loss"]
            and (c["loss"] == 0 or c["retained_length"] >= p["min_retained"])
            and (p.get("max_tso_shift") is None or
                 original_endpoint is not None and original_endpoint >= 0 and
                 abs(c["tso_end"] - original_endpoint) <= p["max_tso_shift"])
        ]
        if not eligible:
            return None, "outside_length_policy"
        signatures = {(c["stem_start"], c["ibar"], c["ibar_length"]) for c in credible}
        if len(signatures) != 1:
            return None, "ambiguous_scaffold_boundary"
        if p.get("ibar_corroboration", "report") == "require" and any(
            c["ibar_state"] not in ("agree", "absent_by_design") for c in credible
        ):
            return None, "ibar_not_corroborated"
        # Length interpretation uses every credible layout, not only the chosen
        # policy. A stricter floor cannot turn uncertainty into physical damage.
        losses = sorted({c["loss"] for c in credible})
        length_state = (
            "intact" if losses == [0] else "shortened_only" if 0 not in losses else "uncertain"
        )
        chosen = min(
            eligible, key=lambda c: (c["anchor_score"], c["tso_distance"], c["tso_end"], c["ibar"])
        )
        pos = chosen["position"]
        profile = self.contract["profiles"][pos]
        observation = dict(
            grna_seq=chosen["spacer"],
            ibar_seq=chosen["ibar"],
            grna_length=chosen["retained_length"],
            ibar_length=chosen["ibar_length"],
            position=pos,
            position_from_length=0,
            position_from_scaffold=pos,
            scaffold_variant=profile["scaffold_variant"],
            tracr_variant=profile["scaffold_variant"],
            structure_status="concordant",
            variant=pos,
            anchor_score=chosen["anchor_score"],
        )
        resolved = self.resolver.resolve(observation, chosen["designed_spacer"])
        if resolved is None:
            return None, "unresolved_reference"
        resolved.update(
            recovery_class="endpoint" if 0 in losses else "short_spacer",
            recovery_length_state=length_state,
            recovery_designed_spacer=chosen["designed_spacer"],
            recovery_possible_losses=json.dumps(losses, separators=(",", ":")),
            recovery_tso_endpoints=json.dumps(
                sorted({c["tso_end"] for c in credible}), separators=(",", ":")
            ),
            recovery_ibar_state=chosen["ibar_state"],
        )
        return resolved, "recovered"

    def match(self, sequence, quality=None, policy=None, original_endpoint=None):
        return self.resolve(self.candidates(sequence, quality), policy, original_endpoint)


def make_policy(contract, options, effective=None):
    selection = effective or dict(
        max_loss=options["max_loss"],
        min_retained=18 if options["min_retained"] == "auto" else options["min_retained"],
        endpoint_policy=options["endpoint_policy"],
        ibar_corroboration=options["ibar_corroboration"],
    )
    if "max_tso_shift" in options and "max_tso_shift" not in selection:
        selection["max_tso_shift"] = options["max_tso_shift"]
    payload = dict(
        schema_version=SCHEMA_VERSION,
        reference_fingerprint=contract["fingerprint"],
        requested={k: v for k, v in options.items() if k not in ("resolved_contract",)},
        effective=selection,
    )
    payload["fingerprint"] = fingerprint(payload)
    return payload


def validate_saved_policy(saved, contract, options):
    body = {k: v for k, v in saved.items() if k not in ("fingerprint", "preflight")}
    if fingerprint(body) != saved.get("fingerprint"):
        raise ValueError("Recovery contract content fingerprint mismatch")
    if saved.get("reference_fingerprint") != contract["fingerprint"]:
        raise ValueError("Stale recovery reference contract")
    requested = {k: v for k, v in options.items() if k != "resolved_contract"}
    if saved.get("requested") != requested:
        raise ValueError("Recovery settings differ from frozen contract")
    return saved


def preflight(rows, contract, config, options, sampling=None, baseline_ids=None):
    """Choose a policy without expression outcomes; validate disjoint families."""
    from scipy.stats import beta

    base = dict(options, min_retained=15)
    engine = RecoveryEngine(
        contract,
        config["grna_structure"]["fixed_seq_1"],
        config["grna_structure"].get("fixed_seq_max_mismatches", 2),
        base,
    )
    truth_engine = RecoveryEngine(contract, engine.tso, 0, dict(base, tso_max_edits=0))
    families = {}
    for row in rows:
        family = (row.get("bc_wells", ""), row.get("umi", ""))
        if not all(family):
            continue
        old = families.get(family)
        if old is None or row["read_id"] < old["read_id"]:
            families[family] = row
    discovery = []
    validation = []
    truth = []
    for family, row in sorted(families.items()):
        fold = (
            int(
                hashlib.sha256(
                    (str(options["seed"]) + "|" + "|".join(family)).encode()
                ).hexdigest()[:8],
                16,
            )
            % 2
        )
        (validation if fold else discovery).append(row)
        if fold:
            cs = truth_engine.candidates(row["r1_seq"], row.get("r1_qual"))
            full = [c for c in cs if c["loss"] == 0]
            if len({(c["position"], c["designed_spacer"], c["tso_end"]) for c in full}) == 1:
                truth.append((row, full[0]))
    floors = [18, 17] if options["min_retained"] == "auto" else [options["min_retained"]]
    if options.get("shared") and options["mode"] == "auto" and options["min_retained"] != "auto":
        floors = sorted({18, 17, options["min_retained"]}, reverse=True)
        floors = [floor for floor in floors if floor >= options["min_retained"]]
    policies = [
        dict(
            max_loss=d,
            min_retained=f,
            endpoint_policy=options["endpoint_policy"],
            ibar_corroboration=options["ibar_corroboration"],
        )
        for f in floors
        for d in range(options["max_loss"] + 1)
    ]
    if "max_tso_shift" in options:
        policies = [dict(policy, max_tso_shift=shift) for policy in policies
                    for shift in range(options["max_tso_shift"] + 1)]
    discovery_candidates = []
    pattern_counts = Counter()
    for row in discovery:
        if baseline_ids is not None and row["read_id"] in baseline_ids:
            continue
        cs = engine.candidates(row["r1_seq"], row.get("r1_qual"))
        endpoint = edlib.align(
            engine.tso,
            row["r1_seq"][: options["search_end"]],
            mode="HW",
            k=options["tso_max_edits"],
            task="locations",
        )
        end = endpoint["locations"][0][1] + 1 if endpoint["editDistance"] >= 0 else -1
        discovery_candidates.append((cs, end))
        # One family contributes once to each plausible interpretation. These
        # counts describe ambiguity; they are not physical damage estimates.
        for pos, loss, shift in {(c["position"], c["loss"], c["tso_end"] - end) for c in cs}:
            pattern_counts[pos, loss, shift] += 1
    # The same molecule stays in one fold; each family counts once per challenge.
    challenges = []
    for row, c in truth:
        seq = row["r1_seq"]
        qual = row["r1_qual"]
        end = c["tso_end"]
        g = c["designed_spacer"]
        for loss in range(options["max_loss"] + 1):
            altered = seq[:end] + ("G" if loss == 0 else "") + seq[end + loss :]
            q = qual[:end] + ("I" if loss == 0 else "") + qual[end + loss :]
            alignment = edlib.align(
                engine.tso,
                altered[: options["search_end"]],
                mode="HW",
                k=options["tso_max_edits"],
                task="locations",
            )
            endpoint = alignment["locations"][0][1] + 1 if alignment["editDistance"] >= 0 else -1
            challenges.append(
                (
                    loss,
                    (c["position"], g),
                    c["ibar_length"],
                    engine.candidates(altered, q),
                    endpoint,
                )
            )
    negatives = []
    for row, c in truth:
        seq = row["r1_seq"]
        end = c["tso_end"]
        stop = end + len(c["designed_spacer"])
        rng = random.Random(fingerprint([options["seed"], row["bc_wells"], row["umi"]]))
        tag = list(seq[end:stop])
        rng.shuffle(tag)
        altered = seq[:end] + "".join(tag) + seq[stop:]
        alignment = edlib.align(
            engine.tso,
            altered[: options["search_end"]],
            mode="HW",
            k=options["tso_max_edits"],
            task="locations",
        )
        endpoint = alignment["locations"][0][1] + 1 if alignment["editDistance"] >= 0 else -1
        negatives.append((engine.candidates(altered, row["r1_qual"]), endpoint))
    layouts = sorted(
        {
            (pos, ilen)
            for pos, profile in contract["profiles"].items()
            for ilen in profile["ibar_lengths"]
        }
    )
    summaries = []
    for policy in policies:
        detected = sum(
            engine.resolve(cs, policy, end)[0] is not None for cs, end in discovery_candidates
        )
        per_class = defaultdict(Counter)
        for loss, identity, ibar_length, cs, end in challenges:
            if loss > policy["max_loss"]:
                continue
            pred, _ = engine.resolve(cs, policy, end)
            label = f"pos{identity[0]}_ibar{ibar_length}_" + (
                "endpoint" if loss == 0 else "loss_" + str(loss)
            )
            per_class[label]["tested"] += 1
            if pred is not None:
                per_class[label]["accepted"] += 1
                per_class[label]["correct"] += (
                    pred["position"],
                    pred["recovery_designed_spacer"],
                ) == identity
        passed = True
        metrics = {}
        for label in [
            f"pos{pos}_ibar{ilen}_" + part
            for pos, ilen in layouts
            for part in ["endpoint"] + ["loss_" + str(d) for d in range(1, policy["max_loss"] + 1)]
        ]:
            counts = per_class[label]
            n = counts["accepted"]
            correct = counts["correct"]
            lower = float(beta.ppf(0.05, correct, n - correct + 1)) if correct else 0.0
            good = n >= options["min_validation_families"] and lower >= options["min_accuracy"]
            passed &= good
            metrics[label] = dict(counts, accuracy_lower_95=lower, passed=good)
        negative_accepted = sum(
            engine.resolve(cs, policy, end)[0] is not None for cs, end in negatives
        )
        upper = (
            float(beta.ppf(0.95, negative_accepted + 1, len(negatives) - negative_accepted))
            if negative_accepted < len(negatives)
            else 1.0
        )
        passed = passed and bool(layouts) and upper <= 1 - options["min_accuracy"]
        summaries.append(
            dict(
                policy=policy,
                discovery_identified_families=detected,
                validation=metrics,
                negative_controls=dict(
                    tested=len(negatives), accepted=negative_accepted, acceptance_upper_95=upper
                ),
                passed=bool(passed),
            )
        )
    passing = [s for s in summaries if s["passed"]]
    if passing:
        best = max(s["discovery_identified_families"] for s in passing)
        close = [
            s
            for s in passing
            if s["discovery_identified_families"] >= best * options["near_optimal_fraction"]
        ]
        selected = min(
            close, key=lambda s: (s["policy"]["max_loss"], s["policy"].get("max_tso_shift", 0), -s["policy"]["min_retained"])
        )["policy"]
        selected = dict(selected, enabled=True)
    else:
        selected = dict(policies[0], enabled=False)
    report = dict(
        sampling=sampling or {},
        sampled_reads=len(rows),
        unique_families=len(families),
        discovery_families=len(discovery),
        validation_families=len(validation),
        truth_families=len(truth),
        settings=summaries,
        status="validated" if passing else "inconclusive",
        note="Synthetic identity challenge precision is not empirical real-read FDR; biological validation is separate.",
    )
    report["structural_patterns"] = [dict(position=p, fiveprime_loss=d, tso_endpoint_shift=t, families=n)
                                     for (p, d, t), n in pattern_counts.most_common()]
    report["suffix_collision_groups"] = engine.reference_safety
    policy = make_policy(contract, options, selected)
    policy["preflight"] = report
    return policy


def prepare_recovery(
    config,
    reference,
    worker_config,
    barcode_path,
    output,
    library,
    raw_fastq=False,
    baseline_processor=None,
):
    """Freeze one validated library policy and attach it to Step 3 workers."""
    import yaml

    options = settings(config)
    worker_config.pop("recovery_contract", None)
    if options["mode"] == "off":
        return None
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    # A missing complete scaffold makes only this route unavailable. Invalid
    # references/configuration still fail compile_design rather than guessing.
    structure = config["grna_structure"]
    single = structure.get("single_guide_mode", False)
    if not structure.get("fixed_seq_4") or (not single and not structure.get("fixed_seq_5")):
        (output / f"{library}_recovery.json").write_text(
            json.dumps({"status": "unavailable_missing_downstream_anchor"}, indent=2) + "\n"
        )
        return None
    contract = recovery_contract(reference, config, worker_config.get("design_contract"))
    if options["resolved_contract"]:
        with open(options["resolved_contract"]) as handle:
            saved = yaml.safe_load(handle)
        policy = validate_saved_policy(saved, contract, options)
    elif options["mode"] == "auto":
        cache = output / f"{library}_recovery_sample.tsv.gz"
        previous = Path(barcode_path).parent / f"{library}_recovery_sample.tsv.gz"
        if previous.exists() and previous.with_suffix(".json").exists():
            metadata = json.loads(previous.with_suffix(".json").read_text())
            if (
                metadata.get("source") == source_stamp(barcode_path)
                and metadata.get("limit") == options["sample_reads"]
                and metadata.get("seed") == options["seed"]
            ):
                cache = previous
        rows, sampling = (
            ([], {"status": "unavailable_raw_fastq_has_no_validated_cell_UMI"})
            if raw_fastq
            else sample_rows(barcode_path, options["sample_reads"], options["seed"], cache)
        )
        baseline_ids = None
        if rows and baseline_processor is not None:
            baseline_ids = {
                r[0]
                for r in baseline_processor(
                    [
                        tuple(
                            row.get(k, "")
                            for k in ("read_id", "bc_wells", "bc_bci", "umi", "r1_seq")
                        )
                        for row in rows
                    ]
                )[0]
            }
        policy = preflight(rows, contract, config, options, sampling, baseline_ids)
    else:
        policy = make_policy(contract, options)
        policy["preflight"] = {
            "status": "audit" if options["mode"] == "audit" else "manual_unvalidated"
        }
    (output / f"{library}_recovery.json").write_text(json.dumps(policy, indent=2) + "\n")
    worker_config.update(
        recovery_contract=contract, recovery_options=options, recovery_policy=policy
    )
    if options["mode"] == "audit":
        worker_config["save_failed_reads"] = True
    print(
        f"  Recovery [{library}]: {options['mode']}, {policy['preflight']['status']}, {policy['effective']}",
        flush=True,
    )
    return policy


class ReservoirWriter:
    """Write Step-2 output unchanged while retaining a bounded read reservoir."""

    def __init__(self, handle, limit, seed):
        self.handle = handle
        self.limit = limit
        self.seed = seed
        self.rng = random.Random(seed)
        self.header = None
        self.rows = []
        self.total = 0

    def write(self, text):
        written = self.handle.write(text)
        for line in text.splitlines():
            if self.header is None:
                self.header = line
                continue
            self.total += 1
            if self.total <= self.limit:
                self.rows.append(line)
            else:
                index = self.rng.randrange(self.total)
                if index < self.limit:
                    self.rows[index] = line
        return written

    def save(self, path, source, coverage=None):
        path = Path(path)
        with gzip.open(path, "wt") as handle:
            handle.write(self.header + "\n")
            for row in self.rows:
                handle.write(row + "\n")
        metadata = dict(
            source=source_stamp(source),
            limit=self.limit,
            seed=self.seed,
            reads_scanned=self.total,
            sampled_reads=len(self.rows),
            sampling="uniform_read_reservoir_during_step2",
            coverage=coverage or {},
        )
        path.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")

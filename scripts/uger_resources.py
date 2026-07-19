#!/usr/bin/env python3
"""Validated UGER resource resolution shared by production launchers."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import yaml


@dataclass(frozen=True)
class Step1aResources:
    slots: int = 16
    memory_per_slot_gb: int = 8
    runtime_hours: int = 48

    @property
    def h_rt(self) -> str:
        return f"{self.runtime_hours:02d}:00:00"

    @property
    def h_vmem(self) -> str:
        return f"{self.memory_per_slot_gb}G"


def _positive_int(value, label: str, maximum: int) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be a positive integer")
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be a positive integer") from exc
    if str(value).strip() != str(parsed):
        raise ValueError(f"{label} must be a positive integer")
    if parsed < 1 or parsed > maximum:
        raise ValueError(f"{label} must be between 1 and {maximum}")
    return parsed


def load_step1a_resources(config_path: str | Path) -> Step1aResources:
    path = Path(config_path)
    with path.open() as handle:
        config = yaml.safe_load(handle) or {}
    uger = config.get("uger", {}) or {}
    if not isinstance(uger, dict):
        raise ValueError("uger must be a mapping")
    raw = uger.get("step1a", {}) or {}
    if not isinstance(raw, dict):
        raise ValueError("uger.step1a must be a mapping")
    return Step1aResources(
        slots=_positive_int(raw.get("slots", 16), "uger.step1a.slots", 256),
        memory_per_slot_gb=_positive_int(
            raw.get("memory_per_slot_gb", 8),
            "uger.step1a.memory_per_slot_gb",
            4096,
        ),
        runtime_hours=_positive_int(
            raw.get("runtime_hours", 48),
            "uger.step1a.runtime_hours",
            744,
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--format", choices=("json", "tsv"), default="json")
    args = parser.parse_args()
    resources = load_step1a_resources(args.config)
    if args.format == "tsv":
        print(
            "\t".join(
                (
                    str(resources.slots),
                    str(resources.memory_per_slot_gb),
                    str(resources.runtime_hours),
                    resources.h_rt,
                )
            )
        )
        return
    payload = asdict(resources)
    payload.update({"h_rt": resources.h_rt, "h_vmem": resources.h_vmem})
    print(json.dumps(payload, sort_keys=True))


if __name__ == "__main__":
    main()

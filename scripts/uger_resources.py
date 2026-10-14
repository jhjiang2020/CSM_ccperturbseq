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


@dataclass(frozen=True)
class Step3Resources(Step1aResources):
    slots: int = 8
    memory_per_slot_gb: int = 4
    runtime_hours: int = 12


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


def _load_resources(config_path, step, resource_type):
    path = Path(config_path)
    with path.open() as handle:
        config = yaml.safe_load(handle) or {}
    uger = config.get("uger", {})
    if uger is None:
        uger = {}
    if not isinstance(uger, dict):
        raise ValueError("uger must be a mapping")
    raw = uger.get(step, {})
    if raw is None:
        raw = {}
    if not isinstance(raw, dict):
        raise ValueError(f"uger.{step} must be a mapping")
    defaults = resource_type()
    return resource_type(
        slots=_positive_int(raw.get("slots", defaults.slots), f"uger.{step}.slots", 256),
        memory_per_slot_gb=_positive_int(
            raw.get("memory_per_slot_gb", defaults.memory_per_slot_gb),
            f"uger.{step}.memory_per_slot_gb",
            4096,
        ),
        runtime_hours=_positive_int(
            raw.get("runtime_hours", defaults.runtime_hours),
            f"uger.{step}.runtime_hours",
            744,
        ),
    )


def load_step1a_resources(config_path: str | Path) -> Step1aResources:
    return _load_resources(config_path, 'step1a', Step1aResources)


def load_step3_resources(config_path: str | Path) -> Step3Resources:
    return _load_resources(config_path, 'step3', Step3Resources)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--step", choices=("step1a", "step3"), default="step1a")
    parser.add_argument("--format", choices=("json", "tsv"), default="json")
    args = parser.parse_args()
    loader = load_step1a_resources if args.step == 'step1a' else load_step3_resources
    resources = loader(args.config)
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

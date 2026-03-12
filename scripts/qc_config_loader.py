#!/usr/bin/env python3
"""QC Config Loader — merges QC configs on top of base configs.

Functions:
    load_qc_config(path)       — merges qc_config.yaml on top of config.yaml
    load_guideseq_config(path) — standalone guideseq config with optional sc_config cross-ref

Usage:
    from scripts.qc_config_loader import load_qc_config, load_guideseq_config
    cfg = load_qc_config('qc_config.yaml')
    gs_cfg = load_guideseq_config('qc_config_guideseq.yaml')
"""

import os
import yaml
import copy


def _deep_merge(base, override):
    """Recursively merge override into base (in-place).

    - Dicts are merged recursively
    - Lists and scalars in override replace base entirely
    """
    for key, value in override.items():
        if (key in base
                and isinstance(base[key], dict)
                and isinstance(value, dict)):
            _deep_merge(base[key], value)
        else:
            base[key] = value
    return base


def load_qc_config(qc_config_path):
    """Load QC config, merged on top of the base config.

    Args:
        qc_config_path: path to qc_config.yaml

    Returns:
        dict — merged configuration (base + QC overrides)
    """
    qc_config_path = os.path.abspath(qc_config_path)
    qc_dir = os.path.dirname(qc_config_path)

    with open(qc_config_path) as f:
        qc_cfg = yaml.safe_load(f) or {}

    # --- Load base config ---
    base_path = qc_cfg.pop('base_config', None)
    if base_path:
        if not os.path.isabs(base_path):
            base_path = os.path.join(qc_dir, base_path)
        with open(base_path) as f:
            base_cfg = yaml.safe_load(f) or {}
    else:
        base_cfg = {}

    # --- Merge: qc overrides base ---
    merged = _deep_merge(copy.deepcopy(base_cfg), qc_cfg)

    # Store resolved paths for convenience
    merged['_qc_config_path'] = qc_config_path
    merged['_base_config_path'] = base_path

    return merged


def load_guideseq_config(gs_config_path):
    """Load standalone guide-seq config with optional single-cell config cross-reference.

    The guideseq config is designed to run INDEPENDENTLY.  If an 'sc_config'
    key is present and the file exists, shared fields (paths, grna_reference)
    are inherited from the single-cell config via deep-merge.  If absent or
    the file does not exist, the config is returned as-is.

    Args:
        gs_config_path: path to qc_config_guideseq.yaml

    Returns:
        dict — guideseq configuration (optionally merged with sc base)
    """
    gs_config_path = os.path.abspath(gs_config_path)
    gs_dir = os.path.dirname(gs_config_path)

    with open(gs_config_path) as f:
        gs_cfg = yaml.safe_load(f) or {}

    # --- Optional single-cell config cross-reference ---
    sc_config_path = gs_cfg.pop('sc_config', None)
    sc_cfg = {}

    if sc_config_path:
        if not os.path.isabs(sc_config_path):
            sc_config_path = os.path.join(gs_dir, sc_config_path)
        if os.path.exists(sc_config_path):
            with open(sc_config_path) as f:
                sc_cfg = yaml.safe_load(f) or {}
            print(f"  [load_guideseq_config] Cross-referencing sc config: {sc_config_path}")
        else:
            print(f"  [load_guideseq_config] sc_config '{sc_config_path}' not found, running standalone.")

    # --- Merge: guideseq config overrides sc config ---
    # Only inherit path/reference fields from sc config, guideseq-specific keys take priority
    if sc_cfg:
        merged = _deep_merge(copy.deepcopy(sc_cfg), gs_cfg)
    else:
        merged = gs_cfg

    # Store resolved paths for convenience
    merged['_guideseq_config_path'] = gs_config_path
    merged['_sc_config_path'] = sc_config_path

    return merged


if __name__ == '__main__':
    import sys
    import json

    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <qc_config.yaml>")
        sys.exit(1)

    cfg = load_qc_config(sys.argv[1])
    print(yaml.dump(cfg, default_flow_style=False, sort_keys=False))


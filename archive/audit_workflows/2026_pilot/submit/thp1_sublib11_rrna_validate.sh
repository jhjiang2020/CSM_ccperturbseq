#!/usr/bin/env bash
# Archived project-specific worker retained for provenance.
#$ -N thp11_rrna_test
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_vmem=8G
#$ -l h_rt=01:00:00

set -euo pipefail

PROJECT_DIR="${CSM_PROJECT_DIR:-/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq}"
CONFIG="${THP_RRNA_CONFIG:-${PROJECT_DIR}/config/THP_1_sublib11_rrna_audit.yaml}"
AUDIT_ROOT="${THP_RRNA_AUDIT_ROOT:-${PROJECT_DIR}/output/THP_1_sublib11/rrna_barcode_diagnostic}"

source /home/unix/jjiahao/my-conda/etc/profile.d/conda.sh
conda activate splitpipe
export PYTHONDONTWRITEBYTECODE=1
cd "${PROJECT_DIR}"

python -m unittest \
    tests.test_barcode_rejection_qc \
    tests.test_external_rrna_audit

python - "${CONFIG}" "${AUDIT_ROOT}/code_validation.json" <<'PY'
import json
import sys
from pathlib import Path

import yaml

from scripts.barcode_contract import expand_well_range, load_runtime_kit_contract
from scripts.barcode_rejection_qc import setup_runtime

config_path = Path(sys.argv[1])
destination = Path(sys.argv[2])
cfg = yaml.safe_load(config_path.read_text())
bcinfo, _kits, kit, chemistry, amplicon = setup_runtime(cfg)
contract = load_runtime_kit_contract(
    bcinfo,
    library="THP_1_sublib11",
    kit=kit,
    chemistry=chemistry,
    amplicon_sequence=amplicon,
    sample_wells=expand_well_range("A1-H12"),
)
counts = {
    str(round_number): int(
        (contract.barcode_table["bc_round"] == round_number).sum()
    )
    for round_number in (1, 2, 3)
}
assert kit == "WT_mega"
assert chemistry == 3
assert len(contract.sample_wells) == 96
assert len(contract.allowed_bc1_winds) == 96
assert counts["1"] == 302
destination.parent.mkdir(parents=True, exist_ok=True)
destination.write_text(
    json.dumps(
        {
            "validated": True,
            "kit": kit,
            "chemistry": chemistry,
            "physical_sample_wells": len(contract.sample_wells),
            "allowed_bc1_winds": len(contract.allowed_bc1_winds),
            "round_sequence_counts": counts,
            "contract_sha256": contract.artifact_sha256,
        },
        indent=2,
        sort_keys=True,
    )
    + "\n"
)
PY

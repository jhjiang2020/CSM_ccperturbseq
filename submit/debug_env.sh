#!/bin/bash
# ==============================================================================
# Debug Script — Check environment compatibility on UGER interactive node
# ==============================================================================
# Run on an interactive UGER node:
#   qrsh -q interactive -l h_vmem=4G
#   bash submit/debug_env.sh
# ==============================================================================

set -uo pipefail

echo "=========================================="
echo "  Environment Debug — $(date)"
echo "  Host: $(hostname)"
echo "=========================================="

# --- 1. System info ---
echo ""
echo "=== 1. SYSTEM INFO ==="
echo "Kernel:  $(uname -r)"
echo "OS:      $(cat /etc/redhat-release 2>/dev/null || cat /etc/os-release 2>/dev/null | head -2 || echo 'unknown')"
echo "Arch:    $(uname -m)"

# glibc version (this is the critical one)
GLIBC_VER=$(ldd --version 2>&1 | head -1)
echo "glibc:   ${GLIBC_VER}"
echo ""
echo "Available GLIBC versions in /lib64/libm.so.6:"
strings /lib64/libm.so.6 2>/dev/null | grep "^GLIBC_" | sort -V | tail -10
echo ""
echo "Available GLIBC versions in /lib64/libc.so.6:"
strings /lib64/libc.so.6 2>/dev/null | grep "^GLIBC_" | sort -V | tail -5

# --- 2. Dotkit / modules ---
echo ""
echo "=== 2. DOTKIT / MODULES ==="
echo "Checking if 'use' command is available..."
if command -v use &> /dev/null; then
    echo "  'use' is available"
    echo "  Trying: use .gcc-13.2.0"
    use .gcc-13.2.0 2>&1 || echo "  (failed or not available)"
    echo "  GCC version after dotkit: $(gcc --version 2>&1 | head -1 || echo 'not found')"
    echo "  NOTE: This loads a newer GCC but does NOT change glibc."
else
    echo "  'use' command not found — dotkit not available in this shell"
    echo "  Try: source /broad/software/scripts/useuse"
fi

# --- 3. Conda environment ---
echo ""
echo "=== 3. CONDA ENVIRONMENT ==="
echo "Activating conda..."
eval "$($HOME/my-conda/condabin/conda shell.bash hook)" 2>&1
conda activate splitpipe 2>&1

echo "Conda env:  $(conda info --envs 2>/dev/null | grep '*' || echo 'unknown')"
echo "Python:     $(which python)"
echo "Python ver: $(python --version)"
echo "Pip:        $(which pip)"
echo ""

# --- 4. Check ALL Python package imports ---
echo "=== 4. PYTHON PACKAGE IMPORTS ==="
echo "Testing each package import individually..."
echo ""

python3 << 'PYEOF'
import sys
import importlib
import traceback

packages = [
    # --- Core (stdlib) ---
    ("sys",           "stdlib"),
    ("os",            "stdlib"),
    ("pathlib",       "stdlib"),
    ("argparse",      "stdlib"),
    ("gzip",          "stdlib"),
    ("multiprocessing", "stdlib"),
    ("collections",   "stdlib"),

    # --- Split-pipe dependencies (from setup.py) ---
    ("numpy",         "split-pipe dep"),
    ("pandas",        "split-pipe dep"),
    ("scipy",         "split-pipe dep"),
    ("matplotlib",    "split-pipe dep"),
    ("natsort",       "split-pipe dep"),
    ("h5py",          "split-pipe dep — KNOWN GLIBC ISSUE"),
    ("scanpy",        "split-pipe dep"),
    ("anndata",       "split-pipe dep"),
    ("leidenalg",     "split-pipe dep"),
    ("igraph",        "split-pipe dep (python-igraph)"),
    ("indexed_gzip",  "split-pipe dep"),
    ("jinja2",        "split-pipe dep"),
    ("psutil",        "split-pipe dep"),
    ("pysam",         "split-pipe dep"),
    ("numba",         "split-pipe dep (conda)"),
    ("Cython",        "split-pipe dep"),

    # --- Pipeline extra dependencies ---
    ("yaml",          "extra dep (pyyaml)"),
    ("edlib",         "extra dep (fast alignment)"),
]

results = {"ok": [], "fail": []}

for pkg_name, category in packages:
    try:
        mod = importlib.import_module(pkg_name)
        ver = getattr(mod, '__version__', getattr(mod, 'version', '?'))
        print(f"  ✓ {pkg_name:<20s} {str(ver):<15s}  [{category}]")
        results["ok"].append(pkg_name)
    except ImportError as e:
        print(f"  ✗ {pkg_name:<20s} FAILED           [{category}]")
        print(f"    Error: {e}")
        results["fail"].append((pkg_name, str(e)))
    except Exception as e:
        print(f"  ✗ {pkg_name:<20s} ERROR            [{category}]")
        print(f"    Error: {e}")
        results["fail"].append((pkg_name, str(e)))

print(f"\n  Summary: {len(results['ok'])} OK, {len(results['fail'])} FAILED")

if results["fail"]:
    print("\n  === FAILED PACKAGES ===")
    for pkg, err in results["fail"]:
        print(f"  {pkg}: {err}")

    # Check if any failures are glibc-related
    glibc_issues = [p for p, e in results["fail"] if "GLIBC" in e]
    if glibc_issues:
        print(f"\n  ⚠ GLIBC-related failures: {glibc_issues}")
        print("  These packages were compiled on a newer OS than this node.")
        print("  Fix options:")
        print("    1. conda install -c conda-forge " + " ".join(glibc_issues))
        print("    2. Request a newer OS node: qrsh -q interactive -l h_vmem=4G,os=RedHat8")
        print("    3. Rebuild the conda env on this node")
PYEOF

# --- 5. Check split-pipe itself ---
echo ""
echo "=== 5. SPLIT-PIPE ==="
SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SUBMIT_DIR}")"
DEBUG_CONFIG="${PROJECT_DIR}/config.yaml"
SPLITPIPE_DIR=$(python3 -c "import yaml; cfg=yaml.safe_load(open('${DEBUG_CONFIG}')); print(cfg['paths']['splitpipe_dir'])" 2>/dev/null || echo "")

if [[ -n "$SPLITPIPE_DIR" && -f "${SPLITPIPE_DIR}/split-pipe" ]]; then
    echo "split-pipe path: ${SPLITPIPE_DIR}/split-pipe"
    echo "Testing split-pipe --version..."
    "${SPLITPIPE_DIR}/split-pipe" --version 2>&1 || echo "  (split-pipe --version failed)"
else
    echo "  Could not locate split-pipe"
fi

# --- 6. Check split-pipe internal imports ---
echo ""
echo "=== 6. SPLIT-PIPE INTERNAL IMPORTS ==="
python3 << PYEOF2
import sys
sys.path.insert(0, "${SPLITPIPE_DIR}")

modules = [
    "splitpipe.bcutils",
    "splitpipe.pb_const",
    "splitpipe.bcinfo",
    "splitpipe.utils",
]

for mod_name in modules:
    try:
        __import__(mod_name)
        print(f"  ✓ {mod_name}")
    except Exception as e:
        print(f"  ✗ {mod_name}: {e}")

# Try cy_utils (Cython)
try:
    from splitpipe import cy_utils
    print(f"  ✓ splitpipe.cy_utils (Cython)")
except ImportError:
    print(f"  ⚠ splitpipe.cy_utils not available (pure Python fallback will be used)")
PYEOF2

# --- 7. Shared library diagnostics for h5py ---
echo ""
echo "=== 7. H5PY SHARED LIBRARY CHECK ==="
H5PY_DIR=$(python3 -c "import os; import h5py; print(os.path.dirname(h5py.__file__))" 2>/dev/null || echo "")
if [[ -n "$H5PY_DIR" ]]; then
    echo "h5py location: ${H5PY_DIR}"
    echo "Checking linked libraries..."
    find "${H5PY_DIR}" -name "*.so" -exec ldd {} \; 2>&1 | grep "not found" | sort -u
    if [[ $? -ne 0 ]]; then
        echo "  All shared libraries resolved OK"
    fi
else
    echo "  h5py import failed — checking .so files in site-packages..."
    SITE=$(python3 -c "import site; print(site.getsitepackages()[0])" 2>/dev/null)
    if [[ -n "$SITE" ]]; then
        echo "  Site-packages: ${SITE}"
        find "${SITE}/h5py" -name "*.so" -exec ldd {} \; 2>&1 | grep "not found" | sort -u || true
        echo ""
        echo "  Checking h5py.libs bundled libraries..."
        find "${SITE}/h5py.libs" -name "*.so*" -exec ldd {} \; 2>&1 | grep "not found" | sort -u || true
    fi
fi

echo ""
echo "=========================================="
echo "  Debug complete: $(date)"
echo "=========================================="

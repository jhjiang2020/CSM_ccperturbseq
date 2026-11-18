#!/bin/bash
# Scheduled worker, read-only base environment plus an isolated assignment venv.
set -euo pipefail
: "${ASSIGN_PREFIX:?}" "${CORE_DIR:?}"
[[ ! -e "$ASSIGN_PREFIX" ]] || { echo 'Refusing to modify existing assignment environment' >&2; exit 1; }
python -m venv --system-site-packages "$ASSIGN_PREFIX"
source "$ASSIGN_PREFIX/bin/activate"
python -m pip install --disable-pip-version-check 'pyro-ppl==1.9.1' 'pyro-api==0.1.2' 'dask==2024.12.1' 'distributed==2024.12.1'
python -m pip install --no-deps 'crispat @ https://github.com/velten-group/crispat/archive/a9372b9cf9b82b7409cc3be4d192056ae57a2673.tar.gz'
python - <<'PY'
import crispat, anndata, scanpy, pandas, numpy, scipy, torch, pyro, yaml
from importlib.metadata import version
assert version('crispat') == '0.9.11'
print('crispat',version('crispat'),'torch',torch.__version__,'pyro',pyro.__version__)
PY
python -m pip freeze > "$ASSIGN_PREFIX/requirements-resolved.txt"
python -c 'import sys; print(sys.version); print(sys.base_prefix)' > "$ASSIGN_PREFIX/base-python.txt"
printf '%s\n' 'crispat 0.9.11 a9372b9cf9b82b7409cc3be4d192056ae57a2673' > "$ASSIGN_PREFIX/SEALED"

#!/bin/bash
# Publish a new dated kallisto tools prefix; existing sealed prefixes are read only.
#$ -q broad
#$ -N csm_shared_kallisto
#$ -j y
#$ -pe smp 4
#$ -binding linear:4
#$ -l h_vmem=8G
#$ -l h_rt=03:00:00
#$ -l os=RedHat8
set -euo pipefail
umask 0027
: "${CSM_PROJECT_DIR:?Absolute project checkout required}"
: "${CSM_TOOLS_PREFIX:?New dated shared tools prefix required}"
: "${CSM_PUBLICATION_DIR:?Writable publication record directory required}"
CSM_PRODUCTION_PREFIX="${CSM_PRODUCTION_PREFIX:-/medpop/rgupta/Jiahao/conda/envs/splitpipe-share}"
CSM_CONDA_EXE="${CSM_CONDA_EXE:-${CSM_PRODUCTION_PREFIX}/condabin/conda}"
CSM_TOOLS_ALIAS="${CSM_TOOLS_ALIAS:-$(dirname "${CSM_TOOLS_PREFIX}")/csm_kallisto-share}"
[[ -x "${CSM_CONDA_EXE}" && -x "${CSM_PRODUCTION_PREFIX}/bin/python" ]]
[[ "${CSM_TOOLS_PREFIX}" != "${CSM_PRODUCTION_PREFIX}" ]] || { echo 'Tools must use a separate prefix' >&2; exit 1; }
if [[ "${CSM_VALIDATE_ONLY:-false}" == true ]]; then
    [[ -x "${CSM_TOOLS_PREFIX}/bin/kallisto" && ! -e "${CSM_TOOLS_PREFIX}/csm_shared_kallisto.complete.json" ]]
else
    [[ ! -e "${CSM_TOOLS_PREFIX}" ]] || { echo "Refusing to modify an existing prefix: ${CSM_TOOLS_PREFIX}" >&2; exit 1; }
fi
mkdir -p "${CSM_PUBLICATION_DIR}"
trap 'touch "${CSM_PUBLICATION_DIR}/FAILED"' ERR
export PYTHONDONTWRITEBYTECODE=1
export PATH="${CSM_PRODUCTION_PREFIX}/bin:${PATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR}/csm-shared-mpl"
export NUMBA_CACHE_DIR="${TMPDIR}/csm-shared-numba"
export XDG_CACHE_HOME="${TMPDIR}/csm-shared-cache"
export CONDA_PKGS_DIRS="${TMPDIR}/csm-shared-conda-pkgs"
export PIP_CACHE_DIR="${TMPDIR}/csm-shared-pip-cache"
echo "Allocated CPUs: ${NSLOTS:-4}; production Python: ${CSM_PRODUCTION_PREFIX}"
cd "${CSM_PROJECT_DIR}"
sha256sum "${CSM_PRODUCTION_PREFIX}/csm_shared_environment.complete.json" \
    "${CSM_PRODUCTION_PREFIX}/conda-explicit.txt" \
    "${CSM_PRODUCTION_PREFIX}/pip-freeze.txt" > "${CSM_PUBLICATION_DIR}/production-before.sha256"
if [[ "${CSM_VALIDATE_ONLY:-false}" != true ]]; then
    "${CSM_CONDA_EXE}" env create --yes --prefix "${CSM_TOOLS_PREFIX}" \
        --file "${CSM_PROJECT_DIR}/environment_kallisto.yml" --solver libmamba
fi
"${CSM_PRODUCTION_PREFIX}/bin/python" -m pip check > "${CSM_PUBLICATION_DIR}/production-pip-check.txt"
"${CSM_TOOLS_PREFIX}/bin/python" -m pip check > "${CSM_PUBLICATION_DIR}/tools-pip-check.txt"
export PYTHONPATH="${CSM_PROJECT_DIR}/tests:${CSM_PROJECT_DIR}"
export CSM_TEST_KALLISTO="${CSM_TOOLS_PREFIX}/bin/kallisto"
export CSM_TEST_BUSTOOLS="${CSM_TOOLS_PREFIX}/bin/bustools"
"${CSM_PRODUCTION_PREFIX}/bin/python" -m unittest \
    test_kallisto_backend test_kallisto_optimization test_splitpipe_gex_wrapper test_uger_resources -v \
    > "${CSM_PUBLICATION_DIR}/runtime-tests.log" 2>&1
"${CSM_CONDA_EXE}" list --prefix "${CSM_TOOLS_PREFIX}" --explicit \
    > "${CSM_TOOLS_PREFIX}/conda-explicit.txt"
"${CSM_TOOLS_PREFIX}/bin/python" -m pip freeze > "${CSM_TOOLS_PREFIX}/pip-freeze.txt"
cp "${CSM_PROJECT_DIR}/environment_kallisto.yml" "${CSM_TOOLS_PREFIX}/environment-source.yml"
sha256sum --check "${CSM_PUBLICATION_DIR}/production-before.sha256"
"${CSM_PRODUCTION_PREFIX}/bin/python" "${CSM_PROJECT_DIR}/scripts/publish_shared_kallisto_environment.py" \
    --prefix "${CSM_TOOLS_PREFIX}" --production-prefix "${CSM_PRODUCTION_PREFIX}" \
    --record-dir "${CSM_PUBLICATION_DIR}" --project-dir "${CSM_PROJECT_DIR}" --group rguptalab
# Touch only the new tools tree. Symlink targets are never followed recursively.
chgrp -hR rguptalab "${CSM_TOOLS_PREFIX}"
chmod -R u+rwX,g+rX,go-w,o-rwx "${CSM_TOOLS_PREFIX}"
"${CSM_PRODUCTION_PREFIX}/bin/python" "${CSM_PROJECT_DIR}/scripts/publish_shared_kallisto_environment.py" \
    --prefix "${CSM_TOOLS_PREFIX}" --production-prefix "${CSM_PRODUCTION_PREFIX}" \
    --record-dir "${CSM_PUBLICATION_DIR}" --project-dir "${CSM_PROJECT_DIR}" --group rguptalab --verify-and-seal \
    --stable-alias "${CSM_TOOLS_ALIAS}"
sha256sum --check "${CSM_PUBLICATION_DIR}/production-before.sha256"
if [[ -f "${CSM_PUBLICATION_DIR}/FAILED" ]]; then
    mv "${CSM_PUBLICATION_DIR}/FAILED" "${CSM_PUBLICATION_DIR}/FAILED.before-${JOB_ID}"
fi
touch "${CSM_PUBLICATION_DIR}/SUCCESS"

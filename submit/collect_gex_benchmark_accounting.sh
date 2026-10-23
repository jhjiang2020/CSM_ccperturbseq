#!/bin/bash
# Login-safe scheduler metadata only. Run after the benchmark jobs finish.
set -eo pipefail
: "${CSM_BENCH_DIR:?absolute benchmark directory}"
source /broad/software/scripts/useuse >/dev/null 2>&1 || true
reuse -q UGER >/dev/null 2>&1 || true
command -v qacct >/dev/null
record="${CSM_BENCH_DIR}/logs/uger_submissions.jsonl"
destination="${CSM_BENCH_DIR}/logs/accounting"
mkdir -p "${destination}"
ids=$({
    if [[ -s "${record}" ]]; then
        sed -n 's/.*"job_id": "\([0-9][0-9]*\)".*/\1/p' "${record}"
    fi
    for sidecar in "${record%.jsonl}.records/"*.json; do
        [[ -f "${sidecar}" ]] || continue
        sed -n 's/.*"job_id": "\([0-9][0-9]*\)".*/\1/p' "${sidecar}"
    done
} | sort -u)
test -n "${ids}"
for job in ${ids}; do
    temporary="${destination}/${job}.qacct.tmp"
    if timeout 20s qacct -j "${job}" > "${temporary}" 2> "${destination}/${job}.qacct.stderr"; then
        test -s "${temporary}"
        mv "${temporary}" "${destination}/${job}.qacct.txt"
    else
        rm -f "${temporary}"
        echo "Accounting not yet available for ${job}" >&2
    fi
done
date -u +%Y-%m-%dT%H:%M:%SZ > "${destination}/collected_utc.txt"

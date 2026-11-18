#!/usr/bin/env bash
# Fetch completed phase reruns and align the 20-step thermo/dump records.
set -euo pipefail
cd "$(dirname "$0")/.."
sock=${1:?active SSH control socket}
host=stu454@sylogin.hpc.sjtu.edu.cn
remote=$(cat remote-path.txt)
export OVITO_THREAD_COUNT=${OVITO_THREAD_COUNT:-4}

while read -r name delay job; do
    [[ "$name" =~ ^(flat|step)-0[12]$ && "$delay" =~ ^(750|1500)$ && "$job" =~ ^[0-9]+$ ]]
    case_name="phase-$name-$delay"
    remote_case="$remote/runs/$case_name"
    record=$(ssh -n -F /dev/null -S "$sock" -o BatchMode=yes "$host" \
        "sacct -X -n -P -j '$job' -o State,ExitCode")
    printf '%s %s %s\n' "$case_name" "$job" "$record"
    if [[ "$record" != 'COMPLETED|0:0' ]]; then
        continue
    fi
    if [[ -s "audit/$case_name/phase_audit.complete" ]]; then
        continue
    fi
    ssh -n -F /dev/null -S "$sock" -o BatchMode=yes "$host" "test -f '$remote_case/phase.complete'"
    case_dir="data/$case_name"
    mkdir -p "$case_dir"
    for file in equilibrated.dump hold.log loading.dump loading.log; do
        target="$case_dir/$file"
        if [[ -s "$target" ]]; then
            continue
        fi
        touch "$target.partial"
        offset=$(stat -c %s "$target.partial")
        ssh -n -F /dev/null -S "$sock" -o BatchMode=yes "$host" \
            "tail -c +$((offset + 1)) '$remote_case/$file'" >> "$target.partial"
        remote_hash=$(ssh -n -F /dev/null -S "$sock" -o BatchMode=yes "$host" \
            "sha256sum '$remote_case/$file'")
        local_hash=$(sha256sum "$target.partial")
        if [[ "${local_hash%% *}" != "${remote_hash%% *}" ]]; then
            printf '%s: check transfer before retrying\n' "$target.partial" >&2
            exit 1
        fi
        mv "$target.partial" "$target"
    done
    ../../.venv/bin/python -B analyze.py "$case_dir" --out audit
    ../../.venv/bin/python -B - "audit/$case_name/loading_audit.csv" <<'PY'
import csv
import sys
from pathlib import Path

with Path(sys.argv[1]).open(newline="") as stream:
    steps = [int(row["step"]) for row in csv.DictReader(stream)]
if steps != list(range(0, 18001, 20)):
    raise SystemExit(f"Unexpected dump sampling: {sys.argv[1]}")
PY
    printf "%s\n" "$job" > "audit/$case_name/phase_audit.complete"
done < phase/active-jobs.txt

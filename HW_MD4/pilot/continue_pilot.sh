#!/usr/bin/env bash
# Audit equilibrated structures, submit accepted cases, and retrieve completed trajectories.
set -euo pipefail
cd "$(dirname "$0")"
exec 9> .continue.lock
flock -n 9 || { echo "Pilot watcher already active" >&2; exit 1; }
sock=${1:?pass active SSH control socket path}
host=stu454@sylogin.hpc.sjtu.edu.cn
remote=$(cat remote-path.txt)
remote_cmd() {
    ssh -n -F /dev/null -S "$sock" -o BatchMode=yes -o ConnectTimeout=10 "$host" "cd '$remote' && $1"
}
fetch() {
    scp -q -F /dev/null -o ControlPath="$sock" -o BatchMode=yes "$@"
}
py=../../.venv/bin/python
export OVITO_THREAD_COUNT=4
mkdir -p data audit
touch load-jobs.txt
while read -r name morph seed job <&3; do
    [[ "$job" =~ ^[0-9]+$ ]] || { echo "Invalid preparation job ID for $name" >&2; exit 1; }
    if awk -v name="$name" '$1==name {found=1} END {exit !found}' load-jobs.txt; then
        continue
    fi
    while true; do
        if remote_cmd "test -f runs/$name/prepare.complete"; then
            break
        fi
        job_states=$(remote_cmd "sacct -X -j $job -P -n -o State")
        job_state=${job_states%%$'\n'*}
        case "$job_state" in
            FAILED*|CANCELLED*|TIMEOUT*|NODE_FAIL*|OUT_OF_MEMORY*|BOOT_FAIL*|DEADLINE*|PREEMPTED*)
                printf '%s: unexpected preparation state %s\n' "$name" "$job_state" >&2
                exit 1
                ;;
        esac
        sleep 15
    done
    mkdir -p "data/$name"
    fetch "$host:$remote/runs/$name/initial.dump" "$host:$remote/runs/$name/minimized.dump" \
        "$host:$remote/runs/$name/equilibrated.dump" "$host:$remote/runs/$name/equil.log" "data/$name/"
    "$py" -B analyze.py "data/$name" --out audit
    "$py" - "$name" "$morph" <<'PY'
import csv
import sys
from pathlib import Path

name, morphology = sys.argv[1:]
with (Path("audit") / name / "initial_audit.csv").open() as stream:
    final = list(csv.DictReader(stream))[-1]
assert final["stage"] == "equilibrated", name
assert int(final["line_count"]) == 0, f"{name}: pre-existing DXA line"
contrast = float(final["terrace_contrast_angstrom"])
assert (contrast > 2.0 if morphology == "stepped" else abs(contrast) < 1.0), name
assert int(final["hcp"]) == 0, f"{name}: HCP atoms need review"
print(f"{name}: geometry and initial defect audit passed")
PY
    remote_cmd "if [[ ! -s runs/$name/load-job.id ]]; then sbatch --parsable --job-name=md4-load-$name run.slurm load $name $morph $seed > runs/$name/load-job.id; fi"
    load_job=$(remote_cmd "cat runs/$name/load-job.id")
    [[ "$load_job" =~ ^[0-9]+$ ]] || { echo "Missing valid load job id for $name" >&2; exit 1; }
    printf '%s %s\n' "$name" "$load_job" | tee -a load-jobs.txt
done 3< equil-jobs.txt
while read -r name job <&3; do
    while true; do
        job_states=$(remote_cmd "sacct -X -j $job -P -n -o State,ExitCode")
        job_record=${job_states%%$'\n'*}
        case "$job_record" in
            'COMPLETED|0:0')
                remote_cmd "test -f runs/$name/load.complete"
                break
                ;;
            COMPLETED*|FAILED*|CANCELLED*|TIMEOUT*|NODE_FAIL*|OUT_OF_MEMORY*|BOOT_FAIL*|DEADLINE*|PREEMPTED*)
                printf '%s: unexpected loading state %s\n' "$name" "$job_record" >&2
                exit 1
                ;;
        esac
        sleep 15
    done
    fetch "$host:$remote/runs/$name/loading.dump" "$host:$remote/runs/$name/loading.log" "data/$name/"
    "$py" -B analyze.py "data/$name" --out audit
done 3< load-jobs.txt

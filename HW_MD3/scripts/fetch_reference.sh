#!/usr/bin/env bash
# Obtain precisely the upstream revision used by this experiment.
set -euo pipefail
project_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
source_dir="$project_root/reference/mmmsfall"
source_url="https://gitee.com/sjtu-liugs/mmmsfall.git"
source_revision="544b5c5a86f9119b3bc77dfd3f06c065c3168743"
potential_digest="c0ae005c4680dc6110f6b3ec135f27ab8a2031efacb4a4b52c83e4ec9b54fb24"

if [[ ! -e "$source_dir" ]]; then
    git clone --no-checkout "$source_url" "$source_dir"
    git -C "$source_dir" checkout --detach "$source_revision"
else
    if [[ ! -d "$source_dir/.git" ]]; then
        printf 'Existing path is not an upstream clone: %s\n' "$source_dir" >&2
        exit 1
    fi
    if [[ "$(git -C "$source_dir" rev-parse HEAD)" != "$source_revision" ]]; then
        printf 'Existing upstream revision differs. Preserve it and check manually: %s\n' "$source_dir" >&2
        exit 1
    fi
    if [[ -n "$(git -C "$source_dir" status --porcelain -- 0Surface)" ]]; then
        printf 'Existing 0Surface files have local modifications; leaving them untouched.\n' >&2
        exit 1
    fi
fi

printf '%s  %s\n' "$potential_digest" "$source_dir/0Surface/CuYM.eam.alloy" | sha256sum --check
printf 'Verified course source: %s\n' "$source_revision"

#!/usr/bin/env bash
# Solve the tool envs into explicit locks for linux-64 and linux-aarch64, or check a built
# image against its platform's locks.
#
# Usage: bash docker/update_locks.sh [<env> ...]
#   Solves yamls/<env>.yaml (default: every yaml) for both platforms in a throwaway
#   condaforge/miniforge3 container (docker/solve_locks.py: mamba --dry-run --platform, so
#   nothing is installed and no emulation is needed), with the same package versions on
#   both, writes docker/locks/linux-64/ and docker/locks/linux-aarch64/ (<env>.lock,
#   <env>.json, and <env>.pip.txt / <env>.r_source.tsv when the yaml asks for them), then
#   runs docker/check_locks.py (fails on any version difference between the platforms) and
#   rewrites docker/locks/tool_versions.json. Review with: git diff -- docker/locks
#   Needs network access (conda-forge, bioconda, PyPI). The solve uses the newest packages
#   the yamls allow on the day it runs: pin a version in the yaml to keep it.
#
# Usage: bash docker/update_locks.sh --check-image <image>
#   Compares every env installed in a built image (conda list, run with --network none)
#   with its platform's lock (name, version and build); prints the differences and exits 1
#   when there are any. Nothing is pushed or changed in the image.
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOCKS="$REPO_ROOT/docker/locks"
# Same base as docker/Dockerfile, so the solver is the one the build's conda uses.
SOLVER_IMAGE=condaforge/miniforge3:26.3.2-3@sha256:532f6ee7a858b009dc895f8313eb6ed875f05a455ec57d832aa6b4c66e2799b9

if [ "${1:-}" = --check-image ]; then
    IMAGE="${2:?usage: update_locks.sh --check-image <image>}"
    WORK="$(mktemp -d)"
    trap 'rm -rf "$WORK"' EXIT
    arch="$(docker run --rm --network none "$IMAGE" uname -m)"
    case "$arch" in
        x86_64) plat=linux-64 ;;
        aarch64) plat=linux-aarch64 ;;
        *) echo "unknown image architecture: $arch" >&2; exit 1 ;;
    esac
    for env in $(docker run --rm --network none "$IMAGE" bash -c 'ls /opt/conda/envs'); do
        docker run --rm --network none "$IMAGE" conda list --json -n "$env" > "$WORK/$env.json"
    done
    python3 - "$LOCKS/$plat" "$WORK" <<'PY'
import json, sys
from pathlib import Path
locks, got = Path(sys.argv[1]), Path(sys.argv[2])
bad = 0
for lock in sorted(locks.glob("*.json")):
    env = lock.stem
    if not (got / f"{env}.json").is_file():
        print(f"[check-image] {env}: not in the image")
        bad += 1
        continue
    want = {p["name"]: (p["version"], p["build"]) for p in json.loads(lock.read_text())}
    have = {p["name"]: (p["version"], p["build_string"]) for p in json.loads((got / f"{env}.json").read_text())
            if p.get("channel") != "pypi"}
    for name in sorted(set(want) | set(have)):
        if want.get(name) != have.get(name):
            print(f"[check-image] {env}: {name} lock {want.get(name)} image {have.get(name)}")
            bad += 1
print(f"[check-image] {locks.name}: {bad} differences")
sys.exit(1 if bad else 0)
PY
    exit $?
fi

ENVS=("$@")
if [ "${#ENVS[@]}" -eq 0 ]; then
    for yaml in "$REPO_ROOT"/yamls/*.yaml; do ENVS+=("$(basename "$yaml" .yaml)"); done
fi
echo "[locks] solving ${ENVS[*]} for linux-64 and linux-aarch64"
docker run --rm --network bridge \
    -v "$REPO_ROOT/yamls:/w/yamls:ro" -v "$REPO_ROOT/docker/solve_locks.py:/w/solve_locks.py:ro" \
    -v "$LOCKS:/w/locks" -e HOME=/tmp --user "$(id -u):$(id -g)" \
    "$SOLVER_IMAGE" python /w/solve_locks.py /w/yamls /w/locks "${ENVS[@]}"
python3 "$REPO_ROOT/docker/check_locks.py" "$LOCKS"
python3 "$REPO_ROOT/docker/make_tool_versions.py" locks "$LOCKS" "$LOCKS/tool_versions.json"
echo "[locks] wrote $LOCKS (review with: git diff -- docker/locks)"

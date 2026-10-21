#!/usr/bin/env bash
# Regenerate docker/locks/<env>.lock and docker/locks/tool_versions.json from a built image.
# Usage: bash docker/update_locks.sh <image>
#
# When to use: after changing a yaml in yamls/. Solve the changed envs from their yamls
# with docker/Dockerfile.solve (`docker build --platform linux/amd64 -f docker/Dockerfile.solve
# --build-arg ENVS="<env> ..." -t bactolomew-solve:tmp .`), run this script against that
# image, review the lock diff, then build the real image with docker/build.sh. Only the
# envs in the image (and base) get new locks; tool_versions.json keeps the other envs.
# Against a full image built by docker/build.sh it rewrites every lock (a check that the
# build matches them, see docker/PUBLISHING.md).
# The image runs with --network none and --rm; nothing is pushed or changed in it.
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMAGE="${1:?usage: update_locks.sh <image>}"
LOCKS="$REPO_ROOT/docker/locks"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

ENVS="$(docker run --rm --network none "$IMAGE" bash -c 'ls /opt/conda/envs')"
for env in base $ENVS; do
    echo "[locks] $env"
    docker run --rm --network none "$IMAGE" conda list --explicit --md5 -n "$env" > "$WORK/$env.lock"
    docker run --rm --network none "$IMAGE" conda list --json -n "$env" > "$WORK/$env.json"
done
for env in base $ENVS; do
    cp "$WORK/$env.lock" "$LOCKS/$env.lock"
done
IMAGE_ID="$(docker image inspect --format '{{.Id}}' "$IMAGE")"
IMAGE_CREATED="$(docker image inspect --format '{{.Created}}' "$IMAGE")"
python3 "$REPO_ROOT/docker/make_tool_versions.py" "$WORK" "$LOCKS/tool_versions.json" \
    --image "$IMAGE" --image-id "$IMAGE_ID" --image-created "$IMAGE_CREATED"
echo "[locks] wrote $LOCKS (review with: git diff -- docker/locks)"

#!/usr/bin/env bash
# Build the BACTolomew unified tool image locally from the repo root.
# Usage: bash docker/build.sh [version] [image-name]
#   version     image tag and org.opencontainers.image.version (default 0.2.0)
#   image-name  repository name (default bactolomew-tools), e.g. ghcr.io/<owner>/bactolomew-tools
# Environment: FASTGEAR_IMAGE=<ref>@sha256:<digest> builds from a mirror of the fastgear
#   base image instead of Docker Hub (see docker/PUBLISHING.md).
# Builds only, always for linux/amd64; nothing is pushed (see docker/PUBLISHING.md).
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VERSION="${1:-0.2.0}"
NAME="${2:-bactolomew-tools}"
IMAGE="$NAME:$VERSION"
cd "$REPO_ROOT"

REVISION="$(git rev-parse HEAD 2>/dev/null || echo unknown)"
if [ "$REVISION" != unknown ] && [ -n "$(git status --porcelain -- docker yamls .dockerignore 2>/dev/null)" ]; then
    # The image content comes from docker/ and yamls/ only; flag local edits there.
    REVISION="$REVISION-dirty"
    echo "[build] WARNING: docker/, yamls/ or .dockerignore has uncommitted changes; the" >&2
    echo "[build] revision label will be $REVISION. Commit before building an image to publish." >&2
fi
CREATED="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
SOURCE="$(git remote get-url origin 2>/dev/null || echo https://github.com/liviurotiul/BACTolomew)"
SOURCE="${SOURCE%.git}"
case "$SOURCE" in
    git@github.com:*) SOURCE="https://github.com/${SOURCE#git@github.com:}" ;;
esac
# The source becomes a public label: drop any user:token@ part of the remote URL.
SOURCE="$(printf '%s' "$SOURCE" | sed -E 's#^([a-z+]+://)[^/@]*@#\1#')"
FASTGEAR_ARGS=()
if [ -n "${FASTGEAR_IMAGE:-}" ]; then
    FASTGEAR_ARGS=(--build-arg FASTGEAR_IMAGE="$FASTGEAR_IMAGE")
fi

echo "[build] building $IMAGE from $REPO_ROOT"
echo "[build] revision $REVISION, created $CREATED, source $SOURCE"
docker build --platform linux/amd64 -f docker/Dockerfile \
    "${FASTGEAR_ARGS[@]}" \
    --build-arg VERSION="$VERSION" \
    --build-arg REVISION="$REVISION" \
    --build-arg CREATED="$CREATED" \
    --build-arg SOURCE="$SOURCE" \
    -t "$IMAGE" .
echo "[build] done: $IMAGE"
echo "[build] image ID: $(docker image inspect --format '{{.Id}}' "$IMAGE")"

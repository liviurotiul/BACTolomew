#!/usr/bin/env bash
# Build the BACTolomew unified tool image locally from the repo root.
# Usage: bash docker/build.sh [--platform linux/amd64|linux/arm64] [version] [image-name]
#   version     image tag and org.opencontainers.image.version (default 0.3.0; 0.2.0 is
#               built from commit 8142861, see PUBLISHING.md)
#   image-name  repository name (default bactolomew), e.g. liviurotiul/bactolomew
# Environment: FASTGEAR_IMAGE=<ref>@sha256:<digest> builds from a mirror of the fastgear
#   base image instead of Docker Hub (see docker/PUBLISHING.md).
#   --platform (or PLATFORM=) picks the platform, default the Docker engine's own (so an
#   Apple Silicon Mac builds linux/arm64 natively, never under emulation); the lock folder
#   follows it (docker/locks/linux-64 or linux-aarch64, the same tool versions). fastGEAR
#   (an x86_64 MATLAB runtime) is included only for linux/amd64; FASTGEAR=none leaves it
#   out there too.
# After the build the image gets the label org.bactolomew.population_methods, the methods
# verify_image.sh found in it. The build fails unless they are the expected ones:
# fastgear,fastbaps for linux/amd64, fastbaps for linux/arm64 or FASTGEAR=none.
# Builds only; nothing is pushed (see docker/PUBLISHING.md).
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [ "${1:-}" = --platform ]; then
    PLATFORM="${2:?--platform needs a value (linux/amd64 or linux/arm64)}"
    shift 2
fi
VERSION="${1:-0.3.0}"
NAME="${2:-bactolomew}"
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
PLATFORM="${PLATFORM:-linux/$(docker version --format '{{.Server.Arch}}')}"
case "$PLATFORM" in
    linux/amd64|linux/arm64) ;;
    *) echo "[build] unsupported platform $PLATFORM (linux/amd64 or linux/arm64)" >&2; exit 2 ;;
esac
FASTGEAR_ARGS=()
EXPECTED_METHODS=fastgear,fastbaps
if [ "$PLATFORM" != linux/amd64 ] || [ "${FASTGEAR:-}" = none ]; then
    FASTGEAR_ARGS=(--build-arg FASTGEAR=none --build-arg FASTGEAR_PROVENANCE=)
    EXPECTED_METHODS=fastbaps
    echo "[build] building without fastGEAR ($PLATFORM${FASTGEAR:+, FASTGEAR=$FASTGEAR})"
elif [ -n "${FASTGEAR_IMAGE:-}" ]; then
    FASTGEAR_ARGS=(--build-arg FASTGEAR_IMAGE="$FASTGEAR_IMAGE")
fi

echo "[build] building $IMAGE for $PLATFORM from $REPO_ROOT"
echo "[build] revision $REVISION, created $CREATED, source $SOURCE"
docker build --platform "$PLATFORM" -f docker/Dockerfile \
    "${FASTGEAR_ARGS[@]}" \
    --build-arg VERSION="$VERSION" \
    --build-arg REVISION="$REVISION" \
    --build-arg CREATED="$CREATED" \
    --build-arg SOURCE="$SOURCE" \
    -t "$IMAGE" .
# The methods are known only inside the built image: add them as a label (a metadata-only
# image on top, no new layer).
METHODS="$(docker run --rm --network none --platform "$PLATFORM" "$IMAGE" /opt/conda/bin/python -c \
    'import json; print(",".join(json.load(open("/opt/bactolomew/image_contents.json"))["population_methods"]))')"
if [ "$METHODS" != "$EXPECTED_METHODS" ]; then
    echo "[build] FAILED: $IMAGE holds the population methods '$METHODS', expected '$EXPECTED_METHODS'" >&2
    exit 1
fi
printf 'FROM %s\nLABEL org.bactolomew.population_methods="%s"\n' "$IMAGE" "$METHODS" \
    | docker build --platform "$PLATFORM" -q -t "$IMAGE" - > /dev/null
echo "[build] done: $IMAGE (population methods: $METHODS)"
echo "[build] image ID: $(docker image inspect --format '{{.Id}}' "$IMAGE")"

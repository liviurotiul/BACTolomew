#!/usr/bin/env bash
# Build the BACTolomew unified tool image locally from the repo root.
# Usage: bash docker/build.sh [image:tag]   (default: bactolomew-tools:0.1.0)
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMAGE="${1:-bactolomew-tools:0.1.0}"
cd "$REPO_ROOT"
echo "[build] building $IMAGE from $REPO_ROOT ..."
docker build -f docker/Dockerfile -t "$IMAGE" .
echo "[build] done: $IMAGE"

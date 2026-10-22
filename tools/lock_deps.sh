#!/usr/bin/env bash
# Pin the host Python dependencies: writes requirements.lock (GUI + pipeline)
# and requirements-pipeline.lock (headless) from requirements.txt and
# requirements-pipeline.txt.
#
# The locks are universal: one file covers every supported Python (3.10 to
# 3.14) on every platform, with environment markers where a package needs a
# different version per Python (numpy, pandas, scipy, pyarrow and matplotlib
# dropped 3.10 in their newest releases), and sha256 hashes for every file.
# install.sh installs them with:
#   pip install --require-hashes --only-binary :all: -r requirements.lock
#
# Uses uv (https://github.com/astral-sh/uv) for the universal resolution: the
# one on PATH, the one named by $UV, or else a temporary copy installed into a
# throwaway venv that is removed afterwards. Only the package index is
# contacted; nothing is uploaded.
#
# Usage: tools/lock_deps.sh [--upgrade]
#   --upgrade   re-resolve to the newest versions within the bounds; without
#               it, versions already in the locks are kept where still allowed.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MIN_PYTHON="3.10"
UPGRADE=()
for arg in "$@"; do
    case "$arg" in
        --upgrade) UPGRADE=(--upgrade) ;;
        -h|--help) sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "lock_deps: unknown option: $arg" >&2; exit 2 ;;
    esac
done

TMP=""
cleanup() { [ -z "$TMP" ] || rm -rf "$TMP"; }
trap cleanup EXIT

uv_bin="${UV:-}"
if [ -z "$uv_bin" ]; then
    uv_bin="$(command -v uv || true)"
fi
if [ -z "$uv_bin" ]; then
    py="$(command -v python3 || command -v python || true)"
    if [ -z "$py" ]; then
        echo "lock_deps: needs uv, or a python3 to install it into a temporary venv" >&2
        exit 1
    fi
    TMP="$(mktemp -d "${TMPDIR:-/tmp}/bactolomew-lock.XXXXXX")"
    echo "lock_deps: installing uv into a temporary venv ($TMP)"
    "$py" -m venv "$TMP/venv"
    "$TMP/venv/bin/pip" install --quiet uv
    uv_bin="$TMP/venv/bin/uv"
fi
echo "lock_deps: $("$uv_bin" --version)"

lock() {
    local src="$1" out="$2"
    echo "lock_deps: $src -> $out"
    (cd "$ROOT" && "$uv_bin" pip compile "$src" \
        --universal \
        --python-version "$MIN_PYTHON" \
        --only-binary :all: \
        --generate-hashes \
        --no-strip-markers \
        --custom-compile-command "tools/lock_deps.sh" \
        ${UPGRADE[@]+"${UPGRADE[@]}"} \
        --quiet \
        --output-file "$out")
}

lock requirements-pipeline.txt requirements-pipeline.lock
lock requirements.txt requirements.lock

echo "lock_deps: pinned versions (requirements.lock):"
grep -E '^[A-Za-z0-9_.-]+==' "$ROOT/requirements.lock" | sed 's/ \\$//; s/^/  /'

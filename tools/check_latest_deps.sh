#!/usr/bin/env bash
# Local check of the newest dependency releases (no CI, nothing uploaded).
#
# Builds two throwaway venvs next to each other:
#   bounded  the newest versions requirements.txt allows
#   latest   the newest versions released, ignoring the upper bounds
# and in each one runs tools/deps_smoke.py (imports, a synthetic run through
# the parquet cache, RunData and QueryEngine) and a short list of test
# modules offscreen. It then prints which packages moved from
# requirements.lock and whether anything broke.
#
# Usage: tools/check_latest_deps.sh [--python PY] [--keep] [--work DIR] [--tests "mod ..."]
#   --python PY  interpreter for the venvs (default: python3); run it once per
#                supported Python (3.10 to 3.14) to cover the range
#   --work DIR   folder for the venvs and logs (default: a new temp folder)
#   --keep       keep the work folder afterwards
#   --tests      test modules to run (default: the list below)
# Exit status: 0 when both venvs pass, 1 when the bounded one fails, 2 when
# only the latest one fails (a bound or the code needs updating before the
# next lock).
#
# uv is used when on PATH (or named by $UV) for speed; otherwise pip.
# The user's QSettings, cache/ and runs/ are never touched (tests isolate
# QSettings; the smoke test writes its cache to the work folder).
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="python3"
WORK=""
KEEP=0
TESTS="test_run_data test_explore_context test_pan_panels test_cgmlst_force_plot test_tree_sources test_annotation_filter_controls"

while [ $# -gt 0 ]; do
    case "$1" in
        --python) PY="$2"; shift 2 ;;
        --work) WORK="$2"; shift 2 ;;
        --keep) KEEP=1; shift ;;
        --tests) TESTS="$2"; shift 2 ;;
        -h|--help) sed -n '2,25p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "check_latest_deps: unknown option: $1" >&2; exit 2 ;;
    esac
done

if [ -z "$WORK" ]; then
    WORK="$(mktemp -d "${TMPDIR:-/tmp}/bactolomew-latest.XXXXXX")"
else
    mkdir -p "$WORK"
fi
if [ "$KEEP" -eq 0 ]; then
    trap 'rm -rf "$WORK"' EXIT
fi
echo "check_latest_deps: work folder $WORK ($("$PY" -V 2>&1))"

uv_bin="${UV:-$(command -v uv || true)}"

make_venv() {  # make_venv DIR
    if [ -n "$uv_bin" ]; then
        "$uv_bin" venv --quiet --python "$PY" "$1"
    else
        "$PY" -m venv "$1" && "$1/bin/pip" install --quiet --upgrade pip
    fi
}

pip_install() {  # pip_install VENV ARGS...
    local venv="$1"; shift
    if [ -n "$uv_bin" ]; then
        VIRTUAL_ENV="$venv" "$uv_bin" pip install --quiet --only-binary :all: "$@"
    else
        "$venv/bin/pip" install --quiet --only-binary :all: "$@"
    fi
}

freeze() {  # freeze VENV -> "name==version" lines, lower case names
    if [ -n "$uv_bin" ]; then
        VIRTUAL_ENV="$1" "$uv_bin" pip freeze 2>/dev/null
    else
        "$1/bin/pip" freeze 2>/dev/null
    fi | tr 'A-Z_' 'a-z-' | sort
}

# The direct dependencies without their version bounds.
unbounded_requirements() {
    cat "$ROOT/requirements-pipeline.txt" "$ROOT/requirements.txt" \
        | sed 's/#.*//' | grep -v '^\s*-r' | grep -v '^\s*$' \
        | sed -E 's/[<>=~!;[:space:]].*$//' | sort -u
}

run_checks() {  # run_checks NAME VENV -> 0 ok, 1 failed
    local name="$1" venv="$2" status=0 log
    log="$WORK/$name.smoke.log"
    if (cd "$ROOT" && QT_QPA_PLATFORM=offscreen PYTHONPATH=src:tests PYTHONDONTWRITEBYTECODE=1 \
            "$venv/bin/python" tools/deps_smoke.py >"$log" 2>&1); then
        grep -E '^(ok|FAIL)' "$log" | sed "s/^/  [$name] /"
    else
        grep -E '^(ok|FAIL)' "$log" | sed "s/^/  [$name] /"
        tail -n 15 "$log" | sed "s/^/  [$name]   /"
        status=1
    fi
    log="$WORK/$name.tests.log"
    # shellcheck disable=SC2086
    if (cd "$ROOT" && QT_QPA_PLATFORM=offscreen PYTHONPATH=src:tests PYTHONDONTWRITEBYTECODE=1 \
            PYTHONFAULTHANDLER=1 timeout 1800 "$venv/bin/python" -m unittest $TESTS >"$log" 2>&1); then
        echo "  [$name] tests: $(grep -E '^Ran ' "$log"), $(grep -E '^OK' "$log" | tail -n 1)"
    else
        echo "  [$name] tests FAILED (log: $log)"
        grep -E '^(ERROR|FAIL):|^Fatal Python error' "$log" | head -n 20 | sed "s/^/  [$name]   /"
        tail -n 3 "$log" | sed "s/^/  [$name]   /"
        status=1
    fi
    return "$status"
}

pinned="$WORK/pinned.txt"
grep -E '^[A-Za-z0-9_.-]+==' "$ROOT/requirements.lock" | sed -E 's/ .*//; s/ \\$//' \
    | tr 'A-Z_' 'a-z-' | sort -u >"$pinned"

declare -A RESULT
for name in bounded latest; do
    venv="$WORK/$name"
    echo
    echo "== $name: installing"
    make_venv "$venv" || { echo "  could not create the venv"; RESULT[$name]=1; continue; }
    if [ "$name" = bounded ]; then
        pip_install "$venv" --upgrade -r "$ROOT/requirements.txt" >"$WORK/$name.install.log" 2>&1
    else
        # shellcheck disable=SC2046
        pip_install "$venv" --upgrade $(unbounded_requirements) >"$WORK/$name.install.log" 2>&1
    fi
    if [ $? -ne 0 ]; then
        echo "  install FAILED:"; tail -n 15 "$WORK/$name.install.log" | sed 's/^/    /'
        RESULT[$name]=1
        continue
    fi
    freeze "$venv" >"$WORK/$name.freeze.txt"
    echo "== $name: versions that differ from requirements.lock"
    moved=0
    while IFS='=' read -r pkg _ ver; do
        locked="$(grep -E "^${pkg}==" "$pinned" | head -n 1 | sed 's/.*==//')"
        # A universal lock pins some packages once per Python: compare with each.
        if [ -n "$locked" ] && ! grep -qxF "${pkg}==${ver}" "$pinned"; then
            echo "    $pkg: $(grep -E "^${pkg}==" "$pinned" | sed 's/.*==//' | paste -sd/ -) -> $ver"
            moved=1
        elif [ -z "$locked" ]; then
            echo "    $pkg: (not in the lock) $ver"
            moved=1
        fi
    done <"$WORK/$name.freeze.txt"
    [ "$moved" -eq 0 ] && echo "    (none)"
    echo "== $name: checks"
    run_checks "$name" "$venv"
    RESULT[$name]=$?
done

echo
echo "== summary ($("$PY" -V 2>&1))"
for name in bounded latest; do
    if [ "${RESULT[$name]:-1}" -eq 0 ]; then
        echo "  $name: OK"
    else
        echo "  $name: BROKEN (logs in $WORK; rerun with --keep to look at them)"
    fi
done
if [ "${RESULT[bounded]:-1}" -ne 0 ]; then
    exit 1
elif [ "${RESULT[latest]:-1}" -ne 0 ]; then
    exit 2
fi
exit 0

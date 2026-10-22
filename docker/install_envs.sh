#!/usr/bin/env bash
# Create conda envs from their explicit locks, then trim, all inside ONE Docker layer.
# Usage (inside the image build): install_envs.sh <env> [<env> ...]
#
# Each env is created from /opt/bactolomew/locks/<env>.lock (the target platform's explicit
# lock from docker/locks/<platform>/, written by docker/update_locks.sh): no solver runs,
# every package is the exact URL and md5 recorded in the lock. Then, when the lock folder
# has them:
#   * <env>.pip.txt: PyPI wheels, each pinned by version and sha256 (pip --require-hashes,
#     wheels only, no dependencies: the lock already holds them);
#   * <env>.r_source.tsv: R packages built from source, one per line (name, version,
#     sha256, then one or more URLs of the same tarball), checked by sha256 and installed
#     with the env's R and compilers, in the listed order.
#
# Size reductions applied here are only the ones that cannot change what a tool does:
#   * envs created in the same layer share their packages: conda hard-links each file from
#     its package cache, and a layer stores a hard-linked file once (a later layer cannot
#     link into an earlier one, which is why the envs are grouped per RUN);
#   * the package cache is emptied in the same layer (`conda clean` keeps copies conda made
#     instead of links, so /opt/conda/pkgs is removed outright; doing it in a later layer
#     would free nothing);
#   * documentation (share/doc, share/man, share/info, share/gtk-doc) and the static
#     libraries directly in lib/ (lib*.a, used only when linking a program; the compilers'
#     own files deeper under lib/ stay) are deleted.
# Not touched on purpose: *.pyc (rebuilding them is impossible on a read-only container
# filesystem), share/locale, include/, test folders, the mob_suite database, prokka and
# barrnap in the panaroo env, and the compilers R pulls into gubbins. See docker/PUBLISHING.md.
set -euo pipefail

LOCKS=/opt/bactolomew/locks
[ "$#" -gt 0 ] || { echo "usage: install_envs.sh <env> [<env> ...]" >&2; exit 2; }

for env in "$@"; do
    lock="$LOCKS/$env.lock"
    [ -f "$lock" ] || { echo "missing lock: $lock" >&2; exit 1; }
    echo "[install_envs] $env <- $lock"
    conda create --yes --quiet --name "$env" --file "$lock"
    if [ -f "$LOCKS/$env.pip.txt" ]; then
        echo "[install_envs] $env <- $LOCKS/$env.pip.txt"
        conda run -n "$env" python -m pip install --no-deps --no-cache-dir --require-hashes \
            --only-binary=:all: --disable-pip-version-check -r "$LOCKS/$env.pip.txt"
    fi
    if [ -f "$LOCKS/$env.r_source.tsv" ]; then
        src="$(mktemp -d)"
        # Packages that ask for C++11 (fastbaps 1.0.8: CXX_STD = CXX11) are compiled as
        # C++14, which the RcppArmadillo 15 headers they include require.
        printf 'CXX11STD = -std=gnu++14\n' > "$src/Makevars"
        export R_MAKEVARS_USER="$src/Makevars"
        while IFS=$'\t' read -r name version sha256 urls; do
            [ -n "$name" ] || continue
            echo "[install_envs] $env <- R source $name $version"
            # shellcheck disable=SC2086
            /opt/conda/bin/python - "$src/$name.tar.gz" "$sha256" $urls <<'PY'
import hashlib, sys, urllib.request
dest, digest, urls = sys.argv[1], sys.argv[2], sys.argv[3:]
for url in urls:
    try:
        data = urllib.request.urlopen(url, timeout=300).read()
    except OSError as exc:
        print(f"[install_envs]   {url}: {exc}", file=sys.stderr)
        continue
    if hashlib.sha256(data).hexdigest() == digest:
        open(dest, "wb").write(data)
        sys.exit(0)
    print(f"[install_envs]   {url}: sha256 differs from the lock", file=sys.stderr)
sys.exit(f"no source with the sha256 {digest}")
PY
            conda run -n "$env" R CMD INSTALL --no-docs --no-html "$src/$name.tar.gz"
            conda run -n "$env" Rscript -e "stopifnot(packageVersion('$name') == '$version')"
        done < "$LOCKS/$env.r_source.tsv"
        unset R_MAKEVARS_USER
        rm -rf "$src"
    fi
done

conda clean --all --force-pkgs-dirs --yes
rm -rf /opt/conda/pkgs/*

for env in "$@"; do
    prefix="/opt/conda/envs/$env"
    rm -rf "$prefix/share/doc" "$prefix/share/man" "$prefix/share/info" "$prefix/share/gtk-doc" "$prefix/man"
    find "$prefix/lib" -maxdepth 1 -type f -name 'lib*.a' -delete
done

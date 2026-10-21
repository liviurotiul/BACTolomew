#!/usr/bin/env bash
# Create conda envs from their explicit locks, then trim, all inside ONE Docker layer.
# Usage (inside the image build): install_envs.sh <env> [<env> ...]
#
# Each env is created from /opt/bactolomew/locks/<env>.lock (`conda list --explicit --md5`
# output): no solver runs, every package is the exact URL and md5 recorded in the lock.
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
done

conda clean --all --force-pkgs-dirs --yes
rm -rf /opt/conda/pkgs/*

for env in "$@"; do
    prefix="/opt/conda/envs/$env"
    rm -rf "$prefix/share/doc" "$prefix/share/man" "$prefix/share/info" "$prefix/share/gtk-doc" "$prefix/man"
    find "$prefix/lib" -maxdepth 1 -type f -name 'lib*.a' -delete
done

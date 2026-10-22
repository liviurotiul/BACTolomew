"""Solve yamls/<env>.yaml into explicit locks for linux-64 and linux-aarch64, with the
same package versions on both platforms.

Runs inside a throwaway condaforge/miniforge3 container (docker/update_locks.sh starts it);
both platforms are solved with `mamba create --dry-run --platform ...`, so nothing is
installed and no emulation is needed.

Per env:
  1. linux-aarch64 (the more restrictive platform) is solved from the yaml's conda
     dependencies (the yaml pins the tools and key dependencies to versions built for
     both platforms);
  2. linux-64 is solved from the same dependencies; then a package the two solves hold
     with different versions is pinned on both (`name==version`, the older version, else
     the newer; platform-named packages such as ld_impl_linux-64 / ld_impl_linux-aarch64
     are pinned under their own name per platform) and a package only one solve holds is
     requested on the other too, and both are solved again until they agree
     (docker/check_locks.py then fails on any difference left);
  3. the packages in the yaml's `x-bactolomew: no_deps` (installed without their
     dependencies, e.g. a recipe whose pins have no aarch64 build) are added to both
     locks from the channel index (the highest build of the requested version);
  4. `pip:` entries (exact `name==version`) become docker/locks/<platform>/<env>.pip.txt
     with the sha256 of the manylinux wheel for that platform (from PyPI);
  5. `x-bactolomew: r_source` entries (R packages built from source, the same tarball on
     every platform) become docker/locks/<platform>/<env>.r_source.tsv.

Writes, per platform: <env>.lock (`@EXPLICIT`, url#md5 lines in install order) and
<env>.json (name, version, build, channel, sha256 of every package, read by
docker/make_tool_versions.py and docker/check_locks.py).

Usage (inside the container): python solve_locks.py <yamls dir> <locks dir> <env> [<env> ...]
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import urllib.request
from pathlib import Path

PLATFORMS = ("linux-aarch64", "linux-64")
CHANNELS = ["--override-channels", "-c", "conda-forge", "-c", "bioconda"]
PIP_WHEEL_TAGS = {"linux-64": "x86_64", "linux-aarch64": "aarch64"}
# The solve must not depend on the machine it runs on: no CUDA builds (tensorflow cpu),
# the glibc of the image's Ubuntu 24.04, generic microarchitecture.
SOLVE_ENV = {"CONDA_OVERRIDE_CUDA": "", "CONDA_OVERRIDE_GLIBC": "2.39"}
ARCHSPEC = {"linux-64": "x86_64", "linux-aarch64": "aarch64"}


def load_yaml(path: Path) -> dict:
    from ruamel.yaml import YAML  # conda's own dependency, in the miniforge base env

    return YAML(typ="safe").load(path.read_text())


def platform_name(name: str, plat: str) -> str:
    """ld_impl_linux-aarch64 -> ld_impl_linux-64 (and back): the same package per platform."""
    other = "linux-64" if plat == "linux-aarch64" else "linux-aarch64"
    return re.sub(rf"_{re.escape(other)}$", f"_{plat}", name)


def solve(plat: str, specs: list[str]) -> list[dict]:
    cmd = ["mamba", "create", "-n", "solve", "--dry-run", "--json", "--yes", "--platform", plat,
           *CHANNELS, *specs]
    env = {**os.environ, **SOLVE_ENV, "CONDA_SUBDIR": plat, "CONDA_OVERRIDE_ARCHSPEC": ARCHSPEC[plat]}
    res = subprocess.run(cmd, capture_output=True, text=True, env=env)
    try:
        doc = json.loads(res.stdout)
    except ValueError:
        doc = {}
    if not doc.get("success"):
        raise RuntimeError(f"{plat}: solve failed for {' '.join(specs)}\n{res.stdout[-4000:]}{res.stderr[-2000:]}")
    return list(doc["actions"]["LINK"])


def no_deps_record(plat: str, spec: str) -> dict:
    cmd = ["mamba", "search", "--json", "--platform", plat, *CHANNELS, spec]
    env = {**os.environ, **SOLVE_ENV, "CONDA_SUBDIR": plat, "CONDA_OVERRIDE_ARCHSPEC": ARCHSPEC[plat]}
    res = subprocess.run(cmd, capture_output=True, text=True, env=env, check=True)
    pkgs = json.loads(res.stdout)["result"]["pkgs"]
    if not pkgs:
        raise RuntimeError(f"{plat}: no package for {spec}")
    return max(pkgs, key=lambda p: (p.get("build_number", 0), p.get("timestamp", 0)))


def pip_lines(plat: str, specs: list[str], python: str) -> list[str]:
    """`name==version --hash=sha256:...` of the manylinux wheel for this platform and the
    env's Python (cpXY)."""
    tag = "cp" + "".join(python.split(".")[:2])
    lines = []
    for spec in specs:
        name, _, version = spec.partition("==")
        if not version:
            raise RuntimeError(f"pip entry {spec!r}: pin an exact version (name==version)")
        meta = json.load(urllib.request.urlopen(f"https://pypi.org/pypi/{name}/{version}/json"))
        arch = PIP_WHEEL_TAGS[plat]
        wheels = [f for f in meta["urls"] if f["packagetype"] == "bdist_wheel" and f"-{tag}-" in f["filename"]
                  and "manylinux" in f["filename"] and f["filename"].endswith(f"{arch}.whl")]
        if not wheels:
            raise RuntimeError(f"{plat}: no manylinux {arch} {tag} wheel of {spec}")
        lines.append(f"{name}=={version} " + " ".join(f"--hash=sha256:{w['digests']['sha256']}" for w in wheels))
    return lines


def versions(pkgs: list[dict], plat: str) -> dict[str, str]:
    """name -> version, platform-named packages under their linux-64 name."""
    return {platform_name(p["name"], "linux-64") if plat == "linux-aarch64" else p["name"]: p["version"]
            for p in pkgs}


def align(env: str, specs: list[str], rounds: int = 6) -> dict[str, list[dict]]:
    """Solve both platforms, then make them agree: a package both solves hold with
    different versions is pinned on both to the older of the two (else the newer), and a
    package only one solve holds is asked for on the other platform too, at its version.
    All the new pins are tried at once, then one at a time (keeping the ones that solve),
    until nothing changes. What still differs is left for docker/check_locks.py."""
    from conda.models.version import VersionOrder

    def solve_both(pins: dict[str, str]) -> dict[str, list[dict]]:
        extra = [f"{n}=={v}" for n, v in sorted(pins.items())]
        out = {}
        for plat in PLATFORMS:
            named = [platform_name(e.split("==")[0], plat) + "==" + e.split("==")[1] for e in extra]
            out[plat] = solve(plat, specs + named)
        return out

    print(f"[solve] {env}", flush=True)
    solved = solve_both({})
    pins: dict[str, str] = {}
    tried: set[tuple[str, str]] = set()
    for _ in range(rounds):
        va, vx = versions(solved["linux-aarch64"], "linux-aarch64"), versions(solved["linux-64"], "linux-64")
        wanted: list[tuple[str, list[str]]] = []
        for name in sorted(set(va) | set(vx)):
            if name in va and name in vx:
                if va[name] != vx[name]:
                    wanted.append((name, sorted({va[name], vx[name]}, key=VersionOrder)))
            else:
                wanted.append((name, [(va | vx)[name]]))
        wanted = [(n, [v for v in vs if (n, v) not in tried]) for n, vs in wanted]
        wanted = [(n, vs) for n, vs in wanted if vs]
        if not wanted:
            break
        print(f"[solve] {env}: pinning {' '.join(f'{n}={vs[0]}' for n, vs in wanted)}", flush=True)
        try:
            solved = solve_both({**pins, **{n: vs[0] for n, vs in wanted}})
            pins.update({n: vs[0] for n, vs in wanted})
            continue
        except RuntimeError:
            pass
        for name, choices in wanted:  # one at a time: keep the pins that solve
            for version in choices:
                tried.add((name, version))
                try:
                    solved = solve_both({**pins, name: version})
                    pins[name] = version
                    break
                except RuntimeError:
                    print(f"[solve] {env}: cannot pin {name}={version} on both platforms", flush=True)
    return solved


def record(pkg: dict) -> dict:
    return {k: pkg.get(k) for k in ("name", "version", "build", "channel", "subdir", "url", "md5", "sha256")}


def main() -> None:
    yamls, locks = Path(sys.argv[1]), Path(sys.argv[2])
    for env in sys.argv[3:]:
        spec = load_yaml(yamls / f"{env}.yaml")
        conda_specs = [str(d) for d in spec.get("dependencies", []) if isinstance(d, str)]
        pip_specs = [str(p) for d in spec.get("dependencies", []) if isinstance(d, dict) for p in d.get("pip", [])]
        extra = spec.get("x-bactolomew") or {}
        solved = align(env, conda_specs)
        for plat in PLATFORMS:
            solved[plat] += [no_deps_record(plat, s) for s in extra.get("no_deps", [])]
        for plat, pkgs in solved.items():
            out = locks / plat
            out.mkdir(parents=True, exist_ok=True)
            body = [f"# {env}: solved from yamls/{env}.yaml by docker/update_locks.sh (docker/solve_locks.py)",
                    f"# platform: {plat}", "@EXPLICIT"]
            body += [f"{p['url']}#{p['md5']}" for p in pkgs]
            (out / f"{env}.lock").write_text("\n".join(body) + "\n")
            (out / f"{env}.json").write_text(json.dumps([record(p) for p in pkgs], indent=1) + "\n")
            python = next((str(p["version"]) for p in pkgs if p["name"] == "python"), "")
            for suffix, lines in ((".pip.txt", pip_lines(plat, pip_specs, python) if pip_specs else []),
                                  (".r_source.tsv", [
                                      "\t".join([r["name"], str(r["version"]), r["sha256"], *r["urls"]])
                                      for r in extra.get("r_source", [])])):
                path = out / f"{env}{suffix}"
                if lines:
                    path.write_text("\n".join(lines) + "\n")
                elif path.exists():
                    path.unlink()
            print(f"[solve] {env} {plat}: {len(pkgs)} packages", flush=True)


if __name__ == "__main__":
    main()

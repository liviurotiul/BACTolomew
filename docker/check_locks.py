"""Check that the linux-64 and linux-aarch64 locks hold the same packages and versions.

Usage: python3 docker/check_locks.py [docker/locks]

Compares docker/locks/linux-64/<env>.json with docker/locks/linux-aarch64/<env>.json (written
by docker/solve_locks.py) by package name and version, and the pip and R source pins of
each env. Exits 1 on any difference other than the platform-only packages listed below.
Packages named after their platform (ld_impl_linux-64 / ld_impl_linux-aarch64,
sysroot_linux-*, gcc_impl_linux-*, ...) count as one package and must have the same version.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

PLATFORMS = ("linux-64", "linux-aarch64")
# Packages that may be in one platform's lock only (with the platform they belong to).
# Each is a platform runtime piece, not a tool or a library a tool's results depend on.
PLATFORM_ONLY: dict[str, str] = {
    # The x86_64 microarchitecture level selector (a metapackage, no files; aarch64 has none).
    "_x86_64-microarch-level": "linux-64",
    # pango 1.58.2 is in both locks; its newer linux-64 build moved the library files into
    # these two packages (same version), the linux-aarch64 build still holds them itself.
    "libpango": "linux-64",
    "libpango-devel": "linux-64",
    # A Perl build helper that turns scripts into Windows .bat files (pulled in by a Perl
    # module's build requirements on linux-64; never run).
    "perl-extutils-pl2bat": "linux-64",
}


def normalize(name: str) -> str:
    return re.sub(r"_linux-(64|aarch64)$", "_linux-<platform>", name)


def load_env(locks: Path, plat: str, env: str) -> dict[str, str]:
    pkgs = json.loads((locks / plat / f"{env}.json").read_text())
    return {normalize(p["name"]): str(p["version"]) for p in pkgs}


def pins(path: Path, column: int = 0) -> list[str]:
    if not path.is_file():
        return []
    return sorted(line.split()[column] if column == 0 else "\t".join(line.split("\t")[:3])
                  for line in path.read_text().splitlines() if line.strip())


def main() -> int:
    locks = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent / "locks"
    envs = [{p.stem for p in (locks / plat).glob("*.json")} for plat in PLATFORMS]
    problems: list[str] = []
    for env in sorted(envs[0] ^ envs[1]):
        problems.append(f"{env}: lock only for {PLATFORMS[0] if env in envs[0] else PLATFORMS[1]}")
    for env in sorted(envs[0] & envs[1]):
        x86, arm = (load_env(locks, plat, env) for plat in PLATFORMS)
        for name in sorted(set(x86) | set(arm)):
            if name in x86 and name in arm:
                if x86[name] != arm[name]:
                    problems.append(f"{env}: {name} {x86[name]} (linux-64) != {arm[name]} (linux-aarch64)")
                continue
            plat = "linux-64" if name in x86 else "linux-aarch64"
            if PLATFORM_ONLY.get(name) != plat:
                problems.append(f"{env}: {name} {(x86 | arm)[name]} only in {plat}")
        for suffix, column in ((".pip.txt", 0), (".r_source.tsv", 1)):
            a, b = (pins(locks / plat / f"{env}{suffix}", column) for plat in PLATFORMS)
            if a != b:
                problems.append(f"{env}{suffix}: {a} (linux-64) != {b} (linux-aarch64)")
    for line in problems:
        print(f"[check_locks] {line}")
    total = len(envs[0] & envs[1])
    print(f"[check_locks] {total} envs, {len(problems)} differences" + ("" if problems else ": the platforms match"))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())

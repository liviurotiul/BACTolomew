"""Write a tool version summary: from the locks (docker/locks/tool_versions.json) or from
what an image really holds (/opt/bactolomew/tool_versions.json, at build time).

Usage:
    python3 docker/make_tool_versions.py locks <docker/locks> <dest>
        (docker/update_locks.sh) the main tools of every env, with their version (the same
        on every platform, see docker/check_locks.py), channel and build per platform;
    python3 make_tool_versions.py image <locks summary> <dest>
        (docker/verify_image.sh, inside the image) only the envs installed in
        /opt/conda/envs, from `conda list --json` (PyPI wheels show as channel pypi), the
        R packages built from source (<env>.r_source.tsv, installed and version-checked by
        install_envs.sh) and fastGEAR only when its files are in the image.

The hand-written "notes" of an env and the "fastgear" section are kept from the existing
locks summary. Every run records the image's "envs" tools as provenance
(pipeline.ContainerRuntime: {env: {tool: "version build"}}).
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import subprocess
from pathlib import Path

PLATFORMS = ("linux-64", "linux-aarch64")
DOCKER_PLATFORMS = {"x86_64": "linux/amd64", "aarch64": "linux/arm64"}
# The tools worth naming per env (the rest of each env is in its lock).
MAIN_TOOLS = {
    "bakta": ["bakta", "ncbi-amrfinderplus", "diamond", "blast", "blast-legacy", "pyrodigal",
              "pyhmmer", "hmmer", "infernal", "trnascan-se", "aragorn", "python"],
    "diamond": ["diamond"],
    "mobsuite": ["mob_suite", "blast", "python"],
    "phylogeny": ["snippy", "raxml-ng", "vcflib", "minimap2", "bwa", "samtools", "bcftools", "freebayes",
                  "snpeff", "python"],
    "gubbins": ["gubbins", "iqtree", "raxml", "raxml-ng", "fasttree", "python"],
    "panaroo": ["panaroo", "intbitset", "mafft", "prank", "cd-hit", "mash", "python"],
    "mafft": ["mafft"],
    "chewbbaca": ["chewbbaca", "blast", "pyrodigal", "mafft", "fasttree", "python"],
    "mlst": ["mlst", "blast", "any2fasta", "perl"],
    "quast": ["quast", "minimap2", "python"],
    "checkm2": ["checkm2", "diamond", "prodigal", "tensorflow", "python"],
    "fastbaps": ["fastbaps", "genie", "r-base", "r-matrix", "r-ape", "r-rcpp", "r-rcpparmadillo"],
    "base": ["conda", "mamba", "python"],  # the base image's own env (not from a lock)
}


def r_sources(path: Path) -> list[dict]:
    out = []
    if path.is_file():
        for line in path.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) >= 4:
                out.append({"name": parts[0], "version": parts[1], "build": f"source sha256:{parts[2]}",
                            "channel": parts[3]})
    return out


def tools_of(env: str, records: dict[str, dict]) -> dict[str, dict]:
    names = MAIN_TOOLS.get(env) or sorted(records)
    return {n: records[n] for n in names if n in records}


def from_locks(locks: Path, old: dict) -> dict:
    old_envs = old.get("envs") if isinstance(old.get("envs"), dict) else {}
    envs: dict[str, dict] = {}
    names = sorted({p.stem for plat in PLATFORMS for p in (locks / plat).glob("*.json")})
    for env in names:
        per_plat = {}
        for plat in PLATFORMS:
            path = locks / plat / f"{env}.json"
            if path.is_file():
                recs = {p["name"]: p for p in json.loads(path.read_text())}
                recs.update({r["name"]: r for r in r_sources(locks / plat / f"{env}.r_source.tsv")})
                pip = locks / plat / f"{env}.pip.txt"
                for line in pip.read_text().splitlines() if pip.is_file() else []:
                    name, _, rest = line.partition("==")
                    recs[name] = {"name": name, "version": rest.split()[0], "build": "wheel", "channel": "pypi"}
                per_plat[plat] = recs
        first = next(iter(per_plat.values()))
        tools = {}
        for name, rec in tools_of(env, first).items():
            tools[name] = {"version": rec["version"], "channel": rec["channel"],
                           "build": {plat: recs[name]["build"] for plat, recs in per_plat.items() if name in recs}}
        entry = {"lock": f"{env}.lock", "packages_total": {plat: len(r) for plat, r in per_plat.items()},
                 "tools": tools}
        notes = (old_envs.get(env) or {}).get("notes")
        if notes:
            entry["notes"] = notes
        envs[env] = entry
    doc = {
        "locks": {
            "solved": dt.date.today().isoformat(),
            "how": "docker/update_locks.sh: docker/solve_locks.py in a throwaway condaforge/miniforge3 "
                   "container solves yamls/<env>.yaml for linux-64 and linux-aarch64 (mamba --dry-run "
                   "--platform, no emulation) with the same package versions on both; docker/check_locks.py "
                   "checks they match",
            "platforms": list(PLATFORMS),
            "note": "The versions the locks name. A built image records what it really holds in its own "
                    "/opt/bactolomew/tool_versions.json (written by docker/verify_image.sh).",
        },
    }
    if "fastgear" in old:
        doc["fastgear"] = old["fastgear"]
    doc["envs"] = envs
    return doc


def from_image(summary: dict) -> dict:
    sum_envs = summary.get("envs") if isinstance(summary.get("envs"), dict) else {}
    envs: dict[str, dict] = {}
    for prefix in [Path("/opt/conda"), *sorted(Path("/opt/conda/envs").iterdir())]:
        env = "base" if prefix == Path("/opt/conda") else prefix.name
        if not (prefix / "conda-meta").is_dir():
            continue
        res = subprocess.run(["/opt/conda/bin/conda", "list", "--json", "-n", env],
                             capture_output=True, text=True, check=True)
        recs = {p["name"]: {"name": p["name"], "version": p["version"], "build": p.get("build_string", ""),
                            "channel": p.get("channel", "")} for p in json.loads(res.stdout)}
        recs.update({r["name"]: r for r in r_sources(Path("/opt/bactolomew/locks") / f"{env}.r_source.tsv")})
        tools = {n: {k: r[k] for k in ("version", "build", "channel")} for n, r in tools_of(env, recs).items()}
        entry = {"lock": f"{env}.lock" if env != "base" else None, "packages_total": len(recs), "tools": tools}
        notes = (sum_envs.get(env) or {}).get("notes")
        if notes:
            entry["notes"] = notes
        envs[env] = entry
    doc = {
        "platform": DOCKER_PLATFORMS.get(os.uname().machine, os.uname().machine),
        "written": "at build time by docker/verify_image.sh from the installed envs",
        "locks": summary.get("locks"),
    }
    if Path("/opt/fastGEAR").is_file() and "fastgear" in summary:
        doc["fastgear"] = summary["fastgear"]
    doc["envs"] = envs
    return doc


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("locks", "image"))
    ap.add_argument("source", type=Path, help="locks: the docker/locks folder; image: the locks summary")
    ap.add_argument("dest", type=Path)
    args = ap.parse_args()
    if args.mode == "locks":
        old = json.loads(args.dest.read_text()) if args.dest.is_file() else {}
        doc = from_locks(args.source, old)
    else:
        doc = from_image(json.loads(args.source.read_text()))
    args.dest.write_text(json.dumps(doc, indent=2) + "\n")


if __name__ == "__main__":
    main()

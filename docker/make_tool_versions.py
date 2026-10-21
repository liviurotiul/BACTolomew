"""Write docker/locks/tool_versions.json from `conda list --json` exports of a built image.

Usage (called by docker/update_locks.sh):
    python3 docker/make_tool_versions.py <dir with <env>.json> <tool_versions.json>
        --image IMAGE --image-id ID --image-created TIME

The "envs" section is rebuilt from the exports (the main tools of each env with version,
build and channel); envs that were not exported (an image from docker/Dockerfile.solve
holds only the envs being changed) keep their existing entry. Hand-written parts of an
existing file are kept: each env's "notes", "build_provenance" and "fastgear". The image copies this file to
/opt/bactolomew/tool_versions.json and every run records its "envs" tools as provenance.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path

# The tools worth naming per env (the rest of each env is in its lock).
MAIN_TOOLS = {
    "bakta": ["bakta", "ncbi-amrfinderplus", "diamond", "blast", "blast-legacy", "pyrodigal",
              "pyhmmer", "hmmer", "infernal", "trnascan-se", "aragorn", "python"],
    "diamond": ["diamond"],
    "mobsuite": ["mob_suite", "blast", "python"],
    "phylogeny": ["snippy", "raxml-ng", "minimap2", "bwa", "samtools", "bcftools", "freebayes", "python"],
    "gubbins": ["gubbins", "iqtree", "raxml", "raxml-ng", "fasttree", "python"],
    "panaroo": ["panaroo", "mafft", "prank", "cd-hit", "diamond", "prodigal", "python"],
    "mafft": ["mafft"],
    "chewbbaca": ["chewbbaca", "blast", "pyrodigal", "mafft", "fasttree", "python"],
    "mlst": ["mlst", "blast", "any2fasta", "python"],
    "assembly_qc": ["quast", "checkm2", "diamond", "prodigal", "tensorflow", "python"],
    "base": ["conda", "mamba", "python"],
}


LOCKS_NOTE = (
    "The image the locks were exported from, not necessarily the image that carries this file: "
    "an image built by docker/Dockerfile installs its envs from these locks, so it holds the "
    "same packages. The carrying image's own identity is in its OCI labels."
)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("exports", type=Path)
    ap.add_argument("dest", type=Path)
    ap.add_argument("--image", required=True)
    ap.add_argument("--image-id", required=True)
    ap.add_argument("--image-created", default="")
    args = ap.parse_args()

    old = json.loads(args.dest.read_text()) if args.dest.is_file() else {}
    old_envs = old.get("envs") if isinstance(old.get("envs"), dict) else {}
    envs: dict[str, dict] = dict(old_envs)
    for export in sorted(args.exports.glob("*.json")):
        env = export.stem
        records = {p["name"]: p for p in json.loads(export.read_text())}
        names = MAIN_TOOLS.get(env) or sorted(records)
        tools = {
            n: {"version": records[n]["version"], "build": records[n]["build_string"],
                "channel": records[n]["channel"]}
            for n in names if n in records
        }
        entry = {"lock": f"{env}.lock", "packages_total": len(records), "tools": tools}
        notes = (old_envs.get(env) or {}).get("notes")
        if notes:
            entry["notes"] = notes
        envs[env] = entry

    doc = {
        "locks_exported_from": {
            "image": args.image,
            "image_id": args.image_id,
            "image_created": args.image_created,
            "exported": f"{dt.date.today().isoformat()}, conda list --explicit --md5 inside the image "
                        "(docker run --rm --network none)",
            "note": LOCKS_NOTE,
        },
        "platform": "linux/amd64",
    }
    for key in ("build_provenance", "fastgear"):
        if key in old:
            doc[key] = old[key]
    doc["envs"] = dict(sorted(envs.items()))
    args.dest.write_text(json.dumps(doc, indent=2) + "\n")


if __name__ == "__main__":
    main()

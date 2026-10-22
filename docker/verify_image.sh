#!/usr/bin/env bash
# Build-time checks and records, run in a layer after the env layers (so adding or changing
# a check never rebuilds an env). The build fails when a check fails.
#
#  1. The population structure methods: FastBAPS (the fastbaps env loads the R package),
#     required on every platform as the default method, and fastGEAR (its launcher, binary
#     and MATLAB runtime under /opt) when the fastgear stage put it there. With fastGEAR,
#     the system libraries its MATLAB runtime needs, which the apt layer installs with
#     fallbacks, must be present under the names Ubuntu 24.04 uses.
#  2. The mob_suite database that the mob_suite package downloads at install time
#     (mob_init: Zenodo record 10304948 over HTTPS without certificate checks, then the
#     NCBI taxonomy through ete3) must match docker/locks/mobsuite_db.sha256: the files
#     of the Zenodo archive and the mash sketch built from them. The taxonomy
#     (taxa.sqlite) is whatever NCBI served on the build date; its checksum is recorded.
#  3. Writes /opt/bactolomew/image_contents.json (mob_suite database date and checksums,
#     "population_methods": the methods present, fastgear first),
#     /opt/bactolomew/apt_packages.tsv (every system package and version) and
#     /opt/bactolomew/tool_versions.json (docker/make_tool_versions.py: the tools of the
#     envs really installed, fastGEAR only when present), which the run provenance reads.
set -euo pipefail
OUT=/opt/bactolomew

methods=()
if [ -x /opt/run_fastGEAR.sh ] && [ -f /opt/fastGEAR ] && [ -d /opt/v901 ]; then
    dpkg -s libglu1-mesa libasound2t64 > /dev/null
    methods+=(fastgear)
fi
if [ -d /opt/conda/envs/fastbaps ] \
    && conda run -n fastbaps Rscript -e 'suppressPackageStartupMessages(library(fastbaps))' > /dev/null; then
    methods+=(fastbaps)
fi
[[ " ${methods[*]} " == *" fastbaps "* ]] || { echo "FastBAPS (the default population structure method) is not usable in the image" >&2; exit 1; }
dpkg-query -W -f='${Package}\t${Version}\n' > "$OUT/apt_packages.tsv"

db=""
for cand in /opt/conda/envs/mobsuite/lib/python3*/site-packages/mob_suite/databases; do
    [ -f "$cand/status.txt" ] && { db="$(realpath "$cand")"; break; }
done
[ -n "$db" ] || { echo "mob_suite database not found (status.txt missing)" >&2; exit 1; }
(cd "$db" && sha256sum --check --strict --quiet "$OUT/locks/mobsuite_db.sha256")

/opt/conda/bin/python - "$db" "$OUT/locks/mobsuite_db.sha256" "$OUT/image_contents.json" "${methods[@]}" <<'EOF'
import hashlib, json, re, sys
from pathlib import Path

db, manifest, dest = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
methods = sys.argv[4:]

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

status = (db / "status.txt").read_text().strip()
m = re.search(r"\d{4}-\d{2}-\d{2}", status)
verified = {}
for line in manifest.read_text().splitlines():
    if line.strip():
        digest, name = line.split(None, 1)
        verified[name.strip()] = digest
doc = {"mobsuite_db": {
    "path": str(db),
    "download_date": m.group(0) if m else status,
    "source": "https://zenodo.org/records/10304948/files/data.tar.gz (mob_init)",
    "verified_sha256": verified,
    "taxa_sqlite_sha256": sha256(db / "taxa.sqlite") if (db / "taxa.sqlite").is_file() else None,
}, "population_methods": methods}
dest.write_text(json.dumps(doc, indent=2) + "\n")
EOF
/opt/conda/bin/python "$OUT/make_tool_versions.py" image "$OUT/locks/tool_versions.json" "$OUT/tool_versions.json"
echo "[verify_image] ok: population methods ${methods[*]}, mob_suite database matches the manifest"

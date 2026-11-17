#!/usr/bin/env bash
# Build vcflib 1.0.9 from source for a baseline CPU, inside the Dockerfile's vcflib stage.
# Usage (inside the stage): build_vcflib.sh <out-dir>
#
# Why: the bioconda vcflib 1.0.9 package (linux-64 build hdcf5f25_8) links a zig static
# library (src/zig, vcfcreatemulti's new code) into libvcflib and every vcf* program. zig
# builds for the CPU of the machine it runs on (the bioconda builder has AVX-512), so its
# compiler runtime (memcpy, memmove, memset and others) uses AVX2, BMI2 and AVX-512
# without any CPU check: vcfstreamsort and vcfuniq (which Snippy runs) die with "Illegal
# instruction" in memset at start on any x86_64 CPU without AVX-512 (most desktop and
# laptop CPUs, every Intel Mac but the Xeon ones, many virtual machines).
# The package's own C++ code is already built without -march=native (its
# march_native.patch).
#
# What: the same source (the release tarball with its submodules, pinned by sha256), the
# same bioconda patches (docker/vcflib/patches/, applied in the recipe's order), the same
# CMake options (zig on, OpenMP on, the env's wfa2-lib, shared libvcflib, pyvcflib), the
# same zig 0.10.1 (pinned by sha256), plus:
#   * zig_baseline.patch: zig builds for the baseline CPU of the architecture (x86-64 v1,
#     generic aarch64) instead of the build machine's;
#   * C and C++ flags -O3 -march=x86-64 -mtune=generic on x86_64 (-march=armv8-a on
#     aarch64), no -march=native;
#   * the glibc 2.17 sysroot (as bioconda's build): the programs linked with conda-forge's
#     2.39 sysroot inherit its crt1.o's "x86 ISA needed: x86-64-v3" note, and glibc then
#     refuses to start them on older CPUs.
# The libraries it links (htslib 1.20, tabixpp 1.1.2, wfa2-lib 2.3.6, zlib, Python 3.12)
# are the phylogeny env's own, created here from the phylogeny lock at the path it has in
# the image; the compilers and build tools come from docker/vcflib/locks/<platform>/
# vcflib_build.lock (solved from vcflib_build.yaml by docker/solve_locks.py; never shipped).
#
# Output: <out-dir>/opt/conda/envs/phylogeny/... holding only the ELF files the bioconda
# vcflib package installs (bin/vcf* programs, lib/libvcflib.so.1, lib/pyvcflib*.so), at
# their paths, plus <out-dir>/opt/bactolomew/source_builds/phylogeny.json (the provenance
# that docker/make_tool_versions.py records). The build fails when a program of the
# package is not rebuilt, or when a rebuilt x86_64 file still holds a VEX (AVX, AVX2,
# FMA), BMI or AVX-512 instruction or is marked as needing more than the x86-64 baseline.
set -euo pipefail
OUT="${1:?usage: build_vcflib.sh <out-dir>}"
HERE=/opt/bactolomew/vcflib
LOCKS=/opt/bactolomew/locks
PREFIX=/opt/conda/envs/phylogeny
BUILD_ENV=/opt/conda/envs/vcflib_build
VCFLIB_VERSION=1.0.9
VCFLIB_URL=https://github.com/vcflib/vcflib/releases/download/v1.0.9/vcflib-1.0.9-src.tar.gz
VCFLIB_SHA256=493204cb8e0515defa1383d56280d00b8efbed55312c6055162f8cb3a4ffb210
ZIG_VERSION=0.10.1
case "$(uname -m)" in
    x86_64)
        ZIG_SHA256=6699f0e7293081b42428f32c9d9c983854094bd15fee5489f12c4cf4518cc380
        PLATFORM=linux-64
        ARCH_FLAGS="-march=x86-64 -mtune=generic" ;;
    aarch64)
        ZIG_SHA256=db0761664f5f22aa5bbd7442a1617dd696c076d5717ddefcc9d8b95278f71f5d
        PLATFORM=linux-aarch64
        ARCH_FLAGS="-march=armv8-a -mtune=generic" ;;
    *) echo "[vcflib] unsupported architecture $(uname -m)" >&2; exit 1 ;;
esac
ZIG_URL="https://ziglang.org/download/$ZIG_VERSION/zig-linux-$(uname -m)-$ZIG_VERSION.tar.xz"

echo "[vcflib] envs: phylogeny (the libraries) and vcflib_build (the tools)"
conda create --yes --quiet --name phylogeny --file "$LOCKS/phylogeny.lock"
conda create --yes --quiet --name vcflib_build --file "$HERE/locks/$PLATFORM/vcflib_build.lock"

WORK="$(mktemp -d)"
fetch() {  # fetch <url> <sha256> <dest-dir>: download, check, unpack
    /opt/conda/bin/python - "$@" <<'PY'
import hashlib, io, sys, tarfile, urllib.request
url, digest, dest = sys.argv[1:4]
data = urllib.request.urlopen(url, timeout=600).read()
got = hashlib.sha256(data).hexdigest()
if got != digest:
    sys.exit(f"[vcflib] {url}: sha256 {got}, expected {digest}")
with tarfile.open(fileobj=io.BytesIO(data)) as tar:
    tar.extractall(dest, filter="tar")
print(f"[vcflib] {url}: sha256 ok")
PY
}
fetch "$VCFLIB_URL" "$VCFLIB_SHA256" "$WORK"
fetch "$ZIG_URL" "$ZIG_SHA256" "$WORK"
SRC="$WORK/vcflib-$VCFLIB_VERSION"
export PATH="$WORK/zig-linux-$(uname -m)-$ZIG_VERSION:$BUILD_ENV/bin:$PATH"
export ZIG_GLOBAL_CACHE_DIR="$WORK/zig-cache" ZIG_LOCAL_CACHE_DIR="$WORK/zig-cache"

# The bioconda recipe's patches in its order, then ours.
for p in pkg-config shared_lib backport_pr_394 tabix march_native zig_baseline; do
    git -C "$SRC" apply -p1 --verbose "$HERE/patches/$p.patch"
done

export CC="$BUILD_ENV/bin/gcc" CXX="$BUILD_ENV/bin/g++"
export PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig"
FLAGS="-O3 $ARCH_FLAGS -D_FILE_OFFSET_BITS=64 -fPIC -I$PREFIX/include"
export CFLAGS="$FLAGS" CXXFLAGS="$FLAGS"
# The libraries as in the bioconda build.sh: libvcflib itself needs wfa2, htslib and tabixpp.
export LDFLAGS="-L$PREFIX/lib -Wl,-rpath-link,$PREFIX/lib -lhts -ltabixpp -lpthread -lz -lm -llzma -lbz2 -fopenmp -lwfa2"
cmake -S "$SRC" -B "$WORK/build" \
    -DZIG=ON -DOPENMP=ON -DWFA_GITMODULE=OFF -DBUILD_DOC=OFF \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_INSTALL_LIBDIR=lib \
    -DCMAKE_INSTALL_PREFIX="$PREFIX" \
    -DCMAKE_PREFIX_PATH="$PREFIX" \
    -DCMAKE_C_FLAGS="$CFLAGS" -DCMAKE_CXX_FLAGS="$CXXFLAGS" \
    -DCMAKE_INSTALL_RPATH='$ORIGIN/../lib' \
    -DPython_EXECUTABLE="$PREFIX/bin/python" -DPYTHON_EXECUTABLE="$PREFIX/bin/python"
cmake --build "$WORK/build" -j "$(nproc)"
DESTDIR="$WORK/stage" cmake --install "$WORK/build"

# Keep only the package's ELF files, at their paths, from the staged install; every ELF
# file of the package must have been rebuilt.
/opt/conda/bin/python - "$PREFIX" "$WORK/stage$PREFIX" "$OUT$PREFIX" <<'PY'
import glob, json, shutil, sys
from pathlib import Path
prefix, stage, out = (Path(a) for a in sys.argv[1:4])
meta = json.loads(Path(glob.glob(f"{prefix}/conda-meta/vcflib-*.json")[0]).read_text())
elf = [f for f in meta["files"] if (prefix / f).is_file() and (prefix / f).read_bytes()[:4] == b"\x7fELF"]
missing = [f for f in elf if not (stage / f).is_file()]
if missing:
    sys.exit(f"[vcflib] not rebuilt: {' '.join(missing)}")
for f in elf:
    (out / f).parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(stage / f, out / f)
print(f"[vcflib] {len(elf)} ELF files of the vcflib package rebuilt")
PY
# conda-forge's gcc driver adds its own env's lib folder to every RPATH: give the copies
# the RPATH the bioconda files have ($ORIGIN/../lib, the phylogeny env's lib).
find "$OUT$PREFIX" -type f -print0 | xargs -0 -I{} cmake -DF={} -P "$HERE/set_rpath.cmake"

if [ "$(uname -m)" = x86_64 ]; then
    # objdump from the build env's binutils: no VEX (v* instructions on xmm/ymm/zmm
    # registers), BMI1/BMI2 or AVX-512 mask register instruction may be left.
    OBJDUMP="$(ls "$BUILD_ENV"/bin/*-objdump | head -n1)"
    [ -x "$OBJDUMP" ] || { echo "[vcflib] no objdump in $BUILD_ENV/bin" >&2; exit 1; }
    bad="$(find "$OUT$PREFIX" -type f -print0 | xargs -0 -n1 bash -c \
        'asm="$("$0" -d --no-show-raw-insn "$1")" || exit 255
         n=$(grep -cE "\s(v[a-z0-9]+\s.*%[xyz]mm|(andn|bextr|blsi|blsmsk|blsr|bzhi|mulx|pdep|pext|rorx|sarx|shlx|shrx)\s|k[a-z]+[bwdq]?\s.*%k[0-7])" <<< "$asm")
         [ "$n" = 0 ] || echo "$1: $n"' "$OBJDUMP")"
    [ -z "$bad" ] || { echo "[vcflib] instructions above the x86-64 baseline in:" >&2; echo "$bad" >&2; exit 1; }
    # Nor a GNU property note that asks for more than the baseline: glibc refuses to load
    # such a file on a lower CPU ("CPU ISA level is lower than required").
    READELF="${OBJDUMP%objdump}readelf"
    bad="$(find "$OUT$PREFIX" -type f -print0 | xargs -0 -n1 bash -c \
        'notes="$("$0" -n "$1")" || exit 255
         ! grep -E "ISA needed:.*x86-64-v[234]" <<< "$notes" > /dev/null || echo "$1"' "$READELF")"
    [ -z "$bad" ] || { echo "[vcflib] marked as needing more than the x86-64 baseline:" >&2; echo "$bad" >&2; exit 1; }
    echo "[vcflib] no AVX, BMI or AVX-512 instruction and no ISA level above the baseline in the rebuilt files"
fi

# The programs Snippy runs, on a tiny VCF: sorted, the duplicate record dropped.
export LD_LIBRARY_PATH="$OUT$PREFIX/lib:$PREFIX/lib"
got="$(printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nc\t5\t.\tA\tG\t9\t.\t.\nc\t2\t.\tC\tT\t9\t.\t.\nc\t2\t.\tC\tT\t9\t.\t.\n' \
    | "$OUT$PREFIX/bin/vcfstreamsort" | "$OUT$PREFIX/bin/vcfuniq" | grep -v '^#' | cut -f2 | tr '\n' ' ')"
unset LD_LIBRARY_PATH
[ "$got" = "2 5 " ] || { echo "[vcflib] vcfstreamsort | vcfuniq gave positions '$got', expected '2 5 '" >&2; exit 1; }
mkdir -p "$OUT/opt/bactolomew/source_builds"
cat > "$OUT/opt/bactolomew/source_builds/phylogeny.json" <<EOF
{
  "vcflib": {
    "version": "$VCFLIB_VERSION",
    "build": "source, $ARCH_FLAGS, zig $ZIG_VERSION baseline CPU",
    "channel": "$VCFLIB_URL",
    "sha256": "$VCFLIB_SHA256",
    "replaces": "the ELF files of the bioconda package (conda-meta keeps its record); docker/vcflib/build_vcflib.sh"
  }
}
EOF
rm -rf "$WORK"
echo "[vcflib] done: $OUT"

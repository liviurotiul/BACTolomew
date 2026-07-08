#!/usr/bin/env bash
# BACTolomew host installer.
#
# Sets up the THIN host environment (PyQt6/WebEngine GUI + pipeline orchestration).
# All bioinformatics tools run in a container (Docker/Singularity) that the app builds
# or pulls on first launch -- this script does NOT install any bioinformatics tools.
#
# Usage:
#   # From a clone (recommended for development):
#   git clone https://github.com/liviurotiul/BACTolomew.git && cd BACTolomew && bash install.sh
#
#   # Bootstrap without cloning first (installs into ~/BACTolomew by default):
#   curl -fsSL https://raw.githubusercontent.com/liviurotiul/BACTolomew/main/install.sh | bash
#
#   # Offline / custom source or target:
#   bash install.sh [--source <repo-dir|tarball>] [--dir <install-dir>] [--ref <branch|tag>]
#
# When run outside a checkout and without --source, the installer clones GITHUB_URL.
set -euo pipefail

say() { printf '\033[1;34m[install]\033[0m %s\n' "$*"; }
die() { printf '\033[1;31m[install] ERROR:\033[0m %s\n' "$*" >&2; exit 1; }

GITHUB_URL="https://github.com/liviurotiul/BACTolomew.git"
GITHUB_REF=""            # branch/tag to clone (empty = repo default branch)

# Self-contained help so it works even over `curl ... | bash -s -- --help` (no file on disk).
usage() {
  cat <<EOF
BACTolomew host installer.

Usage:
  # From a clone (recommended for development):
  git clone $GITHUB_URL && cd BACTolomew && bash install.sh

  # Bootstrap without cloning first (installs into ~/BACTolomew by default):
  curl -fsSL https://raw.githubusercontent.com/liviurotiul/BACTolomew/main/install.sh | bash

  # Offline / custom source or target:
  bash install.sh [--source <repo-dir|tarball>] [--dir <install-dir>] [--ref <branch|tag>]

Outside a checkout and without --source, the installer clones the repo.
EOF
}

SELF="${BASH_SOURCE[0]:-$0}"   # this script's path ("bash" when piped via curl|bash)
SCRIPT_DIR="$(cd "$(dirname "$SELF")" 2>/dev/null && pwd || echo "$PWD")"
SOURCE=""
INSTALL_DIR=""

need_val() { [ "$#" -ge 2 ] || die "$1 requires a value"; }

while [ $# -gt 0 ]; do
  case "$1" in
    --source)       need_val "$@"; SOURCE="$2"; shift 2 ;;
    --dir)          need_val "$@"; INSTALL_DIR="$2"; shift 2 ;;
    --ref|--branch) need_val "$@"; GITHUB_REF="$2"; shift 2 ;;
    -h|--help)      usage; exit 0 ;;
    *) die "unknown argument: $1 (see --help)" ;;
  esac
done

# --- 1. prerequisites -------------------------------------------------------
command -v python3 >/dev/null || die "python3 not found (need >= 3.11)."
python3 -c 'import sys; sys.exit(0 if sys.version_info >= (3,11) else 1)' \
  || die "Python >= 3.11 required (found $(python3 -V 2>&1))."

if command -v docker >/dev/null; then RUNTIME_FOUND="docker"
elif command -v singularity >/dev/null || command -v apptainer >/dev/null; then RUNTIME_FOUND="singularity"
else die "Neither Docker nor Singularity/Apptainer found. Install one:
  Docker:      https://docs.docker.com/engine/install/
  Apptainer:   https://apptainer.org/docs/admin/main/installation.html"
fi
say "Container runtime available: $RUNTIME_FOUND"

# --- 2. obtain the code -----------------------------------------------------
# Precedence: explicit --source  >  running inside a checkout  >  clone from GitHub.
have_checkout() { [ -f "$1/src/annotation_explorer_polars/app.py" ]; }
dir_empty()     { [ ! -e "$1" ] || [ -z "$(ls -A "$1" 2>/dev/null)" ]; }

# Copy a repo directory or tarball into <dest>. For a git checkout, copy only tracked
# files (git archive) so gitignored data -- envs/ DB/ runs/ (tens of GB) -- is skipped.
copy_tree() {  # copy_tree <src> <dest>
  local src="$1" dest="$2"
  mkdir -p "$dest" || die "cannot create $dest"
  case "$src" in
    *.tar.gz|*.tgz|*.tar)
      say "Extracting $src -> $dest"
      tar -xf "$src" -C "$dest" --strip-components=1 || die "failed to extract $src" ;;
    *)
      [ -d "$src" ] || die "--source '$src' is not a directory or a .tar.gz/.tgz/.tar tarball."
      if [ -d "$src/.git" ] && command -v git >/dev/null; then
        say "Copying tracked files $src -> $dest (git archive; skips gitignored data)"
        git -C "$src" archive --format=tar HEAD | tar -x -C "$dest" || die "failed to copy from $src"
      else
        say "Copying $src -> $dest"
        cp -a "$src/." "$dest/" || die "failed to copy from $src"
      fi ;;
  esac
}

# Determine the code source.
if [ -n "$SOURCE" ]; then
  SRC="$SOURCE"; SRC_KIND="source"
elif have_checkout "$SCRIPT_DIR"; then
  SRC="$SCRIPT_DIR"; SRC_KIND="inplace"
else
  SRC=""; SRC_KIND="clone"
fi

# Determine the install dir (default: the checkout itself for in-place; else ~/BACTolomew).
if [ -z "$INSTALL_DIR" ]; then
  if [ "$SRC_KIND" = "inplace" ]; then INSTALL_DIR="$SRC"
  else INSTALL_DIR="${HOME:-$PWD}/BACTolomew"; fi
fi
mkdir -p "$INSTALL_DIR" || die "cannot create install dir: $INSTALL_DIR"
INSTALL_DIR="$(cd "$INSTALL_DIR" && pwd)" || die "cannot access install dir: $INSTALL_DIR"

if [ -n "$GITHUB_REF" ] && [ "$SRC_KIND" != "clone" ]; then
  say "note: --ref '$GITHUB_REF' only applies when cloning; ignoring for $SRC_KIND install."
fi

# Never write into a non-empty directory that is not already a BACTolomew checkout.
if ! dir_empty "$INSTALL_DIR" && ! have_checkout "$INSTALL_DIR"; then
  die "Target $INSTALL_DIR is not empty and is not a BACTolomew checkout.
  Choose an empty --dir <path>, or clear it first."
fi

case "$SRC_KIND" in
  inplace)
    if [ "$SRC" = "$INSTALL_DIR" ]; then say "Installing in-place at $INSTALL_DIR"
    else copy_tree "$SRC" "$INSTALL_DIR"; fi ;;
  source)
    copy_tree "$SRC" "$INSTALL_DIR" ;;
  clone)
    command -v git >/dev/null || die "git not found (needed to clone $GITHUB_URL)."
    if have_checkout "$INSTALL_DIR"; then
      say "Reusing existing checkout at $INSTALL_DIR"
    else
      say "Cloning $GITHUB_URL${GITHUB_REF:+ (ref: $GITHUB_REF)} -> $INSTALL_DIR"
      # Remove a partial clone if the clone fails or is interrupted.
      CLONE_TGT="$INSTALL_DIR"
      trap 'rc=$?; [ -n "${CLONE_TGT:-}" ] && ! have_checkout "$CLONE_TGT" && rm -rf "$CLONE_TGT"; exit $rc' EXIT INT TERM
      if [ -n "$GITHUB_REF" ]; then
        git clone --branch "$GITHUB_REF" "$GITHUB_URL" "$INSTALL_DIR" || die "clone failed (ref: $GITHUB_REF)"
      else
        git clone "$GITHUB_URL" "$INSTALL_DIR" || die "clone failed"
      fi
      CLONE_TGT=""; trap - EXIT INT TERM
    fi ;;
esac

have_checkout "$INSTALL_DIR" \
  || die "Could not find $INSTALL_DIR/src/annotation_explorer_polars/app.py after obtaining the code."
say "Install dir: $INSTALL_DIR"

# --- 3. venv + deps ---------------------------------------------------------
VENV="$INSTALL_DIR/.venv"
say "Creating venv at $VENV"
python3 -m venv "$VENV"
"$VENV/bin/pip" install --upgrade pip >/dev/null
say "Installing host Python dependencies (this can take a minute)..."
"$VENV/bin/pip" install -r "$INSTALL_DIR/requirements.txt"

# --- 4. WebEngine system libraries ------------------------------------------
missing=""
for lib in libnss3.so libxkbcommon.so libGL.so.1 libnspr4.so; do
  ldconfig -p 2>/dev/null | grep -q "$lib" || missing="$missing $lib"
done
if [ -n "$missing" ]; then
  say "WARNING: WebEngine may need system libraries not found:$missing"
  cat <<'HINT'
  If the GUI fails to render, install them (one of):
    Debian/Ubuntu: sudo apt install libnss3 libnspr4 libxkbcommon0 libxcomposite1 \
                     libxdamage1 libxrandr2 libgbm1 libasound2 libxshmfence1 fontconfig libgl1
    Fedora/RHEL:   sudo dnf install nss nspr libxkbcommon libXcomposite libXdamage \
                     libXrandr mesa-libgbm alsa-lib libxshmfence fontconfig mesa-libGL
HINT
fi

# --- 5. launcher ------------------------------------------------------------
LAUNCHER="$INSTALL_DIR/bactolomew"
cat > "$LAUNCHER" <<EOF
#!/usr/bin/env bash
cd "$INSTALL_DIR/src"
exec "$VENV/bin/python" -m annotation_explorer_polars.app "\$@"
EOF
chmod +x "$LAUNCHER"

say "Done. Launch with:  $LAUNCHER"
say "On first launch, BACTolomew will select a container runtime and build/pull the tools image."

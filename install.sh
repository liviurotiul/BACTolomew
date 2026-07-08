#!/usr/bin/env bash
# BACTolomew host installer.
#
# Sets up the THIN host environment (PyQt6/WebEngine GUI + pipeline orchestration).
# All bioinformatics tools run in a container (Docker/Singularity) that the app builds
# or pulls on first launch -- this script does NOT install any bioinformatics tools.
#
# Usage:
#   bash install.sh [--source <repo-dir|tarball>] [--dir <install-dir>]
#
# Today (repo not yet on GitHub): run from inside a checkout, or pass --source <dir>.
# Later (once published): the downloaded bootstrap will `git clone` GITHUB_URL first.
set -euo pipefail

GITHUB_URL="https://github.com/liviurotiul/BACTolomew.git"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE=""
INSTALL_DIR=""

while [ $# -gt 0 ]; do
  case "$1" in
    --source) SOURCE="$2"; shift 2 ;;
    --dir)    INSTALL_DIR="$2"; shift 2 ;;
    -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done

say() { printf '\033[1;34m[install]\033[0m %s\n' "$*"; }
die() { printf '\033[1;31m[install] ERROR:\033[0m %s\n' "$*" >&2; exit 1; }

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
resolve_source() {
  if [ -n "$SOURCE" ]; then echo "$SOURCE"; return; fi
  if [ -f "$SCRIPT_DIR/src/annotation_explorer_polars/app.py" ]; then echo "$SCRIPT_DIR"; return; fi
  die "No --source given and this script is not inside a checkout.
  Pass --source <repo-dir-or-tarball> (the repo is not on GitHub yet)."
}
SRC="$(resolve_source)"

if [ -z "$INSTALL_DIR" ]; then
  if [ -f "$SRC/src/annotation_explorer_polars/app.py" ]; then INSTALL_DIR="$SRC"   # in-place
  else INSTALL_DIR="$HOME/BACTolomew"; fi
fi

if [ ! -f "$INSTALL_DIR/src/annotation_explorer_polars/app.py" ]; then
  say "Copying code from $SRC -> $INSTALL_DIR"
  mkdir -p "$INSTALL_DIR"
  case "$SRC" in
    *.tar.gz|*.tgz) tar -xzf "$SRC" -C "$INSTALL_DIR" --strip-components=1 ;;
    *) cp -a "$SRC/." "$INSTALL_DIR/" ;;
  esac
fi
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

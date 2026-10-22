#!/usr/bin/env bash
# BACTolomew host installer.
#
# Sets up the THIN host environment (PyQt6/WebEngine GUI + pipeline orchestration).
# All bioinformatics tools run in a container (Docker/Podman/Singularity) that setup builds
# or pulls (`./bactolomew setup`, or the first-time setup window): this script does NOT
# install any bioinformatics tools.
#
# Usage:
#   # From a clone (recommended for development):
#   git clone https://github.com/liviurotiul/BACTolomew.git && cd BACTolomew && bash install.sh
#
#   # Bootstrap without cloning first (installs into ~/BACTolomew by default):
#   curl -fsSL -o install.sh https://raw.githubusercontent.com/liviurotiul/BACTolomew/main/install.sh && bash install.sh
#
#   # Custom source, target, Python, a headless (no Qt) or a viewer-only install:
#   bash install.sh [--source <repo-dir|tarball>] [--dir <install-dir>] [--ref <branch|tag>]
#                   [--python <path>] [--headless | --viewer] [--no-lock]
#
# When run outside a checkout and without --source, the installer clones GITHUB_URL.
# Everything runs from main() on the last line, so a partial download does nothing.
set -euo pipefail

say()  { printf '\033[1;34m[install]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[install] WARNING:\033[0m %s\n' "$*" >&2; }
die()  { printf '\033[1;31m[install] ERROR:\033[0m %s\n' "$*" >&2; exit 1; }

GITHUB_URL="https://github.com/liviurotiul/BACTolomew.git"
GITHUB_REF=""            # branch/tag to clone (empty = repo default branch)

# Supported Python (requirements.lock is one universal lock for all of them).
PY_MIN=10
PY_MAX=14
PY_CANDIDATES="python3.13 python3.12 python3.11 python3.10 python3.14 python3"
UV_PYTHON="3.12"
# PyQt6 / PyQt6-WebEngine 6.10+ wheels are manylinux_2_34.
GLIBC_GUI_MIN="2.34"
# The locked pipeline wheels (pyarrow, numpy, scipy, pandas) are manylinux_2_28; biopython has no musllinux wheel.
GLIBC_HEADLESS_MIN="2.28"

# Self-contained help so it works even over `curl ... | bash -s -- --help` (no file on disk).
usage() {
  cat <<EOF
BACTolomew host installer.

Usage:
  # From a clone (recommended for development):
  git clone $GITHUB_URL && cd BACTolomew && bash install.sh

  # Bootstrap without cloning first (installs into ~/BACTolomew by default):
  curl -fsSL -o install.sh https://raw.githubusercontent.com/liviurotiul/BACTolomew/main/install.sh && bash install.sh

Options:
  --source <dir|tarball>  install from a local checkout (tracked and untracked, not ignored files)
                          or a .tar.gz/.tgz/.tar
  --dir <install-dir>     where to install (default: the checkout itself, else ~/BACTolomew)
  --ref <branch|tag>      what to clone (only when cloning)
  --python <path>         the Python to build the venv with (also: PYTHON=<path>); Python 3.$PY_MIN to 3.$PY_MAX
  --headless              pipeline only, no Qt: ./bactolomew then offers setup, run, check and db
  --viewer                the window only, to open runs shared by others: no container runtime,
                          tools image or databases (about 1 GB); \`./bactolomew setup\` adds them later
  --no-lock               install the newest versions within requirements*.txt instead of the tested lock

Outside a checkout and without --source, the installer clones the repo.
EOF
}

# --------------------------------------------------------------------------- #
# The distro, for package hints
# --------------------------------------------------------------------------- #
OS_ID=""; OS_LIKE=""; OS_VERSION=""; PKG=""
detect_distro() {
  if [ -r /etc/os-release ]; then
    OS_ID="$(. /etc/os-release && echo "${ID:-}")"
    OS_LIKE="$(. /etc/os-release && echo "${ID_LIKE:-}")"
    OS_VERSION="$(. /etc/os-release && echo "${VERSION_ID:-}")"
  fi
  case " $OS_ID $OS_LIKE " in
    *" debian "*|*" ubuntu "*) PKG="apt" ;;
    *" fedora "*|*" rhel "*|*" centos "*|*" rocky "*|*" almalinux "*) PKG="dnf" ;;
    *suse*) PKG="zypper" ;;
    *" arch "*) PKG="pacman" ;;
    *) PKG="" ;;
  esac
}

# version_ge A B: A >= B for dotted versions.
version_ge() { [ "$(printf '%s\n%s\n' "$2" "$1" | sort -V | head -n1)" = "$2" ]; }

# Ubuntu 24.04+ and Debian 13+ renamed some libraries for the 64-bit time_t transition.
apt_t64() {
  case "$OS_ID" in
    ubuntu) version_ge "${OS_VERSION:-0}" "24.04" ;;
    debian) version_ge "${OS_VERSION:-0}" "13" ;;
    *) case " $OS_LIKE " in *" ubuntu "*) version_ge "${OS_VERSION:-0}" "24" ;; *) return 1 ;; esac ;;
  esac
}

pkg_install_cmd() {  # pkg_install_cmd <packages...>
  case "$PKG" in
    apt)    echo "sudo apt install $*" ;;
    dnf)    echo "sudo dnf install $*" ;;
    zypper) echo "sudo zypper install $*" ;;
    pacman) echo "sudo pacman -S --needed $*" ;;
    *)      echo "install: $*" ;;
  esac
}

# --------------------------------------------------------------------------- #
# glibc
# --------------------------------------------------------------------------- #
glibc_version() {
  local v=""
  v="$(getconf GNU_LIBC_VERSION 2>/dev/null | awk '{print $2}')" || true
  if [ -z "$v" ] && command -v ldd >/dev/null 2>&1; then
    v="$(ldd --version 2>&1 | head -n1 | grep -i -E 'glibc|gnu libc' | grep -oE '[0-9]+\.[0-9]+' | tail -n1)" || true
  fi
  echo "$v"
}

check_platform() {
  [ "$(uname -s)" = "Linux" ] || warn "BACTolomew is tested on Linux only (this is $(uname -s))."
  [ "$(uname -m)" = "x86_64" ] || warn "BACTolomew is tested on x86_64 only (this is $(uname -m))."
  local glibc; glibc="$(glibc_version)"
  if [ "$HEADLESS" = 1 ]; then
    if [ -z "$glibc" ]; then
      die "This system does not use glibc (musl, as on Alpine?), and the pipeline's Python packages need glibc $GLIBC_HEADLESS_MIN or newer."
    fi
    if ! version_ge "$glibc" "$GLIBC_HEADLESS_MIN"; then
      die "This OS is too old for BACTolomew: it has glibc $glibc, and the pipeline's Python packages (pyarrow, numpy,
  scipy, pandas) need glibc $GLIBC_HEADLESS_MIN or newer, even for a headless install (for example RHEL/Rocky 8 or newer)."
    fi
    say "C library: glibc $glibc (headless install, no Qt needed)"
    return
  fi
  if [ -z "$glibc" ]; then
    die "This system does not use glibc, and the BACTolomew window (PyQt6 wheels) needs glibc $GLIBC_GUI_MIN or newer.
  The pipeline works without the window: rerun with --headless (setup, run and check)."
  fi
  if ! version_ge "$glibc" "$GLIBC_GUI_MIN"; then
    die "This OS is too old for the BACTolomew window: it has glibc $glibc, and PyQt6 needs glibc $GLIBC_GUI_MIN or newer
  (for example Ubuntu 22.04, Debian 12, Rocky/RHEL 9, Fedora 35, openSUSE Leap 15.6 or newer).
  The pipeline works on this OS: rerun with --headless (setup, run and check), and view the
  results on a newer machine."
  fi
  say "C library: glibc $glibc"
}

# --------------------------------------------------------------------------- #
# Python
# --------------------------------------------------------------------------- #
PY=""            # the chosen interpreter
USE_UV=0         # 1: build the venv with `uv venv --python $UV_PYTHON`
PY_SEEN=""       # what was found, for the error message
PY_NOVENV=""     # supported interpreters without a working venv/ensurepip

# py_version <python>: prints "3.X" and succeeds when it is in the supported range.
py_version() {
  "$1" -c "import sys; print('%d.%d' % sys.version_info[:2]); sys.exit(0 if sys.version_info[0] == 3 and $PY_MIN <= sys.version_info[1] <= $PY_MAX else 1)" 2>/dev/null
}

# venv_works <python>: creates a throwaway venv (with pip), the only reliable test for
# Debian/Ubuntu's split-out python3-venv / ensurepip.
venv_works() {
  local tmp rc=0
  tmp="$(mktemp -d "${TMPDIR:-/tmp}/bactolomew-venv-test.XXXXXX")" || return 1
  "$1" -m venv "$tmp/v" >/dev/null 2>&1 && [ -x "$tmp/v/bin/python" ] && "$tmp/v/bin/python" -m pip --version >/dev/null 2>&1 || rc=1
  rm -rf "$tmp"
  return $rc
}

python_hint() {
  case "$PKG" in
    apt)
      if [ "$OS_ID" = "ubuntu" ] && ! version_ge "${OS_VERSION:-0}" "22.04"; then
        echo "  Ubuntu $OS_VERSION: sudo add-apt-repository ppa:deadsnakes/ppa && sudo apt install python3.12 python3.12-venv"
      else
        echo "  $(pkg_install_cmd python3 python3-venv)"
      fi ;;
    dnf)    if [ "$OS_ID" = "fedora" ]; then echo "  $(pkg_install_cmd python3)"
            else echo "  $(pkg_install_cmd python3.12)"; fi ;;
    zypper) echo "  $(pkg_install_cmd python312)" ;;
    pacman) echo "  $(pkg_install_cmd python)" ;;
    *)      echo "  Install Python 3.$PY_MIN to 3.$PY_MAX with your package manager." ;;
  esac
  echo "  Then rerun install.sh (or point it at one: --python /path/to/python3.12)."
  echo "  With uv (https://docs.astral.sh/uv/) on PATH, install.sh builds the venv with \`uv venv --python $UV_PYTHON\` instead."
}

venv_hint() {  # venv_hint <python> <3.X>
  case "$PKG" in
    apt)    echo "  $(pkg_install_cmd "python$2-venv")" ;;
    dnf)    echo "  $(pkg_install_cmd "python$2" "python$2-pip")" ;;
    zypper) echo "  $(pkg_install_cmd "python${2/./}" "python${2/./}-pip")" ;;
    pacman) echo "  $(pkg_install_cmd python python-pip)" ;;
    *)      echo "  Install the venv/ensurepip module for $1." ;;
  esac
}

find_python() {
  local requested="${PYTHON_ARG:-${PYTHON:-}}" c path v
  if [ -n "$requested" ]; then
    path="$(command -v "$requested" 2>/dev/null)" || die "--python $requested: not found."
    v="$(py_version "$path")" || die "$path is Python ${v:-?}; BACTolomew needs Python 3.$PY_MIN to 3.$PY_MAX."
    venv_works "$path" || die "$path (Python $v) cannot create a virtual environment (venv/ensurepip missing). Install it:
$(venv_hint "$path" "$v")"
    PY="$path"; say "Python: $PY ($v)"; return
  fi
  for c in $PY_CANDIDATES; do
    path="$(command -v "$c" 2>/dev/null)" || continue
    if v="$(py_version "$path")"; then
      if venv_works "$path"; then
        PY="$path"; say "Python: $PY ($v)"; return
      fi
      PY_NOVENV="$PY_NOVENV $path:$v"
    fi
    PY_SEEN="$PY_SEEN $c=${v:-?}"
  done
  if command -v uv >/dev/null 2>&1; then
    USE_UV=1
    say "Python: no usable system Python 3.$PY_MIN to 3.$PY_MAX; building the venv with uv (Python $UV_PYTHON)"
    return
  fi
  if [ -n "$PY_NOVENV" ]; then
    local first="${PY_NOVENV# }"; first="${first%% *}"
    die "Python ${first##*:} (${first%%:*}) cannot create a virtual environment: the venv/ensurepip module is missing. Install it:
$(venv_hint "${first%%:*}" "${first##*:}")
  then rerun install.sh."
  fi
  die "No supported Python found (BACTolomew needs Python 3.$PY_MIN to 3.$PY_MAX; found:${PY_SEEN:- none}). Install one:
$(python_hint)"
}

# --------------------------------------------------------------------------- #
# Container runtime (a warning only: a viewer-only install works without one)
# --------------------------------------------------------------------------- #
check_runtime() {
  local usable="" notes="" out rc
  if [ "$VIEWER" = 1 ]; then
    say "Viewer install: no container runtime needed (running predictions later needs one)"
    return 0
  fi
  local to=""; command -v timeout >/dev/null 2>&1 && to="timeout 30"
  if command -v docker >/dev/null 2>&1; then
    rc=0; out="$($to docker info 2>&1 >/dev/null)" || rc=$?
    if [ "$rc" = 0 ]; then
      if docker --version 2>/dev/null | grep -qi podman; then usable="docker (Podman)"; else usable="docker"; fi
    elif printf '%s' "$out" | grep -qi "permission denied"; then
      notes="${notes}
  Docker is installed but this user may not use it (permission denied on the Docker socket).
    Add yourself to the docker group, then log out and back in:  sudo usermod -aG docker \$USER"
    elif printf '%s' "$out" | grep -qi -E "cannot connect|is the docker daemon running|no such file"; then
      notes="${notes}
  Docker is installed but its daemon is not running:  sudo systemctl start docker"
    else
      notes="${notes}
  \`docker info\` failed: $(printf '%s' "$out" | tail -n1)"
    fi
  fi
  if command -v podman >/dev/null 2>&1; then
    if $to podman info >/dev/null 2>&1; then usable="${usable:+$usable, }podman"
    else notes="${notes}
  Podman is installed but \`podman info\` failed."; fi
  fi
  if command -v apptainer >/dev/null 2>&1; then usable="${usable:+$usable, }apptainer"
  elif command -v singularity >/dev/null 2>&1; then usable="${usable:+$usable, }singularity"; fi

  if [ -n "$usable" ]; then
    say "Container runtime: $usable"
    [ -z "$notes" ] || warn "Found but not usable:${notes}"
  else
    warn "No usable container runtime: the window can open existing results, but running the pipeline needs one.${notes}
  $([ -z "$notes" ] && echo "Install one of:" || echo "Fix that, or install one of:")
    Docker     https://docs.docker.com/engine/install/
    Podman     https://podman.io/docs/installation
    Apptainer  https://apptainer.org/docs/admin/main/installation.html
  then run \`./bactolomew setup\`."
  fi
}

# --------------------------------------------------------------------------- #
# Obtaining the code
# --------------------------------------------------------------------------- #
have_checkout() { [ -f "$1/src/annotation_explorer_polars/app.py" ]; }
dir_empty()     { [ ! -e "$1" ] || [ -z "$(ls -A "$1" 2>/dev/null)" ]; }

# Copy a repo directory or tarball into <dest>. For a git checkout, copy the working tree's
# tracked and untracked, not ignored files (so uncommitted work installs too, and gitignored
# data such as envs/ DB/ runs/, tens of GB, is skipped).
copy_tree() {  # copy_tree <src> <dest>
  local src="$1" dest="$2"
  mkdir -p "$dest" || die "cannot create $dest"
  case "$src" in
    *.tar.gz|*.tgz|*.tar)
      say "Extracting $src -> $dest"
      tar -xf "$src" -C "$dest" --strip-components=1 || die "failed to extract $src" ;;
    *)
      [ -d "$src" ] || die "--source '$src' is not a directory or a .tar.gz/.tgz/.tar tarball."
      if [ -e "$src/.git" ] && command -v git >/dev/null 2>&1; then
        say "Copying $src -> $dest (tracked and untracked files; gitignored data skipped)"
        local list; list="$(mktemp "${TMPDIR:-/tmp}/bactolomew-files.XXXXXX")"
        # safe.directory: a checkout owned by another user is still only read.
        git -c safe.directory='*' -C "$src" ls-files -z --cached --others --exclude-standard > "$list" \
          || { rm -f "$list"; die "git ls-files failed in $src"; }
        # Tracked files deleted in the working tree are listed too: copy only what exists.
        local f n=0
        while IFS= read -r -d '' f; do
          if [ -e "$src/$f" ] || [ -L "$src/$f" ]; then
            (cd "$src" && cp -P -p --parents "$f" "$dest/") || { rm -f "$list"; die "failed to copy $src/$f"; }
            n=$((n + 1))
          fi
        done < "$list"
        rm -f "$list"
        say "Copied $n files"
      else
        [ -e "$src/.git" ] && warn "git not found: copying all of $src, ignored data included."
        say "Copying $src -> $dest"
        cp -a "$src/." "$dest/" || die "failed to copy from $src"
      fi ;;
  esac
}

obtain_code() {
  # Precedence: explicit --source  >  running inside a checkout  >  clone from GitHub.
  local src kind
  if [ -n "$SOURCE" ]; then
    src="$SOURCE"; kind="source"
    if [ -d "$src" ]; then src="$(cd "$src" && pwd)"; fi
  elif have_checkout "$SCRIPT_DIR"; then
    src="$SCRIPT_DIR"; kind="inplace"
  else
    src=""; kind="clone"
  fi

  # Install dir (default: the checkout itself for in-place; else ~/BACTolomew).
  if [ -z "$INSTALL_DIR" ]; then
    if [ "$kind" = "inplace" ]; then INSTALL_DIR="$src"
    else INSTALL_DIR="${HOME:-$PWD}/BACTolomew"; fi
  fi
  mkdir -p "$INSTALL_DIR" || die "cannot create install dir: $INSTALL_DIR"
  INSTALL_DIR="$(cd "$INSTALL_DIR" && pwd)" || die "cannot access install dir: $INSTALL_DIR"

  if [ -n "$GITHUB_REF" ] && [ "$kind" != "clone" ]; then
    say "note: --ref '$GITHUB_REF' only applies when cloning; ignoring for $kind install."
  fi

  # Never write into a non-empty directory that is not already a BACTolomew checkout.
  if ! dir_empty "$INSTALL_DIR" && ! have_checkout "$INSTALL_DIR"; then
    die "Target $INSTALL_DIR is not empty and is not a BACTolomew checkout.
  Choose an empty --dir <path>, or clear it first."
  fi

  case "$kind" in
    inplace)
      if [ "$src" = "$INSTALL_DIR" ]; then say "Installing in-place at $INSTALL_DIR"
      else copy_tree "$src" "$INSTALL_DIR"; fi ;;
    source)
      copy_tree "$src" "$INSTALL_DIR" ;;
    clone)
      command -v git >/dev/null 2>&1 || die "git not found (needed to clone $GITHUB_URL). $(pkg_install_cmd git)"
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
}

# --------------------------------------------------------------------------- #
# venv + dependencies
# --------------------------------------------------------------------------- #
install_deps() {
  VENV="$INSTALL_DIR/.venv"
  if [ -x "$VENV/bin/python" ] && py_version "$VENV/bin/python" >/dev/null && "$VENV/bin/python" -m pip --version >/dev/null 2>&1; then
    say "Reusing the venv at $VENV ($("$VENV/bin/python" -V 2>&1))"
  elif [ "$USE_UV" = 1 ]; then
    say "Creating venv at $VENV with uv (Python $UV_PYTHON)"
    # UV_PYTHON_INSTALL_BIN=0: a Python uv downloads stays private (no python3.X links in ~/.local/bin).
    UV_PYTHON_INSTALL_BIN=0 uv venv --seed --python "$UV_PYTHON" "$VENV" || die "uv venv --python $UV_PYTHON failed."
  else
    say "Creating venv at $VENV"
    "$PY" -m venv "$VENV" || die "could not create the venv with $PY.
$(venv_hint "$PY" "$(py_version "$PY" || true)")"
  fi

  local reqs
  if [ "$HEADLESS" = 1 ]; then
    if [ "$USE_LOCK" = 1 ]; then reqs="requirements-pipeline.lock"; else reqs="requirements-pipeline.txt"; fi
  else
    if [ "$USE_LOCK" = 1 ]; then reqs="requirements.lock"; else reqs="requirements.txt"; fi
  fi
  [ -f "$INSTALL_DIR/$reqs" ] || die "$INSTALL_DIR/$reqs is missing."

  "$VENV/bin/python" -m pip install --upgrade pip >/dev/null || die "could not upgrade pip in $VENV (is the network reachable?)"
  say "Installing host Python dependencies from $reqs (this can take a minute)..."
  local log rc=0
  log="$(mktemp "${TMPDIR:-/tmp}/bactolomew-pip.XXXXXX")"
  if [ "$USE_LOCK" = 1 ]; then
    "$VENV/bin/python" -m pip install --require-hashes --only-binary :all: -r "$INSTALL_DIR/$reqs" 2>&1 | tee "$log" || rc=$?
  else
    "$VENV/bin/python" -m pip install --only-binary :all: -r "$INSTALL_DIR/$reqs" 2>&1 | tee "$log" || rc=$?
  fi
  if [ "$rc" != 0 ]; then
    if grep -q -E "No matching distribution|Could not find a version that satisfies" "$log"; then
      rm -f "$log"
      die "pip found no binary package for this system ($("$VENV/bin/python" -V 2>&1), glibc $(glibc_version)).
  BACTolomew needs Python 3.$PY_MIN to 3.$PY_MAX on x86_64 Linux with glibc $([ "$HEADLESS" = 1 ] && echo "$GLIBC_HEADLESS_MIN" || echo "$GLIBC_GUI_MIN") or newer.
  If the network is fine, check the Python version (--python) $([ "$HEADLESS" = 1 ] || echo "or install with --headless (no window)")."
    fi
    rm -f "$log"
    die "pip could not install the dependencies (see above)."
  fi
  rm -f "$log"
}

# --------------------------------------------------------------------------- #
# Qt system libraries
# --------------------------------------------------------------------------- #
# soname -> package, per package manager. apt: "pkg" or "pkg|t64-pkg".
qt_lib_package() {  # qt_lib_package <soname>
  local so="$1" a="" d="" z="" p=""
  case "$so" in
    libglib-2.0.so.*|libgthread-2.0.so.*|libgobject-2.0.so.*|libgio-2.0.so.*)
      a="libglib2.0-0|libglib2.0-0t64"; d="glib2"; p="glib2"
      case "$so" in libgthread*) z="libgthread-2_0-0";; libgobject*) z="libgobject-2_0-0";; libgio*) z="libgio-2_0-0";; *) z="libglib-2_0-0";; esac ;;
    libnss3.so|libnssutil3.so|libsmime3.so) a="libnss3"; d="nss"; z="mozilla-nss"; p="nss" ;;
    libnspr4.so|libplc4.so|libplds4.so) a="libnspr4"; d="nspr"; z="mozilla-nspr"; p="nspr" ;;
    libxkbcommon.so.*) a="libxkbcommon0"; d="libxkbcommon"; z="libxkbcommon0"; p="libxkbcommon" ;;
    libxkbcommon-x11.so.*) a="libxkbcommon-x11-0"; d="libxkbcommon-x11"; z="libxkbcommon-x11-0"; p="libxkbcommon-x11" ;;
    libXcomposite.so.*) a="libxcomposite1"; d="libXcomposite"; z="libXcomposite1"; p="libxcomposite" ;;
    libXdamage.so.*) a="libxdamage1"; d="libXdamage"; z="libXdamage1"; p="libxdamage" ;;
    libXrandr.so.*) a="libxrandr2"; d="libXrandr"; z="libXrandr2"; p="libxrandr" ;;
    libXfixes.so.*) a="libxfixes3"; d="libXfixes"; z="libXfixes3"; p="libxfixes" ;;
    libXtst.so.*) a="libxtst6"; d="libXtst"; z="libXtst6"; p="libxtst" ;;
    libXext.so.*) a="libxext6"; d="libXext"; z="libXext6"; p="libxext" ;;
    libXrender.so.*) a="libxrender1"; d="libXrender"; z="libXrender1"; p="libxrender" ;;
    libX11.so.*) a="libx11-6"; d="libX11"; z="libX11-6"; p="libx11" ;;
    libX11-xcb.so.*) a="libx11-xcb1"; d="libX11-xcb"; z="libX11-xcb1"; p="libx11" ;;
    libxkbfile.so.*) a="libxkbfile1"; d="libxkbfile"; z="libxkbfile1"; p="libxkbfile" ;;
    libxshmfence.so.*) a="libxshmfence1"; d="libxshmfence"; z="libxshmfence1"; p="libxshmfence" ;;
    libSM.so.*) a="libsm6"; d="libSM"; z="libSM6"; p="libsm" ;;
    libICE.so.*) a="libice6"; d="libICE"; z="libICE6"; p="libice" ;;
    libgbm.so.*) a="libgbm1"; d="mesa-libgbm"; z="libgbm1"; p="mesa" ;;
    libdrm.so.*) a="libdrm2"; d="libdrm"; z="libdrm2"; p="libdrm" ;;
    libasound.so.*) a="libasound2|libasound2t64"; d="alsa-lib"; z="libasound2"; p="alsa-lib" ;;
    libfontconfig.so.*) a="libfontconfig1"; d="fontconfig"; z="libfontconfig1"; p="fontconfig" ;;
    libfreetype.so.*) a="libfreetype6"; d="freetype"; z="libfreetype6"; p="freetype2" ;;
    libexpat.so.*) a="libexpat1"; d="expat"; z="libexpat1"; p="expat" ;;
    libGL.so.*|libGLX.so.*) a="libgl1"; d="mesa-libGL"; z="libglvnd"; p="libglvnd" ;;
    libEGL.so.*) a="libegl1"; d="mesa-libEGL"; z="libglvnd"; p="libglvnd" ;;
    libOpenGL.so.*) a="libopengl0"; d="libglvnd-opengl"; z="libglvnd"; p="libglvnd" ;;
    libdbus-1.so.*) a="libdbus-1-3"; d="dbus-libs"; z="libdbus-1-3"; p="dbus" ;;
    libgssapi_krb5.so.*|libkrb5.so.*|libk5crypto.so.*) a="libgssapi-krb5-2"; d="krb5-libs"; z="krb5"; p="krb5" ;;
    libcups.so.*) a="libcups2|libcups2t64"; d="cups-libs"; z="libcups2"; p="libcups" ;;
    libudev.so.*) a="libudev1"; d="systemd-libs"; z="libudev1"; p="systemd-libs" ;;
    libz.so.*) a="zlib1g"; d="zlib"; z="libz1"; p="zlib" ;;
    libxcb.so.*) a="libxcb1"; d="libxcb"; z="libxcb1"; p="libxcb" ;;
    libxcb-cursor.so.*) a="libxcb-cursor0"; d="xcb-util-cursor"; z="libxcb-cursor0"; p="xcb-util-cursor" ;;
    libxcb-icccm.so.*) a="libxcb-icccm4"; d="xcb-util-wm"; z="libxcb-icccm4"; p="xcb-util-wm" ;;
    libxcb-image.so.*) a="libxcb-image0"; d="xcb-util-image"; z="libxcb-image0"; p="xcb-util-image" ;;
    libxcb-keysyms.so.*) a="libxcb-keysyms1"; d="xcb-util-keysyms"; z="libxcb-keysyms1"; p="xcb-util-keysyms" ;;
    libxcb-render-util.so.*) a="libxcb-render-util0"; d="xcb-util-renderutil"; z="libxcb-render-util0"; p="xcb-util-renderutil" ;;
    libxcb-util.so.*) a="libxcb-util1"; d="xcb-util"; z="libxcb-util1"; p="xcb-util" ;;
    libxcb-*.so.*)  # the other xcb extensions: one package each on Debian and SUSE, all in libxcb elsewhere
      local ext="${so#libxcb-}" n="${so##*.so.}"; ext="${ext%%.so*}"; n="${n%%.*}"
      case "$ext" in *[0-9]) a="libxcb-$ext-$n" ;; *) a="libxcb-$ext$n" ;; esac
      z="$a"; d="libxcb"; p="libxcb" ;;
    *) return 1 ;;
  esac
  case "$PKG" in
    apt) if apt_t64; then echo "${a##*|}"; else echo "${a%%|*}"; fi ;;
    dnf) echo "$d" ;;
    zypper) echo "$z" ;;
    pacman) echo "$p" ;;
    *) return 1 ;;
  esac
}

# The full set the window needs on a minimal system, per package manager.
qt_full_list() {
  case "$PKG" in
    apt)
      local asound="libasound2" glib="libglib2.0-0"
      if apt_t64; then asound="libasound2t64"; glib="libglib2.0-0t64"; fi
      echo "libnss3 libnspr4 libxkbcommon0 libxkbcommon-x11-0 libxcomposite1 libxdamage1 libxrandr2" \
           "libxfixes3 libxtst6 libxkbfile1 libgbm1 $asound libxshmfence1 fontconfig libgl1 libegl1 $glib" \
           "libdbus-1-3 libgssapi-krb5-2 libxcb-cursor0 libxcb-icccm4 libxcb-image0 libxcb-keysyms1" \
           "libxcb-render-util0 libxcb-render0 libxcb-shape0 libxcb-util1 libxcb-xkb1" ;;
    dnf)
      echo "nss nspr libxkbcommon libxkbcommon-x11 libXcomposite libXdamage libXrandr libXfixes libXtst" \
           "libxkbfile mesa-libgbm alsa-lib libxshmfence fontconfig mesa-libGL mesa-libEGL glib2 dbus-libs" \
           "krb5-libs xcb-util-cursor xcb-util-wm xcb-util-image xcb-util-keysyms xcb-util-renderutil xcb-util" ;;
    zypper)
      echo "mozilla-nss mozilla-nspr libxkbcommon0 libxkbcommon-x11-0 libXcomposite1 libXdamage1 libXrandr2" \
           "libXfixes3 libXtst6 libxkbfile1 libgbm1 libasound2 libxshmfence1 fontconfig libglvnd" \
           "libglib-2_0-0 libgthread-2_0-0 libdbus-1-3 krb5 libxcb-cursor0 libxcb-icccm4 libxcb-image0" \
           "libxcb-keysyms1 libxcb-render-util0 libxcb-render0 libxcb-shape0 libxcb-util1 libxcb-xkb1" ;;
    pacman)
      echo "nss nspr libxkbcommon libxkbcommon-x11 libxcomposite libxdamage libxrandr libxfixes libxtst" \
           "libxkbfile mesa alsa-lib libxshmfence fontconfig libglvnd glib2 dbus krb5 xcb-util-cursor" \
           "xcb-util-wm xcb-util-image xcb-util-keysyms xcb-util-renderutil xcb-util libxcb" ;;
  esac
}

# After pip: import the Qt modules offscreen and ldd the Qt libraries and the xcb platform
# plugin, then name the missing system libraries and their packages.
check_qt_libs() {
  local import_err="" missing="" qt so pkgs="" unmapped="" f
  import_err="$(QT_QPA_PLATFORM=offscreen "$VENV/bin/python" -c 'import PyQt6.QtWidgets, PyQt6.QtWebEngineWidgets' 2>&1)" || true
  qt="$("$VENV/bin/python" -c 'import os, importlib.util as u; s = u.find_spec("PyQt6"); print(os.path.join(os.path.dirname(s.origin), "Qt6"))' 2>/dev/null)" || qt=""
  if [ -n "$qt" ] && [ -d "$qt" ] && command -v ldd >/dev/null 2>&1; then
    # The PyQt6 modules the app loads (ldd follows their Qt libraries), the xcb plugin and
    # the WebEngine helper process; Qt's optional modules (multimedia, wayland, ...) are left out.
    missing="$(for f in "$qt"/plugins/platforms/libqxcb.so "$qt"/libexec/QtWebEngineProcess \
                        "$qt"/../Qt{Core,Gui,Widgets,Network,Svg,SvgWidgets,OpenGL,OpenGLWidgets,WebChannel,WebEngineCore,WebEngineWidgets}.*so; do
                 [ -e "$f" ] && ldd "$f" 2>/dev/null | awk '/not found/ {print $1}'
               done | grep -v '^libQt6' | sort -u)" || true
  fi
  # The import error names the first missing library even without ldd.
  so="$(printf '%s\n' "$import_err" | grep -oE 'lib[^ :]+\.so[.0-9]*: cannot open shared object' | head -n1 | cut -d: -f1)" || true
  [ -n "$so" ] && missing="$(printf '%s\n%s\n' "$missing" "$so" | grep -v '^$' | sort -u)"

  if [ -z "$missing" ] && [ -z "$import_err" ]; then
    say "Qt libraries: PyQt6 and WebEngine load; the xcb platform plugin has all its libraries"
    return
  fi
  if [ -z "$missing" ]; then
    warn "PyQt6 / WebEngine do not import:
$(printf '%s\n' "$import_err" | tail -n 3 | sed 's/^/    /')"
    return
  fi
  for so in $missing; do
    if f="$(qt_lib_package "$so")"; then
      case " $pkgs " in *" $f "*) ;; *) pkgs="${pkgs:+$pkgs }$f" ;; esac
    else
      unmapped="${unmapped:+$unmapped }$so"
    fi
  done
  warn "The window needs system libraries this machine lacks (pipeline commands are not affected):
$(printf '%s ' $missing | fold -s -w 96 | sed 's/^/    /')"
  if [ -n "$pkgs" ]; then
    printf '  Install them with:\n    %s\n' "$(pkg_install_cmd $pkgs)" >&2
  fi
  [ -z "$unmapped" ] || printf '  No package known for: %s (search your package manager for these files).\n' "$unmapped" >&2
  if [ -n "$PKG" ]; then
    printf '  On a minimal or server install, the full set the window uses:\n    %s\n' "$(pkg_install_cmd $(qt_full_list))" >&2
  fi
  printf '  Then check with:  %s check\n' "$INSTALL_DIR/bactolomew" >&2
}

# --------------------------------------------------------------------------- #
# Launcher
# --------------------------------------------------------------------------- #
write_launcher() {
  LAUNCHER="$INSTALL_DIR/bactolomew"
  if [ "$HEADLESS" = 1 ]; then
    cat > "$LAUNCHER" <<EOF
#!/usr/bin/env bash
# BACTolomew (headless install: no Qt, so no window).
case "\${1:-gui}" in
  setup|run|check|db|--version|-h|--help) ;;
  *) echo "This is a headless BACTolomew install (no Qt): it supports setup, run, check and db." >&2
     echo "For the window, rerun install.sh without --headless on a machine with glibc $GLIBC_GUI_MIN or newer." >&2
     exit 1 ;;
esac
exec "$VENV/bin/python" "$INSTALL_DIR/src/bactolomew_main.py" "\$@"
EOF
  else
    cat > "$LAUNCHER" <<EOF
#!/usr/bin/env bash
exec "$VENV/bin/python" "$INSTALL_DIR/src/bactolomew_main.py" "\$@"
EOF
  fi
  chmod +x "$LAUNCHER"
}

# A viewer install says so in its config, unless a config exists already (that install
# keeps its mode; `./bactolomew setup --viewer` switches it).
write_viewer_config() {
  local config="$INSTALL_DIR/bactolomew_config.json"
  if [ -e "$config" ]; then
    say "Keeping the existing $config"
    return 0
  fi
  printf '{\n  "mode": "viewer"\n}\n' > "$config"
  mkdir -p "$INSTALL_DIR/runs" "$INSTALL_DIR/cache"
}

# --------------------------------------------------------------------------- #
main() {
  SELF="${BASH_SOURCE[0]:-$0}"   # this script's path ("bash" when piped)
  SCRIPT_DIR="$(cd "$(dirname "$SELF")" 2>/dev/null && pwd || echo "$PWD")"
  SOURCE=""; INSTALL_DIR=""; PYTHON_ARG=""; HEADLESS=0; VIEWER=0; USE_LOCK=1

  need_val() { [ "$#" -ge 2 ] || die "$1 requires a value"; }
  while [ $# -gt 0 ]; do
    case "$1" in
      --source)       need_val "$@"; SOURCE="$2"; shift 2 ;;
      --dir)          need_val "$@"; INSTALL_DIR="$2"; shift 2 ;;
      --ref|--branch) need_val "$@"; GITHUB_REF="$2"; shift 2 ;;
      --python)       need_val "$@"; PYTHON_ARG="$2"; shift 2 ;;
      --headless)     HEADLESS=1; shift ;;
      --viewer)       VIEWER=1; shift ;;
      --no-lock)      USE_LOCK=0; shift ;;
      -h|--help)      usage; exit 0 ;;
      *) die "unknown argument: $1 (see --help)" ;;
    esac
  done

  [ "$HEADLESS" = 1 ] && [ "$VIEWER" = 1 ] && die "--viewer needs the window, --headless has none: pick one."

  detect_distro
  # --- 1. prerequisites -----------------------------------------------------
  check_platform
  find_python
  check_runtime
  # --- 2. the code ----------------------------------------------------------
  obtain_code
  # --- 3. venv + deps -------------------------------------------------------
  install_deps
  # --- 4. launcher, then the Qt libraries ----------------------------------
  write_launcher
  [ "$VIEWER" = 1 ] && write_viewer_config
  [ "$HEADLESS" = 1 ] || check_qt_libs

  say "Done. Launch with:  $LAUNCHER"
  if [ "$VIEWER" = 1 ]; then
    say "Viewer install: open shared runs with Import Run... in the run list. To run predictions later: \`$LAUNCHER setup\`."
  elif [ "$HEADLESS" = 1 ]; then
    say "Headless install: \`$LAUNCHER check\`, \`$LAUNCHER setup --yes\` (tool image and databases), \`$LAUNCHER run --genomes DIR --out DIR\`."
  else
    say "The first start (or \`$LAUNCHER setup --yes\` without a window) gets the tools image and the databases."
  fi
}

main "$@"

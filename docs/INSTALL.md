# Installing BACTolomew

BACTolomew has two parts: the **window** (a Python program installed in its own folder) and
the **analysis tools**, which all run in one container image (`liviurotiul/bactolomew` on
Docker Hub) through Docker, Podman or Apptainer. The installer sets up the window and checks
the container runtime; the first start (or `./bactolomew setup`) then downloads the tools
image and the databases.

There are three kinds of install:

| Install | What it does | Needs |
|---|---|---|
| Full (default) | runs analyses and shows the results | Docker (or Podman, Apptainer), about 150 GB free for the first setup |
| Viewer only (`--viewer`) | opens runs someone else made and shared with you | about 1 GB, plus about 0.5 GB per shared run; no Docker |
| Headless (`--headless`) | runs analyses from the command line, no window (servers, HPC) | Docker, Podman or Apptainer, same disk as a full install |

A viewer install can become a full one later with `./bactolomew setup`.

**Contents**

- [Requirements](#requirements)
- [Linux](#linux): [Ubuntu and Debian](#ubuntu-and-debian), [Fedora](#fedora),
  [Rocky, Alma and RHEL](#rocky-alma-and-rhel), [openSUSE Leap and Tumbleweed](#opensuse-leap-and-tumbleweed),
  [Arch](#arch), [Other Linux](#other-linux), [HPC and servers](#hpc-and-servers-headless)
- [macOS](#macos): [Intel](#macos-intel), [Apple Silicon](#macos-apple-silicon-m1-to-m4)
- [Windows 10 and 11](#windows-10-and-11-through-wsl2)
- [Installer options](#installer-options)
- [First launch](#first-launch)
- [Common problems](#common-problems)
- [Uninstall](#uninstall)

## Requirements

| | Full or headless install | Viewer only |
|---|---|---|
| System | Linux x86_64 (glibc 2.34 or newer for the window, 2.28 for headless); macOS 13 or newer; Windows 10 (2004) or 11 through WSL2 | the same |
| Python | 3.10 to 3.14 (the installer says how to get one) | the same |
| Container runtime | Docker (tested), Podman or Apptainer/Singularity | none |
| Disk | about 150 GB free for the first setup (databases about 125 GB while unpacking, the tools image about 20 GB); about 100 GB stay afterwards, plus your runs | about 1 GB, plus about 0.5 GB per shared run of 150 genomes |
| Memory | 32 GB recommended; 16 GB works, more slowly | 8 GB is enough |
| CPU | x86_64 with AVX2 (Intel Haswell, AMD Zen or newer), or Apple Silicon | the same |

What was tested: the full and the viewer install in clean containers of Ubuntu 22.04 and 24.04,
Debian 12, Fedora 42, Rocky 9, openSUSE Leap 15.6 and Tumbleweed and Arch Linux, and the
headless install on Rocky 8. macOS and Windows are not yet tested on a real machine.

## Linux

On every distribution the steps are the same:

1. Install git and a Python with `venv` (commands below).
2. Install Docker and let your user run it (not needed for a viewer install).
3. Run the installer:

   ```bash
   git clone https://github.com/liviurotiul/BACTolomew.git
   cd BACTolomew
   bash install.sh              # or: bash install.sh --viewer
   ```

   Without git you can download the installer alone; it then clones BACTolomew into
   `~/BACTolomew`:

   ```bash
   curl -fsSL -o install.sh https://raw.githubusercontent.com/liviurotiul/BACTolomew/main/install.sh && bash install.sh
   ```

4. Start it with `./bactolomew` (see [First launch](#first-launch)).

The installer stops with the exact command to run when Python or a system library is missing:
run that command and start `install.sh` again. On a desktop system the window's libraries are
usually present already; the lists below are for minimal and server installs.

### Ubuntu and Debian

Ubuntu 22.04 or newer, Debian 12 or newer (older releases: headless only).

```bash
sudo apt update
sudo apt install git curl python3 python3-venv
```

Docker, from the distribution (the simplest) or from
[Docker's own repository](https://docs.docker.com/engine/install/ubuntu/):

```bash
sudo apt install docker.io
sudo usermod -aG docker $USER      # then log out and back in
```

The window's libraries on a minimal or server install (Ubuntu 24.04; on Ubuntu 22.04 and
Debian 12 write `libasound2 libglib2.0-0` instead of `libasound2t64 libglib2.0-0t64`):

```bash
sudo apt install libegl1 libgl1 libx11-xcb1 libx11-6 libxcomposite1 libxdamage1 libxext6 \
  libxfixes3 libxrandr2 libxtst6 libasound2t64 libdbus-1-3 libfontconfig1 libfreetype6 libgbm1 \
  libglib2.0-0t64 libgssapi-krb5-2 libnspr4 libnss3 libxcb-cursor0 libxcb-dri3-0 libxcb-icccm4 \
  libxcb-image0 libxcb-keysyms1 libxcb-randr0 libxcb-render-util0 libxcb-render0 libxcb-shape0 \
  libxcb-shm0 libxcb-sync1 libxcb-util1 libxcb-xfixes0 libxcb-xkb1 libxcb1 libxkbcommon-x11-0 \
  libxkbcommon0 libxkbfile1
```

### Fedora

```bash
sudo dnf install git python3
```

Docker, from Fedora (`moby-engine`) or from
[Docker's own repository](https://docs.docker.com/engine/install/fedora/):

```bash
sudo dnf install moby-engine
sudo systemctl enable --now docker
sudo usermod -aG docker $USER      # then log out and back in
```

Podman, which Fedora usually has already, works too.

The window's libraries on a minimal install:

```bash
sudo dnf install mesa-libEGL mesa-libGL libX11-xcb libX11 libXcomposite libXdamage libXext \
  libXfixes libXrandr libXtst alsa-lib dbus-libs fontconfig freetype mesa-libgbm nspr nss \
  xcb-util-cursor libxcb xcb-util-wm xcb-util-image xcb-util-keysyms xcb-util-renderutil xcb-util \
  libxkbcommon-x11 libxkbcommon libxkbfile
```

### Rocky, Alma and RHEL

Version 9 or newer for the window; version 8 works headless only (`bash install.sh --headless`).
The system Python of version 9 is 3.9, too old, so install 3.12 beside it (the installer finds
it by itself):

```bash
sudo dnf install git python3.12
```

Docker, from [Docker's repository](https://docs.docker.com/engine/install/centos/) (or use
Podman: `sudo dnf install podman`):

```bash
sudo dnf install dnf-plugins-core
sudo dnf config-manager --add-repo https://download.docker.com/linux/centos/docker-ce.repo
sudo dnf install docker-ce docker-ce-cli containerd.io
sudo systemctl enable --now docker
sudo usermod -aG docker $USER      # then log out and back in
```

The window's libraries on a minimal install are the same list as for [Fedora](#fedora).

### openSUSE Leap and Tumbleweed

Leap 15.6 or newer, or Tumbleweed. Leap's system Python is 3.6, so install 3.12 beside it:

```bash
sudo zypper install git python312
```

Docker:

```bash
sudo zypper install docker
sudo systemctl enable --now docker
sudo usermod -aG docker $USER      # then log out and back in
```

The window's libraries and a font on a minimal install (on Leap leave out `krb5`):

```bash
sudo zypper install dejavu-fonts libglvnd libX11-xcb1 libX11-6 libXcomposite1 libXdamage1 \
  libXext6 libXfixes3 libXrandr2 libXtst6 libasound2 libdbus-1-3 libfontconfig1 libfreetype6 \
  libgbm1 krb5 libgthread-2_0-0 mozilla-nspr mozilla-nss libxcb-cursor0 libxcb-dri3-0 \
  libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-randr0 libxcb-render-util0 libxcb-render0 \
  libxcb-shape0 libxcb-shm0 libxcb-sync1 libxcb-util1 libxcb-xfixes0 libxcb-xkb1 libxcb1 \
  libxkbcommon-x11-0 libxkbcommon0 libxkbfile1
```

### Arch

```bash
sudo pacman -S --needed git python
sudo pacman -S --needed docker
sudo systemctl enable --now docker
sudo usermod -aG docker $USER      # then log out and back in
```

The window's libraries and a font on a minimal install:

```bash
sudo pacman -S --needed ttf-dejavu libglvnd libx11 libxcomposite libxdamage libxext libxfixes \
  libxrandr libxtst alsa-lib fontconfig freetype2 mesa nspr nss xcb-util-cursor libxcb \
  xcb-util-wm xcb-util-image xcb-util-keysyms xcb-util-renderutil xcb-util libxkbcommon-x11 \
  libxkbcommon libxkbfile
```

### Other Linux

Any x86_64 distribution with glibc 2.34 or newer can run the window (check with
`ldd --version`); glibc 2.28 to 2.33 can run the pipeline with `--headless`. You need git,
Python 3.10 to 3.14 with `venv`, a container runtime and, for the window, the X11/xcb, GL, NSS
and font libraries listed above. If your distribution has no suitable Python, the installer
uses [uv](https://docs.astral.sh/uv/) when it is on your PATH (`curl -LsSf https://astral.sh/uv/install.sh | sh`),
or point it at any Python with `--python /path/to/python3.12`.

Not supported: musl systems such as Alpine, and the window on 64-bit ARM Linux (no PyQt6Ads
build for it).

### HPC and servers (headless)

On a cluster or a server without a display, install the pipeline only:

```bash
git clone https://github.com/liviurotiul/BACTolomew.git
cd BACTolomew
bash install.sh --headless
./bactolomew check
./bactolomew setup --yes --runtime apptainer      # or docker, podman; --sif PATH for where the image goes
./bactolomew run --genomes /path/to/assemblies --out /path/to/run
```

Most clusters have Apptainer or Singularity rather than Docker; `setup` then converts the tools
image to a `.sif` file (about 15 GB more on the same disk). The Apptainer path has not been
tested on a real cluster yet; please try one small run first. Inside a Slurm job the run uses
the job's cores and memory (`--cpus-per-task`, `--mem`). Copy the run folder to a machine with
the window, or share it with Export for Sharing, to look at the results.

## macOS

Not yet tested on a real Mac. Needs macOS 12 or newer (13 or newer for the window).

### macOS Intel

1. Install the command line tools (git) and Python 3.12, with [Homebrew](https://brew.sh):

   ```bash
   xcode-select --install
   brew install python@3.12
   ```

   or with the installer from [python.org](https://www.python.org/downloads/macos/) (3.12 or 3.13).

2. Install [Docker Desktop for Mac](https://docs.docker.com/desktop/setup/install/mac-install/)
   (Intel chip), start it, and in Settings > Resources set the memory to at least 16 GB and
   the disk limit to at least 150 GB; Apply & restart. Not needed for a viewer install.
   Docker Desktop is free for personal use, education and small companies.

3. In Terminal:

   ```bash
   git clone https://github.com/liviurotiul/BACTolomew.git
   cd BACTolomew
   bash install.sh              # or: bash install.sh --viewer
   ./bactolomew
   ```

There is no app bundle yet: start BACTolomew from Terminal with `./bactolomew` in its folder.

**Keep runs on a case-sensitive volume.** The Mac disk ignores the case of file names by
default, and some tool outputs differ only by case. Make a separate volume for your runs: Disk
Utility > select your disk > Add Volume (+), format "APFS (Case-sensitive)". It shares free
space with the rest of the disk. Sample names that differ only by case (`S2` and `s2`) are
refused for the same reason.

### macOS Apple Silicon (M1 to M4)

The same steps as for [Intel](#macos-intel), with Docker Desktop for Apple silicon. Run the
installer in a native Terminal, not one opened with Rosetta (the installer warns if it is).

**The tools image.** The published tools image is for Intel (x86_64) only; an Apple Silicon
image is not published yet. Two ways to get it:

- **Build it on your Mac (what setup does by itself).** The first `./bactolomew setup` (or the
  first-time setup window) sees that Docker runs arm64, skips the published Intel image (it
  does not download it) and builds the image from `docker/Dockerfile`. This takes about an hour and runs at full native
  speed. FastBAPS (the default lineage method) is included; fastGEAR is not, because it exists
  for x86_64 only. Docker Desktop needs at least 16 GB of memory for the build.
- **Run the Intel image under Rosetta (untested).** In Docker Desktop, Settings > General, turn
  on "Use Rosetta for x86_64/amd64 emulation on Apple Silicon", then pull the image before
  setup:

  ```bash
  docker pull --platform linux/amd64 "liviurotiul/bactolomew:0.3.0@sha256:7ed7408c88aac76f6c00b2417862e24f1dd3d91fe9ea89c5d62d57df7770ffde"
  ```

  Setup then finds it and does not build (it notes that the image runs emulated). If you later
  build a native image as well, runs use the native one. It is the same image as on Linux, but the tools run
  emulated, more slowly, and some may not run at all (CheckM2 most likely needs macOS 15 or
  newer under Rosetta).

We recommend the first way until an Apple Silicon image is published.

## Windows 10 and 11 (through WSL2)

Not yet tested on a real Windows PC.

BACTolomew is a Linux program; on Windows it runs inside WSL2 (the Windows Subsystem for
Linux). Its window opens on the Windows desktop through WSLg and its tools run in Docker. You
need Windows 10 version 2004 (build 19041) or newer, or Windows 11; for the window, Windows 10
21H2 or newer or Windows 11 with the current WSL (`wsl --update`).

In a PowerShell window (not as administrator; the script asks for administrator rights only
to turn WSL on):

```powershell
Invoke-WebRequest -Uri https://raw.githubusercontent.com/liviurotiul/BACTolomew/main/install_windows.ps1 -OutFile install_windows.ps1
powershell -ExecutionPolicy Bypass -File .\install_windows.ps1
```

For a viewer install add `-Viewer`. Other options: `-Distro <name>` (default Ubuntu-24.04,
installed when missing), `-Dir <linux path>` (default `~/BACTolomew`), `-Ref <tag or branch>`
(default main), `-Docker Desktop|Engine|Skip` (default: ask), `-DesktopShortcut`.

The script checks and, where needed, installs WSL2 and the Linux distribution, the Linux
packages the window needs, and Docker, then runs `install.sh` inside the distribution and adds
BACTolomew to the Start Menu. When it says Windows must restart (the first time WSL is turned
on), restart and run it again with the same options; it continues where it stopped. It is safe
to run again at any time, for example after installing Docker.

**Docker.** Two choices:

- Docker Desktop for Windows (recommended, the simplest): Settings > General, "Use the WSL 2
  based engine" on; Settings > Resources > WSL integration, your distribution on. Docker
  Desktop is free for personal use, education and small companies; larger companies need a
  paid subscription.
- Docker Engine inside the distribution (free for everyone, command line only): the script
  turns on systemd in the distribution and installs it. Do not use both for one distribution.

**Starting it:** the Start Menu entry, or in PowerShell (with your distribution and install
folder):

```powershell
wsl -d Ubuntu-24.04 --cd ~ -- ~/BACTolomew/bactolomew
```

**Keep your data inside Linux.** Genomes and runs on the Windows drive (C:\, seen in Linux as
`/mnt/c`) work but are many times slower, because every file access crosses from Linux to
Windows. Copy inputs into your Linux home first, for example
`cp -r /mnt/c/Users/<you>/genomes ~/genomes`. In File Explorer the Linux files are under
`\\wsl.localhost\Ubuntu-24.04\home\<user>` (also "Linux" in the navigation pane).

**Memory.** WSL2 gets half of the PC's memory by default, and Docker Desktop's containers share
that limit (with the WSL 2 engine its own Resources sliders do not apply). BACTolomew analyses
want about 32 GB. To give WSL more, create `%USERPROFILE%\.wslconfig` with

```ini
[wsl2]
memory=48GB
```

(leave a few GB for Windows), then run `wsl --shutdown` and start BACTolomew again.

**Disk.** First-time setup needs about 150 GB free. The distribution's files live in a virtual
disk (ext4.vhdx) that grows on the drive where the distribution is installed (C: by default)
and does not shrink by itself when files are deleted; Docker Desktop keeps its images in a
virtual disk of its own (Settings > Resources > Advanced > Disk image location). On a recent
WSL, `wsl --manage Ubuntu-24.04 --set-sparse true` lets the disk give space back, and
`wsl --manage Ubuntu-24.04 --move D:\WSL` moves it to another drive.

Windows on ARM is not tested: there is no ARM tools image yet, so setup would build one itself
(several hours, without fastGEAR); a `-Viewer` install is the safer choice there.

## Installer options

`bash install.sh [options]` (Linux, macOS and inside WSL):

| Option | What it does |
|---|---|
| `--viewer` | the window only, to open shared runs; no container runtime, tools image or databases |
| `--headless` | the pipeline only, no window (`setup`, `run`, `check`); the only choice on glibc 2.28 to 2.33 |
| `--dir <folder>` | where to install (default: the checkout it runs from, or `~/BACTolomew`) |
| `--ref <tag or branch>` | which version to clone, for example `v0.3.0` |
| `--python <path>` | which Python to use (also `PYTHON=...`) |
| `--source <folder or .tar.gz>` | install from a local copy instead of cloning |
| `--no-lock` | the newest package versions instead of the tested ones |

The installer creates a virtual environment (`.venv`) in the install folder from the tested,
hash-pinned package list (`requirements.lock`, or `requirements-pipeline.lock` for
`--headless`) and writes the `./bactolomew` launcher. It installs no bioinformatics tools and
changes nothing outside the install folder. It is safe to run again, for example after installing a missing package.

## First launch

```bash
./bactolomew
```

The first start asks what the install is for: **Full (run analyses)** or **Viewer only (open
results from others)**. A full install then runs first-time setup: the tools image (about 20 GB,
pulled from Docker Hub), the Bakta database (a 39 GB download that unpacks to 76 GB, 1 to 2
hours) and the CheckM2 database. Parts already installed are kept, so an interrupted setup can
simply be started again.

Without a window: `./bactolomew setup --yes` does the same, `./bactolomew check` shows what is
ready (Python, Qt, container runtime, image, databases, CPU and memory budget).

## Common problems

**`permission denied` on the Docker socket.** Your user is not in the `docker` group yet:
`sudo usermod -aG docker $USER`, then log out and back in (or reboot). `./bactolomew check`
should then show Docker as usable.

**`Cannot connect to the Docker daemon`.** Docker is not running: `sudo systemctl enable --now docker`
on Linux; on macOS and Windows start Docker Desktop and wait until it says the engine is running.

**The installer says a system library is missing.** Install the packages it prints (or the
full list for your distribution above), then run `./bactolomew check`.

**The window shows boxes instead of text.** No font is installed (minimal Arch and openSUSE
systems, some servers and WSL distributions): `sudo pacman -S ttf-dejavu`,
`sudo zypper install dejavu-fonts`, `sudo apt install fonts-dejavu-core` or
`sudo dnf install dejavu-sans-fonts`. The installer warns about this and prints the command.

**No suitable Python.** The installer prints the command for your system. Any Python 3.10 to
3.14 works: `bash install.sh --python /path/to/python3.12`.

**"glibc ... too old".** The window needs glibc 2.34 or newer (Ubuntu 22.04, Debian 12,
RHEL 9 or newer). On an older system use `--headless` and view the results elsewhere.

**Older CPU without AVX2.** `./bactolomew check` says so; in the install folder run
`.venv/bin/pip install 'polars[rtcompat]'`.

**Docker Desktop has too little memory or disk (macOS, Windows).** macOS: Settings >
Resources, memory 16 GB or more and the disk limit 150 GB or more. Windows: see
[Memory](#windows-10-and-11-through-wsl2) above (`.wslconfig`).

**Windows: the window does not appear.** WSLg is missing or out of date: in PowerShell run
`wsl --update`, then `wsl --shutdown`, and start BACTolomew again. Windows 10 needs 21H2 or
newer for WSLg.

**Setup is slow on the tools image.** The image is about 20 GB (about 7 GB to download). If the download fails, setup
builds the image locally instead, which takes about an hour; check the network and run setup
again to download it.

## Uninstall

Everything BACTolomew installs is in its install folder (the program, `.venv`, the databases
`bakta_db/`, `checkm2_db/` and `DB/`, your `runs/`, `cache/` and `bactolomew_config.json`),
plus the container image and a few per-user files. Deleting the install folder also deletes
your runs, so move any you want to keep first.

**Linux:**

```bash
rm -rf ~/BACTolomew                                        # or your --dir or checkout
docker image rm liviurotiul/bactolomew:0.3.0               # podman image rm with Podman
rm -rf ~/.config/BACTolomew ~/.cache/bactolomew ~/.cache/bactolomew_tiles
```

With Apptainer, the `.sif` file is in the install folder unless you chose another path with `--sif`.

**macOS:**

```bash
rm -rf ~/BACTolomew                                        # or your checkout
docker image rm liviurotiul/bactolomew:0.3.0               # or the locally built image (docker image ls)
rm -rf ~/.cache/bactolomew ~/.cache/bactolomew_tiles
rm -f ~/Library/Preferences/app.bactolomew.BACTolomew.plist
```

Docker Desktop is removed from its own menu (Troubleshoot > Uninstall).

**Windows:** delete the BACTolomew shortcut from the Start Menu (and the desktop), then either
remove only BACTolomew inside the distribution:

```powershell
wsl -d Ubuntu-24.04 -- rm -rf ~/BACTolomew ~/.config/BACTolomew ~/.cache/bactolomew ~/.cache/bactolomew_tiles
wsl -d Ubuntu-24.04 -- docker image rm liviurotiul/bactolomew:0.3.0
```

or, if the distribution was installed only for BACTolomew, remove all of it (this deletes every
file in it): `wsl --unregister Ubuntu-24.04`.

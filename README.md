# BACTolomew

BACTolomew is a bacterial comparative-genomics platform: a desktop **annotation explorer**
(PyQt6 / WebEngine) backed by a reproducible, containerized **analysis pipeline**. It takes
bacterial assemblies from QC through annotation, typing, pangenome and phylogenetic analysis,
and lets you explore the results interactively: tables, trees, pangenome views, plasmid
views and geographic maps.

![BACTolomew explore window, cgMLST minimum spanning tree, pangenome clustering, geographic map and gene neighborhoods](docs/images/bactolomew_overview.png)

A run of 153 *Acinetobacter baumannii* genomes: the explore window (phylogeny, gene
presence/absence, map, plasmidome), the cgMLST minimum spanning tree, pangenome clustering,
the geographic map and the *bla*<sub>OXA-23</sub> gene neighborhoods.

## Features

- **Pipeline** (one command or one button): assembly QC, Bakta annotation, MLST,
  AMRFinderPlus (acquired genes and point mutations), MOB-suite plasmids, gene categories
  from the run's own evidence, Panaroo pangenome, a recombination-aware core-genome
  phylogeny (core alignment, FastBAPS lineages, Gubbins per lineage, IQ-TREE), a Snippy
  SNP tree, chewBBACA cgMLST and searches of the databases you add yourself. Every step is
  resumable and reused only when its inputs match.
- **Explorer**: dockable, multi-instance panels in one or several windows, all following one
  filter set and one sample set: Sample tree, Annotations, Genome track, Map, Phylogeny,
  Gene presence/absence, Metadata summary, Plasmidome, Pangenome clustering, cgMLST MST,
  cgMLST MST (compact), Gene neighborhoods and Sequence alignment. Built-in and saved
  layouts (Window > Layouts).
- **Metadata**: a CSV with a `Sample` column, or "Use run data only" (species, MLST,
  fastGEAR clades and assembly QC); More > Change Metadata... on the run's card in the
  run list swaps it later.
- **Shared runs**: More > Export for Sharing... packs what the explorer needs into one file;
  Import Run... in the run list opens it on another install, without the pipeline.
- **Runs on a laptop, a workstation or an HPC node**: the same `bactolomew` command opens
  the window or runs headless, with Docker, Podman (rootful or rootless) or
  Apptainer/Singularity.

## Architecture

BACTolomew is split into a **thin host layer** and a **fat container layer**:

- **Host (this repo)**: a lightweight Python venv running the PyQt6/WebEngine GUI and the
  in-process pipeline orchestration. No bioinformatics tools are installed on the host.
- **Container (`liviurotiul/bactolomew` on Docker Hub)**: a single image carrying every bioinformatics tool.
  Each tool call is dispatched as `conda run -n <env>` inside the container, with host
  directories bind-mounted at identical paths so absolute paths resolve unchanged. See
  [`docker/Dockerfile`](docker/Dockerfile).

This keeps the GUI fast and portable while making the analysis environment fully reproducible.

## Pipeline tools

The conda environments are defined in [`yamls/`](yamls/), pinned exactly in
[`docker/locks/`](docker/locks/) and baked into the tools image. Versions from
[`docker/locks/tool_versions.json`](docker/locks/tool_versions.json):

| Stage              | Tools |
| ------------------ | ----- |
| Assembly QC        | QUAST 5.3.0, CheckM2 1.1.0 |
| Annotation         | Bakta 1.10.4 |
| AMR                | AMRFinderPlus 4.0.23 (`--plus`, point mutations with `--organism` from the MLST species) |
| Typing             | mlst 2.35.0 (PubMLST data of 2025-12-29), chewBBACA 3.5.4 (cgMLST) |
| Plasmids           | MOB-suite 3.1.9 |
| Pangenome          | Panaroo 1.8.0 |
| Core phylogeny     | FastBAPS 1.0.8 (the default lineages; fastGEAR from Docker Hub `pgcbioinfo/fastgear`, pinned by digest, MATLAB Runtime R2016a, optional and x86_64 only), Gubbins 3.4.3, IQ-TREE 2.4.0 |
| SNP phylogeny      | Snippy 4.6.0, RAxML-NG 1.2.2 |
| Database searches  | DIAMOND 2.2.3 |
| Sequence alignment | MAFFT 7.526 (explorer) |

## Installation

On every system the installer checks what is missing and prints the exact command to fix it.
Add `--viewer` (`-Viewer` on Windows) to install only the viewer, for looking at runs others
shared with you: no Docker, no databases, about 1 GB.

**Requirements** for a full install: about 150 GB free for the first setup (about 100 GB stay),
32 GB of RAM recommended, and Docker, Podman or Apptainer. Linux needs glibc 2.34 or newer for
the window (2.28 or newer headless); macOS needs version 13 or newer.

### Ubuntu and Debian

Ubuntu 22.04 or newer, Debian 12 or newer.

```bash
sudo apt install git python3 python3-venv docker.io
sudo usermod -aG docker $USER          # then log out and back in
git clone https://github.com/liviurotiul/BACTolomew.git && cd BACTolomew
bash install.sh
```

### Fedora

```bash
sudo dnf install git python3 moby-engine
sudo systemctl enable --now docker
sudo usermod -aG docker $USER          # then log out and back in
git clone https://github.com/liviurotiul/BACTolomew.git && cd BACTolomew
bash install.sh
```

Podman, which Fedora usually has, works instead of Docker.

### Rocky, Alma and RHEL 9

```bash
sudo dnf install git python3.12 podman
git clone https://github.com/liviurotiul/BACTolomew.git && cd BACTolomew
bash install.sh
```

Docker from [Docker's repository](https://docs.docker.com/engine/install/centos/) works too.
Version 8 runs headless only (`bash install.sh --headless`).

### openSUSE Leap and Tumbleweed

```bash
sudo zypper install git python312 docker dejavu-fonts
sudo systemctl enable --now docker
sudo usermod -aG docker $USER          # then log out and back in
git clone https://github.com/liviurotiul/BACTolomew.git && cd BACTolomew
bash install.sh
```

### Arch

```bash
sudo pacman -S --needed git python docker ttf-dejavu
sudo systemctl enable --now docker
sudo usermod -aG docker $USER          # then log out and back in
git clone https://github.com/liviurotiul/BACTolomew.git && cd BACTolomew
bash install.sh
```

### macOS (Intel)

1. Install Python 3.12 with [Homebrew](https://brew.sh) (`brew install python@3.12`) or from
   [python.org](https://www.python.org/downloads/macos/).
2. Install [Docker Desktop](https://docs.docker.com/desktop/setup/install/mac-install/), start
   it, and in Settings > Resources give it at least 16 GB of memory and 150 GB of disk.
3. In Terminal:

   ```bash
   git clone https://github.com/liviurotiul/BACTolomew.git && cd BACTolomew
   bash install.sh
   ```

Keep your runs on a case-sensitive APFS volume (Disk Utility > Add Volume > "APFS
(Case-sensitive)").

### macOS (Apple Silicon, M1 to M4)

The same steps as for Intel, with Docker Desktop for Apple silicon. The published tools image is
for Intel only, so the first `./bactolomew setup` builds it on your Mac (about an hour, at native
speed). fastGEAR is not available on Apple Silicon; FastBAPS, the default, is.

### Windows 10 and 11

BACTolomew runs inside WSL2, Windows' built-in Linux; its window opens on the Windows desktop.
In PowerShell (not as administrator):

```powershell
Invoke-WebRequest -Uri https://raw.githubusercontent.com/liviurotiul/BACTolomew/main/install_windows.ps1 -OutFile install_windows.ps1
powershell -ExecutionPolicy Bypass -File .\install_windows.ps1
```

The script sets up WSL2 and Ubuntu, Docker (Docker Desktop or Docker inside Linux, it asks),
installs BACTolomew and adds it to the Start Menu. If it asks you to restart Windows, restart and
run it again; it continues where it stopped.

### HPC and servers

```bash
git clone https://github.com/liviurotiul/BACTolomew.git && cd BACTolomew
bash install.sh --headless
./bactolomew setup --yes --runtime apptainer      # or docker, podman
./bactolomew run --genomes /path/to/assemblies --out /path/to/run
```

### After installing

Start BACTolomew with `./bactolomew` in its folder (on Windows, from the Start Menu). The first
start downloads the tools image and the databases.

The Linux installs above are tested; macOS and Windows are not yet tested on a real machine.
Minimal systems may need more packages (the installer names them), and the full guide with
common problems and how to uninstall is in [docs/INSTALL.md](docs/INSTALL.md).

## First launch and setup

```bash
./bactolomew
```

On first launch BACTolomew asks what the install is for: **Full (run analyses)** or **Viewer
only (open results from others)**. A viewer downloads nothing: open runs shared with you through
Import Run... in the run list; features that need the tools (such as aligning sequences from a
gene card, adding cgMLST schemas or databases, starting or resuming a run) are greyed out with
the reason. Settings > "Set Up Pipeline…" (or `./bactolomew setup`) makes it a full install
later; `./bactolomew setup --viewer` sets up a viewer without a window, and `./bactolomew check`
shows the mode.

For a full install, a **first-time setup** window gets everything ready by itself: the tools image,
the full Bakta database (with its AMRFinderPlus database, which the AMRFinderPlus step uses) and
the CheckM2 database, then a quick check that makes the `DB/` folder; parts already installed are
kept. The full Bakta database is a 39.4 GB download that unpacks to 76 GB and takes about 1 to 2
hours; before it starts anything, setup checks for about 125 GB free (download and unpacked
database side by side, 115.7 GB, plus CheckM2 4.8 GB and a 5 GB margin; about 15 GB more for a
Singularity `.sif` on the same disk). The tools image is pulled when `container.image` names a
registry; until a published image is configured it is built locally from `docker/Dockerfile`
(about an hour, needs network access to conda-forge, bioconda and Docker Hub). Setup writes a
machine-local `bactolomew_config.json` (git-ignored).

Then the runs window lists your runs (View, Monitor, Resume, Stop) and starts a new one from a
folder of assembled genomes ("Add Assemblies…", "Start New Run"). Its gear button opens
Settings: compute limits, cgMLST schemes and databases. Opening a run asks for its metadata
(a CSV with a `Sample` column, or "Use run data only").

## Running

### Annotation explorer (GUI)

`./bactolomew` (or `./bactolomew gui`) runs `src/bactolomew_main.py`, which starts the app. From
a checkout without the wrapper:

```bash
cd src
../.venv/bin/python -m annotation_explorer_polars.app
```

### Command line, headless and HPC

The same entry point works without a display (only the window needs Qt):

```bash
./bactolomew check                                   # Python, Qt, runtime, image, databases, CPU/RAM budget
./bactolomew setup --yes [--runtime docker|podman|singularity|apptainer] [--sif PATH] [--dry-run]
./bactolomew run --genomes DIR --out DIR [--threads N] [--resume]
./bactolomew db --list | --add FASTA --name NAME --category CAT | --remove NAME
./bactolomew --version
```

`setup` is first-time setup with text progress; a config without a runtime gets the first
usable one found. `run` reads the databases from `bactolomew_config.json` (no `--bakta_db`
or catalog paths needed; `run --help` lists the overrides) and refuses an output folder that
already holds a run unless `--resume` is given. The core and RAM budget follows the CPU
affinity, the cgroup limits (v2 `cpu.max` / `memory.max`, v1 fallbacks) and the Slurm
allocation (`SLURM_CPUS_PER_TASK`, `SLURM_MEM_PER_NODE` / `SLURM_MEM_PER_CPU`), capped by
`resources.max_threads` / `max_ram_gb`, and is logged at the start of each run. On a cluster
without Docker, use Apptainer/Singularity (`setup --runtime apptainer [--sif PATH]`).

## Databases

BACTolomew ships no marker or compound databases. The gene categories come from what the
pipeline already computes (step `bakta_evidence`):

- **AMR**, point mutations, stress/metal/biocide (BMRG) and AMRFinderPlus "plus" virulence:
  AMRFinderPlus, with the AMRFinderPlus database that comes with the Bakta database. A call on
  no Bakta gene is kept as a gene of its own on its contig coordinates.
- **Virulence** and **mobile elements**: Bakta's VFDB and IS hits (and, for mobile elements and
  phages, Bakta product names, marked as product name matches).
- **Plasmid**: MOB-suite's replicon and relaxase (and MPF) hits.

**Your own databases**: Settings > Databases ("Add database…") or `./bactolomew db --add`
take a protein FASTA (plain or `.gz`) with a category (AMR, BMRG, VIRULENCE, MOBILE, PHAGE,
PLASMID, or COMPOUND for a reference locus). They are searched with DIAMOND after Bakta from the
next run on, and their hits are shown as "<name> (U)". You are responsible for the license of
any database you add (the optional note can record it).

**cgMLST schemes** are installed from Settings ("Add schema…" for files you downloaded, one or
many; "Download from cgMLST.org…") or with `PYTHONPATH=src .venv/bin/python -m pipeline.CgmlstInstall`
(`--list`, `--status`, `--install <species>`, `--install-all`, `--import <paths>`; see `--help`).
First-time setup makes `DB/cgmlst` with the scheme species catalog
(`src/pipeline/data/scheme_species_catalog.json`). cgMLST.org schemes are for non-commercial use
only and may not be redistributed; the terms are shown once before the first download. A run
whose species has a cgMLST.org scheme that is not installed installs it during its chewBBACA
step (the GUI asks for the terms before such a run starts; without them a command-line run stops
at that step until they are accepted with `python -m pipeline.CgmlstInstall --accept-terms`). A
species without any scheme stops the run at the chewBBACA step until one is imported ("Add
schema…" or `--import`); Resume then continues. An installed scheme grows with every run, so it
is only rebuilt with `--replace`, and the previous version is kept under `.replaced/`.

The Bakta, AMRFinderPlus, CheckM2, MOB-suite and PubMLST data keep their own licenses and
citations; cite them with the tools above.

## Outputs

A run folder (by default under `runs/`) holds one folder per step (`assembly_qc/`, `bakta/`,
`mlst/`, `amrfinder/`, `mobsuite/`, `evidence/`, `panaroo/`, `core_alignment/`, `fastbaps/` or `fastgear/`,
`gubbins/`, `snippy/`, `chewbbaca/`, `dmnd/`), the trees with their index (`trees/`,
`trees/trees.json`) and three records: `run_config.json` (settings and the environment: image,
tool versions, databases), `run_state.json` (step progress, warnings) and `run_report.json`
(the provenance of the finished run). The explorer keeps its parquet cache under `cache/`.

## Development

Run straight from a local checkout (no `install.sh`):

```bash
git clone https://github.com/liviurotiul/BACTolomew.git
cd BACTolomew
python3 -m venv .venv                        # Python 3.10 to 3.14
.venv/bin/pip install --require-hashes --only-binary :all: -r requirements.lock
cd src && ../.venv/bin/python -m annotation_explorer_polars.app
```

The module must be launched from `src/` (its parent is the import root). Tests are unittest
style and run offscreen:

```bash
QT_QPA_PLATFORM=offscreen PYTHONPATH=src .venv/bin/python -m unittest discover -s tests -p 'test_*.py'
```

The tools image is built with `bash docker/build.sh [version]`; publishing it is described in
[`docker/PUBLISHING.md`](docker/PUBLISHING.md).

## Configuration

First-time setup writes `bactolomew_config.json` with the paths to the databases, the
container runtime and image, and resource limits. Most of it can be changed later in
**Settings**.

## Citation

If you use BACTolomew, please cite the release you used, e.g.
"BACTolomew v0.2.0, https://github.com/liviurotiul/BACTolomew/releases/tag/v0.2.0",
together with the tools listed above.

## License

BACTolomew is free software, licensed under the GNU General Public License version 3. See
[`LICENSE`](LICENSE).

# BACTolomew

BACTolomew is a bacterial comparative-genomics platform: a desktop **annotation explorer**
(PyQt6 / WebEngine) backed by a reproducible, containerized **analysis pipeline**. It takes
bacterial assemblies from QC through annotation, typing, pangenome and phylogenetic analysis,
and lets you explore the results interactively — tables, trees, pangenome graphs, plasmid
views and geographic maps.

> Formerly developed under the name *Contexar*.

## Architecture

BACTolomew is split into a **thin host layer** and a **fat container layer**:

- **Host (this repo)** — a lightweight Python venv running the PyQt6/WebEngine GUI and the
  in-process pipeline orchestration. No bioinformatics tools are installed on the host.
- **Container (`bactolomew-tools`)** — a single Docker/Singularity image carrying every
  bioinformatics tool. Each tool call is dispatched as `conda run -n bactolomew-<tool>`
  inside the container, with host directories bind-mounted at identical paths so absolute
  paths resolve unchanged. See [`docker/Dockerfile`](docker/Dockerfile).

This keeps the GUI fast and portable while making the analysis environment fully reproducible.

## Pipeline tools

Conda environments are defined in [`yamls/`](yamls/) and baked into the container image:

| Stage            | Tools |
| ---------------- | ----- |
| Assembly QC      | QUAST, CheckM2 |
| Annotation       | Bakta |
| Typing           | mlst, chewBBACA (cgMLST) |
| Pangenome        | Panaroo, DIAMOND |
| Phylogenetics    | Snippy, RAxML-NG, Gubbins, FastTree, MAFFT |
| Plasmids         | MOB-suite |

## Requirements

- **Python ≥ 3.11** on the host
- **Docker** (or **Singularity/Apptainer**) for the analysis tools
- Host Python dependencies (see [`requirements.txt`](requirements.txt)): PyQt6 + PyQt6-WebEngine,
  Polars, pandas, NumPy, PyArrow, SciPy, Matplotlib, PyQtGraph, Biopython

## Installation

```bash
git clone https://github.com/liviurotiul/BACTolomew.git
cd BACTolomew
bash install.sh
```

`install.sh` verifies prerequisites, creates a thin virtual environment (`.venv`) and installs
the host dependencies. The container image with the bioinformatics tools is built or pulled on
first launch.

To set up the host venv manually:

```bash
python -m venv .venv
.venv/bin/pip install -r requirements.txt
```

## Usage

### Annotation explorer (GUI)

```bash
.venv/bin/python -m annotation_explorer_polars.app
```

On first launch a **first-time setup** dialog seeds the reference databases and verifies the
tool environments, then writes a machine-local `bactolomew_config.json` (git-ignored).

### Pipeline (headless)

```bash
.venv/bin/python src/main.py --help
```

## Data & databases

Large inputs, reference databases and results live **outside** version control and are
git-ignored:

- `DB/` — reference databases (cgMLST schemes, marker/compound catalogs, Bakta DB)
- `runs/` — per-sample pipeline outputs
- `envs/` — prebuilt conda tool environments
- `cache/` — intermediate caches

Databases are seeded during first-time setup from a source root (default `/Data/BACTolomewDB`),
configurable in the setup dialog.

## Repository layout

```
src/
  main.py                        # pipeline entry point
  pipeline/                      # workflow orchestration, containerized tool dispatch (ToolEnv)
  annotation_explorer_polars/    # PyQt6 GUI: tables, trees, pangenome, plasmid & geo views
  scripts/                       # standalone analysis / plotting utilities
yamls/                           # conda environment definitions (one per tool group)
docker/                          # Dockerfile + build script for the tools image
map_tests/                       # geographic-map rendering backend experiments
install.sh                       # host installer
requirements.txt                 # host Python dependencies
```

## Configuration

`bactolomew_config.json` is **generated per-machine** by first-time setup (absolute paths to
envs, databases and the container runtime) and is intentionally not tracked in git.

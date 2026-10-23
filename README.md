# BACTolomew

BACTolomew takes a set of bacterial genome assemblies through QC, annotation, typing, AMR,
plasmid, pangenome and phylogenetic analysis, then lets you explore the results in a desktop
app: linked tables, trees, pangenome and plasmid views, gene neighborhoods and maps.

Every analysis tool runs from one container image, so the host needs only a small Python
app and a container runtime. The same `bactolomew` command opens the window on a
laptop or runs headless on an HPC node. Finished runs can be packed into one file and opened
by colleagues who only install the viewer.

![BACTolomew explore window, cgMLST minimum spanning tree, pangenome clustering, geographic map and gene neighborhoods](docs/images/bactolomew_overview.png)

A run of 153 *Acinetobacter baumannii* genomes: the explore window (phylogeny, gene
presence/absence, map, plasmidome), the cgMLST minimum spanning tree, pangenome clustering,
the geographic map and the *bla*<sub>OXA-23</sub> gene neighborhoods.

## What it does

The pipeline runs assembly QC, Bakta annotation, MLST, AMRFinderPlus (acquired genes and point
mutations), MOB-suite plasmids, a Panaroo pangenome, a recombination-aware core-genome
phylogeny (core alignment, FastBAPS lineages, Gubbins per lineage, IQ-TREE), a Snippy SNP
tree, chewBBACA cgMLST, and searches of databases you add yourself. An interrupted run can be
resumed from where it stopped (Resume in the run list, or `run --resume`), and a step is rerun
only when its inputs change.

The explorer is a set of dockable panels: an annotation table, a genome track, gene
neighborhoods and sequence alignment for looking at genes; a core-genome phylogeny, cgMLST
minimum spanning trees and pangenome clustering for population structure; a presence/absence
heatmap and a plasmidome view for gene content; and a map and summary charts for metadata. Every
panel follows the same filters and sample selection. You can spread the panels over several
windows and arrange them however you like; Window > Layouts has a few built-in layouts to start
from, and you can save your own there.

Metadata comes from a CSV with a `Sample` column, or from the run itself (species, MLST,
lineages from FastBAPS or fastGEAR, and assembly QC).

## Installation

BACTolomew can be installed in two ways, and a viewer can be turned into a full install later:

- **Full pipeline** (run analyses and view them): a container runtime (Docker, Podman or
  Apptainer/Singularity), 16 GB RAM (32 GB recommended) and the disk space described below.
- **Viewer only** (open runs others shared): about 1 GB of disk, no container runtime.

Installing takes two steps. `bash install.sh` installs the app (about 1 GB). The first time
you start it, unless you choose Viewer only, setup downloads about 50 GB: the tools image (about
7 GB compressed, about 20 GB once unpacked), the Bakta database (39.4 GB, unpacked to 76 GB)
and the CheckM2 database (1.7 GB). Before it starts, setup checks that the disk for the databases has
about 125 GB free, since Bakta needs room for the archive while it unpacks; the tools image
takes its 20 GB in the container runtime's storage on top of that. Afterwards about 100 GB stay
in use (Bakta 76 GB, the tools image 20 GB, CheckM2 3 GB), plus your runs.

`bash install.sh --viewer` (`-Viewer` on Windows) skips the Full or Viewer only choice at first
start and installs a viewer straight away. The installer checks what is missing and prints the exact
command to fix it.

On Linux (x86_64) the desktop window needs glibc 2.34 or newer (check with `ldd --version`).
Older systems, down to glibc 2.28, can still run analyses from the command line: install with
`bash install.sh --headless`, which leaves out the window (so it can't be combined with
`--viewer`). On 64-bit ARM Linux the window is not supported, and the command line has not
been tested there. macOS needs version 13 or newer for the window (12 for the command line
only).

### Ubuntu and Debian

Ubuntu 22.04 or newer, Debian 12 or newer.

```bash
sudo apt install git python3 python3-venv docker.io
sudo usermod -aG docker $USER
```

### Fedora

```bash
sudo dnf install git python3 moby-engine
sudo systemctl enable --now docker
sudo usermod -aG docker $USER
```

Fedora usually comes with Podman, which works just as well: if you use it, install only
`git python3` and skip the `systemctl` and `usermod` lines.

### Rocky, Alma and RHEL (8 and 9)

```bash
sudo dnf install git python3.12 podman
```

Docker from [Docker's repository](https://docs.docker.com/engine/install/centos/) works too.
On version 8 the window does not run; install with `bash install.sh --headless` and use the
command line.

### openSUSE Leap and Tumbleweed

```bash
sudo zypper install git python312 docker dejavu-fonts
sudo systemctl enable --now docker
sudo usermod -aG docker $USER
```

### Arch

```bash
sudo pacman -S --needed git python docker ttf-dejavu
sudo systemctl enable --now docker
sudo usermod -aG docker $USER
```

### Then, on any Linux

Log out and back in so the docker group takes effect (not needed with Podman), then:

```bash
git clone https://github.com/liviurotiul/BACTolomew.git && cd BACTolomew
bash install.sh
```

### macOS

1. Install Python 3.12 with [Homebrew](https://brew.sh) (`brew install python@3.12`) or from
   [python.org](https://www.python.org/downloads/macos/).
2. Install [Docker Desktop](https://docs.docker.com/desktop/setup/install/mac-install/) (the
   Intel or Apple silicon build to match your Mac), start it, and in Settings > Resources give
   it at least 16 GB of memory (32 GB if you can) and a disk limit of at least 150 GB.
3. Create a case-sensitive APFS volume (Disk Utility > Add Volume > "APFS (Case-sensitive)").
   Runs are saved inside the BACTolomew folder, and some tool outputs have names that differ
   only in case, so BACTolomew goes on that volume.
4. In Terminal:

   ```bash
   cd /Volumes/<volume name>
   git clone https://github.com/liviurotiul/BACTolomew.git && cd BACTolomew
   bash install.sh
   ```

On Apple Silicon (M1 to M4) the published tools image is for Intel only, so the first setup
builds a native one on your Mac (about an hour). fastGEAR, an optional alternative for lineage
detection, does not run there, but FastBAPS (the default) does. Running the Intel image under
Rosetta is also possible but untested; see [docs/INSTALL.md](docs/INSTALL.md).

### Windows 10 and 11

BACTolomew runs inside WSL2 (the Linux system built into Windows), and its window opens on the
Windows desktop. In PowerShell, not as administrator:

```powershell
Invoke-WebRequest -Uri https://raw.githubusercontent.com/liviurotiul/BACTolomew/main/install_windows.ps1 -OutFile install_windows.ps1
powershell -ExecutionPolicy Bypass -File .\install_windows.ps1
```

The script sets up WSL2 and Ubuntu, then Docker (it asks whether to use Docker Desktop or Docker
inside Linux), installs BACTolomew and adds it to the Start Menu. If it asks you to restart
Windows, restart and run it again; it carries on.

### HPC and servers

```bash
git clone https://github.com/liviurotiul/BACTolomew.git && cd BACTolomew
bash install.sh --headless
./bactolomew setup --yes --runtime apptainer      # or docker, podman
```

With Apptainer, setup converts the tools image to a `.sif` file (about 15 GB more). The
Apptainer route has not been tested on a real cluster yet, so try a small run first. See [Command line and HPC](#command-line-and-hpc) below for running jobs.

### After installing

Start BACTolomew with `./bactolomew` in its folder, or from the Start Menu on Windows, then
follow [Getting started](#getting-started) below.

The Linux installs above are tested; macOS and Windows have not been tested on a real machine
yet. Minimal systems may need extra packages (the installer names them).
[docs/INSTALL.md](docs/INSTALL.md) has the full guide, common problems and how to uninstall.

## Getting started

### Set up the full pipeline

1. Install as above (`bash install.sh`).
2. Run `./bactolomew` and choose **Full (run analyses)**. On a server without a screen, run
   `./bactolomew setup --yes` instead.
3. Setup first checks that there is enough free disk space (see
   [Installation](#installation) for the sizes), then downloads each of these and shows its
   progress:
   - the tools image, `liviurotiul/bactolomew` from Docker Hub (built locally on Apple Silicon);
   - the full Bakta database (the largest, about 1 to 2 hours), which includes the
     AMRFinderPlus database;
   - the CheckM2 database.

   It finishes with a quick check that everything works. If something fails, click Retry and
   setup carries on from there.
4. When setup finishes, the [run list](#the-run-list) opens. (From a terminal,
   `./bactolomew check` now prints "Ready to run.")
5. Click **Start New Run**, then **Add Assemblies…** and pick one FASTA file per genome
   (`.fasta`, `.fa` or `.fna`, gzipped or not). Each file name, minus the suffix, becomes the
   sample name. Names may only contain letters, digits and `_ . - |`, and two names can't
   differ only in case.
6. Click OK to start. The run is saved under `runs/` in the BACTolomew folder. You don't need
   to give the species: MLST works it out, and if cgMLST.org has a scheme for it, the scheme is
   downloaded during the run (the first time, you'll be asked to accept the cgMLST.org terms;
   see [cgMLST schemes](#cgmlst-schemes) for species without one). **Monitor** shows each step
   as it runs.
7. When the run is done, click **View**. BACTolomew asks for metadata: pick a CSV whose
   `Sample` column matches the sample names, or choose **Use run data only**. The other
   columns become filters and colors; a `Country` column (or `Latitude` and `Longitude`)
   places samples on the map. **Change Metadata...** in the run's **More** menu swaps the CSV
   later.

### Set up a viewer

1. Install with `bash install.sh --viewer` (on Windows:
   `powershell -ExecutionPolicy Bypass -File .\install_windows.ps1 -Viewer`).

   If you already ran plain `bash install.sh`, skip the reinstall: at first start choose
   **Viewer only (open results from others)**, or run `./bactolomew setup --viewer` on a
   machine without a screen.
2. Run `./bactolomew` (on Windows, BACTolomew in the Start Menu). It opens on the run list.
   (To confirm from a terminal, `./bactolomew check` prints "Ready to view runs.")
3. Click **Import Run...** and pick the `.zip` you were sent (or a run folder). A run made with
   a newer version shows a warning but can still be opened (**Open Anyway**).
4. Click **View** and pick a metadata CSV whose `Sample` column matches the sample names, or
   choose **Use run data only** (see step 7 of the full pipeline setup).

In a viewer, the few features that need the tools are grayed out, with a note saying why:
aligning sequences from a gene card (unless MAFFT is installed on the computer), adding
cgMLST schemes or databases, and starting or resuming runs.

### Share a run

On the machine that made the run, open its **More** menu in the run list and choose **Export
for Sharing...**. This writes one `.zip` with what the viewer needs and no paths from your
computer: about 0.5 GB for 150 genomes, plus about 4 MB per genome for the assemblies. The
assemblies and the metadata CSV are included by default; untick the assemblies to make the
file smaller. The other person can import the `.zip` as described above.

### Turn a viewer into a full install

Click **Set Up Pipeline…** in Settings (the gear button in the run list), or run
`./bactolomew setup`. It installs everything listed for the full pipeline.
`./bactolomew check` always tells you which mode you are in.

## The run list

The run list is the first window. Each run there has View, Monitor and Resume buttons and a
**More** menu. From here you start, monitor, resume and stop runs, open
finished ones, and import runs that others have shared with you. Each run's **More** menu has
the less common actions: run details, changing its metadata, exporting it and showing it in the
file manager.

The gear button opens Settings: compute limits, the lineage method (FastBAPS by default, or
fastGEAR, which is not available on Apple Silicon; `phylo.population_method` in the config),
cgMLST schemes and your own databases.

## Command line and HPC

`./bactolomew` (or `./bactolomew gui`) opens the window. Everything else works without a
display:

```bash
./bactolomew check                                   # checks the install, shows the cores and RAM a run would use
./bactolomew setup --yes [--runtime docker|podman|singularity|apptainer] [--sif PATH] [--viewer] [--dry-run]
./bactolomew run --genomes DIR --out DIR [--threads N] [--resume]
./bactolomew db --list
./bactolomew db --add proteins.faa --name MyDB --category AMR
./bactolomew db --remove MyDB
./bactolomew --version
```

- `setup` does the same first-time setup as the window and prints its progress as text.
  `--yes` answers yes to every prompt. Without `--runtime`, it uses the runtime already in the
  config, or else the first of Docker, Podman and Apptainer/Singularity that it finds. `--sif PATH` sets where the
  Apptainer/Singularity image lives (and implies `--runtime singularity`): if there is already
  an image at that path, setup uses it, and if not, setup pulls the tools image from Docker Hub
  and saves it there as a `.sif`.
- `run` takes a `--genomes` folder of FASTA files, one per genome, named as in
  [step 5](#set-up-the-full-pipeline) above. It takes the database paths from
  `bactolomew_config.json` (`run --help` lists the overrides) and refuses an output folder that
  already holds a run unless you pass `--resume`.
- Before the first run that needs a new cgMLST scheme, accept the cgMLST.org terms once by running, from the
  BACTolomew folder, `PYTHONPATH=src .venv/bin/python -m pipeline.CgmlstInstall --accept-terms`.
  Otherwise the run stops at the chewBBACA step; accept the terms, then continue with `run --resume`.
- The number of cores and the RAM a run uses follow the CPU affinity, cgroup limits and the
  Slurm allocation (`SLURM_CPUS_PER_TASK`, `SLURM_MEM_PER_NODE` or `SLURM_MEM_PER_CPU`), capped
  by `resources.max_threads` and `resources.max_ram_gb` in the config. The budget is logged at
  the start of each run.
- On a cluster without Docker, use Apptainer/Singularity (`setup --runtime apptainer`; see
  [HPC and servers](#hpc-and-servers)).

Setup writes `bactolomew_config.json` in the BACTolomew folder: database paths, container
runtime and image, and resource limits. Settings in the window edits the same file.

On Slurm, for example (submitted from the BACTolomew folder):

```bash
sbatch --cpus-per-task=32 --mem=64G --wrap "./bactolomew run --genomes /path/to/assemblies --out /path/to/run"
```

To look at a headless run, copy its run folder to any computer with BACTolomew (a viewer is
enough) and open it with **Import Run...**.

## Gene categories

Genes are put into categories (AMR, virulence, mobile elements, plasmid) using the tools the
pipeline already runs.

- **AMR**, point mutations, stress/metal/biocide resistance (BMRG) and the virulence genes
  that AMRFinderPlus reports with its `--plus` option come from AMRFinderPlus, with the AMRFinderPlus database that ships inside
  the Bakta database.
- **Virulence** and **mobile elements** come from Bakta's VFDB (virulence factor) and
  insertion sequence (IS) hits, plus Bakta product names for mobile elements and phages
  (marked as product name matches).
- MOB-suite's replicon, relaxase and mating pair formation (MPF) hits mark genes as
  **Plasmid**.

## Your own databases

BACTolomew comes with no marker databases of its own; the categories above come from the tools
the pipeline runs. To add one, in Settings > Databases click "Add database…", or run `./bactolomew db --add`. Give a protein
FASTA (plain or `.gz`) and a category: AMR, BMRG, VIRULENCE, MOBILE, PHAGE, PLASMID, or
COMPOUND for a locus of several genes (an operon or a resistance island, for example) that is
matched as one element. From the next run on, it is searched with DIAMOND and its hits appear
in the tables as "<name> (U)", where U marks a database you added. Check the license of any
database you add; the Add database dialog has a note field for it.

## cgMLST schemes

In Settings, "Download from cgMLST.org…" fetches a scheme directly, and "Add schema…" imports
one or more schemes you downloaded yourself. There is no `bactolomew` subcommand for schemes
yet; from the BACTolomew folder run:

```bash
PYTHONPATH=src .venv/bin/python -m pipeline.CgmlstInstall --list    # also --status, --install <species>, --install-all, --import <paths>
```

- cgMLST.org schemes are for non-commercial use only and may not be redistributed. The terms
  are shown once, before the first download.
- If your species has a cgMLST.org scheme that is not installed yet, the run installs it at
  its chewBBACA step. The window asks you to accept the terms before such a run starts; a
  command-line run stops at that step until you accept them (see
  [Command line and HPC](#command-line-and-hpc)).
- If your species has no scheme at all, the run stops at the chewBBACA step until you import
  one ("Add schema…" or `--import`); Resume then continues.
- chewBBACA adds the new alleles it finds to the installed scheme, so allele numbers stay
  comparable between runs. To start over from a fresh download, use
  `--install <species> --replace`; the old scheme is moved to `.replaced/` in the schemes
  folder.

## Outputs

A run folder (under `runs/` by default) has one subfolder per step (mostly named after its tool), and
the trees in `trees/`. The files most people look for:

- Annotations (GFF3, GenBank, protein FASTA): `bakta/<sample>/`
- AMR calls: `amrfinder/<sample>.tsv`
- Gene presence/absence: `panaroo/out/gene_presence_absence.csv`
- Trees (Newick): `trees/`, listed in `trees/trees.json`
- cgMLST allele profiles: `chewbbaca/out/results_alleles.tsv`

Three JSON files describe the run:

- `run_config.json`: the settings and the environment (image, tool versions, databases)
- `run_state.json`: step progress and warnings
- `run_report.json`: the provenance of the finished run

## Tools and versions

All tools come from the `liviurotiul/bactolomew` image on Docker Hub. The versions below are from
[`docker/locks/tool_versions.json`](docker/locks/tool_versions.json).

| Stage              | Tools |
| ------------------ | ----- |
| Assembly QC        | QUAST 5.3.0, CheckM2 1.1.0 |
| Annotation         | Bakta 1.10.4 |
| AMR                | AMRFinderPlus 4.0.23 (`--plus`, point mutations with `--organism` from the MLST species) |
| Typing             | mlst 2.35.0 (PubMLST data of 2025-12-29), chewBBACA 3.5.4 (cgMLST) |
| Plasmids           | MOB-suite 3.1.9 |
| Pangenome          | Panaroo 1.8.0 |
| Core phylogeny     | FastBAPS 1.0.8 (lineages, default), Gubbins 3.4.3, IQ-TREE 2.4.0; fastGEAR optional (x86_64 only) |
| SNP phylogeny      | Snippy 4.6.0, RAxML-NG 1.2.2 |
| Database searches  | DIAMOND 2.2.3 |
| Sequence alignment | MAFFT 7.526 (explorer) |

## Development

To run from a checkout without `install.sh`:

```bash
git clone https://github.com/liviurotiul/BACTolomew.git
cd BACTolomew
python3 -m venv .venv                        # Python 3.10 to 3.14
.venv/bin/pip install --require-hashes --only-binary :all: -r requirements.lock
cd src && ../.venv/bin/python -m annotation_explorer_polars.app
```

The tests use unittest and run offscreen:

```bash
QT_QPA_PLATFORM=offscreen PYTHONPATH=src .venv/bin/python -m unittest discover -s tests -p 'test_*.py'
```

The tool environments are defined in [`yamls/`](yamls/) and pinned exactly in
[`docker/locks/`](docker/locks/). Build the tools image with `bash docker/build.sh [version]`;
publishing it is described in [`docker/PUBLISHING.md`](docker/PUBLISHING.md).

## Citation

If you use BACTolomew, cite the release you used (`./bactolomew --version` prints it), for
example "BACTolomew v0.3.0, https://github.com/liviurotiul/BACTolomew/releases/tag/v0.3.0",
together with the tools listed in [Tools and versions](#tools-and-versions) and the databases
they use (Bakta, AMRFinderPlus, CheckM2, MOB-suite, PubMLST), which have their own licenses and
citations.

## License

BACTolomew is free software under the GNU General Public License version 3. See
[`LICENSE`](LICENSE).

# Publishing the tool image

This is the checklist for publishing the tools image `liviurotiul/bactolomew` to Docker Hub, and (last section)
for minting the code DOI. None of it has been run yet: no image is pushed, no registry login
exists, and this machine's `bactolomew_config.json` still names the local `bactolomew-tools:0.1.0`.

Nothing else is published: there is no database bundle. First-time setup downloads the full
Bakta database (with its AMRFinderPlus database) and CheckM2 from their own sources, and
makes `DB/` with empty catalogs and the shipped cgMLST scheme species catalog.
(`tools/zenodo_upload.py` was written for a database bundle and is not needed for this.)

## Version and tag

Decided 2026-10-21: the first public release is 0.2.0. At commit `d0951e9`, `__version__`
in `src/pipeline/__init__.py`, `TOOLS_IMAGE` in `src/pipeline/ToolEnv.py` and
`src/scripts/demo_install.py` name 0.2.0; tag the GitHub release `v0.2.0`. The working tree
names 0.3.0.dev0 and `bactolomew:0.3.0`.

**The 0.2.0 image is built from commit `d0951e9`** (pushed to `origin/main`; `docker/` last
changed in `8142861`). The working tree's `docker/` now holds the 0.3.0 changes (FastBAPS,
`population_methods`, the per-platform locks with the common tool configuration below), so a
build from it is not the tested 0.2.0: `docker/build.sh` defaults to `0.3.0-dev`, and steps
2 and 3 below run in a separate worktree of `d0951e9`. The in-app fallback build
(`container.build_on_missing`, `ContainerRuntime.build_image_local`) builds the working
tree's Dockerfile under the configured tag; the working tree's `TOOLS_IMAGE` is
`bactolomew:0.3.0`, so it never puts new content under the 0.2.0 tag unless a
config names 0.2.0 and that image is missing locally.
Two local images exist: `bactolomew-tools:0.1.0` (29.3 GB, the image the reference run used)
and `bactolomew-tools:0.2.0` (22.2 GB, the tested rebuild from the locks). Until a published
image is configured, setup builds the image named by `TOOLS_IMAGE` locally from the current
Dockerfile (`container.build_on_missing` defaults to true).

## What is already in place

| Piece | Where |
| --- | --- |
| Exact package locks per platform and env (from 0.3.0; explicit `url#md5` lists, plus PyPI wheels and R sources pinned by sha256 where an env needs them) | `docker/locks/linux-64/`, `docker/locks/linux-aarch64/` |
| Solving the yamls for both platforms with the same versions (throwaway miniforge container, `mamba --dry-run --platform`, no emulation), the lock check, the summary | `docker/update_locks.sh`, `docker/solve_locks.py`, `docker/check_locks.py`, `docker/make_tool_versions.py` |
| Readable tool version summary of the locks (its notes are copied into the image) | `docker/locks/tool_versions.json` |
| Dockerfile installing from the target platform's locks (picked from `TARGETARCH`), base images pinned by digest, OCI labels | `docker/Dockerfile` |
| Env install and trim step (one layer per env group; pip wheels and R sources after the conda packages) | `docker/install_envs.sh` |
| Build-time checks: the population structure methods present (FastBAPS required, fastGEAR when included), the MATLAB runtime system libraries with fastGEAR, mob_suite database files equal to recorded sha256s; writes `/opt/bactolomew/image_contents.json` (with `population_methods`), `apt_packages.tsv` and `tool_versions.json` (the envs really installed) | `docker/verify_image.sh`, `docker/locks/mobsuite_db.sha256` |
| Build script (`--platform`, default the Docker engine's own, version argument, labels, fastGEAR only on amd64 or `FASTGEAR=none`, the label `org.bactolomew.population_methods`, failing unless the methods are `fastgear,fastbaps` on amd64 and `fastbaps` otherwise, no credentials in labels, fastgear mirror override, prints the image ID) | `docker/build.sh` |
| Checking a built image against its platform's locks | `docker/update_locks.sh --check-image <image>` |
| Tool runs on a tiny input set in two images, outputs compared (both here, or one image per machine and a compare) | `docker/smoke_compare.sh` |
| Build context limited to the files the Dockerfile copies (about 260 KB) | `.dockerignore` |
| Digest pinning in `container.image` (`name:tag@sha256:...`), with a warning on mismatch | `src/pipeline/ContainerRuntime.py` |
| Image acquisition: pull any reference that names a registry (Docker and Singularity), build locally only when `container.build_on_missing` is true; pulls and builds use the engine's own platform (`linux/amd64` or `linux/arm64`, `engine_platform`), and only an amd64 build includes fastGEAR | `src/pipeline/ContainerRuntime.py` |
| Run provenance: `environment` (image, tool versions, databases) and `environment_history` in `run_config.json` and `run_report.json` | `src/pipeline/Runner.py` |

Why a new version instead of pushing 0.1.0: the local 0.1.0 was built from an older
Dockerfile (label `contexar-tools`, yamls under `/opt/contexar/yamls`, no `libncurses.so.5`
link), so the current Dockerfile does not describe it, and it has no
`/opt/bactolomew/tool_versions.json`. Version 0.2.0 is built from the locks, so it carries
the same tool versions as 0.1.0 (see `tool_versions.json`), with a matching Dockerfile.

## Population structure methods (from 0.3.0)

From 0.3.0 the image always holds FastBAPS, and on linux/amd64 fastGEAR too; the pipeline picks one per run
(`phylo.population_method`: fastbaps (the default), fastgear or auto, which takes FastBAPS
when the image has it, else fastGEAR; see `src/pipeline/PopulationStructure.py`). Left
unset, the default falls back to fastGEAR on an image without FastBAPS (0.2.0 and older);
a method set in the config or for the run is required. Runs made before FastBAPS, and runs
that recorded fastGEAR, keep fastGEAR when resumed.

- FastBAPS is the `fastbaps` env on every platform: R 4.4.3 and fastbaps's dependencies
  from conda, `genie` 1.0.7 (CRAN) and `fastbaps` 1.0.8 (GitHub tag `v1.0.8`) built from
  source by `install_envs.sh` (`docker/locks/<platform>/fastbaps.r_source.tsv`, checked by
  sha256; each tarball lists a second URL, a CRAN mirror or GitHub's codeload, in case
  one host is down), because bioconda's r-fastbaps 1.0.8 and r-genie are linux-64 only. As
  the default method it is required: a failed download or source build fails the build,
  so no image without FastBAPS is built from these files. Its layer comes after the other
  envs.
- fastGEAR comes from the fastgear image (MATLAB runtime, x86_64 only), so it is included
  only for `linux/amd64` (`FASTGEAR=none` leaves it out there too). It is the only part of
  the image that differs by platform.
- `verify_image.sh` records the methods present in `/opt/bactolomew/image_contents.json`
  (`"population_methods": ["fastgear", "fastbaps"]`) and fails without FastBAPS;
  `build.sh` fails unless they are `fastgear,fastbaps` (linux/amd64) or `fastbaps`
  (linux/arm64, `FASTGEAR=none`), then adds them as the label
  `org.bactolomew.population_methods` (a metadata-only image on top of the build). Before
  releasing an amd64 image, check that the build ended with
  `[build] done: ... (population methods: fastgear,fastbaps)`. The pipeline reads the label, else the file,
  else (images before 0.3.0) takes fastGEAR when its files are in the image.

## One tool configuration on every platform (from 0.3.0)

Decided 2026-10-22: linux/amd64 and linux/arm64 use the same configuration wherever
possible; where linux-aarch64 is more restrictive, its configuration is used on both.
`docker/check_locks.py` fails on any version difference between `docker/locks/linux-64` and
`docker/locks/linux-aarch64` other than the platform-only packages it lists by name (a
metapackage and packaging splits, no tool). What this changes against 0.2.0 (on x86 too):

- **QUAST and CheckM2 in two envs**, `quast` (QUAST 5.3.0, Python 3.13) and `checkm2`
  (CheckM2 1.1.0 with TensorFlow 2.18.0, Python 3.12): on aarch64 quast 5.3.0 needs Python
  3.13 and checkm2's TensorFlow builds stop at Python 3.12. scikit-learn is pinned at
  1.6.1 (as in 0.2.0), the version CheckM2's pickled MinMaxScaler was saved with: a newer
  one warns (`InconsistentVersionWarning`) that the model may give invalid results. The
  pipeline calls each tool in its own env and falls back to the old `assembly_qc` env in
  images before 0.3.0 (`ToolEnv.LEGACY_TOOL_ENVS`).
- **BLAST 2.16.0 everywhere** (bakta, mobsuite, chewbbaca, mlst, quast): the only version
  built for aarch64. trnascan-se stays 2.0.12, raxml-ng stays 1.2.2 (phylogeny) and 2.0.2
  (gubbins), all built for both.
- **Panaroo 1.8.0 without its dependencies**: intbitset has no aarch64 conda build, so
  the env holds the recipe's other dependencies from conda, `intbitset` 4.1.2 from PyPI
  (manylinux wheels for both, sha256 in `panaroo.pip.txt`) and the panaroo package itself
  (`x-bactolomew: no_deps` in `yamls/panaroo.yaml`). prokka is no longer in the env (the
  pipeline never used it; panaroo 1.8.0's current build does not ask for it).
- **Snippy 4.6.0 without its dependencies**: its recipe pins vcflib <= 1.0.2 and tabixpp
  1.1.0, not built for aarch64. The env holds vcflib 1.0.9 (it still has the vcf programs
  snippy calls) and `perl-bioperl-core` in place of the perl-bioperl bundle (whose aarch64
  build pulls samtools 0.1.19), with samtools 1.20, bwa 0.7.18 and minimap2 2.28.
- **mlst 2.35.0 without its dependencies**: its recipe asks for perl-bioperl, which mlst
  never loads and which pulls samtools 0.1.19 built against different zlib and OpenSSL
  versions per platform.
- **FastBAPS built from source** on both (above), with its R dependencies genie 1.0.7,
  genieclust 1.2.0 and quitefastmst 0.9.2 (CRAN sources checked by sha256, see
  `yamls/fastbaps.yaml`), because their bioconda/conda-forge builds are linux-64 only.

`bactolomew:0.3.0` (2026-10-22, local only, linux/amd64, image ID
`sha256:b918a55d6c3e`, 19.5 GB, fastGEAR and FastBAPS) is built from these locks
(`update_locks.sh --check-image`: 0 differences). Against 0.2.0 on the smoke inputs (step 4)
every output is the same: Snippy, RAxML-NG, IQ-TREE, Gubbins (tree and recombination),
Panaroo, chewBBACA, mlst and QUAST; the pipeline's assembly QC step (QUAST and CheckM2,
TensorFlow 2.18.0 against 2.17.0, scikit-learn 1.6.1 in both, so no `InconsistentVersionWarning`)
gave the same summary on two genomes. FastBAPS built
from source gave the same clusters as bioconda's linux-64 binary, on the smoke alignment
and on the 153-genome reference core alignment. The image is 2.7 GB smaller than 0.2.0
(prokka and its dependencies left the panaroo env). linux/arm64 is not built yet (below);
its locks were checked only by solving and by fetching every package URL (all served).

The yamls hold the intent (tools and key versions, the `no_deps` packages, the pip and R
source pins in an `x-bactolomew` block that conda ignores). Native (non-container)
installs create envs from the yamls with `conda env create`, which does not read that
block: they lack panaroo, snippy, mlst and fastbaps. The container image is the supported
way to run the pipeline.

## Size

0.1.0 is 29.3 GB. The Dockerfile applies only the reductions that cannot change what a
tool does:

- envs grouped per layer so shared packages are stored once (about -5 GB):
  `quast checkm2 bakta gubbins phylogeny mafft` (about 6.8 GB), `panaroo chewbbaca mlst diamond`
  (about 7.2 GB), `mobsuite` alone (about 5.1 GB with its database). Every layer stays
  under 10 GB uncompressed;
- the conda package cache removed in the same layer (about -2 GB);
- docs (`share/doc`, `share/man`, `share/info`, `share/gtk-doc`) and the static libraries
  directly in each env's `lib/` (about -0.5 GB before sharing).

Expected result: roughly 21 to 22 GB uncompressed. Measure it after the build
(`docker image ls`, `docker history`).

Deferred on purpose, each needs a smoke test or a behavior decision first:

- mob_suite database (2.9 GB): the mob_suite package downloads it at install time
  (`mob_init`, see "The mob_suite database" below). Moving it out of the image means
  `mobsuite.database_directory` must be set (route B).
- prokka and barrnap in the panaroo env (about 2.2 GB): unused by the pipeline, but
  removing them needs a panaroo run on a small set.
- compilers and sysroot that R and k8 pull into gubbins, quast, phylogeny and fastbaps
  (about 1 GB): needs a gubbins and quast/minimap2 smoke test.
- `share/locale`, `include/`, test folders inside site-packages. Never delete `*.pyc`.

## Steps to publish (in order)

### 1. The registry

Docker Hub, repository `liviurotiul/bactolomew` (public; anonymous pulls are rate limited,
about 100 per 6 hours per IP, which is plenty for installs). The fastGEAR base image is
mirrored as `liviurotiul/fastgear:2022-12-01`. Published so far: `0.2.0` (release v0.2.0,
commit 8142861) and `0.3.0`.

The steps below use `IMAGE=liviurotiul/bactolomew`.

### 2. Build 0.2.0 from its commit

The local test image `bactolomew-tools:0.2.0` was built from an uncommitted tree: its
`org.opencontainers.image.revision` label is `d11b88d...-dirty`, which names no commit. Do
not publish it. Its `docker/` was then committed and pushed as part of `8142861`
(`origin/main` is at `d0951e9`). Do not build 0.2.0 from the working tree: it holds the
0.3.0 changes (FastBAPS, see "Population structure methods"). Build it from a detached
worktree of `d0951e9` (no branch is created; build.sh labels the image with that clean
commit), **without pruning the Docker build cache first** (no `docker builder prune`, no
`docker system prune`). BuildKit then reuses the cached apt and env layers of the tested
0.2.0, so the tools and the mob_suite database are byte-identical to what was smoke-tested.
A pruned cache fetches the apt packages and the mob_suite database again, with different
content (see below).

```bash
git worktree add --detach ../BACTolomew-0.2.0 d0951e9
cd ../BACTolomew-0.2.0
bash docker/build.sh 0.2.0 liviurotiul/bactolomew
docker image ls liviurotiul/bactolomew
docker image inspect --format '{{json .Config.Labels}}' liviurotiul/bactolomew:0.2.0
```

The revision label must read `d0951e9...` without `-dirty`, and there must be no
`org.bactolomew.population_methods` label (that comes from 0.3.0's build.sh).

The Dockerfile at `d0951e9` differs from the one the tested 0.2.0 was built with only after
the env layers, so the rebuild redoes only the layers from the `tool_versions.json` copy on
(that file gained a note; the fastGEAR copy that follows is the same content): the
image-wide `LD_LIBRARY_PATH` is gone (the MATLAB runtime's 2016 libstdc++, libgcc_s, libssl
and others no longer shadow the system and conda libraries; fastGEAR's launcher sets its
own path, and fastGEAR on its example data gave identical output with and without the
variable), and `docker/verify_image.sh` runs as a new layer. Check the build log shows
`CACHED` for the `install_envs.sh` steps; if it does not, the env layers were rebuilt and
step 4 must be repeated in full.

The build needs network access to conda-forge, bioconda and Docker Hub (and Zenodo and
NCBI when the mob_suite layer is rebuilt). It always builds `linux/amd64` (the locks and
the MATLAB runtime are x86_64 only), also on an arm64 host. It does not touch
`bactolomew-tools:0.1.0`.

### 3. Check the build matches the locks

In the same worktree (its `update_locks.sh` and `make_tool_versions.py` are the 0.2.0 ones;
the working tree's version would treat 0.2.0 as a partial export, since 0.2.0 has no
fastbaps env, and would not give the result below):

```bash
bash docker/update_locks.sh liviurotiul/bactolomew:0.2.0
git diff --stat -- docker/locks
```

The `.lock` files must be unchanged; only `locks_exported_from` in `tool_versions.json`
changes. Discard that change in the worktree (`git checkout -- docker/locks/tool_versions.json`);
step 6 records the published image where it belongs.

This check reads `conda-meta`, so it holds by construction (the envs were installed from
those locks) and cannot see the files `install_envs.sh` deleted after install (static
libraries in `lib/`, docs). It shows the right packages went in, not that the tools still
work: step 4 is the evidence for that.

Steps 4 to 7 use the image only, so they can run from either tree. Remove the worktree once
0.2.0 is pushed: `git worktree remove ../BACTolomew-0.2.0`.

### 4. Smoke test

```bash
I=liviurotiul/bactolomew:0.2.0
for c in "bakta bakta --version" "bakta amrfinder --version" "diamond diamond version" \
         "mobsuite mob_recon --version" "phylogeny snippy --version" "phylogeny raxml-ng --version" \
         "gubbins run_gubbins.py --version" "gubbins iqtree2 --version" "panaroo panaroo --version" \
         "chewbbaca chewBBACA.py --version" "mlst mlst --version" "assembly_qc quast.py --version" \
         "assembly_qc checkm2 --version" "mafft mafft --version"; do
    set -- $c; docker run --rm --network none "$I" conda run -n "$1" "${@:2}" 2>&1 | head -1
done
docker run --rm "$I" cat /opt/bactolomew/tool_versions.json | head
```

Then run the heavier tools on a tiny input set in 0.1.0 and the new image and compare the
outputs (Snippy, RAxML-NG, IQ-TREE, Gubbins, Panaroo, chewBBACA CreateSchema and AlleleCall;
each container runs with `--network none`):

```bash
bash docker/smoke_compare.sh bactolomew-tools:0.1.0 "$I" <inputs-dir> <empty work-dir>
```

The script's header says what the inputs folder holds (four Bakta GFF3s and their contigs,
a Snippy reference, a small core alignment). Every line should read `SAME`.

Done on 2026-10-13 for 0.1.0 against the tested 0.2.0, on four ST492 genomes (Bakta GFF3s
and contigs copied from the reference run, Snippy reference GCA_028622775.1) and a
12-sequence, 400 kb core alignment with about 10,700 variable sites: all seven outputs
identical (6,090 Snippy SNPs, 6 Gubbins recombination blocks, 3,838 Panaroo gene families,
chewBBACA alleles for 4 genomes). fastGEAR on its bundled example (`/opt/data3_many`) gave
identical lineage and recombination files with and without the old image-wide
`LD_LIBRARY_PATH`. Not yet done: a pipeline run on a few genomes with the new image
(below), QUAST/CheckM2 and Bakta beyond `--version`.

Then run the pipeline on a few genomes with the new image: back up
`bactolomew_config.json`, set `container.image` to the new tag in a copy or temporarily,
run into a fresh output folder (never into the reference run), check
`run_config.json` -> `environment` names the image and its tool versions, and restore the
config.

### 5. Log in and push

```bash
# A Docker Hub personal access token with Read & Write access
# (hub.docker.com > Account settings > Personal access tokens)
docker login -u liviurotiul          # paste the token as the password
docker push liviurotiul/bactolomew:0.3.0
```

A new Docker Hub repository is public by default; check it under Repositories.

### 6. Record the digest

```bash
docker image inspect --format '{{index .RepoDigests 0}}' liviurotiul/bactolomew:0.2.0
# liviurotiul/bactolomew@sha256:<digest>
```

Write the digest into the README and the release notes, and into the main tree's
`docker/locks/tool_versions.json` (`locks_exported_from`: the locks are the 0.2.0 ones;
the fastbaps entry keeps its own `exported_from`). Optionally tag the commit the image was
built from (`git tag tools-0.2.0 d0951e9`).

### 7. Check a clean pull

On another machine, or after `docker image rm` of the pushed tag only:

```bash
docker pull liviurotiul/bactolomew:0.2.0@sha256:<digest>
```

### 8. Singularity / Apptainer (HPC)

With `runtime: "singularity"`, `container.image` set to the full reference (step 9) and
`container.sif_path` set to where the `.sif` should go, the app pulls it on first start
(`singularity pull <sif_path> docker://<image>`); no Docker daemon is needed. By hand:

```bash
singularity pull bactolomew_0.3.0.sif \
    docker://liviurotiul/bactolomew:0.3.0@sha256:<digest>
# or, from a local Docker image:
singularity build bactolomew_0.3.0.sif docker-daemon://liviurotiul/bactolomew:0.3.0
```

The `.sif` is about the compressed size of the image; build it on a filesystem with room
for twice that (`SINGULARITY_TMPDIR`). Set `container.sif_path` to it. The Singularity
path is untested here (no Singularity on this machine): try it once on the HPC before
announcing it.

### 9. Point the app at the published image

`bactolomew_config.json` (per machine, not in git):

```json
"container": {
  "image": "liviurotiul/bactolomew:0.2.0@sha256:<digest>",
  "registry": "",
  "sif_path": "",
  "build_on_missing": false
}
```

Use the full reference in `image` and leave `registry` empty; `registry` only prefixes a
bare name (one without `/`). Any reference with a `/` is pulled as it is, by Docker and by
Singularity. With a digest in `image`, the app pulls exactly that image.

`build_on_missing` decides what happens when the image cannot be pulled (network down, a
package still private, a typo in the digest) or names no registry: `false` stops with an
error that says why; `true` builds the image locally from `docker/Dockerfile` (about an
hour, downloads the mob_suite database again, and the result does not match the pinned
digest; the app then warns and uses it). Use `false` for installs of the published image,
`true` only for development. With `true` and a local build of the same name and tag, the
app warns that the digest differs and uses the local one.

In the code, bump the default image for new installs:

- `TOOLS_IMAGE` in `src/pipeline/ToolEnv.py` (the full registry reference, so new installs
  pull it), and `__version__` in `src/pipeline/__init__.py` if the app version follows the
  image;
- the image in `src/scripts/demo_install.py`;
- the README ("First launch and setup": the image is pulled, the local build is only for
  development; "Uninstall": the image name).

## Changing a tool version later

1. Edit `yamls/<env>.yaml` (pin the tool and the dependencies whose version matters).
2. `bash docker/update_locks.sh <env> ...` (no argument: every env). It solves the yamls
   for linux-64 and linux-aarch64 in a throwaway `condaforge/miniforge3` container (the
   same base as the Dockerfile, `mamba create --dry-run --platform`: nothing is installed,
   no emulation), pins both platforms to the same versions, writes both lock folders, runs
   `docker/check_locks.py` and rewrites `tool_versions.json`. A difference the solver
   cannot remove makes the check fail: pin a version both platforms have in the yaml (or
   use `no_deps` as panaroo, snippy and mlst do), or, for a package that only exists on one
   platform and holds no tool, add it to `PLATFORM_ONLY` in `check_locks.py` with the reason.
3. Review `git diff -- docker/locks yamls`.
4. Build with `docker/build.sh <new version>`, check it with
   `bash docker/update_locks.sh --check-image <image>` (the installed envs equal the
   platform's locks), run the smoke test (step 4) against the previous version, and
   publish as above.

## Building on Apple Silicon (arm64)

Not done yet: this machine is x86_64 and does not emulate arm64. The arm64 image is
built natively on an Apple Silicon Mac with Docker Desktop, from the same commit as the
x86 image, and its outputs are compared with the x86 image's on the same tiny inputs.

**Docker Desktop settings** (Settings > Resources): memory at least 16 GB (Docker Desktop
gives its Linux VM only part of the Mac's RAM by default, often half; the build needs about
8 GB, Bakta and Panaroo on real data more), swap 2 GB or more, and a virtual disk limit of
at least 120 GB (the image is about 20 GB and the build cache as much again; the smoke test
needs little). Settings > General: "Use Virtualization framework" on (the default); Rosetta is
not needed for the arm64 build (it would only run the x86 image).

**Build** (in a clone of the repository at the commit to build; nothing is pushed):

```bash
docker version --format '{{.Server.Arch}}'        # arm64
bash docker/build.sh --platform linux/arm64 0.3.0-dev
# the same: bash docker/build.sh 0.3.0-dev (the default is the engine's own platform)
docker image ls bactolomew:0.3.0
docker image inspect --format '{{.Architecture}} {{json .Config.Labels}}' bactolomew:0.3.0
```

The build installs `docker/locks/linux-aarch64/` (the Dockerfile picks it from
`TARGETARCH`), skips fastGEAR, and must end with
`[build] done: bactolomew:0.3.0 (population methods: fastbaps)`. Expect about an
hour (the conda packages, the mob_suite database and the FastBAPS source build).

**Check the build and the tools**:

```bash
bash docker/update_locks.sh --check-image bactolomew:0.3.0    # 0 differences
I=bactolomew:0.3.0
for c in "bakta bakta --version" "bakta amrfinder --version" "diamond diamond version" \
         "mobsuite mob_recon --version" "phylogeny snippy --version" "phylogeny raxml-ng --version" \
         "gubbins run_gubbins.py --version" "gubbins iqtree2 --version" "panaroo panaroo --version" \
         "chewbbaca chewBBACA.py --version" "mlst mlst --version" "quast quast.py --version" \
         "checkm2 checkm2 --version" "mafft mafft --version"; do
    set -- $c; docker run --rm --network none "$I" conda run -n "$1" "${@:2}" 2>&1 | tail -1
done
docker run --rm --network none "$I" conda run -n fastbaps Rscript -e 'packageVersion("fastbaps")'
docker run --rm --network none "$I" conda run -n checkm2 python -c 'import tensorflow as tf; print(tf.__version__)'
docker run --rm --network none "$I" cat /opt/bactolomew/image_contents.json
docker run --rm --network none "$I" cat /opt/bactolomew/tool_versions.json > tool_versions_arm64.json
```

**Compare the outputs with the x86 image**: copy the smoke inputs (four Bakta GFF3s and
their contigs, the Snippy reference and the core alignment; about 45 MB, the folder that
`smoke_compare.sh` describes) from the Linux machine, run the arm64 image on them, and
compare with the x86 image's outputs from the same inputs:

```bash
# on the Mac
bash docker/smoke_compare.sh --run bactolomew:0.3.0 smoke_inputs smoke_arm64
# on the x86 machine, the same inputs and the x86 image of the same commit
bash docker/smoke_compare.sh --run bactolomew:0.3.0 smoke_inputs smoke_amd64
# either machine, with both output folders
bash docker/smoke_compare.sh --compare smoke_amd64 smoke_arm64
```

Every line should read `SAME`. A `DIFFERS` line is not necessarily a bug (floating point
or thread scheduling can change a tree's last digits or a tie between equally good
trees), but each needs a look: diff the two files and compare the trees' topologies.

**Report back**: the build log's last 50 lines and its duration, the image ID and size,
the `--check-image` result, the version lines above, `tool_versions_arm64.json`, the
`--compare` output, and for each `DIFFERS` the diff. If the build fails, the failing step
and its log; a FastBAPS download or source build error fails the `install_envs.sh fastbaps`
step, with the R error in the log.

## The mob_suite database

The mob_suite package's post-link script runs `mob_init` during `install_envs.sh mobsuite`.
`mob_init` downloads `data.tar.gz` from Zenodo record 10304948 with certificate checks
turned off and no checksum check (its `check_hash_or_size` is never called), then builds
the BLAST and mash databases and downloads the NCBI taxonomy through ete3. Conda offers no
way to skip a post-link script, so the image does not replace this step; it checks it:

- `docker/verify_image.sh` (a layer after the envs) runs `sha256sum --check` against
  `docker/locks/mobsuite_db.sha256`: the eight files of the Zenodo archive and the mash
  sketch built from them. The build fails on any difference, so a tampered or truncated
  download cannot reach a published image. These checksums are the same in 0.1.0
  (database of 2026-07-08) and 0.2.0 (2026-10-13): the Zenodo record is fixed.
- The NCBI taxonomy (`taxa.sqlite`) is whatever NCBI served on the build date and differs
  between builds (0.1.0 and 0.2.0 differ). It is not checked, only recorded: its sha256 and
  the download date go into `/opt/bactolomew/image_contents.json`, which every run copies
  into its provenance.

This is why the cached mob_suite layer matters in step 2. Moving the database out of the
image (route B, `mobsuite.database_directory`) would make it a versioned download instead.

## What a rebuild depends on

The published image is the archive; the Dockerfile is the recipe, and a rebuild from it
needs services outside this repository:

- **`pgcbioinfo/fastgear`** (MATLAB runtime and fastGEAR binaries) is a third-party image,
  last pushed in 2022. If it is deleted from Docker Hub, the Dockerfile cannot be built,
  even with the digest pinned. A local `docker save` copy keeps it buildable here, and a
  mirror in our own registry is the fix for everyone else (see "Mirroring fastGEAR" below).
  A build uses the mirror with
  `FASTGEAR_IMAGE=liviurotiul/fastgear:2022-12-01@sha256:<digest> bash docker/build.sh`.
- **conda.anaconda.org** must keep serving every URL in the locks.
- **Ubuntu's archive**: the apt layer installs what Ubuntu 24.04 serves on the build date
  (unpinned; 27 packages already differ between 0.1.0 and 0.2.0). Every installed package
  and version is recorded in `/opt/bactolomew/apt_packages.tsv`, and `verify_image.sh` fails
  the build when `libglu1-mesa` or `libasound2t64` is missing (the apt step tries fallback
  names). The apt layer is left unchanged for 0.2.0 so the rebuild in step 2 reuses the
  tested layer; for the next version, pin it with Ubuntu's snapshot service
  (`apt-get update --snapshot <timestamp>`) and install those two names without fallbacks.
- **Zenodo and NCBI** for the mob_suite database (above).

## Mirroring fastGEAR to our own registry (later)

**Not done yet** (it needs a registry login, so it waits for publishing). What is in place
without uploading anything:

- The app's default fastGEAR reference is pinned by digest (`ToolEnv.DEFAULT_FASTGEAR_IMAGE`,
  the same digest as the Dockerfile); `container.fastgear_image` names another one.
- A local copy outside the repository (in the maintainer's backups folder,
  `BACTolomew_backups/fastgear_image/fastgear_2022-12-01.tar.gz`: a `docker save` of image ID
  `sha256:ac1746a455a6...`, with its `.sha256` and a README). With
  `container.fastgear_archive` set to that file, a local build of the tools image whose
  fastGEAR image is missing loads it (`docker load` / `podman load`, checked against the
  `.sha256` file) instead of pulling it.
- Images built from the current Dockerfile carry the label `org.bactolomew.fastgear.image`,
  which the run provenance records as `environment.fastgear_image` (`null` for images
  built before the label: they used the Dockerfile's pinned digest).

The mirror, when there is a registry login:

```bash
F=pgcbioinfo/fastgear:latest@sha256:afb299e7fb327129ebf2acc1bcae00e0a54cd6fbe71d114d0a9f1642a2216b79
docker pull "$F"     # or: docker load -i fastgear_2022-12-01.tar.gz; F=pgcbioinfo/fastgear:latest
docker tag "$F" liviurotiul/fastgear:2022-12-01
docker push liviurotiul/fastgear:2022-12-01
docker image inspect --format '{{json .RepoDigests}}' liviurotiul/fastgear:2022-12-01
```

Then pin the mirror (`liviurotiul/fastgear:2022-12-01@sha256:<digest>`) as the
Dockerfile's `FASTGEAR_IMAGE` default and as `DEFAULT_FASTGEAR_IMAGE`, and make the package
public.

## Run provenance

Every run writes `environment` into `run_config.json` (at start) and `run_report.json`:
runtime, image reference (and the configured one when they differ), local image ID, repo
digests, the image's OCI version/revision/created/source labels, the pinned digest and
whether it matches, and the tool versions read from `/opt/bactolomew/tool_versions.json`.
Images without that file (0.1.0) record `tool_versions: null` and point to the probes in
`run_report.json` `tools.versions`. For Singularity it records the SIF path, its sha256
(computed once per file version and cached in `~/.cache/bactolomew/sif_sha256.json`) and
the labels `singularity inspect` reports (a SIF pulled from `docker://` keeps the OCI
labels).

`environment.databases` names the databases the tools used: the Bakta DB (`version.json`:
date, version, type, DOI, and the AMRFinderPlus DB release it bundles), the CheckM2 DIAMOND
file (path, size, modification time), and the mob_suite database (in the image: download
date and checksums from `image_contents.json`, or only the date for older images; or the
external `mobsuite.database_directory` and its date).

When a run folder is resumed with another image or other databases, the earlier
`environment` is kept in `environment_history` (with `recorded_at`, the earlier run's
start), because reused outputs (for example Bakta annotations) were computed by it. A folder
first run before provenance existed gets an entry with `environment: null`.

Cost: one `docker image inspect` and one short container start per run.

## The code DOI (GitHub to Zenodo)

Decided: the code DOI comes from Zenodo's GitHub integration, not from an upload by hand.
Not done yet (it needs the public repository and a Zenodo login):

0. Make the GitHub repository public (Zenodo archives only public repositories).
1. On zenodo.org, log in with GitHub, open Account > GitHub and switch the
   `liviurotiul/BACTolomew` repository on (before the release; Zenodo only archives releases
   made after the switch).
2. Optionally add a `CITATION.cff` (or `.zenodo.json`) at the repository root, so the record
   gets the right authors, title and license (GPL-3.0) instead of what Zenodo guesses.
3. Publish a GitHub release (a tag such as `v0.2.0`). Zenodo archives that tag's source and
   mints a version DOI plus a concept DOI that always resolves to the newest version.
4. Put the concept DOI into the README ("Citation"), and the version DOI into the release
   notes. Every later GitHub release becomes a new Zenodo version by itself.

The tools image is not part of that record: it is cited by its registry reference and digest
(step 6), which every run also records in its provenance.

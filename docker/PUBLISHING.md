# Publishing the tool image

This is the checklist for publishing `bactolomew-tools` to a registry, and (last section)
for minting the code DOI. None of it has been run yet: no image is pushed, no registry login
exists, and this machine's `bactolomew_config.json` still names the local `bactolomew-tools:0.1.0`.

Nothing else is published: there is no database bundle. First-time setup downloads the full
Bakta database (with its AMRFinderPlus database) and CheckM2 from their own sources, and
makes `DB/` with empty catalogs and the shipped cgMLST scheme species catalog.
(`tools/zenodo_upload.py` was written for a database bundle and is not needed for this.)

## Version and tag

Decided 2026-10-21: the first public release is 0.2.0. `__version__` in
`src/pipeline/__init__.py`, `TOOLS_IMAGE` in `src/pipeline/ToolEnv.py`, `docker/build.sh`
and `src/scripts/demo_install.py` all name 0.2.0; tag the GitHub release `v0.2.0`.
Two local images exist: `bactolomew-tools:0.1.0` (29.3 GB, the image the reference run used)
and `bactolomew-tools:0.2.0` (22.2 GB, the tested rebuild from the locks). Until a published
image is configured, setup builds the image named by `TOOLS_IMAGE` locally from the current
Dockerfile (`container.build_on_missing` defaults to true).

## What is already in place

| Piece | Where |
| --- | --- |
| Exact package locks, one per env (`conda list --explicit --md5`, exported from 0.1.0) | `docker/locks/<env>.lock` |
| Readable tool version summary (copied into the image as `/opt/bactolomew/tool_versions.json`) | `docker/locks/tool_versions.json` |
| Dockerfile installing from the locks, base images pinned by digest, OCI labels | `docker/Dockerfile` |
| Env install and trim step (one layer per env group) | `docker/install_envs.sh` |
| Build-time checks: MATLAB runtime system libraries present, mob_suite database files equal to recorded sha256s; writes `/opt/bactolomew/image_contents.json` and `apt_packages.tsv` | `docker/verify_image.sh`, `docker/locks/mobsuite_db.sha256` |
| Build script (version argument, labels, `linux/amd64`, no credentials in labels, fastgear mirror override, prints the image ID) | `docker/build.sh` |
| Solving envs from their yamls (new locks after a tool change) | `docker/Dockerfile.solve` |
| Lock regeneration from a built image | `docker/update_locks.sh`, `docker/make_tool_versions.py` |
| Tool runs on a tiny input set in two images, outputs compared | `docker/smoke_compare.sh` |
| Build context limited to the files the Dockerfile copies (about 260 KB) | `.dockerignore` |
| Digest pinning in `container.image` (`name:tag@sha256:...`), with a warning on mismatch | `src/pipeline/ContainerRuntime.py` |
| Image acquisition: pull any reference that names a registry (Docker and Singularity), build locally only when `container.build_on_missing` is true | `src/pipeline/ContainerRuntime.py` |
| Run provenance: `environment` (image, tool versions, databases) and `environment_history` in `run_config.json` and `run_report.json` | `src/pipeline/Runner.py` |

Why a new version instead of pushing 0.1.0: the local 0.1.0 was built from an older
Dockerfile (label `contexar-tools`, yamls under `/opt/contexar/yamls`, no `libncurses.so.5`
link), so the current Dockerfile does not describe it, and it has no
`/opt/bactolomew/tool_versions.json`. Version 0.2.0 is built from the locks, so it carries
the same tool versions as 0.1.0 (see `tool_versions.json`), with a matching Dockerfile.

## Size

0.1.0 is 29.3 GB. The Dockerfile applies only the reductions that cannot change what a
tool does:

- envs grouped per layer so shared packages are stored once (about -5 GB):
  `assembly_qc bakta gubbins phylogeny mafft` (about 6.8 GB), `panaroo chewbbaca mlst diamond`
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
- compilers and sysroot that R and k8 pull into gubbins, assembly_qc and phylogeny
  (about 1 GB): needs a gubbins and quast/minimap2 smoke test.
- `share/locale`, `include/`, test folders inside site-packages. Never delete `*.pyc`.

## Steps to publish (in order)

### 1. Choose the registry

- **GHCR** (`ghcr.io/liviurotiul/bactolomew-tools`): free for public images, tied to the
  GitHub repo (the `org.opencontainers.image.source` label links the package to it). Check
  the current per-layer size limit before pushing (reported as about 10 GB; the largest
  layer here is about 7 GB uncompressed).
- Docker Hub or Quay.io work the same way; Docker Hub has pull rate limits for anonymous
  users.
- The code DOI comes from the GitHub to Zenodo link (last section). Archiving the image
  itself on Zenodo (`docker save | gzip`, about the compressed image size) is optional.

The steps below assume GHCR; `IMAGE=ghcr.io/liviurotiul/bactolomew-tools`.

### 2. Commit and push the source, then build

The local test image `bactolomew-tools:0.2.0` was built from an uncommitted tree: its
`org.opencontainers.image.revision` label is `d11b88d...-dirty`, which names no commit. Do
not publish it. Instead:

1. Commit `docker/`, `yamls/` and `.dockerignore` (build.sh appends `-dirty`, and warns,
   when those paths have local changes).
2. Push that commit to `origin`, so the `revision` label resolves on GitHub.
3. Rebuild with the same build.sh, **without pruning the Docker build cache first**
   (no `docker builder prune`, no `docker system prune`). BuildKit then reuses the cached
   apt and env layers of the tested 0.2.0, so the tools and the mob_suite database are
   byte-identical to what was smoke-tested. A pruned cache fetches the apt packages and the
   mob_suite database again, with different content (see below).

```bash
bash docker/build.sh 0.2.0 ghcr.io/liviurotiul/bactolomew-tools
docker image ls ghcr.io/liviurotiul/bactolomew-tools
docker image inspect --format '{{json .Config.Labels}}' ghcr.io/liviurotiul/bactolomew-tools:0.2.0
```

The Dockerfile changed after the tested 0.2.0 was built, only after the env layers, so the
rebuild redoes only the layers from the `tool_versions.json` copy on (that file gained a
note; the fastGEAR copy that follows is the same content): the image-wide
`LD_LIBRARY_PATH` is gone (the MATLAB runtime's 2016 libstdc++, libgcc_s, libssl and others
no longer shadow the system and conda libraries; fastGEAR's launcher sets its own path,
and fastGEAR on its example data gave identical output with and without the variable), and
`docker/verify_image.sh` runs as a new layer. Check the build log shows `CACHED` for the
`install_envs.sh` steps; if it does not, the env layers were rebuilt and step 4 must be
repeated in full.

The build needs network access to conda-forge, bioconda and Docker Hub (and Zenodo and
NCBI when the mob_suite layer is rebuilt). It always builds `linux/amd64` (the locks and
the MATLAB runtime are x86_64 only), also on an arm64 host. It does not touch
`bactolomew-tools:0.1.0`.

### 3. Check the build matches the locks

```bash
bash docker/update_locks.sh ghcr.io/liviurotiul/bactolomew-tools:0.2.0
git diff --stat -- docker/locks
```

The `.lock` files must be unchanged; only `locks_exported_from` in `tool_versions.json`
changes. Keep that change (it now names the published image) or discard it with
`git checkout -- docker/locks/tool_versions.json`.

This check reads `conda-meta`, so it holds by construction (the envs were installed from
those locks) and cannot see the files `install_envs.sh` deleted after install (static
libraries in `lib/`, docs). It shows the right packages went in, not that the tools still
work: step 4 is the evidence for that.

### 4. Smoke test

```bash
I=ghcr.io/liviurotiul/bactolomew-tools:0.2.0
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
# A classic personal access token with write:packages (and read:packages)
echo "$GHCR_TOKEN" | docker login ghcr.io -u liviurotiul --password-stdin
docker push ghcr.io/liviurotiul/bactolomew-tools:0.2.0
```

On GitHub, Packages > bactolomew-tools > Package settings: set the visibility to public.

### 6. Record the digest

```bash
docker image inspect --format '{{index .RepoDigests 0}}' ghcr.io/liviurotiul/bactolomew-tools:0.2.0
# ghcr.io/liviurotiul/bactolomew-tools@sha256:<digest>
```

Write the digest into `docker/locks/tool_versions.json` (`locks_exported_from`), the README
and the release notes. Optionally tag the commit (`git tag tools-0.2.0`).

### 7. Check a clean pull

On another machine, or after `docker image rm` of the pushed tag only:

```bash
docker pull ghcr.io/liviurotiul/bactolomew-tools:0.2.0@sha256:<digest>
```

### 8. Singularity / Apptainer (HPC)

With `runtime: "singularity"`, `container.image` set to the full reference (step 9) and
`container.sif_path` set to where the `.sif` should go, the app pulls it on first start
(`singularity pull <sif_path> docker://<image>`); no Docker daemon is needed. By hand:

```bash
singularity pull bactolomew-tools_0.2.0.sif \
    docker://ghcr.io/liviurotiul/bactolomew-tools:0.2.0@sha256:<digest>
# or, from a local Docker image:
singularity build bactolomew-tools_0.2.0.sif docker-daemon://ghcr.io/liviurotiul/bactolomew-tools:0.2.0
```

The `.sif` is about the compressed size of the image; build it on a filesystem with room
for twice that (`SINGULARITY_TMPDIR`). Set `container.sif_path` to it. The Singularity
path is untested here (no Singularity on this machine): try it once on the HPC before
announcing it.

### 9. Point the app at the published image

`bactolomew_config.json` (per machine, not in git):

```json
"container": {
  "image": "ghcr.io/liviurotiul/bactolomew-tools:0.2.0@sha256:<digest>",
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

1. Edit `yamls/<env>.yaml` (the intent).
2. Solve the changed envs from their yamls:
   `docker build --platform linux/amd64 -f docker/Dockerfile.solve --build-arg ENVS="<env> ..." -t bactolomew-solve:tmp .`
3. `bash docker/update_locks.sh bactolomew-solve:tmp` and review `git diff -- docker/locks`
   (only the solved envs and `base` get new locks; `tool_versions.json` keeps the others).
4. Build the real image with `docker/build.sh <new version>`, run step 4 against the
   previous version, and publish as above.

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
  `FASTGEAR_IMAGE=ghcr.io/liviurotiul/fastgear:2022-12-01@sha256:<digest> bash docker/build.sh`.
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
docker tag "$F" ghcr.io/liviurotiul/fastgear:2022-12-01
docker push ghcr.io/liviurotiul/fastgear:2022-12-01
docker image inspect --format '{{json .RepoDigests}}' ghcr.io/liviurotiul/fastgear:2022-12-01
```

Then pin the mirror (`ghcr.io/liviurotiul/fastgear:2022-12-01@sha256:<digest>`) as the
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

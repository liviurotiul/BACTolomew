# BACTolomew unified tool image: every bioinformatics tool the pipeline runs, each in
# its own conda env (env name == ToolEnv env key), plus fastGEAR (MATLAB MCR) folded in.
#
# Build from the repo root:   bash docker/build.sh [version]     (default version 0.2.0)
#   or directly:              docker build --platform linux/amd64 -f docker/Dockerfile -t bactolomew-tools:0.2.0 .
# The image is linux/amd64 only (the locks are linux-64 packages, the MATLAB runtime is
# x86_64): build.sh passes --platform linux/amd64, so an arm64 host builds amd64 too.
#
# Run a tool:                 docker run --rm bactolomew-tools:0.2.0 conda run -n bakta bakta --version
#
# Reproducibility: every env is installed from docker/locks/<env>.lock, an explicit
# `conda list --explicit --md5` export (exact package URLs and md5s, no solver). The
# yamls/*.yaml files stay the human-facing intent (which tools an env should hold); after
# changing a yaml, solve it once with docker/Dockerfile.solve and export the new lock with
# docker/update_locks.sh (see docker/PUBLISHING.md). Base images are pinned by digest; the
# fastgear image is third party and unmaintained, so PUBLISHING.md says how to mirror it.

# --- stage: fastgear (MATLAB MCR runtime + fastGEAR launcher/specs) ---
# Digest of the :latest tag the 0.1.0 image was built from (unchanged on Docker Hub since
# 2022-12-01). If it disappears from Docker Hub, build from a mirror of the same image
# (PUBLISHING.md): FASTGEAR_IMAGE=<mirror>@sha256:<digest> bash docker/build.sh.
ARG FASTGEAR_IMAGE=pgcbioinfo/fastgear:latest@sha256:afb299e7fb327129ebf2acc1bcae00e0a54cd6fbe71d114d0a9f1642a2216b79
FROM ${FASTGEAR_IMAGE} AS fastgear

# --- stage: main ---
# Miniforge3 26.3.2-3 (conda 26.3.2, mamba 2.5.0, Ubuntu 24.04), the release 0.1.0 was built on.
FROM condaforge/miniforge3:26.3.2-3@sha256:532f6ee7a858b009dc895f8313eb6ed875f05a455ec57d832aa6b4c66e2799b9
SHELL ["/bin/bash", "-lc"]

# System libraries the MATLAB MCR runtime needs (fastGEAR); harmless for the conda tools.
# Core X libs are stable names; the audio/GLU libs got t64 renames on newer Debian, so
# the fallbacks below try both names. docker/verify_image.sh (a later layer) fails the
# build when the Ubuntu 24.04 names (libglu1-mesa, libasound2t64) are missing, and records
# every installed package in /opt/bactolomew/apt_packages.tsv. The packages are whatever
# Ubuntu serves on the build date; this layer is kept as is so a rebuild reuses the cached
# layer of the tested image (see docker/PUBLISHING.md, "Rebuilding").
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates procps \
        libx11-6 libxext6 libxt6 libxmu6 libxpm4 libxrandr2 libxcomposite1 libxcursor1 \
        libxdamage1 libxfixes3 libxft2 libxi6 libxrender1 libxtst6 libgl1 \
    && (apt-get install -y --no-install-recommends libglu1-mesa || apt-get install -y --no-install-recommends libglu1 || true) \
    && (apt-get install -y --no-install-recommends libasound2t64 || apt-get install -y --no-install-recommends libasound2 || true) \
    && rm -rf /var/lib/apt/lists/*

# The locks (what is installed) travel with the image, at /opt/bactolomew/locks.
COPY docker/locks/*.lock /opt/bactolomew/locks/
COPY docker/install_envs.sh /opt/bactolomew/install_envs.sh

# One conda env per lock, env name == lock basename == ToolEnv env key. The envs are
# grouped so the ones sharing the most packages are created in the same layer (shared
# files are stored once) while every layer stays well under the per-layer size limits of
# common registries (largest group about 7 GB uncompressed). mobsuite is alone: its
# package downloads the mob_suite database (about 2.9 GB) at install time; the files
# are checked against docker/locks/mobsuite_db.sha256 by docker/verify_image.sh below.
RUN /opt/bactolomew/install_envs.sh assembly_qc bakta gubbins phylogeny mafft
RUN /opt/bactolomew/install_envs.sh panaroo chewbbaca mlst diamond
RUN /opt/bactolomew/install_envs.sh mobsuite

# The yamls (intent) and the tool version summary, copied after the envs so editing them
# does not rebuild an env; the run provenance reads /opt/bactolomew/tool_versions.json.
COPY yamls/ /opt/bactolomew/yamls/
COPY docker/locks/tool_versions.json /opt/bactolomew/tool_versions.json

# fastGEAR: MATLAB MCR + fastGEAR/plot binaries + launcher + specs, kept at /opt/* so the
# launcher's hardcoded /opt/v901 and /opt/fastGEAR references still resolve. (The fastgear
# image has no /opt/conda or /opt/bactolomew, so merging into our /opt is safe.)
COPY --from=fastgear /opt /opt
# The R2016a runtime links libncurses.so.5, which Ubuntu 24.04 no longer ships; ncurses 6
# serves it (the pipeline links the same file in fastGEAR's working directory for older images).
RUN ln -s /usr/lib/x86_64-linux-gnu/libncursesw.so.6 /opt/v901/bin/glnxa64/libncurses.so.5
# No image-wide LD_LIBRARY_PATH: /opt/v901 ships old libstdc++, libgcc_s, libssl and
# others that would shadow the system and conda libraries of every other tool. The
# launcher run_fastGEAR.sh sets the MATLAB runtime path for fastGEAR alone.
ENV MCR_ROOT=/opt/v901

# Build-time checks and records (system libraries, mob_suite database checksums), in a
# layer of their own after the envs so a new check never rebuilds an env.
COPY docker/locks/mobsuite_db.sha256 /opt/bactolomew/locks/mobsuite_db.sha256
COPY docker/verify_image.sh /opt/bactolomew/verify_image.sh
RUN /opt/bactolomew/verify_image.sh

# OCI labels. docker/build.sh passes the version, the git commit, the build time and the
# repository URL; they are last so a new commit or date only rebuilds this metadata.
ARG VERSION=0.2.0
ARG REVISION=unknown
ARG CREATED=unknown
ARG SOURCE=https://github.com/liviurotiul/BACTolomew
# The fastGEAR image the fastgear stage came from, read into the run provenance.
# FASTGEAR_PROVENANCE: the pinned ref a local copy was checked against (a `docker load`
# copy has no registry digest, so FASTGEAR_IMAGE is then its bare name).
ARG FASTGEAR_IMAGE
ARG FASTGEAR_PROVENANCE=${FASTGEAR_IMAGE}
LABEL org.opencontainers.image.title="bactolomew-tools" \
      org.opencontainers.image.description="BACTolomew bioinformatics tool suite (bakta, amrfinderplus, mlst, mob_suite, gubbins, panaroo, chewbbaca, diamond, quast, checkm2, snippy, raxml-ng, iqtree, mafft, fastGEAR), one conda env per tool, installed from explicit locks" \
      org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.revision="${REVISION}" \
      org.opencontainers.image.created="${CREATED}" \
      org.opencontainers.image.source="${SOURCE}" \
      org.opencontainers.image.url="${SOURCE}" \
      org.opencontainers.image.base.name="docker.io/condaforge/miniforge3:26.3.2-3" \
      org.opencontainers.image.base.digest="sha256:532f6ee7a858b009dc895f8313eb6ed875f05a455ec57d832aa6b4c66e2799b9" \
      org.bactolomew.fastgear.image="${FASTGEAR_PROVENANCE}"

# No ENTRYPOINT: the host invokes `... <image> conda run -n <env> <exe> <args>`.
CMD ["bash"]

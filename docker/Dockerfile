# BACTolomew unified tool image: every bioinformatics tool the pipeline runs, each in
# its own conda env (env name == ToolEnv env key), plus fastGEAR (MATLAB MCR) folded in.
#
# Build from the repo root:   docker build -f docker/Dockerfile -t bactolomew-tools:0.1.0 .
# Run a tool:                 docker run --rm bactolomew-tools:0.1.0 conda run -n bakta bakta --version

# --- stage: fastgear (MATLAB MCR runtime + fastGEAR launcher/specs) ---
FROM pgcbioinfo/fastgear:latest AS fastgear

# --- stage: main ---
FROM condaforge/miniforge3:latest
SHELL ["/bin/bash", "-lc"]

# System libraries the MATLAB MCR runtime needs (fastGEAR); harmless for the conda tools.
# Core X libs are stable names; the audio/GLU libs got t64 renames on newer Debian, so
# install those best-effort (fastGEAR compute doesn't need them).
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates procps \
        libx11-6 libxext6 libxt6 libxmu6 libxpm4 libxrandr2 libxcomposite1 libxcursor1 \
        libxdamage1 libxfixes3 libxft2 libxi6 libxrender1 libxtst6 libgl1 \
    && (apt-get install -y --no-install-recommends libglu1-mesa || apt-get install -y --no-install-recommends libglu1 || true) \
    && (apt-get install -y --no-install-recommends libasound2t64 || apt-get install -y --no-install-recommends libasound2 || true) \
    && rm -rf /var/lib/apt/lists/*

# Speed up/soothe conda solves.
RUN conda config --system --set channel_priority flexible

# One conda env per yaml. Separate RUN layers so a single yaml change doesn't rebuild all
# envs, and so a failing env is easy to spot. Env name == basename == ToolEnv env key.
COPY yamls/ /opt/bactolomew/yamls/
RUN mamba env create -n bakta       -f /opt/bactolomew/yamls/bakta.yaml       && mamba clean -a -y
RUN mamba env create -n diamond     -f /opt/bactolomew/yamls/diamond.yaml     && mamba clean -a -y
RUN mamba env create -n mobsuite    -f /opt/bactolomew/yamls/mobsuite.yaml    && mamba clean -a -y
RUN mamba env create -n phylogeny   -f /opt/bactolomew/yamls/phylogeny.yaml   && mamba clean -a -y
RUN mamba env create -n gubbins     -f /opt/bactolomew/yamls/gubbins.yaml     && mamba clean -a -y
RUN mamba env create -n panaroo     -f /opt/bactolomew/yamls/panaroo.yaml     && mamba clean -a -y
RUN mamba env create -n chewbbaca   -f /opt/bactolomew/yamls/chewbbaca.yaml   && mamba clean -a -y
RUN mamba env create -n mlst        -f /opt/bactolomew/yamls/mlst.yaml        && mamba clean -a -y
RUN mamba env create -n assembly_qc -f /opt/bactolomew/yamls/assembly_qc.yaml && mamba clean -a -y
RUN mamba env create -n mafft       -f /opt/bactolomew/yamls/mafft.yaml       && mamba clean -a -y

# fastGEAR: MATLAB MCR + fastGEAR/plot binaries + launcher + specs, kept at /opt/* so the
# launcher's hardcoded /opt/v901 and /opt/fastGEAR references still resolve. (The fastgear
# image has no /opt/conda, so merging into our /opt is safe.)
COPY --from=fastgear /opt /opt
ENV MCR_ROOT=/opt/v901 \
    LD_LIBRARY_PATH=/opt/v901/runtime/glnxa64:/opt/v901/bin/glnxa64:/opt/v901/sys/os/glnxa64:/opt/v901/sys/opengl/lib/glnxa64

LABEL org.opencontainers.image.title="bactolomew-tools" \
      org.opencontainers.image.description="BACTolomew bioinformatics tool suite (bakta, mlst, mobsuite, gubbins, panaroo, chewbbaca, diamond, quast, checkm2, snippy, raxml-ng, mafft, fastGEAR)"

# No ENTRYPOINT: the host invokes `... <image> conda run -n <env> <exe> <args>`.
CMD ["bash"]

# Minimal Augustus container built from Bioconda using micromamba
# Build locally, then push to your registry (Quay/GHCR), or convert to SIF.

FROM mambaorg/micromamba:1.5.8
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Install augustus and common helpers
# - Use conda-forge + bioconda, pin versions as needed
# - Clean caches to keep image small
RUN micromamba install -y -n base -c bioconda -c conda-forge \
    augustus=3.* \
    gffread \
    bedtools \
    samtools \
    && micromamba clean -a -y

# Augustus expects AUGUSTUS_CONFIG_PATH; in Bioconda this lives under $CONDA_PREFIX/share/augustus/config
ENV AUGUSTUS_CONFIG_PATH=/opt/conda/envs/base/share/augustus/config
ENV PATH=/opt/conda/envs/base/bin:/opt/conda/bin:$PATH
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV TZ=UTC

# Default to bash
SHELL ["/bin/bash", "-lc"]
ENTRYPOINT ["/bin/bash"]

RUN which augustus || (echo "WARN: augustus not found in PATH at build time; it may exist under /opt/conda/envs/base/bin at runtime" && ls -R /opt/conda/envs/base/share 2>/dev/null || true)

# Note: Config location may vary across builds; AUGUSTUS_CONFIG_PATH is set and
# runtime scripts detect/override if needed.

# Minimal Augustus container built from Bioconda using micromamba
# Build locally, then push to your registry (Quay/GHCR), or convert to SIF.

FROM mambaorg/micromamba:1.5.8
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Install augustus and common helpers
# - Use conda-forge + bioconda, pin versions as needed
# - Clean caches to keep image small
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    augustus=3.* \
    gffread \
    bedtools \
    samtools \
    && micromamba clean -a -y

# Augustus expects AUGUSTUS_CONFIG_PATH; in Bioconda this lives under $CONDA_PREFIX/share/augustus/config
ENV AUGUSTUS_CONFIG_PATH=/opt/conda/share/augustus/config
ENV PATH=/opt/conda/bin:$PATH
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV TZ=UTC

# Default to bash
SHELL ["/bin/bash", "-lc"]
ENTRYPOINT ["/bin/bash"]

# Build-time sanity checks (fail early if configs or binary missing)
RUN echo "AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH" && \
    ls -d "$AUGUSTUS_CONFIG_PATH" && \
    augustus --version && \
    augustus --printSpecies | head -n 20

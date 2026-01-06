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
ENV AUGUSTUS_CONFIG_PATH=/opt/conda/share/augustus/config
ENV PATH=/opt/conda/bin:$PATH
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV TZ=UTC

# Default to bash
SHELL ["/bin/bash", "-lc"]
ENTRYPOINT ["/bin/bash"]

RUN if ! command -v augustus >/dev/null 2>&1; then \
            echo "ERROR: augustus not found in PATH"; \
            echo "Searching under /opt/conda for the binary..."; \
            find /opt/conda -maxdepth 4 -type f -name augustus -print; \
            exit 1; \
        fi

# Note: Config location may vary across builds; AUGUSTUS_CONFIG_PATH is set and
# runtime scripts detect/override if needed.

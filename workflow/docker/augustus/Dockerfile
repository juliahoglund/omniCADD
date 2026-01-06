# Minimal Augustus container built from Bioconda using micromamba
# Build locally, then push to your registry (Quay/GHCR), or convert to SIF.

FROM ghcr.io/conda-incubator/miniforge:latest

# Install augustus and common helpers
# - Use conda-forge + bioconda, pin versions as needed
# - Clean caches to keep image small
RUN conda create -y -n augustus -c bioconda -c conda-forge \
    augustus=3.* \
    gffread \
    bedtools \
    samtools \
    && conda clean -a -y

# Augustus expects AUGUSTUS_CONFIG_PATH; in Bioconda this lives under $CONDA_PREFIX/share/augustus/config
ENV PATH=/opt/conda/envs/augustus/bin:/opt/conda/bin:$PATH
ENV AUGUSTUS_CONFIG_PATH=/opt/conda/envs/augustus/share/augustus/config
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV TZ=UTC

# Default to bash
SHELL ["/bin/bash", "-lc"]
ENTRYPOINT ["/bin/bash"]

RUN which augustus && echo "Configs at: $AUGUSTUS_CONFIG_PATH" && ls -d "$AUGUSTUS_CONFIG_PATH" || (echo "ERROR: augustus not found or configs missing" && exit 1)

# Note: Config location may vary across builds; AUGUSTUS_CONFIG_PATH is set and
# runtime scripts detect/override if needed.

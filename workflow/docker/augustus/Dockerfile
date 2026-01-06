# Minimal Augustus container built from Bioconda using micromamba
# Build locally, then push to your registry (Quay/GHCR), or convert to SIF.

FROM mambaorg/micromamba:1.5.8

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

# Default to bash
SHELL ["/bin/bash", "-lc"]
ENTRYPOINT ["/bin/bash"]

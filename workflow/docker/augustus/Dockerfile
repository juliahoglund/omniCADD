# Minimal Augustus container built from Bioconda using Miniforge (Docker Hub)
FROM docker.io/condaforge/miniforge3:23.11.0-0

# Install augustus and common helpers
# - Use conda-forge + bioconda, pin versions as needed
# - Clean caches to keep image small
RUN conda create -y -n augustus -c conda-forge -c bioconda \
            augustus=3.* \
            snpeff=5.* \
            openjdk=17 \
            gffread \
            bedtools \
            samtools \
        && conda clean -a -y

# Normalize AUGUSTUS config location and PATH
ENV PATH=/opt/conda/envs/augustus/bin:/opt/conda/bin:$PATH
SHELL ["/bin/bash", "-lc"]

# Ensure AUGUSTUS config path is present at share/augustus/config and export AUGUSTUS_CONFIG_PATH
RUN set -euxo pipefail; \
    PREFIX="/opt/conda/envs/augustus"; \
    if [ -d "$PREFIX/config" ]; then CFG="$PREFIX/config"; \
    elif [ -d "$PREFIX/share/augustus/config" ]; then CFG="$PREFIX/share/augustus/config"; \
    else echo "ERROR: Augustus config directory not found in $PREFIX" >&2; exit 1; fi; \
    mkdir -p "$PREFIX/share/augustus"; \
    if [ ! -e "$PREFIX/share/augustus/config" ]; then ln -s "$CFG" "$PREFIX/share/augustus/config"; fi; \
    echo "Detected Augustus config at: $CFG"; \
    true

ENV AUGUSTUS_CONFIG_PATH=/opt/conda/envs/augustus/share/augustus/config
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV TZ=UTC

ENTRYPOINT ["/bin/bash"]

RUN which augustus && echo "Configs at: $AUGUSTUS_CONFIG_PATH" && ls -d "$AUGUSTUS_CONFIG_PATH" || (echo "ERROR: augustus not found or configs missing" && exit 1)
RUN which snpEff && snpEff -version || (echo "ERROR: snpEff not found" && exit 1)

# Note: Config location may vary across builds; AUGUSTUS_CONFIG_PATH is set and
# runtime scripts detect/override if needed.

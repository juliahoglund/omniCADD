# Minimal SNPEff container built from Bioconda using Miniforge
FROM docker.io/condaforge/miniforge3:23.11.0-0

# Install snpeff and openjdk
# - Use conda-forge + bioconda
# - Clean caches to keep image small
RUN conda create -y -n snpeff -c conda-forge -c bioconda \
            snpeff=5.1 \
            openjdk=21 \
        && conda clean -a -y

# Add conda environment to path
ENV PATH=/opt/conda/envs/snpeff/bin:/opt/conda/bin:$PATH

ENTRYPOINT ["/bin/bash"]

# Verify installation
RUN java -version && which snpEff && snpEff -version || (echo "ERROR: snpEff not found or incompatible Java" && exit 1)

# Use the official Ensembl VEP image as base
# PINNED to release 111.0 for reproducibility
# To use an absolute freeze, replace the tag with the SHA256 digest:
# FROM ensemblorg/ensembl-vep@sha256:4b4f...
FROM ensemblorg/ensembl-vep:release_111.0

################## METADATA ######################
LABEL maintainer="Julia Höglund"
LABEL version=111.0
LABEL description="Custom VEP 111.0 image for omniCADD with fixed dependencies"

################## DEPENDENCIES ##################
# Switch to root to install system-level Perl dependencies
USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    liblist-moreutils-perl \
    libdbd-mysql-perl \
    libarchive-extract-perl \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Switch back to the 'vep' user provided by the base image
USER vep

# Ensure the entrypoint is correct
WORKDIR /opt/vep/src/ensembl-vep
ENTRYPOINT ["vep"]

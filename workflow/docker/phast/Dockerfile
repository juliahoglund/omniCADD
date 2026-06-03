FROM ubuntu:jammy

### ENVIRONMENT ###
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH=/opt/phast/bin:/opt/kent/bin:$PATH
ENV DEBIAN_FRONTEND=noninteractive

### CORE DEPENDENCIES ###
RUN apt-get update && \
    apt-get -y upgrade && \
    apt-get -y install \
      build-essential \
      wget \
      curl \
      git \
      unzip \
      make \
      gcc \
      g++ \
      zlib1g-dev \
      libbz2-dev \
      liblzma-dev \
      libcurl4-openssl-dev \
      libssl-dev \
      libncurses5-dev \
      liblapack-dev \
      pkg-config \
      rsync \
      && apt-get clean && \
      rm -rf /var/lib/apt/lists/*

### INSTALL PHAST ###
WORKDIR /opt
RUN wget https://github.com/CshlSiepelLab/phast/archive/refs/tags/v1.5.tar.gz && \
    tar -xzf v1.5.tar.gz && \
    rm v1.5.tar.gz && \
    mv phast-1.5 phast && \
    cd phast/src && \
    make

# Test phast installation
RUN phyloFit 2>&1 | head -n 1 && \
    phastCons 2>&1 | head -n 1 && \
    phyloP 2>&1 | head -n 1

### INSTALL UCSC KENT TOOLS (for wig2bed) ###
WORKDIR /opt
RUN mkdir -p kent/bin && \
    cd kent/bin && \
    rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/wigToBigWig . && \
    chmod +x wigToBigWig

# Download bedops for wig2bed (alternative since UCSC doesn't have wig2bed directly)
RUN wget https://github.com/bedops/bedops/releases/download/v2.4.41/bedops_linux_x86_64-v2.4.41.tar.bz2 && \
    tar -xjf bedops_linux_x86_64-v2.4.41.tar.bz2 && \
    mv bin/* /opt/kent/bin/ && \
    rm -rf bedops_linux_x86_64-v2.4.41.tar.bz2 bin

# Test wig2bed
RUN wig2bed --help 2>&1 | head -n 1

WORKDIR /data

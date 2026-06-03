FROM ubuntu:jammy

### ENVIRONMENT ###
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH=/opt/gerp/bin:$PATH
ENV DEBIAN_FRONTEND=noninteractive

### CORE DEPENDENCIES ###
RUN apt-get update && \
    apt-get -y upgrade && \
    apt-get -y install \
      build-essential \
      wget \
      make \
      gcc \
      g++ \
      && apt-get clean && \
      rm -rf /var/lib/apt/lists/*

### INSTALL GERP++ ###
WORKDIR /opt
RUN wget http://mendel.stanford.edu/SidowLab/downloads/gerp/gerp++.tar.gz && \
    tar -xzf gerp++.tar.gz && \
    rm gerp++.tar.gz && \
    mv gerp++ gerp && \
    cd gerp/src && \
    make

# Test GERP installation
RUN gerpcol 2>&1 | head -n 1

WORKDIR /data


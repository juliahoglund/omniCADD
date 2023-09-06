FROM ubuntu:jammy

### ENVIRONMENT ###
ENV  LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV  PATH=/opt/src:$PATH
ENV  LC_ALL=C
ENV  DEBIAN_FRONTEND=noninteractive

### CORE DEPENDENCIES ###
WORKDIR /opt/src/

RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get clean

RUN apt-get -y install \
      build-essential \
      software-properties-common \
      libatlas-base-dev \
      wget \
      curl \
      git \
      unzip \
      libbz2-dev \
      libc6-dev \
      libcurl4-openssl-dev \
      libdbd-mysql-perl \
      libffi-dev \
      libgdbm-dev \
      libgsl0-dev \
      liblzma-dev \
      libmysqlclient-dev \
      libncurses5-dev \
      libncursesw5-dev \
      libnss3-dev \
      libpcap-dev \
      libperl-dev \
      libreadline-dev \
      libsqlite3-dev \
      libtk8.6 \
      libz-dev \
      openssl \
      zlib1g-dev \
      tk-dev \
      perl \
      cpanminus  \
      mysql-server \
      python3 \
      pip \
        && rm -rf /var/lib/apt/lists/*

### PERL & ENSEMBL API ###
# set perl environment variables
ENV PERL_PATH=/opt/src/perl/
ENV PERL5LIB=$PERL_PATH:$PERL_PATH/lib/perl5:$PERL5LIB
ENV PERL_MM_OPT="INSTALL_BASE=$PERL_PATH"
ENV PERL_MB_OPT="--install_base $PERL_PATH"
ENV PATH="$PERL_PATH/bin:$PATH"

### BioPerl
WORKDIR /opt/src/
RUN wget https://github.com/bioperl/bioperl-live/archive/refs/tags/release-1-7-2.zip
RUN unzip release-1-7-2.zip
RUN mv bioperl-live-release-1-7-2 bioperl-1.7.2

### Ensembl API
RUN git config --global http.sslverify false
RUN git clone https://github.com/Ensembl/ensembl.git
RUN git clone https://github.com/Ensembl/ensembl-variation.git
RUN git clone https://github.com/Ensembl/ensembl-compara.git
RUN git clone https://github.com/Ensembl/ensembl-funcgen.git
RUN git clone https://github.com/Ensembl/ensembl-tools.git

# RUN wget ftp://ftp.ensembl.org/pub/ensembl-api.tar.gz
# RUN tar -zxf ensembl-api.tar.gz

ENV PERL5LIB ${PERL5LIB}:/opt/src/bioperl-1.7.2
ENV PERL5LIB=${PERL5LIB}:/opt/src/ensembl/modules
ENV PERL5LIB=${PERL5LIB}:/opt/src/ensembl-variation/modules
ENV PERL5LIB=${PERL5LIB}:/opt/src/ensembl-compara/modules
ENV PERL5LIB=${PERL5LIB}:/opt/src/ensembl-funcgen/modules
ENV PERL5LIB=${PERL5LIB}:/opt/src/ensembl-tools/modules

# Install perl dependencies
RUN cpanm --force -l $PERL_PATH \
    Archive::Extract \
    Archive::Tar \
    Archive::Zip \
    DBI \
    DBD \
    DBD::mysql

# install python dependencies
RUN pip3 install --upgrade pip
RUN pip3 install pysam

### ping to test ###
WORKDIR /opt/src/ensembl/misc-scripts
RUN ./ping_ensembl.pl

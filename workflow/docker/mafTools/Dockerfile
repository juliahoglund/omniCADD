FROM python:2.7.18-stretch

ENV  LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV  PATH=/opt/mafTools/bin:$PATH
ENV  DEBIAN_FRONTEND=noninteractive
ENV  LC_ALL=C

# Debian Stretch is EOL - use archive repositories
RUN sed -i 's|http://deb.debian.org|http://archive.debian.org|g' /etc/apt/sources.list && \
    sed -i 's|http://security.debian.org|http://archive.debian.org|g' /etc/apt/sources.list && \
    sed -i '/stretch-updates/d' /etc/apt/sources.list

RUN apt-get update && \
   apt-get -y --no-install-recommends install --fix-missing \
   wget curl git build-essential libz-dev valgrind libatlas-base-dev gfortran \
   liblz4-tool gzip

# Install numpy first (scipy depends on it during build)
RUN pip install --no-cache 'numpy>=1.8.2,<1.17'

# Then install scipy
RUN pip install --no-cache scipy==1.2

WORKDIR /opt
RUN  git clone --recursive --depth 1 --shallow-submodules https://github.com/benedictpaten/sonLib/
RUN  git clone --recursive --depth 1 --shallow-submodules https://github.com/benedictpaten/pinchesAndCacti
RUN  git clone --recursive --depth 1 --shallow-submodules https://github.com/dentearl/mafTools/

# Capture mafTools version
WORKDIR /opt/mafTools
RUN  git rev-parse HEAD > /opt/maftools_version.txt && \
     echo "mafTools version: $(cat /opt/maftools_version.txt)"

WORKDIR /opt/sonLib
RUN  export PYTHONPATH=${PYTHONPATH}:/opt
RUN  make && make test || true

WORKDIR /opt/pinchesAndCacti
RUN  make && make test || true

WORKDIR /opt/mafTools
RUN  make
RUN  make

CMD exec "$@"
# docker run mafTools bablabalba

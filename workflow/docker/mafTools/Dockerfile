FROM python:2.7.18-stretch

ENV  LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV  PATH=/opt/mafTools/bin:$PATH
ENV  DEBIAN_FRONTEND=noninteractive
ENV  LC_ALL=C

RUN apt-get update && \
   apt-get -y --no-install-recommends install --fix-missing \
   wget curl git build-essential libz-dev valgrind libatlas-base-dev

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

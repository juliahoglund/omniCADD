FROM r-base

RUN export DOCKER_BUILDKIT=0
RUN export COMPOSE_DOCKER_CLI_BUILD=0

RUN apt-get update && \
   apt-get -y --no-install-recommends install --fix-missing \
   wget curl git build-essential libz-dev libatlas-base-dev \
   libcurl4-openssl-dev \
   libssl-dev libxml++2.6-dev librsvg2-dev cmake \
   libharfbuzz-dev libfribidi-dev

RUN R -e "install.packages('renv')"

WORKDIR /omniCADD
# change to from config later
COPY config/renv.lock renv.lock

ENV RENV_PATHS_LIBRARY config/renv/library

RUN R -e "renv::restore()"

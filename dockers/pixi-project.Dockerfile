# Dockerfile that containerizes a pixi project.
ARG ENV_NAME

# Use an ubuntu (release noble/2026.04) with pixi 0.82.1 pre-installed for the build image
# https://github.com/prefix-dev/pixi-docker/pkgs/container/pixi/1421279891?tag=0.82.1-noble
FROM ghcr.io/prefix-dev/pixi@sha256:f83a5fc6d1bb9a2ced27faf2b26c04c5eb618eed22fe13eb60f9225db4817ce0 AS build
ARG ENV_NAME

# git is needed to resolve git+https pypi dependencies during pixi install
RUN apt-get update && apt-get install -y --no-install-recommends git && rm -rf /var/lib/apt/lists/*

# create hack version of date to fix nextflow bug https://github.com/nextflow-io/nextflow/issues/7114
RUN echo '#!/usr/bin/env bash' > /usr/local/bin/date && \
    echo 'if [ "$1" = "+%s%3N" ]; then /usr/bin/date +%s%N | cut -c1-13 ; else /usr/bin/date "$@" ; fi' >> /usr/local/bin/date && \
    chmod a+x /usr/local/bin/date

# Install the project into the ENV_NAME
COPY . "/opt/$ENV_NAME"
WORKDIR "/opt/$ENV_NAME"

# install the packages into this environment
RUN pixi install && ls -lah

# create shell hook
RUN pixi shell-hook > /opt/shell-hook.sh
# extend the shell-hook script to run the command passed to the container
RUN echo 'exec "${@}"' >> /opt/shell-hook.sh

# ubuntu 26.04
# https://hub.docker.com/layers/library/ubuntu/26.04/images/
FROM ubuntu@sha256:c098726dd9b9733f3e885f8de0ef375206a311b0a317eadc40d347b1115fb944
ARG ENV_NAME
WORKDIR "/opt/$ENV_NAME"

# Setup a non-root user
RUN groupadd --system --gid 999 nonroot \
    && useradd --system --gid 999 --uid 999 --create-home nonroot

# use the nonroot user
USER nonroot

# for apptainer, need to update path
ENV PATH="/opt/$ENV_NAME/.pixi/envs/default/bin:$PATH"

# copy temporary fix to bug with nextflow date on rust coreutils
COPY --from=build /usr/local/bin/date /usr/local/bin/date
# copy the shell-hook so the environment can auto-activate
COPY --from=build --chown=nonroot:nonroot /opt/shell-hook.sh /opt/shell-hook.sh

# set the entrypoint to the shell-hook script (activate the environment and run the command)
ENTRYPOINT ["/bin/bash", "/opt/shell-hook.sh"]

# copy the environments into prod container
COPY --from=build --chown=nonroot:nonroot "/opt/$ENV_NAME/" "/opt/$ENV_NAME/"

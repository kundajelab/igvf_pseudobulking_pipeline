# Dockerfile that uses pixi to containerize a conda environment. In addition to being only slightly
# more complicated than using conda, it has the advantage of using .lock files for cross-platform
# reproducibility.
# This docker expects an external script to be maintaining the .toml and .lock files using pixi and
# a base .yaml
ARG ENV_NAME

# Use an ubuntu (release noble/2026.04) with pixi 0.70.1 pre-installed for the build image
# https://github.com/prefix-dev/pixi-docker/pkgs/container/pixi/914261862?tag=0.70.1-noble
FROM ghcr.io/prefix-dev/pixi@sha256:2537738f8b7e2c7a7f070f56928ab959c4559a8d7e04f71eb16b0f779f0588f6 AS build
ARG ENV_NAME

# Install the project into the ENV_NAME
WORKDIR "/opt/$ENV_NAME"

# install the packages into this environment
RUN --mount=type=cache,target=/root/.cache/pixi \
    --mount=type=bind,source="$ENV_NAME.lock",target="/opt/$ENV_NAME/pixi.lock" \
    --mount=type=bind,source="$ENV_NAME.toml",target="/opt/$ENV_NAME/pixi.toml" \
    pixi install

# create shell hook
RUN --mount=type=cache,target=/root/.cache/pixi \
    --mount=type=bind,source="$ENV_NAME.toml",target="/opt/$ENV_NAME/pixi.toml" \
    pixi shell-hook > /opt/shell-hook.sh

# extend the shell-hook script to run the command passed to the container
RUN echo 'exec "$@"' >> /opt/shell-hook.sh

# create hack version of date to fix nextflow bug https://github.com/nextflow-io/nextflow/issues/7114
RUN echo '#!/usr/bin/env bash' > /usr/local/bin/date && \
    echo 'if [ "$1" = "+%s%3N" ]; then /usr/bin/date +%s%N | cut -c1-13 ; else /usr/bin/date "$@" ; fi' >> /usr/local/bin/date && \
    chmod a+x /usr/local/bin/date


# ubuntu 26.04
# https://hub.docker.com/layers/library/ubuntu/26.04/images/sha256-37a0633b900e99d0937986022b4b4018908e1361704425313cc1c348bebd7230
FROM ubuntu@sha256:f3d28607ddd78734bb7f71f117f3c6706c666b8b76cbff7c9ff6e5718d46ff64 AS production
ARG ENV_NAME
WORKDIR "/opt/$ENV_NAME"

# Setup a non-root user
RUN groupadd --system --gid 999 nonroot \
    && useradd --system --gid 999 --uid 999 --create-home nonroot

# copy the environments into prod container
COPY --from=build --chown=nonroot:nonroot "/opt/$ENV_NAME/.pixi/envs" "/opt/$ENV_NAME/.pixi/envs"
# copy the shell-hook so the environment can auto-activate
COPY --from=build --chown=nonroot:nonroot /opt/shell-hook.sh /opt/shell-hook.sh
# copy temporary fix to bug with nextflow date on rust coreutils
COPY --from=build /usr/local/bin/date /usr/local/bin/date

# use the nonroot user
USER nonroot

# set the entrypoint to the shell-hook script (activate the environment and run the command)
# no more pixi needed in the prod container
ENTRYPOINT ["/bin/bash", "/opt/shell-hook.sh"]

# for apptainer, need to update path
ENV PATH="/opt/$ENV_NAME/.pixi/envs/default/bin:$PATH"

# Dockerfile that containerizes a uv project.
ARG ENV_NAME

# Use a debian image with uv version 0.13.4 pre-installed
# https://github.com/astral-sh/uv/pkgs/container/uv/1416269754?tag=0.13.4-trixie-slim
FROM ghcr.io/astral-sh/uv@sha256:2d5298090a57b670458ccae215742352f753f3d583db23f3cd35081d26859ada AS build
ARG ENV_NAME

# create hack version of date to fix nextflow bug https://github.com/nextflow-io/nextflow/issues/7114
RUN echo '#!/usr/bin/env bash' > /usr/local/bin/date && \
    echo 'if [ "$1" = "+%s%3N" ]; then /usr/bin/date +%s%N | cut -c1-13 ; else /usr/bin/date "$@" ; fi' >> /usr/local/bin/date && \
    chmod a+x /usr/local/bin/date

RUN apt-get update && apt-get install -y --no-install-recommends build-essential && rm -rf /var/lib/apt/lists/*

# Install the project into the ENV_NAME
WORKDIR "/opt/$ENV_NAME"

# Copy from the cache instead of linking since it's a mounted volume
ENV UV_LINK_MODE=copy

# Omit development dependencies
ENV UV_NO_DEV=1

# Ensure installed tools can be executed out of the box
ENV UV_TOOL_BIN_DIR=/usr/local/bin

# Install the project's dependencies using the lockfile and settings
COPY uv.lock /opt/$ENV_NAME/
COPY pyproject.toml /opt/$ENV_NAME/
COPY .python-version /opt/$ENV_NAME/
RUN  uv sync --locked --no-dev --no-install-project

# Then, add the rest of the project source code and install it
# Installing separately from its dependencies allows optimal layer caching
COPY --exclude=uv.lock --exclude=pyproject.toml . "/opt/$ENV_NAME"
RUN uv build && \
    uv sync --locked --no-dev

# create a shell-hook to activate
RUN echo source "/opt/$ENV_NAME/.venv/bin/activate" > /opt/shell-hook.sh
# extend the shell-hook script to run the command passed to the container
RUN echo 'exec "$@"' >> /opt/shell-hook.sh

# make uv python symlink point to /opt/local
RUN PYTHON_LOC=$(readlink -f "/opt/$ENV_NAME/.venv/bin/python" | sed 's,^/root/.local,/opt/local,') && \
    ln -sf "$PYTHON_LOC" "/opt/$ENV_NAME/.venv/bin/python"

# ubuntu 26.04
# https://hub.docker.com/layers/library/ubuntu/26.04/images/
FROM ubuntu@sha256:c098726dd9b9733f3e885f8de0ef375206a311b0a317eadc40d347b1115fb944 AS production
ARG ENV_NAME
WORKDIR "/opt/$ENV_NAME"

# Setup a non-root user
RUN groupadd --system --gid 999 nonroot \
    && useradd --system --gid 999 --uid 999 --create-home nonroot

# use the nonroot user
USER nonroot

# for apptainer, need to update path
ENV PATH="/opt/$ENV_NAME/.venv/bin:$PATH"

# env variables for packages that want to write to / or /root by default:
ENV NUMBA_CACHE_DIR=/tmp/numba
ENV MPLCONFIGDIR=/tmp/matplotlib

# copy temporary fix to bug with nextflow date on rust coreutils
COPY --from=build /usr/local/bin/date /usr/local/bin/date
# Copy shell-hook so that the environment can auto-activate
COPY --from=build --chown=nonroot:nonroot /opt/shell-hook.sh /opt/shell-hook.sh

# set the entrypoint to the shell-hook script (activate the environment and run the command)
ENTRYPOINT ["/bin/bash", "/opt/shell-hook.sh"]

# Copy the python environment, to /opt/local
COPY --from=build --chown=nonroot:nonroot /root/.local /opt/local

# Copy whole env into the production image. We're doing science so debug-ability is more important
# than maximum hardening
COPY --from=build --chown=nonroot:nonroot "/opt/$ENV_NAME" "/opt/$ENV_NAME"

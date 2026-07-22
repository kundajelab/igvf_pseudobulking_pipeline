#!/usr/bin/env bash
export JAVA_HOME="$PIXI_PROJECT_ROOT/.pixi/envs/default/lib/jvm"
"alias conda=mamba"
# Note: quoting the eval causes errors, it must remain unquoted.
# shellcheck disable=SC2046
eval $(pixi shell-hook)
if [[ ! -f "${CONDA_PREFIX}/bin/conda" ]]; then
    pushd "${CONDA_PREFIX}/bin" &> /dev/null || exit 1
    ln -s mamba conda
fi

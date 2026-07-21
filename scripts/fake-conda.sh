#!/usr/bin/env bash
"alias conda=mamba"
eval $(pixi shell-hook)
pushd "${CONDA_PREFIX}/bin" &> /dev/null
if [[ ! -f conda ]]; then
    ln -s mamba conda
fi
# maybe causes crashes?
#export LD_LIBRARY_PATH="$PIXI_PROJECT_ROOT/.pixi/envs/$PIXI_ENVIRONMENT_NAME/lib:${LD_LIBRARY_PATH}"

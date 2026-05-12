#!/usr/bin/env bash
if [[ -n "$SCRATCH" ]]; then
    echo "$SCRATCH/workspace"
else
    echo "$PIXI_PROJECT_ROOT/workspace"
fi

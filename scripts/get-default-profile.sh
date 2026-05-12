#!/usr/bin/env bash


if command -v sbatch &> /dev/null; then
    # slurm is available, use slurm and apptainer
    profile="slurm,apptainer"
else
    # run locally: prefer docker but can fall back to conda
    if command -v docker &> /dev/null; then
        # run local docker, check if we're on an ARM mac
        if [[ $(/usr/bin/arch) == "arm64" ]]; then
            profile="local,docker,arm"
        else
            profile="local,docker"
        fi
    else
        profile="local,conda"
    fi
fi
echo "$profile"

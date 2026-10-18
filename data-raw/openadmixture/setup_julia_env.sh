#!/usr/bin/env bash
# Create the conda environment used by data-raw/openadmixture/regenerate_julia_fixtures.R to
# run OpenADMIXTURE.jl and store its results for tests/testthat/test_openadmixture_vs_julia.R
# (the tests themselves do not need Julia).
#
# Julia and all its packages live inside the conda environment (the
# conda-forge julia package points JULIA_DEPOT_PATH at
# $CONDA_PREFIX/share/julia), so nothing is written to ~/.julia.
#
# Usage (from the root of the package):
#   bash data-raw/openadmixture/setup_julia_env.sh
#
# OpenADMIXTURE.jl is pinned to the commit used for the port. The Julia
# version is pinned too: conda-forge's Julia 1.10.4 aborts inside libuv when
# downloading packages on recent Linux kernels, and conda-forge has no 1.11
# builds, so we use 1.12.7.
set -euo pipefail

ENV_NAME="${OPENADMIXTURE_JULIA_ENV:-openadmixturer_julia}"
OPENADMIXTURE_JL_URL="https://github.com/OpenMendel/OpenADMIXTURE.jl"
OPENADMIXTURE_JL_REV="401ac5c4056f0f6f2a7ea6ffc3af69c2b1e56fdc"

source "$(conda info --base)/etc/profile.d/conda.sh"
if ! conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  conda create -y -n "$ENV_NAME" -c conda-forge "julia=1.12.7"
fi

conda run --no-capture-output -n "$ENV_NAME" julia -e "
  using Pkg
  Pkg.add(url = \"$OPENADMIXTURE_JL_URL\", rev = \"$OPENADMIXTURE_JL_REV\")
  Pkg.add([\"SnpArrays\", \"DelimitedFiles\", \"Random\"])
  Pkg.precompile()
  using OpenADMIXTURE
  println(\"OpenADMIXTURE.jl installed in conda environment $ENV_NAME\")
"

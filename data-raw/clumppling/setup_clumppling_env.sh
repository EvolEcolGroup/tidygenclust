#!/usr/bin/env bash
# Create the conda environment with the python implementation of clumppling
# (version 2.0.10), used by data-raw/clumppling/regenerate_python_fixtures.R
# and data-raw/clumppling/validate_vs_python.R to create the reference results
# of tests/testthat/test_*clumppling*.R (the tests themselves do not need
# python).
#
# Usage (from the root of the package):
#   bash data-raw/clumppling/setup_clumppling_env.sh
#
# The python packages are pinned in data-raw/clumppling/requirements.txt.
# Note that the 'infomap' package is not installed, so that cdlib uses the
# implementation of infomap in igraph, as the native port does.
set -euo pipefail

ENV_NAME="${CLUMPPLING_REFERENCE_ENV:-clumppling_reference}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "The conda environment '${ENV_NAME}' already exists; remove it with"
  echo "  conda env remove -n ${ENV_NAME}"
  echo "to create it again."
  exit 1
fi

conda create -y -n "${ENV_NAME}" -c conda-forge python=3.11 pip
conda run -n "${ENV_NAME}" python -m pip install \
  -r "${SCRIPT_DIR}/requirements.txt"
conda run -n "${ENV_NAME}" python -c \
  "import importlib.metadata as m; print('clumppling', m.version('clumppling'))"

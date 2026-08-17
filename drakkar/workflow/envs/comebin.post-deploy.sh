#!/usr/bin/env bash
# Snakemake post-deployment script for comebin.yaml.
# Runs once, inside the freshly created environment, right after conda has
# resolved it. Replaces the conda-resolved PyTorch with the exact CUDA-enabled
# build that COMEBin is known to work with on Mjolnir's GPU nodes. COMEBin 1.0.4
# runs on Python 3.7, so the last CUDA wheel available for that interpreter
# (1.10.2+cu113) is pinned by URL.

set -euo pipefail

# Avoid contamination from HPC software modules / user site-packages
unset PYTHONPATH || true
unset PYTHONHOME || true
export PYTHONNOUSERSITE=1

PYTHON="${CONDA_PREFIX}/bin/python"

echo "Installing CUDA-enabled PyTorch for COMEBin"

# --no-deps keeps pip from resolving the dependency tree against PyPI, which on
# Python 3.7 would pull incompatible replacements for conda-installed packages.
# The dependencies of torch 1.10.2 are already satisfied by the conda solve.
"${PYTHON}" -m pip install \
    --ignore-installed \
    --no-deps \
    "https://download.pytorch.org/whl/cu113/torch-1.10.2%2Bcu113-cp37-cp37m-linux_x86_64.whl"

echo "Checking PyTorch installation"

"${PYTHON}" - <<'PY'
import torch

print("PyTorch:", torch.__version__)
print("Torch location:", torch.__file__)
print("CUDA build:", torch.version.cuda)

assert torch.__version__ == "1.10.2+cu113", \
    "Unexpected PyTorch version: {}".format(torch.__version__)

assert torch.version.cuda == "11.3", \
    "Expected CUDA 11.3 build, got: {}".format(torch.version.cuda)

print("CUDA-enabled PyTorch installation OK")
PY

#!/usr/bin/env bash
# Snakemake post-deployment script for semibin.yaml.
# Runs once, inside the freshly created environment, right after conda has
# resolved it. Replaces the conda-resolved PyTorch with the exact CUDA-enabled
# build that SemiBin2 is known to work with on Mjolnir's A100 nodes.

set -euo pipefail

# Avoid contamination from HPC software modules / user site-packages
unset PYTHONPATH || true
unset PYTHONHOME || true
export PYTHONNOUSERSITE=1

PYTHON="${CONDA_PREFIX}/bin/python"

echo "Installing CUDA-enabled PyTorch for SemiBin2"

"${PYTHON}" -m pip install \
    --force-reinstall \
    torch==2.5.1 \
    --index-url https://download.pytorch.org/whl/cu124

echo "Checking PyTorch installation"

"${PYTHON}" - <<'PY'
import torch

print("PyTorch:", torch.__version__)
print("Torch location:", torch.__file__)
print("CUDA build:", torch.version.cuda)

assert torch.__version__.startswith("2.5.1"), \
    "Unexpected PyTorch version: {}".format(torch.__version__)

assert torch.version.cuda == "12.4", \
    "Expected CUDA 12.4 build, got: {}".format(torch.version.cuda)

print("CUDA-enabled PyTorch installation OK")
PY

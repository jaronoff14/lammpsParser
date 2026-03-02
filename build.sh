#!/usr/bin/env bash
set -euo pipefail
if [ -z "${CONDA_PREFIX:-}" ]; then
  echo "Warning: CONDA_PREFIX not set. Make sure the right Python is selected."
fi
mkdir -p build
cd build
cmake -DPython_EXECUTABLE="$(which python)" ..
cmake --build . --config Release
echo "Built. The module file will be in the build directory; to test run Python from this dir."

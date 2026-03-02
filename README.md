# lammpsParser

Starter repo for the LAMMPS dump parser C++ extension exposed to Python via pybind11.

This minimal repo includes:
- C++ parser + pybind11 bindings
- Column-major NumPy API for fast numeric access
- Convenience functions (install into active env)
- `pyproject.toml` for pip/scikit-build-based builds

See build instructions in this README.

---

## Quick build (local dev)

Create and activate a Python environment (conda recommended):

```bash
conda create -n pybind_tes python=3.10
conda activate pybind_tes
conda install -c conda-forge cmake pybind11 scikit-build-core numpy
# or use pip:
# pip install --upgrade pip
# pip install scikit-build-core cmake pybind11 numpy

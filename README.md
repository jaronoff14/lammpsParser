# lammpsParser

LAMMPS dump parser C++ extension exposed to Python via pybind11.



See build instructions in this README.

---

## Quick build (local dev)

Create and activate a Python environment (conda recommended):

```bash
conda create -n pybind_test python=3.10
conda activate pybind_test
conda install -c conda-forge cmake pybind11 scikit-build-core numpy
# or use pip:
# pip install --upgrade pip
# pip install scikit-build-core cmake pybind11 numpy


# MMTO

Multi-spectral source detection tool

[![DOI](https://img.shields.io/badge/DOI-10.1515/ACCESS-27cff5.svg)](https://doi.org/10.1109/ACCESS.2024.3403309)
[![Higra](https://img.shields.io/badge/Powered%20by-Higra-c98a22.svg)](https://higra.readthedocs.io/)
[![Astropy](https://img.shields.io/badge/powered%20by-Astropy-D84315.svg)](https://www.astropy.org/)

Installation
------------

**Requires a C++ 14 compiler and cmake**


```bash
python -m venv mmto                  # Create a virtual environment
source mmto/bin/activate             # Activate the virtual environment

pip install -U pip setuptools wheel scikit-build cmake ninja

pip install "numpy>=2.3" "higra>=0.6.12" "pybind11>=2.12"

# PIP_NO_BUILD_ISOLATION=1 pip install -v ./MMTO
pip install ./MMTO

```


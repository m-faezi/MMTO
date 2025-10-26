# MMTO

Multi-spectral source detection tool

[![DOI](https://img.shields.io/badge/DOI-10.1515/ACCESS-000.svg)](https://doi.org/10.1109/ACCESS.2024.3403309)
[![Higra](https://img.shields.io/badge/Powered%20by-Higra-000.svg)](https://higra.readthedocs.io/)
[![Astropy](https://img.shields.io/badge/powered%20by-Astropy-000.svg)](https://www.astropy.org/)

## Installation

```bash
python -m venv mmto
source mmto/bin/activate
pip install -r ./requirements/requirements_base.txt
pip install -r ./requirements/requirements_torch.txt || pip install -r ./requirements/requirements_torch_fallback.txt
pip install -U pip setuptools wheel scikit-build cmake ninja
pip install ./mmtolib
```

<p align="center">
    <img src="./assets/MMTO-pipeline.svg" alt="MMTO-pipeline" width="93.6%">
</p>

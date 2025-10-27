# <span style="color:#FF0000">M</span>MTO - <span  style="color:#666666">Multi-spectral Astronomical Source Detection Tool</span>

[![DOI](https://img.shields.io/badge/DOI-10.1515/ACCESS-000.svg)](https://doi.org/10.1109/ACCESS.2024.3403309)
[![Higra](https://img.shields.io/badge/Powered%20by-Higra-000.svg)](https://higra.readthedocs.io/)
[![Astropy](https://img.shields.io/badge/powered%20by-Astropy-000.svg)](https://www.astropy.org/)

MMTO ([Faezi et al.](#1)) is a multi-spectral photometric object detection and color extraction software, representing and processing on max-tree ([Salembier et al.](#2)) data structure across multiple spectral bands, built on [![MTO2](https://img.shields.io/badge/📡-MTO2-purple?style=for-the-badge&logo=github)](https://github.com/m-faezi/MTO2).

- **Multi-band processing**: Simultaneous analysis across multiple spectral bands
- **Cross-band matching**: Intelligent source association between bands using spatial and similarity metrics


## Installation

### Dependencies

The dependencies are listed in the [./requirements](requirements) directory.

```bash
python -m venv mmto
source mmto/bin/activate
pip install -r ./requirements/requirements_base.txt
pip install -r ./requirements/requirements_torch.txt || pip install -r ./requirements/requirements_torch_fallback.txt
pip install -U pip setuptools wheel scikit-build cmake ninja
pip install --no-build-isolation ./mmtolib
```

### Overview
MMTO extends multiple max-trees by integrating semantically meaningful node partitions, derived from statistical tests, into a structured graph. This integration enables the exploration of correlations among cross-band emissions, enhancing segmentation accuracy.
<p align="center">
    <img src="./assets/MMTO-pipeline.svg" alt="MMTO-pipeline" width="93.6%">
</p>

### Tuned run
```bash
python main.py --co_sim 0.9 --pix_dist 3.0 
```

### Command line arguments

| Option              | Description                                 | Type      | Default | Range/Values |
|---------------------|---------------------------------------------|-----------|---------|--------------|
| `--pix_dist`        | Cross-band central distance (pixel)         | float     | 3.00    | ≥ 0          |
| `--co_sim`          | Cross-band cosine similarity threshold      | float     | 0.90    | [0.0, 1.0)   |
| `-h`, `--help`      | Show the help message and exit              | flag      | -       | -            |


## Bibliography

- <a id="1">Faezi M. H., Peletier R., & Wilkinson M. H. (2024). “Multi-Spectral Source-Segmentation Using Semantically-Informed Max-Trees”. In: *IEEE Access* 12, pp. 72288 - 72302. DOI: [10.1109/ACCESS.2024.3403309](https://doi.org/10.1109/ACCESS.2024.3403309).</a>
- <a id="2">Salembier P., Oliveras A., & Garrido L. (1998). “Antiextensive connected operators for image and sequence processing”. In: *IEEE Transactions on Image Processing* 7.4, pp. 555–570. DOI: [10.1109/83.663500](https://doi.org/10.1109/83.663500).</a>

# <span style="color:#FF0000">M</span>MTO - <span  style="color:#666666">Multi-spectral Astronomical Source Detection Tool</span>

[![DOI](https://img.shields.io/badge/DOI-10.1515/ACCESS-000.svg)](https://doi.org/10.1109/ACCESS.2024.3403309)
[![Higra](https://img.shields.io/badge/Powered%20by-Higra-000.svg)](https://higra.readthedocs.io/)
[![Astropy](https://img.shields.io/badge/powered%20by-Astropy-000.svg)](https://www.astropy.org/)
[![contributions welcome](https://img.shields.io/badge/contributions-welcome-66ff00.svg)](https://github.com/m-faezi/MMTO/blob/main/CONTRIBUTING.md)
[![Matrix](https://img.shields.io/badge/Matrix-Join%20Chat-brightgreen?logo=matrix&logoColor=white)](https://matrix.to/#/%23MMTO:matrix.org)


MMTO ([Faezi et al.](#1)) is a multi-spectral photometric object detection and color extraction software, representing and processing on max-tree ([Salembier et al.](#2)) data structure across multiple spectral bands, built on [![MTO2](https://img.shields.io/badge/📡-MTO2-purple?style=flat&logo=github)](https://github.com/m-faezi/MTO2).
- **Multi-band processing**: Simultaneous analysis across multiple spectral bands
- **Cross-band matching**: Intelligent source association between bands using spatial and similarity metrics

<!-- omit in toc -->
## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
  - [Dependencies](#dependencies)
- [Usage](#usage)
  - [Configuration](#configuration)
  - [Tuned run](#tuned-run)
  - [Command line arguments](#command-line-arguments)
- [Citation](#citation)
- [Bibliography](#bibliography)
- [License](#license)

## Overview
MMTO extends multiple max-trees by integrating semantically meaningful node partitions, derived from statistical tests, into a structured graph. This integration enables the exploration of correlations among cross-band emissions, enhancing segmentation accuracy.
<p align="center">
    <img src="./assets/MMTO-pipeline.svg" alt="MMTO-pipeline" width="93.6%">
    <br>
    <em>MMTO</em> processing pipeline workflow.
</p>

## Installation

### Dependencies

The dependencies are listed in the [./requirements](requirements) directory.

```bash
python -m venv ./venvs/mmto
source ./venvs/mmto/bin/activate
pip install -r ./requirements/requirements_base.txt
pip install -r ./requirements/requirements_torch.txt || pip install -r ./requirements/requirements_torch_fallback.txt
pip install -U pip setuptools wheel scikit-build cmake ninja
pip install --no-build-isolation ./mmtolib
```

## Usage
> [!TIP]
> MMTO supports any number of bands to be processed.
> In this version, input images should have the same resolution (pixels) and sly coordinates.

### Configuration
Configure the bands and processing parameters (see [MTO2 documentation](https://github.com/m-faezi/MTO2)) in config.yaml file.
Here is an instance of 4-band setting:

```yaml
bands:
  
  band_1:
    file_path: "./1st.fits"
    background_mode: "const"
    move_factor: 5
    area_ratio: 0.90
    s_sigma: 4.0
    G_fit: false
    skip_reduction: true
  
  band_2:
    file_path: "./2nd.fits"
    background_mode: "morph"
    move_factor: 8
    area_ratio: 0.90
    s_sigma: 4.0
    G_fit: true
    skip_reduction: false
  
  band_3:
    file_path: "./3rd.fits"
    background_mode: "const"
    move_factor: 8
    area_ratio: 0.93
    s_sigma: 2.7
    G_fit: true
    skip_reduction: true
  
  band_4:
    file_path: "./4th.fits"
    background_mode: "const"
    move_factor: 3
    area_ratio: 0.93
    s_sigma: 3.13
    G_fit: false
    skip_reduction: true
```


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


## Citation

If you use MMTO in your research, please cite the following paper:

```bibtex
@ARTICLE{10535192,
  author={Hashem Faezi, Mohammad and Peletier, Reynier and Wilkinson, Michael H. F.},
  journal={IEEE Access}, 
  title={Multi-Spectral Source-Segmentation Using Semantically-Informed Max-Trees}, 
  year={2024},
  volume={12},
  number={},
  pages={72288-72302},
  doi={10.1109/ACCESS.2024.3403309}
}
```


## Bibliography

- <a id="1">Faezi M. H., Peletier R., & Wilkinson M. H. (2024). “Multi-Spectral Source-Segmentation Using Semantically-Informed Max-Trees”. In: *IEEE Access* 12, pp. 72288 - 72302. DOI: [10.1109/ACCESS.2024.3403309](https://doi.org/10.1109/ACCESS.2024.3403309).</a>
- <a id="2">Salembier P., Oliveras A., & Garrido L. (1998). “Antiextensive connected operators for image and sequence processing”. In: *IEEE Transactions on Image Processing* 7.4, pp. 555–570. DOI: [10.1109/83.663500](https://doi.org/10.1109/83.663500).</a>

## License

This project is licensed under the [MIT License](https://opensource.org/license/MIT) - see the [LICENSE](LICENSE) file for details.


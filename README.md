# COMPASS

Main status:
[![Main status](https://gitlab.obspm.fr/cosmic-rtc/compass/badges/main/pipeline.svg)](https://gitlab.obspm.fr/cosmic-rtc/compass/commits/main)

Develop status:
[![Develop status](https://gitlab.obspm.fr/cosmic-rtc/compass/badges/develop/pipeline.svg)](https://gitlab.obspm.fr/cosmic-rtc/compass/commits/develop)
[![coverage report](https://gitlab.obspm.fr/cosmic-rtc/compass/badges/develop/coverage.svg)](https://cosmic-rtc.pages.obspm.fr/compass/coverage/index.html)

- [COMPASS](#compass)
  - [Citations](#citations)
  - [Overview](#overview)
    - [Hardware requirements](#hardware-requirements)
    - [Environment requirements](#environment-requirements)
  - [Installation](#installation)
    - [Install Miniforge3 with python3](#install-miniforge3-with-python3)
      - [setup .bashrc](#setup-bashrc)
      - [Download and installation](#download-and-installation)
    - [Install the platform](#install-the-platform)
      - [Download sources](#download-sources)
      - [Install dependencies (if not already done)](#install-dependencies-if-not-already-done)
      - [Install COMPASS](#install-compass)
  - [Contributing](#contributing)
  - [License](#license)

## Citations
If you use COMPASS in your research, please cite one of the following papers:

- [Ferreira, F. et al, “COMPASS: an efficient GPU-based simulation software for adaptive optics systems”, HPCS 2018](https://doi.org/10.1109/HPCS.2018.00043)
- [Ferreira, F. et al., “Real-time end-to-end AO simulations at ELT scale on multiple GPUs with the COMPASS platform”, SPIE 2018](https://doi.org/10.1117/12.2312593)
- [Gratadour, D. et al, "COMPASS: an efficient, scalable and versatile numerical platform for the development of ELT AO systems"](https://doi.org/10.1117/12.2056358)
- [Gratadour, D. et al, “GPUs for adaptive optics: simulations and real-time control”, SPIE, 2012](https://doi.org/10.1117/12.925723)

## Overview

The COMPASS platform is distributed as a single bundle of CArMA and SuTrA C++ / Cuda libraries and their Python extensions NAGA & SHESHA.

### Hardware requirements

The system must contain at least an x86 CPU and a CUDA capable GPU. list of compatible GPUs can be found here <http://www.nvidia.com/object/cuda_gpus.html>. Specific requirements apply to clusters (to be updated).

### Environment requirements

The system must be running a 64 bit distribution of Linux with the latest NVIDIA drivers and CUDA toolkit. The following installation instructions are valid if the default installation paths have been selected for these components.

## Installation

### Install Miniforge3 with python3

more info: <https://github.com/mamba-org/mamba>

#### setup .bashrc

```bashrc
export MAMBA_ROOT=$HOME/miniforge3
export PATH=$MAMBA_ROOT/bin:$PATH
```

#### Download and installation

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p $MAMBA_ROOT
mamba init
```

### Install the platform

#### Download sources

First check out the latest version from the svn repository :

```bash
git clone https://gitlab.obspm.fr/cosmic-rtc/compass.git
```

once there, you need to modify system variables in our .bashrc :

```bash
## CUDA default definitions
export CUDA_ROOT=/usr/local/cuda
export CUDA_INC_PATH=$CUDA_ROOT/include
export CUDA_LIB_PATH=$CUDA_ROOT/lib
export CUDA_LIB_PATH_64=$CUDA_ROOT/lib64
export PATH=$CUDA_ROOT/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_LIB_PATH_64:$CUDA_LIB_PATH:$LD_LIBRARY_PATH

#COMPASS default definitions
export COMPASS_ROOT=$HOME/compass
export COMPASS_INSTALL_ROOT=$COMPASS_ROOT/local
export SHESHA_ROOT=$COMPASS_ROOT/shesha
export LD_LIBRARY_PATH=$COMPASS_INSTALL_ROOT/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$NAGA_ROOT:$SHESHA_ROOT:$COMPASS_INSTALL_ROOT/python:$PYTHONPATH
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$COMPASS_INSTALL_ROOT/lib/pkgconfig
```

#### Install dependencies (if not already done)

```bash
cd $COMPASS_ROOT
mamba env create --file environment.yml
mamba activate compass
export VCPKG_ROOT=$HOME/vcpkg
export PATH=$VCPKG_ROOT:$PATH
./script/install_vcpkg.sh $VCPKG_ROOT
```

#### Install COMPASS

```bash
cd $COMPASS_ROOT
./compile_vcpkg.py
```

## Contributing

If you want to contribute to COMPASS, please read the [CONTRIBUTING.md](CONTRIBUTING.md) file.
Contributions must go through merge requests.

## License

COMPASS is distributed under the LGPLv3 license. See [LICENSE](LICENSE) for more information.
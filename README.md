# COMPASS

Master status:
[![Master status](https://gitlab.obspm.fr/compass/compass/badges/master/pipeline.svg)](https://gitlab.obspm.fr/compass/compass/commits/master)

Develop status:
[![Develop status](https://gitlab.obspm.fr/compass/compass/badges/develop/pipeline.svg)](https://gitlab.obspm.fr/compass/compass/commits/develop)
[![coverage report](https://gitlab.obspm.fr/compass/compass/badges/develop/coverage.svg)](https://compass.pages.obspm.fr/compass/coverage/index.html)

- [COMPASS](#compass)
  - [Overview](#overview)
    - [Hardware requirements](#hardware-requirements)
    - [Environment requirements](#environment-requirements)
  - [Install Anaconda with python3](#install-anaconda-with-python3)
    - [setup .bashrc](#setup-bashrc)
    - [Download and installation](#download-and-installation)
  - [Install the platform](#install-the-platform)
    - [Download sources](#download-sources)
    - [Install dependencies (if not already done)](#install-dependencies-if-not-already-done)
      - [Conan dependencies](#conan-dependencies)
    - [Install COMPASS](#install-compass)

## Overview

The COMPASS platform is distributed as a single bundle of CArMA and SuTrA C++ / Cuda libraries and their Python extensions NAGA & SHESHA.

### Hardware requirements

The system must contain at least an x86 CPU and a CUDA capable GPU. list of compatible GPUs can be found here <http://www.nvidia.com/object/cuda_gpus.html>. Specific requirements apply to clusters (to be updated).

### Environment requirements

The system must be running a 64 bit distribution of Linux with the latest NVIDIA drivers and CUDA toolkit. The following installation instructions are valid if the default installation paths have been selected for these components.

## Install Anaconda with python3

more info: <https://www.continuum.io/downloads#linux>

### setup .bashrc

```bashrc
export CONDA_ROOT=$HOME/miniconda3
export PATH=$CONDA_ROOT/bin:$PATH
```

### Download and installation

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $CONDA_ROOT
```
## Install the platform

### Download sources

First check out the latest version from the svn repository :

```bash
git clone https://gitlab.obspm.fr/compass/compass

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
export COMPASS_DO_HALF="OFF"  # set to ON if you want to use half precision RTC (needs SM>=60)
export NAGA_ROOT=$COMPASS_ROOT/naga
export SHESHA_ROOT=$COMPASS_ROOT/shesha
export LD_LIBRARY_PATH=$COMPASS_INSTALL_ROOT/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$NAGA_ROOT:$SHESHA_ROOT:$COMPASS_INSTALL_ROOT/python:$PYTHONPATH
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$COMPASS_INSTALL_ROOT/lib/pkgconfig
```

### Install dependencies (if not already done)

```bash
pip install -r requirements.txt
cd $COMPASS_ROOT
./install_dependencies.sh
```

### Install COMPASS

```bash
cd $COMPASS_ROOT
./compile.sh
```

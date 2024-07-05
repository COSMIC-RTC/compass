# Advanced Installation Guide

This document provides a detailed guide for setting up the project environment and executing build and test operations based on the `.gitlab-ci.yml` configuration.

## Prerequisites

### Packages required

For debian based systems, the following packages are required:

```bash
apt install -y git wget bzip2 pkg-config gfortran make vim doxygen  bash-completion curl zip unzip tar autoconf automake autoconf-archive
```	

For RedHat based systems, the following packages are required:

```bash
dnf install -y dnf-plugins-core epel-release
dnf config-manager --set-enabled powertools
dnf install -y git wget vim xorg-x11-xauth zip gvim clang curl zip unzip tar autoconf automake autoconf-archive perl-IPC-Cmd
```

### Conda/Mamba Installation

For conda installation, we recommend to use Mamba, a faster and more efficient package manager. To install Mamba, run the following commands:

```bash
export CONDA_ROOT=${HOME}/mambaforge
export PATH=${CONDA_ROOT}/bin:${PATH}
wget -O Mambaforge-Linux-x86_64.sh https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh -b -p ${CONDA_ROOT}
rm Mambaforge-Linux-x86_64.sh
```

Note that the `CONDA_ROOT` path can be adjusted as needed and `wget -O` can be replaced with `curl --output`.

### Conda/Mamba Environment Setup

Create a new environment named `compass` using the provided `environment.yml` file:

```bash
mamba env create -f environment.yml
```

### Vcpkg Installation

To install Vcpkg, run the following commands:

```bash
export VCPKG_ROOT=${HOME}/vcpkg
export PATH=${VCPKG_ROOT}:${PATH}
install_vcpkg.sh
```

### bashrc Configuration

Add the following lines to the `.bashrc` file:

```bash
export COMPASS_ROOT=${HOME}/compass
export SHESHA_ROOT=${COMPASS_ROOT}/shesha
export NAGA_ROOT=${COMPASS_ROOT}/naga
export VCPKG_ROOT=${HOME}/vcpkg
export COMPASS_DO_HALF=OFF
export PYTHONPATH=${SHESHA_ROOT}:${COMPASS_ROOT}/local/
export LD_LIBRARY_PATH=${COMPASS_ROOT}/local/lib
export PATH=${VCPKG_ROOT}:${PATH}
```

## Configuration Overview

### Environment Variables Setup

make sure to set the previous environment variables (defined in the `.bashrc` file):

### Build Stage

1. **Environment Activation**: Activate the `compass` environment: `conda activate compass`.
    
2. **Dependency Installation**: Compile Vcpkg dependencies with `compile_vcpkg.py`.

NOTE: if you update the dependencies, you need to remove the `build` directory and run `compile_vcpkg.py` again.

### Tests

1. **Unit Tests**: Run the unit tests with:
    - `pytest $COMPASS_ROOT/python_module/carmaWrap/test`
    - `pytest $COMPASS_ROOT/shesha/tests/pytest/rtc`
    - `pytest $COMPASS_ROOT/shesha/tests/pytest/supervisor`

2. **Integration Tests**: Run the integration tests with: `$COMPASS_ROOT/shesha/tests/checkCompass.sh`.
# COMPASS

Master status: [![Master status](https://gitlab.obspm.fr/compass/compass/badges/master/pipeline.svg)](https://gitlab.obspm.fr/compass/compass/commits/master)

Develop status: [![Develop status](https://gitlab.obspm.fr/compass/compass/badges/develop/pipeline.svg)](https://gitlab.obspm.fr/compass/compass/commits/develop)

- [COMPASS](#compass)
  - [Overview](#overview)
    - [Hardware requirements](#hardware-requirements)
    - [Environment requirements](#environment-requirements)
  - [Install Anaconda with python3](#install-anaconda-with-python3)
    - [setup .bashrc](#setup-bashrc)
    - [Download and installation](#download-and-installation)
  - [Install MAGMA](#install-magma)
    - [Why MAGMA?](#why-magma)
    - [Extraction](#extraction)
    - [Configure MAGMA with MKL & installation](#configure-magma-with-mkl--installation)
    - [Configure with Makefile](#configure-with-makefile)
    - [Configure with CMake (not working fine...)](#configure-with-cmake-not-working-fine)
    - [Tuning (not tested)](#tuning-not-tested)
  - [Install the platform](#install-the-platform)
    - [Download sources](#download-sources)
    - [Install dependencies (if not already done)](#install-dependencies-if-not-already-done)
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

## Install MAGMA

### Why MAGMA?

The MAGMA project aims to develop a dense linear algebra library similar to LAPACK but for heterogeneous/hybrid architectures, starting with current "Multicore+GPU" systems.

Unlike CULA, MAGMA propose a dense linear algebra library handling double for free.

But MAGMA needs a LAPACK and a BLAS implementation. Actually, we try two options : openBLAS (free, easy to install) and MKL (free, need a registration but better optimized on Intel processors)

### Extraction

MAGMA is available here : <http://icl.cs.utk.edu/magma/software/index.html>

extract the tgz file and go into the new directory

```bash
wget http://icl.cs.utk.edu/projectsfiles/magma/downloads/magma-2.5.0.tar.gz -O - | tar xz
cd magma-2.5.0
```

### Configure MAGMA with MKL & installation

Installation of dependencies using anaconda

```bash
conda install -y numpy mkl-include pyqtgraph ipython pyqt qt matplotlib astropy blaze h5py hdf5 pytest-html pandas scipy docopt tqdm tabulate
```

### Configure with Makefile

You have to create your own make.inc based on make.inc.openblas:

```bash
cp make.inc-examples/make.inc.mkl-gcc make.inc
sed -i -e 's:/intel64: -Wl,-rpath=$(CUDADIR)/lib64 -Wl,-rpath=$(MKLROOT)/lib:' make.inc
```

just compile the shared target (and test if you want)

```bash
export MKLROOT=$CONDA_ROOT
export CUDA_ROOT=/usr/local/cuda
export NCPUS=8
GPU_TARGET=sm_XX MKLROOT=$MKLROOT CUDADIR=$CUDA_ROOT make -j $NCPUS shared sparse-shared
```

Where:

- sm_XX is compatible with the [compute capability](http://www.nvidia.com/object/cuda_gpus.html). For example, sm_60 for Tesla Tesla P100
- NCPUS is the number of CPUs in your system

To install libraries and include files in a given prefix, run:

```bash
GPU_TARGET=sm_XX MKLROOT=$MKLROOT CUDADIR=$CUDA_ROOT make install prefix=$HOME/local/magma
```

### Configure with CMake (not working fine...)

```bash
wget https://gitlab.obspm.fr/snippets/30/raw -O CMakeLists.txt
wget https://gitlab.obspm.fr/snippets/31/raw -O magma.pc.in
mkdir build
cd build
export CUDA_ROOT=/usr/local/cuda
export NCPUS=8
cmake .. -DGPU_TARGET=sm_XX -DLAPACK_LIBRARIES=$CONDA_ROOT/lib/libmkl_core.so -DMKLROOT=$CONDA_ROOT -DCMAKE_INSTALL_PREFIX=$HOME/local/magma
make -j $NCPUS
make install
```

Where:

- sm_XX is compatible with the [compute capability](http://www.nvidia.com/object/cuda_gpus.html). For example, sm_60 for Tesla Tesla P100
- NCPUS is the number of CPUs in your system

Note: If your gcc/g++ is too recent, please specify

```
-DCMAKE_CXX_COMPILER=$CUDA_ROOT/bin/g++ -DCMAKE_C_COMPILER=$CUDA_ROOT/bin/gcc
```

### Tuning (not tested)

For multi-GPU functions, set $MAGMA_NUM_GPUS to set the number of GPUs to use.

For multi-core BLAS libraries, set $OMP_NUM_THREADS or $MKL_NUM_THREADS or $VECLIB_MAXIMUM_THREADS to set the number of CPU threads, depending on your BLAS library.

## Install the platform

### Download sources

First check out the latest version from the svn repository :

```bash
git clone https://gitlab.obspm.fr/compass/compass --recurse-submodules
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

#MAGMA definitions
export MAGMA_ROOT=$HOME/local/magma
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MAGMA_ROOT/lib
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$MAGMA_ROOT/lib/pkgconfig

#third party lib path
export CUB_ROOT=$COMPASS_ROOT/tplib/cub
export WYRM_ROOT=$COMPASS_ROOT/tplib/wyrm

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

conda install -y numpy pyqtgraph ipython pyqt qt matplotlib astropy blaze h5py hdf5 pytest-html pandas scipy docopt tqdm tabulate

```

### Install COMPASS

```bash
cd $COMPASS_ROOT
./compile.sh
```

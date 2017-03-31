source: <https://projets-lesia.obspm.fr/projects/compass/wiki/Install_the_platform>

# Install Anaconda with python2

more info: <https://www.continuum.io/downloads#linux>

## Download and installation

```bash
wget https://repo.continuum.io/archive/Anaconda2-4.2.0-Linux-x86_64.sh
bash ./Anaconda2-4.2.0-Linux-x86_64.sh
```

update your .batchrc,
add anaconda2/bin into the $PATH

## Add more packets

To avoid any incompatibility this python modules, it's highly recommended to use the gcc provided with anaconda:

```bash
conda config --add channels compass
conda install compass
```

# Install MAGMA

## Why MAGMA ?

The MAGMA project aims to develop a dense linear algebra library similar to LAPACK but for heterogeneous/hybrid architectures, starting with current "Multicore+GPU" systems.

Unlike CULA, MAGMA propose a dense linear algebra library handling double for free.

But MAGMA needs a LAPACK and a BLAS implementation. Actually, we try two options : openBLAS (free, easy to install) and MKL (free, need a registration but better optimized on Intel processors)

## Configure MAGMA with openBLAS

### Dependencies : openblas (<http://www.openblas.net>)

First, clone the GIT repository:

```bash
git clone https://github.com/xianyi/OpenBLAS.git
```

compile it:

```bash
cd OpenBLAS/
make
```

install it:

```bash
make install PREFIX=$HOME/local/openblas
```

add to you .bashrc:

```bash
export OPENBLAS_ROOT=$HOME/local/openblas
```

### extraction

MAGMA is available here : <http://icl.cs.utk.edu/magma/software/index.html>

extract the tgz file and go into the new directory

```bash
wget http://icl.cs.utk.edu/projectsfiles/magma/downloads/magma-2.2.0.tar.gz
tar xf magma-2.2.0.tar.gz
cd magma-2.2.0
```

### configuration

You have to create your own make.inc based on make.inc.openblas:

```bash
cp make.inc-examples/make.inc.openblas make.inc
```

example : please verify GPU_TARGET, OPENBLASDIR, CUDADIR

```Makefile
#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA (version 2.2.0) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date November 2016
#//////////////////////////////////////////////////////////////////////////////

# GPU_TARGET contains one or more of Fermi, Kepler, or Maxwell,
# to specify for which GPUs you want to compile MAGMA:
#     Fermi   - NVIDIA compute capability 2.x cards
#     Kepler  - NVIDIA compute capability 3.x cards
#     Maxwell - NVIDIA compute capability 5.x cards
#     Pascal  - NVIDIA compute capability 6.x cards
# The default is "Fermi Kepler".
# Note that NVIDIA no longer supports 1.x cards, as of CUDA 6.5.
# See http://developer.nvidia.com/cuda-gpus
#
GPU_TARGET ?= Pascal

# --------------------
# programs

CC        = gcc
CXX       = g++
NVCC      = nvcc
FORT      = gfortran

ARCH      = ar
ARCHFLAGS = cr
RANLIB    = ranlib

# --------------------
# flags

# Use -fPIC to make shared (.so) and static (.a) library;
# can be commented out if making only static library.
FPIC      = -fPIC

CFLAGS    = -O3 $(FPIC) -DNDEBUG -DADD_ -Wall -fopenmp
FFLAGS    = -O3 $(FPIC) -DNDEBUG -DADD_ -Wall -Wno-unused-dummy-argument
F90FLAGS  = -O3 $(FPIC) -DNDEBUG -DADD_ -Wall -Wno-unused-dummy-argument -x f95-cpp-input
NVCCFLAGS = -O3         -DNDEBUG -DADD_       -Xcompiler "$(FPIC)"
LDFLAGS   =     $(FPIC)                       -fopenmp

# C++11 (gcc >= 4.7) is not required, but has benefits like atomic operations
CXXFLAGS := $(CFLAGS) -std=c++11
CFLAGS   += -std=c99

# --------------------
# libraries

# gcc with OpenBLAS (includes LAPACK)
LIB       = -lopenblas

LIB      += -lcublas -lcusparse -lcudart -lcudadevrt

# --------------------
# directories

# define library directories preferably in your environment, or here.
OPENBLASDIR ?= $(HOME)/local/openblas
CUDADIR ?= /usr/local/cuda
-include make.check-openblas
-include make.check-cuda

LIBDIR    = -L$(CUDADIR)/lib64 \
            -L$(OPENBLASDIR)/lib

INC       = -I$(CUDADIR)/include \
            -I$(OPENBLASDIR)/include
```

## Configure MAGMA with MKL

### extraction

To download MKL, you have to create a account here : <https://registrationcenter.intel.com/RegCenter/NComForm.aspx?ProductID=1517>

extract l_ccompxe_2013_sp1.1.106.tgz and go into l_ccompxe_2013_sp1.1.106

install it with ./install_GUI.sh and add IPP stuff to default choices configuration

You have to create your own make.inc based on make.inc.mkl-gcc-ilp64:

example: please verify GPU_TARGET, MKLROOT, CUDADIR

```Makefile
#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA (version 2.1.0) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date August 2016
#//////////////////////////////////////////////////////////////////////////////

# GPU_TARGET contains one or more of Fermi, Kepler, or Maxwell,
# to specify for which GPUs you want to compile MAGMA:
#     Fermi   - NVIDIA compute capability 2.x cards
#     Kepler  - NVIDIA compute capability 3.x cards
#     Maxwell - NVIDIA compute capability 5.x cards
#     Pascal  - NVIDIA compute capability 6.x cards
# The default is "Fermi Kepler".
# Note that NVIDIA no longer supports 1.x cards, as of CUDA 6.5.
# See http://developer.nvidia.com/cuda-gpus
#
#GPU_TARGET ?= Fermi Kepler

# --------------------
# programs

CC        = icc
CXX       = icpc
NVCC      = nvcc
FORT      = ifort

ARCH      = ar
ARCHFLAGS = cr
RANLIB    = ranlib

# --------------------
# flags

# Use -fPIC to make shared (.so) and static (.a) library;
# can be commented out if making only static library.
FPIC      = -fPIC

CFLAGS    = -O3 $(FPIC) -openmp -DADD_ -Wall -Wshadow -DMAGMA_WITH_MKL
FFLAGS    = -O3 $(FPIC)         -DADD_ -warn all -warn nounused -nogen-interfaces
F90FLAGS  = -O3 $(FPIC)         -DADD_ -warn all -warn nounused
NVCCFLAGS = -O3                 -DADD_ -Xcompiler "$(FPIC) -Wall -Wno-unused-function"
LDFLAGS   =     $(FPIC) -openmp

# Defining MAGMA_ILP64 or MKL_ILP64 changes magma_int_t to int64_t in include/magma_types.h
CFLAGS    += -DMKL_ILP64
FFLAGS    += -integer-size 64
F90FLAGS  += -integer-size 64
NVCCFLAGS += -DMKL_ILP64

# Options to do extra checks for non-standard things like variable length arrays;
# it is safe to disable all these
CFLAGS   += -pedantic -Wno-long-long
#CFLAGS   += -Werror  # uncomment to ensure all warnings are dealt with

# C++11 (icc >= 13) is not required, but has benefits like atomic operations
CXXFLAGS := $(CFLAGS) -std=c++11
CFLAGS   += -std=c99

# --------------------
# libraries

# IMPORTANT: these link lines are for 64-bit int !!!!
# For regular 64-bit builds using 64-bit pointers and 32-bit int,
# use the lp64 library, not the ilp64 library. See make.inc.mkl-gcc or make.inc.mkl-icc.

# see MKL Link Advisor at http://software.intel.com/sites/products/mkl/
# icc with MKL 10.3, Intel OpenMP threads, 64-bit int
# note -DMAGMA_ILP64 or -DMKL_ILP64, and -integer-size 64 in FLAGS above
LIB       = -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lpthread -lstdc++ -lm

LIB      += -lcublas -lcusparse -lcudart

# --------------------
# directories

# define library directories preferably in your environment, or here.
# for MKL run, e.g.: source /opt/intel/composerxe/mkl/bin/mklvars.sh intel64
#MKLROOT ?= /opt/intel/composerxe/mkl
#CUDADIR ?= /usr/local/cuda
-include make.check-mkl
-include make.check-cuda

LIBDIR    = -L$(CUDADIR)/lib64 \
            -L$(MKLROOT)/lib/intel64

INC       = -I$(CUDADIR)/include \
            -I$(MKLROOT)/include
```

In this example, I use gcc but with MKL, you can use icc instead of gcc. In this case, you have to compile yorick with icc. For this, you have to change the CC flag in Make.cfg

## compilation and installation

### compilation

just compile the shared target (and test if you want)

```bash
make -j 8 shared sparse-shared
```

### installation

To install libraries and include files in a given prefix, run:

```bash
make install prefix=$HOME/local/magma
```

The default prefix is /usr/local/magma. You can also set prefix in make.inc. tuning (not tested)

For multi-GPU functions, set $MAGMA_NUM_GPUS to set the number of GPUs to use. For multi-core BLAS libraries, set $OMP_NUM_THREADS or $MKL_NUM_THREADS or $VECLIB_MAXIMUM_THREADS to set the number of CPU threads, depending on your BLAS library.

# Install the platform

The COMPASS platform is distributed as a single bundle of CArMA and SuTrA C++ / Cuda libraries and their Python extensions NAGA & SHESHA.

## Hardware requirements

The system must contain at least an x86 CPU and a CUDA capable GPU. list of compatible GPUs can be found here <http://www.nvidia.com/object/cuda_gpus.html>. Specific requirements apply to clusters (to be updated).

## Environment requirements

The system must be running a 64 bit distribution of Linux with the latest NVIDIA drivers and CUDA toolkit. The following installation instructions are valid if the default installation paths have been selected for these components.

Additionally, to benefit from the user-oriented features of the platform, Anaconda2 should be installed (<https://www.continuum.io/downloads#_unix>). In the last versions of compass (r608+), Yorick is no more supported.

## Installation process

First check out the latest version from the svn repository :

```bash
git clone https://gitlab.obspm.fr/compass/compass
```

then go in the newly created directory and then go into the branch "develop":

```bash
cd compass
git checkout develop
```

once there, you need to modify system variables in our .bashrc :

```bash
# CUDA default definitions
export CUDA_ROOT=$CUDA_ROOT #/usr/local/cuda
export CUDA_INC_PATH=$CUDA_ROOT/include
export CUDA_LIB_PATH=$CUDA_ROOT/lib
export CUDA_LIB_PATH_64=$CUDA_ROOT/lib64
export PATH=$CUDA_ROOT/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_LIB_PATH_64:$CUDA_LIB_PATH:$LD_LIBRARY_PATH
```

in this file, you also have to indicate the proper architecture of your GPU so as the compiler will generate the appropriate code.

```bash
export GENCODE="arch=compute_52,code=sm_52"
```

and change both 52 to your architecture : for instance a Tesla Fermi will have 2.0 computing capabilities so change 52 to 20, a Kepler GPU will have 3.0 or 3.5 (K20) computing capabilities, change 52 to 30 (or 35), a Maxwell GPU have 5.2 (M6000), a Pascal have 6.0 (P100). (more informations here: <https://developer.nvidia.com/cuda-gpus>)

If you are using CULA, you have to specify it:

```bash
# CULA default definitions
export CULA_ROOT= /usr/local/cula
export CULA_INC_PATH= $CULA_ROOT/include
export CULA_LIB_PATH= $CULA_ROOT/lib
export CULA_LIB_PATH_64= $CULA_ROOT/lib64
export LD_LIBRARY_PATH=$CULA_LIB_PATH_64:$CULA_LIB_PATH:$LD_LIBRARY_PATH
```

If you are using MAGMA, you have to specify it:

```bash
# MAGMA definitions (uncomment this line if MAGMA is installed)
export MAGMA_ROOT=$HOME/local/magma
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MAGMA_ROOT/lib
export PKG_CONFIG_PATH=$MAGMA_ROOT/lib/pkgconfig
```

Last variables to define:

```bash
export COMPASS_ROOT=/path/to/compass/trunk
export NAGA_ROOT=$COMPASS_ROOT/naga
export SHESHA_ROOT=$COMPASS_ROOT/shesha
export LD_LIBRARY_PATH=$COMPASS_ROOT/libcarma:$COMPASS_ROOT/libsutra:$LD_LIBRARY_PATH
```

At the end, you .bashrc shoud containts all those informations:

```bash
# conda default definitions
export CONDA_ROOT=/your/path/anaconda2
export PATH=$CONDA_ROOT/bin:$PATH

# CUDA default definitions
export CUDA_INC_PATH=$CUDA_ROOT/include
export CUDA_LIB_PATH=$CUDA_ROOT/lib
export CUDA_LIB_PATH_64=$CUDA_ROOT/lib64
export PATH=$CUDA_ROOT/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_LIB_PATH_64:$CUDA_LIB_PATH:$LD_LIBRARY_PATH
export GENCODE="arch=compute_52,code=sm_52"

# OPENBLAS definitions
export OPENBLAS_ROOT=$HOME/local/openblas
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OPENBLAS_ROOT/lib
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$OPENBLAS_ROOT/lib/pkgconfig

# MAGMA definitions
export MAGMA_ROOT=$HOME/local/magma
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MAGMA_ROOT/lib
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$MAGMA_ROOT/lib/pkgconfig

# COMPASS default definitions
export COMPASS_ROOT=/your/path/compass
export NAGA_ROOT=$COMPASS_ROOT/naga
export SHESHA_ROOT=$COMPASS_ROOT/shesha
export LD_LIBRARY_PATH=$COMPASS_ROOT/libcarma:$COMPASS_ROOT/libsutra:$LD_LIBRARY_PATH
```

Once this is done, you're ready to compile the whole library:

```bash
make clean all
```

If you did not get any error, CArMA, SuTrA, NAGA and SHESHA are now installed on your machine. You can check that everything is working by launching a GUI to test a simulation:

```bash
ipython -i $SHESHA_ROOT/widgets/widget_ao.py
```

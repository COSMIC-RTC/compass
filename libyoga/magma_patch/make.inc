#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA (version 1.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of	Colorado, Denver
#      November 2011
#
#      Contributed by: Allan Menezes (Ontario, Canada)
#//////////////////////////////////////////////////////////////////////////////

#
# GPU_TARGET specifies for which GPU you want to compile MAGMA
#      0: Tesla family
#      1: Fermi Family
#
GPU_TARGET = 0

CC        = gcc
NVCC      = nvcc
FORT      = gfortran

ARCH      = ar
ARCHFLAGS = cr
RANLIB    = ranlib

OPTS      = -O3 -DADD_ -fPIC
FOPTS     = -O3 -DADD_ -x f95-cpp-input -fPIC
NVOPTS    = --compiler-options -fno-strict-aliasing -DUNIX -O3 -DADD_ -m64 -Xcompiler -fPIC
LDOPTS    = -fPIC -Xlinker -zmuldefs

LIB       = -latlas -lcblas -lpthread -lcublas -lcudart -llapack -lm 

CUDADIR   = /usr/local/cuda

LIBDIR    = -L/usr/local/cuda/lib64 -L/usr/lib64
INC       = -I$(CUDADIR)/include

#LIBMAGMA     = $(MAGMA_DIR)/lib/magma.a
#LIBMAGMABLAS = $(MAGMA_DIR)/lib/magmablas.a

CC = gcc
NVCC = nvcc

LIB_NAME = libkp_kalman.a

LAM_KP_DIR = $(COMPASS_ROOT_DIR)/lam/kalman_CPU_GPU
LAM_KP_INTERFACE_DIR = $(LAM_KP_DIR)/compass_interface

#MLIB    = icc -O2 -Wall -mkl -c -I$(MAGMA_ROOT)/include -DHAVE_CUBLAS

CFLAGS = -O3 -Wall  -m64 -I$(MAGMA_ROOT)/include -DMAGMA_SETAFFINITY  -DUSE_MAGMA_1_4 -DADD_ -DUSE_MAGMA_WITH_MKL -DHAVE_CUBLAS  -I$(LAM_KP_DIR) -I$(LAM_KP_INTERFACE_DIR) -I$(COMPASS_ROOT_DIR)/libcarma/include.h -I$(COMPASS_ROOT_DIR)/libsutra/include.h -fPIC -I$(CUDA_ROOT)/include  #-I/home/usr.off/include

# We have some temporal problems with cuda thrust
CUFLAGS = -O3 -m64 -Xcompiler -fPIC -use_fast_math -arch=compute_30 -code=sm_30 -I$(CUDA_ROOT)/include -I $(COMPASS_ROOT_DIR)/libcarma/include.h  -DHAVE_CUBLAS -I$(MAGMA_ROOT)/include  #-I/home/usr.off/include 

ifeq ($(COMPILATION_LAM),double)
  CUFLAGS += -DKP_DOUBLE
  CFLAGS += -DKP_DOUBLE
endif

SRC_CPP := $(wildcard $(LAM_KP_DIR)/*.cpp)

SRC_CU := $(wildcard $(LAM_KP_DIR)/*.cu)


OBJ_CPP := $(patsubst %.cpp,%.o,$(SRC_CPP))
OBJ_CU := $(patsubst %.cu,%.cuo,$(SRC_CU))


LDFLAGS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -lcublas -lcudart -lstdc++ -lgfortran

all: $(LIB_NAME)


$(LAM_KP_DIR)/%.cuo: $(LAM_KP_DIR)/%.cu  
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(NVCC) $(CUFLAGS) -c -o $@ $^
	
$(LAM_KP_DIR)/%.o: $(LAM_KP_DIR)/%.cpp  
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(CC) $(CFLAGS) -c -o $@ $^

$(LIB_NAME): $(OBJ_CU) $(OBJ_CPP)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Linking     " $@
	#@$(CC) -shared -o $@ $^
	@ar -rv $@ $^
clean:
	@rm -f $(OBJ_CPP) $(OBJ_CU) $(LIB_NAME)

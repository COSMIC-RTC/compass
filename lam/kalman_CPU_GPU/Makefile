CC = gcc
NVCC = nvcc

LIB_NAME = libkp_kalman.a

LAM_KP_DIR = $(COMPASS_ROOT_DIR)/lam/kalman_CPU_GPU
LAM_KP_INTERFACE_DIR = $(LAM_KP_DIR)/compass_interface

#MLIB    = icc -O2 -Wall -mkl -c -I$(MAGMA_ROOT)/include -DHAVE_CUBLAS 
CFLAGS = -O2 -Wall  -m64 -I$(MAGMA_ROOT)/include -DHAVE_CUBLAS  -I$(LAM_KP_DIR) -I$(LAM_KP_INTERFACE_DIR) -I$(COMPASS_ROOT_DIR)/libcarma/include.h -I$(COMPASS_ROOT_DIR)/libsutra/include.h -fPIC -I/home/tgautrais/soft/cuda-6.5/include #-I/home/usr.off/include 
 

# We have some temporal problems with cuda thrust
CUFLAGS = -m64 -Xcompiler -fPIC -use_fast_math -arch=sm_35 -I/home/tgautrais/soft/cuda-6.5/include #-I/home/usr.off/include 

SRC_CPP := $(wildcard $(LAM_KP_DIR)/*.cpp)
SRC_CPP += $(LAM_KP_INTERFACE_DIR)/kp_carma_tools.cpp

SRC_CU := $(wildcard $(LAM_KP_DIR)/*.cu)


OBJ_CPP := $(patsubst %.cpp,%.o,$(SRC_CPP))
OBJ_CU := $(patsubst %.cu,%.cuo,$(SRC_CU))


LDFLAGS =  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm

all: $(LIB_NAME)


$(LAM_KP_DIR)/%.cuo: $(LAM_KP_DIR)/%.cu  
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(NVCC) $(CUFLAGS) $(LDFLAGS) -c -o $@ $^
	
$(LAM_KP_DIR)/%.o: $(LAM_KP_DIR)/%.cpp  
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(CC) $(CFLAGS) $(LDFLAGS) -c -o $@ $^

$(LIB_NAME): $(OBJ_CU) $(OBJ_CPP)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Linking     " $@
	#@$(CC) -shared -o $@ $^
	@ar -rv $@ $^
clean:
	@rm -f $(OBJ_CPP) $(OBJ_CU) $(LIB_NAME)



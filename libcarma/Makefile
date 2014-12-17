
CU_DIR=src.cu
CPP_DIR=src.cpp
CUH_DIR=include.cuh
H_DIR=include.h
OBJ_DIR=obj

LIB_NAME=libcarma.so

#listing sources and futur objects
SRC_CPP = $(wildcard $(CPP_DIR)/*.cpp) 
SRC_CU  = $(wildcard $(CU_DIR)/*.cu) 

OBJS := $(patsubst $(CPP_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_CPP)) 
OBJS += $(patsubst $(CU_DIR)/%.cu,$(OBJ_DIR)/%.cuo,$(SRC_CU)) 

LIB_DEPLIBS  = -fPIC -Wall -lstdc++ 
LIB_DEPLIBS += -L$(CUDA_LIB_PATH_64) -lcuda -lcudart -lcublas -lcufft -lcurand -lcusparse -lnvToolsExt  -lpthread

# set compiler (or rarely loader) flags specific to this package
LIB_CPPFLAGS = -c -fPIC -Wall -fno-common -I$(H_DIR) -I$(CUDA_INC_PATH)

# ------------------------------------- targets and rules for this package
# color for printf
#31=red, 32=green, 33=yellow,34=blue, 35=pink, 36=cyan, 37=white

MPICXX           = $(shell which mpicxx 2>/dev/null)
ifeq ($(MPICXX),)
  GCC ?= $(CC)
else
  GCC = $(MPICXX)
endif
GCC = gcc

NVCC = $(CUDA_ROOT)/bin/nvcc
LIB_NVCCFLAGS  = -m64 -Xcompiler -fPIC -c --compiler-options 
LIB_NVCCFLAGS += -I$(CUH_DIR)  -I$(H_DIR) -I$(Y_MAKEDIR)/include
LIB_NVCCFLAGS += -gencode $(GENCODE)

#-gencode arch=compute_11,code=sm_11 -gencode arch=compute_13,code=sm_13 
#-gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 
#-gencode arch=compute_35,code=sm_35 

# --ptxas-options=-v # #-arch sm_13 -fpermissive 
ifneq ($(DEBUG),)
 LIB_CPPFLAGS  += -g -DDEBUG
 LIB_NVCCFLAGS += -g -G -lineinfo -DDEBUG
else
 LIB_CPPFLAGS  += -O3
 LIB_NVCCFLAGS += -O3 
endif

#ifneq ($(CULA_ROOT),)
# LIB_CPPFLAGS  += -I$(CULA_INC_PATH) -DUSE_CULA
# LIB_NVCCFLAGS += -I$(CULA_INC_PATH) -DUSE_CULA
# LIB_DEPLIBS   += -L$(CULA_LIB_PATH_64) -lcula_lapack_basic 
#endif

ifneq ($(MAGMA_ROOT),)
 #use the same keyword as in make.inc
 LIB_CPPFLAGS  += -I$(MAGMA_ROOT)/include -DUSE_MAGMA -DUSE_MAGMA_WITH_MKL
 LIB_CPPFLAGS  += -DADD_ -DHAVE_CUBLAS -DMAGMA_SETAFFINITY #-DMAGMA_ILP64 -DMKL_ILP64
 #comment this line if you are using magma 1.5
# LIB_CPPFLAGS  += -DUSE_MAGMA_1_4
 LIB_NVCCFLAGS += -Xcompiler -fno-strict-aliasing
 LIB_DEPLIBS   += -L$(MAGMA_ROOT)/lib -lmagma 
endif

ifneq ($(COMPILATION_LAM),)
  LIB_DEPLIBS += -L$(NVIDIA_CURRENT_LIB_PATH_64) 
endif

all: $(LIB_NAME)

clean:
	@rm -f $(OBJS)
	@rm -f *.so

$(OBJ_DIR)/%.cuo: $(CU_DIR)/%.cu
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(NVCC) $(LIB_NVCCFLAGS) -o $@ $< 

$(OBJ_DIR)/%.o: $(CPP_DIR)/%.cpp
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(GCC) $(CFLAGS) $(LIB_CPPFLAGS) -c -o $@ $<

$(LIB_NAME): $(OBJS)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Linking       " $(LIB_NAME)
	@$(GCC) -shared $(OBJS) -o $@ $(LIB_DEPLIBS)


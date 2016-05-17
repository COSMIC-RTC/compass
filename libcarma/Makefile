
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

LIB_DEPLIBS  = -fPIC -shared -Wl,-rpath=$(CUDA_LIB_PATH_64) -L$(CUDA_LIB_PATH_64)
LIB_DEPLIBS += -lcuda -lcudart -lcublas -lcufft -lcurand -lcusparse -lnvToolsExt  -lpthread

# set compiler (or rarely loader) flags specific to this package
LIB_CPPFLAGS = -c -fPIC -Wall -fno-common -std=c++11 #-DUSE_UVA
LIB_CPPFLAGS+= -I$(H_DIR) -I$(CUDA_INC_PATH)

# ------------------------------------- targets and rules for this package
# color for printf
#31=red, 32=green, 33=yellow,34=blue, 35=pink, 36=cyan, 37=white

#MPICXX           = $(shell which mpicxx 2>/dev/null)
ICC              = $(shell which icc 2>/dev/null)

ifeq ($(MPICXX),)
ifeq ($(ICC),)
  GCC ?= $(CC)
else
  GCC ?= $(ICC)
endif
else
  GCC = $(MPICXX)
endif

#force GCC if needed
#GCC := gcc

NVCC = $(CUDA_ROOT)/bin/nvcc
LIB_NVCCFLAGS  = -m64 -Xcompiler -fPIC -use_fast_math -c -Xptxas -dlcm=ca --maxrregcount 32 -std=c++11
LIB_NVCCFLAGS += -gencode $(GENCODE)
LIB_NVCCFLAGS += -I$(CUH_DIR) -I$(H_DIR) -I$(CUDA_INC_PATH)

#-gencode arch=compute_11,code=sm_11 -gencode arch=compute_13,code=sm_13 
#-gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 
#-gencode arch=compute_35,code=sm_35 

# --ptxas-options=-v # #-arch sm_13 -fpermissive 
ifneq ($(DEBUG),)
 LIB_CPPFLAGS  += -g -DDEBUG
 LIB_NVCCFLAGS += -g -DDEBUG -G -lineinfo
else
 LIB_CPPFLAGS  += -O3
 LIB_NVCCFLAGS += -O3 
endif

#ifneq ($(CULA_ROOT),)
# LIB_CPPFLAGS  += -I$(CULA_INC_PATH) -DUSE_CULA
# LIB_NVCCFLAGS += -I$(CULA_INC_PATH) -DUSE_CULA
# LIB_DEPLIBS   += -L$(CULA_LIB_PATH_64) -lcula_lapack_basic 
#endif

MAGMA = $(shell pkg-config --exists magma && echo "1" || echo "0")
#$(info $(MAGMA))
ifeq ($(MAGMA),1)
 #use the same keyword as in make.inc
 LIB_CPPFLAGS  += -DUSE_MAGMA $(shell pkg-config --cflags magma)

 LIB_NVCCFLAGS += -Xcompiler -fno-strict-aliasing
 LIB_DEPLIBS   += -Wl,-rpath=$(MAGMA_ROOT)/lib
 LIB_DEPLIBS   += $(shell pkg-config --libs magma)
 ifneq ($(USE_MAGMA_SPARSE),)
  LIB_CPPFLAGS  += -DUSE_MAGMA_SPARSE
  LIB_DEPLIBS   += -lmagma_sparse
 endif
else
 $(warning "WARNING MAGMA not found, did you add $(MAGMA_ROOT)/lib/pkgconfig to PKG_CONFIG_PATH")
endif

ifneq ($(KSPARSE_ROOT),)
 #use the same keyword as in make.inc
 LIB_CPPFLAGS  += -I$(KSPARSE_ROOT)/include -DUSE_KSPARSE
 LIB_DEPLIBS   += -Wl,-rpath=$(KSPARSE_ROOT)/lib -L$(KSPARSE_ROOT)/lib -lksparse 
endif

#ifneq ($(COMPILATION_LAM),)
#  LIB_DEPLIBS += -L$(NVIDIA_CURRENT_LIB_PATH_64) 
#endif

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
	$(GCC) $(OBJS) -o $@ $(LIB_DEPLIBS)


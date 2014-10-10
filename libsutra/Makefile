
CU_DIR=src.cu
CPP_DIR=src.cpp
CUH_DIR=include.cuh
H_DIR=include.h
I_DIR=yorick
OBJ_DIR=obj

LIB_NAME=libsutra.so


SRC_CU  = $(wildcard $(CU_DIR)/*.cu) 
SRC_CPP = $(wildcard $(CPP_DIR)/*.cpp) 

ifneq ($(COMPILATION_LAM),) 
  LAM_KP_DIR = $(COMPASS_ROOT_DIR)/lam/kalman_CPU_GPU
  LAM_KP_INTERFACE_DIR =  $(LAM_KP_DIR)/compass_interface
#else
#  SRC_CPP = $(filter-out $(CPP_DIR)/sutra_controller_kalman.cpp, $(wildcard $(CPP_DIR)/*.cpp)) 
endif


CURED_DIR = $(COMPASS_ROOT_DIR)/aao/cured
CURED_INC_PATH = $(CURED_DIR)
CURED_CPP = $(wildcard $(CURED_DIR)/*.cpp)

OBJS := $(patsubst $(CPP_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_CPP)) 
OBJS += $(patsubst $(CU_DIR)/%.cu,$(OBJ_DIR)/%.cuo,$(SRC_CU)) 
OBJS += $(patsubst %.cpp,%.o,$(CURED_CPP))


LIB_DEPLIBS= -L$(CUDA_LIB_PATH_64) -lcudart -lcublas -lcufft -lstdc++ -lcurand
LIB_DEPLIBS += -L$(COMPASS_ROOT_DIR)/libcarma -lcarma


# set compiler (or rarely loader) flags specific to this package


LIB_CPPFLAGS = -c -fPIC -Wall -DUNIX -I$(H_DIR) -I$(COMPASS_ROOT_DIR)/libcarma/$(H_DIR) -I$(CUDA_INC_PATH) -I$(CURED_INC_PATH)

ifneq ($(DEBUG),)
  LIB_CPPFLAGS+= -g -DDEBUG
  NVCCFLAGS+= -g -DDEBUG -G -lineinfo
else
  LIB_CPPFLAGS+= -O3
  NVCCFLAGS+= -O3
endif

MPICXX           = $(shell hich mpicxx 2>/dev/null)
ifeq ($(MPICXX),)
  GCC ?= $(CC)
else
  GCC = $(MPICXX)
endif
GCC = $(CC)

NVCC = $(CUDA_ROOT)/bin/nvcc
NVCCFLAGS = -m64 -Xcompiler -fPIC -use_fast_math -c -I$(CUH_DIR)  -I$(H_DIR) \
	-I$(COMPASS_ROOT_DIR)/libcarma/$(CUH_DIR)  -I$(COMPASS_ROOT_DIR)/libcarma/$(H_DIR) -I$(CURED_INC_PATH) -gencode $(GENCODE)

ifneq ($(COMPILATION_LAM),)
  CFLAGS = -DCOMPILATION_LAM
  LIB_DEPLIBS += -L$(NVIDIA_CURRENT_LIB_PATH_64) -lcuda -lgsl
  LIB_DEPLIBS += -L$(MAGMA_ROOT)/lib -lmagma
  LIB_DEPLIBS += $(LAM_KP_DIR)/libkp_kalman.a
  LIB_CPPFLAGS += -m64 -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -I$(MAGMA_ROOT)/include -DHAVE_CUBLAS -I$(LAM_KP_DIR) -I$(LAM_KP_INTERFACE_DIR) 
endif



all: $(LIB_NAME)

clean:
	@rm -rf $(OBJS)
	@rm -f	*.so

$(OBJ_DIR)/%.cuo: $(CU_DIR)/%.cu
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(NVCC) $(NVCCFLAGS) -o $@ $< 

$(OBJ_DIR)/%.o: $(CPP_DIR)/%.cpp
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(GCC) $(CFLAGS) -DADD_ $(LIB_CPPFLAGS) -c -o $@ $< 

$(CURED_DIR)/%.o: $(CURED_DIR)/%.cpp
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(GCC) $(CFLAGS) -DADD_ $(LIB_CPPFLAGS) -c -o $@ $< 

$(LIB_NAME): $(OBJS)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Linking       " $(LIB_NAME)
	@$(GCC) -shared -o $@ $^ $(LIB_DEPLIBS)



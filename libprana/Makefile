CU_DIR=src.cu
CPP_DIR=src.cpp
CUH_DIR=include.cuh
H_DIR=include.h
I_DIR=yorick
OBJ_DIR=obj

LIB_NAME=libprana.so

SRC_CPP = $(wildcard $(CPP_DIR)/*.cpp) 
SRC_CU  = $(wildcard $(CU_DIR)/*.cu) 

OBJS := $(patsubst $(CPP_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_CPP)) 
OBJS += $(patsubst $(CU_DIR)/%.cu,$(OBJ_DIR)/%.cuo,$(SRC_CU)) 
OBJS += $(patsubst %.cpp,%.o,$(CURED_CPP))

LIB_DEPLIBS= -L$(CUDA_LIB_PATH_64) -lcudart -lcublas -lcufft -lstdc++ -lcurand
LIB_DEPLIBS += -L$(COMPASS_ROOT_DIR)/libcarma -lcarma

# set compiler (or rarely loader) flags specific to this package

LIB_CPPFLAGS = -c -fPIC -Wall -DUNIX -I$(H_DIR) -I$(COMPASS_ROOT_DIR)/libcarma/$(H_DIR) -I$(CUDA_INC_PATH)

ifneq ($(DEBUG),)
  LIB_CPPFLAGS+= -O0 -g -DDEBUG
else
  LIB_CPPFLAGS+= -O2
endif

MPICXX           = $(shell which mpicxx 2>/dev/null)
ifeq ($(MPICXX),)
  GCC ?= $(CC)
else
  GCC = $(MPICXX)
endif
GCC = $(CC)

NVCC = $(CUDA_ROOT)/bin/nvcc
NVCCFLAGS = -m64 -Xcompiler -fPIC -use_fast_math -c -I$(CUH_DIR)  -I$(H_DIR) \
	-I$(COMPASS_ROOT_DIR)/libcarma/$(CUH_DIR)  -I$(COMPASS_ROOT_DIR)/libcarma/$(H_DIR) -I$(CURED_INC_PATH) #-gencode $(GENCODE)

ifneq ($(DEBUG),)
  NVCCFLAGS+= -O0 -g -DDEBUG # -G -lineinfo
else
  NVCCFLAGS+= -O2
endif

all: $(LIB_NAME)

clean:
	@rm -rf $(OBJS)

$(OBJ_DIR)/%.cuo: $(CU_DIR)/%.cu
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(NVCC) $(NVCCFLAGS) -o $@ $< 

$(OBJ_DIR)/%.o: $(CPP_DIR)/%.cpp
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(GCC) $(CFLAGS) -DADD_ $(LIB_CPPFLAGS) -c -o $@ $< 

$(LIB_NAME): $(OBJS)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Linking       " $(LIB_NAME)
	@$(GCC) -shared -o $@ $^ $(LIB_DEPLIBS)



CU_DIR=src.cu
CPP_DIR=src.cpp
CUH_DIR=include.cuh
H_DIR=include.h
I_DIR=yorick
OBJ_DIR=obj

DEBUG=1

LIB_NAME=libsutra.so

LIB_OBJS= sutra_ao_utils.cuo sutra_phase.o sutra_turbu.o sutra_turbu.cuo sutra_lgs.cuo sutra_lgs.o sutra_target.cuo sutra_target.o sutra_wfs.cuo sutra_wfs.o sutra_centroider.cuo sutra_centroider.o sutra_controler.cuo sutra_controler.o sutra_rtc.o sutra_dm.o sutra_dm.cuo sutra_telemetry.o sutra_aotemplate.o sutra_kl.o sutra_kl.cuo sutra_aotemplate.cuo sutra_acquisim.o sutra_acquisim.cuo 
OBJS := $(addprefix $(OBJ_DIR)/,$(LIB_OBJS))

LIB_DEPLIBS= -L$(CUDA_LIB_PATH_64) -lcudart -lcublas -lcufft -lstdc++ -lcurand #-lcula -L/usr/local/cula/lib 
LIB_DEPLIBS += -L$(COMPASS_ROOT_DIR)/libcarma -lcarma

# set compiler (or rarely loader) flags specific to this package

LIB_CPPFLAGS = -c -fPIC -Wall -DUNIX -I$(H_DIR) -I$(COMPASS_ROOT_DIR)/libcarma/$(H_DIR) -I$(CUDA_INC_PATH) #-I/usr/local/cula/include
ifneq ($(DEBUG),)
LIB_CPPFLAGS+= -g -DDEBUG
endif

MPICXX           = $(shell which mpicxx 2>/dev/null)
ifeq ($(MPICXX),)
  GCC ?= $(CC)
else
  GCC = $(MPICXX)
endif

NVCC = $(CUDA_ROOT)/bin/nvcc
NVCCFLAGS = -m64 -Xcompiler -fPIC -use_fast_math -c -I$(CUH_DIR)  -I$(H_DIR) -I$(COMPASS_ROOT_DIR)/libcarma/$(CUH_DIR)  -I$(COMPASS_ROOT_DIR)/libcarma/$(H_DIR) #-gencode $(GENCODE)

#-gencode arch=compute_11,code=sm_11 -gencode arch=compute_13,code=sm_13 -gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 #-I/usr/local/cula/include # --ptxas-options=-v #-arch sm_13 --compiler-options -fpermissive 

all: $(LIB_NAME)

clean:
	@rm -rf $(OBJS)

$(OBJ_DIR)/%.cuo: $(CU_DIR)/%.cu
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(NVCC) $(NVCCFLAGS) -o $@ $< 

$(OBJ_DIR)/%.o: $(CPP_DIR)/%.cpp
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(GCC) $(CFLAGS) -DADD_ $(LIB_CPPFLAGS) -c -o $@ $< #@$(GCC) $(CFLAGS) -DADD_ $(LIB_CPPFLAGS) -c -o $@ $<

$(LIB_NAME): $(OBJS)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Linking       " $(LIB_NAME)
	@$(GCC) -shared $(LIB_DEPLIBS) -o $@ $^ #@$(GCC) -shared $(LIB_DEPLIBS) -o $@ $^ 


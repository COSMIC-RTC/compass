CU_DIR=src.cu
CPP_DIR=src.cpp
CUH_DIR=include.cuh
H_DIR=include.h
I_DIR=yorick
OBJ_DIR=obj

DEBUG=1

LIB_NAME=libyoga_ao.so

LIB_OBJS= yoga_ao_utils.cuo yoga_phase.o yoga_turbu.o yoga_turbu.cuo yoga_lgs.cuo yoga_lgs.o yoga_target.cuo yoga_target.o yoga_wfs.cuo yoga_wfs.o yoga_centroider.cuo yoga_centroider.o yoga_controler.cuo yoga_controler.o yoga_rtc.o yoga_dm.o yoga_dm.cuo yoga_telemetry.o yoga_aotemplate.o yoga_kl.o yoga_kl.cuo yoga_aotemplate.cuo #yoga_turbu.cuo 
OBJS := $(addprefix $(OBJ_DIR)/,$(LIB_OBJS))

LIB_DEPLIBS= -L$(CUDA_LIB_DIR) -lcudart -lcublas -lcufft -lstdc++ -lcurand #-lcula -L/usr/local/cula/lib 
#PKG_DEPLIBS += ../cudpp/lib/libcudpp_x86_64.a
LIB_DEPLIBS += -L$(YOGA_DIR)/libyoga -lyoga

# set compiler (or rarely loader) flags specific to this package

LIB_CPPFLAGS = -c -fPIC -Wall -DUNIX -I$(H_DIR) -I$(YOGA_DIR)/libyoga/$(H_DIR) -I$(CUDA_PATH)/include #-I/usr/local/cula/include
ifneq ($(DEBUG),)
LIB_CPPFLAGS+= -g -DDEBUG
endif

GCC ?= $(CC)
MPICXX           = $(shell which mpiCC)
ifeq ($(MPICXX),)
MPICXX = $(GCC)
endif

NVCC = /usr/local/cuda/bin/nvcc
NVCCFLAGS = -m64 -Xcompiler -fPIC -use_fast_math -c -I$(CUH_DIR)  -I$(H_DIR) -I../../libyoga/$(CUH_DIR)  -I../../libyoga/$(H_DIR) -gencode arch=compute_12,code=sm_12
#-gencode arch=compute_11,code=sm_11 -gencode arch=compute_13,code=sm_13 -gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 #-I/usr/local/cula/include #-I../cudpp/include # --ptxas-options=-v #-arch sm_13 --compiler-options -fpermissive 

all: $(LIB_NAME)

clean:
	@rm -rf $(OBJS)

$(OBJ_DIR)/%.cuo: $(CU_DIR)/%.cu
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(NVCC) $(NVCCFLAGS) -o $@ $< 

$(OBJ_DIR)/%.o: $(CPP_DIR)/%.cpp
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(MPICXX) $(CFLAGS) -DADD_ $(LIB_CPPFLAGS) -c -o $@ $< #@$(GCC) $(CFLAGS) -DADD_ $(LIB_CPPFLAGS) -c -o $@ $<

$(LIB_NAME): $(OBJS)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Linking       " $(LIB_NAME)
	@$(MPICXX) -shared $(LIB_DEPLIBS) -o $@ $^ #@$(GCC) -shared $(LIB_DEPLIBS) -o $@ $^ 


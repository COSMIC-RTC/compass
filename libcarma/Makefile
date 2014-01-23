
CU_DIR=src.cu
CPP_DIR=src.cpp
CUH_DIR=include.cuh
H_DIR=include.h
OBJ_DIR=obj

#DEBUG=1

LIB_NAME=libcarma.so

LIB_OBJS=carma_host_obj.o carma_obj.o carma_obj.cuo carma_utils.o carma_fft.o carma_rng.o carma_rng.cuo carma_transpose.cuo carma_sum.cuo carma_cublas.o convolutionFFT2D.cuo carma_fft_conv.o carma_fft_conv.cuo carma_context.o carma_streams.o carma_cula.o #shrUtils.o 
ifneq ($(MAGMA_PATH),)
 LIB_OBJS += carma_magma.o
endif
OBJS := $(addprefix $(OBJ_DIR)/,$(LIB_OBJS))

LIB_DEPLIBS= -fPIC -Wall -O0 -lstdc++ -L$(CUDA_LIB_PATH_64) -L$(CULA_LIB_PATH_64) -lcudart -lcublas -lcufft -lcurand -lcula_lapack_basic #-lcula_lapack_link -lcula_core 

ifneq ($(MAGMA_PATH),)
 LIB_DEPLIBS += -L$(MAGMA_PATH)/lib -lmagma 
 #LIB_DEPLIBS += -DADD_ -DGPUSHMEM=130 -Xlinker -zmuldefs
endif

# set compiler (or rarely loader) flags specific to this package
#LIB_CPPFLAGS= -Wall -g -O3 -DUNIX -DADD_ -DGPUSHMEM=130 -fPIC -Xlinker -zmuldefs -fno-common -I$(H_DIR) -I/usr/local/cuda/include 
LIB_CPPFLAGS= -c -fPIC -Wall -O0 -fno-common -I$(H_DIR) -I$(CUDA_INC_PATH) -I$(CULA_INC_PATH) -DUNIX
ifneq ($(DEBUG),)
LIB_CPPFLAGS+= -g -DDEBUG
endif
ifneq ($(MAGMA_PATH),)
 LIB_CPPFLAGS +=-I$(MAGMA_PATH)/include -DUSE_MAGMA -DHAVE_CUBLAS -DADD_ -DMAGMA_SETAFFINITY
endif

# ------------------------------------- targets and rules for this package
# color for printf
#31=red, 32=green, 33=yellow,34=blue, 35=pink, 36=cyan, 37=white

MPICXX           = $(shell which mpicxx 2>/dev/null)
ifeq ($(MPICXX),)
  GCC ?= $(CC)
else
  GCC = $(MPICXX)
endif

NVCC = $(CUDA_ROOT)/bin/nvcc
NVCCFLAGS = -m64 -Xcompiler -fPIC -c --compiler-options -I$(CUH_DIR)  -I$(H_DIR) -I$(Y_MAKEDIR)/include -I$(CULA_INC_PATH) -gencode $(GENCODE)

#-gencode arch=compute_11,code=sm_11 -gencode arch=compute_13,code=sm_13 -gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 

# --ptxas-options=-v # #-arch sm_13 -fpermissive 
ifneq ($(DEBUG),)
NVCCFLAGS+= -g -DDEBUG
endif
ifneq ($(MAGMA_PATH),)
 NVCCFLAGS += -DADD_ -I$(MAGMA_PATH)/include -Xcompiler -fno-strict-aliasing
 NVCCFLAGS +=-DGPUSHMEM=200
endif

all: $(LIB_NAME)

clean:
	@rm -f $(OBJS)

cleanmagma: 
	@printf '\033[36m%s\033[31m%s\033[m\n' "cleaning magma install"
	@rm -rf $(MAGMA_PATH)

magma: $(MAGMA_PATH)

$(MAGMA_PATH): magma_1.1.0
	@printf '\033[36m%s\033[31m%s\033[m\n' "Patching     " $<
	@cp magma_patch/make.inc $<
#	@cd $</src && patch -ts <../../magma_patch/clatrd.cpp.patch
#	@cd $</src && patch -ts <../../magma_patch/zlatrd.cpp.patch
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $<
	@cd $< && make lib
#	@cd $< && prefix=$(MAGMA_PATH) make install
	@mkdir -p $(MAGMA_PATH)/include $(MAGMA_PATH)/lib
	@cp -Rf $</include $(MAGMA_PATH)
	@cp -Rf $</lib $(MAGMA_PATH)
	@cp $</control/common.o $(MAGMA_PATH)/lib
	@cp $</COPYRIGHT $(MAGMA_PATH)
	@cp $</README $(MAGMA_PATH)
	@cp $</ReleaseNotes $(MAGMA_PATH)
	@rm -rf $<
	@rm -f $<.tar.gz

magma_%: 
	@printf '\033[36m%s\033[31m%s\033[m\n' "Retreiving      " $@
	@wget http://icl.cs.utk.edu/projectsfiles/magma/$@.tar.gz
	@printf '\033[36m%s\033[31m%s\033[m\n' "Unziping         " $@
	@tar xf $@.tar.gz

$(OBJ_DIR)/%.cuo: $(CU_DIR)/%.cu
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(NVCC) $(NVCCFLAGS) -o $@ $< 

$(OBJ_DIR)/%.o: $(CPP_DIR)/%.cpp
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(GCC) $(CFLAGS) $(LIB_CPPFLAGS) -c -o $@ $<

$(LIB_NAME): $(OBJS)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Linking       " $(LIB_NAME)
	@$(GCC) -shared $(LIB_DEPLIBS) $(OBJS) -o $@ 


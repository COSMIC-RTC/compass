
CU_DIR=src.cu
CPP_DIR=src.cpp
CUH_DIR=include.cuh
H_DIR=include.h
OBJ_DIR=obj

CUDPP_PATH= #cudpp
CUDPP_LIB= $(CUDPP_PATH)/lib/libcudpp_x86_64.a
MAGMA_PATH= #magma
MAGMA_LIB= $(MAGMA_PATH)/lib/libmagma.a $(MAGMA_PATH)/lib/libmagmablas.a

DEBUG=1
DEBUG=0

LIB_NAME=libyoga.so

LIB_OBJS=yoga_host_obj.o yoga_obj.o yoga_obj.cuo yoga_utils.o yoga_fft.o shrUtils.o yoga_rng.o yoga_rng.cuo yoga_transpose.cuo yoga_sum.cuo yoga_cublas.o convolutionFFT2D.cuo yoga_fft_conv.o yoga_fft_conv.cuo yoga_context.o yoga_streams.o yoga_cula.o yoga_multithread.o
ifneq ($(CUDPP_PATH),)
 LIB_OBJS += yoga_cudpp.o
endif
ifneq ($(MAGMA_PATH),)
 LIB_OBJS += yoga_magma.o
endif
OBJS := $(addprefix $(OBJ_DIR)/,$(LIB_OBJS))


LIB_DEPLIBS= -fPIC -Wall -O0 -lstdc++ -L/usr/local/cuda/lib64 -L$(CULA_LIB_PATH_64) -lcuda -lcudart -lcublas -lcufft -lcurand -lcula_lapack_basic #-lcula_lapack_link -lcula_core 
ifneq ($(CUDPP_PATH),)
 LIB_DEPLIBS +=-L$(CUDPP_PATH)/lib -lcudpp_x86_64 
endif
#-lcurand -lcula -L/usr/local/cula/lib 
ifneq ($(MAGMA_PATH),)
 LIB_DEPLIBS += -L$(MAGMA_PATH)/lib -lmagma 
 LIB_DEPLIBS += -latlas -lcblas -lpthread -llapack 
 LIB_DEPLIBS +=  -lmagmablas #-llapacktest
 LIB_DEPLIBS +=  $(MAGMA_PATH)/lib/common.o
 #LIB_DEPLIBS +=  $(MAGMA_PATH)/lib/libmagmablas.a $(MAGMA_PATH)/lib/liblapacktest.a
 LIB_DEPLIBS += -DADD_ -DGPUSHMEM=130 -Xlinker -zmuldefs
endif

# set compiler (or rarely loader) flags specific to this package
#LIB_CPPFLAGS= -Wall -g -O3 -DUNIX -DADD_ -DGPUSHMEM=130 -fPIC -Xlinker -zmuldefs -fno-common -I$(H_DIR) -I/usr/local/cuda/include 
LIB_CPPFLAGS= -c -fPIC -Wall -O0 -fno-common -I$(H_DIR) -I$(CUDA_PATH)/include -I$(CULA_INC_PATH) -DUNIX
ifneq ($(DEBUG),)
LIB_CPPFLAGS+= -g -DDEBUG
endif
ifneq ($(CUDPP_PATH),)
 LIB_CPPFLAGS +=-I$(CUDPP_PATH)/include -DUSE_CUDPP
endif
ifneq ($(MAGMA_PATH),)
 LIB_CPPFLAGS +=-I$(MAGMA_PATH)/include -DUSE_MAGMA
 LIB_CPPFLAGS +=-DADD_ -DGPUSHMEM=130
endif

# ------------------------------------- targets and rules for this package
# color for printf
#31=red, 32=green, 33=yellow,34=blue, 35=pink, 36=cyan, 37=white

GCC ?= $(CC)
NVCC = /usr/local/cuda/bin/nvcc
NVCCFLAGS = -m64 -Xcompiler -fPIC -c --compiler-options -I$(CUH_DIR)  -I$(H_DIR) -I$(Y_MAKEDIR)/include -I/usr/local/cula/include -gencode arch=compute_12,code=sm_12

#-gencode arch=compute_11,code=sm_11 -gencode arch=compute_13,code=sm_13 -gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 

# --ptxas-options=-v # #-arch sm_13 -fpermissive 
ifneq ($(DEBUG),)
NVCCFLAGS+= -g -DDEBUG
endif
ifneq ($(CUDPP_PATH),)
 NVCCFLAGS +=-I$(CUDPP_PATH)/include 
endif
ifneq ($(MAGMA_PATH),)
 NVCCFLAGS +=-I$(MAGMA_PATH)/include 
 NVCCFLAGS +=-DGPUSHMEM=130
endif

MPICXX           = $(shell which mpicxx)
ifeq ($(MPICXX),)
  MPICXX = GCC
  MPI=1
endif

all: $(LIB_NAME)

clean:
	@rm -f $(OBJS)

cudpp:$(CUDPP_LIB) 

$(CUDPP_LIB): cudpp_src_1.1.1 #cudpp_src_2.0 
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " "cupdd"
	@cd $</cudpp && CXX_ARCH_FLAGS="-fPIC" NVCCFLAGS="--compiler-options -Xcompiler -fPIC" make
#	@cd $< && cmake . && CXX_ARCH_FLAGS="-fPIC" NVCCFLAGS="--compiler-options -Xcompiler -fPIC" make
	@mkdir -p $(CUDPP_PATH)/include $(CUDPP_PATH)/lib
	@cp $</cudpp/include/* $(CUDPP_PATH)/include
#	@cp $</include/* $(CUDPP_PATH)/include
	@cp $</lib/* $(CUDPP_PATH)/lib
	@rm -rf $<
	@rm -f $<.zip

cudpp_src_%: 
	@printf '\033[36m%s\033[31m%s\033[m\n' "Retreiving     " $@
	@wget http://cudpp.googlecode.com/files/$@.zip
	@printf '\033[36m%s\033[31m%s\033[m\n' "Unziping       " $@
	@unzip $@.zip

cleanmagma: 
	@printf '\033[36m%s\033[31m%s\033[m\n' "cleaning magma install"
	@rm -rf $(MAGMA_PATH)

magma: $(MAGMA_LIB)

$(MAGMA_LIB): magma_1.1.0
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
	@$(MPICXX) $(CFLAGS) $(LIB_CPPFLAGS) -c -o $@ $<

$(LIB_NAME): $(OBJS)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Linking       " $(LIB_NAME)
	@$(MPICXX) -shared $(LIB_DEPLIBS) $(OBJS) -o $@ 


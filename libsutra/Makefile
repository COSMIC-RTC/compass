CU_DIR=src.cu
CPP_DIR=src.cpp
CUH_DIR=include.cuh
H_DIR=include.h
OBJ_DIR=obj

LIB_NAME=libsutra.so


SRC_CU  = $(wildcard $(CU_DIR)/*.cu)
SRC_CPP = $(wildcard $(CPP_DIR)/*.cpp)

ifneq ($(COMPILATION_LAM),)
  LAM_KALMAN_DIR = $(COMPASS_ROOT)/lam/kalman_CPU_GPU
  LAM_KP_CLASSES_DIR =  $(COMPASS_ROOT)/lam/kp_classes
  LAM_KP_INTERFACE_DIR =  $(COMPASS_ROOT)/lam/compass_interface
#else
#  SRC_CPP = $(filter-out $(CPP_DIR)/sutra_controller_kalman.cpp, $(wildcard $(CPP_DIR)/*.cpp))
endif


CURED_DIR = $(COMPASS_ROOT)/aao/cured
CURED_INC_PATH = $(CURED_DIR)
CURED_CPP = $(wildcard $(CURED_DIR)/*.cpp)

OBJS := $(patsubst $(CPP_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_CPP))
OBJS += $(patsubst $(CU_DIR)/%.cu,$(OBJ_DIR)/%.cuo,$(SRC_CU))
OBJS += $(patsubst %.cpp,%.o,$(CURED_CPP))

OS_SIZE = $(shell uname -m | sed -e "s/i.86/32/" -e "s/x86_64/64/" -e "s/armv7l/32/")
OS_ARCH = $(shell uname -m | sed -e "s/i386/i686/")

ifeq ($(OS_ARCH),armv7l)
  ifeq ($(abi),gnueabi)
    CCFLAGS += -mfloat-abi=softfp
  else
    # default to gnueabihf
    override abi := gnueabihf
    #      LDFLAGSS += --dynamic-linker=/lib/ld-linux-armhf.so.3
    CCFLAGS   += -mfloat-abi=hard
    NVCCFLAGS := -target-cpu-arch ARM -Xlinker --dynamic-linker=/lib/ld-linux-armhf.so.3 -Xcompiler -mfloat-abi=hard
  endif
else
  CCFLAGS += -m64
endif


LIB_DEPLIBS  = -Wl,-rpath=$(CUDA_LIB_PATH_64) -L$(CUDA_LIB_PATH_64)
LIB_DEPLIBS += -lcudart -lcublas -lcufft -lstdc++ -lcurand
LIB_DEPLIBS += -Wl,-rpath=$(COMPASS_ROOT)/libcarma -L$(COMPASS_ROOT)/libcarma -lcarma


# set compiler (or rarely loader) flags specific to this package


LIB_CPPFLAGS = -c -fPIC -Wall -Wno-unused-result -DUNIX -std=c++11 $(CCFLAGS)
LIB_CPPFLAGS+= -I$(H_DIR) -I$(COMPASS_ROOT)/libcarma/$(H_DIR) -I$(CUDA_INC_PATH) -I$(CURED_INC_PATH)

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

#GCC := icc #/usr/local/mpe/bin/mpecc -mpilog

NVCC = $(CUDA_ROOT)/bin/nvcc
LIB_NVCCFLAGS  = -Xcompiler -fPIC -use_fast_math -c -Xptxas -dlcm=ca --maxrregcount 32 -std=c++11 $(NVCCFLAGS)
LIB_NVCCFLAGS += -I$(CUH_DIR)  -I$(H_DIR) -I$(COMPASS_ROOT)/libcarma/$(CUH_DIR)
LIB_NVCCFLAGS += -I$(COMPASS_ROOT)/libcarma/$(H_DIR) -I$(CURED_INC_PATH) -gencode $(GENCODE)

ifneq ($(DEBUG),)
  LIB_CPPFLAGS+= -g -DDEBUG
  LIB_NVCCFLAGS+= -g -DDEBUG -G -lineinfo
else
  LIB_CPPFLAGS+= -O3
  LIB_NVCCFLAGS+= -O3
endif

ifneq ($(COMPILATION_LAM),)
  CFLAGS = -DCOMPILATION_LAM -O3 -DHAVE_CUBLAS -I$(MAGMA_ROOT)/include -I$(LAM_KP_CLASSES_DIR) -I$(LAM_KP_INTERFACE_DIR) -I$(LAM_KALMAN_DIR)
  LIB_DEPLIBS += -Wl,-rpath=$(NVIDIA_CURRENT_LIB_PATH_64) -L$(NVIDIA_CURRENT_LIB_PATH_64) -lcuda -lgsl
  LIB_DEPLIBS += -Wl,-rpath=$(MAGMA_ROOT)/lib -L$(MAGMA_ROOT)/lib -lmagma
  LIB_DEPLIBS += $(LAM_KP_CLASSES_DIR)/lib_kp_classes.a
  LIB_DEPLIBS += $(LAM_KALMAN_DIR)/lib_kalman.a
  LIB_CPPFLAGS += -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
  OBJS += $(LAM_KP_CLASSES_DIR)/lib_kp_classes.a  $(LAM_KALMAN_DIR)/lib_kalman.a
endif

ifeq ($(filter $(COMPILATION_LAM),double),double)
  CFLAGS += -DKP_DOUBLE
endif

ifneq ($(BRAMA_ROOT),)
	LIB_CPPFLAGS  += -I$(BRAMA_ROOT) -I$(ACE_ROOT) -I$(TAO_ROOT) -I$(TAO_ROOT)/orbsvcs -I$(DDS_ROOT) -DUSE_BRAMA
  LIB_DEPLIBS   += -Wl,-rpath=$(ACE_ROOT)/lib,-rpath=$(TAO_ROOT)/tao,-rpath=$(DDS_ROOT)/lib
# LIB_DEPLIBS   += -L$(ACE_ROOT)/lib -L$(TAO_ROOT)/tao -L$(DDS_ROOT)/lib
	LIB_DEPLIBS   += -Wl,-rpath=$(BRAMA_ROOT) -L$(BRAMA_ROOT) -lBRAMACommon
# LIB_DEPLIBS   += -Wl,-rpath=$(ACE_ROOT)/lib -L$(ACE_ROOT)/lib
#	LIB_DEPLIBS 	+= -lTAO_Svc_Utils -lOpenDDS_InfoRepoDiscovery -lOpenDDS_Dcps -lTAO_PortableServer -lTAO_AnyTypeCode -lTAO -lACE
  LIB_CPPFLAGS  += -fvisibility-inlines-hidden -fno-strict-aliasing -Wno-deprecated-declarations -pthread -pipe -D_GNU_SOURCE -D__ACE_INLINE__
  LIB_DEPLIBS   += -Wl,-E -fvisibility-inlines-hidden -fno-strict-aliasing -pthread -pipe -L$(ACE_ROOT)/lib -L$(DDS_ROOT)/lib -L$(TAO_ROOT)/tao -lBRAMACommon -lOpenDDS_InfoRepoDiscovery -lOpenDDS_Dcps -lTAO_PortableServer -lTAO_AnyTypeCode -lTAO -lACE -ldl -lrt
endif

.PHONY: ${OBJ_DIR}

all: ${OBJ_DIR} $(LIB_NAME)

${OBJ_DIR}:
	mkdir -p ${OBJ_DIR}

clean:
	@rm -rf $(OBJS)
	@rm -f	*.so

$(OBJ_DIR)/%.cuo: $(CU_DIR)/%.cu
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(NVCC) $(LIB_NVCCFLAGS) -o $@ $<

$(OBJ_DIR)/%.o: $(CPP_DIR)/%.cpp
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(GCC) $(CFLAGS) -DADD_ $(LIB_CPPFLAGS) -c -o $@ $<

$(CURED_DIR)/%.o: $(CURED_DIR)/%.cpp
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(GCC) $(CFLAGS) -DADD_ $(LIB_CPPFLAGS) -c -o $@ $<

$(LIB_NAME): $(OBJS)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Linking       " $@
	@$(GCC) -shared -o $@ $^ $(LIB_DEPLIBS)

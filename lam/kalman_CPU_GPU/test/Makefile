

CXX=g++


SRC_CPP = kp_init4matio.cpp  

KALMAN_DIR = $(COMPASS_ROOT_DIR)/lam/kalman_CPU_GPU
KP_DIR = $(COMPASS_ROOT_DIR)/lam/kp_classes

LIB_KALMAN =  $(KALMAN_DIR)/lib_kalman.a
LIB_KP =  $(KP_DIR)/lib_kp_classes.a


CARMA_INC =  $(COMPASS_ROOT_DIR)/libcarma/include.h

ifeq ($(CXX),icc)
	CXXFLAGS_MKL = -mkl
else
	CXXFLAGS_MKL = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
endif

CXXFLAGS_CPU = -m64 -O2 -Wall $(CXXFLAGS_MKL) $(LIBS) -I$(KALMAN_DIR) -I$(KP_DIR) -g -I$(CUDA_INC_PATH) -I$(CARMA_INC)
ifeq ($(filter $(COMPILATION_LAM),double),double)
  CXXFLAGS_CPU += -DKP_DOUBLE
endif

CXXFLAGS_GPU = $(CXXFLAGS_CPU)  -I$(MAGMA_ROOT)/include -DHAVE_CUBLAS

LDFLAGS_CPU =  -L$(CUDA_LIB_PATH_64) -lmatio -lgsl -lcudart -lcublas -lcusparse -L$(MAGMA_ROOT)/lib -lmagma
LDFLAGS_GPU = $(LDFLAGS_CPU) -L$(MAGMA_ROOT)/lib -lcublas -lcusparse -lcudart -lmagma


all:kalman_CPU kalman_GPU





kalman_CPU: kalman_CPU.cpp $(SRC_CPP) $(LIB_KALMAN) $(LIB_KP)
	$(CXX) $^ $(CXXFLAGS_CPU) $(LDFLAGS_CPU) -o $@

kalman_GPU: kalman_GPU.cpp $(SRC_CPP) $(LIB_KALMAN) $(LIB_KP)
	$(CXX) $^ $(CXXFLAGS_GPU) $(LDFLAGS_GPU) -o $@

clean:
	rm -f kalman_GPU kalman_CPU

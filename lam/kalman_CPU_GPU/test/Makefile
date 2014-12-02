
CXX=g++


SRC_CPP = kp_init4matio.cpp  sp_parmap.cpp

CARMA_INC =  $(COMPASS_ROOT_DIR)/libcarma/include.h

ifeq ($(CXX),icc)
	CXXFLAGS_MKL = -mkl
else
	CXXFLAGS_MKL = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
endif

CXXFLAGS_CPU = -m64 -O2 -Wall $(CXXFLAGS_MKL) $(LIBS) -I.. -g -I$(CUDA_INC_PATH) -I$(CARMA_INC)
ifeq ($(COMPILATION_LAM),double)
  CXXFLAGS_CPU += -DKP_DOUBLE
endif

CXXFLAGS_GPU = $(CXXFLAGS_CPU) -I$(CUDA_INC_PATH)  -I$(MAGMA_ROOT)/include -DHAVE_CUBLAS

LDFLAGS_CPU =  -L$(CUDA_LIB_PATH_64) -lmatio -lgsl -lcudart -lcublas -lcusparse -L$(MAGMA_ROOT)/lib -lmagma
LDFLAGS_GPU = $(LDFLAGS_CPU) -L$(MAGMA_ROOT)/lib -lcublas -lcusparse -lcudart -lmagma


all:kalman_CPU kalman_GPU





kalman_CPU: kalman_CPU.cpp $(SRC_CPP) ../libkp_kalman.a
	$(CXX) $(SRC_CPP) $< ../libkp_kalman.a $(CXXFLAGS_CPU) $(LDFLAGS_CPU) -o $@

kalman_GPU: kalman_GPU.cpp $(SRC_CPP) ../libkp_kalman.a
	$(CXX) $(SRC_CPP) $< ../libkp_kalman.a $(CXXFLAGS_GPU) $(LDFLAGS_GPU) -o $@

clean:
	rm -f kalman_GPU kalman_CPU





#LIBS= -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm \
     -lmatio -lgsl -lcublas -lcusparse -lcudart -lmagma -L$(MAGMA_ROOT)/lib #-L/data/local/intel/composer_xe_2013_sp1.1.106/composer_xe_2013_sp1.1.106/mkl/lib/intel64 
#INC=-I../ -I$(MAGMA_ROOT)/include -DHAVE_CUBLAS -fPIC
#COMP=gcc -m64 -O2 $(LIBS) $(INC)
#LLIBS=kp_init4matio.cpp  kp_multif_phasetur.cpp sp_parmap.cpp

#all:kalman_CPU kalman_GPU


#kalman_CPU: kalman_CPU.cpp
#	$(COMP) $(LLIBS) kalman_CPU.cpp -o kalman_CPU ../libkp_kalman.a

#kalman_GPU: kalman_GPU.cpp
#	$(COMP) $(LLIBS) kalman_GPU.cpp -o kalman_GPU ../lib_kp_kalman.a -lmagma
	
#clean:
#	rm -f kalman_GPU kalman_CPU

CC = gcc
NVCC = nvcc

LAM_DIR = $(COMPASS_ROOT_DIR)/lam
LAM_KP_CLASSES = $(LAM_DIR)/kp_classes
LAM_KALMAN_DIR = $(LAM_DIR)/kalman_CPU_GPU
LAM_COMPASS_INTERFACE_DIR = $(LAM_DIR)/compass_interface

LIB_KP_CLASSES = $(LAM_KP_CLASSES)/lib_kp_classes.a
LIB_KALMAN = $(LAM_KALMAN_DIR)/lib_kalman.a


CFLAGS = -O3 -Wall -m64 -fPIC -fopenmp  -DMAGMA_SETAFFINITY  -DUSE_MAGMA_1_4 -DADD_ -DUSE_MAGMA_WITH_MKL -DHAVE_CUBLAS -I$(MAGMA_ROOT)/include -I$(LAM_KP_CLASSES) -I$(LAM_COMPASS_INTERFACE_DIR) -I$(COMPASS_ROOT_DIR)/libcarma/include.h -I$(COMPASS_ROOT_DIR)/libsutra/include.h -I$(CUDA_INC_PATH)

CUFLAGS = -O3 -m64 -Xcompiler -fPIC -use_fast_math -gencode $(GENCODE) -DHAVE_CUBLAS -I$(CUDA_INC_PATH) -I $(COMPASS_ROOT_DIR)/libcarma/include.h -I$(COMPASS_ROOT_DIR)/libsutra/include.h -I$(MAGMA_ROOT)/include 

ifeq ($(filter $(COMPILATION_LAM),double),double)
  CUFLAGS += -DKP_DOUBLE
  CFLAGS += -DKP_DOUBLE
endif


CFLAGS_KALMAN = $(CFLAGS) -I$(LAM_KALMAN_DIR) 
ifeq ($(filter $(COMPILATION_LAM),standalone),standalone)
  CFLAGS_KALMAN += -D__STANDALONE__
endif
ifeq ($(filter $(COMPILATION_LAM),temps_kalman),temps_kalman)
  CFLAGS_KALMAN += -D__TEMPS_BOUCLE_KALMAN__
endif

SRC_KP_CPP := $(wildcard $(LAM_KP_CLASSES)/*.cpp)
SRC_KP_CU := $(wildcard $(LAM_KP_CLASSES)/*.cu)
SRC_KP_HPP := $(wildcard $(LAM_KP_CLASSES)/*.hpp)
SRC_KP_HU := $(wildcard $(LAM_KP_CLASSES)/*.hu)



SRC_KALMAN_CPP := $(wildcard $(LAM_KALMAN_DIR)/*.cpp)



OBJ_KP_CU := $(patsubst $(LAM_KP_CLASSES)/%.cu,$(LAM_KP_CLASSES)/%.cuo,$(SRC_KP_CU))
OBJ_KP_CPP := $(patsubst $(LAM_KP_CLASSES)/%.cpp,$(LAM_KP_CLASSES)/%.o,$(SRC_KP_CPP))
OBJ_KALMAN_CPP := $(patsubst $(LAM_KALMAN_DIR)/%.cpp,$(LAM_KALMAN_DIR)/%.o,$(SRC_KALMAN_CPP))

OBJ_CU := $(OBJ_KP_CU)
OBJ_CPP := $(OBJ_KP_CPP) $(OBJ_KALMAN_CPP)

LDFLAGS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -lcublas -lcudart -lstdc++ -lgfortran -lz 

all:  $(LIB_KALMAN)
kalman: $(LIB_KALMAN)



$(LAM_KP_CLASSES)/%.cuo:  $(LAM_KP_CLASSES)/%.cu $(SRC_KP_HU)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(NVCC) $(CUFLAGS) $(LDFLAGS) -c -o $@ $<
	
$(LAM_KP_CLASSES)/%.o: $(LAM_KP_CLASSES)/%.cpp $(SRC_KP_HPP)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(CC) $(CFLAGS) $(LDFLAGS) -c -o $@ $< 

$(LIB_KP_CLASSES): $(OBJ_KP_CPP) $(OBJ_KP_CU) $(SRC_KP_HPP)
	@ar -rv $@ $(OBJ_KP_CPP) $(OBJ_KP_CU)


$(LAM_KALMAN_DIR)/%.o: $(LAM_KALMAN_DIR)/%.cpp
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(CC) $(CFLAGS_KALMAN) $(LDFLAGS) -c -o $@ $^

$(LIB_KALMAN): $(LIB_KP_CLASSES) $(OBJ_KALMAN_CPP)
	@ar -rv $@ $^



clean:
	@rm -f $(OBJ_KP_CU) $(OBJ_KP_CPP) $(LIB_KP_CLASSES) $(OBJ_KALMAN_CPP) $(LIB_KALMAN)


cleankalman:
	@rm -f $(OBJ_KP_CU) $(OBJ_KP_CPP) $(LIB_KP_CLASSES) $(OBJ_KALMAN_CPP) $(LIB_KALMAN)

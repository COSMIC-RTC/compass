# CUDA default definitions
export CUDA_ROOT?=/usr/local/cuda
export CUDA_INC_PATH?=$(CUDA_ROOT)/include
export CUDA_LIB_PATH?=$(CUDA_ROOT)/lib
export CUDA_LIB_PATH_64?=$(CUDA_ROOT)/lib64

# CULA default definitions
export CULA_ROOT?=/usr/local/cula
export CULA_INC_PATH?=$(CULA_ROOT)/include
export CULA_LIB_PATH?=$(CULA_ROOT)/lib
export CULA_LIB_PATH_64?=$(CULA_ROOT)/lib64

# YOGA default definitions
export COMPASS_ROOT_DIR?=/home/brujo/compass
export YOGA_DIR?=$COMPASS_ROOT_DIR/yoga_ao
export YOGA_AO_DIR?=$(YOGA_DIR)/yoga_ao
export YOGA_AO_TOP?=$(YOGA_AO_DIR)

export GENCODE?="arch=compute_12,code=sm_12"

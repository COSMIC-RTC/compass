
function check_var() {
    if [[ "$2"!="a" ]]
    then
	echo $1 defined to $2
    else
	echo $1 not defined, set to $3
	export $1=$3
    fi
}

# CUDA default definitions
check_var CUDA_ROOT a$CUDA_ROOT /usr/local/cuda
export CUDA_ROOT=/usr/local/cuda


#export CUDA_ROOT?=/usr/local/cuda
#export CUDA_INC_PATH?=$(CUDA_ROOT)/include
#export CUDA_LIB_PATH?=$(CUDA_ROOT)/lib
#export CUDA_LIB_PATH_64?=$(CUDA_ROOT)/lib64

# CULA default definitions
#export CULA_ROOT?=/usr/local/cula
#export CULA_INC_PATH?=$(CULA_ROOT)/include
#export CULA_LIB_PATH?=$(CULA_ROOT)/lib
#export CULA_LIB_PATH_64?=$(CULA_ROOT)/lib64

# YOGA default definitions
#export YOGA_DIR?=/obs/sevin/compass
#export YOGA_AO_DIR?=$(YOGA_DIR)/yoga_ao
#export YOGA_AO_TOP?=$(YOGA_AO_DIR)


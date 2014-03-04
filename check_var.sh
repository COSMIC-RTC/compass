#!/bin/bash

function check_var() {
    if [[ $2 != "" ]]
    then
	echo $1 defined to $2
    else
	echo $1 not defined !!!
    fi
}

# CUDA default definitions
check_var CUDA_ROOT $CUDA_ROOT
check_var CUDA_INC_PATH $CUDA_INC_PATH
check_var CUDA_LIB_PATH $CUDA_LIB_PATH
check_var CUDA_LIB_PATH_64 $CUDA_LIB_PATH_64

# MAGMA default definitions
check_var MAGMA_ROOT $MAGMA_ROOT

# CULA default definitions
check_var CULA_ROOT $CULA_ROOT
check_var CULA_INC_PATH $CULA_INC_PATH
check_var CULA_LIB_PATH $CULA_LIB_PATH
check_var CULA_LIB_PATH_64 $CULA_LIB_PATH_64

# YOGA default definitions
check_var YOGA_DIR $YOGA_DIR
check_var YOGA_AO_DIR $YOGA_AO_DIR
check_var YOGA_AO_TOP $YOGA_AO_TOP

